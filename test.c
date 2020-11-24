/**************************
* Kryzstof Kudlak         *
* Image processing part 2 *
* Finished 11/24/2020     *
**************************/

#include <stdio.h> // standard lib
#include <unistd.h> // for forks and pipes
#include <sys/wait.h> // for waiting for children
#include <errno.h> // also for waiting for children
#include <stdlib.h> // not certain
#include <math.h> // for calculating standard deviation
#include "iplib2New.c" // for connecting to library functions for images

/* mask matrices to be used in edge detection */
int horz_mask[3][3] = { { -1, -2, -1 }, {  0,  0,  0 }, {  1,  2,  1 } };
int vert_mask[3][3] = { { -1,  0,  1 }, { -2,  0,  2 }, { -1,  0,  1 } };

/* function definitions as included in iplib2 */
image_ptr read_pnm(char *filename, int *rows, int *cols, int *type);
int getnum(FILE *fp);
void write_pnm(image_ptr ptr, char *filename, int rows, int cols, int type);

/* function definitons from this file */
int* getnxn(image_ptr, int, int, int, int);
int prod(int*, int mask[3][3]);
double threshold(image_ptr, int, int);
double mean(image_ptr, int, int);
double stdev(image_ptr, int);

int cmpfunc(const void * a, const void * b) { // comparison function for qsort used later, wasn't sure where else to put it
	return ( *(int*)a < *(int*)b );
}


int main(int argc, char **argv) {


	/***
	* References to be used later
	***/
	int rows, cols, type;      // metadata about the image being read
	image_ptr origin;          // original image to be read from file
	
	image_ptr himg, vimg, img; // edge detection image pointers
	image_ptr hbin, vbin, bin; // binarized image pointers



	/***
	* Handle command line args
	***/
	char *file_names[8] = {    // default naming if no args specified
		"a.exe",               // executeable file
		"in.pgm",              // input file
		"hout.pgm",            // output for horizontal edge detection
		"vout.pgm",            // output for vertical edge detection
		"out.pgm",             // output for combined
		"hbin.pgm",            // horizontal binary threshold
		"vbin.pgm",            // vertical binary threshold
		"bin.pgm"              // combined binary threshold image
	};

	/* basic error checking to make sure we have the right # of params */
	if (argc > 8) {
		printf("ERR: too many parameters");
		exit(1);
	}

	/* override defaults for each in/out file specified */
	for (int i = 0; i < argc; i++)
		file_names[i] = argv[i];

	/* read the image to initialize the origin image pointer */
	printf("reading input image ... \n");
	origin = read_pnm(file_names[1], &rows, &cols, &type);
	printf("image read successfully \n");
	printf("rows=%d, cols=%d, type=%d \n", rows, cols, type);

	/***
	* Conclude handling of command line args
	***/



	/***
	* Create pipes and processes
	***/
	int proc_num = 0;
	while (proc_num < 1) {
		printf("Input desired number of processes: \n");
		scanf("%d", &proc_num);
	}

	int starting_index = 0;          // each process will start reading origin from a certain index
	int size = (rows*cols)/proc_num; // this is the # of elements to be processed by each process, same for each child


	// creates a pipe for each process
	int hpipe[proc_num][2]; // for horizontal response
	int vpipe[proc_num][2]; // for vertical response
	int tpipe[proc_num][2]; // for combined response
	int apipe[proc_num][2]; // for the local means to be sent
	for (int i = 0; i < proc_num; i++) {
		if (pipe(hpipe[i]) < 0 || pipe(vpipe[i]) < 0 || pipe(tpipe[i]) < 0 || pipe(apipe[i]) < 0) // where pipe creation occurs
			return 1; // if it errors out, exit program
	}


	// forking a bunch of children
	int pid = fork();
	int pipe_id = 0;
	for (int i = 1; i < proc_num; i++) {
		if (pid != 0) {
			pid = fork();
			pipe_id = i;
			starting_index += size; // sets starting indexes of processes, different for each child
		}
	}


	// close the unused ends of each pipe
	// 0 is read, 1 is write
	// children write, parent reads
	for (int i = 0; i < proc_num; i++) {
		// for children
		if (pid == 0) {
			close(hpipe[i][0]); // close the read channel of every pipe
			close(vpipe[i][0]);
			close(tpipe[i][0]);
			close(apipe[i][0]);
			if (i != pipe_id) {
				close(hpipe[i][1]); // close the write channel if this isnt the appropriate pipe id
				close(vpipe[i][1]); 
				close(tpipe[i][1]); 
				close(apipe[i][1]); 
			}
		}
		// for the parent
		else {
			close(hpipe[i][1]); // close the writing channel for every pipe, it reads from every process
			close(vpipe[i][1]); 
			close(tpipe[i][1]); 
			close(apipe[i][1]); 
		}

	}
	/***
	* K the pipes and process are created
	***/


	/* process edge detection in each child process & send to parent */ 
	if (pid == 0) {

		himg = (image_ptr)malloc(size*sizeof(unsigned char));
		vimg = (image_ptr)malloc(size*sizeof(unsigned char));
		 img = (image_ptr)malloc(size*sizeof(unsigned char));

		if (himg == NULL || vimg == NULL || img == NULL) {
			printf("\nUnable to allocate memory in process %d\n", pipe_id);
			return 2;
		}

		// edge detection
		for (int i = starting_index; i < starting_index + size; i++) {

			int index = i - starting_index;
			int *nxn = getnxn(origin, rows, cols, i/rows, i%rows);
			himg[index] = prod(nxn, horz_mask);
			vimg[index] = prod(nxn, vert_mask);
			free(nxn);
			 img[index] = himg[index] + vimg[index];

		}

		// send responses through the pipe
		if (write(hpipe[pipe_id][1], himg, sizeof(unsigned char)*size) < 0){
			printf("ERR: Process %d failed to write himg to pipe\n", pipe_id);
			return 1;
		}

		if (write(vpipe[pipe_id][1], vimg, sizeof(unsigned char)*size) < 0){
			printf("ERR: Process %d failed to write vimg to pipe\n", pipe_id);
			return 1;
		}

		if (write(tpipe[pipe_id][1],  img, sizeof(unsigned char)*size) < 0){
			printf("ERR: Process %d failed to write img to pipe\n", pipe_id);
			return 1;
		}

		double local_mean = mean(origin, starting_index, size);
		if (write(apipe[pipe_id][1],  &local_mean, sizeof(double)) < 0){
			printf("ERR: Process %d failed to write img to pipe\n", pipe_id);
			return 1;
		}

		printf("PROCESS: %d SUCCESSFULLY WROTE TO PIPE\n", pipe_id);
	}

	//  wait for all the parent's children to finish
	while (wait(NULL) != -1 || errno != ECHILD) {}


	/***
	* MAGIC TO RECEIVE MESSAGES FROM PIPES
	***/
	double mean_of_means;
	double stdev_of_means;
	double median_of_medians; // variables to hold statistics for output later

	if (pid != 0) {

		printf("\n\nParent reading in from children ... \n");

		himg = (image_ptr) malloc(rows*cols*(sizeof(unsigned char)));
		vimg = (image_ptr) malloc(rows*cols*(sizeof(unsigned char)));
		img  = (image_ptr) malloc(rows*cols*(sizeof(unsigned char)));

		if (himg == NULL || vimg == NULL || img == NULL)	{
			printf("Unable to allocate memory for image pointers in parent process");
			exit(1);
		}

		double means[proc_num];

		for (int i = 0; i < proc_num; i++) {

			image_ptr htemp = malloc(size*(sizeof(unsigned char)));
			image_ptr vtemp = malloc(size*(sizeof(unsigned char)));
			image_ptr ttemp = malloc(size*(sizeof(unsigned char)));

			if (read(hpipe[i][0], htemp, size*(sizeof(unsigned char))) < 0) { // error checking
				printf("Parent failed to receive data from process %d \n", i);
				return 2;
			}

			if (read(vpipe[i][0], vtemp, size*(sizeof(unsigned char))) < 0) { // error checking
				printf("Parent failed to receive data from process %d \n", i);
				return 2;
			}

			if (read(tpipe[i][0], ttemp, size*(sizeof(unsigned char))) < 0) { // error checking
				printf("Parent failed to receive data from process %d \n", i);
				return 2;
			}

			double temp;
			if (read(apipe[i][0], &temp, sizeof(double)) < 0) { // error checking
				printf("Parent failed to receive average from process %d \n", i);
				return 2;
			}
			means[i] = temp;

			// piece together images
			for (int j = 0; j < size; j++) {
				himg[i*size + j] = htemp[j];
				vimg[i*size + j] = vtemp[j];
				 img[i*size + j] = ttemp[j];
			}

		}

		// figure out the average of the means
		int sum = 0;
		for (int i = 0; i < proc_num; i++) {
			sum += means[i];
		}
		mean_of_means = (double)(((double)sum)/((double)proc_num));


		// figure out std deviation of means
		int stdsum = 0;
		for (int i = 0; i < proc_num; i++)
			stdsum += ((img[i] - mean_of_means) * (img[i] - mean_of_means));

		stdev_of_means = sqrt(((double)stdsum)/((double)(proc_num-1)));

		// bubble sort the means (later i learned that there's a qsort function for this)
		for (int i = 0; i < proc_num; i++) {
			for (int j = 0; j < proc_num-i-1; j++) {
				if (means[j] > means[j+1]) {
					double temp = means[j];
					means[j] = means[j+1];
					means[j+1] = temp;
				}
			}
		}

		// get their median
		if (proc_num%2 == 1)
			median_of_medians = means[proc_num/2]; // integer division to get to middle of list
		else 
			median_of_medians = (means[proc_num/2] + means[proc_num/2-1])/2; 
		
	}
	/***
	* PIPE READING/DATA PROCESSING IS COMPLETE FOR THE PARENT PROCESS
	***/


	// close the remaining pipe channels
	if (pid == 0) { // children
		close(hpipe[pipe_id][1]);
		close(vpipe[pipe_id][1]);
		close(tpipe[pipe_id][1]);
		close(apipe[pipe_id][1]);
		return 0; // kill the children processes, F in chat
	}
	else { // parent has to close all its receiving channels
		for (int i = 1; i < proc_num; i++) {
			close(hpipe[i][0]);
			close(vpipe[i][0]);
			close(tpipe[i][0]);
			close(apipe[i][0]); // notice the parent stays alive to finish
		}
	}


	/**
	* This part is mostly old code from assignment 2 with the statistics added at the end
	**/

	// calculate threshold values
	printf("\ndetermining threshold values ... \n");
	double th = threshold(himg, rows, cols);
	double tv = threshold(vimg, rows, cols);
	double tt = threshold(img , rows, cols);
	printf("threshold values determined \n");



	/** OG Code to binarize and write **/

	/* Binarize the images */
	printf("\nbinarizing images ... \n");
	hbin = (image_ptr) malloc(rows*cols*(sizeof(unsigned char)));
	vbin = (image_ptr) malloc(rows*cols*(sizeof(unsigned char)));
	bin  = (image_ptr) malloc(rows*cols*(sizeof(unsigned char)));

	if(hbin == NULL || vbin == NULL || bin == NULL)	{
		printf("Unable to allocate memory for binarized image pointers");
		exit(1);
	}

	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < cols; j++) {
			hbin[i*cols + j] = (himg[i*cols + j] > th ? 255 : 0);
			vbin[i*cols + j] = (vimg[i*cols + j] > tv ? 255 : 0);
			 bin[i*cols + j] = ( img[i*cols + j] > tt ? 255 : 0);
		}
	}
	printf("binarizing complete \n");

	/* write all images into their respective files */
	printf("\nwriting images ... \n");
	write_pnm(himg, file_names[2], rows, cols, type);
	write_pnm(vimg, file_names[3], rows, cols, type);
	write_pnm( img, file_names[4], rows, cols, type);
	write_pnm(hbin, file_names[5], rows, cols, type);
	write_pnm(vbin, file_names[6], rows, cols, type);
	write_pnm( bin, file_names[7], rows, cols, type);
	printf("\nsuccess! \n");

	// calculating median without a separate function because im lazy
	size = rows*cols-proc_num-1;
	image_ptr sort = (image_ptr) malloc(size*(sizeof(unsigned char)));
	for (int i = 0; i < size; i++)
		sort[i] = origin[i]; // copy img into sort so that the original isn't mutated
	qsort(sort, size, sizeof(image_ptr), cmpfunc); // you have no idea how happy i am that this exists

	double median;
	if (size%2 == 1)
		median = (double)sort[size/2]; // integer division to get to middle of list
	else 
		median = (double)(((double)sort[size/2] + (double)sort[size/2-1])/2); // when the list is even, take average of 2 middle elements

	// ah finally, statistics can be printed
	printf("\n\nMetrics as required in the project handout:\n");
	printf("\nGlobal Image Mean: %lf\n",                mean(origin, 0, rows*cols));
	printf("Global Standard Deviation: %lf\n",         stdev(origin,    rows*cols));
	printf("Global Median: %lf\n", median);

	printf("\nHorizontal Edge Response Average: %lf\n", mean(himg  , 0, rows*cols));
	printf("Vertical Edge Response Average: %lf\n",     mean(vimg  , 0, rows*cols));
	printf("Combined Edge Response Average: %lf\n",     mean( img  , 0, rows*cols));

	printf("\nOverall image mean: %lf\n",               mean_of_means);
	printf("Overall image standard deviation: %lf\n",   stdev_of_means);
	printf("Overall image median: %lf\n\n\n",           median_of_medians);

	// levi if you're reading this, I'm sorry for making my program so horribly monolithic

	return 0;
	
}


/*********************************************************************
* Description: this method returns the 9x9 matrix at the given point *
*         - origin: ptr to the image being processed                 *
*         - rows: number of rows in the image                        *
*         - cols: number of cols in the image                        *
*         - x: x coord of the point for the matrix being found       *
*         - y: y coord of the point for the matrix being found       *
**********************************************************************/ 
int* getnxn(image_ptr origin, int rows, int cols, int x, int y) {
	int* ptr = malloc(sizeof(int)*9);

	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			if ( (x - 1 < 0 && i == 0) || (x + 1 >= cols && i == 2)  || (y - 1 < 0 && j == 0) || (y + 1 >= rows && j == 2) )
				ptr[i*3 + j] = origin[x*cols + y]; // if we have edges out of bounds, set their value = to middle val
			else
				ptr[i*3 + j] = origin[(x+(i-1))*cols + (y+(j-1))];  // base case, nxn matrix is well within its bounds
 		}
	}

	return ptr;
}

/*********************************************************************
* Description: returns the product of two 3x3 matrices               *
*         - nxn: the 3x3 array of the target element & its neighbors *
*         - mask: the 3x3 array representing the mask                *
**********************************************************************/ 
int prod(int* nxn, int mask[3][3]) {
	int sum = 0;
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			sum += nxn[i*3 + j] * mask[i][j];
	return (int)abs((double)sum);
}


/*********************************************************************
* Description: returns the threshold calculation for an image        *
*         - img: the image for the threshold to be calculated from   *
*         - rows: the number of rows in the image                    *
*         - cols: the number of cols in the image                    *
**********************************************************************/
double threshold(image_ptr img, int rows, int cols) {

	double m = mean(img, 0, rows*cols);
	double std = stdev(img, rows*cols);

	return (m + std);
}


/*********************************************************************
* Description: returns the mean calculation for an image             *
*         - img: the image for the mean to be calculated from        *
*         - size: the length of the given image array                *
**********************************************************************/
double mean(image_ptr img, int starting_index, int size) {
	int sum = 0;
	for (int i = starting_index; i < starting_index+size; i++)
		sum += (img[i] > 255 ? 255 : img[i]);

	return (double)( ((double)sum)/((double)size) );
}

/*********************************************************************
* Description: returns the stdev calculation for an image            *
*         - img: the image for the stdev to be calculated from       *
*         - size: the length of the given image array                *
**********************************************************************/
double stdev(image_ptr img, int size) {

	double m = mean(img, 0, size);
	int stdsum = 0;

	for (int i = 0; i < size; i++) {
		int elem = (img[i] > 255 ? 255 : img[i]);
		stdsum += ((elem - m) * (elem - m));
	}

	return sqrt(((double)stdsum)/((double)(size-1)));
}

