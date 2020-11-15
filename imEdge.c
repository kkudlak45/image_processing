/*******************
* Kryzstof Kudlak  *
* Image Processing *
* 10/27/2020       *
********************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "iplib2New.c"

/* mask matrixes to be used in edge detection */
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

int main(int argc, char **argv)	{

	int rows, cols, type;      // metadata about the image being read
	image_ptr origin;          // original image to be read from file
	
	image_ptr himg, vimg, img; // edge detection image pointers
	image_ptr hbin, vbin, bin; // binarized image pointers

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

	/* memory allocation for each new image to be created */
	himg = (image_ptr) malloc(rows*cols*(sizeof(unsigned char)));
	vimg = (image_ptr) malloc(rows*cols*(sizeof(unsigned char)));
	img  = (image_ptr) malloc(rows*cols*(sizeof(unsigned char)));
	hbin = (image_ptr) malloc(rows*cols*(sizeof(unsigned char)));
	vbin = (image_ptr) malloc(rows*cols*(sizeof(unsigned char)));
	bin  = (image_ptr) malloc(rows*cols*(sizeof(unsigned char)));

	if(himg == NULL || vimg == NULL || img == NULL
		 || hbin == NULL || vbin == NULL || bin == NULL)	{
		printf("Unable to allocate memory for image pointers");
		exit(1);
	}

	/* process edge detection */ 
	printf("\nprocessing edge detection ... \n");
	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < cols; j++) {
			int *nxn = getnxn(origin, rows, cols, i, j);
			himg[i*cols + j] = prod(nxn, horz_mask);
			vimg[i*cols + j] = prod(nxn, vert_mask);
			free(nxn);
			img[i*cols + j]  = himg[i*cols + j] + vimg[i*cols + j];
 		}
	}
	printf("edge detection processing complete \n");

	free(origin); // original image reference is no longer necessary

	/* determine threshold calculations */
	printf("\ndetermining threshold values ... \n");
	double th = threshold(himg, rows, cols);
	double tv = threshold(vimg, rows, cols);
	double tt = threshold(img , rows, cols);
	printf("threshold values determined \n");

	/* Binarize the images */
	printf("\nbinarizing images ... \n");
	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < cols; j++) {
			hbin[i*cols + j] = (himg[i*cols + j] > th ? 255 : 0);
			vbin[i*cols + j] = (vimg[i*cols + j] > tv ? 255 : 0);
			 bin[i*cols + j] = ( img[i*cols + j] > tt ? 255 : 0);
		}
	}
	printf("binarizing complete \n");

	/* write all images into their respective files */
	printf("\nwriting image ... \n");
	write_pnm(himg, file_names[2], rows, cols, type);
	write_pnm(vimg, file_names[3], rows, cols, type);
	write_pnm(img , file_names[4], rows, cols, type);
	write_pnm(hbin, file_names[5], rows, cols, type);
	write_pnm(vbin, file_names[6], rows, cols, type);
	write_pnm(bin , file_names[7], rows, cols, type);
	printf("\nsuccess! \n");

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

	/* calculate mean */
	int sum = 0;
	for (int i = 0; i < rows; i++)
		for (int j = 0; j < cols; j++)
			sum += img[i*rows + j];
	double mean = sum/(rows*cols);

	/* calculate standard deviation */
	int stdsum;
	for (int i = 0; i < rows; i++)
		for (int j = 0; j < cols; j++)
			stdsum += ((img[i*rows + j] - mean) * (img[i*rows + j] - mean));
	double stdev = sqrt(((double)stdsum/((double)(rows*cols-1))));

	printf("Mean: %lf | Std Dev: %lf\n", mean, stdev);

	return (mean+stdev);
}
