/* source: http://marathon.csee.usf.edu/edge/edge_detection.html */
/* URL: ftp://figment.csee.usf.edu/pub/Edge_Comparison/source_code/canny.src */

/* ECPS 203 Assignment 4 solution */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

#define VERBOSE 0

#define NOEDGE 255
#define POSSIBLE_EDGE 128
#define EDGE 0
#define BOOSTBLURFACTOR 90.0
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#define SIGMA 0.6
#define TLOW  0.3
#define THIGH 0.8

#define COLS 2704
#define ROWS 1520
#define SIZE COLS*ROWS
#define VIDEONAME "Engineering"
#define IMG_IN    "video/" VIDEONAME "%03d.pgm"
#define IMG_OUT   VIDEONAME "%03d_edges.pgm"
#define IMG_NUM   30 /* number of images processed (1 or more) */
#define AVAIL_IMG 30 /* number of different image frames (1 or more) */

/* upper bound for the size of the gaussian kernel
 * SIGMA must be less than 4.0
 * check for 'windowsize' below
 */
int mag_times = 0;
double mag_totaltime = 0.0;
int der_times = 0;
double der_totaltime = 0.0;
int smo_times = 0;
double smo_totaltime = 0.0;
double x_totaltime = 0.0;
double y_totaltime = 0.0;
int ker_times = 0;
double ker_totaltime = 0.0;
int fol_times = 0;
double fol_totaltime = 0.0;
int app_times = 0;
double app_totaltime = 0.0;
int non_times = 0;
double non_totaltime = 0.0;
int read_times = 0;
double read_totaltime = 0.0;
int write_times = 0;
double write_totaltime = 0.0;
int can_times = 0;
double can_totaltime = 0.0;

#define WINSIZE 21

/* Function Declarations */
int read_pgm_image(const char *infilename, unsigned char *image, int rows, int cols);
int write_pgm_image(const char *outfilename, unsigned char *image, int rows,
    int cols, const char *comment, int maxval);

void canny(unsigned char *image, int rows, int cols, float sigma,
         float tlow, float thigh, unsigned char *edge);
void gaussian_smooth(unsigned char *image, int rows, int cols, float sigma,
        short int *smoothedim);
void make_gaussian_kernel(float sigma, float *kernel, int *windowsize);
void derrivative_x_y(short int *smoothedim, int rows, int cols,
        short int *delta_x, short int *delta_y);
void magnitude_x_y(short int *delta_x, short int *delta_y, int rows, int cols,
        short int *magnitude);
void apply_hysteresis(short int *mag, unsigned char *nms, int rows, int cols,
        float tlow, float thigh, unsigned char *edge);

void non_max_supp(short *mag, short *gradx, short *grady, int nrows, int ncols,
    unsigned char *result);

int main(int argc, char* argv[])
{
   unsigned char image[SIZE];
   unsigned char edge[SIZE];
   char infilename[70];
   char outfilename[128];    /* Name of the output "edge" image */
   int i=0, n=0;

   for(i=0; i<IMG_NUM; i++)
   {
      n = i % AVAIL_IMG;
      sprintf(infilename, IMG_IN, n+1);
      /****************************************************************************
      * Read in the image. This read function allocates memory for the image.
      ****************************************************************************/
      if(VERBOSE) printf("Reading the image %s.\n", infilename);
      if(read_pgm_image(infilename, image, ROWS, COLS) == 0){
         fprintf(stderr, "Error reading the input image, %s.\n", infilename);
         exit(1);
      }

      /****************************************************************************
      * Perform the edge detection. All of the work takes place here.
      ****************************************************************************/
      if(VERBOSE) printf("Starting Canny edge detection.\n");
      canny(image, ROWS, COLS, SIGMA, TLOW, THIGH, edge);

      /****************************************************************************
      * Write out the edge image to a file.
      ****************************************************************************/
      n = i % AVAIL_IMG;
      sprintf(outfilename, IMG_OUT, n+1);
      if(VERBOSE) printf("Writing the edge iname in the file %s.\n", outfilename);
      if(write_pgm_image(outfilename, edge, ROWS, COLS, "", 255) == 0){
         fprintf(stderr, "Error writing the edge image, %s.\n", outfilename);
         exit(1);
      }
   }
   //printf("canny run time %f secs\n", can_totaltime / can_times); 
   //printf("canny runs %d times\n", can_times);
   printf("magnitude_x_y run time %f secs\n", mag_totaltime / mag_times); 
   printf("mag runs %d times\n", mag_times);
   printf("derrivative_x_y run time %f secs\n", der_totaltime / der_times); 
   printf("der runs %d times\n", der_times);
   //printf("gaussian_smooth run time %f secs\n", smo_totaltime / smo_times); 
   //printf("gau runs %d times\n", smo_times);
   printf("blur x run time %f secs\n", x_totaltime / smo_times); 
   printf("blur y run time %f secs\n", y_totaltime / smo_times); 
   printf("make_gaussian_kernel run time %f secs\n", ker_totaltime / ker_times); 
   printf("ker runs %d times\n", ker_times);
   printf("follow_edges run time %f secs\n", fol_totaltime / fol_times); 
   printf("fol runs %d times\n", fol_times);
   printf("apply_hysteresis run time %f secs\n", app_totaltime / app_times); 
   printf("app runs %d times\n", app_times);
   printf("non_max_supp run time %f secs\n", non_totaltime / non_times); 
   printf("non runs %d times\n", non_times);
   printf("read_pgm_image run time %f secs\n", read_totaltime / read_times); 
   printf("read runs %d times\n", read_times);
   printf("write_pgm_image run time %f secs\n", write_totaltime / write_times); 
   printf("write runs %d times\n", write_times);

   return(0); /* exit cleanly */
}

/*******************************************************************************
* PROCEDURE: canny
* PURPOSE: To perform canny edge detection.
* NAME: Mike Heath
* DATE: 2/15/96
*******************************************************************************/
void canny(unsigned char *image, int rows, int cols, float sigma,
         float tlow, float thigh, unsigned char *edge)
{
   can_times++;
   //clock_t start = clock();
   unsigned char nms[SIZE]    /* Points that are local maximal magnitude. */
			= {0};
   short int smoothedim[SIZE] /* The image after gaussian smoothing.      */
			= {0},
             delta_x[SIZE]    /* The first devivative image, x-direction. */
			= {0},
             delta_y[SIZE]    /* The first derivative image, y-direction. */
			= {0},
             magnitude[SIZE]  /* The magnitude of the gadient image.      */
			= {0};

   /****************************************************************************
   * Perform gaussian smoothing on the image using the input standard
   * deviation.
   ****************************************************************************/
   if(VERBOSE) printf("Smoothing the image using a gaussian kernel.\n");
   gaussian_smooth(image, rows, cols, sigma, smoothedim);

   /****************************************************************************
   * Compute the first derivative in the x and y directions.
   ****************************************************************************/
   if(VERBOSE) printf("Computing the X and Y first derivatives.\n");
   clock_t start_der = clock();
   derrivative_x_y(smoothedim, rows, cols, delta_x, delta_y);
   clock_t finish_der = clock();
   der_totaltime += (double)(finish_der - start_der)/CLOCKS_PER_SEC;

   /****************************************************************************
   * Compute the magnitude of the gradient.
   ****************************************************************************/
   if(VERBOSE) printf("Computing the magnitude of the gradient.\n");
   clock_t start_mag = clock();
   magnitude_x_y(delta_x, delta_y, rows, cols, magnitude);
   clock_t finish_mag = clock();
   mag_totaltime += (double)(finish_mag - start_mag)/CLOCKS_PER_SEC; 

   /****************************************************************************
   * Perform non-maximal suppression.
   ****************************************************************************/
   if(VERBOSE) printf("Doing the non-maximal suppression.\n");
   clock_t start_non = clock();
   non_max_supp(magnitude, delta_x, delta_y, rows, cols, nms);
   clock_t finish_non = clock();
   non_totaltime += (double)(finish_non - start_non)/CLOCKS_PER_SEC;

   /****************************************************************************
   * Use hysteresis to mark the edge pixels.
   ****************************************************************************/
   if(VERBOSE) printf("Doing hysteresis thresholding.\n");
   clock_t start_app = clock();
   apply_hysteresis(magnitude, nms, rows, cols, tlow, thigh, edge);
   clock_t finish_app = clock();
   app_totaltime += (double)(finish_app - start_app)/CLOCKS_PER_SEC;
   //can_totaltime += (double)(finish-start)/CLOCKS_PER_SEC; 
}

/*******************************************************************************
* PROCEDURE: magnitude_x_y
* PURPOSE: Compute the magnitude of the gradient. This is the square root of
* the sum of the squared derivative values.
* NAME: Mike Heath
* DATE: 2/15/96
*******************************************************************************/
void magnitude_x_y(short int *delta_x, short int *delta_y, int rows, int cols,
        short int *magnitude)
{
   mag_times++;
   //clock_t start = clock();
   int r, c, pos, sq1, sq2;

   for(r=0,pos=0;r<rows;r++){
      for(c=0;c<cols;c++,pos++){
         sq1 = (int)delta_x[pos] * (int)delta_x[pos];
         sq2 = (int)delta_y[pos] * (int)delta_y[pos];
         magnitude[pos] = (short)(0.5 + sqrt((float)sq1 + (float)sq2));
      }
   }
   //clock_t finish = clock();
}

/*******************************************************************************
* PROCEDURE: derrivative_x_y
* PURPOSE: Compute the first derivative of the image in both the x any y
* directions. The differential filters that are used are:
*
*                                          -1
*         dx =  -1 0 +1     and       dy =  0
*                                          +1
*
* NAME: Mike Heath
* DATE: 2/15/96
*******************************************************************************/

void derrivative_x_y(short int *smoothedim, int rows, int cols,
        short int *delta_x, short int *delta_y)
{
   der_times++;
   //clock_t start = clock();
   int r, c, pos;

   /****************************************************************************
   * Compute the x-derivative. Adjust the derivative at the borders to avoid
   * losing pixels.
   ****************************************************************************/
   if(VERBOSE) printf("   Computing the X-direction derivative.\n");
   for(r=0;r<rows;r++){
      pos = r * cols;
      delta_x[pos] = smoothedim[pos+1] - smoothedim[pos];
      pos++;
      for(c=1;c<(cols-1);c++,pos++){
         delta_x[pos] = smoothedim[pos+1] - smoothedim[pos-1];
      }
      delta_x[pos] = smoothedim[pos] - smoothedim[pos-1];
   }

   /****************************************************************************
   * Compute the y-derivative. Adjust the derivative at the borders to avoid
   * losing pixels.
   ****************************************************************************/
   if(VERBOSE) printf("   Computing the Y-direction derivative.\n");
   for(c=0;c<cols;c++){
      pos = c;
      delta_y[pos] = smoothedim[pos+cols] - smoothedim[pos];
      pos += cols;
      for(r=1;r<(rows-1);r++,pos+=cols){
         delta_y[pos] = smoothedim[pos+cols] - smoothedim[pos-cols];
      }
      delta_y[pos] = smoothedim[pos] - smoothedim[pos-cols];
   }
   //clock_t finish = clock();
   //der_totaltime += (double)(finish-start)/CLOCKS_PER_SEC; 
}

/*******************************************************************************
* PROCEDURE: gaussian_smooth
* PURPOSE: Blur an image with a gaussian filter.
* NAME: Mike Heath
* DATE: 2/15/96
*******************************************************************************/

void gaussian_smooth(unsigned char *image, int rows, int cols, float sigma,
        short int *smoothedim)
{
   smo_times++;
   //clock_t start = clock();
   int r, c, rr, cc,     /* Counter variables. */
      windowsize,        /* Dimension of the gaussian kernel. */
      center;            /* Half of the windowsize. */
   float tempim[SIZE]    /* Buffer for separable filter gaussian smoothing. */
		= {0.0},
         kernel[WINSIZE] /* A one dimensional gaussian kernel. */
		= {0.0},
         dot,            /* Dot product summing variable. */
         sum;            /* Sum of the kernel weights variable. */

   /****************************************************************************
   * Create a 1-dimensional gaussian smoothing kernel.
   ****************************************************************************/
   if(VERBOSE) printf("   Computing the gaussian smoothing kernel.\n");
   clock_t start_ker = clock();
   make_gaussian_kernel(sigma, kernel, &windowsize);
   clock_t finish_ker = clock();
   ker_totaltime += (double)(finish_ker - start_ker)/CLOCKS_PER_SEC;
   center = windowsize / 2; 

   /****************************************************************************
   * Blur in the x - direction.
   ****************************************************************************/
   clock_t start_x = clock();
   if(VERBOSE) printf("   Bluring the image in the X-direction.\n");
   for(r=0;r<rows;r++){
      for(c=0;c<cols;c++){
         dot = 0.0;
         sum = 0.0;
         for(cc=(-center);cc<=center;cc++){
            if(((c+cc) >= 0) && ((c+cc) < cols)){
               dot += (float)image[r*cols+(c+cc)] * kernel[center+cc];
               sum += kernel[center+cc];
            }
         }
         tempim[r*cols+c] = dot/sum;
      }
   }
   clock_t finish_x = clock();
   x_totaltime += (double)(finish_x - start_x)/CLOCKS_PER_SEC; 

   /****************************************************************************
   * Blur in the y - direction.
   ****************************************************************************/
   clock_t start_y = clock();
   if(VERBOSE) printf("   Bluring the image in the Y-direction.\n");
   for(c=0;c<cols;c++){
      for(r=0;r<rows;r++){
         sum = 0.0;
         dot = 0.0;
         for(rr=(-center);rr<=center;rr++){
            if(((r+rr) >= 0) && ((r+rr) < rows)){
               dot += tempim[(r+rr)*cols+c] * kernel[center+rr];
               sum += kernel[center+rr];
            }
         }
         smoothedim[r*cols+c] = (short int)(dot*BOOSTBLURFACTOR/sum + 0.5);
      }
   }
   clock_t finish_y = clock();
   y_totaltime += (double)(finish_y-start_y)/CLOCKS_PER_SEC;

   //clock_t finish = clock();
   //smo_totaltime += (double)(finish-start)/CLOCKS_PER_SEC; 
}

/*******************************************************************************
* PROCEDURE: make_gaussian_kernel
* PURPOSE: Create a one dimensional gaussian kernel.
* NAME: Mike Heath
* DATE: 2/15/96
*******************************************************************************/

void make_gaussian_kernel(float sigma, float *kernel, int *windowsize)
{
   ker_times++;
   //clock_t start = clock();
   int i, center;
   float x, fx, sum=0.0;

   *windowsize = 1 + 2 * ceil(2.5 * sigma);
   center = (*windowsize) / 2;

   if(VERBOSE) printf("      The kernel has %d elements.\n", *windowsize);

   for(i=0;i<(*windowsize);i++){
      x = (float)(i - center);
      fx = pow(2.71828, -0.5*x*x/(sigma*sigma)) / (sigma * sqrt(6.2831853));
      kernel[i] = fx;
      sum += fx;
   }

   for(i=0;i<(*windowsize);i++) kernel[i] /= sum;

   if(VERBOSE){
      printf("The filter coefficients are:\n");
      for(i=0;i<(*windowsize);i++)
         printf("kernel[%d] = %f\n", i, kernel[i]);
   }
   //clock_t finish = clock();
   //ker_totaltime += (double)(finish-start)/CLOCKS_PER_SEC; 
}

/*******************************************************************************
* PROCEDURE: follow_edges
* PURPOSE: This procedure edges is a recursive routine that traces edgs along
* all paths whose magnitude values remain above some specifyable lower
* threshhold.
* NAME: Mike Heath
* DATE: 2/15/96
*******************************************************************************/

void follow_edges(unsigned char *edgemapptr, short *edgemagptr, short lowval,
   int cols)
{
   fol_times++;
   //clock_t start = clock();
   short *tempmagptr;
   unsigned char *tempmapptr;
   int i;
   int x[8] = {1,1,0,-1,-1,-1,0,1},
       y[8] = {0,1,1,1,0,-1,-1,-1};

   for(i=0;i<8;i++){
      tempmapptr = edgemapptr - y[i]*cols + x[i];
      tempmagptr = edgemagptr - y[i]*cols + x[i];

      if((*tempmapptr == POSSIBLE_EDGE) && (*tempmagptr > lowval)){
         *tempmapptr = (unsigned char) EDGE;
         follow_edges(tempmapptr,tempmagptr, lowval, cols);
      }
   }
   //clock_t finish = clock();
   //fol_totaltime += (double)(finish-start)/CLOCKS_PER_SEC; 
}

/*******************************************************************************
* PROCEDURE: apply_hysteresis
* PURPOSE: This routine finds edges that are above some high threshhold or
* are connected to a high pixel by a path of pixels greater than a low
* threshold.
* NAME: Mike Heath
* DATE: 2/15/96
*******************************************************************************/

void apply_hysteresis(short int *mag, unsigned char *nms, int rows, int cols,
	float tlow, float thigh, unsigned char *edge)
{
   app_times++;
//   clock_t start = clock();
   int r, c, pos, numedges, highcount, lowthreshold, highthreshold, hist[32768];
   short int maximum_mag;

   /****************************************************************************
   * Initialize the edge map to possible edges everywhere the non-maximal
   * suppression suggested there could be an edge except for the border. At
   * the border we say there can not be an edge because it makes the
   * follow_edges algorithm more efficient to not worry about tracking an
   * edge off the side of the image.
   ****************************************************************************/
   for(r=0,pos=0;r<rows;r++){
      for(c=0;c<cols;c++,pos++){
	 if(nms[pos] == POSSIBLE_EDGE) edge[pos] = POSSIBLE_EDGE;
	 else edge[pos] = NOEDGE;
      }
   }

   for(r=0,pos=0;r<rows;r++,pos+=cols){
      edge[pos] = NOEDGE;
      edge[pos+cols-1] = NOEDGE;
   }
   pos = (rows-1) * cols;
   for(c=0;c<cols;c++,pos++){
      edge[c] = NOEDGE;
      edge[pos] = NOEDGE;
   }

   /****************************************************************************
   * Compute the histogram of the magnitude image. Then use the histogram to
   * compute hysteresis thresholds.
   ****************************************************************************/
   for(r=0;r<32768;r++) hist[r] = 0;
   for(r=0,pos=0;r<rows;r++){
      for(c=0;c<cols;c++,pos++){
	 if(edge[pos] == POSSIBLE_EDGE) hist[mag[pos]]++;
      }
   }

   /****************************************************************************
   * Compute the number of pixels that passed the nonmaximal suppression.
   ****************************************************************************/
   for(r=1,numedges=0;r<32768;r++){
      if(hist[r] != 0) maximum_mag = r;
      numedges += hist[r];
   }

   highcount = (int)(numedges * thigh + 0.5);

   /****************************************************************************
   * Compute the high threshold value as the (100 * thigh) percentage point
   * in the magnitude of the gradient histogram of all the pixels that passes
   * non-maximal suppression. Then calculate the low threshold as a fraction
   * of the computed high threshold value. John Canny said in his paper
   * "A Computational Approach to Edge Detection" that "The ratio of the
   * high to low threshold in the implementation is in the range two or three
   * to one." That means that in terms of this implementation, we should
   * choose tlow ~= 0.5 or 0.33333.
   ****************************************************************************/
   r = 1;
   numedges = hist[1];
   while((r<(maximum_mag-1)) && (numedges < highcount)){
      r++;
      numedges += hist[r];
   }
   highthreshold = r;
   lowthreshold = (int)(highthreshold * tlow + 0.5);

   if(VERBOSE){
      printf("The input low and high fractions of %f and %f computed to\n",
	 tlow, thigh);
      printf("magnitude of the gradient threshold values of: %d %d\n",
	 lowthreshold, highthreshold);
   }

   /****************************************************************************
   * This loop looks for pixels above the highthreshold to locate edges and
   * then calls follow_edges to continue the edge.
   ****************************************************************************/
   for(r=0,pos=0;r<rows;r++){
      for(c=0;c<cols;c++,pos++){
	 if((edge[pos] == POSSIBLE_EDGE) && (mag[pos] >= highthreshold)){
            edge[pos] = EDGE;
            clock_t start_fol = clock();
            follow_edges((edge+pos), (mag+pos), lowthreshold, cols);
            clock_t finish_fol = clock();
            fol_totaltime += (double)(finish_fol - start_fol)/CLOCKS_PER_SEC;

	 }
      }
   }

   /****************************************************************************
   * Set all the remaining possible edges to non-edges.
   ****************************************************************************/
   for(r=0,pos=0;r<rows;r++){
      for(c=0;c<cols;c++,pos++) if(edge[pos] != EDGE) edge[pos] = NOEDGE;
   }

  // clock_t finish = clock();
  // app_totaltime += (double)(finish-start)/CLOCKS_PER_SEC; 
}

/*******************************************************************************
* PROCEDURE: non_max_supp
* PURPOSE: This routine applies non-maximal suppression to the magnitude of
* the gradient image.
* NAME: Mike Heath
* DATE: 2/15/96
*******************************************************************************/

void non_max_supp(short *mag, short *gradx, short *grady, int nrows, int ncols,
    unsigned char *result)
{
    non_times++;
    //clock_t start = clock();

    int rowcount, colcount,count;
    short *magrowptr,*magptr;
    short *gxrowptr,*gxptr;
    short *gyrowptr,*gyptr,z1,z2;
    short m00,gx,gy;
    float mag1,mag2,xperp,yperp;
    unsigned char *resultrowptr, *resultptr;

   /****************************************************************************
   * Zero the edges of the result image.
   ****************************************************************************/
    for(count=0,resultrowptr=result,resultptr=result+ncols*(nrows-1);
        count<ncols; resultptr++,resultrowptr++,count++){
        *resultrowptr = *resultptr = (unsigned char) 0;
    }

    for(count=0,resultptr=result,resultrowptr=result+ncols-1;
        count<nrows; count++,resultptr+=ncols,resultrowptr+=ncols){
        *resultptr = *resultrowptr = (unsigned char) 0;
    }

   /****************************************************************************
   * Suppress non-maximum points.
   ****************************************************************************/
   for(rowcount=1,magrowptr=mag+ncols+1,gxrowptr=gradx+ncols+1,
      gyrowptr=grady+ncols+1,resultrowptr=result+ncols+1;
      rowcount<=nrows-2;	// bug fix 3/29/17, RD
      rowcount++,magrowptr+=ncols,gyrowptr+=ncols,gxrowptr+=ncols,
      resultrowptr+=ncols){
      for(colcount=1,magptr=magrowptr,gxptr=gxrowptr,gyptr=gyrowptr,
         resultptr=resultrowptr;colcount<=ncols-2;	// bug fix 3/29/17, RD
         colcount++,magptr++,gxptr++,gyptr++,resultptr++){
         m00 = *magptr;
         if(m00 == 0){
            *resultptr = (unsigned char) NOEDGE;
         }
         else{
            xperp = -(gx = *gxptr)/((float)m00);
            yperp = (gy = *gyptr)/((float)m00);
           // gx = *gxptr;
			//gy = *gyptr;
			//xperp = -(gx<<16)/m00;
			//yperp = (gy<<16)/m00; 
         }

         if(gx >= 0){
            if(gy >= 0){
                    if (gx >= gy)
                    {
                        /* 111 */
                        /* Left point */
                        z1 = *(magptr - 1);
                        z2 = *(magptr - ncols - 1);

                        mag1 = (m00 - z1)*xperp + (z2 - z1)*yperp;

                        /* Right point */
                        z1 = *(magptr + 1);
                        z2 = *(magptr + ncols + 1);

                        mag2 = (m00 - z1)*xperp + (z2 - z1)*yperp;
                    }
                    else
                    {
                        /* 110 */
                        /* Left point */
                        z1 = *(magptr - ncols);
                        z2 = *(magptr - ncols - 1);

                        mag1 = (z1 - z2)*xperp + (z1 - m00)*yperp;

                        /* Right point */
                        z1 = *(magptr + ncols);
                        z2 = *(magptr + ncols + 1);

                        mag2 = (z1 - z2)*xperp + (z1 - m00)*yperp;
                    }
                }
                else
                {
                    if (gx >= -gy)
                    {
                        /* 101 */
                        /* Left point */
                        z1 = *(magptr - 1);
                        z2 = *(magptr + ncols - 1);

                        mag1 = (m00 - z1)*xperp + (z1 - z2)*yperp;

                        /* Right point */
                        z1 = *(magptr + 1);
                        z2 = *(magptr - ncols + 1);

                        mag2 = (m00 - z1)*xperp + (z1 - z2)*yperp;
                    }
                    else
                    {
                        /* 100 */
                        /* Left point */
                        z1 = *(magptr + ncols);
                        z2 = *(magptr + ncols - 1);

                        mag1 = (z1 - z2)*xperp + (m00 - z1)*yperp;

                        /* Right point */
                        z1 = *(magptr - ncols);
                        z2 = *(magptr - ncols + 1);

                        mag2 = (z1 - z2)*xperp  + (m00 - z1)*yperp;
                    }
                }
            }
            else
            {
                if ((gy = *gyptr) >= 0)
                {
                    if (-gx >= gy)
                    {
                        /* 011 */
                        /* Left point */
                        z1 = *(magptr + 1);
                        z2 = *(magptr - ncols + 1);

                        mag1 = (z1 - m00)*xperp + (z2 - z1)*yperp;

                        /* Right point */
                        z1 = *(magptr - 1);
                        z2 = *(magptr + ncols - 1);

                        mag2 = (z1 - m00)*xperp + (z2 - z1)*yperp;
                    }
                    else
                    {
                        /* 010 */
                        /* Left point */
                        z1 = *(magptr - ncols);
                        z2 = *(magptr - ncols + 1);

                        mag1 = (z2 - z1)*xperp + (z1 - m00)*yperp;

                        /* Right point */
                        z1 = *(magptr + ncols);
                        z2 = *(magptr + ncols - 1);

                        mag2 = (z2 - z1)*xperp + (z1 - m00)*yperp;
                    }
                }
                else
                {
                    if (-gx > -gy)
                    {
                        /* 001 */
                        /* Left point */
                        z1 = *(magptr + 1);
                        z2 = *(magptr + ncols + 1);

                        mag1 = (z1 - m00)*xperp + (z1 - z2)*yperp;

                        /* Right point */
                        z1 = *(magptr - 1);
                        z2 = *(magptr - ncols - 1);

                        mag2 = (z1 - m00)*xperp + (z1 - z2)*yperp;
                    }
                    else
                    {
                        /* 000 */
                        /* Left point */
                        z1 = *(magptr + ncols);
                        z2 = *(magptr + ncols + 1);

                        mag1 = (z2 - z1)*xperp + (m00 - z1)*yperp;

                        /* Right point */
                        z1 = *(magptr - ncols);
                        z2 = *(magptr - ncols - 1);

                        mag2 = (z2 - z1)*xperp + (m00 - z1)*yperp;
                    }
                }
            }

            /* Now determine if the current point is a maximum point */

            if ((mag1 > 0.0) || (mag2 > 0.0))
            {
                *resultptr = (unsigned char) NOEDGE;
            }
            else
            {
                if (mag2 == 0.0)
                    *resultptr = (unsigned char) NOEDGE;
                else
                    *resultptr = (unsigned char) POSSIBLE_EDGE;
            }
        }
    }
    //clock_t finish = clock();
    //non_totaltime += (double)(finish-start)/CLOCKS_PER_SEC; 
}

/******************************************************************************
* Function: read_pgm_image
* Purpose: This function reads in an image in PGM format. The image can be
* read in from either a file or from standard input. The image is only read
* from standard input when infilename = NULL. Because the PGM format includes
* the number of columns and the number of rows in the image, these are read
* from the file. Memory to store the image is allocated OUTSIDE this function.
* The found image size is checked against the expected rows and cols.
* All comments in the header are discarded in the process of reading the
* image. Upon failure, this function returns 0, upon sucess it returns 1.
******************************************************************************/

int read_pgm_image(const char *infilename, unsigned char *image, int rows, int cols)
{
   read_times++;
   clock_t start = clock();
   FILE *fp;
   char buf[71];
   int r, c;

   /***************************************************************************
   * Open the input image file for reading if a filename was given. If no
   * filename was provided, set fp to read from standard input.
   ***************************************************************************/
   if(infilename == NULL) fp = stdin;
   else{
      if((fp = fopen(infilename, "r")) == NULL){
         fprintf(stderr, "Error reading the file %s in read_pgm_image().\n",
            infilename);
         return(0);
      }
   }

   /***************************************************************************
   * Verify that the image is in PGM format, read in the number of columns
   * and rows in the image and scan past all of the header information.
   ***************************************************************************/
   fgets(buf, 70, fp);
   if(strncmp(buf,"P5",2) != 0){
      fprintf(stderr, "The file %s is not in PGM format in ", infilename);
      fprintf(stderr, "read_pgm_image().\n");
      if(fp != stdin) fclose(fp);
      return(0);
   }
   do{ fgets(buf, 70, fp); }while(buf[0] == '#');  /* skip all comment lines */
   sscanf(buf, "%d %d", &c, &r);
   if(c != cols || r != rows){
      fprintf(stderr, "The file %s is not a %d by %d image in ", infilename, cols, rows);
      fprintf(stderr, "read_pgm_image().\n");
      if(fp != stdin) fclose(fp);
      return(0);
   }
   do{ fgets(buf, 70, fp); }while(buf[0] == '#');  /* skip all comment lines */

   /***************************************************************************
   * Read the image from the file.
   ***************************************************************************/
   if((unsigned)rows != fread(image, cols, rows, fp)){
      fprintf(stderr, "Error reading the image data in read_pgm_image().\n");
      if(fp != stdin) fclose(fp);
      return(0);
   }

   if(fp != stdin) fclose(fp);
   clock_t finish = clock();
   read_totaltime += (double)(finish-start)/CLOCKS_PER_SEC; 
   return(1);
}

/******************************************************************************
* Function: write_pgm_image
* Purpose: This function writes an image in PGM format. The file is either
* written to the file specified by outfilename or to standard output if
* outfilename = NULL. A comment can be written to the header if coment != NULL.
******************************************************************************/

int write_pgm_image(const char *outfilename, unsigned char *image, int rows,
    int cols, const char *comment, int maxval)
{
   write_times++;
   clock_t start = clock();
   FILE *fp;

   /***************************************************************************
   * Open the output image file for writing if a filename was given. If no
   * filename was provided, set fp to write to standard output.
   ***************************************************************************/
   if(outfilename == NULL) fp = stdout;
   else{
      if((fp = fopen(outfilename, "w")) == NULL){
         fprintf(stderr, "Error writing the file %s in write_pgm_image().\n",
            outfilename);
         return(0);
      }
   }

   /***************************************************************************
   * Write the header information to the PGM file.
   ***************************************************************************/
   fprintf(fp, "P5\n%d %d\n", cols, rows);
   if(comment != NULL)
      if(strlen(comment) <= 70) fprintf(fp, "# %s\n", comment);
   fprintf(fp, "%d\n", maxval);

   /***************************************************************************
   * Write the image data to the file.
   ***************************************************************************/
   if((unsigned)rows != fwrite(image, cols, rows, fp)){
      fprintf(stderr, "Error writing the image data in write_pgm_image().\n");
      if(fp != stdout) fclose(fp);
      return(0);
   }

   if(fp != stdout) fclose(fp);
   clock_t finish = clock();
   write_totaltime += (double)(finish-start)/CLOCKS_PER_SEC; 
   return(1);
}

