/* source: http://marathon.csee.usf.edu/edge/edge_detection.html */
/* URL: ftp://figment.csee.usf.edu/pub/Edge_Comparison/source_code/canny.src */

/* ECPS 203 Assignment 5 solution */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "systemc.h"

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
#define SET_STACK_SIZE set_stack_size(128*1024*1024);
/* upper bound for the size of the gaussian kernel
 * SIGMA must be less than 4.0
 * check for 'windowsize' below
 */
#define WINSIZE 21

typedef struct Float_Image
{
	float img[SIZE];

	Float_Image(void)
	{
	   for (int i=0; i<SIZE; i++)
	   { 
	      img[i] = 0;
	   }
	}

	Float_Image& operator=(const Float_Image& copy)
	{
	   for (int i=0; i<SIZE; i++)
	   { 
	      img[i] = copy.img[i];
	   }            
	   return *this;
	}

	operator float*()
	{
	   return img;
	}

	float& operator[](const int index)
	{
	   return img[index];
	}
} FIMAGE;

typedef struct Float_Kernel
{
	float img[WINSIZE];

	Float_Kernel(void)
	{
	   for (int i=0; i<WINSIZE; i++)
	   { 
	      img[i] = 0;
	   }
	}

	Float_Kernel& operator=(const Float_Kernel& copy)
	{
	   for (int i=0; i<WINSIZE; i++)
	   { 
	      img[i] = copy.img[i];
	   }            
	   return *this;
	}

	operator float*()
	{
	   return img;
	}

	float& operator[](const int index)
	{
	   return img[index];
	}
} FKERNEL;


typedef struct Image_s
{
	unsigned char img[SIZE];

	Image_s(void)
	{
	   for (int i=0; i<SIZE; i++)
	   { 
	      img[i] = 0;
	   }
	}

	Image_s& operator=(const Image_s& copy)
	{
	   for (int i=0; i<SIZE; i++)
	   { 
	      img[i] = copy.img[i];
	   }            
	   return *this;
	}

	operator unsigned char*()
	{
	   return img;
	}

	unsigned char& operator[](const int index)
	{
	   return img[index];
	}
} IMAGE;

typedef struct Image_Short
{
	short int img[SIZE];

	Image_Short(void)
	{
	   for (int i=0; i<SIZE; i++)
	   { 
	      img[i] = 0;
	   }
	}

	Image_Short& operator=(const Image_Short& copy)
	{
	   for (int i=0; i<SIZE; i++)
	   { 
	      img[i] = copy.img[i];
	   }            
	   return *this;
	}

	operator short int*()
	{
	   return img;
	}

	short int& operator[](const int index)
	{
	   return img[index];
	}
} SIMAGE;

SC_MODULE(Stimulus){
	IMAGE imageout;
	sc_fifo_out<IMAGE> ImgOut;

	int read_pgm_image(const char *infilename, unsigned char *image, int rows, int cols)
	{
	   FILE *fp;
	   char buf[71];
	   int r, c;

	   if(infilename == NULL) fp = stdin;
	   else{
	      if((fp = fopen(infilename, "r")) == NULL){
	         fprintf(stderr, "Error reading the file %s in read_pgm_image().\n",
	            infilename);
	         return(0);
	      }
	   }

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

	   if((unsigned)rows != fread(image, cols, rows, fp)){
	      fprintf(stderr, "Error reading the image data in read_pgm_image().\n");
	      if(fp != stdin) fclose(fp);
	      return(0);
	   }

	   if(fp != stdin) fclose(fp);
	   return(1);
	}

	void main(void)
	{
	   int i=0, n=0;
	   char infilename[40];

	   for(i=0; i<IMG_NUM; i++)
	   {
	      n = i % AVAIL_IMG;
	      sprintf(infilename, IMG_IN, n+1);

	      if(VERBOSE) printf("Reading the image %s.\n", infilename);
	      if(read_pgm_image(infilename, imageout, ROWS, COLS) == 0){
	         fprintf(stderr, "Error reading the input image, %s.\n", infilename);
	         exit(1);
	      }
	      ImgOut.write(imageout);
	   }
	}

	SC_CTOR(Stimulus)
	{
	   SC_THREAD(main);
	   SET_STACK_SIZE
	}
};

SC_MODULE(Monitor){
	IMAGE imagein;
	sc_fifo_in<IMAGE>  ImgIn;

	int write_pgm_image(const char *outfilename, unsigned char *image, int rows,
	    int cols, const char *comment, int maxval)
	{
	   FILE *fp;

	   if(outfilename == NULL) fp = stdout;
	   else{
	      if((fp = fopen(outfilename, "w")) == NULL){
	         fprintf(stderr, "Error writing the file %s in write_pgm_image().\n",
	            outfilename);
	         return(0);
	      }
	   }

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
	   return(1);
	}

	void main(void)
	{
	   char outfilename[128];    /* Name of the output "edge" image */
	   int i, n;

	   for(i=0; i<IMG_NUM; i++)
	   {
	      ImgIn.read(imagein);

	      /****************************************************************************
	      * Write out the edge image to a file.
	      ****************************************************************************/
	      n = i % AVAIL_IMG;
	      sprintf(outfilename, IMG_OUT, n+1);
	      if(VERBOSE) printf("Writing the edge image in the file %s.\n", outfilename);
	      if(write_pgm_image(outfilename, imagein, ROWS, COLS,"", 255) == 0){
	         fprintf(stderr, "Error writing the edge image, %s.\n", outfilename);
	         exit(1);
	      }
	   }
	   if(VERBOSE) printf("Monitor exits simulation.\n");
	      sc_stop();	// done testing, quit the simulation
	}
	
	SC_CTOR(Monitor)
	{
	   SC_THREAD(main);
	   SET_STACK_SIZE
	}
};

SC_MODULE(DataIn)
{
	sc_fifo_in<IMAGE> ImgIn;
	sc_fifo_out<IMAGE> ImgOut;
	IMAGE Image;

	void main()
	{
	   while(1)
	   {
	      ImgIn.read(Image);
	      ImgOut.write(Image);
	   }
	}
	
	SC_CTOR(DataIn)
	{
	   SC_THREAD(main);
	   SET_STACK_SIZE
	}
};

SC_MODULE(DataOut)
{
	sc_fifo_in<IMAGE> ImgIn;
	sc_fifo_out<IMAGE> ImgOut;
	IMAGE Image;

	void main()
	{
	   while(1)
	   {
	      ImgIn.read(Image);
	      ImgOut.write(Image);
	   }
	}
	
	SC_CTOR(DataOut)
	{
	   SC_THREAD(main);
	   SET_STACK_SIZE
	}
};

SC_MODULE(Receive_Image){
	sc_fifo_in<IMAGE> Receive_In;
	sc_fifo_out<IMAGE> Receive_Out;
	IMAGE img;
	void main(){
		while(true){
			Receive_In.read(img);
			Receive_Out.write(img);
		}
	}
	SC_CTOR(Receive_Image){
		SC_THREAD(main);
		SET_STACK_SIZE
	}
};

SC_MODULE(BlurX){
	//image from gaus
	IMAGE imgin;
	//image to blury
	FIMAGE tempim;
	//kernel to blur x
	FKERNEL kernel;
	//window size
	int windowsize;
	sc_fifo_in<FKERNEL> BX_In_Kernel;
	sc_fifo_in<IMAGE> BX_In_Img;
	sc_fifo_in<int> BX_In_WINDOW;
	sc_fifo_out<FIMAGE> BX_OUT_BY;

	void blurx(unsigned char *imgin, float *kernel, int rows, int cols, float *tempim, int *windowsize){
		int r, c, rr, cc, center;            /* Half of the windowsize. */
	    float dot, sum;            /* Sum of the kernel weights variable. */
	    if(VERBOSE) printf("   Bluring the image in the X-direction.\n");
	    center = *windowsize / 2;
	    //blurx
		for(r=0;r<rows;r++){
		   for(c=0;c<cols;c++){
		   	  dot = 0.0;
		      sum = 0.0;
		      for(cc=(-center);cc<=center;cc++){
		         if(((c+cc) >= 0) && ((c+cc) < cols)){
		            dot += (float)imgin[r*cols+(c+cc)] * kernel[center+cc];
		            sum += kernel[center+cc];
		         }
		      }
	          tempim[r*cols+c] = dot/sum;
		   }
		}
	}

	void main(){
		while(true){
			BX_In_Img.read(imgin);
			BX_In_Kernel.read(kernel);
			BX_In_WINDOW.read(windowsize);
			blurx(imgin, kernel, ROWS, COLS, tempim, &windowsize);
			BX_OUT_BY.write(tempim);
		}
	}

	SC_CTOR(BlurX){
	   SC_THREAD(main);
	   SET_STACK_SIZE
	}
};

SC_MODULE(BlurY){
	//image from blurx
	FIMAGE tempim;
	//output image
	SIMAGE smoothedim;
	//from kernel
	FKERNEL kernel;
	//windowsize
	int windowsize;
	sc_fifo_in<FKERNEL> BY_In_Kernel;
	sc_fifo_in<FIMAGE> BY_In_Tempim;
	sc_fifo_in<int> BY_In_WINDOW;
	sc_fifo_out<SIMAGE> BY_OUT_SmoothedIm;

	void blury(float *tempim, float *kernel, int rows, int cols, short int *smoothedim, int *windowsize){
	   int r, c, rr, cc, center;     /* Counter variables. */
	   float dot, sum;            /* Sum of the kernel weights variable. */
	   center = *windowsize / 2;
	   //blury
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
	}
	void main(){
		while(true){
			BY_In_Tempim.read(tempim);
			BY_In_WINDOW.read(windowsize);
			BY_In_Kernel.read(kernel);
			blury(tempim, kernel, ROWS, COLS, smoothedim, &windowsize);
			BY_OUT_SmoothedIm.write(smoothedim);
		}
	}
	SC_CTOR(BlurY){
	   SC_THREAD(main);
	   SET_STACK_SIZE
	}
};


SC_MODULE(Gaussian_Kernel){
	FKERNEL kernel;
	int windowsize;
	sc_fifo_out<FKERNEL> Kernel_Out_Kernel;
	sc_fifo_out<FKERNEL> Kernel_Out_Kernel2;

	sc_fifo_out<int> Kernel_Out_Windowsize;
	sc_fifo_out<int> Kernel_Out_Windowsize2;
	void make_gaussian_kernel(float sigma, float *kernel, int *windowsize);
	void main(){
		while(true){
			make_gaussian_kernel(SIGMA, kernel, &windowsize);
			Kernel_Out_Kernel.write(kernel);
			Kernel_Out_Kernel2.write(kernel);

			Kernel_Out_Windowsize.write(windowsize);
			Kernel_Out_Windowsize2.write(windowsize);
		}
	}

	SC_CTOR(Gaussian_Kernel){
		SC_THREAD(main);
		SET_STACK_SIZE
	}
};

void Gaussian_Kernel::make_gaussian_kernel(float sigma, float *kernel, int *windowsize)
	{
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
	}

SC_MODULE(Gaussian_Smooth){
	//modules
	Receive_Image receive_image;
	Gaussian_Kernel gaussian_kernel;
	BlurX blurx;
	BlurY blury;
	//channels
	sc_fifo<IMAGE> REC_BX;
	sc_fifo<FKERNEL> KER_BX;
	sc_fifo<FKERNEL> KER_BY;
	sc_fifo<int> WIN_BX;
	sc_fifo<int> WIN_BY;
	sc_fifo<FIMAGE> BX_BY;
	//ports
	sc_fifo_in<IMAGE> Gaussian_In;
	sc_fifo_out<SIMAGE> Gaussian_Out;
	//data
	IMAGE img;

	void before_end_of_elaboration(){
		//img in
		receive_image.Receive_In.bind(Gaussian_In);
		//receive and blurx
		receive_image.Receive_Out.bind(REC_BX);
		blurx.BX_In_Img.bind(REC_BX);
		//kernel and blurx
		gaussian_kernel.Kernel_Out_Kernel.bind(KER_BX);
		blurx.BX_In_Kernel.bind(KER_BX);
		//kernel and blury
		gaussian_kernel.Kernel_Out_Kernel2.bind(KER_BY);
		blury.BY_In_Kernel.bind(KER_BY);
		//windowsize and blurx
		gaussian_kernel.Kernel_Out_Windowsize.bind(WIN_BX);
		blurx.BX_In_WINDOW.bind(WIN_BX);
		//windowsize and blury
		gaussian_kernel.Kernel_Out_Windowsize2.bind(WIN_BY);
		blury.BY_In_WINDOW.bind(WIN_BY);
		//blurx and blury
		blurx.BX_OUT_BY.bind(BX_BY);
		blury.BY_In_Tempim.bind(BX_BY);
		//blury out
		blury.BY_OUT_SmoothedIm.bind(Gaussian_Out);
	}

	SC_CTOR(Gaussian_Smooth)
	:REC_BX("REC_BX", 1)
	,KER_BX("KER_BX", 1)
	,KER_BY("KER_BY", 1)
	,WIN_BX("WIN_BX", 1)
	,WIN_BY("WIN_BY", 1)
	,BX_BY("BX_BY", 1)
	,receive_image("receive_image")
	,gaussian_kernel("gaussian_kernel")
	,blurx("blurx")
	,blury("blury")
	{}
};

SC_MODULE(Derivative_X_Y){
	SIMAGE smoothedim;
	SIMAGE delta_x;
	SIMAGE delta_y;
	sc_fifo_in<SIMAGE> Derivative_In;
	sc_fifo_out<SIMAGE> Derivative_Out_X;
	sc_fifo_out<SIMAGE> Derivative_Out_X2;
	sc_fifo_out<SIMAGE> Derivative_Out_Y;
	sc_fifo_out<SIMAGE> Derivative_Out_Y2;

	//derivate_x_y
	void derivative_x_y(short int *smoothedim, int rows, int cols,
	        short int *delta_x, short int *delta_y)
	{
	   int r, c, pos;
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
	}
	
	void main(){
		while(true){
			Derivative_In.read(smoothedim);
			derivative_x_y(smoothedim, ROWS, COLS, delta_x, delta_y);
			Derivative_Out_X.write(delta_x);
			Derivative_Out_X2.write(delta_x);
			Derivative_Out_Y.write(delta_y);
			Derivative_Out_Y2.write(delta_y);

		}
	}
	
	SC_CTOR(Derivative_X_Y)
	{
	   SC_THREAD(main);
	   SET_STACK_SIZE
	}
};

SC_MODULE(Magnitude_X_Y){
	SIMAGE delta_x;
	SIMAGE delta_y;
	SIMAGE magnitude;
	sc_fifo_in<SIMAGE> Magnitude_In_X;
	sc_fifo_in<SIMAGE> Magnitude_In_Y;

	sc_fifo_out<SIMAGE> Magnitude_Out;
	sc_fifo_out<SIMAGE> Magnitude_Out2;


	void magnitude_x_y(short int *delta_x, short int *delta_y, int rows, int cols,
	        short int *magnitude)
	{
	   int r, c, pos, sq1, sq2;

	   for(r=0,pos=0;r<rows;r++){
	      for(c=0;c<cols;c++,pos++){
	         sq1 = (int)delta_x[pos] * (int)delta_x[pos];
	         sq2 = (int)delta_y[pos] * (int)delta_y[pos];
	         magnitude[pos] = (short)(0.5 + sqrt((float)sq1 + (float)sq2));
	      }
	   }
	}
	void main(){
		while(true){
			Magnitude_In_X.read(delta_x);
			Magnitude_In_Y.read(delta_y);
			magnitude_x_y(delta_x, delta_y, ROWS, COLS, magnitude);
			Magnitude_Out.write(magnitude);
			Magnitude_Out2.write(magnitude);

		}
	}

	SC_CTOR(Magnitude_X_Y)
	{
	   SC_THREAD(main);
	   SET_STACK_SIZE
	}
};

SC_MODULE(Non_Max_Supp){
	SIMAGE mag;
	SIMAGE gradx;
	SIMAGE grady;
	IMAGE result;
	sc_fifo_in<SIMAGE> Max_In_D;
	sc_fifo_in<SIMAGE> Max_In_X;
	sc_fifo_in<SIMAGE> Max_In_Y;
	sc_fifo_out<IMAGE> Max_Out;

	void non_max_supp(short *mag, short *gradx, short *grady, int nrows, int ncols,
	    unsigned char *result)
	{
	    int rowcount, colcount,count;
	    short *magrowptr,*magptr;
	    short *gxrowptr,*gxptr;
	    short *gyrowptr,*gyptr,z1,z2;
	    short m00,gx,gy;
	    float mag1,mag2,xperp,yperp;
	    unsigned char *resultrowptr, *resultptr;

	    for(count=0,resultrowptr=result,resultptr=result+ncols*(nrows-1);
	        count<ncols; resultptr++,resultrowptr++,count++){
	        *resultrowptr = *resultptr = (unsigned char) 0;
	    }

	    for(count=0,resultptr=result,resultrowptr=result+ncols-1;
	        count<nrows; count++,resultptr+=ncols,resultrowptr+=ncols){
	        *resultptr = *resultrowptr = (unsigned char) 0;
	    }

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
	}

	void main(){
		while(true){
			Max_In_D.read(mag);
			Max_In_X.read(gradx);
			Max_In_Y.read(grady);
			non_max_supp(mag, gradx, grady, ROWS, COLS, result);
			Max_Out.write(result);
		}
	}

	SC_CTOR(Non_Max_Supp)
	{
	   SC_THREAD(main);
	   SET_STACK_SIZE
	}
};

SC_MODULE(Apply_Hysteresis){
	SIMAGE mag;
	IMAGE nms;
	IMAGE edge;

	sc_fifo_in<SIMAGE> Apply_In_Mag;
	sc_fifo_in<IMAGE> Apply_In_Non;
	sc_fifo_out<IMAGE> Apply_Out;

	void follow_edges(unsigned char *edgemapptr, short *edgemagptr, short lowval,
	   int cols)
	{
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
	}

	void apply_hysteresis(short int *mag, unsigned char *nms, int rows, int cols,
		float tlow, float thigh, unsigned char *edge)
	{
	   int r, c, pos, numedges, highcount, lowthreshold, highthreshold, hist[32768];
	   short int maximum_mag;

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

	   for(r=0;r<32768;r++) hist[r] = 0;
	   for(r=0,pos=0;r<rows;r++){
	      for(c=0;c<cols;c++,pos++){
		 if(edge[pos] == POSSIBLE_EDGE) hist[mag[pos]]++;
	      }
	   }

	   for(r=1,numedges=0;r<32768;r++){
	      if(hist[r] != 0) maximum_mag = r;
	      numedges += hist[r];
	   }

	   highcount = (int)(numedges * thigh + 0.5);

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

	   for(r=0,pos=0;r<rows;r++){
	      for(c=0;c<cols;c++,pos++){
		 if((edge[pos] == POSSIBLE_EDGE) && (mag[pos] >= highthreshold)){
	            edge[pos] = EDGE;
	            follow_edges((edge+pos), (mag+pos), lowthreshold, cols);
		 }
	      }
	   }

	   for(r=0,pos=0;r<rows;r++){
	      for(c=0;c<cols;c++,pos++) if(edge[pos] != EDGE) edge[pos] = NOEDGE;
	   }
	}

	void main(){
		while(true){
			Apply_In_Mag.read(mag);
			Apply_In_Non.read(nms);
			apply_hysteresis(mag, nms, ROWS, COLS, TLOW, THIGH, edge);
			Apply_Out.write(edge);
		}
	}

	SC_CTOR(Apply_Hysteresis)
	{
	   SC_THREAD(main);
	   SET_STACK_SIZE
	}
};


SC_MODULE(DUT)
{
	//modules
	Gaussian_Smooth gaussian_smooth;
	Derivative_X_Y derivative_x_y;
	Magnitude_X_Y magnitude_x_y;
	Non_Max_Supp non_max_supp;
	Apply_Hysteresis apply_hysteresis;
	//channels
	sc_fifo<SIMAGE> GAUS_DERI;
	sc_fifo<SIMAGE> DERI_NON_MAX_X;
	sc_fifo<SIMAGE> DERI_NON_MAX_Y;
	sc_fifo<SIMAGE> DERI_MAG_X;
	sc_fifo<SIMAGE> DERI_MAG_Y;
	sc_fifo<SIMAGE> MAG_APP;
	sc_fifo<SIMAGE> MAG_NON_MAX;
	sc_fifo<IMAGE> NON_MAX_APP;
	//ports
	sc_fifo_in<IMAGE> ImgIn;
	sc_fifo_out<IMAGE> ImgOut;
	//data
	IMAGE imageout;
	IMAGE imagein;

	void before_end_of_elaboration(){
		//image to gaussian smooth
		gaussian_smooth.Gaussian_In.bind(ImgIn);
		//gaussian smooth and derivatexy
		gaussian_smooth.Gaussian_Out.bind(GAUS_DERI);
	    derivative_x_y.Derivative_In.bind(GAUS_DERI);
	    //derivativexy and nonmaxsupp
	    derivative_x_y.Derivative_Out_X.bind(DERI_NON_MAX_X);
	    derivative_x_y.Derivative_Out_Y.bind(DERI_NON_MAX_Y);
	    non_max_supp.Max_In_X.bind(DERI_NON_MAX_X);
	    non_max_supp.Max_In_Y.bind(DERI_NON_MAX_Y);
	    //derivativexy and magnitudexy
	    derivative_x_y.Derivative_Out_X2.bind(DERI_MAG_X);
	    derivative_x_y.Derivative_Out_Y2.bind(DERI_MAG_Y);
	    magnitude_x_y.Magnitude_In_X.bind(DERI_MAG_X);
	    magnitude_x_y.Magnitude_In_Y.bind(DERI_MAG_Y);
	    //magnitude and nonmax
	    magnitude_x_y.Magnitude_Out.bind(MAG_NON_MAX);
	    non_max_supp.Max_In_D.bind(MAG_NON_MAX);
	    //magnitude and apply
	    magnitude_x_y.Magnitude_Out2.bind(MAG_APP);
	    apply_hysteresis.Apply_In_Mag.bind(MAG_APP);
	    //nonmax and apply
	    non_max_supp.Max_Out.bind(NON_MAX_APP);
	    apply_hysteresis.Apply_In_Non.bind(NON_MAX_APP);
	    //image out
	    apply_hysteresis.Apply_Out.bind(ImgOut);
	}

	SC_CTOR(DUT)
	:GAUS_DERI("GAUS_DERI", 1)
	,DERI_NON_MAX_X("DERI_NON_MAX_X", 1)
	,DERI_NON_MAX_Y("DERI_NON_MAX_Y", 1)
	,DERI_MAG_X("DERI_MAG_X", 1)
	,DERI_MAG_Y("DERI_MAG_Y", 1)
	,MAG_APP("MAG_APP", 1)
	,MAG_NON_MAX("MAG_NON_MAX", 1)
	,NON_MAX_APP("NON_MAX_APP", 1)
	,gaussian_smooth("gaussian_smooth")
	,derivative_x_y("derivative_x_y")
	,magnitude_x_y("magnitude_x_y")
	,non_max_supp("non_max_supp")
	,apply_hysteresis("apply_hysteresis")
	{}
};

SC_MODULE(Platform)
{
	sc_fifo_in<IMAGE> ImgIn;
	sc_fifo_out<IMAGE> ImgOut;
	sc_fifo<IMAGE> q1;
	sc_fifo<IMAGE> q2;

	DataIn din;
	DUT canny;
	DataOut dout;

	void before_end_of_elaboration(){
	   din.ImgIn.bind(ImgIn);
	   din.ImgOut.bind(q1);
	   canny.ImgIn.bind(q1);
	   canny.ImgOut.bind(q2);
	   dout.ImgIn.bind(q2);
	   dout.ImgOut.bind(ImgOut);
	}

	SC_CTOR(Platform)
	:q1("q1",1)
	,q2("q2",1)
	,din("din")
	,canny("canny")
	,dout("dout")
	{}
};

SC_MODULE(Top)
{
	sc_fifo<IMAGE> q1;
	sc_fifo<IMAGE> q2;
	Stimulus stimulus;
	Platform platform;
	Monitor monitor;

	void before_end_of_elaboration(){
	   stimulus.ImgOut.bind(q1);
	   platform.ImgIn.bind(q1);
	   platform.ImgOut.bind(q2);
	   monitor.ImgIn.bind(q2);
	}

	SC_CTOR(Top)
	:q1("q1",1)
	,q2("q2",1)
	,stimulus("stimulus")
	,platform("platform")
	,monitor("monitor")
	{}
};

Top top("top");

int sc_main(int argc, char* argv[])
{
	sc_start();
	return 0;
}
