#include <highgui.h>
#include <cv.h>

using namespace cv;
void my_HomoFilter(Mat srcImg, Mat &dst);  
void my_Doublepla(Mat src,Mat &dst);

 double* CreateKernel(double sigma);
 int* CreateFastKernel(double sigma);
 void FilterGaussian(IplImage* img, double sigma);
 void FastFilter(IplImage *img, double sigma);

 void my_Retinex
(IplImage *img, double sigma, int gain = 128, int offset = 128);
 void MultiScaleRetinex
(IplImage *img, int scales, double *weights, double *sigmas, int gain = 128, int offset = 128);

//extern void MultiScaleRetinexCR
//(IplImage *img, int scales, double *weights, double *sigmas, int gain = 128, int offset = 128,
// double restoration_factor = 6, double color_gain = 2);