#include "stdafx.h"
#include <IOSTREAM> // write for test you can neglect it
#include "cv.h"
#include "math.h"
#include "cxcore.h"
#include "highgui.h"
#include <vector>
//#include "homo.h"
 
using namespace cv;

//以下为同态滤波代码
void my_HomoFilter(Mat srcImg, Mat &dst)  
{  
      
    srcImg.convertTo(srcImg, CV_64FC1);  
    dst.convertTo(dst, CV_64FC1);  
    //1. ln  
    for (int i = 0; i < srcImg.rows; i++)  
    {  
        double* srcdata = srcImg.ptr<double>(i);  
        double* logdata = srcImg.ptr<double>(i);  
        for (int j = 0; j < srcImg.cols; j++)  
        {  
            logdata[j] = log(srcdata[j]+0.0001);  
        }  
    }  
      
    //spectrum  
    //2. dct  
    Mat mat_dct = Mat::zeros(srcImg.size(), CV_64FC1); 
	int qwe=srcImg.channels();
	int qw=mat_dct.channels();
	//mat_dct.ere();hh
    dct(srcImg, mat_dct);  
    //imshow("dct", mat_dct);  
      
    //3. linear filter  
    Mat H_u_v;  
    double gammaH = 1;//1.5;  
    double gammaL = 0.8;//0.5;  
    double C = 1;  
    double d0 = (srcImg.rows/2)*(srcImg.rows/2) + (srcImg.cols/2)*(srcImg.cols/2);  
    double d2 = 0;  
    H_u_v = Mat::zeros(srcImg.rows, srcImg.cols, CV_64FC1);  
      
    double totalWeight = 0.0;  
    for (int i = 0; i < srcImg.rows; i++)  
    {  
        double * dataH_u_v = H_u_v.ptr<double>(i);  
        for (int j = 0; j < srcImg.cols; j++)  
        {  
            d2 = pow((i), 2.0) + pow((j), 2.0);  
            dataH_u_v[j] =  (gammaH - gammaL)*(1 - exp(-C*d2/d0)) + gammaL;  
            totalWeight += dataH_u_v[j];  
        }  
    }  
    H_u_v.ptr<double>(0)[0] = 1.1;  
      
    //H_u_v = Mat::ones(srcImg.rows, srcImg.cols, CV_64FC1);  
    //imshow("H_u_v", H_u_v);  
  
  
    //imshow("before filter", mat_dct);  
  
    mat_dct = mat_dct.mul(H_u_v);  
    //Mat tmp = mat_dct.mul(H_u_v);  
    //tmp.copyTo(mat_dct);  
    //4. idct  
    idct(mat_dct, dst);  
      
#if 0  
    //spatial high high pass filter  
    Mat tmp = Mat::zeros(srcImg.rows, srcImg.cols, CV_64FC1);  
    GaussianBlur(srcImg, tmp, Size(9, 9), 1.5, 1.5);  
    const double alpha = 0.5;  
      
    for (int i = 0; i < srcImg.rows; i++)  
    {  
        double* srcdata = srcImg.ptr<double>(i);  
        double* blurdata = tmp.ptr<double>(i);  
        double* dstdata = dst.ptr<double>(i);  
        for (int j = 0; j < srcImg.cols; j++)  
        {  
            dstdata[j] = (1+alpha)*srcdata[j] - alpha*blurdata[j];  
            //dstdata[j] = blurdata[j];  
              
        }  
    }  
       
      
      
#endif  
    //5. exp  
    for (int i = 0; i < srcImg.rows; i++)  
    {  
        double* srcdata = dst.ptr<double>(i);  
        double* dstdata = dst.ptr<double>(i);  
        for (int j = 0; j < srcImg.cols; j++)  
        {  
            dstdata[j] = exp(srcdata[j]);  
        }  
    }  
      
    //imshow("dst", dst);  
    dst.convertTo(dst, CV_8UC1);  
  
}  

//以下为双平台直方图代码
void my_Doublepla(Mat src,Mat &dst)
{
	//Mat src;
//src=Mat(img);
//Mat dest;  
//Src.copyTo(src);  
  
if (src.channels() > 1)  
{  
    cvtColor(src, src, CV_BGR2GRAY);  
}  
  
  
  
MatND  hist;  
const int histSize = 256;  
float range[] = { 0, 255 };  
const float *ranges[] = { range };  
const int channels = 0;  
cv::calcHist(&src, 1, &channels, Mat(), hist, 1, &histSize, ranges);  
float total = src.size().width* src.size().height;  
  
float bins[histSize] = { 0 };  
float binsAcc[histSize] = { 0 };  
Mat lut(1, 256, CV_8U);  
vector<float> vectorBins;  
vector<float> maxBins;  
float sumBins = 0.0;  
int countMax = 0;  
float maxValue=0;
float TValueMax = 0;  
float TValueMin = 0;  
//float TValue = 0;  
  
// Find the mapping table  
for (int i = 0; i<histSize; i++)  
{  
    float bin_val = hist.at<float>(i); // 第i灰度级上的数  
    bins[i] = bin_val / total;  
	/*int a;
		a=cv::max(bins[]);*/
    if (bins[i] > 0)  
    {  
        vectorBins.push_back(bins[i]);  
    }  
	for (int i = 0; i<histSize; i++)  {
		if (bins[i]>maxValue)
			maxValue=bins[i];
	}
      
  
}  
  

//TValue = sumBins / countMax;  
TValueMax=25*maxValue/100;
TValueMin=0.05*maxValue/100;
  
  
// Find the mapping table  
for (int i = 0; i<histSize; i++)  
{  
  
    if (bins[i] > TValueMax)  
    {  
        bins[i] = TValueMax;  
    }  
	if (bins[i] < TValueMin)  
    {  
        bins[i] = TValueMin;  
    }  
  
    if (i>0)  
    {  
        binsAcc[i] = binsAcc[i - 1] + bins[i];  
    }  
    else  
    {  
        binsAcc[0] = bins[0];  
    }  
}  
  
for (int i = 0; i < histSize; i++)  
{  
    lut.at<uchar>(i) = static_cast<uchar>(cvRound(binsAcc[i] * 255 / binsAcc[255]));  
}  
  
LUT(src, lut, dst);  
//imshow("src", src);  
//imshow("equlization", dst);  
}


//以下为RETINEX代码
//#define USE_EXACT_SIGMA


#define pc(image, x, y, c) image->imageData[(image->widthStep * y) + (image->nChannels * x) + c]

#define INT_PREC 1024.0
#define INT_PREC_BITS 10

inline double int2double(int x) { return (double)x / INT_PREC; }
inline int double2int(double x) { return (int)(x * INT_PREC + 0.5); }

inline int int2smallint(int x) { return (x >> INT_PREC_BITS); }
inline int int2bigint(int x) { return (x << INT_PREC_BITS); }


double*
CreateKernel(double sigma)
{
	int i, x, filter_size;
	double* filter;
	double sum;

	// Reject unreasonable demands
	if ( sigma > 200 ) sigma = 200;

	// get needed filter size (enforce oddness)
	filter_size = (int)floor(sigma*6) / 2;
	filter_size = filter_size * 2 + 1;

	// Allocate kernel space
	filter = new double[filter_size];

	// Calculate exponential
	sum = 0;
	for (i = 0; i < filter_size; i++) {
		x = i - (filter_size / 2);
		filter[i] = exp( -(x*x) / (2*sigma*sigma) );

		sum += filter[i];
	}

	// Normalize
	for (i = 0, x; i < filter_size; i++)
		filter[i] /= sum;

	return filter;
}


int*
CreateFastKernel(double sigma)
{
	double* fp_kernel;
	int* kernel;
	int i, filter_size;

	// Reject unreasonable demands
	if ( sigma > 200 ) sigma = 200;

	// get needed filter size (enforce oddness)
	filter_size = (int)floor(sigma*6) / 2;
	filter_size = filter_size * 2 + 1;

	// Allocate kernel space
	kernel = new int[filter_size];

	fp_kernel = CreateKernel(sigma);

	for (i = 0; i < filter_size; i++)
		kernel[i] = double2int(fp_kernel[i]);

	delete fp_kernel;

	return kernel;
}



void
FilterGaussian(IplImage* img, double sigma)
{
	int i, j, k, source, filter_size;
	int* kernel;
	IplImage* temp;
	int v1, v2, v3;

	// Reject unreasonable demands
	if ( sigma > 200 ) sigma = 200;

	// get needed filter size (enforce oddness)
	filter_size = (int)floor(sigma*6) / 2;
	filter_size = filter_size * 2 + 1;

	kernel = CreateFastKernel(sigma);

	temp = cvCreateImage(cvSize(img->width, img->height), img->depth, img->nChannels);

	// filter x axis
	for (j = 0; j < temp->height; j++)
		for (i = 0; i < temp->width; i++) {

			// inner loop has been unrolled

			v1 = v2 = v3 = 0;
			for (k = 0; k < filter_size; k++) {

				source = i + filter_size / 2 - k;

				if (source < 0) source *= -1;
				if (source > img->width - 1) source = 2*(img->width - 1) - source;

				v1 += kernel[k] * (unsigned char)pc(img, source, j, 0);
				if (img->nChannels == 1) continue;
				v2 += kernel[k] * (unsigned char)pc(img, source, j, 1);
				v3 += kernel[k] * (unsigned char)pc(img, source, j, 2);

			}

			// set value and move on
			pc(temp, i, j, 0) = (char)int2smallint(v1);
			if (img->nChannels == 1) continue;
			pc(temp, i, j, 1) = (char)int2smallint(v2);
			pc(temp, i, j, 2) = (char)int2smallint(v3);

		}

		// filter y axis
		for (j = 0; j < img->height; j++)
			for (i = 0; i < img->width; i++) {

				v1 = v2 = v3 = 0;
				for (k = 0; k < filter_size; k++) {

					source = j + filter_size / 2 - k;

					if (source < 0) source *= -1;
					if (source > temp->height - 1) source = 2*(temp->height - 1) - source;

					v1 += kernel[k] * (unsigned char)pc(temp, i, source, 0);
					if (img->nChannels == 1) continue;
					v2 += kernel[k] * (unsigned char)pc(temp, i, source, 1);
					v3 += kernel[k] * (unsigned char)pc(temp, i, source, 2);

				}

				// set value and move on
				pc(img, i, j, 0) = (char)int2smallint(v1);
				if (img->nChannels == 1) continue;
				pc(img, i, j, 1) = (char)int2smallint(v2);
				pc(img, i, j, 2) = (char)int2smallint(v3);

			}


			cvReleaseImage( &temp );

			delete kernel;

}


void
FastFilter(IplImage *img, double sigma)
{
	int filter_size;

	// Reject unreasonable demands
	if ( sigma > 200 ) sigma = 200;

	// get needed filter size (enforce oddness)
	filter_size = (int)floor(sigma*6) / 2;
	filter_size = filter_size * 2 + 1;

	// If 3 sigma is less than a pixel, why bother (ie sigma < 2/3)
	if(filter_size < 3) return;

	// Filter, or downsample and recurse
	if (filter_size < 10) {

#ifdef USE_EXACT_SIGMA
		FilterGaussian(img, sigma)
#else
		cvSmooth( img, img, CV_GAUSSIAN, filter_size, filter_size );
#endif

	}
	else {
		if (img->width < 2 || img->height < 2) return;

		IplImage* sub_img = cvCreateImage(cvSize(img->width / 2, img->height / 2), img->depth, img->nChannels);

		cvPyrDown( img, sub_img );

		FastFilter( sub_img, sigma / 2.0 );

		cvResize( sub_img, img, CV_INTER_LINEAR );

		cvReleaseImage( &sub_img );
	}

}


void
my_Retinex(IplImage *img, double sigma, int gain, int offset)
{
	IplImage *A, *fA, *fB, *fC;

	// Initialize temp images
	fA = cvCreateImage(cvSize(img->width, img->height), IPL_DEPTH_32F, img->nChannels);
	fB = cvCreateImage(cvSize(img->width, img->height), IPL_DEPTH_32F, img->nChannels);
	fC = cvCreateImage(cvSize(img->width, img->height), IPL_DEPTH_32F, img->nChannels);//原为IPL_DEPTH_32F

	// Compute log image
	cvConvert( img, fA );
	cvLog( fA, fB );

	// Compute log of blured image
	A = cvCloneImage( img );
	FastFilter( A, sigma );
	cvConvert( A, fA );
	cvLog( fA, fC );

	// Compute difference
	cvSub( fB, fC, img );
	// Restore
	/*cvConvertScale( fA, img, gain, offset);*/

	// Release temp images
	cvReleaseImage( &A );
	cvReleaseImage( &fA );
	cvReleaseImage( &fB );
	cvReleaseImage( &fC );

}


void
MultiScaleRetinex(IplImage *img, int scales, double *weights, double *sigmas, int gain, int offset)
{
	int i;
	double weight;
	IplImage *A, *fA, *fB, *fC;

	// Initialize temp images
	fA = cvCreateImage(cvSize(img->width, img->height), IPL_DEPTH_32F, img->nChannels);
	fB = cvCreateImage(cvSize(img->width, img->height), IPL_DEPTH_32F, img->nChannels);
	fC = cvCreateImage(cvSize(img->width, img->height), IPL_DEPTH_32F, img->nChannels);


	// Compute log image
	cvConvert( img, fA );
	cvLog( fA, fB );

	// Normalize according to given weights
	for (i = 0, weight = 0; i < scales; i++)
		weight += weights[i];

	if (weight != 1.0) cvScale( fB, fB, weight );

	// Filter at each scale
	for (i = 0; i < scales; i++) {
		A = cvCloneImage( img );
		FastFilter( A, sigmas[i] );

		cvConvert( A, fA );
		cvLog( fA, fC );
		cvReleaseImage( &A );

		// Compute weighted difference
		cvScale( fC, fC, weights[i] );
		cvSub( fB, fC, fB );
	}

	// Restore
	cvConvertScale( fB, img, gain, offset);

	// Release temp images
	cvReleaseImage( &fA );
	cvReleaseImage( &fB );
	cvReleaseImage( &fC );
}


//void
//MultiScaleRetinexCR(IplImage *img, int scales, double *weights, double *sigmas,
//					int gain, int offset, double restoration_factor, double color_gain)
//{
//	int i;
//	double weight;
//	IplImage *A, *B, *C, *fA, *fB, *fC, *fsA, *fsB, *fsC, *fsD, *fsE, *fsF;
//
//	// Initialize temp images
//	fA = cvCreateImage(cvSize(img->width, img->height), IPL_DEPTH_32F, img->nChannels);
//	fB = cvCreateImage(cvSize(img->width, img->height), IPL_DEPTH_32F, img->nChannels);
//	fC = cvCreateImage(cvSize(img->width, img->height), IPL_DEPTH_32F, img->nChannels);
//	fsA = cvCreateImage(cvSize(img->width, img->height), IPL_DEPTH_32F, 1);
//	fsB = cvCreateImage(cvSize(img->width, img->height), IPL_DEPTH_32F, 1);
//	fsC = cvCreateImage(cvSize(img->width, img->height), IPL_DEPTH_32F, 1);
//	fsD = cvCreateImage(cvSize(img->width, img->height), IPL_DEPTH_32F, 1);
//	fsE = cvCreateImage(cvSize(img->width, img->height), IPL_DEPTH_32F, 1);
//	fsF = cvCreateImage(cvSize(img->width, img->height), IPL_DEPTH_32F, 1);
//
//	// Compute log image
//	cvConvert( img, fB );
//	cvLog( fB, fA );
//
//	// Normalize according to given weights
//	for (i = 0, weight = 0; i < scales; i++)
//		weight += weights[i];
//
//	if (weight != 1.0) cvScale( fA, fA, weight );
//
//	// Filter at each scale
//	for (i = 0; i < scales; i++) {
//		A = cvCloneImage( img );
//		FastFilter( A, sigmas[i] );
//
//		cvConvert( A, fB );
//		cvLog( fB, fC );
//		cvReleaseImage( &A );
//
//		// Compute weighted difference
//		cvScale( fC, fC, weights[i] );
//		cvSub( fA, fC, fA );
//	}
//
//	// Color restoration
//	if (img->nChannels > 1) {
//		A = cvCreateImage(cvSize(img->width, img->height), img->depth, 1);
//		B = cvCreateImage(cvSize(img->width, img->height), img->depth, 1);
//		C = cvCreateImage(cvSize(img->width, img->height), img->depth, 1);
//
//		// Divide image into channels, convert and store sum
//		cvCvtPixToPlane( img, A, B, C, NULL );
//
//		cvConvert( A, fsA );
//		cvConvert( B, fsB );
//		cvConvert( C, fsC );
//
//		cvReleaseImage( &A );
//		cvReleaseImage( &B );
//		cvReleaseImage( &C );
//
//		// Sum components
//		cvAdd( fsA, fsB, fsD );
//		cvAdd( fsD, fsC, fsD );
//
//		// Normalize weights
//		cvDiv( fsA, fsD, fsA, restoration_factor);
//		cvDiv( fsB, fsD, fsB, restoration_factor);
//		cvDiv( fsC, fsD, fsC, restoration_factor);
//
//		cvConvertScale( fsA, fsA, 1, 1 );
//		cvConvertScale( fsB, fsB, 1, 1 );
//		cvConvertScale( fsC, fsC, 1, 1 );
//
//		// Log weights
//		cvLog( fsA, fsA );
//		cvLog( fsB, fsB );
//		cvLog( fsC, fsC );
//
//		// Divide retinex image, weight accordingly and recombine
//		cvCvtPixToPlane( fA, fsD, fsE, fsF, NULL );
//
//		cvMul( fsD, fsA, fsD, color_gain);
//		cvMul( fsE, fsB, fsE, color_gain );
//		cvMul( fsF, fsC, fsF, color_gain );
//
//		cvCvtPlaneToPix( fsD, fsE, fsF, NULL, fA );
//	}
//
//	// Restore
//	cvConvertScale( fA, img, gain, offset);
//
//	// Release temp images
//	cvReleaseImage( &fA );
//	cvReleaseImage( &fB );
//	cvReleaseImage( &fC );
//	cvReleaseImage( &fsA );
//	cvReleaseImage( &fsB );
//	cvReleaseImage( &fsC );
//	cvReleaseImage( &fsD );
//	cvReleaseImage( &fsE );
//	cvReleaseImage( &fsF );
//}
