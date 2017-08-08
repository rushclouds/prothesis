#include <cv.h>
#include <highgui.h>
#include <cxcore.h>
void generatePhosphene( IplImage* img,int PhosRadius,double Sigma);
void drawPhosphene( IplImage* img, IplImage* phos, CvPoint center);
void generatePhospheneGreytwo( IplImage* img,int PhosRadius,double Sigma);
void generatePhospheneGreythree( IplImage* img,int PhosRadius,double Sigma);
void generatePhospheneGreyfour( IplImage* img,int PhosRadius,double Sigma);
void generatePhospheneGreyfive( IplImage* img,int PhosRadius,double Sigma);
void generatePhospheneGreysix( IplImage* img,int PhosRadius,double Sigma);
void generatePhospheneGreyseven( IplImage* img,int PhosRadius,double Sigma);
void setROI(IplImage* src,int length);
void cvThin( IplImage* src, IplImage* dst, int iterations);
void AdaptiveFindThreshold(const CvArr* image, double *low, double *high, int aperture_size=3);
void _AdaptiveFindThreshold(CvMat *dx, CvMat *dy, double *low, double *high);
CvMat* genDisArr(int rsize,int PhosRadius,double Sigma1);
CvMat* genArr(int rsize,int PhosRadius);
CvMat* genNeiArr(CvMat* disLoc,double para,int rsize,int PhosRadius);
double getAverage(IplImage* img,CvRect rect);
void TowThred(IplImage* img0,IplImage* img1,int thredl,int thredh);
void lhMorpRemoveBoderObj(const IplImage* src, IplImage* dst);
void lhMorpRDilate(const IplImage* src, const IplImage* msk, IplImage* dst, IplConvKernel* se = NULL, int iterations=-1);
int  lhImageCmp(const IplImage* img1, const IplImage* img2);