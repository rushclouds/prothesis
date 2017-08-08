#include "StdAfx.h"  
#include "Reading.h"
#include "math.h"
#include <cv.h>
#include <highgui.h>
#include <cxcore.h>
#include <fstream>


/********产生一个光幻视点********/
void generatePhosphene( IplImage* img,int PhosRadius,double Sigma)
{

	double r;
for(int y=0;y<img->height;y++)
	
	{
		uchar* ptr=(uchar*)(img->imageData+y*img->widthStep);
		for(int x=0;x<img->width;x++)
		{
			r=pow((double)(y-PhosRadius),2)+pow((double)(x-PhosRadius),2);
			if (r<=pow((double)PhosRadius,2))
			ptr[x]=cvRound(255*exp(-r/(2*pow(Sigma,2))));
		}
	}
}

/********画一个光幻视点********/
void drawPhosphene( IplImage* img, IplImage* phos, CvPoint center)
{
	if (center.x<0||center.y<0)
		return;
	int x0,y0;
	int temp;
	x0=center.x-phos->width/2;
	y0=center.y-phos->height/2;
	for(int y=0;y<phos->height;y++)
	{
		uchar* ptri=(uchar*)(img->imageData+(y+y0)*img->widthStep);
		uchar* ptrp=(uchar*)(phos->imageData+y*phos->widthStep);
		for(int x=0;x<phos->width;x++)
		{
			temp=ptri[x0+x]+ptrp[x];
			if (temp>255)
				temp=255;
			ptri[x0+x]=temp;
		}
	}
}
void generatePhospheneGreytwo( IplImage* img,int PhosRadius,double Sigma)
{
	double r;
	
	for(int y=0;y<img->height;y++)
	{
		uchar* ptr=(uchar*)(img->imageData+y*img->widthStep);
		for(int x=0;x<img->width;x++)
		{
			r=pow((double)(y-PhosRadius),2)+pow((double)(x-PhosRadius),2);
			if (r<=pow((double)PhosRadius,2))
			ptr[x]=cvRound(219*exp(-r/(2*pow(Sigma,2))));
		}
	}
}


void generatePhospheneGreythree( IplImage* img,int PhosRadius,double Sigma)
{
	double r;
	
	for(int y=0;y<img->height;y++)
	{
		uchar* ptr=(uchar*)(img->imageData+y*img->widthStep);
		for(int x=0;x<img->width;x++)
		{
			r=pow((double)(y-PhosRadius),2)+pow((double)(x-PhosRadius),2);
			if (r<=pow((double)PhosRadius,2))
			ptr[x]=cvRound(183*exp(-r/(2*pow(Sigma,2))));
		}
	}
}

void generatePhospheneGreyfour( IplImage* img,int PhosRadius,double Sigma)
{
	double r;
	
	for(int y=0;y<img->height;y++)
	{
		uchar* ptr=(uchar*)(img->imageData+y*img->widthStep);
		for(int x=0;x<img->width;x++)
		{
			r=pow((double)(y-PhosRadius),2)+pow((double)(x-PhosRadius),2);
			if (r<=pow((double)PhosRadius,2))
			ptr[x]=cvRound(146*exp(-r/(2*pow(Sigma,2))));
		}
	}
}

void generatePhospheneGreyfive( IplImage * img,int PhosRadius,double Sigma)
{
	double r;
	
	for(int y=0;y<img->height;y++)
	{
		uchar* ptr=(uchar*)(img->imageData+y*img->widthStep);
		for(int x=0;x<img->width;x++)
		{
			r=pow((double)(y-PhosRadius),2)+pow((double)(x-PhosRadius),2);
			if (r<=pow((double)PhosRadius,2))
			ptr[x]=cvRound(110*exp(-r/(2*pow(Sigma,2))));
		}
	}
}

void generatePhospheneGreysix( IplImage* img,int PhosRadius,double Sigma)
{
	double r;
	
	for(int y=0;y<img->height;y++)
	{
		uchar* ptr=(uchar*)(img->imageData+y*img->widthStep);
		for(int x=0;x<img->width;x++)
		{
			r=pow((double)(y-PhosRadius),2)+pow((double)(x-PhosRadius),2);
			if (r<=pow((double)PhosRadius,2))
			ptr[x]=cvRound(73*exp(-r/(2*pow(Sigma,2))));
		}
	}
}

/*void generatePhospheneGreyseven( IplImage* img,int PhosRadius,double Sigma)
{
	double r;
	
	for(int y=0;y<img->height;y++)
	{
		uchar* ptr=(uchar*)(img->imageData+y*img->widthStep);
		for(int x=0;x<img->width;x++)
		{
			r=pow((double)(y-4),2)+pow((double)(x-4),2);
			if (r<=pow((double)4,2))
			ptr[x]=cvRound(37*exp(-r/(2*pow(Sigma,2))));
		}
	}
}
*/
void generatePhospheneGreyseven( IplImage* img,int PhosRadius,double Sigma)
{
	double r;
	
	for(int y=0;y<img->height;y++)
	{
		uchar* ptr=(uchar*)(img->imageData+y*img->widthStep);
		for(int x=0;x<img->width;x++)
		{
			r=pow((double)(y-PhosRadius),2)+pow((double)(x-PhosRadius),2);
			if (r<=pow((double)PhosRadius,2))
			ptr[x]=cvRound(37*exp(-r/(2*pow(Sigma,2))));
		}
	}
}

/******从图像中心设置方形ROI区域******/
void setROI(IplImage* src,int length)
{
	int x,y;
	y=(src->height-length)/2+1;
	x=(src->width-length)/2+1;
	cvSetImageROI(src,cvRect(x,y,length,length));
}
/******骨架提取******/
void cvThin( IplImage* src, IplImage* dst, int iterations)//将IPL_DEPTH_8U型二值图像进行细化
{
 CvSize size = cvGetSize(src);
 cvCopy(src, dst);//拷贝一个数组给另一个数组
    int n = 0,i = 0,j = 0;
 for(n=0; n<iterations; n++)
 {
//IplImage* t_image = cvCloneImage(dst);
IplImage* t_image=cvCreateImage(size,src->depth,src->nChannels);
cvCopy(dst,t_image);
// t_image = cvCloneImage(dst);

  for(i=0; i<size.height;  i++)
  {
   for(int j=0; j<size.width; j++)
   {
    if(CV_IMAGE_ELEM(t_image,byte,i,j)==1)
    {
     int ap=0;
     int p2 = (i==0)?0:CV_IMAGE_ELEM(t_image,byte, i-1, j);
     int p3 = (i==0 || j==size.width-1)?0:CV_IMAGE_ELEM(t_image,byte, i-1, j+1);
     if (p2==0 && p3==1)
     {
      ap++;
     }
     int p4 = (j==size.width-1)?0:CV_IMAGE_ELEM(t_image,byte,i,j+1);
     if(p3==0 && p4==1)
     {
      ap++;
     }
     int p5 = (i==size.height-1 || j==size.width-1)?0:CV_IMAGE_ELEM(t_image,byte,i+1,j+1);
     if(p4==0 && p5==1)
     {
      ap++;
     }
     int p6 = (i==size.height-1)?0:CV_IMAGE_ELEM(t_image,byte,i+1,j);
     if(p5==0 && p6==1)
     {
      ap++;
     }
     int p7 = (i==size.height-1 || j==0)?0:CV_IMAGE_ELEM(t_image,byte,i+1,j-1);
     if(p6==0 && p7==1)
     {
      ap++;
     }
     int p8 = (j==0)?0:CV_IMAGE_ELEM(t_image,byte,i,j-1);
     if(p7==0 && p8==1)
     {
      ap++;
     }
     int p9 = (i==0 || j==0)?0:CV_IMAGE_ELEM(t_image,byte,i-1,j-1);
     if(p8==0 && p9==1)
     {
      ap++;
     }
     if(p9==0 && p2==1)
     {
      ap++;
     }
     if((p2+p3+p4+p5+p6+p7+p8+p9)>1 && (p2+p3+p4+p5+p6+p7+p8+p9)<7)
     {
      if(ap==1)
      {
       if(p2*p4*p8==0)
       {
        if(p2*p6*p8==0)
        {
         CV_IMAGE_ELEM(dst, byte,i,j)=0;
        }
       }
      }
     }                    
    }

   }

  }            
  cvReleaseImage(&t_image);

}}








/******产生不规则阵列坐标********/
CvMat* genDisArr(int rsize,int PhosRadius,double Sigma1)
{
	int PhosDiameter=2*PhosRadius;
	//int Cd=2*PhosDiameter;
	int Cd=PhosDiameter;
	CvMat* disLoc=cvCreateMat(rsize,rsize,CV_16SC2);

	CvRandState rng_state;
//	CvRNG rng_state = cvRNG(cvGetTickCount());
	cvRandInit(&rng_state,0,Sigma1,(int)cvGetTickCount(),CV_RAND_NORMAL);

//	cvRandInit(rng_state.state,0,Sigma1,(int)cvGetTickCount(),CV_RAND_NORMAL);
	cvRandArr( &rng_state.state, disLoc,CV_RAND_NORMAL,cvScalar(0,0,0,0),cvScalar(Sigma1,Sigma1,0,0));
//	cvRandArr( &rng_state, disLoc,CV_RAND_NORMAL,cvScalar(0,0),cvScalar(Sigma1,Sigma1));
	for (int y=0; y<rsize; y++)
	{
		short* ptr=(short*)(disLoc->data.s+y*(disLoc->step/2));
		for (int x=0; x<rsize; x++)
		{
//			ptr[2*x]=Cd*x+PhosDiameter;
//			ptr[2*x+1]=Cd*y+PhosDiameter;
//			CvPoint* ptr=(CvPoint*)cvPtr2D(disLoc,y,x,0);

			int dy=Cd*y+PhosDiameter+ptr[2*x+1];
			int dx=Cd*x+PhosDiameter+ptr[2*x];

			int max=Cd*rsize-PhosRadius;
			if (dy<PhosRadius)
				ptr[2*x+1]=PhosRadius;
			else if (dy>max)
				ptr[2*x+1]=max;
			else
				ptr[2*x+1]=dy;
			if (dx<PhosRadius)
				ptr[2*x]=PhosRadius;
			else if (dx>max)
				ptr[2*x]=max;
			else
				ptr[2*x]=dx;
		}
	}

	return disLoc;
}


/******产生规则阵列坐标********/
CvMat* genArr(int rsize,int PhosRadius)
{
	int PhosDiameter=2*PhosRadius;
	int Cd=cvFloor(1.5*PhosDiameter);
	CvMat* disLoc=cvCreateMat(rsize,rsize,CV_16SC2);
	for (int y=0; y<rsize; y++)
	{
		for (int x=0; x<rsize; x++)
		{
			short* ptr=(short*)(disLoc->data.s+y*(disLoc->step/2));
				ptr[2*x]=Cd*x+PhosDiameter;
				ptr[2*x+1]=Cd*y+PhosDiameter;
//			CvPoint pt=*(CvPoint*)cvPtr2D(disLoc,y,x,0);
//			printf("%d %d; ",pt.y, pt.x);

		}
	}

	return disLoc;
}


/******搜索邻近点阵列********/
CvMat* genNeiArr(CvMat* disLoc,double para,int rsize,int PhosRadius)
{
	int PhosDiameter=2*PhosRadius;
	int Cd=2*PhosDiameter;
	int tx=-1;
	int ty=-1;
	int ti=-1;
	int tj=-1;
	double min;
	double temp;
	bool sign;
//	min=pow((double)PhosNum,2);
	CvMat* neiLoc=cvCreateMat(rsize,rsize,CV_16SC2);
	CvMat* distemp=cvCloneMat(disLoc);
	for (int y=0; y<rsize; y++)
	{
		short* ptrn=(short*)(neiLoc->data.s+y*(neiLoc->step/2));
		for (int x=0; x<rsize; x++)
		{
			sign=false;
			min=pow(para*Cd,2);
			for(int i=0;i<rsize;i++)
			{
				short* ptrd=(short*)(distemp->data.s+i*(distemp->step/2));
				for(int j=0;j<rsize;j++)
				{
					if (ptrd[2*j+1]!=-1&&ptrd[2*j]!=-1)
					{
						temp=pow((double)ptrd[2*j+1]-Cd*y-PhosDiameter,2)+pow((double)ptrd[2*j]-Cd*x-PhosDiameter,2);
						if (temp<=min)
						{
							tx=ptrd[2*j];
							ty=ptrd[2*j+1];
							ti=i;
							tj=j;
							min=temp;
							sign=true;
						}
					}
				}
			}
			if (sign)
			{
				ptrn[2*x]=tx;
				ptrn[2*x+1]=ty;
				((short*)(distemp->data.s+ti*(distemp->step/2)))[2*tj]=-1;
				((short*)(distemp->data.s+ti*(distemp->step/2)))[2*tj+1]=-1;

			}
			else
			{
				ptrn[2*x]=-1;
				ptrn[2*x+1]=-1;
			}
		}
	}
	cvReleaseMat(&distemp);
	return neiLoc;
}

/*******求二值图像矩形区域平均灰度值*********/
double getAverage(IplImage* img,CvRect rect)
{	
	float sum=0;
	double average;
	int step=img->widthStep/sizeof(float);
float *data=(float*)img->imageData;
float temp1;
	for(int y=0;y<rect.height;y++)
	{
		//float* ptr=(float*)(img->imageData+(y+rect.y)*img->widthStep);
		for(int x=0;x<rect.width;x++)
		{
			temp1=data[x+rect.x+(y+rect.y)*step];
			sum+=temp1;
		}
	}
	average=(double)sum/(rect.width*rect.height);
	return average;

}
//double getAverage(IplImage* img,CvRect rect)
//{	
//	int sum=0;
//	double average;
//	for(int y=0;y<rect.height;y++)
//	{
//		uchar* ptr=(uchar*)(img->imageData+(y+rect.y)*img->widthStep);
//		for(int x=0;x<rect.width;x++)
//		{
//			sum+=ptr[x+rect.x];
//		}
//	}
//	average=(double)sum/(rect.width*rect.height);
//	return average;
//
//}
/*********双阈值二值化************/
void TowThred(IplImage* img0,IplImage* img1,int thredl,int thredh)
{
	for(int y=0;y<img1->height;y++)
	{ 
		uchar* ptr0=(uchar*)(img0->imageData+y*img0->widthStep);
		uchar* ptr1=(uchar*)(img1->imageData+y*img1->widthStep);
		for(int x=0;x<img1->width;x++)
		{
			if (ptr1[x]<thredl)
				ptr0[x]=1;
			else if(ptr1[x]>thredh)
				ptr0[x]=0;
		}
	}
	cvCopy(img0,img1);
}
void AdaptiveFindThreshold(const CvArr* image, double *low, double *high, int aperture_size)  
{                                                                                
    cv::Mat src = cv::cvarrToMat(image);                                     
    const int cn = src.channels();                                           
    cv::Mat dx(src.rows, src.cols, CV_16SC(cn));                             
    cv::Mat dy(src.rows, src.cols, CV_16SC(cn));                             
                                                                                 
   /* cv::Sobel(src, dx, CV_16S, 1, 0, aperture_size, 1, 0, cv::BORDER_REPLIC);  
    cv::Sobel(src, dy, CV_16S, 0, 1, aperture_size, 1, 0, cv::BORDER_REPLIC);  
                                                                                 */
	 cv::Sobel(src, dx, CV_16S, 1, 0, aperture_size, 1, 0, 4);  
    cv::Sobel(src, dy, CV_16S, 0, 1, aperture_size, 1, 0, 4);  
                                                                                 
    CvMat _dx = dx, _dy = dy;                                                
    _AdaptiveFindThreshold(&_dx, &_dy, low, high);                           
                                                                                 
}                                                                                
                                                                                 
// 仿照matlab，自适应求高低两个门限                                              
void _AdaptiveFindThreshold(CvMat *dx, CvMat *dy, double *low, double *high)     
{                                                                                
    CvSize size;                                                             
    IplImage *imge=0;                                                        
    int i,j;                                                                 
    CvHistogram *hist;                                                       
    int hist_size = 255;                                                     
    float range_0[]={0,256};                                                 
    float* ranges[] = { range_0 };                                           
    double PercentOfPixelsNotEdges = 0.7;                                    
    size = cvGetSize(dx);                                                    
    imge = cvCreateImage(size, IPL_DEPTH_32F, 1);                            
    // 计算边缘的强度, 并存于图像中                                          
    float maxv = 0;                                                          
    for(i = 0; i < size.height; i++ )                                        
    {                                                                        
        const short* _dx = (short*)(dx->data.ptr + dx->step*i);          
        const short* _dy = (short*)(dy->data.ptr + dy->step*i);          
        float* _image = (float *)(imge->imageData + imge->widthStep*i);  
        for(j = 0; j < size.width; j++)                                  
        {                                                                
            _image[j] = (float)(abs(_dx[j]) + abs(_dy[j]));          
            maxv = maxv < _image[j] ? _image[j]: maxv;               
                                                                             
        }                                                                
    }                                                                        
    if(maxv == 0){                                                           
        *high = 0;                                                       
        *low = 0;                                                        
        cvReleaseImage( &imge );                                         
        return;                                                          
    }                                                                        
                                                                                 
    // 计算直方图                                                            
    range_0[1] = maxv;                                                       
    hist_size = (int)(hist_size > maxv ? maxv:hist_size);                    
    hist = cvCreateHist(1, &hist_size, CV_HIST_ARRAY, ranges, 1);            
    cvCalcHist( &imge, hist, 0, NULL );                                      
    int total = (int)(size.height * size.width * PercentOfPixelsNotEdges);   
    float sum=0;                                                             
    int icount = hist->mat.dim[0].size;                                      
                                                                                 
    float *h = (float*)cvPtr1D( hist->bins, 0 );                             
    for(i = 0; i < icount; i++)                                              
    {                                                                        
        sum += h[i];                                                     
        if( sum > total )                                                
            break;                                                   
    }                                                                        
    // 计算高低门限                                                          
    *high = (i+1) * maxv / hist_size ;                                       
    *low = *high * 0.4;                                                      
    cvReleaseImage( &imge );                                                 
    cvReleaseHist(&hist);                                                    
}   

void lhMorpRemoveBoderObj(const IplImage* src, IplImage* dst)
{
	IplImage *temp = cvCloneImage(src);
	//double min, max;
	//cvMinMaxLoc(src, &min, &max);
	
	//标记图像
	cvRectangle(temp, cvPoint(3,3), cvPoint(temp->width-7, temp->height-7), CV_RGB(0,0,0), -1);
	//cvRectangle(temp, cvPoint(1,1), cvPoint(temp->width, temp->height), CV_RGB(0,0,0), -1);
	//将原图像作为掩模图像
	 IplConvKernel *element1 = cvCreateStructuringElementEx( 3, 3, 0, 0, CV_SHAPE_RECT, 0);
	lhMorpRDilate(temp, src, dst,NULL,-1);
	cvReleaseImage(&temp);
	cvReleaseStructuringElement(&element1);
	//cvSet((CvArr*)src, cvScalar(min), dst);
	//cvCopy(src, dst);
	cvSub(src, dst, dst);
}

void lhMorpRDilate(const IplImage* src, const IplImage* msk, IplImage* dst, IplConvKernel* se , int iterations)
{

	assert(src != NULL && msk != NULL && dst != NULL && src != dst );

	//if(iterations < 0)
	//{
		//膨胀重建
		cvMin(src, msk, dst);
		cvDilate(dst, dst, se);
		cvMin(dst, msk, dst);
		IplImage* temp1=cvCreateImage(cvGetSize(src),src->depth,src->nChannels);
		
		ofstream htt; //write只是个名字 你可以定义为任何其他的名字
htt.open("D:\\htt\\study\\material\\htt1.txt"); //表示你要把内容输出到“text.txt"这个文件里 如果没有这个文件，会自动创建这个文件
 clock_t start1 = clock();      
int i=0;
	do
		{
			cvCopy(dst, temp1);
			cvDilate(dst, dst, se);
			cvMin(dst, msk, dst);	
			//std::swap(dst, temp1);
			/*cvCopy(dst, temp1);
			cvDilate(dst, temp1, se);
			cvMin(temp1, msk, temp1);*/
			//i=i+1;
		}
		while(lhImageCmp(temp1, dst)!= 0);	

 //for(int i=0;i<100;i++)
 //{
	// 	cvCopy(dst, temp1);
	//   cvDilate(dst, dst, se);
	//	cvMin(dst, msk, dst);
 //}
 clock_t end1   = clock();  
cout<<end1 - start1<<endl;
htt << end1-start1; //这里是你想要输出的内容，这里是输出了一个string abc
//htt << i;
htt.close(); // 输出完毕后关闭这个文件
 
		//cvReleaseImage(&temp1);
		//cvReleaseImage(&temp2);

		return;	
}


////////////////比较两个图像是否相同， 0 相同

int  lhImageCmp(const IplImage* img1, const IplImage* img2)
{
	assert(img1->width == img2->width && img1->height == img2->height && img1->imageSize == img2->imageSize );
	return memcmp(img1->imageData, img2->imageData, img1->imageSize);
}

