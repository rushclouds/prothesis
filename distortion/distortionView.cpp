// distortionView.cpp : CdistortionView 类的实现
//

#include "stdafx.h"
#include "distortion.h"
#include "distortionView.h"
#include "CvvImage.h"
#include "Reading.h"
#include "math.h"
#include <cv.h>
#include <highgui.h>
#include <cxcore.h>
#include "time.h"
#include <fstream>
#include<iostream>
#include <functional>
#include <stdio.h>
#include <opencv2/opencv.hpp>
#include "Infrare.h"
using namespace std;
#ifdef _DEBUG
#define new DEBUG_NEW
#endif
extern CdistortionApp theApp;
double a[25]={-1,0.5,0.8,0.6,0.7,0.5,0.6,-1,0.7,0.8,0.6,0.7,0.5,0.8,-1,0.7,0.8,0.6,-1,0.5,0.8,-1,0.7,0.5,0.6};  
double irr=0.6;
double para=a[0];
int ParagraphNum=25;
int index=0;
int PhosDiameter=12;
//int PhosDiameter=8;
int PhosRadius=PhosDiameter/2;
//int Cd=2*PhosDiameter;     //相邻点中心距离
int Cd=cvFloor(1.5*PhosDiameter);     //相邻点中心距离
int PhosNum=24;
//int PhosNum=16;
//double Sigma=PhosDiameter/3.0;
double Sigma=3.85;
double Sigma1=irr*Cd;
int FinalSize=442;
int cutImgLen=173;//108
double cellLen=((double)cutImgLen)/PhosNum;
double BrightProp=1; //缺失率为1-BrightProp
CvRandState rng_state;
CvRandState rng_st;
// CdistortionView
VideoCapture cap;
Mat frame;


CdistortionView::CdistortionView()       //??
{
}

CdistortionView::~CdistortionView()
{
}


BEGIN_MESSAGE_MAP(CdistortionView, CWnd)   //宏开始消息映射
	ON_WM_PAINT()   //消息映射入口
	ON_WM_CREATE()
	ON_WM_TIMER()
	ON_WM_LBUTTONDOWN()
	ON_WM_CHAR()
END_MESSAGE_MAP()



// CdistortionView 消息处理程序

BOOL CdistortionView::PreCreateWindow(CREATESTRUCT& cs) //在创建窗口前，可以设置一些窗口的属性（窗口大小，位置，工具栏，状态栏）,在OnCreate之前调用PreCreate
{
	if (!CWnd::PreCreateWindow(cs))
		return FALSE;

//	cs.dwExStyle |= WS_EX_CLIENTEDGE;
	cs.style &= ~WS_BORDER;
	cs.lpszClass = AfxRegisterWndClass(CS_HREDRAW|CS_VREDRAW|CS_DBLCLKS, 
		::LoadCursor(NULL, IDC_ARROW), reinterpret_cast<HBRUSH>(COLOR_WINDOW+1), NULL);

	return TRUE;
}



void CdistortionView::OnPaint() 
{
	CPaintDC dc(this); // 用于绘制的设备上下文
	
	// TODO: 在此处添加消息处理程序代码
	CRect   rect;  
	GetClientRect(rect); 
	dc.FillSolidRect(rect,RGB(0,0,0));

	initcapture();

	// 不要为绘制消息而调用 CWnd::OnPaint()
}


void CdistortionView::initcapture(void)
{
	theApp.capture = NULL;  
	theApp.tempFrame=cvCreateImage(cvSize(cutImgLen,cutImgLen),IPL_DEPTH_8U,1);
    theApp.outFrame=cvCreateImage(cvSize(Cd*PhosNum,Cd*PhosNum),IPL_DEPTH_8U,1);
	theApp.currentFrame=cvCreateImage(cvSize(FinalSize,FinalSize),IPL_DEPTH_8U,1);
	theApp.inFrame=cvCreateImage(cvSize(FinalSize,FinalSize),IPL_DEPTH_32F,1);
	theApp.phosphene=cvCreateImage(cvSize(PhosDiameter,PhosDiameter),IPL_DEPTH_8U,1);           
	theApp.phospheneGreytwo=cvCreateImage(cvSize(PhosDiameter,PhosDiameter),IPL_DEPTH_8U,1);
	theApp.phospheneGreythree=cvCreateImage(cvSize(PhosDiameter,PhosDiameter),IPL_DEPTH_8U,1);
	theApp.phospheneGreyfour=cvCreateImage(cvSize(PhosDiameter,PhosDiameter),IPL_DEPTH_8U,1);
	theApp.phospheneGreyfive=cvCreateImage(cvSize(PhosDiameter,PhosDiameter),IPL_DEPTH_8U,1);
	theApp.phospheneGreysix=cvCreateImage(cvSize(PhosDiameter,PhosDiameter),IPL_DEPTH_8U,1);
	theApp.phospheneGreyseven=cvCreateImage(cvSize(PhosDiameter,PhosDiameter),IPL_DEPTH_8U,1);
	//theApp.phospheneGreysix=cvCreateImage(cvSize(8,8),IPL_DEPTH_8U,1);
    //theApp.phospheneGreyseven=cvCreateImage(cvSize(8,8),IPL_DEPTH_8U,1);
	cvZero(theApp.phosphene);                                                                         
	cvZero(theApp.phospheneGreytwo);
	cvZero(theApp.phospheneGreythree);
	cvZero(theApp.phospheneGreyfour);
	cvZero(theApp.phospheneGreyfive);
	cvZero(theApp.phospheneGreysix);
	cvZero(theApp.phospheneGreyseven);
    generatePhosphene(theApp.phosphene,PhosRadius,Sigma);                                       ////生成7个灰度级
	generatePhospheneGreytwo( theApp.phospheneGreytwo,PhosRadius, Sigma);
	generatePhospheneGreythree( theApp.phospheneGreythree,PhosRadius,Sigma);
	generatePhospheneGreyfour( theApp.phospheneGreyfour,PhosRadius,Sigma);
	generatePhospheneGreyfive( theApp.phospheneGreyfive,PhosRadius,Sigma);
	generatePhospheneGreysix( theApp.phospheneGreysix,PhosRadius,Sigma);
	generatePhospheneGreyseven( theApp.phospheneGreyseven,PhosRadius,Sigma);

	theApp.disLoc=genArr(PhosNum,PhosRadius);     //生成阵列
	theApp.ElpLoc=cvCreateMat(PhosNum,PhosNum,CV_32F);
	cvRandInit( &rng_state,0,1,CV_RAND_UNI);    //(&rng_state,0,1,CV_RAND_UNI);         //2013.3.19;3.25修改
	cvRandArr(&rng_state.state,theApp.ElpLoc,CV_RAND_UNI,cvRealScalar(0),cvRealScalar(1) );
	cap.open("rtsp://192.168.1.12:554/user=admin&password=&channel=1&stream=1.sdp?");
	//cap.open(0);
	cap.read(frame);
   /////// theApp.inFrame = cvCloneImage(&(IplImage)frame); //读取一帧图像
	//theApp.inFrame = cvQueryFrame(theApp.capture);
	//setROI(frame,cutImgLen);              //视角转换，截取中间cutImgLen×cutImgLen大小的图像
	//cvCvtColor(theApp.inFrame,theApp.tempFrame0,CV_BGR2GRAY);   //灰度转换   
	theApp.hDC =::GetWindowDC(m_hWnd);             //??返回窗口设备环境获得的设备环境。覆盖了整个窗口（包括非客户区），例如标题栏、菜单、滚动条，以及边框。这使得程序能够在非客户区域实现自定义图形，例如自定义标题或者边框。
	GetClientRect( &theApp.rect );                 //取得窗口客户区(不包括非客户区)在客户区坐标系下的RECT坐标,可以得到窗口的大小
	int rw = 640;                                  // 求出图片控件的宽和高
    int rh = 480;
    int iw = FinalSize;                            // 读取图片的宽和高
    int ih = FinalSize;
    int tx = (int)(rw - iw)/2;                     // 使图片的显示位置正好在控件的正中
    int ty = (int)(rh - ih)/2;
    SetRect( theApp.rect, tx, ty, tx+iw, ty+ih );
	theApp.hDC =::GetWindowDC(m_hWnd);
	SetTimer(1,30,NULL);
}

void CdistortionView::showimage(void)
{
    theApp.m_CvvImage.CopyOf( theApp.currentFrame );                            // 复制图片
    theApp.m_CvvImage.DrawToHDC( theApp.hDC, &theApp.rect );                    // 将图片绘制到显示控件的指定区域内

}

int CdistortionView::OnCreate(LPCREATESTRUCT lpCreateStruct)
{
	if (CWnd::OnCreate(lpCreateStruct) == -1)
		return -1;
	return 0;
}

void CdistortionView::OnTimer(UINT_PTR nIDEvent)
{
	// TODO: 在此添加消息处理程序代码和/或调用默认值
	double average;
	//double average2;
	short* ptr;
	int i,j;
	// int nFrmNum = 0;
	//CString strSaveImagePath = "D:\\htt\\result\\";   //!保存图片的目录名
    //double g,M;
	 cap.read(frame);
	cvtColor(frame, frame, CV_RGB2GRAY);	
	//GaussianBlur(frame, frame, Size(3, 3), 0, 0);
	
	//下面为红外处理
	
	
	//1.同态

	//Mat frame_pro(frame);
	//my_HomoFilter(frame, frame_pro); //运行同态得把上面转灰度的语句注释
	//frame_pro = frame_pro(Rect(58, 90, 173 , 173));
	//frame_pro.convertTo(frame_pro, CV_32FC3, 1.0 / 255);
	//IplImage* m=cvCloneImage(&(IplImage) frame_pro);

	//2.双平台

	//frame = frame(Rect(58, 90, 173 , 173));
	//Mat frame_pro(frame);
	//my_Doublepla(frame, frame_pro);
	//frame_pro.convertTo(frame_pro, CV_32FC3, 1.0 / 255);
	//IplImage* m=cvCloneImage(&(IplImage) frame_pro);

		
    //3.retinex
	//当不需要进行增强而直接DP时，执行下面三句即可
	frame = frame(Rect(58, 90, 173 , 173));
	frame.convertTo(frame, CV_32FC3, 1.0 / 255);
	IplImage* m=cvCloneImage(&(IplImage) frame);
	//my_Retinex(m,15,128,128);//retinex处理，需要CV_32FC3格式才能处理
	
			cvZero(theApp.outFrame);

			for (i=0; i<PhosNum; i++)
			{
			ptr=(short*)(theApp.disLoc->data.s+i*(theApp.disLoc->step/2));
			for(j=0; j<PhosNum; j++)
			{

				average = getAverage(m, cvRect(j*cellLen, i*cellLen, cellLen, cellLen)) * 255;
			float a=cvGetReal2D(theApp.ElpLoc,i,j);
			if ((average>238) &&(a<BrightProp))
			drawPhosphene(theApp.outFrame,theApp.phosphene,cvPoint(ptr[2*j],ptr[2*j+1]));
			else if ((average>202)&&(a<BrightProp))
			drawPhosphene(theApp.outFrame,theApp.phospheneGreytwo,cvPoint(ptr[2*j],ptr[2*j+1]));
			else if((average>166)&&(a<BrightProp))
			drawPhosphene(theApp.outFrame,theApp.phospheneGreythree,cvPoint(ptr[2*j],ptr[2*j+1]));
			else if ((average>129)&&(a<BrightProp))
			drawPhosphene(theApp.outFrame,theApp.phospheneGreyfour,cvPoint(ptr[2*j],ptr[2*j+1]));
			else if ((average>92)&&(a<BrightProp))
			drawPhosphene(theApp.outFrame,theApp.phospheneGreyfive,cvPoint(ptr[2*j],ptr[2*j+1]));
			else if ((average>56)&&(a<BrightProp))
			drawPhosphene(theApp.outFrame,theApp.phospheneGreysix,cvPoint(ptr[2*j],ptr[2*j+1]));
			else if ((average>18)&&(a<BrightProp))
			drawPhosphene(theApp.outFrame,theApp.phospheneGreyseven,cvPoint(ptr[2*j],ptr[2*j+1]));
			}
			}
			cvResize(theApp.outFrame,theApp.currentFrame,CV_INTER_LINEAR);      //调整像素数
			//showimage();//显示图像
			theApp.m_CvvImage.CopyOf(theApp.currentFrame);        //theApp.currentFrame                    // 复制图片
		theApp.m_CvvImage.DrawToHDC(theApp.hDC, &theApp.rect);
	
   /////////////////释放内存////////////////////
	/*cvReleaseImage(&img1);
	cvReleaseImage(&img2);
	cvReleaseImage(&img3);
	cvReleaseImage(&img4);
	cvReleaseImage(&img5);
    cvReleaseImage(&img6);
	cvReleaseImage(&img7);
	cvReleaseImage(&temp);
	cvReleaseImage(&temp0);*/
	cvReleaseImage(&m);
	CWnd::OnTimer(nIDEvent);
	
	}

void CdistortionView::OnLButtonDown(UINT nFlags, CPoint point)
{
	// TODO: 在此添加消息处理程序代码和/或调用默认值
//	theApp.disLoc=genDisArr(PhosNum,PhosRadius,Sigma1); 
//	theApp.neiLoc=genNeiArr(theApp.disLoc,para,PhosNum,PhosRadius);
//	CWnd::OnLButtonDown(nFlags, point);
}

void CdistortionView::OnChar(UINT nChar, UINT nRepCnt, UINT nFlags)
{
	// TODO: 在此添加消息处理程序代码和/或调用默认值
	if(nChar==VK_ESCAPE)
	{
		::PostMessage(AfxGetMainWnd()->GetSafeHwnd(),WM_CLOSE,0,0);
	}
	else if (nChar==13)
	{
		index++;
		if (index<ParagraphNum)
		{
			theApp.disLoc=genDisArr(PhosNum,PhosRadius,Sigma1);
			para=a[index];
			if (para>0)
				theApp.neiLoc=genNeiArr(theApp.disLoc,para,PhosNum,PhosRadius);
			else
				theApp.neiLoc=theApp.disLoc;
		}
		else
		{
			::PostMessage(AfxGetMainWnd()->GetSafeHwnd(),WM_CLOSE,0,0);
		}
	}
	else if (nChar==8)
	{
		if(index>0)
		{
			theApp.disLoc=genDisArr(PhosNum,PhosRadius,Sigma1);
			index--;
			para=a[index];
			if (para>0)
				theApp.neiLoc=genNeiArr(theApp.disLoc,para,PhosNum,PhosRadius);
			else
				theApp.neiLoc=theApp.disLoc;
		}
	}
	CWnd::OnChar(nChar, nRepCnt, nFlags);
}
