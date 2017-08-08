// distortion.h : distortion 应用程序的主头文件
//
#pragma once

#ifndef __AFXWIN_H__
	#error "在包含此文件之前包含“stdafx.h”以生成 PCH 文件"
#endif

#include "resource.h"       // 主符号
#include "CvvImage.h"
#include <cv.h>
#include <highgui.h>
#include <cxcore.h>

// CdistortionApp:
// 有关此类的实现，请参阅 distortion.cpp
//

class CdistortionApp : public CWinApp
{
public:
	CdistortionApp();


// 重写
public:
	virtual BOOL InitInstance();
	virtual int ExitInstance();

// 实现

public:
	afx_msg void OnAppAbout();
	DECLARE_MESSAGE_MAP()

//变量
public:
	int SM_CXSCREEN_x;
	int SM_CYSCREEN_y;

	CvvImage m_CvvImage;  
	CRect rect;   
	HDC hDC;      
	CvCapture *capture;  
	IplImage *currentFrame; 
	IplImage* phosphene;
	IplImage* phospheneGreytwo;
	IplImage* phospheneGreythree;
	IplImage* phospheneGreyfour;
	IplImage* phospheneGreyfive;
	IplImage* phospheneGreysix;
	IplImage* phospheneGreyseven;
	IplImage* outFrame;
	IplImage* inFrame;
	IplImage* tempFrame;


	//IplImage* tempFrame0;

	CvMat* disLoc;
	CvMat* neiLoc;
	CvMat* ElpLoc;
};

extern CdistortionApp theApp;