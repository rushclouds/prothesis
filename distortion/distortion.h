// distortion.h : distortion Ӧ�ó������ͷ�ļ�
//
#pragma once

#ifndef __AFXWIN_H__
	#error "�ڰ������ļ�֮ǰ������stdafx.h�������� PCH �ļ�"
#endif

#include "resource.h"       // ������
#include "CvvImage.h"
#include <cv.h>
#include <highgui.h>
#include <cxcore.h>

// CdistortionApp:
// �йش����ʵ�֣������ distortion.cpp
//

class CdistortionApp : public CWinApp
{
public:
	CdistortionApp();


// ��д
public:
	virtual BOOL InitInstance();
	virtual int ExitInstance();

// ʵ��

public:
	afx_msg void OnAppAbout();
	DECLARE_MESSAGE_MAP()

//����
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