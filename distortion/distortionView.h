// distortionView.h : CdistortionView ��Ľӿ�
//


#pragma once


// CdistortionView ����

class CdistortionView : public CWnd
{
// ����
public:
	CdistortionView();

// ����
public:

// ����
public:

// ��д
	protected:
	virtual BOOL PreCreateWindow(CREATESTRUCT& cs);

// ʵ��
public:
	virtual ~CdistortionView();

	// ���ɵ���Ϣӳ�亯��
protected:
	afx_msg void OnPaint();
	DECLARE_MESSAGE_MAP()
public:
	void initcapture(void);
	void showimage(void);
	afx_msg int OnCreate(LPCREATESTRUCT lpCreateStruct);
	afx_msg void OnTimer(UINT_PTR nIDEvent);
	afx_msg void OnLButtonDown(UINT nFlags, CPoint point);
	afx_msg void OnChar(UINT nChar, UINT nRepCnt, UINT nFlags);
		IplImage* img1;
	IplImage* img2;
	IplImage* img3;
	IplImage* img4;
	IplImage* img5;
	IplImage* img6;
	IplImage* img7;
	IplImage* m;
	IplImage* temp;
};

