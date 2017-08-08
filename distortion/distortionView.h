// distortionView.h : CdistortionView 类的接口
//


#pragma once


// CdistortionView 窗口

class CdistortionView : public CWnd
{
// 构造
public:
	CdistortionView();

// 属性
public:

// 操作
public:

// 重写
	protected:
	virtual BOOL PreCreateWindow(CREATESTRUCT& cs);

// 实现
public:
	virtual ~CdistortionView();

	// 生成的消息映射函数
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

