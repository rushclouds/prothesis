// distortionView.cpp : CdistortionView 类的实现
//

#include "stdafx.h"
#include "distortion.h"
#include "distortionView.h"
#include "CvvImage.h"
#include "Reading.h"
#include "Saliency.h"
#include "math.h"
#include <cv.h>
#include <highgui.h>
#include <cxcore.h>
#include "time.h"
#include <fstream>



#ifdef _DEBUG
#define new DEBUG_NEW
#endif

extern CdistortionApp theApp;
const char* Saliency::SAL_TYPE_DES[SAL_TYPE_NUM] = {"_HC", "_HC", "_HC", "_HC", "_HC"};
//const char* Saliency::SAL_TYPE_DES[SAL_TYPE_NUM] = {"_RC", "_RC", "_RC", "_RC", "_RC"};
const Saliency::GET_SAL_FUNC Saliency::gFuns[SAL_TYPE_NUM] ={GetHC, GetHC, GetHC, GetHC, GetHC};
double a[25]={-1,0.5,0.8,0.6,0.7,0.5,0.6,-1,0.7,0.8,0.6,0.7,0.5,0.8,-1,0.7,0.8,0.6,-1,0.5,0.8,-1,0.7,0.5,0.6};  //??
double irr=0.6;
double para=a[0];
int ParagraphNum=25;
int index=0;
int PhosDiameter=12;
//int PhosDiameter=8;
int PhosRadius=PhosDiameter/2;
//int Cd=2*PhosDiameter;     //相邻点中心距离
int Cd=cvFloor(1.5*PhosDiameter);     //相邻点中心距离
int PhosNum=32;
//double Sigma=PhosDiameter/3.0;
double Sigma=3.85;
double Sigma1=irr*Cd;
int FinalSize=462;
int cutImgLen=250;//108
double cellLen=((double)cutImgLen)/PhosNum;
double BrightProp=1; //缺失率为1-BrightProp
CvRandState rng_state;
CvRandState rng_st;
// CdistortionView




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
//	if (para>0)
//		theApp.neiLoc=genNeiArr(theApp.disLoc,para,PhosNum,PhosRadius);     //生成近邻搜索阵列
//	else
//		theApp.neiLoc=theApp.disLoc;
	
	theApp.capture=cvCreateCameraCapture(0);      //捕获视频  
	//theApp.capture =  cvCreateFileCapture("I:\\htt.wmv");//播放视频 就换成这句
       if(!theApp.capture)  
        {  
            AfxMessageBox("读取摄像头失败！");  
            return;  
        }

	theApp.inFrame = cvQueryFrame(theApp.capture); //读取一帧图像
	setROI(theApp.inFrame,cutImgLen);              //视角转换，截取中间cutImgLen×cutImgLen大小的图像
	
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
	SetTimer(1,40,NULL);
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
	double average2;
	short* ptr;
	int i,j;
	 int nFrmNum = 0;
	CString strSaveImagePath = "D:\\htt\\result\\";   //!保存图片的目录名

	double g,M;
	//for(;;)
 //{
	
	theApp.inFrame = cvQueryFrame(theApp.capture);
	//theApp.inFrame = cvLoadImage("C:\\htt.wmv");
	 
cvSmooth(theApp.inFrame,theApp.inFrame,CV_GAUSSIAN,3,3,0,0);//3x3
	IplImage* temp0=cvCloneImage(theApp.inFrame);
	IplImage* img2=cvCreateImage(cvGetSize(temp0),IPL_DEPTH_8U,1);
    IplImage* img3=cvCreateImage(cvGetSize(temp0),IPL_DEPTH_8U,1);
	IplImage* temp=cvCreateImage(cvGetSize(theApp.tempFrame),theApp.tempFrame->depth,theApp.tempFrame->nChannels);

cvCvtColor(temp0,img2,CV_BGR2GRAY);
Mat sal, img3f (theApp.inFrame,0) ;
img3f.convertTo(img3f, CV_32FC3, 1.0/255);
sal =Saliency::gFuns[1](img3f);

IplImage* m=cvCloneImage(&(IplImage) sal);
IplImage* img1=cvCreateImage(cvGetSize(temp0),IPL_DEPTH_8U,1);
   
IplImage* img4=cvCreateImage(cvGetSize(m),m->depth,m->nChannels);
IplImage* img5=cvCreateImage(cvGetSize(m),m->depth,m->nChannels);
IplImage* img6=cvCreateImage(cvGetSize(m),m->depth,m->nChannels);
IplImage* img7=cvCreateImage(cvGetSize(m),m->depth,m->nChannels);		

////////////改进的canny///////////////////////
double low = 0.0, high = 0.0;  
AdaptiveFindThreshold(img2, &low, &high,3);
 
cvErode(img2,img3,NULL,1);
cvSub(img2,img3,img3);
cvSmooth(img3,img3,CV_GAUSSIAN,3,3,0,0);//3x3
 cvScale(img3,img3,20);

 ////////////数据类型转换/////////////////



cvMinMaxLoc(img3, &g, &M, NULL, NULL, NULL);  
cvScale(img3, img7, 1.0/(M-g), 1.0*(-g)/(M-g));//图像数据转换到[0,1]区间
cvSmooth(img7,img7,CV_GAUSSIAN,3,3,0,0);//3x3 
float *pchar;

 int width=m->width;
 int heigh=m->height;
 int step=m->widthStep/sizeof(float);
 int step1=img7->widthStep/sizeof(float);
  int step2=img6->widthStep/sizeof(float);
float *data=(float*)m->imageData;
float *data1=(float*)img7->imageData;
float *data2=(float*)img6->imageData;
float temp1;
float temp2;
 for (i=0;i<heigh;i++)
 {
 for(j=0;j<width;j++)
	 {
	 temp1=data1[i*step1+j];
		if (temp1>0.55)
		//temp1=temp1*temp2;
data1[i*step1+j]=temp1;
		else
		{data1[i*step1+j]=0;}
 temp2=data[i*step+j]+data1[i*step1+j];
	  data2[i*step2+j]=temp2;
		 
	 } 
	 }

 cvMorphologyEx(img6,img4,temp,NULL,CV_MOP_OPEN,1);
cvMorphologyEx(img4,img5,temp,NULL,CV_MOP_CLOSE,1);
  cvSmooth(img5,img5,CV_GAUSSIAN,3,3,0,0);//3x3

	cvZero(theApp.outFrame);
	
	for (i=0; i<PhosNum; i++)
	{
		ptr=(short*)(theApp.disLoc->data.s+i*(theApp.disLoc->step/2));
		for(j=0; j<PhosNum; j++)
		{
			
			average=getAverage(m,cvRect(j*cellLen,i*cellLen,cellLen,cellLen))*255;
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
theApp.m_CvvImage.CopyOf(theApp.currentFrame);                            // 复制图片
theApp.m_CvvImage.DrawToHDC( theApp.hDC, &theApp.rect );	

	/////////////////保存每一帧图像为图片形式///////////////
	  /* CString strSaveImageName;
       CString strTmp;
       strTmp.Format("SaveImage_%d.jpg",nFrmNum);
       strSaveImageName = strSaveImagePath + strTmp;
       cvSaveImage (strSaveImageName,img3,0);
		nFrmNum++;*/
   /////////////////释放内存////////////////////
	cvReleaseImage(&img1);
	cvReleaseImage(&img2);
	cvReleaseImage(&img3);
	cvReleaseImage(&img4);
	cvReleaseImage(&img5);
    cvReleaseImage(&img6);
	cvReleaseImage(&img7);
	cvReleaseImage(&temp);
	cvReleaseImage(&temp0);
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
Mat Saliency::GetHC(const Mat &img3f)
{
	// Quantize colors and
	Mat idx1i, binColor3f, colorNums1i, weight1f, _colorSal;
	Quantize(img3f, idx1i, binColor3f, colorNums1i);
	cvtColor(binColor3f, binColor3f, CV_BGR2Lab);

	normalize(colorNums1i, weight1f, 1, 0, NORM_L1, CV_32F);
	GetHC(binColor3f, weight1f, _colorSal);
	float* colorSal = (float*)(_colorSal.data);
	Mat salHC1f(img3f.size(), CV_32F);
	for (int r = 0; r < img3f.rows; r++)
	{
		float* salV = salHC1f.ptr<float>(r);
		int* _idx = idx1i.ptr<int>(r);
		for (int c = 0; c < img3f.cols; c++)
			salV[c] = colorSal[_idx[c]];
	}
	GaussianBlur(salHC1f, salHC1f, Size(3, 3), 0);
	normalize(salHC1f, salHC1f, 0, 1, NORM_MINMAX);
	return salHC1f;
}


 void Saliency::GetHC(const Mat &binColor3f, const Mat &weight1f, Mat &_colorSal)
{
	int binN = binColor3f.cols; 
	_colorSal = Mat::zeros(1, binN, CV_32F);
	float* colorSal = (float*)(_colorSal.data);
	vector<vector<CostfIdx>> similar(binN); // Similar color: how similar and their index
	Vec3f* color = (Vec3f*)(binColor3f.data);
	float *w = (float*)(weight1f.data);
	for (int i = 0; i < binN; i++)
	{
		vector<CostfIdx> &similari = similar[i];
		similari.push_back(make_pair(0.f, i));
		for (int j = 0; j < binN; j++)
		{
			if (i == j)
				continue;
			float dij = vecDist3<float>(color[i], color[j]);
			similari.push_back(make_pair(dij, j));
			colorSal[i] += w[j] * dij;
		}
		sort(similari.begin(), similari.end());
	}

	SmoothSaliency(binColor3f, _colorSal, 4.0f, similar);
}

void Saliency::SmoothSaliency(const Mat &binColor3f, Mat &sal1d, float delta, const vector<vector<CostfIdx>> &similar)
{
	if (sal1d.cols < 2)
		return;
	CV_Assert(binColor3f.size() == sal1d.size() && sal1d.rows == 1);
	int binN = binColor3f.cols;
	Vec3f* color = (Vec3f*)(binColor3f.data);
	Mat tmpSal;
	sal1d.copyTo(tmpSal);
	float *sal = (float*)(tmpSal.data);
	float *nSal = (float*)(sal1d.data);

	//* Distance based smooth
	int n = max(cvRound(binN/delta), 2);
	vecF dist(n, 0), val(n);
	for (int i = 0; i < binN; i++)
	{
		const vector<CostfIdx> &similari = similar[i];
		float totalDist = 0;

		val[0] = sal[i];
		for (int j = 1; j < n; j++)
		{
			int ithIdx =similari[j].second;
			dist[j] = similari[j].first;
			val[j] = sal[ithIdx];
			totalDist += dist[j];
		}
		float valCrnt = 0;
		for (int j = 0; j < n; j++)
			valCrnt += val[j] * (totalDist - dist[j]);

		nSal[i] =  valCrnt / ((n-1) * totalDist);
	}	
	

	  /* Gaussian smooth  */
	const float guassCoeff = -0.5f/(delta*delta);
	for (int i = 0; i < binN; i++)
	{
		const vector<CostfIdx> &similari = similar[i];
		float saliencyI = sal[i], totalW = 1;

		for (int j = 1; j < binN; j++)
		{
			float w = expf(sqr(similari[j].first)*guassCoeff);
			if (w < 1e-8f)
				break;
			saliencyI += w * sal[similari[j].second];
			totalW += w;
		}
		nSal[i] = saliencyI / totalW;
	}
	//
}


int Saliency::Quantize(const Mat& img3f, Mat &idx1i, Mat &_color3f, Mat &_colorNum, double ratio)
{
	static const int clrNums[3] = {12, 12, 12};
	static const float clrTmp[3] = {clrNums[0] - 0.0001f, clrNums[1] - 0.0001f, clrNums[2] - 0.0001f};
	static const int w[3] = {clrNums[1] * clrNums[2], clrNums[2], 1};

	CV_Assert(img3f.data != NULL);
	idx1i = Mat::zeros(img3f.size(), CV_32S);
	int rows = img3f.rows, cols = img3f.cols;
	if (img3f.isContinuous() && idx1i.isContinuous())//  如果内存连续存放为一行 则rows=1
	{
		cols *= rows;
		rows = 1;
	}

	// Build color pallet
	map<int, int> pallet;//Map是c++的一个标准容器
	for (int y = 0; y < rows; y++)
	{
		const float* imgData = img3f.ptr<float>(y);//ptr是opencv的智能指针
		int* idx = idx1i.ptr<int>(y);
		for (int x = 0; x < cols; x++, imgData += 3)
		{
			idx[x] = (int)(imgData[0]*clrTmp[0])*w[0] + (int)(imgData[1]*clrTmp[1])*w[1] + (int)(imgData[2]*clrTmp[2]);
			pallet[idx[x]] ++;
		}
	}

	// Fine significant colors
	int maxNum = 0;
	{
		int count = 0;
		vector<pair<int, int>> num; // (num, color) pairs in num
		num.reserve(pallet.size());
		for (map<int, int>::iterator it = pallet.begin(); it != pallet.end(); it++)
			num.push_back(pair<int, int>(it->second, it->first)); // (color, num) pairs in pallet
		sort(num.begin(), num.end(), std::greater<pair<int, int>>());

		maxNum = (int)num.size();
		int maxDropNum = cvRound(rows * cols * (1-ratio));
		for (int crnt = num[maxNum-1].first; crnt < maxDropNum && maxNum > 1; maxNum--)
			crnt += num[maxNum - 2].first;
		maxNum = min(maxNum, 256); // To avoid very rarely case
		if (maxNum < 10)
			maxNum = min((int)pallet.size(), 100);
		pallet.clear();
		for (int i = 0; i < maxNum; i++)
			pallet[num[i].second] = i; 

		vector<Vec3i> color3i(num.size());
		for (unsigned int i = 0; i < num.size(); i++)
		{
			color3i[i][0] = num[i].second / w[0];
			color3i[i][1] = num[i].second % w[0] / w[1];
			color3i[i][2] = num[i].second % w[1];
		}

		for (unsigned int i = maxNum; i < num.size(); i++)//剩余5%用最近距离代替
		{
			int simIdx = 0, simVal = INT_MAX;
			for (int j = 0; j < maxNum; j++)
			{
				int d_ij = vecSqrDist3(color3i[i], color3i[j]);
				if (d_ij < simVal)
					simVal = d_ij, simIdx = j;
			}
			pallet[num[i].second] = pallet[num[simIdx].second];
		}
	}

	_color3f = Mat::zeros(1, maxNum, CV_32FC3);
	_colorNum = Mat::zeros(_color3f.size(), CV_32S);

	Vec3f* color = (Vec3f*)(_color3f.data);
	int* colorNum = (int*)(_colorNum.data);
	for (int y = 0; y < rows; y++) 
	{
		const Vec3f* imgData = img3f.ptr<Vec3f>(y);
		int* idx = idx1i.ptr<int>(y);
		for (int x = 0; x < cols; x++)
		{
			idx[x] = pallet[idx[x]];
			color[idx[x]] += imgData[x];
			colorNum[idx[x]] ++;
		}
	}
	for (int i = 0; i < _color3f.cols; i++)
		//color[i] /= colorNum[i];
		color[i]=color[i]/colorNum[i];

	return _color3f.cols;
}

