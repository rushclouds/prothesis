// distortion.cpp : 定义应用程序的类行为。
//

#include "stdafx.h"
#include "distortion.h"
#include "MainFrm.h"


#ifdef _DEBUG
#define new DEBUG_NEW
#endif

// CdistortionApp

BEGIN_MESSAGE_MAP(CdistortionApp, CWinApp)                   //使用BEGIN_MESSAGE_MAP宏开始消息映射，
	ON_COMMAND(ID_APP_ABOUT, &CdistortionApp::OnAppAbout)    //然后为每个消息处理函数加入一个入口，
END_MESSAGE_MAP()                                            //最后用END_MESSAGE_MAP宏结束消息映射
                                                             //都是宏定义，不是函数

// CdistortionApp 构造

CdistortionApp::CdistortionApp()
{
	// TODO: 在此处添加构造代码，
	// 将所有重要的初始化放置在 InitInstance 中
}


// 唯一的一个 CdistortionApp 对象

CdistortionApp theApp;


// CdistortionApp 初始化

BOOL CdistortionApp::InitInstance()
{
	// 如果一个运行在 Windows XP 上的应用程序清单指定要
	// 使用 ComCtl32.dll 版本 6 或更高版本来启用可视化方式，
	//则需要 InitCommonControlsEx()。否则，将无法创建窗口。
	INITCOMMONCONTROLSEX InitCtrls;
	InitCtrls.dwSize = sizeof(InitCtrls);
	// 将它设置为包括所有要在应用程序中使用的
	// 公共控件类。
	InitCtrls.dwICC = ICC_WIN95_CLASSES;
	InitCommonControlsEx(&InitCtrls);

	CWinApp::InitInstance();

	// 标准初始化
	// 如果未使用这些功能并希望减小
	// 最终可执行文件的大小，则应移除下列
	// 不需要的特定初始化例程
	// 更改用于存储设置的注册表项
	// TODO: 应适当修改该字符串，
	// 例如修改为公司或组织名
	SetRegistryKey(_T("应用程序向导生成的本地应用程序"));
	// 若要创建主窗口，此代码将创建新的框架窗口
	// 对象，然后将其设置为应用程序的主窗口对象
	CMainFrame* pFrame = new CMainFrame;
	if (!pFrame)
		return FALSE;
	m_pMainWnd = pFrame;
	// 创建并加载框架及其资源
	pFrame->LoadFrame(IDR_MAINFRAME,
		WS_OVERLAPPEDWINDOW | FWS_ADDTOTITLE, NULL,
		NULL);






	// 唯一的一个窗口已初始化，因此显示它并对其进行更新
	pFrame->ShowWindow(SW_SHOW);
	pFrame->UpdateWindow();
	// 仅当具有后缀时才调用 DragAcceptFiles
	//  在 SDI 应用程序中，这应在 ProcessShellCommand 之后发生

	/**************设置分辨率**************/
	//保存原始屏幕分辨率，便于退出后还原
	SM_CXSCREEN_x=GetSystemMetrics(SM_CXSCREEN);
	SM_CYSCREEN_y=GetSystemMetrics(SM_CYSCREEN);
	//将屏幕设置成640*480
	DEVMODE DevMode;  //屏幕信息结构体
	EnumDisplaySettings(NULL,ENUM_CURRENT_SETTINGS,&DevMode);  //获取当前的数据
    DevMode.dmPelsWidth = 640;  //修改成你想要的分辨率
	DevMode.dmPelsHeight = 480;
	//DevMode.dmPelsWidth = 1920;  //修改成你想要的分辨率
	//DevMode.dmPelsHeight = 1080;
	ChangeDisplaySettings(&DevMode, CDS_UPDATEREGISTRY);  //设置生效

	return TRUE;
}



int CdistortionApp::ExitInstance()          //退出应用程序？
{
	//将屏幕设置还原
	DEVMODE DevMode;  //屏幕信息结构体
	EnumDisplaySettings(NULL,ENUM_CURRENT_SETTINGS,&DevMode);  //获取当前的数据
	DevMode.dmPelsWidth = SM_CXSCREEN_x;  //修改成你想要的分辨率
	DevMode.dmPelsHeight =SM_CYSCREEN_y;
	ChangeDisplaySettings(&DevMode, CDS_UPDATEREGISTRY);  //设置生效

	//释放变量
	cvReleaseCapture(&capture);
	cvReleaseImage(&inFrame);//recite in 7_7
	cvReleaseImage(&currentFrame);
	cvReleaseImage(&tempFrame);
	//cvReleaseImage(&tempFrame0);
	cvReleaseImage(&outFrame);
	cvReleaseImage(&phosphene);
	cvReleaseMat(&disLoc);
	cvReleaseMat(&neiLoc);
	return CWinApp::ExitInstance();
}

// CdistortionApp 消息处理程序




// 用于应用程序“关于”菜单项的 CAboutDlg 对话框

class CAboutDlg : public CDialog
{
public:
	CAboutDlg();

// 对话框数据
	enum { IDD = IDD_ABOUTBOX };

protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV 支持

// 实现
protected:
	DECLARE_MESSAGE_MAP()
};

CAboutDlg::CAboutDlg() : CDialog(CAboutDlg::IDD)
{
}

void CAboutDlg::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
}

BEGIN_MESSAGE_MAP(CAboutDlg, CDialog)
END_MESSAGE_MAP()

// 用于运行对话框的应用程序命令
void CdistortionApp::OnAppAbout()
{
	CAboutDlg aboutDlg;
	aboutDlg.DoModal();
}


// CdistortionApp 消息处理程序

