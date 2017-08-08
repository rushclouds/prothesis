// distortion.cpp : ����Ӧ�ó��������Ϊ��
//

#include "stdafx.h"
#include "distortion.h"
#include "MainFrm.h"


#ifdef _DEBUG
#define new DEBUG_NEW
#endif

// CdistortionApp

BEGIN_MESSAGE_MAP(CdistortionApp, CWinApp)                   //ʹ��BEGIN_MESSAGE_MAP�꿪ʼ��Ϣӳ�䣬
	ON_COMMAND(ID_APP_ABOUT, &CdistortionApp::OnAppAbout)    //Ȼ��Ϊÿ����Ϣ����������һ����ڣ�
END_MESSAGE_MAP()                                            //�����END_MESSAGE_MAP�������Ϣӳ��
                                                             //���Ǻ궨�壬���Ǻ���

// CdistortionApp ����

CdistortionApp::CdistortionApp()
{
	// TODO: �ڴ˴���ӹ�����룬
	// ��������Ҫ�ĳ�ʼ�������� InitInstance ��
}


// Ψһ��һ�� CdistortionApp ����

CdistortionApp theApp;


// CdistortionApp ��ʼ��

BOOL CdistortionApp::InitInstance()
{
	// ���һ�������� Windows XP �ϵ�Ӧ�ó����嵥ָ��Ҫ
	// ʹ�� ComCtl32.dll �汾 6 ����߰汾�����ÿ��ӻ���ʽ��
	//����Ҫ InitCommonControlsEx()�����򣬽��޷��������ڡ�
	INITCOMMONCONTROLSEX InitCtrls;
	InitCtrls.dwSize = sizeof(InitCtrls);
	// ��������Ϊ��������Ҫ��Ӧ�ó�����ʹ�õ�
	// �����ؼ��ࡣ
	InitCtrls.dwICC = ICC_WIN95_CLASSES;
	InitCommonControlsEx(&InitCtrls);

	CWinApp::InitInstance();

	// ��׼��ʼ��
	// ���δʹ����Щ���ܲ�ϣ����С
	// ���տ�ִ���ļ��Ĵ�С����Ӧ�Ƴ�����
	// ����Ҫ���ض���ʼ������
	// �������ڴ洢���õ�ע�����
	// TODO: Ӧ�ʵ��޸ĸ��ַ�����
	// �����޸�Ϊ��˾����֯��
	SetRegistryKey(_T("Ӧ�ó��������ɵı���Ӧ�ó���"));
	// ��Ҫ���������ڣ��˴��뽫�����µĿ�ܴ���
	// ����Ȼ��������ΪӦ�ó���������ڶ���
	CMainFrame* pFrame = new CMainFrame;
	if (!pFrame)
		return FALSE;
	m_pMainWnd = pFrame;
	// ���������ؿ�ܼ�����Դ
	pFrame->LoadFrame(IDR_MAINFRAME,
		WS_OVERLAPPEDWINDOW | FWS_ADDTOTITLE, NULL,
		NULL);






	// Ψһ��һ�������ѳ�ʼ���������ʾ����������и���
	pFrame->ShowWindow(SW_SHOW);
	pFrame->UpdateWindow();
	// �������к�׺ʱ�ŵ��� DragAcceptFiles
	//  �� SDI Ӧ�ó����У���Ӧ�� ProcessShellCommand ֮����

	/**************���÷ֱ���**************/
	//����ԭʼ��Ļ�ֱ��ʣ������˳���ԭ
	SM_CXSCREEN_x=GetSystemMetrics(SM_CXSCREEN);
	SM_CYSCREEN_y=GetSystemMetrics(SM_CYSCREEN);
	//����Ļ���ó�640*480
	DEVMODE DevMode;  //��Ļ��Ϣ�ṹ��
	EnumDisplaySettings(NULL,ENUM_CURRENT_SETTINGS,&DevMode);  //��ȡ��ǰ������
    DevMode.dmPelsWidth = 640;  //�޸ĳ�����Ҫ�ķֱ���
	DevMode.dmPelsHeight = 480;
	//DevMode.dmPelsWidth = 1920;  //�޸ĳ�����Ҫ�ķֱ���
	//DevMode.dmPelsHeight = 1080;
	ChangeDisplaySettings(&DevMode, CDS_UPDATEREGISTRY);  //������Ч

	return TRUE;
}



int CdistortionApp::ExitInstance()          //�˳�Ӧ�ó���
{
	//����Ļ���û�ԭ
	DEVMODE DevMode;  //��Ļ��Ϣ�ṹ��
	EnumDisplaySettings(NULL,ENUM_CURRENT_SETTINGS,&DevMode);  //��ȡ��ǰ������
	DevMode.dmPelsWidth = SM_CXSCREEN_x;  //�޸ĳ�����Ҫ�ķֱ���
	DevMode.dmPelsHeight =SM_CYSCREEN_y;
	ChangeDisplaySettings(&DevMode, CDS_UPDATEREGISTRY);  //������Ч

	//�ͷű���
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

// CdistortionApp ��Ϣ�������




// ����Ӧ�ó��򡰹��ڡ��˵���� CAboutDlg �Ի���

class CAboutDlg : public CDialog
{
public:
	CAboutDlg();

// �Ի�������
	enum { IDD = IDD_ABOUTBOX };

protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV ֧��

// ʵ��
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

// �������жԻ����Ӧ�ó�������
void CdistortionApp::OnAppAbout()
{
	CAboutDlg aboutDlg;
	aboutDlg.DoModal();
}


// CdistortionApp ��Ϣ�������

