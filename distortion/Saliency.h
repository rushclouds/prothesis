#pragma once

struct Saliency

{
	// Get saliency values of a group of images.
	// Input image names and directory name for saving saliency maps.
	static void Get(const string &imgNameW, const string &salDir);
	// Histogram Contrast of [3]
	static Mat GetHC(const Mat &img3f);
	// Color quantization
	static int Quantize(const Mat& img3f, Mat &idx1i, Mat &_color3f, Mat &_colorNum, double ratio = 0.95);
//private:
	static const int SAL_TYPE_NUM = 5;
	typedef Mat (*GET_SAL_FUNC)(const Mat &);
	static const char* SAL_TYPE_DES[SAL_TYPE_NUM];
	static const GET_SAL_FUNC gFuns[SAL_TYPE_NUM];

	// Histogram based Contrast
	static void GetHC(const Mat &binColor3f, const Mat &weights1f, Mat &colorSaliency);
	static void SmoothSaliency(const Mat &binColor3f, Mat &sal1d, float delta, const vector<vector<CostfIdx>> &similar);

	
};

/************************************************************************/
/*[1]R. Achanta, S. Hemami, F. Estrada and S. Susstrunk, Frequency-tuned*/
/*   Salient Region Detection, IEEE International Conference on Computer*/
/*	 Vision and Pattern Recognition (CVPR), 2009.						*/
/*[2]Y. Zhai and M. Shah. Visual attention detection in video sequences */
/*   using spatiotemporal cues. In ACM Multimedia, pages 815¨C824. ACM, */
/*   2006.																*/
/*[3]Our paper															*/
/*[4]X. Hou and L. Zhang. Saliency detection: A spectral residual		*/
/*   approach. In IEEE Conference on Computer Vision and Pattern		*/
/*	 Recognition, 2007. CVPR¡¯07, pages 1¨C8, 2007.						*/
/************************************************************************/
