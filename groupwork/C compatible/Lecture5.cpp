#include <stdio.h>
#include <conio.h>
#include "cv.h" 
#include "cxcore.h"
#include "highgui.h"
#include "math.h"

//Offset of histogram
#define XOFFSET 100
#define YOFFSET 450
#define PIXEL_MAX 255
#define PIXEL_MIN 0
#define PI 3.14159265358979323846


void hist(const char* histStr, IplImage *img, int height, int width, int scale) {
	unsigned char *channel; // pre-allocated array holding image data for the color channel or the grayscale image.
	unsigned char value = 0;
	int k = 256;
	int histogram[256]; // histogram array - remember to set to zero initially

	//resize window
	cvNamedWindow(histStr);
	cvMoveWindow(histStr, 900, 100);
	cvResizeWindow(histStr, width, height);

	//extract image data
	channel = (unsigned char*)img->imageData;

	while (k-- > 0)
		histogram[k] = 0; // reset histogram entry
	for (int i = 0; i < (img->width)*(img->height); i++)
	{
		value = channel[i];
		histogram[value] += 1;
	}

	//create new image for window
	img = cvCreateImage(cvSize(width, height), IPL_DEPTH_8U, 1);
	cvRectangle(img, cvPoint(0, 0), cvPoint(width - 1, height - 1), cvScalar(0), CV_FILLED);

	//Draw color lines below histogram
	for (int i = 0; i < 256; i += 64) {
		cvLine(img, cvPoint(XOFFSET + i, YOFFSET + 10), cvPoint(XOFFSET + i + 64, YOFFSET + 10), cvScalar(i + 20), 3);
	}

	//Draw histogram
	for (int i = 0; i < sizeof(histogram) / sizeof(int); i++) {
		cvLine(img, cvPoint(XOFFSET + i, YOFFSET), cvPoint(XOFFSET + i, YOFFSET - histogram[i] / scale), cvScalar(255), 0);
	}
	cvShowImage(histStr, img);

}

void thresholdImg(IplImage *img, unsigned char threshold) {
	for (uint32_t i = 0; i < (uint32_t)(img->height) * (img->width); i++) {
		//printf("%d, %d, ", i%9 ,(unsigned char) img->imageData[i]);
		if ((unsigned char)img->imageData[i] < threshold)
			img->imageData[i] = (char)PIXEL_MAX;
		else
			img->imageData[i] = (char)PIXEL_MIN;
		//printf("%d\n", (unsigned char)img->imageData[i]);
	}
	printf("%d,\n", (img->height) * (img->width));
}





int main(int argc, char* argv[]) {
	int width = 480;
	int height = 480;
	int scale = 10;
	uint32_t avg[3][3] = {
	{1,1,1},
	{1,1,1},
	{1,1,1}
	};
	uint32_t lap[3][3] = {
	{0,1,0},
	{1,-4,1},
	{0,1,0}
	};

	// Names for plots and images
	const char* str = "1";
	const char* strCopy = "2";
	const char* strTest = "gradiant";
	const char* histStr = "22";
	const char* histStrCopy = "gradiant2";

	unsigned char threshold = 90;
	//	unsigned long int center[2];
	//	unsigned long int moment[6];

		//Define image streams
	IplImage *img, *imgLP, *imgLAP, *imgGY, *imgGX;

	img = cvLoadImage("C:\\Users\\30330\\Documents\\Exersizes\\grand-canyon-colorado-river.jpg", 0);

	//init copy of image
	imgLP = cvCreateImage(cvSize(img->width, img->height), IPL_DEPTH_8U, img->nChannels);
	imgLAP = cvCreateImage(cvSize(img->width, img->height), IPL_DEPTH_8U, img->nChannels);
	imgGY = cvCreateImage(cvSize(img->width, img->height), IPL_DEPTH_8U, img->nChannels);
	imgGX = cvCreateImage(cvSize(img->width, img->height), IPL_DEPTH_8U, img->nChannels);
	//copy img to imgCopy
	cvCopy(img, imgLP, 0);


	
	//float K1[] = { 1,1,1,1,1, 1,1,1,1};
	//float K1[] = { 1,2,1, 2,4,2, 1,2,1 };
	float K1[] = {	2,	4,	5,	9,	5,	4,	2,
					4,	9,	12, 15,	12,	9,	4, 
					5,	12,	15,	18,	15,	12,	5,
					9,	15,	18,	25,	18,	15,	12,
					5,	12,	15,	18,	15,	12,	5,
					4,	9,	12, 15,	12,	9,	4,
					2,	4,	5,	9,	5,	4,	2,};
	int32_t K2[] = { 0,1,0, 1,-4,1, 0,1,0 };
	//int32_t K2[] = { 1,1,1, 1,-8,1, 1,1,1 };
	int32_t Ky[] = { 1,0,-1, 2,0,-2, 1,0,-1 };
	int32_t Kx[] = { 1,2,1, 0,0,0, -1,-2,-1 };
	int ksize = sizeof(K1) / sizeof(int32_t);
	int dsize = int(sqrt(ksize));
	int temp = 0;
	for (int i = 0; i < ksize; i++) {
		temp  += K1[i];
	}
	for (int i = 0; i < ksize; i++) {
		K1[i] /= temp;
	}

	CvMat kernel1 = cvMat(dsize, dsize, CV_32FC1, K1);
	CvMat kernel2 = cvMat(3, 3, CV_32SC1, K2);
	CvMat kernel3 = cvMat(3, 3, CV_32SC1, Ky);
	CvMat kernel4 = cvMat(3, 3, CV_32SC1, Kx);

	cvFilter2D(img, imgLP, &kernel1);
	cvFilter2D(imgLP, imgLAP, &kernel2);
	cvFilter2D(imgLP, imgGY, &kernel3);
	cvFilter2D(imgLP, imgGX, &kernel4);

	for (int i = 0; i < imgGY->height *imgGY->width; i++) {
		// test if it is working 
		imgGX->imageData[i] = sqrt(pow(imgGY->imageData[i],2) + pow(imgGX->imageData[i],2))/sqrt(2);
		// old piece of code
		//imgGX->imageData[i] = (imgGY->imageData[i] + imgGX->imageData[i])/2;
	}

	//for (int i = 0; i < imgGY->height *imgGY->width; i++) {
	//	imgLAP->imageData[i] = (imgLAP->imageData[i] + imgGX->imageData[i] / 2);
	//}

	printf("%d, %d\n",imgGY->height,img->height);
	//cvFilter2D(imgCopy, imgCopy, &kernel4);
	//cvFilter2D(img, test, &kernel2);
	//cvFilter2D(imgCopy, test, &kernel3);

	hist(histStr, imgLAP, height, width, scale);
	hist(histStrCopy, imgGX, height, width, scale);
	//hist(histStrCopy, imgLP, height, width, scale);

	thresholdImg(imgLAP, 5);
	thresholdImg(imgGX, 40);
	
	//thresholdImg(imgLP, 120);

	cvShowImage(str, img);
	cvShowImage(strCopy, imgLAP);
	cvShowImage(strTest, imgGX);
	//cvShowImage(strTest, imgLP);

	while (true) { //Makes the window stay open
		if (cvWaitKey(5) > 0)
			break;
	}

	cvReleaseImage(&img);
	cvDestroyWindow(str);
	return 0;
}