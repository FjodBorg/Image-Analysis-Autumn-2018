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

int64_t ticks = 0;


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


// Check om type casts er rigtige
void filter(IplImage *img, IplImage *filtered, int32_t filter[3][3]){
	uint8_t divider = 0;
	int32_t temp;
	uint32_t height = img->height; 
	uint32_t width = img->width; 
	uint32_t index = 0;
	int8_t value = 0;
	
	for (int k = 0; k < 3; k++) {
		for (int l = 0; l < 3; l++) {
			divider += filter[k][l];
			printf("%d ", filter[k][l]);
		}
	}
	if (divider == 0) {
		divider = 1;
	}
	printf("%d ", divider);

	for (uint32_t i = 1; height-1 > i; i++) {
		for (uint32_t j = 1; width-1 > j; j++) {
			temp = 0;
			for (int k = 0; k < 3; k++) {
				for (int l = 0; l < 3; l++) {
					index = j + (l - 1) + (i + (k - 1)) * width;
					value = img->imageData[index];
					temp += (int32_t)(value) * (filter[k][l]);
				
				}
			} 
			if (temp < 0)
				temp = 0;
			filtered->imageData[j + i * width]= (uint8_t)(temp/divider);
		}
	}
}

int main(int argc, char* argv[]) {
	int width = 480;
	int height = 480;
	int scale = 10;
	int32_t avg[3][3] = {
	{1,1,1},
	{1,1,1},
	{1,1,1}
	};
	int32_t lap[3][3] = {
	{0,1,0},
	{1,-4,1},
	{0,1,0}
	};

	// Names for plots and images
	const char* str = "rocket";
	const char* strCopy = "spaceMorten";
	const char* strTest = "paceMorten";
	const char* histStr = "Morteb";
	const char* histStrCopy = "spaceMorteb";

	unsigned char threshold = 90;
//	unsigned long int center[2];
//	unsigned long int moment[6];

	//Define image streams
	IplImage *img, *imgCopy, *test;

	img = cvLoadImage("C:\\Users\\30330\\Documents\\Exersizes\\PEN.png", 0);

	//init copy of image
	imgCopy = cvCreateImage(cvSize(img->width, img->height), IPL_DEPTH_8U, img->nChannels);
	test = cvCreateImage(cvSize(img->width, img->height), IPL_DEPTH_8U, img->nChannels);
	//copy img to imgCopy
	cvCopy(img, imgCopy, 0);


	hist(histStr, img, height, width, scale);

	ticks = cvGetTickCount();
	filter(img, imgCopy,avg);
	//filter(img, imgCopy, lap);
	printf("Tickcount %ld\n ", cvGetTickCount() - ticks);
	

	float K[] = {1,1,1, 1,1,1, 1,1,1};
	for (int i = 0; i < 9; i++) {
		K[i] = K[i]/9;
	}

	ticks = cvGetTickCount();
	CvMat kernel = cvMat(3, 3, CV_32FC1, K);

	cvFilter2D(img, test, &kernel);
	printf("Tickcount %ld\n ", cvGetTickCount() - ticks);

	cvShowImage(str, img);
	cvShowImage(strTest, test);
	cvShowImage(strCopy, imgCopy);
	while (true) { //Makes the window stay open
		if (cvWaitKey(5) > 0)
			break;
	}

	cvReleaseImage(&img);
	cvDestroyWindow(str);
	return 0;
}
