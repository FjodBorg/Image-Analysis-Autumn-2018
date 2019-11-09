#include <stdio.h>
#include <conio.h>
#include "cxcore.h"
#include "highgui.h"

#define XOFFSET 100
#define YOFFSET 450

int main2(int argc, char* argv[]) {
	IplImage *img, *img_copy;
	const char* str = "rocket";
	const char* str2 = "Morteb";
	int hist[256]; 
	//cvInvert

	unsigned char *channel; // pre-allocated array holding image data for the color channel or the grayscale image.
	unsigned char value = 0; // index value for the histogram (not really needed)
	int histogram[256]; // histogram array - remember to set to zero initially
	int width = 480; // say, 320
	int height = 480; // say, 240 
	int k = 256;
	int scale = 5;

	//frame = cvLoadImage("C:\\Users\\30330\\Documents\\Exersizes\\ariane5_1b.jpg", 0);
	img = cvLoadImage("C:\\Users\\30330\\Desktop\\ariane5_1b.jpg", 1);
	cvNamedWindow(str);
	cvMoveWindow(str, 100, 100);
	//cvResizeWindow(str, width, height);
	
	printf("\nPreeeeeeeeeeeeeeess any key\n");

	cvNamedWindow(str2);
	cvMoveWindow(str2, 900, 100);
	cvResizeWindow(str2, width, height);
	
	channel = (unsigned char*)img->imageData;
	
	while (k-- > 0)
		histogram[k] = 0; // reset histogram entry
	for (int i = 0; i < (img->width)*(img->height); i++)
	{
		value = channel[i];
		histogram[value] += 1;
	}
	for (int i = 0; i < 256; i++) {
		printf("%d:  %d, %d, %d\n", (img->imageSize), histogram[i], img->height, img->width);
	}

	cvShowImage(str, img);
	
	
	img = cvCreateImage(cvSize(width, height),IPL_DEPTH_8U,1);
	cvRectangle(img, cvPoint(0,0), cvPoint(width-1, height-1), cvScalar(0), CV_FILLED);

	for (int i = 0; i < 256; i+=64) {
		cvLine(img, cvPoint(XOFFSET + i, YOFFSET+10), cvPoint(XOFFSET + i+64, YOFFSET + 10), cvScalar(i+20), 3);
	}

	for (int i = 0; i < sizeof(histogram)/sizeof(int); i++) {
		cvLine(img, cvPoint(XOFFSET+i, YOFFSET), cvPoint(XOFFSET+i, YOFFSET-histogram[i]/scale), cvScalar(255), 0);
	}
	
	
	cvShowImage(str2, img);
	
	//cvLine(img);
	


	
	while (true) {
		if (cvWaitKey(5) > 0)
			break;
	} 


	cvReleaseImage(&img);
	cvDestroyWindow(str);
	return 0;
}