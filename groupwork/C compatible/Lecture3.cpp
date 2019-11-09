#include <stdio.h>
#include <conio.h>
#include "cv.h" 
#include "cxcore.h"
#include "highgui.h"
#include "math.h"

//pictures height and width have to be dividable by 32bit (4bytes or else it will read artifacts)

//Offset of histogram
#define XOFFSET 100
#define YOFFSET 450
#define PIXEL_MAX 255
#define PIXEL_MIN 0
#define PI 3.14159265358979323846


void hist(const char* histStr, IplImage *img, int height, int width, int scale){
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
	img = cvCreateImage(cvSize(width, height),IPL_DEPTH_8U,1);
	cvRectangle(img, cvPoint(0,0), cvPoint(width - 1, height - 1), cvScalar(0), CV_FILLED);

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

double getAngle(double moment[6]) {
	double mu11 = (moment[3] - moment[2] / moment[0] * moment[1]);
	double mu20 = (moment[4] - moment[1] / moment[0] * moment[1]);
	double mu02 = (moment[5] - moment[2] / moment[0] * moment[2]);
	printf("%lf %lf %lf\n", mu11, mu20, mu02);

	double mu11p = (double)(mu11 / moment[0]);
	double mu20p = (double)(mu20 / moment[0]);
	double mu02p = (double)(mu02 / moment[0]);
	return atan2((mu20p - mu02p),2 * mu11p) / 2 * 180 / PI;//degrees not needed
	//not working with 45° 90°
}


void printMoment(double moment[6]) {
	printf("M00 :%lf\nM10 :%lf\nM01 :%lf\nM11 :%lf\nM20 :%lf\nM02 :%lf\n", 
		moment[0], moment[1], moment[2] , moment[3], moment[4] , moment[5]);
}

double calcInvariantMoments(double mu, double mu00, double j, double i){
	double eta;
	double bot = pow(mu00, 1+(i+j)/2); // maybe not power
	eta = mu/bot;
	return eta;
}

void printInvariantMoment(double moment[6]){
	double mu00 = moment[0];
	double mu01 = 0;
	double mu10 = 0;
	double mu11 = (moment[3] - moment[2] / mu00 * moment[1]);
	double mu20 = (moment[4] - moment[1] / mu00 * moment[1]);
	double mu02 = (moment[5] - moment[2] / mu00 * moment[2]);

	double temp[7];
	temp[0] = calcInvariantMoments(mu00, mu00, 0, 0);
	temp[1] = calcInvariantMoments(mu10, mu00, 1, 0);
	temp[2] = calcInvariantMoments(mu01, mu00, 0, 1);
	temp[3] = calcInvariantMoments(mu11, mu00, 1, 1);
	temp[4] = calcInvariantMoments(mu20, mu00, 2, 0);
	temp[5] = calcInvariantMoments(mu02, mu00, 0, 2);
	temp[6] = temp[4]+temp[5];

	printf("eta00 :%lf\neta10 :%lf\neta01 :%lf\neta11 :%lf\neta20 :%lf\neta02 :%lf\nI1 :%lf\n", 
		temp[0], temp[1], temp[2] , temp[3], temp[4] , temp[5], temp[6]);
}



void getMoment(IplImage *binaryImg, double moment[6]) {
	unsigned long int M00 = 0; //area
	unsigned long int M01 = 0; //center of mass y
	unsigned long int M10 = 0; //center of mass x
	unsigned long int M11 = 0; //Correlation between M02 and M20
	unsigned long int M02 = 0; //distributions of mass in Y (moment of inertia)
	unsigned long int M20 = 0; //distributions of mass in X (Moment of inertia)
	unsigned long int x, y;

	unsigned int height = binaryImg->height;
	unsigned int width = binaryImg->width;

	for (unsigned int i = 0; i < height; i++) {
		for (unsigned int j = 0; j < width; j++) {
			if ((unsigned char)binaryImg->imageData[j + i * width] == PIXEL_MAX) {
				x = (j + 1);
				y = (i + 1);
				M00++;
				M10 += x;
				M01 += y;
				M11 += x * y; 
				M20 += x * x;
				M02 += y * y;
			}
		}
		printf("%d, %d, %d\n", M20, M02, M11);
	}
	moment[0] = M00;
	moment[1] = M10;
	moment[2] = M01;
	moment[3] = M11;
	moment[4] = M20;
	moment[5] = M02;
}

void centerOfMass(IplImage *binaryImg, unsigned long int center[2]) {
	unsigned long int x = 0;
	unsigned long int y = 0;
	unsigned long int area = 0;
	unsigned long int temp=0, temp2=0;

	uint32_t height = binaryImg->height;
	uint32_t width = binaryImg->width;

	for (uint32_t i = 0; i < height; i++) {
		for (uint32_t j = 0; j < width; j++) {
			if ((unsigned char)binaryImg->imageData[j+i*width] == PIXEL_MAX) {
				area++;
				x += j + 1;
				y += i + 1;
			}
		}
		if ((area-temp2))
			printf("%d\n", (x-temp)/(area-temp2));
		temp = x;
		temp2 = area;
		//printf("%d, %d, %d, %d,\n",width*height, width + i * width, width, height);
	}
	
	center[0] = x / area; 
	center[1] = y / area;
	printf("%d, %d\n", center[0], center[1]);
}





int main(int argc, char* argv[]){
	// Histogram scale and size
	int width = 480; 
	int height = 480;
	int scale = 10;

	// Names for plots and images
	const char* str = "rocket";
	const char* strCopy = "spaceMorten";
	const char* histStr = "Morteb";
	const char* histStrCopy = "spaceMorteb";

	unsigned char threshold = 90;
	unsigned long int center[2];
	double moment[6];

	//Define image streams
	IplImage *img, *imgCopy;

	img = cvLoadImage("C:\\Users\\30330\\Documents\\Exersizes\\test5.png", 0);

	//init copy of image
	imgCopy = cvCreateImage(cvSize(img->width, img->height), IPL_DEPTH_8U, img->nChannels);
	//copy img to imgCopy
	cvCopy(img, imgCopy, 0);

	hist(histStr, img, height, width, scale);


	//apply threshold to image copy
	thresholdImg(imgCopy, threshold);

	centerOfMass(imgCopy, center);
	getMoment(imgCopy, moment);
	printMoment(moment);
	printInvariantMoment(moment);


	cvShowImage(strCopy, imgCopy);

	// dot plz
	cvLine(img, cvPoint(center[0], center[1]), cvPoint(center[0], center[1]), cvScalar(150), 5);

	//show img on window "str"
	cvShowImage(str, img);

	printf("\n%lf",getAngle(moment));

	cvSaveImage("C:\\Users\\30330\\Documents\\Exersizes\\test6.png", imgCopy);

	while (true) { //Makes the window stay open
		if (cvWaitKey(5) > 0)
			break;
	}

	cvReleaseImage(&img);
	cvDestroyWindow(str);
	return 0;
}
