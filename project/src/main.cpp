#include "opencv2/core.hpp"
#include "opencv2/highgui.hpp"
#include "opencv2/imgproc.hpp"

#include "preprocess.hpp"

#include <stdio.h>
#include <math.h>

int main(int argc, char *argv[]) {
	cv::Mat frame, frameLP, frameGauss, filteredX, filteredY, filteredLap, frameOG, frameGradMag, frameGradAng, frameOut, frameCircle;
	Circle foundCircle;
	Objects* map = NULL;
	uint32_t numObjects;
	float imageScale = 0.5f;
	const char* image = "original";
	const char* imageCopy = "copy";
	const char* imageX = "filtered x";
	const char* imageY = "filtered y";
	const char* imageLap = "laplace";
	const char* imageHist = "hist";

	int width = 380; // say, 320
	int height = 520; // say, 240
	int scale = 10;
	uint8_t threshold = 90;
	uint16_t center[2] = { 0, 0 };
	uint32_t moment[6] = { 0 };

	int32_t lowpass[3][3] = {
		{1,1,1},
		{1,1,1},
		{1,1,1}
	};
	int32_t laplace[3][3] = {
		{0,1,0},
		{1,-4,1},
		{0,1,0}
	};
	float divisorLP = 0;
	float divisorGrad = 255;
	float divisorGauss = 159;
	float kLP[] = { 1,1,1,1,1,1,1,1,1 };
	float kX[] = {1,0,-1,2,0,-2,1,0,-1};
	float kY[] = {1,2,1,0,0,0,-1,-2,-1};
	float kLap[] = { 0,1,0,1,-4,1,0,1,0 };
	float kGauss[] = { 2,4,5,4,2 , 4,9,12,9,4 , 5,12,15,12,5 , 4,9,12,9,4 , 2,4,5,4,2 };
	
//	myfile.open("../data/otsuCalc.txt");
//	myfile << "test\n";
//	myfile.close();
	
	for (uint8_t i = 0; i < 9; i++) {
		divisorLP += kLP[i];
	}
	
	for (uint8_t i = 0; i < 9; i++) {
		kLP[i] /= divisorLP;
		//kX[i] *= divisorGrad;
		//kY[i] *= divisorGrad;
		//kLap[i] *= divisorGrad;
	}
	for (uint8_t i = 0; i < 25; i++) {
		kGauss[i] /= divisorGauss;
	}
	
	cv::Mat kernelLP = cv::Mat(3, 3, CV_32FC1, kLP);
	cv::Mat kernelGauss = cv::Mat(5, 5, CV_32FC1, kGauss);
	cv::Mat kernelX = cv::Mat(3, 3, CV_32FC1, kX);
	cv::Mat kernelY = cv::Mat(3, 3, CV_32FC1, kY);
	cv::Mat kernelLap = cv::Mat(3, 3, CV_32FC1, kLap);

//	frame = cv::imread("../fig/envelope/cropped/long2.png",0);
	frameOG = cv::imread("../fig/envelope/cropped/long1.png", 0); // (argv[1], 0);  //ariane5_1b.jpg
	cv::resize(frameOG,frame,cv::Size(),imageScale,imageScale);
	
	frameOut = cv::Mat(cv::Size(frame.cols,frame.rows), CV_8UC1);
	frameCircle = cv::Mat(cv::Size(frame.cols,frame.rows), CV_8UC3);
	cv:cvtColor(frame, frameCircle, CV_GRAY2RGB);

	//findEllipse(frameLP, 50, 150, 10);
	cv::medianBlur(frame,frame,11);

//	cv::GaussianBlur(frame, frameGauss, cv::Size(13,13),0,0);
//	cv::Canny(frameGauss,frameGauss,30,60,3,1);
//	uint8_t thresh = 80;
//	cannyDetector(frame,frameOut,7,thresh/2,thresh);
	cannyDetector(frame,frameOut,11,0.03f);
	cv::imshow("test",frameOut);
	//morphologyEx(frameGradMag,frameOut,MORPH_OPEN,getStructuringElement(MORPH_RECT,Size(3,3)));

	numObjects = labeling(frameOut,&map,10,0);
	thresholdImage(frameOut,120);

//	labeling(frameOut, 15);

//	cv::imshow("before",frameOut);

//	cv::dilate(frameOut,frameOut, cv::getStructuringElement(0, cv::Size(3,3)));

//	cv::imshow("after erosion",frameOut);

	foundCircle = ransac(frameOut,map,numObjects,150*imageScale,900*imageScale,0.5f); //150 to 900
//	ransacVoting(frameOut,map,numObjects,75,450,0.5f);
	cv::cvtColor(frameOut,frameOut,CV_GRAY2RGB);
	cv::circle(frameOut, cv::Point(foundCircle.x,foundCircle.y), foundCircle.radius, cv::Scalar(255,0,0), 1);
	cv::circle(frameCircle, cv::Point(foundCircle.x,foundCircle.y), foundCircle.radius, cv::Scalar(255,0,0), 2);
	
	//histogram("gradmag hist",frame,height,width,scale*5);
	//histogram("geadmag 2",frame,height,width,0);
	//histogram("canny",frameOut,height,width,scale);

//	histogram("median threshold",frameGradMag,520,380,medT);
//	histogram("mean threshold",frameGradMag,520,380,meanT);
//	histogram("procent threshold",frameGradMag,520,380,procT);

	cv::imshow("out",frameOut);
	cv::moveWindow("out",400,0);
//	cv::imwrite("../fig/circleFound/9gauss_2perc_50_200_50perc/flash2Circles_edge.png", frameOut);
//	cv::imshow("opencv",frameGauss);
	cv::imshow("circle",frameCircle);
	cv::moveWindow("circle", 1000,0);
//	cv::imwrite("../fig/circleFound/9gauss_2perc_50_200_50perc/flash2Circles_circle.png", frameCircle);
	
//	cv::imwrite("../fig/test/otsu/gradMagOtsu.png",frameOut);
//	cv::imwrite("../fig/filtered/flashM038_gradMagCircle_9_30.png",frameCircle);
	
	for (;;) {
		if (cv::waitKey(5) > 0) {
			break;
		}
	}
	
	return 0;
}
