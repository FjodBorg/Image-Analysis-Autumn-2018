#include <stdio.h> 
#include "Exercise1.h" 
#include "Exercise2.h"
#include "Exercise3.h"
#include "Exercise4.h"
#include <conio.h> 
#include <string>
#include <iostream>

using namespace cv;
using namespace std;


int main(int argc, char* argv[])
{
	
	string imageName("../Debug/Exersices/ariane5_1b.jpg"); // by default
	/*
	Exercise1 exe1(imageName);
	exe1.Display("Original");
	exe1.Histogram();
	exe1.CameraCapture();
	Mat copied_image;
	copied_image=exe1.getImage().clone();
	exe1.setImage(copied_image);
	exe1.Display("Copy");
	exe1.displayWithThreshold("AfterThresholding",90);
	
	Exercise2 exe2;
	exe2.setImage(exe1.getImage());
	exe2.computeMoments();
	exe2.displayCenterOfMass("Center of mass");

	

	Exercise3 exe3(imageName);
	exe3.lowPassFilter();
	exe3.highPassFilter();
	exe3.laplaceGaussianFilter();
	exe3.fractileFilter(0.8);
	*/

	Exercise4 exe4(imageName);
	exe4.gradient();
	exe4.laplacian();
	
	//cv::Mat nextPicture = cv::imread("../Debug/Exersices/PEN.png", IMREAD_GRAYSCALE);
	//exe4.setImage(nextPicture);
	//exe4.contourSearch(125);

	waitKey(0); // Wait for a keystroke in the window	
	return 0;
}