#include "cv.h" 
#include "cxcore.h" 
#include "highgui.h" 
#include <conio.h> 
#include <string>
#include <iostream>

#pragma once
class Exercise1
{
public:
	Exercise1(std::string);
	~Exercise1() {};
	void Display(std::string windowName);
	void CameraCapture();
	void Histogram();
	void setImage(cv::Mat image) { m_Image = image; }
	cv::Mat getImage() { return m_Image; }
	void displayWithThreshold(std::string windowName,int t);
private:
	cv::VideoCapture m_Cap;
	cv::Mat m_Image;
	double m_Histogram[256];
};

