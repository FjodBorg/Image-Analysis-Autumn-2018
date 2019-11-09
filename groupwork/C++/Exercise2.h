#include "cv.h" 
#include "cxcore.h" 
#include "highgui.h" 
#include <conio.h> 
#include <string>
#include <iostream>
#include <numeric>

#pragma once
class Exercise2
{
public:
	Exercise2();
	~Exercise2() {};
	void setImage(cv::Mat t_image) { m_Image = t_image.clone(); }
	void computeMoments();
	void displayCenterOfMass(std::string windowName);

	double m_M00;
	double m_M01;
	double m_M10;
	double m_M11;
	double m_M02;
	double m_M20;
	double mu00, mu11 ,mu20 ,mu02;

private:
	cv::Mat m_Image;
	std::vector<int> m_Vertical;
	std::vector<int> m_Horizontal;
	cv::Moments m_Moments;
};



