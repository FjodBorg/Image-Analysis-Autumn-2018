#include "cv.h" 
#include "cxcore.h" 
#include "highgui.h" 
#include <conio.h> 
#include <string>
#include <iostream>
#include <numeric>
#include <map>
#include <chrono>
#include <thread>
#pragma once
class Exercise4
{
public:
	Exercise4(std::string path);
	~Exercise4() {};
	void gradient();
	void laplacian();
	void contourSearch(double threshold);
	void setImage(cv::Mat image){ originalPicture = image.clone(); }
private:
	cv::Mat originalPicture;
	cv::Mat filteredPicture;
};
