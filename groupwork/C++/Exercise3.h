#include "cv.h" 
#include "cxcore.h" 
#include "highgui.h" 
#include <conio.h> 
#include <string>
#include <iostream>
#include <numeric>
#include <array>
#include <functional>
#include <algorithm>
#include <math.h>

#pragma once
class Exercise3
{
public:
	Exercise3(std::string path);
	~Exercise3() {};
	void lowPassFilter();
	void highPassFilter();
	void fractileFilter(double percentage);
	void laplaceGaussianFilter();
private:
	cv::Mat originalPicture;
	cv::Mat filteredPicture;
};

