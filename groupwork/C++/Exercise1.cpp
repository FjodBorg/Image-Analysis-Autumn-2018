#include "Exercise1.h"



Exercise1::Exercise1(std::string pathToFile)
{
	m_Image = cv::imread(pathToFile.c_str(), cv::IMREAD_GRAYSCALE); // Read the file
	if (m_Image.empty())                      // Check for invalid input
	{
		std::cout << "Could not open or find the image" << std::endl;
	}
}

void Exercise1::Display(std::string windowName)
{
	cv::namedWindow(windowName, cv::WINDOW_AUTOSIZE); // Create a window for display.
	cv::imshow(windowName, m_Image);                // Show our image inside it
}

void Exercise1::displayWithThreshold(std::string windowName,int t) {
	cv::namedWindow(windowName, cv::WINDOW_AUTOSIZE); // Create a window for display.
	cv::threshold(m_Image, m_Image, t, 255, cv::THRESH_BINARY);
	cv::imshow(windowName, m_Image);                // Show our image inside it
}



void Exercise1::Histogram()
{
	int max;
	int min;
	int k = 256;
	while (k-- > 0)
		m_Histogram[k] = 0; // reset histogram entry 

	for (int i = 0; i < m_Image.rows; i++)
		for (int j = 0; j < m_Image.cols; j++)
			m_Histogram[m_Image.at<uchar>(i, j)] += 1;

	min = *std::min_element(m_Histogram, m_Histogram + 256);
	max = *std::max_element(m_Histogram, m_Histogram + 256);

	for (int i = 0; i < 256; i++)
		m_Histogram[i] = 255 * ((m_Histogram[i] - min) / (max - min));

	cv::namedWindow("Histogram", cv::WINDOW_AUTOSIZE); // Create a window for display.
	cv::Mat histImg = cv::Mat::zeros(cv::Size(510, 255), CV_8U);

	for (int i = 0; i < 255; i++)
		line(histImg, cv::Point(i * 2, 255), cv::Point(i * 2, 255 - m_Histogram[i]), cv::Scalar(255, 255, 0), 2, 8, 0);
	imshow("Histogram", histImg);
}

void Exercise1::CameraCapture()
{
	if (!m_Cap.open(0))
		std::cout << "Cannot access camera device" << std::endl;
	for (;;)
	{
		m_Cap >> m_Image;
		if (m_Image.empty()) break; // end of video stream
		cv::imshow("this is you, smile! :)", m_Image);
		if (cv::waitKey(10) == 27) break; // stop capturing by pressing ESC 
	}
	cv::imwrite("../Debug/frame.jpg", m_Image);
}
