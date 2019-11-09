#include "Exercise2.h"



Exercise2::Exercise2()
{
}

void Exercise2::computeMoments()
{
	m_Horizontal.resize(m_Image.cols,0);
	m_Vertical.resize(m_Image.rows,0);

	for (int y = 0; y < m_Image.rows; y++)
	{
		for (int x = 0; x < m_Image.cols; x++)
		{
			m_M11 += x * y*m_Image.at<uchar>(y, x);
			if (!m_Image.at<uchar>(y, x))
			{
				m_Vertical.at(y) += 1;
				m_Horizontal.at(x) += 1;
			}
		}
	}

	for (int i = 0; i < m_Horizontal.size(); i++)
	{ 
		m_M00 += m_Horizontal.at(i);
		m_M10 += i * m_Horizontal.at(i);
		m_M20 += i * i * m_Horizontal.at(i);
	}

	for (int i = 0; i < m_Vertical.size(); i++)
	{
		m_M01 += i * m_Vertical.at(i);
		m_M02 += i * i * m_Vertical.at(i);
	}

	int X = m_M10 / m_M00;
	int Y = m_M01 / m_M00;

	mu00 = m_M00;
	mu11 = m_M11 - X * m_M01;
	mu20 = m_M20 - X * m_M10;
	mu02 = m_M02 - Y * m_M01;

}

void Exercise2::displayCenterOfMass(std::string windowName)
{
	cv::namedWindow(windowName, cv::WINDOW_AUTOSIZE); // Create a window for display.
	int X = m_M10/ m_M00;
	int Y = m_M01/ m_M00;

	int X1 = m_Moments.m10 / m_Moments.m00;
	int Y1=  m_Moments.m01 / m_Moments.m00;
	double angle = tan((2 * mu11) / (mu20 - mu02));
	double length = 60;

	cv::line(m_Image, cv::Point(X,Y), cv::Point(X + length * cos(angle), Y + length * sin(angle)), cv::Scalar(133, 0, 0));
	cv::drawMarker(m_Image, cv::Point(X, Y), cv::Scalar(133,0,0));

	cv::imshow(windowName, m_Image);                // Show our image inside it	
}