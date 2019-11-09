#include "Exercise3.h"



Exercise3::Exercise3(std::string path)
{
	originalPicture = cv::imread(path.c_str(), cv::IMREAD_GRAYSCALE); // Read the file
	filteredPicture = originalPicture.clone();
}


void Exercise3::lowPassFilter()
{


	float Kernel[3][3] = {
					  {1 / 9.0, 1 / 9.0, 1 / 9.0},
					  {1 / 9.0, 1 / 9.0, 1 / 9.0},
					  {1 / 9.0, 1 / 9.0, 1 / 9.0}};

	cv::Mat data = cv::Mat(3, 3, CV_32F, &Kernel);
	cv::filter2D(originalPicture, filteredPicture,-1,data);

	cv::namedWindow("lowPass", cv::WINDOW_AUTOSIZE); // Create a window for display.
	cv::imshow("lowPass", filteredPicture);                // Show our image inside it
	cv::namedWindow("original", cv::WINDOW_AUTOSIZE); // Create a window for display.
	cv::imshow("original",originalPicture);                // Show our image inside it

}

void Exercise3::highPassFilter()
{

	float Kernel[3][3] = {
		{ 0, 1 , 0 },
		{ 1, -4, 1 },
		{ 0, 1, 0 } };

	cv::Mat data = cv::Mat(3, 3, CV_32F, &Kernel);
	cv::filter2D(originalPicture, filteredPicture, -1, data);

	cv::namedWindow("highPass", cv::WINDOW_AUTOSIZE); // Create a window for display.
	cv::imshow("highPass", filteredPicture);                // Show our image inside it
	cv::namedWindow("original", cv::WINDOW_AUTOSIZE); // Create a window for display.
	cv::imshow("original", originalPicture);                // Show our image inside it

}

void Exercise3::laplaceGaussianFilter()
{

	float Kernel[9][9] = {
		{ 0, 0, 1, 2, 2, 2, 1, 0, 0 },
		{ 0, 1, 5, 10, 12, 10, 5, 1, 0 },
		{ 1, 5, 15, 19, 16, 19, 15, 5, 1 },
		{ 2, 10, 19, -19,-64,-19,19,10,2 },
		{ 2, 12, 16, -64, -148,-64,16,12,2 },
		{ 2, 10, 19,-19,-64,-19,19,10,2 },
		{ 1, 5, 15,19,16,19,15,5,1 },
		{ 0,1,5,10,12,10,5,1,0 },
		{ 0,0,1,2,2,2,1,0,0 },
	};

	cv::GaussianBlur(originalPicture, originalPicture, cv::Size(3, 3), 0, 0, cv::BORDER_DEFAULT);

	cv::Mat data = cv::Mat(9, 9, CV_32F, &Kernel);
	cv::filter2D(originalPicture, filteredPicture, -1, data);

	cv::namedWindow("LAPLACE", cv::WINDOW_AUTOSIZE); // Create a window for display.
	cv::imshow("LAPLACE", filteredPicture);                // Show our image inside it
	cv::namedWindow("original", cv::WINDOW_AUTOSIZE); // Create a window for display.
	cv::imshow("original", originalPicture);                // Show our image inside it

}

void Exercise3::fractileFilter(double percentage)
{
	std::array<uchar, 9>window;
	int tmp = ceil(9*percentage);
	for (int i = 1; i < originalPicture.rows-1; i++)
	{
		for (int j = 1; j < originalPicture.cols - 1; j++)
		{
			for (int y = -1; y < 2; y++)
			{
				for (int x = -1; x < 2; x++)
				{
					window.at((y+1)*3+x+1)= originalPicture.at<uchar>(i + y, j + x);
				}
			}
			std::sort(window.begin(), window.end());
			filteredPicture.at<uchar>(i, j) = window.at(tmp);
		}
	}
	cv::namedWindow("fractal", cv::WINDOW_AUTOSIZE); // Create a window for display.
	cv::imshow("fractal", filteredPicture);                // Show our image inside 
}
