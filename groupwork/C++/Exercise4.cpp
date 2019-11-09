#include "Exercise4.h"



Exercise4::Exercise4(std::string path)
{
	originalPicture = cv::imread(path.c_str(), cv::IMREAD_GRAYSCALE); // Read the file
	filteredPicture = originalPicture.clone();
}

void Exercise4::gradient()
{
	float KernelX[3][3] = {{-1,0,1},
	                      {-1,0,1},
						  {-1,0,1}};


	cv::Mat data = cv::Mat(3, 3, CV_32F, &KernelX);
	cv::filter2D(originalPicture, filteredPicture, -1, data);
	cv::namedWindow("gradient in X", cv::WINDOW_AUTOSIZE); // Create a window for display.
	cv::imshow("gradient in X", filteredPicture);                // Show our image inside it
	
	float KernelY[3][3] = {{-1,-1,-1},
					       {0,0,0},
					       {1,1,1}};

	data = cv::Mat(3, 3, CV_32F, &KernelY);
	cv::filter2D(originalPicture, filteredPicture, -1, data);
	cv::namedWindow("gradient in y", cv::WINDOW_AUTOSIZE); // Create a window for display.
	cv::imshow("gradient in y", filteredPicture);                // Show our image inside it
	
	cv::namedWindow("original", cv::WINDOW_AUTOSIZE); // Create a window for display.
	cv::imshow("original", originalPicture);                // Show our image inside it
}

void Exercise4::laplacian()
{
	float Kernel[3][3] = {{0,1,0},
						   {1,-4,1},
						   {0,1,0}};

	cv::Mat data = cv::Mat(3, 3, CV_32F, &Kernel);
	cv::filter2D(originalPicture, filteredPicture, -1, data);
	cv::namedWindow("laplace", cv::WINDOW_AUTOSIZE); // Create a window for display.
	cv::imshow("laplace", filteredPicture);                // Show our image inside it

}

void Exercise4::contourSearch(double threshold) {
	cv::threshold(originalPicture, originalPicture, threshold,255, cv::THRESH_BINARY_INV);
	cv::Size size=cv::Size(originalPicture.cols, originalPicture.rows);
	filteredPicture = cv::Mat::zeros(size,CV_8U);

	int searchingDirection;

	std::vector<std::vector<int>> values = {{0,-1},{1,-1},{1,0},{1,1},{0,1},{-1,1},{-1,0},{-1,-1}};
	std::map<int, std::vector<int>> map;
	for (int i = 0; i < values.size(); i++)
		map.insert(std::pair<int, std::vector<int>>(i, values.at(i)));
	
	int direction = 2;
	int xIn, yIn;
	bool found=FALSE;

	for (int x = 0; x < originalPicture.cols; x++)
	{
		for (int y = 0; y < originalPicture.rows; y++)
		{
			if (originalPicture.at<uchar>(y, x) > threshold)
			{
				xIn = x;
				yIn = y;
				found = TRUE;
				break;
			}
		}
		if(found)
			break;
	}

	int xOut = INT16_MAX;
	int yOut = INT16_MAX;
	int xTmp = xIn;
	int yTmp = yIn;

	int nextDirection;
	int maxIterations = 20000;

	cv::namedWindow("blank", cv::WINDOW_AUTOSIZE); // Create a window for display.
	while (xOut != xIn && yOut != yIn)
	{
		
		filteredPicture.at<uchar>(yTmp, xTmp) = 255;
		nextDirection = (direction > 1) ? direction - 2 : direction + 6;
		std::rotate(values.begin(), values.begin() + nextDirection, values.end());
		for (int i = 0; i < values.size(); i++)
		{
			if (originalPicture.at<uchar>(yTmp + values.at(i)[1], xTmp + values.at(i)[0]) > 0)
			{
				xTmp += values.at(i)[0];
				yTmp += values.at(i)[1];
				direction = i;
				break;
			}
		}
		maxIterations--;
		if(!maxIterations)
			break;
	}
	cv::imshow("blank", filteredPicture);                // Show our image inside it

	


	int count = 0;
	const int MAX_RAND = 500;  //max iterations
	int pos = 0;
	//int B = 1;  // what does this do? y*width maybe?
	int B = width*y;
	int end = MAX_RAND;

	int rimx[MAX_RAND];
	int rimy[MAX_RAND];

	unsigned char *pic; // placeholder for image data 
	int newpos, oldpos, local_tresh, draw_type;
	pic = originalPicture.data;
	

	draw_type = 0;
	newpos = pos; // pos equals the starting position in the image ( =  y*Width + x)
	
	while (newpos >= 0L && newpos < end)
	{
		oldpos = newpos;
		rimx[count] = newpos % B; // save current position in list 
		rimy[count] = newpos / B;
		count++;
	
		// my fix
		draw_type = (draw_type +6) & 7; // AND with 0x00000111 (faster than modulus %8) // Select next search direction always -90deg from last move to ensure always the right contour
		for (int i = 0; i<8 && oldpos == newpos; i++){
			draw_type = (draw_type + i) & 7;
			switch (draw_type){
				case 0: if (pic[newpos + 1] > local_tresh) { 
					newpos += 1; break; }
				case 1: if (pic[newpos + B + 1] > local_tresh){
					newpos += B + 1; break;}
				case 2: if (pic[newpos + B] > local_tresh) { 
					newpos += B; break; }
				case 3: if (pic[newpos + B - 1] > local_tresh) {
					newpos += B -	1; break;}
				case 4: if (pic[newpos - 1] > local_tresh) { 
					newpos -= 1; break; }
				case 5: if (pic[newpos - B - 1] > local_tresh) {
					newpos -= B + 1; break;}
				case 6: if (pic[newpos - B] > local_tresh) { 
					newpos -= B; break; }
				case 7: if (pic[newpos - B + 1] > local_tresh) {
					newpos -= B -1; break;}
			}
		}
		/* doesn't work?
		draw_type = (draw_type + 6) % 8; // Select next search direction 
		switch (draw_type)
		{
		case 0: if (pic[newpos + 1] > local_tresh) { 
			newpos += 1; 
			draw_type = 0; 
			break; }
		case 1: if (pic[newpos + B + 1] > local_tresh){
			newpos += B + 1; 
			draw_type = 1; 
			break;}
		case 2: if (pic[newpos + B] > local_tresh) { 
			newpos += B; 
			draw_type = 2; 
			break; }
		case 3: if (pic[newpos + B - 1] > local_tresh) {
			newpos += B -	1; 
			draw_type = 3; 
			break;}
		case 4: if (pic[newpos - 1] > local_tresh) { 
			newpos -= 1; 
			draw_type = 4; 
			break; }
		case 5: if (pic[newpos - B - 1] > local_tresh) {
			newpos -= B + 1; 
			draw_type = 5; 
			break;}
		case 6: if (pic[newpos - B] > local_tresh) { 
			newpos -= B; 
			draw_type = 6; 
			break; }
		case 7: if (pic[newpos - B + 1] > local_tresh) {
			newpos -= B -1; 
			draw_type = 7; 
			break;}



		// case 8: if (pic[newpos + 1] > local_tresh) { 
		// 	newpos += 1; 
		// 	draw_type = 0; 
		// 	break; }
		// case 9: if (pic[newpos + B + 1] > local_tresh){
		// 	newpos += B + 1; 
		// 	draw_type = 1; 
		// 	break;}
		// case 10: if (pic[newpos + B] > local_tresh) { 
		// 	newpos += B; 
		// 	draw_type = 2; 
		// 	break; }
		// case 11: if (pic[newpos + B - 1] > local_tresh) {
		// 	newpos += B -1; 
		// 	draw_type = 3; 
		// 	break;}
		// case 12: if (pic[newpos - 1] > local_tresh) { 
		// 	newpos -= 1; 
		// 	draw_type = 4; 
		// 	break; }
		// case 13: if (pic[newpos - B - 1] > local_tresh) {
		// 	newpos -= B + 1; 
		// 	draw_type = 5; 
		// 	break;}
		// case 14: if (pic[newpos - B] > local_tresh) { 
		// 	newpos -= B; 
		// 	draw_type = 6; 
		// 	break; }
					*/
		}
		// If we are back at the beginning, we declare success 
		if (newpos == pos)
			break;
		// Abort if the contour is too complex. 
		if (count >= MAX_RAND)
			break;
	}
	
	// draw found line
}
