#include "include/preprocess.hpp"

cv::Mat histogram(const char* histStr, cv::Mat img, int height, int width, int scale) {
	//unsigned char* channel; // pre-allocated array holding image data for the color channel or the grayscale image.
	unsigned char value = 0;
	uint16_t k = 256;
	uint16_t histogram[256]; // histogram array - remember to set to zero initially
	uint16_t max = 0;
	uint32_t size = img.cols*img.rows;
//	cv::namedWindow(histStr, CV_WINDOW_NORMAL);
//	cv::resizeWindow(histStr, width, height);

	//channel = img.data;

	while (k-- > 0)
		histogram[k] = 0; // reset histogram entry
	for (int i = 0; i < (img.cols*img.rows); i++)
	{
		value = img.data[i];
		histogram[value] += 1;
		if(histogram[value]>max){
			max = histogram[value];
		}
	}

	img = cv::Mat(cv::Size(width,height), CV_8UC1);
	cv::rectangle(img, cv::Point(0, 0), cv::Point(width - 1, height - 1), cv::Scalar(0), CV_FILLED);

	for (int i = 0; i < 256; i += 64) {
		cv::line(img, cv::Point(XOFFSET + i, YOFFSET + 10), cv::Point(XOFFSET + i + 64, YOFFSET + 10), cv::Scalar(i + 20), 3);
	}
	
	if(scale==-2){
		max = 20*log10((float)max);
		for (uint32_t i = 0; i < 256; i++) {
			if(histogram[i]>0){
				cv::line(img, cv::Point(XOFFSET+i, YOFFSET), cv::Point(XOFFSET+i, (uint16_t)(YOFFSET-20*log10((float)histogram[i])/max*NORMALIZED_SCALE)), cv::Scalar(255), 1);		
			}
		}
	}else if(scale==-1){ // normalized by max(histogram)
		for (uint32_t i = 0; i < 256; i++) {
			cv::line(img, cv::Point(XOFFSET+i, YOFFSET), cv::Point(XOFFSET+i, (uint16_t)(YOFFSET-(float)histogram[i]/max*NORMALIZED_SCALE)), cv::Scalar(255), 1);
		}
	}else if(scale == 0){ // normalized by length(img)
		for (uint32_t i = 0; i < 256; i++) {
			cv::line(img, cv::Point(XOFFSET+i, YOFFSET), cv::Point(XOFFSET+i, (uint16_t)(YOFFSET-(float)histogram[i]/size*NORMALIZED_SCALE)), cv::Scalar(255), 1);
		}
	}else if(scale > 0){
		for (uint32_t i = 0; i < 256; i++) {
			cv::line(img, cv::Point(XOFFSET+i, YOFFSET), cv::Point(XOFFSET+i, YOFFSET-histogram[i]/scale), cv::Scalar(255), 1);
		}
	}
	
//	for (uint32_t i = 0; i < 256; i++) {
//			line(img, Point(XOFFSET+i, YOFFSET), Point(XOFFSET+i, YOFFSET-histogram[i]/scale), Scalar(255), 1);
//	}

	for (int i = 0; i <= 256; i += 32) {
		cv::line(img, cv::Point(XOFFSET+i, YOFFSET), cv::Point(XOFFSET+i, YOFFSET+20), cv::Scalar(255), 1);
	}
	cv::imshow(histStr, img);
	return img;
//	cv::imwrite("../fig/histogram/gradMagMed.png",img);
}

cv::Mat histogram(const char* histStr, cv::Mat img, int height, int width, uint8_t threshold) {
	cv::Mat output;
	unsigned char value = 0;
	uint16_t k = 256;
	uint16_t histogram[256]; // histogram array - remember to set to zero initially
	uint16_t max = 0;
	uint32_t size = img.cols*img.rows;
//	cv::namedWindow(histStr, CV_WINDOW_NORMAL);
//	cv::resizeWindow(histStr, width, height);

	while (k-- > 0)
		histogram[k] = 0; // reset histogram entry
	for (int i = 0; i < (img.cols*img.rows); i++)
	{
		value = img.data[i];
		histogram[value] += 1;
		if(histogram[value]>max){
			max = histogram[value];
		}
	}

	img = cv::Mat(cv::Size(width,height), CV_8UC3);
	cv::rectangle(img, cv::Point(0, 0), cv::Point(width - 1, height - 1), cv::Scalar(0), CV_FILLED);

	for (int i = 0; i < 256; i += 64) {
		cv::line(img, cv::Point(XOFFSET + i, YOFFSET + 10), cv::Point(XOFFSET + i + 64, YOFFSET + 10), cv::Scalar(i + 20,i+20,i+20), 3);
	}

	for (uint32_t i = 0; i < 256; i++) {
		cv::line(img, cv::Point(XOFFSET+i, YOFFSET), cv::Point(XOFFSET+i, (uint16_t)(YOFFSET-(float)histogram[i]/max*NORMALIZED_SCALE)), cv::Scalar(255,255,255), 1);
	}

	for (int i = 0; i <= 256; i += 32) {
		cv::line(img, cv::Point(XOFFSET+i, YOFFSET), cv::Point(XOFFSET+i, YOFFSET+20), cv::Scalar(255,255,255), 1);
	}
	cv::line(img,cv::Point(XOFFSET+threshold,YOFFSET-20),cv::Point(XOFFSET+threshold,YOFFSET+40),cv::Scalar(0,0,255),2);
	cv::imshow(histStr, img);
//	cv::imwrite("../fig/histogram/gradMagMean.png",img);
	return img;
}

uint16_t histogram(cv::Mat img, float histP[], uint8_t normalized){
	uint8_t temp;
	uint16_t grayMax = 0;
	uint16_t k = 256;
//	uint16_t length;
	//float histogram[k];
	uint32_t nonZero = 0;
	uint32_t size = img.cols*img.rows;
//	float tempHist[k];
	
	while(k-- > 0){
		histP[k] = 0.0;
	}
	for(uint32_t i = 0; i < size; i++){
		temp = img.data[i];
		histP[temp] += 1.0;
		grayMax = temp>grayMax ? temp:grayMax;
	}
	for(uint16_t i = 0; i < 256; i++){
		histP[i] = histP[i] / size;
	}

//	std::ofstream myfile;
//	myfile.open("../fig/test/otsu/histogramData.txt");
//	myfile << "number\t" << "percentage\n";
//	for(uint16_t i = 0; i < grayMax; i++){
//		myfile << i << "\t" << histP[i] << "\n";
//	}
//	myfile.close();;
	return grayMax;
}

void thresholdImage(cv::Mat img, uint8_t threshold) {
	uint8_t *channel;
	channel = (unsigned char*)img.data;
	for (uint32_t i = 0; i < (img.cols)*(img.rows); i++) {
		if ((unsigned char)img.data[i] < threshold) {
			img.data[i] = PIXEL_MIN;
		} else {
			img.data[i] = PIXEL_MAX;
		}
	}
}

void thresholdImage(cv::Mat img, uint8_t lowerT, uint8_t upperT) {
	for (uint32_t i = 0; i < (img.cols)*(img.rows); i++) {
		if (img.data[i] >= upperT) {
			img.data[i] = EDGE_STRONG;
		} else if(img.data[i] >=lowerT) {
			img.data[i] = EDGE_WEAK;
		} else {
			img.data[i] = PIXEL_MIN;
		}
	}
}

void getMoments(cv::Mat img, uint32_t moment[6]) {
	//moment[6] M00, M10, M01, M11, M20, M02
	uint32_t height = img.rows;
	uint32_t width = img.cols;

	for (uint32_t i = 0; i < height; i++) {
		for (uint32_t j = 0; j < width; j++) {
			if ((unsigned char)img.data[j + i * width] == PIXEL_MAX) {
				uint32_t tempJ = j + 1;
				uint32_t tempI = i + 1;
				moment[1] += tempJ;
				moment[2] += tempI;
				moment[3] += tempI * tempJ;
				moment[4] += tempJ * tempJ;
				moment[5] += tempI * tempI;
				moment[0]++;
			}
		}
	}
	printf("%d\n", moment[0]);
}

void printMoments(uint32_t moment[6]) {
	printf("M00 :%d\nM01 :%d\nM10 :%d\nM11 :%d\nM02 :%d\nM20 :%d\n",
		moment[0], moment[1], moment[2], moment[3], moment[4], moment[5]);
}

double getAngle(uint32_t moment[6]) {
	double mu11 = (double) (moment[3] - moment[1] * moment[2] / moment[0]);
	double mu20 = (double) (moment[4] - moment[1] * moment[1] / moment[0]);
	double mu02 = (double) (moment[5] - moment[2] * moment[2] / moment[0]);
	double mu11R = (double) (mu11 / moment[0]);
	double mu20R = (double) (mu20 / moment[0]);
	double mu02R = (double) (mu02 / moment[0]);

	return atan(2 * mu11R / (mu20R - mu02R)) / 2 * 180 / PI;
}

void filter3x3(cv::Mat original, cv::Mat filter, int32_t stencil[3][3]) {
	uint8_t temp = 0;
	uint32_t height = original.rows;
	uint32_t width = original.cols;
	int32_t sum = 0;
	int8_t divisor = 0;

	
	for (uint8_t m = 0; m < 3; m++) {
		for (uint8_t n = 0; n < 3; n++) {
			divisor += stencil[n][m];
		}
	}
	if (divisor == 0) {
		divisor = 1;
	}

	for (uint32_t i = 1; i < height - 1; i++) {
		for (uint32_t j = 1; j < width - 1; j++) {
			sum = 0;
			for (uint32_t k = 0; k < 3; k++) {
				for (uint32_t l = 0; l < 3; l++) {
					temp = original.data[(j + (l - 1)) + (i + (k - 1)) * width];
					sum += (int32_t) (temp) * stencil[l][k];
				}
			}
			filter.data[j + i * width] = (sum<0) ?  0:(uint8_t) (sum / divisor);
		}
	}
}

float median(float m[], int n){
    int         i, less, greater, equal;
    float  min, max, guess, maxltguess, mingtguess;

    min = max = m[0];
    for (i=1 ; i<n ; i++) {
        if (m[i]<min) min=m[i];
        if (m[i]>max) max=m[i];
    }

    while (1) {
        guess = (min+max)/2;
        less = 0; greater = 0; equal = 0;
        maxltguess = min ;
        mingtguess = max ;
        for (i=0; i<n; i++) {
            if (m[i]<guess) {
                less++;
                if (m[i]>maxltguess) maxltguess = m[i] ;
            } else if (m[i]>guess) {
                greater++;
                if (m[i]<mingtguess) mingtguess = m[i] ;
            } else equal++;
        }
        if (less <= (n+1)/2 && greater <= (n+1)/2) break ; 
        else if (less>greater) max = maxltguess ;
        else min = mingtguess;
    }
    if (less >= (n+1)/2) return maxltguess;
    else if (less+equal >= (n+1)/2) return guess;
    else return mingtguess;
}

uint8_t otsu(cv::Mat img) {
	uint16_t grayMax;
	uint32_t length = img.rows*img.cols;
	//float max_between, mean[256], prob[256], between[256];
	
	/*
	prob = class probability
	mean = class mean
	between = between class variance
	*/
	uint8_t threshold;
	float w0, w1, mu0, mu1, sum0, muT, between, max_between;
	float histP[256];
	
	grayMax = histogram(img, histP, 1); // probability distribution, i.e. normalized
	muT = 0.0;
	for(uint16_t i = 0; i < 256; i++){
		muT += i*histP[i];
	}
	//printf("%f\n",sum);
	w0 = 0.0;
	w1 = 0.0;
	sum0 = 0.0;
	max_between = 0.0;
	for(uint16_t i = 0; i < grayMax; i++){
		w0 += histP[i]; // probability class 0
		if(w0==0){
			//printf("continue\n");
			continue;
		}
		w1 = 1-w0; 			// probability class 1
		if(w1==0){
			//printf("break\n");
			break;
		}
		
		sum0 += i*histP[i];
		mu0 = sum0/w0;
		float sum1 = 0.0;
//		for(uint16_t j = i+1; j < graymax; j++){
//			sum1 += i*histp[i];
//		}
		mu1 = (muT-sum0)/w1;
//		between = w0*pow(mu0-muT,2) + w1*pow(mu1-muT,2);
		between = w0*w1*pow(mu0-mu1,2);
		if(between > max_between){
			max_between = between;
			threshold = i;
		}
	}
//	threshold = 0;
//	max_between = 0;
//	for(uint8_t i = 0; i < 256; i++){
//		prob[i] = 0.0;
//		mean[i] = 0.0;
//		between[i] = 0.0;
//	}
//	prob[0] = histogram[0];
//	for(uint8_t i=0; i < 256; i++){
//		prob[i] = prob[i-1] + histogram[i]; // updating w_i
//		mean[i] = mean[i-1] + i*histogram[i]; // updating numerator part of mu_i; mu_T = mean[255]
//	}
//	for(uint8_t i = 0; i < 256; i++){
//		if(prob[i]!=0.0 && prob[i]!=1.0){
//			between[i] = prob[i]*(1-prob[i])*pow( // very wrong; try instead:
//			/*
//			w0*mu0+w1*mu1=muT AND w0+w1=1 => w0*mu0+(1-w0)*mu1=w0*mu0[255]
//			w0*mu0=w0*mu0[255]+(w0-1)*mu1 => mu0=mu0[255]+(1-1/w0)*mu1 => mu0=mu0[255]+mu1-mu1/w0 OR mu1=(w0*mu0[255]-w0*mu0)/(1-w0)
//			sig^2=w0*w1*(mu0-mu1)^2 => sig^2=w0*w1*(mu0[255]-mu1/w0)^2
//			sig^2=(w0^2*w1^2*mu0[255]-w0^2*w1^2*mu1/w0)^2
//			*/
//			//between[i] = pow(mean[i]-muT,2)
//		}else{
//			between[i] = 0.0;
//		}
//		if(between[i] >= max_between){
//			max_between = [between[i];
//			threshold = i;
//		}
	//histogram("mean threshold",img,520,380,threshold);
	return threshold;
}

uint8_t otsuMed(cv::Mat img){
	uint8_t threshold;
	uint16_t grayMax;
	uint32_t length = img.rows*img.cols;
	float prob, probTotal, w0, w1, mT, m0, m1, between, max_between;
	float histP[256];
	w0 = 0.0;
	threshold = 0;
	max_between = 0.0;
	grayMax = histogram(img, histP, 1); // probability distribution
//	for(uint16_t i = 0; i < sizeof(histF)/sizeof(float); i++){
//		histP[i] = histF[i]/length;
//	}
//	std::ofstream myfile;
//	myfile.open("../data/otsuCalc.txt");
//	myfile << "index\t" << "prob\nt" << "med0\t" << "med1\n";
	
	prob = 0.0;
	for(uint16_t i = 0; i < grayMax; i++){
		prob += histP[i];
		if(prob>0.5){
			mT = i;
			break;
		}else if(prob==0.5){
			mT = (2*i+1)/2;
			break;
		}
//		printf("%d\t%f\t%f\n",i,histP[i],prob);
	}
//	printf("i\tw0\t\tw1\t\tm0\t\tm1\t\tprobT\t\tbetween\n");
	for(uint16_t i = 0; i < grayMax; i++){
		w0 += histP[i];
		if(w0 == 0){
			continue;
		}
		w1 = 1-w0;
		if(w1 == 0){
			break;
		}
		
		prob = 0.0;
		for(uint16_t j = 0; j <= i; j++){
			prob += histP[j];
			if(prob > w0/2){
				m0 = j;
				break;
			}else if(prob == w0/2){
				m0 = j+0.5; // average of j and j+1
				break;
			}
		}
//		printf("%d\t%f\t%f\t%f\t",i,w0,w1,m0);
		prob = 0.0;
		for(uint16_t j = i+1; j < grayMax; j++){
			prob += histP[j];
			if(prob > w1/2){
				m1 = j;
				break;
			}else if(prob == w1/2){
				m1 = j-0.5; // average of j and j-1
				break;
			}
		}
//		myfile << i << "\t" << histP[i] << "\t" << m0 << "\t" << m1 << "\n";
		between = w0*pow(m0-mT,2) + w1*pow(m1-mT,2);
//		printf("%f\t%f\t%f\n",m1,mT,between);
		//printf("%f\t", m0);
		if(between > max_between){
			max_between = between;
			threshold = i;
		}
	}
//	myfile << 255 << "\t" << histP[255] << "\t" << NULL << "\t" << NULL << "\n";
//	printf("medT:%f\n",mT);
//	myfile.close();
	printf("%f\n",mT);
	return threshold;
}

uint8_t percThreshold(cv::Mat img, float desiredPerc){
	float currentPerc;
	float histP[256];
	histogram(img, histP, 1);
	currentPerc = 0;
	for(uint16_t i = 255; i >= 0; i--){
		currentPerc += histP[i];
//		printf("currentP:%f\tthreshold:%d\n",currentperc,i);
		if(currentPerc>=desiredPerc){
			return i;
		}
	}
	return 0;
}

void cannyDetector(cv::Mat original, cv::Mat output, uint8_t gaussDim, float percent){//, uint8_t lowerT, uint8_t upperT){
	cv::Mat frameGauss, frameGradMag, frameGradAng;
	uint32_t height = original.rows;
	uint32_t width = original.cols;
	uint8_t meanT, medT, percT;
	
	frameGradMag = cv::Mat(cv::Size(width,height), CV_8UC1);
	frameGradAng = cv::Mat(cv::Size(width,height), CV_8UC1);
	
	/* Gauss LP Filtering */
	cv::GaussianBlur(original, frameGauss, cv::Size(gaussDim,gaussDim),0,0);
	
	gradientMagAngle(frameGauss, frameGradMag, frameGradAng);
	
	meanT = otsu(frameGradMag); // not sure if this will work particularly well
	medT = otsuMed(frameGradMag); // think there's an error when calculating median
	percT = percThreshold(frameGradMag, percent);
	printf("mean:%d\tmedian:%d\tpercent:%d\n",meanT,medT,percT);
//	histogram("log",frameGradMag,520,380,-2);
	
	//cv::imwrite("../data/test.png",frameGradMag);
//	histogram("min suppresion",frameGradMag,520,380,-1);
//	histogram("median threshold",frameGradMag,520,380,medT);
//	histogram("mean threshold",frameGradMag,520,380,meanT);
//	histogram("percent threshold",frameGradMag,520,380,percT);
	
	nonMaxSuppression(frameGradMag, frameGradAng, output);
	
	thresholdImage(output,percT/2,percT);
	//thresholdImage(output, 30);
//	cv::imwrite("../fig/report/Implementation/thresholdMean.png",histogram("mean threshold",frameGradMag,520,380,meanT));
//	cv::imwrite("../fig/report/Implementation/thresholdMedian.png",histogram("median threshold",frameGradMag,520,380,medT));
//	cv::imwrite("../fig/report/Implementation/thresholdPercent.png",histogram("min suppresion",frameGradMag,520,380,percT));
//	cv::namedWindow("gradmag");
//	cv::imshow("gradmag",frameGradMag);
}


void cannyDetector(cv::Mat original, cv::Mat output, uint8_t gaussDim, uint8_t lowerT, uint8_t upperT){//, uint8_t lowerT, uint8_t upperT){
	cv::Mat frameGauss, frameGradMag, frameGradAng;
	uint32_t height = original.rows;
	uint32_t width = original.cols;
	uint8_t meanT, medT, percT;
	
	frameGradMag = cv::Mat(cv::Size(width,height), CV_8UC1);
	frameGradAng = cv::Mat(cv::Size(width,height), CV_8UC1);
	
	/* Gauss LP Filtering */
	cv::GaussianBlur(original, frameGauss, cv::Size(gaussDim,gaussDim),0,0);
	
	gradientMagAngle(frameGauss, frameGradMag, frameGradAng);

	//cv::imwrite("../data/test.png",frameGradMag);
//	histogram("min suppresion",frameGradMag,520,380,-1);
//	histogram("median threshold",frameGradMag,520,380,medT);
//	histogram("mean threshold",frameGradMag,520,380,meanT);
//	histogram("percent threshold",frameGradMag,520,380,percT);
	//histogram("upper threshold",frameGradMag,520,380,upperT);

	nonMaxSuppression(frameGradMag, frameGradAng, output);
	
	thresholdImage(output,lowerT,upperT);
	//thresholdImage(output, 30);
	histogram("thresh",frameGradMag,520,380,upperT);
}

void printImageData(const char* fileLocation, cv::Mat frame){
	std::ofstream myfile;
	myfile.open(fileLocation);
	myfile << "number\t" << "grayscaleValue\n";
	for(uint32_t i = 0; i < frame.rows*frame.cols; i++){
		myfile << i << "\t" << (int)frame.data[i] << "\n";
	}
	myfile.close();
}

void gradientMagAngle(cv::Mat input, cv::Mat frameGradMag, cv::Mat frameGradAng){
	cv::Mat frameX, frameX_abs, frameY, frameY_abs, tempGradMag;
	uint32_t height = input.rows;
	uint32_t width = input.cols;
	float tempAng, tempMag, maxMag;
	int16_t x, y;
	int16_t kX[] = {1,0,-1, 2,0,-2, 1,0,-1};//{3,0,-3,10,0,-10,3,0,-3};//
	int16_t kY[] = {1,2,1, 0,0,0, -1,-2,-1};//{3,10,3,0,0,0,-3,-10,-3};//
	int16_t kX5[] = {1,2,0,-2,-1, 	4,8,0,-8,-4, 	6,12,0,-12,-6, 	4,8,0,-8,-4, 	1,2,0,-2,-1};
	int16_t kY5[] = {1,4,6,4,1, 	2,8,12,8,2, 	0,0,0,0,0,0, 	-2,-8,-12,-8,-2, 	-1,-4,-6,-4,-1};
	cv::Mat sobelX = cv::Mat(3,3,CV_16SC1, kX);
	cv::Mat sobelY = cv::Mat(3,3,CV_16SC1, kY);
	cv::Mat sobelX5 = cv::Mat(5,5,CV_16SC1, kX5);
	cv::Mat sobelY5 = cv::Mat(5,5,CV_16SC1, kY5);

	frameX = cv::Mat(cv::Size(width,height), CV_16SC1);
	frameY = cv::Mat(cv::Size(width,height), CV_16SC1);
	tempGradMag = cv::Mat(cv::Size(width,height), CV_16SC1);
	cv::filter2D(input, frameX, CV_16S, sobelX);
	cv::filter2D(input, frameY, CV_16S, sobelY);
	cv::convertScaleAbs(frameX,frameX_abs);
	cv::convertScaleAbs(frameY,frameY_abs);
	
//	cv::addWeighted(frameX_abs,1,frameY_abs,1,0,frameGradMag);
	
//	cv::phase(frameX,frameY,frameGradAng,1); // inDegrees true
//	for(uint16_t i = 0; i < 10000; i++){
//		frameGradAng.data[i] = 255;
//	}
//	float tempMag[width*height]
	maxMag = 0;
	for (uint32_t i = 0; i < width*height; i++) {
		memcpy(&x, &(frameX.data[i*sizeof(x)]), sizeof(x));
		memcpy(&y, &(frameY.data[i*sizeof(y)]), sizeof(y));

		tempMag = hypot(x,y);
		tempGradMag.data[i*sizeof(int16_t)] = (int16_t)tempMag;
		maxMag = tempMag > maxMag ? tempMag : maxMag;		
//		frameGradMag.data[i] = tempMag>255 ? 255 : (uint8_t)tempMag;

		tempAng = atan2(y,x)*180/PI;		
//		maxMag = tempMag[i] > maxMag ? tempMag[i] : maxMag;
		tempAng = tempAng < 0 ? tempAng+180 : tempAng;
		/* 
		using positive y axis as down on the screen
		therefore North is down on the screen and south is up.
		Positive x axis is to the right on the screen
		therefore east is to the right and west to the left 
		
		angle goes from positiv x asxis towards positive y
		there clock wise
		*/
		if(tempAng <= 22.5 || tempAng >= 157.5){
			frameGradAng.data[i] = 0;
//			printf("%d\n",frameGradAng.data[i]);
		}else if(tempAng < 67.5){
			frameGradAng.data[i] = 45;
//			printf("%d\n",frameGradAng.data[i]);
		}else if(tempAng <= 112.5){
			frameGradAng.data[i] = 90;
//			printf("%d\n",frameGradAng.data[i]);
		}else{
			frameGradAng.data[i] = 135;
//			printf("%d\n",frameGradAng.data[i]);
		
		}
	}
//	cv::convertScaleAbs(tempGradMag,frameGradMag,(double)255/maxMag);
	for(uint32_t i = 0; i < width*height; i++){ // the same as convert scale abs by inspection
		frameGradMag.data[i] = (uint8_t)abs(tempGradMag.data[i*sizeof(int16_t)]*255/maxMag);
	}
//	for(uint32_t i = 0; i < width*height; i++){
//		frameGradMag.data[i] = (uint8_t)(255*tempMag[i]/maxMag);
//	}
//	histogram("grad mag", frameGradMag, 520, 380, (uint8_t)32);
//	printImageData("../fig/test/otsu/gradMagOtsu.txt",frameGradMag);
//	cv::imwrite("../fig/test/otsu/gradMagOtsu.png",frameGradMag);
//	cv::imshow("frameGradMag", frameGradMag);
//	cv::imwrite("../fig/report/Implementation/gradMagWeightedOverflow_1.png", frameGradMag);

}

void nonMaxSuppression(cv::Mat frameGradMag, cv::Mat frameGradAng, cv::Mat output){
	uint32_t height = frameGradMag.rows;
	uint32_t width = frameGradMag.cols;

	for (uint32_t i = 1; i < height - 1; i++) {
		for (uint32_t j = 1; j < width - 1; j++) {
			switch(frameGradAng.data[j + i * width]){
				case 0:
					if(frameGradMag.data[j+i*width] >= frameGradMag.data[(j-1)+i*width] 
					&& frameGradMag.data[j+i*width] >= frameGradMag.data[(j+1)+i*width]){
						output.data[j+i*width] = frameGradMag.data[j+i*width];
					}else{
						output.data[j+i*width] = 0;
					}
					break;
				case 90:
					if(frameGradMag.data[j+i*width] >= frameGradMag.data[(j)+(i-1)*width] 
					&& frameGradMag.data[j+i*width] >= frameGradMag.data[j+(i+1)*width]){
						output.data[j+i*width] = frameGradMag.data[j+i*width];
					}else{
						output.data[j+i*width] = 0;
					}
					break;
				case 135:
					if(frameGradMag.data[j+i*width] >= frameGradMag.data[(j+1)+(i-1)*width] 
					&& frameGradMag.data[j+i*width] >= frameGradMag.data[(j-1)+(i+1)*width]){
						output.data[j+i*width] = frameGradMag.data[j+i*width];
					}else{
						output.data[j+i*width] = 0;
					}
					break;
				case 45:
					if(frameGradMag.data[j+i*width] >= frameGradMag.data[(j+1)+(i+1)*width] 
					&& frameGradMag.data[j+i*width] >= frameGradMag.data[(j-1)+(i-1)*width]){
						output.data[j+i*width] = frameGradMag.data[j+i*width];
					}else{
						output.data[j+i*width] = 0;
					}
					break;
				default:
					return;
			}
		}
	}
	//cv::imwrite("../fig/report/Implementation/nonMaxSuppression.png", output);

}

uint32_t findMaxLabels(cv::Mat img, uint32_t *labelMap){
	uint32_t maxLabels = 0;
	uint32_t height = img.rows;
	uint32_t width = img.cols;
	for (uint32_t i = 1; i < height; i++) {
		for (uint32_t j = 1; j < width-1; j++) {
			if ((uint8_t)img.data[j + i * width] > PIXEL_MIN) {
				//printf("index %d has address %p and points at %p\n",labelPointer[0].index, &labelPointer[0], labelPointer[0].point);
				labelMap[j + i * width] = 0;
				if ((uint8_t)img.data[(j-1) + i * width] != PIXEL_MIN){
					labelMap[j + i * width] = labelMap[(j-1) + i * width]; 
				} else if ((uint8_t)img.data[(j-1) + (i-1) * width] != PIXEL_MIN){
					labelMap[j + i * width] = labelMap[(j-1) + (i-1) * width];
				} else if ((uint8_t)img.data[(j) + (i-1) * width] != PIXEL_MIN){
					labelMap[j + i * width] = labelMap[(j) + (i-1) * width];
				} else if ((uint8_t)img.data[(j+1) + (i-1) * width] != PIXEL_MIN){
					labelMap[j + i * width] = labelMap[(j+1) + (i-1) * width];
				}  
				if(labelMap[j + i * width] == PIXEL_MIN) {
					maxLabels++;
					labelMap[j + i * width] = maxLabels;
				}
			} else{
				labelMap[j + i * width]=0;
			}
		}
	}
	return maxLabels;
}

void removeObjects(cv::Mat out, uint32_t *labelMap, LabelIndices *labelPointer, uint32_t minObjectSize, uint32_t maxObjectSize){
	uint32_t height = out.rows;
	uint32_t width = out.cols;
	uint32_t isNotZero;
	uint32_t isNotStrongEdge;
	uint32_t objectSize; 
	maxObjectSize = (maxObjectSize == 0) ? -1 : maxObjectSize;
	for (uint32_t i = 1; i < height; i++) {
		for (uint32_t j = 1; j < width-1; j++) {
			isNotZero = (labelMap[j + i * width] != 0);	
	
			if(isNotZero){
				isNotStrongEdge = !(labelPointer[labelMap[j + i * width]-1].hasStrongEdge);
				objectSize = labelPointer[labelMap[j + i * width]-1].size;	

				if(isNotStrongEdge){
					labelMap[j + i * width] = 0;
				} else if (objectSize < minObjectSize){
					labelMap[j + i * width] = 0; 
				} else if (objectSize > maxObjectSize){
					labelMap[j + i * width] = 0;
				}			
			}
		}
	}
}


void applyMap(cv::Mat out, uint32_t *labelMap){
	uint32_t height = out.rows;
	uint32_t width = out.cols;
	for (uint32_t i = 0; i < height; i++) {
		for (uint32_t j = 0; j < width; j++) {
			if (labelMap[j + i * width] == 0){
				out.data[j + i * width]=0; 
			}
		}
	}
}


uint32_t labeling(cv::Mat img, Objects **objectMap, uint32_t minObjectSize, uint32_t maxObjectSize) {
	uint32_t height = img.rows;
	uint32_t width = img.cols;
	uint32_t *labelMap = (uint32_t*) malloc(sizeof(uint32_t*) * height*width); //size problems if used on stack thus pointers to put it on heap
	uint32_t currentMemoryScale = 1;
	uint32_t currentLabelIndex = 0;
	uint32_t prevMemoryScale;

	//remove edges
	for (uint32_t j = 0; j < width; j++) {		
		img.data[j] = 0 ;
		labelMap[j] = 0;
	}
	for (uint32_t i = 0; i < height; i++) {
		img.data[i * width] = 0 ;
		img.data[(width-1) + i * width] = 0 ;
		labelMap[i * width] = 0;
		labelMap[i * width + (width-1)] = 0;
	}

	uint32_t maxLabels = findMaxLabels(img, labelMap);
	LabelIndices *labelPointer = (LabelIndices*) malloc(maxLabels * sizeof(LabelIndices)); 
	for (uint32_t k = 0; k < maxLabels; k++){
		labelPointer[k].index = k+1;
		labelPointer[k].point = NULL;
		labelPointer[k].size = 0; 
		labelPointer[k].hasStrongEdge = 0;
	}	




    // the interesting areas will be within the area of 1->height and 1->width
	//labeling
	for (uint32_t i = 1; i < height; i++) {
		for (uint32_t j = 1; j < width-1; j++) {
			if ((uint8_t)img.data[j + i * width] > PIXEL_MIN) {
				//printf("index %d has address %p and points at %p\n",labelPointer[0].index, &labelPointer[0], labelPointer[0].point);
				labelMap[j + i * width] = 0;
				if ((uint8_t)img.data[(j-1) + i * width] != PIXEL_MIN){
					labelMap[j + i * width] = labelMap[(j-1) + i * width]; 
				} if ((uint8_t)img.data[(j-1) + (i-1) * width] != PIXEL_MIN){
					if (labelMap[j + i * width] != 0 && labelMap[(j-1) + (i-1) * width] != labelMap[j + i * width]){
						labelPointer[labelMap[(j-1) + (i-1) * width]-1].point = &labelPointer[labelMap[j + i * width]-1];
					} else
						labelMap[j + i * width] = labelMap[(j-1) + (i-1) * width];
				} 
				if ((uint8_t)img.data[(j) + (i-1) * width] != PIXEL_MIN){
					if (labelMap[j + i * width] != 0 && labelMap[(j) + (i-1) * width] != labelMap[j + i * width]){
						labelPointer[labelMap[(j) + (i-1) * width]-1].point = &labelPointer[labelMap[j + i * width]-1];
					} else
						labelMap[j + i * width] = labelMap[(j) + (i-1) * width];
				} 
				if ((uint8_t)img.data[(j+1) + (i-1) * width] != PIXEL_MIN){
					//printf("\t\tpoints at %p\n", labelPointer[0].point);
					if (labelMap[j + i * width] != 0 && labelMap[(j+1) + (i-1) * width] != labelMap[j + i * width]){
						labelPointer[labelMap[(j+1) + (i-1) * width]-1].point = &labelPointer[labelMap[j + i * width]-1]; 
					} else 
						labelMap[j + i * width] = labelMap[(j+1) + (i-1) * width];
				}  
				
				if(labelMap[j + i * width] == 0) {
					currentLabelIndex++;
					labelMap[j + i * width]=currentLabelIndex;
				}
			} else{
				labelMap[j + i * width]=0;
			}
		}
	}
	for (uint32_t k = 0; k < maxLabels; k++){
	//	printf("index %d has address %p and points at %p\n",labelPointer[k].index, &labelPointer[k], labelPointer[k].point);
	}


	for (uint32_t i = 0; i < height; i++) {
		for (uint32_t j = 0; j < width; j++) {
			//printf("%d\t",labelMap[j + i * width]);
		}
		//printf("\n");
	}

	//printf("\n\n");
//set neighbors equal
	LabelIndices *tempLabelPointer;
	for (uint32_t i = 1; i < height; i++) {
		for (uint32_t j = 1; j < width-1; j++) {
			
			if (labelMap[j + i * width] != 0){
				
				tempLabelPointer = &labelPointer[labelMap[j + i * width]-1];
				for( uint32_t k = 0; tempLabelPointer->point != NULL; k++){
					tempLabelPointer = tempLabelPointer->point;
					if (k == 99)
						tempLabelPointer->point = NULL;
				}
				labelMap[j + i * width] = tempLabelPointer->index; 
				labelPointer[labelMap[j + i * width]-1].size++;
			}
			if((uint8_t)img.data[j + i * width] == EDGE_STRONG){
				labelPointer[labelMap[j + i * width]-1].hasStrongEdge++;
			}	
		}		
	}

	for (uint32_t i = 0; i < height; i++) {
		for (uint32_t j = 0; j < width; j++) {
			//printf("%d\t",labelMap[j + i * width]);
		}
		//printf("\n");
	}
	//printf("\n\n");
	
	//Make conneced edges map
	removeObjects(img, labelMap, labelPointer, minObjectSize, maxObjectSize);
	//getEdgeMap();	
	

	//allocate memory for object map
	uint32_t maxObjects = 0;
	uint32_t newObjectIndex = 0;
	uint32_t objectSize = 0;
	uint32_t oldObjectIndex = 0;
	for (uint32_t k = 0; k <maxLabels; k++){
		if (labelPointer[k].size){
			maxObjects++; 				
		}	
	}
	
	uint32_t tempArray[maxLabels];
//	uint32_t tempSize[maxObjects];
//	if (maxObjects){
			*objectMap = (Objects*) realloc(*objectMap , sizeof(Objects)*(maxObjects));
//		if (!objectMap){
//			printf("NULL pointer returned\n");
//		}
//	}
//	else{
//		applyMap(img, labelMap);
//		free(labelMap);
//		free(labelPointer);
//		return 0;
//	} 
	
	//printf("pp %p\n", objectMap);
	//define vector size and new indecies	
	for (uint32_t k = 0; k <maxLabels; k++){
		tempArray[k] = -1;
		if (labelPointer[k].size){	
			objectSize = labelPointer[k].size;
//			tempSize[newObjectIndex] = objectSize;		
			(*objectMap)[newObjectIndex].size = 0;
			(*objectMap)[newObjectIndex].pos = (Vector*) malloc(sizeof(Vector) * objectSize);	
			tempArray[k] = newObjectIndex;
			//printf("%d, %d\n",tempArray[k], newObjectIndex);
			newObjectIndex++;			
		}	
//		printf("%d -> %d     %d \n", k, tempArray[k], objectSize);	
	}

	uint32_t newLabel;
	uint32_t currentPos;
	for (uint32_t i = 0; i < height; i++){
		for (uint32_t j = 0; j < width; j++){
			uint32_t oldLabel = labelMap[j + i * width]-1;		
			if (oldLabel != -1){
				newLabel = tempArray[oldLabel];
				currentPos = (*objectMap)[newLabel].size; 
				(*objectMap)[newLabel].size++;
				(*objectMap)[newLabel].pos[currentPos].x = j;
				(*objectMap)[newLabel].pos[currentPos].y = i;
				//printf("%d,%d, %d, %d, %d\n",labelMap[j + i * width],newLabel, currentPos, (*objectMap)[newLabel].pos[currentPos].x, (*objectMap)[newLabel].pos[currentPos].y);
				labelMap[j + i * width] = newLabel+1;
				if (newLabel == 0) cv::circle(img, cv::Point(j,i), 1, cv::Scalar(255,0,0), 1);
				else if (newLabel == 1) cv::circle(img, cv::Point(j,i), 1, cv::Scalar(50,0,0), 1);
				else if (newLabel == 2) cv::circle(img, cv::Point(j,i), 1, cv::Scalar(150,0,0), 1);
				else if (newLabel == 3) cv::circle(img, cv::Point(j,i), 1, cv::Scalar(255,0,0), 1);
			}
		}
	}
	for (uint32_t i = 0; i < height; i++) {
		for (uint32_t j = 0; j < width; j++) {
			//printf("%d\t",labelMap[j + i * width]);
		}
		//printf("\n");
	}
	
	for (uint32_t k = 0; k <maxObjects; k++){
		//printf("%d,%d\n",k ,(*objectMap)[k].size);
	}
	uint32_t labelSize = maxObjects-1;
	applyMap(img, labelMap);
	free(labelMap);
	free(labelPointer);
	printf("%d\n",labelSize);
	return labelSize;
}



void getCircle(Vector &p1,Vector &p2,Vector &p3, Vector &center, float &radius){
  	float x1 = p1.x;
  	float x2 = p2.x;
  	float x3 = p3.x;

  	float y1 = p1.y;
  	float y2 = p2.y;
  	float y3 = p3.y;
	

  	// PLEASE CHECK FOR TYPOS IN THE FORMULA :)
  	center.x = (x1*x1+y1*y1)*(y2-y3) + (x2*x2+y2*y2)*(y3-y1) + (x3*x3+y3*y3)*(y1-y2);
  	center.x /= ( 2*(x1*(y2-y3) - y1*(x2-x3) + x2*y3 - x3*y2) );

  	center.y = (x1*x1 + y1*y1)*(x3-x2) + (x2*x2+y2*y2)*(x1-x3) + (x3*x3 + y3*y3)*(x2-x1);
  	center.y /= ( 2*(x1*(y2-y3) - y1*(x2-x3) + x2*y3 - x3*y2) );
	
  	radius = sqrt((center.x-x1)*(center.x-x1) + (center.y-y1)*(center.y-y1));
}

uint8_t dynamicRealloc(Vector *point, uint32_t *prev, uint32_t *current){
	if ((*prev<<1) >= (*current)-2){
		*prev = *current;
		point = (Vector*) realloc(point, sizeof(Vector)*((*current)<<1));	
		//printf("%p %lld\n ",point, (*current)<<2);			
		if (point == NULL)
			return 1;		
	}
	return 0;
}


uint32_t getPointPositions(cv::Mat img, Vector *pointPositions){
	uint32_t height = img.rows;
	uint32_t width = img.cols;
	uint32_t counter = 1;
	uint32_t prevCounter = 0;
 	for(uint32_t y=0; y<height; y++){
     	for(uint32_t x=0; x<width; x++){
         	if(img.data[x + y * width] > 0) {
				counter++;
				//dynamicRealloc(pointPositions, &prevCounter, &counter);			
				pointPositions[counter].x = x;
				pointPositions[counter].y = y;
			}
     	}
 	}
	
 	return counter;
}
 //vote implementation
uint8_t verifyRadius(cv::Mat img, Vector &center, float radius, float t, float maxInlierDist, float radiusTolerance){
	float isInsideCircle = 0;

	for (float i = -radiusTolerance; i <= radiusTolerance; i++){
	 	float cX = (radius-i)*cos(t) + center.x;
	 	float cY = (radius-i)*sin(t) + center.y; 		
		isInsideCircle += img.at<float>(cY,cX) < maxInlierDist;
	}
	return isInsideCircle;
}

float verifyCircle(cv::Mat img, Vector &center, float radius, Vector *inlierSet){
	uint32_t height = img.rows;
	uint32_t width = img.cols;	
	uint32_t prevInlier = 1;
 	uint32_t inlier = 0;
	float counter = 0;
 	float minInlierDist = 2.0f;
 	float maxInlierDistMax = 100.0f;
 	float maxInlierDist = radius/25.0f;
	
	if(maxInlierDist<minInlierDist) 
		maxInlierDist = minInlierDist;
 	if(maxInlierDist>maxInlierDistMax) 
		maxInlierDist = maxInlierDistMax;

	// if no cirlces are present look for half circles!
 	for(float t =0; t<2*PI; t+= PI/100){
	 	counter+=0.1; //place ocunter there to count outside image as not inlier
	 	float cX = radius*cos(t) + center.x;
	 	float cY = radius*sin(t) + center.y;
		uint32_t tol = 4; //tolerance each side
	 	if(cX < img.cols-tol && cX >= tol && cY < img.rows-tol && cY >= tol){
			counter++; //count ignore pixels outside the image 
			if(verifyRadius(img, center, radius, t, maxInlierDist, tol)){//img.at<float>(cY,cX) < maxInlierDist	
				inlier++;				  // why doesnt this give almost 100%???
				inlierSet[inlier].x = (uint32_t) cX;
				inlierSet[inlier].y = (uint32_t) cY;
			}
		}
	} 
 	return (float)inlier/float(counter);
}

/*
float verifyCircle(cv::Mat img, Vector &center, float radius, Vector *inlierSet){
	uint32_t height = img.rows;
	uint32_t width = img.cols;	
	uint32_t prevInlier = 1;
 	uint32_t inlier = 0;
	uint32_t counter = 0;
 	float minInlierDist = 2.0f;
 	float maxInlierDistMax = 100.0f;
 	float maxInlierDist = radius/25.0f;
	
	if(maxInlierDist<minInlierDist) 
		maxInlierDist = minInlierDist;
 	if(maxInlierDist>maxInlierDistMax) 
		maxInlierDist = maxInlierDistMax;

 	// choose samples along the circle and count inlier percentage
 	for(float t =0; t<2*PI; t+= PI/100){
     	counter++;
     	float cX = radius*cos(t) + center.x;
     	float cY = radius*sin(t) + center.y;

     	if(cX < img.cols && cX >= 0 && cY < img.rows && cY >= 0){
			if (center.x == 64 && center.y==65)			
				printf("%d, %f, %f\n",counter, img.at<float>(cY,cX),maxInlierDist);
 			if(img.at<float>(cY,cX) < maxInlierDist){ // this might be wrong
				
				inlier++;				
    			inlierSet[inlier].x = (uint32_t) cX;
				inlierSet[inlier].y = (uint32_t) cY;
			}
		}
 	}
 	return (float)inlier/float(counter);
}
*/

Circle ransac(cv::Mat img, float minRadius, float maxRadius, float minCirclePercentage){

	uint32_t height = img.rows;
	uint32_t width = img.cols;	
	float radius;
	int32_t maxNrOfIterations;
	uint32_t maxSize = 2;
	Circle bestCircle;
	for (uint32_t i = 0; i < height; i++){		
		for (uint32_t j = 0; j < width; j++){
			if (img.data[j + i * width]){
				maxSize++;
			}
		}
	}
	
    Vector *pointPositions = (Vector*) malloc(sizeof(Vector)*maxSize);
	Vector *inlierSet = (Vector*) malloc(sizeof(Vector)*maxSize);
	Vector center; 

    maxNrOfIterations = getPointPositions(img, pointPositions); //returns high pixels

    // create distance transform to efficiently evaluate distance to nearest edge
    cv::Mat dt;
    cv::distanceTransform(255-img, dt,CV_DIST_L1, 3);

    //TODO: maybe seed random variable for real random numbers.

    uint32_t nIterations = 0; // TODO: adjust this parameter or include some real ransac criteria with inlier/outlier percentages to decide when to stop

    Vector bestCircleCenter;
    float bestCircleRadius;
    float bestCirclePercentage = 0;
	//float maxRadius = 100;
    //float minRadius = 50;   // TODO: ADJUST THIS PARAMETER TO YOUR NEEDS, otherwise smaller circles wont be detected or "small noise circles" will have a high percentage of completion


    for(uint32_t its=0; its< maxNrOfIterations*1; its++){
        //RANSAC: randomly choose 3 point and create a circle:
        //TODO: choose randomly but more intelligent, 
        //so that it is more likely to choose three points of a circle. 
        //For example if there are many small circles, it is unlikely to randomly choose 3 points of the same circle.
        uint32_t idx1 = rand()%maxNrOfIterations;
        uint32_t idx2 = rand()%maxNrOfIterations;
        uint32_t idx3 = rand()%maxNrOfIterations;

        // we need 3 different samples:
        if(idx1 == idx2) 
			continue;
        if(idx1 == idx3) 
			continue;
        if(idx3 == idx2) 
			continue;

        // create circle from 3 points:
        
        getCircle(pointPositions[idx1],pointPositions[idx2],pointPositions[idx3],center,radius);

        // inlier set unused at the moment but could be used to approximate a (more robust) circle from alle inlier
       

        //verify or falsify the circle by inlier counting:
        float cPerc = verifyCircle(dt,center,radius,inlierSet);

        // update best circle information if necessary
        if(cPerc >= bestCirclePercentage)
        if(radius >= minRadius && radius <= maxRadius){
            bestCirclePercentage = cPerc;
            bestCircleRadius = radius;
            bestCircleCenter = center;
			//printf("circle my dood");

//			printf("%d %d %f %f\n",bestCircleCenter.x,bestCircleCenter.y,bestCircleRadius,cPerc);		
			
			//printf("\n percent %f",cPerc);
        }
    }

	if(bestCirclePercentage >= minCirclePercentage && bestCircleRadius >= minRadius && bestCircleRadius <= maxRadius)	{
//		cv::circle(img, cv::Point(bestCircleCenter.x, bestCircleCenter.y),bestCircleRadius, cv::Scalar(100,100,0),1);
		bestCircle.x = bestCircleCenter.x;
		bestCircle.y = bestCircleCenter.y;
		bestCircle.radius = bestCircleRadius;
//		printf("\n yee %d %d %f",bestCircleCenter.x,bestCircleCenter.y,bestCircleRadius);
		//printf("\n yee %d %d %f %f\n",bestCircle.x,bestCircle.y,bestCircle.radius, bestCirclePercentage);
//		cv::circle(img, cv::Point(bestCircleCenter.x, bestCircleCenter.y),bestCircleRadius, cv::Scalar(100,100,0),1);
	} else
		bestCircle.radius = 0;
	

	
	free(pointPositions);
    free(inlierSet);
    //cv::imshow("output",img);
    //cv::imshow("mask",mask);
    //cv::waitKey(0);
	//imshow("circles", img );
    return bestCircle;
}



Circle ransac(cv::Mat img, Objects *map, uint32_t objects, float minRadius, float maxRadius, float minCirclePercentage){

	uint32_t height = img.rows;
	uint32_t width = img.cols;	
	float radius;
	int32_t maxNrOfIterations;
	uint32_t maxSize = 2;
	uint32_t maxCircles = 4000;
	Circle bestCircle;
	for (uint32_t i = 0; i < height; i++){		
		for (uint32_t j = 0; j < width; j++){
			if (img.data[j + i * width]){
				maxSize++;
			}
		}
	}
	Circle *circles = (Circle*) calloc(maxCircles, sizeof(Circle));
    Vector *inlierSet = (Vector*) malloc(sizeof(Vector)*maxSize);
	Vector center; 


	
    // create distance transform to efficiently evaluate distance to nearest edge
    cv::Mat dt;
    cv::distanceTransform(255-img, dt,CV_DIST_L1, 3);

    //TODO: maybe seed random variable for real random numbers.

    uint32_t nIterations = 0; // TODO: adjust this parameter or include some real ransac criteria with inlier/outlier percentages to decide when to stop

    Vector bestCircleCenter;
    float bestCircleRadius;
    float bestCirclePercentage = 0;
	
	//float maxRadius = 100;
    //float minRadius = 50;   // TODO: ADJUST THIS PARAMETER TO YOUR NEEDS, otherwise smaller circles wont be detected or "small noise circles" will have a high percentage of completion

	uint32_t circleIndex = 0;
	// maybe add the +1
	uint32_t offset = 3; //3, 9, 17
	for(uint32_t i=0; i< objects; i++){
		//printf("%d %d %d\n",i , objects,map[i].size);
				

		for(uint32_t j=0; j< map[i].size*1 ; j++){
			
			uint32_t idx1 = rand()%map[i].size;
			uint32_t idx2 = rand()%map[i].size;
			uint32_t idx3 = rand()%map[i].size;

			// we need 3 different samples:
			if(idx1 == idx2) 
				continue;
			if(idx1 == idx3) 
				continue;
			if(idx3 == idx2) 
				continue;
			//printf("%d,%d,%d\n", idx1, idx2,idx3);
			
			// create circle from 3 points:
			getCircle(map[i].pos[idx1],map[i].pos[idx2],map[i].pos[idx3],center,radius);

			// looking for circle bug
			//cv::circle(img, cv::Point(map[i].pos[idx1].x,map[i].pos[idx1].y), 0, cv::Scalar(a,b,c), 2);
			///

			// inlier set unused at the moment but could be used to approximate a (more robust) circle from alle inlier
		   

			//verify or falsify the circle by inlier counting:
			float cPerc = verifyCircle(dt,center,radius,inlierSet);
					
			//
			// it only checks within the object
			//


			//printf("%u %u %f %f\n",center.x,center.y,radius,cPerc);		
			// update best circle information if necessary
			if(radius >= minRadius && radius <= maxRadius){
				//printf("%d %d %f %f\n",center.x,center.y,radius,cPerc);
				if (cPerc > minCirclePercentage && circleIndex < maxCircles){
					//printf("%d\n",circleIndex);
					circles[circleIndex].x = center.x;
					circles[circleIndex].y = center.y;
					circles[circleIndex].radius = radius;
					circles[circleIndex].perc = cPerc;
					circleIndex++;		
					//cv::circle(img, cv::Point(center.x, center.y),radius, cv::Scalar(100,100,0),1);			
				}	
				if(cPerc >= bestCirclePercentage){
					bestCirclePercentage = cPerc;
					bestCircleRadius = radius;
					bestCircleCenter = center;
					//printf("circle my dood");
					//printf("%d %d %f %f\n",bestCircleCenter.x,bestCircleCenter.y,bestCircleRadius,cPerc);		
					//cv::circle(img, cv::Point(center.x, center.y),radius, cv::Scalar(255,100,0),1);	
					//printf("\n percent %f",cPerc);
				}
			}
		}
	}
	


	CircleOccurrence similarCircles[circleIndex];
	uint32_t similarCirclesIndex = 1;
	uint32_t pScale = 6;
	uint32_t rScale = 5;
	similarCircles[0].x = circles[0].x;
	similarCircles[0].y = circles[0].y;
	similarCircles[0].radius = circles[0].radius;
	similarCircles[0].occurrence = 1;
	for (uint32_t i = 1; i < circleIndex; i++){
		for (uint32_t j = 0; j <= similarCirclesIndex; j++){
			//printf("%d,%d %d\n",j, i, similarCirclesIndex);	
			if ((uint32_t)circles[i].x>>pScale == (uint32_t)similarCircles[j].x>>pScale
				&& (uint32_t)circles[i].y>>pScale == (uint32_t)similarCircles[j].y>>pScale
				&& (uint32_t)circles[i].radius>>rScale == (uint32_t)similarCircles[j].radius>>rScale) {
				similarCircles[j].occurrence++;
				break;			
			} else if (j == similarCirclesIndex ){
				similarCirclesIndex++;
				similarCircles[j].occurrence =1;
				similarCircles[j].perc = circles[i].perc;//0;//
				similarCircles[j].x = circles[i].x;//0;//
				similarCircles[j].y = circles[i].y;//0;//
				similarCircles[j].radius = circles[i].radius;//0;//
				//printf("%d, %d, %d, %f %d\n", (uint32_t)similarCircles[j].x>>5, (uint32_t)similarCircles[j].y>>5, ((uint32_t)similarCircles[j].radius)>>4, similarCircles[j].perc, similarCircles[j].occurrence);
				break;
			}

		}
		
		
		//printf("%d, %d, %d, %f, %d\n", (uint32_t)circles[i].x>>5, (uint32_t)circles[i].y>>5, ((uint32_t)circles[i].radius)>>4, circles[i].perc, i);
	}	
	
	uint32_t maxOccurrences = 0;
	CircleOccurrence finalCircles[similarCirclesIndex];
		
	for (uint32_t i = 0; i < similarCirclesIndex; i++){

		//printf("%d, %d, %d, %f %d\n", (uint32_t)similarCircles[i].x>>pScale, (uint32_t)similarCircles[i].y>>pScale, ((uint32_t)similarCircles[i].radius)>>rScale, similarCircles[i].perc, similarCircles[i].occurrence);
		finalCircles[i].perc = 0;//circles[i].perc;//
		finalCircles[i].x = 0;//circles[i].x;//
		finalCircles[i].y = 0;//circles[i].y;//
		finalCircles[i].radius = 0;//circles[i].radius;
		finalCircles[i].occurrence = similarCircles[i].occurrence;
		if (maxOccurrences < finalCircles[i].occurrence)
			maxOccurrences = finalCircles[i].occurrence;
	}

	for (uint32_t i = 0; i < circleIndex; i++){
		for (uint32_t j = 0; j <= similarCirclesIndex; j++){
			//printf("%d,%d %d\n",j, i, similarCirclesIndex);	
			if ((uint32_t)circles[i].x>>pScale == (uint32_t)similarCircles[j].x>>pScale
				&& (uint32_t)circles[i].y>>pScale == (uint32_t)similarCircles[j].y>>pScale
				&& (uint32_t)circles[i].radius>>rScale == (uint32_t)similarCircles[j].radius>>rScale) {
				finalCircles[j].perc += circles[i].perc;
				finalCircles[j].x += circles[i].x;
				finalCircles[j].y += circles[i].y;
				finalCircles[j].radius += circles[i].radius;
				//printf("%d, %d, %d, %f %d\n", (uint32_t)finalCircles[j].x>>5, (uint32_t)finalCircles[j].y>>5, ((uint32_t)finalCircles[j].radius)>>4, finalCircles[j].perc, finalCircles[j].occurrence);
				break;
			}
		}
	}

	//cv::cvtColor(frame, frameCircle, CV_GRAY2RGB);
 	
	//printf("\n\n");


	
	
	//remove or rework this
	if(bestCirclePercentage >= minCirclePercentage && bestCircleRadius >= minRadius && bestCircleRadius <= maxRadius)	{
//		cv::circle(img, cv::Point(bestCircleCenter.x, bestCircleCenter.y),bestCircleRadius, cv::Scalar(100,100,0),1);
		bestCircle.x = bestCircleCenter.x;
		bestCircle.y = bestCircleCenter.y;
		bestCircle.radius = bestCircleRadius;
//		cv::circle(img, cv::Point(bestCircleCenter.x, bestCircleCenter.y),bestCircleRadius, cv::Scalar(100,100,0),1);
	} else{
		bestCircle.x = -1;
		bestCircle.y = -1;
		bestCircle.radius = 0;
		printf("No acceptable circle found\nbest percentage: %d\n ",bestCirclePercentage);
		return bestCircle; 
	}

	



	
	printf("x: %f pixels \n y: %f pixels\n radius: %f pixels\n percentage: %f\n",bestCircle.x,bestCircle.y,bestCircle.radius, 100*bestCirclePercentage);

    free(inlierSet);
    //cv::imshow("output",img);
    //cv::imshow("mask",mask);
    //cv::waitKey(0);
	//imshow("circles", img );
    return bestCircle;
}



Circle ransacVoting(cv::Mat img, Objects *map, uint32_t objects, float minRadius, float maxRadius, float minCirclePercentage){

	uint32_t height = img.rows;
	uint32_t width = img.cols;	
	float radius;
	int32_t maxNrOfIterations;
	uint32_t maxSize = 2;
	uint32_t maxCircles = 4000;
	Circle bestCircle;
	for (uint32_t i = 0; i < height; i++){		
		for (uint32_t j = 0; j < width; j++){
			if (img.data[j + i * width]){
				maxSize++;
			}
		}
	}
	Circle *circles = (Circle*) calloc(maxCircles, sizeof(Circle));
    Vector *inlierSet = (Vector*) malloc(sizeof(Vector)*maxSize);
	Vector center; 


	
    // create distance transform to efficiently evaluate distance to nearest edge
    cv::Mat dt;
    cv::distanceTransform(255-img, dt,CV_DIST_L1, 3);

    //TODO: maybe seed random variable for real random numbers.

    uint32_t nIterations = 0; // TODO: adjust this parameter or include some real ransac criteria with inlier/outlier percentages to decide when to stop

    Vector bestCircleCenter;
    float bestCircleRadius;
    float bestCirclePercentage = 0;
	
	//float maxRadius = 100;
    //float minRadius = 50;   // TODO: ADJUST THIS PARAMETER TO YOUR NEEDS, otherwise smaller circles wont be detected or "small noise circles" will have a high percentage of completion

	uint32_t circleIndex = 0;
	// maybe add the +1
	uint32_t offset = 3; //3, 9, 17
	for(uint32_t i=0; i< objects; i++){
		//printf("%d %d %d\n",i , objects,map[i].size);
				

		for(uint32_t j=0; j< map[i].size*1 ; j++){
			
			uint32_t idx1 = rand()%map[i].size;
			uint32_t idx2 = rand()%map[i].size;
			uint32_t idx3 = rand()%map[i].size;

			// we need 3 different samples:
			if(idx1 == idx2) 
				continue;
			if(idx1 == idx3) 
				continue;
			if(idx3 == idx2) 
				continue;
			//printf("%d,%d,%d\n", idx1, idx2,idx3);
			
			// create circle from 3 points:
			getCircle(map[i].pos[idx1],map[i].pos[idx2],map[i].pos[idx3],center,radius);

			// looking for circle bug
			//cv::circle(img, cv::Point(map[i].pos[idx1].x,map[i].pos[idx1].y), 0, cv::Scalar(a,b,c), 2);
			///

			// inlier set unused at the moment but could be used to approximate a (more robust) circle from alle inlier
		   

			//verify or falsify the circle by inlier counting:
			float cPerc = verifyCircle(dt,center,radius,inlierSet);
					
			//
			// it only checks within the object
			//


			//printf("%u %u %f %f\n",center.x,center.y,radius,cPerc);		
			// update best circle information if necessary
			if(radius >= minRadius && radius <= maxRadius){
				//printf("%d %d %f %f\n",center.x,center.y,radius,cPerc);
				if (cPerc > minCirclePercentage && circleIndex < maxCircles){
					//printf("%d\n",circleIndex);
					circles[circleIndex].x = center.x;
					circles[circleIndex].y = center.y;
					circles[circleIndex].radius = radius;
					circles[circleIndex].perc = cPerc;
					circleIndex++;		
					//cv::circle(img, cv::Point(center.x, center.y),radius, cv::Scalar(100,100,0),1);			
				}	
				if(cPerc >= bestCirclePercentage){
					bestCirclePercentage = cPerc;
					bestCircleRadius = radius;
					bestCircleCenter = center;
					//printf("circle my dood");
					//printf("%d %d %f %f\n",bestCircleCenter.x,bestCircleCenter.y,bestCircleRadius,cPerc);		
					//cv::circle(img, cv::Point(center.x, center.y),radius, cv::Scalar(255,100,0),1);	
					//printf("\n percent %f",cPerc);
				}
			}
		}
	}
	


	CircleOccurrence similarCircles[circleIndex];
	uint32_t similarCirclesIndex = 1;
	uint32_t pScale = 6;
	uint32_t rScale = 5;
	similarCircles[0].x = circles[0].x;
	similarCircles[0].y = circles[0].y;
	similarCircles[0].radius = circles[0].radius;
	similarCircles[0].occurrence = 1;
	for (uint32_t i = 1; i < circleIndex; i++){
		for (uint32_t j = 0; j <= similarCirclesIndex; j++){
			//printf("%d,%d %d\n",j, i, similarCirclesIndex);	
			if ((uint32_t)circles[i].x>>pScale == (uint32_t)similarCircles[j].x>>pScale
				&& (uint32_t)circles[i].y>>pScale == (uint32_t)similarCircles[j].y>>pScale
				&& (uint32_t)circles[i].radius>>rScale == (uint32_t)similarCircles[j].radius>>rScale) {
				similarCircles[j].occurrence++;
				break;			
			} else if (j == similarCirclesIndex ){
				similarCirclesIndex++;
				similarCircles[j].occurrence =1;
				similarCircles[j].perc = circles[i].perc;//0;//
				similarCircles[j].x = circles[i].x;//0;//
				similarCircles[j].y = circles[i].y;//0;//
				similarCircles[j].radius = circles[i].radius;//0;//
				//printf("%d, %d, %d, %f %d\n", (uint32_t)similarCircles[j].x>>5, (uint32_t)similarCircles[j].y>>5, ((uint32_t)similarCircles[j].radius)>>4, similarCircles[j].perc, similarCircles[j].occurrence);
				break;
			}

		}
		
		
		//printf("%d, %d, %d, %f, %d\n", (uint32_t)circles[i].x>>5, (uint32_t)circles[i].y>>5, ((uint32_t)circles[i].radius)>>4, circles[i].perc, i);
	}	
	
	uint32_t maxOccurrences = 0;
	CircleOccurrence finalCircles[similarCirclesIndex];
		
	for (uint32_t i = 0; i < similarCirclesIndex; i++){

		//printf("%d, %d, %d, %f %d\n", (uint32_t)similarCircles[i].x>>pScale, (uint32_t)similarCircles[i].y>>pScale, ((uint32_t)similarCircles[i].radius)>>rScale, similarCircles[i].perc, similarCircles[i].occurrence);
		finalCircles[i].perc = 0;//circles[i].perc;//
		finalCircles[i].x = 0;//circles[i].x;//
		finalCircles[i].y = 0;//circles[i].y;//
		finalCircles[i].radius = 0;//circles[i].radius;
		finalCircles[i].occurrence = similarCircles[i].occurrence;
		if (maxOccurrences < finalCircles[i].occurrence)
			maxOccurrences = finalCircles[i].occurrence;
	}

	for (uint32_t i = 0; i < circleIndex; i++){
		for (uint32_t j = 0; j <= similarCirclesIndex; j++){
			//printf("%d,%d %d\n",j, i, similarCirclesIndex);	
			if ((uint32_t)circles[i].x>>pScale == (uint32_t)similarCircles[j].x>>pScale
				&& (uint32_t)circles[i].y>>pScale == (uint32_t)similarCircles[j].y>>pScale
				&& (uint32_t)circles[i].radius>>rScale == (uint32_t)similarCircles[j].radius>>rScale) {
				finalCircles[j].perc += circles[i].perc;
				finalCircles[j].x += circles[i].x;
				finalCircles[j].y += circles[i].y;
				finalCircles[j].radius += circles[i].radius;
				//printf("%d, %d, %d, %f %d\n", (uint32_t)finalCircles[j].x>>5, (uint32_t)finalCircles[j].y>>5, ((uint32_t)finalCircles[j].radius)>>4, finalCircles[j].perc, finalCircles[j].occurrence);
				break;
			}
		}
	}

	//cv::cvtColor(frame, frameCircle, CV_GRAY2RGB);
 	
	//printf("\n\n");


	
	uint32_t mostOccurrence = 0;
	float bestPercentage = 0;
	uint32_t minOccurrence = 12;
		
	for (uint32_t i = 0; i < similarCirclesIndex; i++){
		
		
		finalCircles[i].perc /= finalCircles[i].occurrence; // finalCircles[i].occurrence;//Confidence based on occurence and percentage fit
		finalCircles[i].x /= finalCircles[i].occurrence;
		finalCircles[i].y /= finalCircles[i].occurrence;
		finalCircles[i].radius /= finalCircles[i].occurrence;		
		if (finalCircles[i].perc>minCirclePercentage){
			if (finalCircles[i].occurrence > mostOccurrence && finalCircles[i].occurrence > minOccurrence){
				mostOccurrence = finalCircles[i].occurrence;
				bestCircle.x = finalCircles[i].x;
				bestCircle.y = finalCircles[i].y;
				bestCircle.radius = finalCircles[i].radius;
				bestPercentage = finalCircles[i].perc;
				//cv::circle(img, cv::Point(finalCircles[i].x, finalCircles[i].y),finalCircles[i].radius, cv::Scalar(100,100,0),1);	
				
			}
			//printf("%d, %d, %d, %f %d\n", (uint32_t)finalCircles[i].x, (uint32_t)finalCircles[i].y, ((uint32_t)finalCircles[i].radius), finalCircles[i].perc, finalCircles[i].occurrence);
			//cv::circle(img, cv::Point(finalCircles[i].x, finalCircles[i].y),finalCircles[i].radius, cv::Scalar(100,100,0),1);						
		}
	}
	
	if (mostOccurrence < minOccurrence){
		bestCircle.x = -1;
		bestCircle.y = -1;
		bestCircle.radius = 0;
		printf("No acceptable circle found\nmax occurrences: %d\n ",mostOccurrence);
		return bestCircle;
	}
		
	
	


	
	printf("x: %f pixels \n y: %f pixels\n radius: %f pixels\n percentage: %f\n occurrence: %d times\n",bestCircle.x,bestCircle.y,bestCircle.radius, 100*bestPercentage, mostOccurrence);

    free(inlierSet);
    //cv::imshow("output",img);
    //cv::imshow("mask",mask);
    //cv::waitKey(0);
	//imshow("circles", img );
    return bestCircle;
}


void findEllipse(cv::Mat img, uint16_t minDist, uint16_t maxDist, uint16_t accumThresh){
/* Based on:
/* A New Efficient Ellipse Detection Method */
/* Xie and Ji */
	cv::Mat output;
	uint32_t width = img.cols;
	uint32_t height = img.rows;
	/*
	ellipseInfo* potEllipse;
	EdgePointer* edgeMap; // = (ellipseInfo*) malloc(sizeof(ellipseInfo)); 	
//	edgeMap[0].point = NULL;
	edgeMap[0].x = 0;
	edgeMap[0].y = 0;
	edgeMap[0].remove = 0;
	*/
	uint32_t totalEdge = 0;
	float x0, y0, a; // center coordinates and half length of major axis
	float alpha; // orientation
	float bSquared, cosTau, f; // half length of minor axis approximation, cosTau, and dist
	float tempDist, xDiff, yDiff;
	uint16_t index; // accumulator index
	uint16_t max = 0;
	uint16_t maxIndex;
	
	output = cv::Mat(cv::Size(img.cols,img.rows),CV_8UC3);
	
	for(uint32_t i = 0; i < height; i++){
		for(uint32_t j = 0; j < width; j++){
			if(img.data[j+i*height] == PIXEL_MAX){
//				edgeMap[totalEdge].x = j;
//				edgeMap[totalEdge].y = i;
				totalEdge++;
			}	
		}
	}
	printf("%d\n",totalEdge);
	EdgePointer* edgeMap = (EdgePointer*) malloc(sizeof(EdgePointer)*totalEdge);
	uint16_t l = 0;
	for(uint16_t i = 0; i < height; i++){
		for(uint16_t j = 0; j < width; j++){
			if(img.data[j+i*height] == PIXEL_MAX){
				edgeMap[l].x = j;
				edgeMap[l].y = i;
				edgeMap[l].remove = 0;
				l++;
			}	
		}
	}
	Ellipse potEllipse;
	potEllipse.votes = 0;
	l = 0;
	for(uint16_t i = 0; i < totalEdge; i++){
		/* DECLARE ACCUMULATOR */
		printf("new accum\n");
		uint16_t accumulator[maxDist-minDist];
		if(edgeMap[i].remove) {
			printf("remove 1\n");
			continue;
		}
		for(uint16_t j = 0; j < totalEdge; j++){
			if(edgeMap[j].remove) {
				printf("remove 2\n");
				continue;
			}
			//printf("%d\n",j);
			xDiff = (float)edgeMap[i].x-edgeMap[j].x;
			yDiff = (float)edgeMap[i].y-edgeMap[j].y;
			tempDist = hypot(xDiff,yDiff);
			//printf("%f\t%f\t%f\n",xDiff,yDiff,tempDist);
			if(tempDist > maxDist){
				continue;
			}
			if(tempDist > minDist){
				x0 = ((float)edgeMap[j].x+edgeMap[i].x)/2;
				y0 = ((float)edgeMap[j].y+edgeMap[i].y)/2;
				a = tempDist;
				alpha = atan2(yDiff,xDiff);
				for(uint16_t k = 0; k < totalEdge; k++){	
					xDiff = x0-edgeMap[k].x;
					yDiff = y0-edgeMap[k].y;
					tempDist = hypot(xDiff,yDiff); // maybe replce with function
					if(tempDist > minDist){
						xDiff = (float)edgeMap[k].x - edgeMap[j].x;
						yDiff = (float)edgeMap[k].y - edgeMap[j].y;
						f = hypot(xDiff,yDiff); // maybe replace with function
						cosTau = pow((pow(a,2)+pow(tempDist,2)-pow(f,2))/(2*a*tempDist),2);
						printf("a%f\td%f\tf%f\t%f\n",a,tempDist,f,cosTau);
						bSquared = (pow(a,2)*pow(tempDist,2)*(1-cosTau))/(pow(a,2)-pow(tempDist,2)*cosTau); // WHY IS THIS NEGATIVE????????????????????????????????????????????
						printf("%f\n",bSquared);
						index = (uint16_t)roundf(bSquared-minDist);
						//printf("%d\n",index);
						accumulator[index]++;
						if(accumulator[index] > max){
							max = accumulator[index];
							maxIndex = index;
						}
					}
				}
				if(accumulator[maxIndex] > accumThresh && accumulator[maxIndex] > potEllipse.votes){
//					potEllipse[l].x = x0;
//					potEllipse[l].y = y0;
//					potEllipse[l].a = a;
//					potEllipse[l].alpha = alpha;
//					potEllipse[l].votes = accumulator[maxIndex];
//					l++;
					printf("found ellipse\n");
					potEllipse.x = x0;
					potEllipse.y = y0;
					potEllipse.a = a;
					potEllipse.alpha = alpha;
					potEllipse.votes = accumulator[maxIndex];
					for(uint16_t k = 0; k < totalEdge; k++){
						if(hypot(x0-edgeMap[k].x,y0-edgeMap[k].y) == a){
							edgeMap[k].remove = 1;
						}
					}
					break;
				}
			}
		}
	}
	printf("drawing\n");
	cv::circle(img,cv::Point(potEllipse.x,potEllipse.y),potEllipse.a,cv::Scalar(255,0,0),1);
	cv::namedWindow("with ellipse");
	cv::imshow("with ellipse", img);
}
