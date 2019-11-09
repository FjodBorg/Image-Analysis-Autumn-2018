#ifndef PREPROCESS_INCLUDE
#define PREPROCESS_INCLUDE

#include "opencv2/core.hpp"
#include "opencv2/highgui.hpp"
#include "opencv2/imgproc.hpp"

#include <stdio.h>
#include <math.h>
//#include <curses.h>

#include <iostream>
#include <fstream>

#define XOFFSET 60
#define YOFFSET 500
#define NORMALIZED_SCALE 450
#define PIXEL_MAX 255
#define PIXEL_MIN 0
#define EDGE_WEAK 128
#define EDGE_STRONG 255
#define PI 3.14159265358979323846

typedef struct LabelIndices{
	uint32_t index;
	uint32_t size;
	uint8_t hasStrongEdge;
	LabelIndices *point;
} LabelIndices;

typedef struct Circle{
	float x;
	float y;	
	float radius;	
	float perc; 
} Circle;

typedef struct CircleOccurrence{
	float x;
	float y;	
	float radius;
	float perc;	
	uint32_t occurrence; 
} CircleOccurrence;

typedef struct Vector{
	float x;
	float y;
} Vector;

typedef struct Objects{
	uint32_t size;
	Vector* pos;
} Objects;

typedef struct EdgePointer{
	float x;
	float y;
	uint8_t remove;
} EdgePointer;

typedef struct Ellipse{
    uint32_t x;
    uint32_t y;
    uint32_t a;
    uint32_t alpha;
    uint16_t votes;
} Ellipse;

cv::Mat histogram(const char* histStr, cv::Mat img, int height, int width, int scale);
cv::Mat histogram(const char* histStr, cv::Mat img, int height, int width, uint8_t threshold);
uint16_t histogram(cv::Mat img, float* hist);
void thresholdImage(cv::Mat img, uint8_t threshold);
void thresholdImage(cv::Mat img, uint8_t lowerT, uint8_t upperT);
void getMoments(cv::Mat img, uint32_t moment[6]);
double getAngle(uint32_t moment[6]);
void filter3x3(cv::Mat original, cv::Mat filter, int32_t stencil[3][3]);
float median(float m[], int n);

uint8_t otsu(cv::Mat img);
uint8_t otsuMed(cv::Mat img);
uint8_t percThreshold(cv::Mat img, float desiredPerc);
void cannyDetector(cv::Mat original, cv::Mat output, uint8_t gaussDim, float percent);
void cannyDetector(cv::Mat original, cv::Mat output, uint8_t gaussDim, uint8_t lowerT, uint8_t upperT);
void printImageData(const char* fileLocation, cv::Mat frame);
void gradientMagAngle(cv::Mat input, cv::Mat frameGradMag, cv::Mat frameGradAng);
void nonMaxSuppression(cv::Mat frameGradMag, cv::Mat frameGradAng, cv::Mat output);

uint32_t findMaxLabels(cv::Mat img, uint32_t *labelMap);
void removeObjects(cv::Mat out, uint32_t *labelMap, LabelIndices *labelPointer, uint32_t minObjectSize, uint32_t maxObjectSize);
void applyMap(cv::Mat out, uint32_t *labelMap);
uint32_t labeling(cv::Mat img, Objects **objectMap, uint32_t minObjectSize, uint32_t maxObjectSize);
void findEllipse(cv::Mat img, uint16_t minDist, uint16_t maxDist, uint16_t accumThresh);
float verifyCircle(cv::Mat img,Vector &center, float radius, Vector *inlierSet);
uint8_t dynamicRealloc(Vector *point, uint32_t *prev, uint32_t *current);
void getCircle(Vector &p1,Vector &p2,Vector &p3, Vector &center, float& radiusMinor);
uint32_t getPointPositions(cv::Mat img, Vector *pointPositions);
Circle ransac(cv::Mat img, float maxRadius, float minRadius, float minCirclePercentage);
Circle ransac(cv::Mat img, Objects *map, uint32_t objects, float maxRadius, float minRadius, float minCirclePercentage);
Circle ransacVoting(cv::Mat img, Objects *map, uint32_t objects, float maxRadius, float minRadius, float minCirclePercentage);

#endif
