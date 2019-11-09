#include <conio.h>
#include <stdio.h>
#include "cv.h"
#include "cxcore.h"
#include "highgui.h"
#include "stdlib.h"
#include "math.h"
 
// Offset of histogram
#define XOFFSET 100
#define YOFFSET 450
#define PIXEL_MAX 255
#define PIXEL_MIN 0
#define PI 3.14159265358979323846
 

// readme! 
// pick elements in each corner!

typedef struct cor {
    int x;
    int y;
    int xp;
    int yp;
    int sum;
} cor;

void *test = NULL;
 
 

int findSmallAreaCorrelation(cor* area, IplImage* img1, IplImage* img2, int x, int y,
	int xp, int yp, int Msize, int margin) {
	if (x < img1->width - Msize && y < img1->height - Msize && area != NULL) {
		int sum = 0;
		int temp, a1, a2;

		for (int i = 0; i < Msize; i++) {
			for (int j = 0; j < Msize; j++) {
				a1 = (unsigned char)img1->imageData[(x + j) + (y + i) * img1->width];  // area 1
				a2 = (unsigned char)img2->imageData[(xp + j) + (yp + i) * img2->width];  // area 2
				temp = (int)a1 - int(a2);
				if (temp < 0) {
					temp = -temp;
				}
				sum += temp;
				if (sum > margin) {
					return NULL;
				}
			}
		}
		// why do i need this????
		//area = (cor*)test; 
		//
		area->x = x;
		area->y = y;
		area->xp = xp;
		area->yp = yp;
		area->sum = sum;
		return 1;
	}
	else
		return NULL;
}

int md_comparator(const void *v1, const void *v2)
{
    const cor *p1 = (cor *)v1;
    const cor *p2 = (cor *)v2;
    if (p1->sum < p2->sum)
        return -1;
    else if (p1->sum > p2->sum)
        return +1;
    else
        return 0;
}
 
int main(int argc, char* argv[]) {
    int width = 480;
    int height = 480;
 
    // Names for plots and images
	const char* str1 = "original 1";
	const char* str2 = "original 2";
	const char* str3 = "drawn 1";
	const char* str4 = "drawn 2";
	const char* str5 = "transposed 1";

    IplImage *img1, *img2, *img1Temp, *imgTrans, *img2Temp;
 
    img1 = cvLoadImage("C:\\Users\\30330\\Documents\\Exersizes\\PIC1_R.PNG", 0);
    img2 = cvLoadImage("C:\\Users\\30330\\Documents\\Exersizes\\PIC1_L.PNG", 0);
	img1Temp = cvLoadImage("C:\\Users\\30330\\Documents\\Exersizes\\PIC1_R.PNG", 0); // make a copy
	img2Temp = cvLoadImage("C:\\Users\\30330\\Documents\\Exersizes\\PIC1_L.PNG", 0); // make a copy
	imgTrans = cvLoadImage("C:\\Users\\30330\\Documents\\Exersizes\\PIC1_L.PNG", 0); // make a copy

    int Msize = 10;
    int Mcounter = 0;
    int const max = 50;
	int const maxMatch = 30;
	int const HFtol = 150; //
    int const margin = 600;  // intensity
    int xp = 0;
    int yp = 0;
    int xn[max]; // first pic
    int yn[max]; // fisrt pic
    int tempMin;
    int matchCounter = 0;
	int noError = 0;
 
	
	cor *temp = (cor*)malloc(sizeof(cor*));
	test = temp; 
	cor match[max];
	int a1 = 0;
	// we dont want to see the borders
    for (int y = Msize; y < img1->height - Msize*2 && Mcounter < max; y += Msize*2) {
        for (int x = Msize; x < img1->width - Msize*2 && Mcounter < max; x += Msize*2) {
            xn[Mcounter] = x;
            yn[Mcounter] = y;
			int HF=0, a1;
			for (int i = 0; i < Msize; i++) {
				for (int j = 0; j < Msize; j++) {
					a1 = (unsigned char)
						img1->imageData[(x + j) +
						(y + i) * img1->width]; //area 1
					// alternating to see if same color
					if (i & 1 == true) {

						HF -= a1;

					}					
					else {
						HF += a1;					
					}
				}
			}
			if (HF < 0) {
				HF = -HF;
			}
			if (HF > HFtol) {
				Mcounter++; // try with a new value
			}
        }
    }
    for (int z = 0; z < Mcounter && matchCounter <= maxMatch; z++) {
        tempMin = margin + 1;
		int yp;
		int xp;



        for (yp = 0; yp < img1->height - Msize; yp += 1) {
            for (xp = 0; xp < img1->width - Msize; xp += 1) {

                noError = findSmallAreaCorrelation(temp, 
                    img1, img2, xn[z], yn[z], xp, yp, Msize,
                    margin);  // remeber this is not the center pixel but the
                              // corner one!


                if (noError == NULL) {
                    //printf("\n No match or error :) \n");
                } 
				else if (temp->sum <= tempMin && matchCounter<maxMatch) {
                    if (tempMin >margin) {  // 
											//if first time assigning this value
                        matchCounter++;
                    }
					tempMin = temp->sum;
                    match[matchCounter - 1] = *temp;  // assign match values
                }
            }
        }
		cvDrawRect(img1Temp,
			cvPoint(match[matchCounter - 1].x, match[matchCounter - 1].y),
			cvPoint(match[matchCounter - 1].x + Msize, match[matchCounter - 1].y + Msize),
			cvScalar(0), 2);		
		cvDrawRect(img1Temp,
				cvPoint(match[matchCounter - 1].x, match[matchCounter - 1].y),
				cvPoint(match[matchCounter - 1].x + Msize, match[matchCounter - 1].y + Msize),
				cvScalar(255), 1);
		cvDrawRect(img2Temp,
			cvPoint(match[matchCounter - 1].xp, match[matchCounter - 1].yp),
			cvPoint(match[matchCounter - 1].xp + Msize, match[matchCounter - 1].yp + Msize),
			cvScalar(0), 2);
		cvDrawRect(img2Temp,
			cvPoint(match[matchCounter - 1].xp, match[matchCounter - 1].yp),
			cvPoint(match[matchCounter - 1].xp + Msize, match[matchCounter - 1].yp + Msize),
			cvScalar(255), 1);

    }
 
    qsort(match,matchCounter,sizeof(cor), md_comparator);
 
 
    double p[3];
	double q[3];
	double x1 =  (double) match[0].x;
	double x2 =  (double) match[1].x;
	double x3 =  (double) match[2].x;
	double y1 =  (double) match[0].y;
	double y2 =  (double) match[1].y;
	double y3 =  (double) match[2].y;
	double xp1 = (double) match[0].xp;
	double xp2 = (double) match[1].xp;
	double xp3 = (double) match[2].xp;
	double yp1 = (double) match[0].yp;
	double yp2 = (double) match[1].yp;
	double yp3 = (double) match[2].yp;
    p[0] = -(x1*xp2*y3 - x1*xp3*y2 - x2*xp1*y3 + x2*xp3*y1 + x3*xp1*y2 - x3*xp2*y1)/(x1*y2 - x2*y1 - x1*y3 + x3*y1 + x2*y3 - x3*y2);
    q[0] =  (x1*y2*yp3 - x1*y3*yp2 - x2*y1*yp3 + x2*y3*yp1 + x3*y1*yp2 - x3*y2*yp1)/(x1*y2 - x2*y1 - x1*y3 + x3*y1 + x2*y3 - x3*y2);
    p[1] =  (xp1*y2 - xp2*y1 - xp1*y3 + xp3*y1 + xp2*y3 - xp3*y2)/(x1*y2 - x2*y1 - x1*y3 + x3*y1 + x2*y3 - x3*y2);
    q[1] = -(y1*yp2 - y2*yp1 - y1*yp3 + y3*yp1 + y2*yp3 - y3*yp2)/(x1*y2 - x2*y1 - x1*y3 + x3*y1 + x2*y3 - x3*y2);
    p[2] =  (x1*xp2 - x2*xp1 - x1*xp3 + x3*xp1 + x2*xp3 - x3*xp2)/(x1*y2 - x2*y1 - x1*y3 + x3*y1 + x2*y3 - x3*y2);
    q[2] =  (x1*yp2 - x2*yp1 - x1*yp3 + x3*yp1 + x2*yp3 - x3*yp2)/(x1*y2 - x2*y1 - x1*y3 + x3*y1 + x2*y3 - x3*y2);
 
    for (int y = 0; y < img1->height; y += 1) {
        for (int x = 0; x < img1->width; x += 1) {
            xp = (int)(p[0]+p[1]*x+p[2]*y);
            yp = (int)(q[0]+q[1]*x+q[2]*y);
			if (xp < img1->width && yp < img1->height && xp >= 0 && yp >= 0)
				imgTrans->imageData[x + y * img1Temp->width] = img2->imageData[xp + yp * img1->width] ; // see if it is zero!
            else
                imgTrans->imageData[x + y*img1->width] = (char)255;
        }
    }
 
    cvShowImage(str1, img1);
    cvShowImage(str2, img2);
	cvShowImage(str3, img1Temp);
	cvShowImage(str4, img2Temp);	
	cvShowImage(str5, imgTrans);
 
    while (true) {  // Makes the window stay open
        if (cvWaitKey(5) > 0) break;
    }
 
    cvReleaseImage(&img1);
    cvReleaseImage(&img2);
    cvDestroyWindow(str1);
    cvDestroyWindow(str2);
	free(temp); 
    return 0;
}