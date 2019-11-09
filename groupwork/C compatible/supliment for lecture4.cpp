#include <math.h>
#include <stdio.h>
#include <algorithm>
#include <array>
#include <iostream>

void filter(unsigned char *img, unsigned char *filtered, int32_t filter[3][3]) {
   uint8_t divider = 0;
   int32_t temp;
   uint32_t height = 7;
   uint32_t width = 7;
   uint32_t index = 0;
   int8_t value = 0;

   for (int k = 0; k < 3; k++) {
      for (int l = 0; l < 3; l++) {
         divider += filter[k][l];
         printf("%2d ", filter[k][l]);
      }
      printf("\n");
   }
   if (divider == 0) {
      divider = 1;
   }

   for (uint32_t i = 1; height - 1 > i; i++) {
      for (uint32_t j = 1; width - 1 > j; j++) {
         temp = 0;
         for (int k = 0; k < 3; k++) {
            for (int l = 0; l < 3; l++) {
               index = j + (l - 1) + (i + (k - 1)) * width;
               value = img[index];
               temp += (int32_t)(value) * (filter[k][l]);
            }
         }
         if (temp < 0) temp = 0;
         filtered[j + i * width] = (uint8_t)(temp / divider);
         printf("%3d  ", (uint8_t)(temp / divider));
      }
      printf("\n");
   }
}

// unsigned char originalPicture[49] = {
//     1, 2, 3,  4,  5,  6, 7, 1, 5,  3, 4, 51, 6, 7,  1, 2, 3,
//     4, 5, 56, 7,  12, 2, 3, 4, 65, 6, 7, 1,  9, 38, 4, 5, 26,
//     7, 1, 22, 35, 4,  5, 6, 7, 1,  2, 3, 4,  5, 6,  7};
    unsigned char originalPicture[49] = {
    0,0,0,0,0,0,0, 
    5,8,8,8,8,8,5,  
    0,5,0,0,0,0,5,  
    0,5,5,5,5,20,5,  
    0,0,0,6,6,0,0,  
    0,0,0,6,5,0,0, 
    0,0,0,5,6,0,0};
unsigned char filteredPicture[49];

int main(void) {
   float percentage = 0.5;
   int height = 7;
   int width = 7;

   for (int i = 0; i < height ; i++) {
      for (int j = 0; j < width ; j++) {
         printf("%3d  ", originalPicture[j + i * width]);
      }
      printf("\n");
   }
printf("\n\n\n\n\n");
   unsigned char window[9];
   // std::array<unsigned char, 9> window;
   int tmp = ceil((9 - 1) * percentage);  // 0 indexed
   for (int i = 1; i < height - 1; i++) {
      for (int j = 1; j < width - 1; j++) {
         for (int y = -1; y <= 1; y++) {
            for (int x = -1; x <= 1; x++) {
               window[(y + 1) * 3 + x + 1] =
                   originalPicture[(i + y) * width + j + x];
            }
         }
         std::sort(window, window + 9);
         filteredPicture[i * width + j] = window[tmp];
         printf("%3d  ", window[tmp]);
      }
      printf("\n");
   }
   printf("\n\n\n\n");
   int32_t kernel[3][3] = {{-1, -2, -1}, {0, 0, 0}, {1, 2, 1}};
   filter(originalPicture, filteredPicture, kernel);
   printf("\n\n\n\n");
   int32_t kernel2[3][3] = {{-4, 0, 0}, {0, -4, 0}, {12, 0, -4}};
   filter(originalPicture, filteredPicture, kernel2);
   printf("ree");
   return 0;
}