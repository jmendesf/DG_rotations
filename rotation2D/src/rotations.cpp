#include "../include/rotations.h"
#include "../include/tools.h"

// assigns the pixel value corresponding to pixel with clamped coordinates (x,y) to val
void getClampedPixelValue(Image image, int x, int y, float& val)
{
  int upperBoundX = image.domain().upperBound()[0];
  int upperBoundY = image.domain().upperBound()[1];

  x = clampInt(x, 0, upperBoundX - 1);
  y = clampInt(y, 0, upperBoundY - 1);

  val = image.operator()({x,y});
}

// Bilinear interpolation function
// x, y the coordinates calculated by backward rotation
float computeBilinearInterpolation(Image image, float x, float y)
{
  // compute floor coordinates
  int x1 = floor(x);
  int x2 = x1 + 1;
  int y1 = floor(y);
  int y2 = y1 + 1;

  // values of the 4 adjacent pixels
  float Q11 = image.operator()({x1, y1});
  float Q12 = image.operator()({x1, y2});
  float Q21 = image.operator()({x2, y1});
  float Q22 = image.operator()({x2, y2});

  // interpolation computation
  float res = (1/((x2 - x1) * (y2 - y1))) 
        * (Q11 * (x2 - x) * (y2 - y) 
        + Q21 * (x - x1) * (y2 - y) 
        + Q12 * (x2 - x) * (y- y1) 
        + Q22 * (x - x1) * (y - y1));

  return res;
}

// computes the cubic hermic interpolation 
float cubicHermite(float A, float B, float C, float D, float t)
{
  float a = -A / 2.0f + (3.0f*B) / 2.0f - (3.0f*C) / 2.0f + D / 2.0f;
  float b = A - (5.0f*B) / 2.0f + 2.0f*C - D / 2.0f;
  float c = -A / 2.0f + C / 2.0f;
  float d = B;

  return a*t*t*t + b*t*t + c*t + d;
}

// computes the bicubic interpolation of a pixel 
float computeBicubicInterpolation(Image image, float x, float y)
{
  int xInt = (int) x;
  float xFract = x - floor(x);

  int yInt = (int) y;
  float yFract = y - floor(y);

  // 16 pixel values
  float p00, p10, p20, p30;
  float p01, p11, p21, p31;
  float p02, p12, p22, p32;
  float p03, p13, p23, p33;

  // clamp to avoid having out of bounds pixels
  // first row
  getClampedPixelValue(image, xInt - 1, yInt - 1, p00);
  getClampedPixelValue(image, xInt + 0, yInt - 1, p10);
  getClampedPixelValue(image, xInt + 1, yInt - 1, p20);
  getClampedPixelValue(image, xInt + 2, yInt - 1, p30);

  // second row
  getClampedPixelValue(image, xInt - 1, yInt, p01);
  getClampedPixelValue(image, xInt + 0, yInt, p11);
  getClampedPixelValue(image, xInt + 1, yInt, p21);
  getClampedPixelValue(image, xInt + 2, yInt, p31);

  // third row
  getClampedPixelValue(image, xInt - 1, yInt + 1, p02);
  getClampedPixelValue(image, xInt + 0, yInt + 1, p12);
  getClampedPixelValue(image, xInt + 1, yInt + 1, p22);
  getClampedPixelValue(image, xInt + 2, yInt + 1, p32);

  // fourth row
  getClampedPixelValue(image, xInt - 1, yInt + 2, p03);
  getClampedPixelValue(image, xInt + 0, yInt + 2, p13);
  getClampedPixelValue(image, xInt + 1, yInt + 2, p23);
  getClampedPixelValue(image, xInt + 2, yInt + 2, p33);

  // col interpolation
  float col0 = cubicHermite(p00, p10, p20, p30, xFract);
  float col1 = cubicHermite(p01, p11, p21, p31, xFract);
  float col2 = cubicHermite(p02, p12, p22, p32, xFract);
  float col3 = cubicHermite(p03, p13, p23, p33, xFract);

  // return the complete interpolation
  return cubicHermite(col0, col1, col2, col3, yFract);
}

// Compute rotation with nearest neighbor choice of pixel value
Image rotateBackward(Image image, float angle, INTERPOLATION_METHOD method)
{
  // domain's extrema
  int maxX = image.domain().upperBound()[0];
  int maxY = image.domain().upperBound()[1];
  int minX, minY;

  // Compute rotation point
  PointVector<2, int> center(maxX / 2, maxY/ 2);
    
  int x1, y1;
  int x2, y2;
  int x3, y3;
  int x4, y4;  

  // Rotation for extremas will be done with forward rotation
  angle = -angle;

  // Compute position of the 4 corners of the domain
  x1 = center[0] + (0 - center[0]) * cos(angle) - (0 - center[1]) * sin(angle);
  y1 = center[1] + (0 - center[0]) * sin(angle) + (0 - center[1]) * cos(angle);

  x2 = center[0] + (maxX - center[0]) * cos(angle) - (0 - center[1]) * sin(angle);
  y2 = center[1] + (maxX - center[0]) * sin(angle) + (0 - center[1]) * cos(angle);

  x3 = center[0] + (maxX - center[0]) * cos(angle) - (maxY - center[1]) * sin(angle);
  y3 = center[1] + (maxX - center[0]) * sin(angle) + (maxY - center[1]) * cos(angle);

  x4 = center[0] + (0 - center[0]) * cos(angle) - (maxY - center[1]) * sin(angle);
  y4 = center[1] + (0 - center[0]) * sin(angle) + (maxY - center[1]) * cos(angle);

  // Compute extrema values
  minX = std::min(std::min(x1,x2), std::min(x3,x4));
  minY = std::min(std::min(y1,y2), std::min(y3,y4));
  maxX = std::max(std::max(x1,x2), std::max(x3,x4));
  maxY = std::max(std::max(y1,y2), std::max(y3,y4));

  // Create the corresponding domain and image
  Z2i::Domain domain(Z2i::Point(minX,minY), Z2i::Point(maxX,maxY));
  Image rotIm(domain);

  angle = -angle;

  float backX, backY;

  for(int y = minY; y < rotIm.domain().upperBound()[1]; ++y)
  {
    for(int x = minX; x < rotIm.domain().upperBound()[0]; ++x)
    {
      // Compute backward rotation
      backX = center[0] + (x - center[0]) * cos(angle) - (y - center[1]) * sin(angle);
      backY = center[1] + (x - center[0]) * sin(angle) + (y - center[1]) * cos(angle);

      // Rounding: Nearest neighbor
      if(method == NEAREST_NEIGHBOR)
      {
        backX = round(backX);
        backY = round(backY);
      }
      
      // Ensure position is valid
      if((backX >= image.domain().upperBound()[0]) || (backX < 0))
        continue;

      if((backY >= image.domain().upperBound()[1]) || (backY < 0))
        continue;  

      // Set the value in the rotated Image
      // Nearest neighbor
      if(method == NEAREST_NEIGHBOR)
        rotIm.setValue({x,y}, image.operator()({(int)backX, (int)backY}));
      else if(method == BILINEAR_INTERPOLATION)
        rotIm.setValue({x,y}, computeBilinearInterpolation(image, backX, backY));
      else if(method == BICUBIC_INTERPOLATION)
        rotIm.setValue({x,y}, computeBicubicInterpolation(image, backX, backY));
    }
  }
  return rotIm;
}