#include <iostream>
#include <fstream>
#include <algorithm>
#include <string>
#include <boost/algorithm/minmax_element.hpp>
#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/base/BasicFunctors.h"
#include "DGtal/kernel/BasicPointPredicates.h"
#include "DGtal/kernel/sets/DigitalSetInserter.h"
#include "DGtal/images/ImageContainerBySTLVector.h"
#include "DGtal/images/ImageHelper.h"
#include "DGtal/geometry/volumes/distance/DistanceTransformation.h"
#include "DGtal/images/IntervalForegroundPredicate.h"
#include "DGtal/io/boards/Board2D.h"
#include "DGtal/io/readers/PGMReader.h"
#include "DGtal/io/colormaps/HueShadeColorMap.h"
#include "DGtal/io/colormaps/GrayscaleColorMap.h"
#include "DGtal/images/RigidTransformation2D.h"
#include "Board/Point.h"

using namespace std;
using namespace DGtal;
using namespace functors;

// Enum on the interpolation method
enum INTERPOLATION_METHOD {
  NEAREST_NEIGHBOR,
  BILINEAR_INTERPOLATION,
  BICUBIC_INTERPOLATION,
  ALL
};

// 2D image definition
typedef ImageContainerBySTLVector<Z2i::Domain, float>  Image;
// Grayscale mapping
typedef GrayscaleColorMap<float> Gray;

// Create a binarizer
typedef functors::IntervalForegroundPredicate<Image> Binarizer;  

// Create and apply a distance transform to the binary image with the L2 norm
// Def L1 norm aswell
typedef DistanceTransformation<Z2i::Space, Binarizer, Z2i::L2Metric> DTL2;
typedef DistanceTransformation<Z2i::Space, Binarizer, Z2i::L1Metric> DTL1;

// Print if wrong args
void usage(bool help)
{
  cout << endl;
  cout << "usage: ./rotation2D <path_to_pgm_file> <angle> <method> <minThresh> <maxThresh>" << endl;
  cout << "method: nn (nearest neighbor), bli (bilinear), bic (bicubic), all" << endl;
  
  if(help)
  {
    cout << endl << "Threshold values:" << endl;
    cout << "ContourS.pgm: 1, 135" << endl;
    cout << "key.pgm: 150, 255" << endl << endl;
    cout << "example: ./rotation2D ../samples/contourS.pgm 2.5 all 1 135" << endl;
  } else
  {
    cout << "./rotation2D help for more" << endl;
  }
  
  cout << endl;
}

// Inverse the given image
void inverseImage(Image& image)
{
  for(int y = 0; y < image.domain().myUpperBound[1]; y++)
    for(int x = 0; x < image.domain().myUpperBound[0]; x++)
      image.setValue({x,y}, 255 - int(image.operator()({x,y})));
}

// Returns an image corresponding to the given DT
Image createImageFromDT(DTL1 dtl1, int maxValue, bool toGS)
{
  Image im(dtl1.domain());
  int step;
  if(toGS)
    step = maxValue != 0 ? 255 / maxValue : 0;

  int value;

  for(int y = 0; y < im.domain().upperBound()[1]; ++y)
    for(int x = 0; x < im.domain().upperBound()[0]; ++x)
    {
      value = dtl1.operator()({x,y});
      im.setValue({x,y}, toGS ? step * value : value);
    }
  return im;
}

Image createImageFromDT(DTL2 dtl2, int maxValue, bool toGS)
{
  Image im(dtl2.domain());
  int step;
  if(toGS)
    step = maxValue != 0 ? 255 / maxValue : 0;

  float value = 0;

  for(int y = 0; y < im.domain().upperBound()[1]; ++y)
    for(int x = 0; x < im.domain().upperBound()[0]; ++x)
    {
      value = dtl2.operator()({x,y});
      im.setValue({x,y}, toGS ? step * value : value);
    }
  return im;
}

// Returns the addition of two images
Image addImages(Image im1, Image im2)
{
  // Takes the size of the biggest image
  Image im((im1.domain().size() > im2.domain().size() ? im1.domain() : im2.domain()));
  float value;

  for(int y = 0; y < im.domain().upperBound()[1]; ++y)
    for(int x = 0; x < im.domain().upperBound()[0]; ++x)
    {
      value = im1.operator()({x,y}) + im2.operator()({x,y});
      if(value == 0.)
      {
        value = -0.5;
      }

      // image value cannot be more than 255
      if(value > 255.)
        value = 255.;

      im.setValue({x, y}, value);
    }
  return im;
}

// Returns the addition of two DT in the form of an image
Image addImages(DTL2 dtl2Im1, DTL2 dtl2Im2)
{
  // Takes the size of the biggest image
  Image im((dtl2Im1.domain().size() > dtl2Im2.domain().size() ? dtl2Im1.domain() : dtl2Im2.domain()));

  // Max values of DT for each DT
  DTL1::Value maxDT2 = (*std::max_element(dtl2Im1.constRange().begin(), dtl2Im1.constRange().end()));
  DTL1::Value maxDT1 = (*std::max_element(dtl2Im2.constRange().begin(), dtl2Im2.constRange().end()));

  // Compute grayscale coefficient depending on the max DT value
  int step = maxDT2 != 0 ? 255 / maxDT2 : 1;
  int step2 = maxDT1 != 0 ? 255 / maxDT1 : 1;

  int value, value2;

  for(int y = 0; y < im.domain().upperBound()[1]; ++y)
  {
    for(int x = 0; x < im.domain().upperBound()[0]; ++x)
    {
      value = dtl2Im1.operator()({x,y});
      value2 = dtl2Im2.operator()({x,y});
      
      // the point at x,y takes the maximum value between the two
      im.setValue({x,y}, max(value * step, value2 * step2));
    }
  }
  return im;
}

// Converts a DT image into its grayscale equivalent (values from 0 to 255)
void imDTToGS(Image& imDT, int minValue, int maxValue)
{
  // avoid division by 0
  int step = (maxValue != minValue) ? 255 / (maxValue - minValue) : 1;
  float value;

  // compute equivalent for DT in GS
  for(int y = 0; y < imDT.domain().upperBound()[1]; ++y)
  {
    for(int x = 0; x < imDT.domain().upperBound()[0]; ++x)
    {
      value = imDT.operator()({x,y});

      if(imDT.operator()({x,y}) < 0)
        value = 255;

      imDT.setValue({x,y}, abs(step * value));
    }
  }
}

// if pixel is positive: 0 (background)
// if pixel is negative: 255 (foreground)
void thresholdDTImage(Image src, Image& dst)
{
  for(int y = src.domain().lowerBound()[1]; y < src.domain().upperBound()[1]; y++)
  {
    for(int x = src.domain().lowerBound()[0]; x < src.domain().upperBound()[0]; x++)
    {
      dst.setValue({x,y}, src.operator()({x,y}) < 0 ? 255 : 0);
    }
  }
}

// Prepares the DT for rotation. 
// Contour pixels take the value of 0.5 if belonging to the BG, -0.5 if FG.
void processDT(Image& imDT, bool isInterior)
{
  for(int y = 0; y < imDT.domain().upperBound()[1]; y++)
  {
    for(int x = 0; x < imDT.domain().upperBound()[0]; x++)
    {
      // if(imDT.operator()({x,y}) == 1.){
      //   imDT.setValue({x,y}, isInterior ? -0.5 : 0.5);
      // } else
      // {
      if(isInterior && (imDT.operator()({x,y}) != 0))
        imDT.setValue({x,y}, -imDT.operator()({x,y}));
      // }
    }
  }
}


// Compute rotation with nearest neighbor choice of pixel value
Image rotateBackwardNearestNeighbor(Image image, float angle)
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
  minX = min(min(x1,x2), min(x3,x4));
  minY = min(min(y1,y2), min(y3,y4));
  maxX = max(max(x1,x2), max(x3,x4));
  maxY = max(max(y1,y2), max(y3,y4));

  // Create the corresponding domain and image
  Z2i::Domain domain(Z2i::Point(minX,minY), Z2i::Point(maxX,maxY));
  Image rotIm(domain);

  angle = -angle;

  double backX, backY;

  for(int y = minY; y < rotIm.domain().upperBound()[1]; ++y)
  {
    for(int x = minX; x < rotIm.domain().upperBound()[0]; ++x)
    {
      // Compute backward rotation
      backX = center[0] + (x - center[0]) * cos(angle) - (y - center[1]) * sin(angle);
      backY = center[1] + (x - center[0]) * sin(angle) + (y - center[1]) * cos(angle);

      // Rounding: Nearest neighbor
      backX = round(backX);
      backY = round(backY);

      // Ensure position is valid
      if((backX >= image.domain().upperBound()[0]) || (backX < 0))
        continue;

      if((backY >= image.domain().upperBound()[1]) || (backY < 0))
        continue;  

      // Set the value in the rotated Image
      rotIm.setValue({x,y}, image.operator()({(int)backX, (int)backY}));
    }
  }
  return rotIm;
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

// clamps an integer value
int clampInt(int value, int low, int high)
{
  return value < low ? low : value > high ? high : value; 
}

// assigns the pixel value corresponding to pixel with clamped coordinates (x,y) to val
void getClampedPixelValue(Image image, int x, int y, float& val)
{
  int upperBoundX = image.domain().upperBound()[0];
  int upperBoundY = image.domain().upperBound()[1];

  x = clampInt(x, 0, upperBoundX - 1);
  y = clampInt(y, 0, upperBoundY - 1);

  val = image.operator()({x,y});
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

// Compute the backward rotation of all pixels of an image using bicubic interpolation
Image rotateBackwardBicubicInterpolation(Image image, float angle)
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
  minX = min(min(x1,x2), min(x3,x4));
  minY = min(min(y1,y2), min(y3,y4));
  maxX = max(max(x1,x2), max(x3,x4));
  maxY = max(max(y1,y2), max(y3,y4));

  // Create the corresponding domain and image
  Z2i::Domain domain(Z2i::Point(minX,minY), Z2i::Point(maxX,maxY));
  Image rotIm(domain);

  angle = -angle;

  double backX, backY;

  for(int y = minY; y < rotIm.domain().upperBound()[1]; ++y)
  {
    for(int x = minX; x < rotIm.domain().upperBound()[0]; ++x)
    {
      // Compute backward rotation
      backX = center[0] + (x - center[0]) * cos(angle) - (y - center[1]) * sin(angle);
      backY = center[1] + (x - center[0]) * sin(angle) + (y - center[1]) * cos(angle);

      // Ensure position is valid
      if((backX > image.domain().upperBound()[0]) || (backX < 0))
        continue;
      if((backY > image.domain().upperBound()[1]) || (backY < 0))
        continue;

      // set value to the value of the bicubic interpolation
      rotIm.setValue({x,y}, computeBicubicInterpolation(image, backX, backY));
    }
  }

  return rotIm;
}

Image rotateBackwardBilinearInterpolation(Image image, float angle)
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
  minX = min(min(x1,x2), min(x3,x4));
  minY = min(min(y1,y2), min(y3,y4));
  maxX = max(max(x1,x2), max(x3,x4));
  maxY = max(max(y1,y2), max(y3,y4));

  // Create the corresponding domain and image
  Z2i::Domain domain(Z2i::Point(minX,minY), Z2i::Point(maxX,maxY));
  Image rotIm(domain);

  angle = -angle;

  double backX, backY;

  for(int y = minY; y < rotIm.domain().upperBound()[1]; ++y)
  {
    for(int x = minX; x < rotIm.domain().upperBound()[0]; ++x)
    {
      // Compute backward rotation
      backX = center[0] + (x - center[0]) * cos(angle) - (y - center[1]) * sin(angle);
      backY = center[1] + (x - center[0]) * sin(angle) + (y - center[1]) * cos(angle);

      // Ensure position is valid
      if((backX > image.domain().upperBound()[0]) || (backX < 0))
        continue;
      if((backY > image.domain().upperBound()[1]) || (backY < 0))
        continue;

      // Set the interpolated value 
      rotIm.setValue({x,y}, computeBilinearInterpolation(image, backX, backY));
    }
  }
  return rotIm;
}

// Save board from image (
// Save to path
void saveImage(Board2D board, Image image, int minVal, int maxVal, string path)
{
  board.clear();
  Display2DFactory::drawImage<Gray>(board, image, minVal, maxVal);
  board.saveEPS(path.c_str());
}

// Save board from set
// Save to path
void saveSet(Board2D board, Z2i::DigitalSet set, string path)
{
  board.clear();
  board << set.domain() << set;
  board.saveEPS(path.c_str());
}

void processImage(Image& image, float angle, INTERPOLATION_METHOD method, int minThresh, int maxThresh)
{
  cout << endl;

  // Create board and set its size
  Board2D board;
  board << image.domain();
  board.clear();

  // copy the image then inverse it 
  Image imInv = image;
  inverseImage(imInv);

  // Import image to the board
  Display2DFactory::drawImage<Gray>(board, image, (unsigned char)0, (unsigned char)255);

  // Threshold the image
  Binarizer b(image, minThresh, maxThresh);
  Binarizer bInv(imInv, minThresh, maxThresh);

  // DTL2 dtl2(&image.domain(), &b, &Z2i::l2Metric);
  DTL2 dtl2(&image.domain(), &b, &Z2i::l2Metric);
  DTL2 dtl2Inv(&imInv.domain(), &bInv, &Z2i::l2Metric);  

  // Compute the max value of the DT
  DTL2::Value maxDT2 = (*std::max_element(dtl2.constRange().begin(), dtl2.constRange().end()));
  DTL2::Value maxDT1 = (*std::max_element(dtl2Inv.constRange().begin(), dtl2Inv.constRange().end()));

  // Create GS DT images
  Image imGS = createImageFromDT(dtl2, maxDT2, false);
  Image imInvGS = createImageFromDT(dtl2Inv, maxDT1, false);

  processDT(imGS, true);
  processDT(imInvGS, false);

  // Add both DT
  Image imAddDTL = addImages(imInvGS, imGS);

  // Compute rotations
  // Nearest neighbor
  if((method == NEAREST_NEIGHBOR) || (method == ALL))
  {
    cout << "Computing rotation using nearest neighbor -";
    Image rotIm = rotateBackwardNearestNeighbor(imAddDTL, angle);
    thresholdDTImage(rotIm, rotIm);
    saveImage(board, rotIm, 0, 255, "../output/rot_NN.eps");
    cout << " done." << endl;
    cout << "Output save as ../output/rot_NN.eps" << endl << endl;
  }
  
  // Bilinear interpolation
  if((method == BILINEAR_INTERPOLATION) || (method == ALL))
  {
    cout << "Computing rotation using bilinear interpolation -";
    Image rotIm = rotateBackwardBilinearInterpolation(imAddDTL, angle);
    thresholdDTImage(rotIm, rotIm);
    saveImage(board, rotIm, 0, 255, "../output/rot_BLI.eps");
    cout << " done." << endl;
    cout << "Output save as ../output/rot_BLI.eps" << endl << endl;
  }

  if((method == BICUBIC_INTERPOLATION) || (method == ALL))
  {
    cout << "Computing rotation using bicubic interpolation -";
    Image rotIm = rotateBackwardBicubicInterpolation(imAddDTL, angle);
    thresholdDTImage(rotIm, rotIm);
    saveImage(board, rotIm, 0, 255, "../output/rot_BIC.eps");
    cout << " done." << endl;
    cout << "Output save as ../output/rot_BIC.eps" << endl << endl;
  }

  imDTToGS(imAddDTL, -maxDT2, maxDT1);

  // Create a set for the image and its inverse
  Z2i::DigitalSet set(image.domain());
  Z2i::DigitalSet setInv(imInv.domain());
  DigitalSetInserter<Z2i::DigitalSet> inserter(set);
  DigitalSetInserter<Z2i::DigitalSet> inserterInv(setInv);
  DGtal::setFromImage(dtl2, inserter, 1, 135);
  DGtal::setFromImage(imInv, inserterInv, 1, 135);

  // Save output 
  saveImage(board, imAddDTL, 0, 255, "../output/im_add_DT.eps");
  saveImage(board, imGS, 0, 255, "../output/im_GS_DT.eps");
  saveImage(board, imInvGS, 0, 255, "../output/im_inv_GS_DT");
  saveSet(board, set, "../output/set.eps");
  saveSet(board, setInv, "../output/set_inv.eps");

  cout << endl;
}

int main(int argc, char** argv)
{
  if(argc == 2)
    if(strcmp(argv[1], "help") == 0)
    {
      usage(true);
      return 0;
    }

  // Check args
  if(argc != 6)
  {
    usage(false);
    return 0;
  }

  // Load image
  Image image = PGMReader<Image>::importPGM(argv[1]);
  
  // process depending on the user's choice
  if(strcmp(argv[3], "bli") == 0)
    processImage(image, -stof(argv[2]), BILINEAR_INTERPOLATION, stoi(argv[4]), stoi(argv[5]));
  else if(strcmp(argv[3],"nn") == 0)
    processImage(image, -stof(argv[2]), NEAREST_NEIGHBOR, stoi(argv[4]), stoi(argv[5]));
  else if(strcmp(argv[3], "bic") == 0)
    processImage(image, -stof(argv[2]), BICUBIC_INTERPOLATION, stoi(argv[4]), stoi(argv[5]));
  else if(strcmp(argv[3], "all") == 0)
    processImage(image, -stof(argv[2]), ALL, stoi(argv[4]), stoi(argv[5]));
  else 
    usage(false); 

  return 0;
}
