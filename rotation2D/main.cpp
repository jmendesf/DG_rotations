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

enum METHOD {
  NEAREST_NEIGHBOR,
  BILINEAR_INTERPOLATION,
  ALL
};

// 2D image definition
typedef ImageContainerBySTLVector<Z2i::Domain, float>  Image;
// Grayscale mapping
typedef GrayscaleColorMap<float> Gray;

// Rigid transformations
typedef ForwardRigidTransformation2D < Z2i::Space > ForwardTrans;
typedef BackwardRigidTransformation2D < Z2i::Space > BackwardTrans;
typedef ConstImageAdapter<Image, Z2i::Domain, BackwardTrans, Image::Value, Identity > MyImageBackwardAdapter;
typedef DomainRigidTransformation2D < Z2i::Domain, ForwardTrans > MyDomainTransformer;
typedef MyDomainTransformer::Bounds Bounds;

// Create a binarizer
typedef functors::IntervalForegroundPredicate<Image> Binarizer;  

// Create and apply a distance transform to the binary image with the L2 norm
// Def L1 norm aswell
typedef DistanceTransformation<Z2i::Space, Binarizer, Z2i::L2Metric> DTL2;
typedef DistanceTransformation<Z2i::Space, Binarizer, Z2i::L1Metric> DTL1;

// Print if wrong args
void usage()
{
  cout << endl;
  cout << "usage: ./rotation <path_to_pgm_file> <angle> <method>" << endl;
  cout << "method: bli (bilinear interpolation), nn (nearest neighbor), all" << endl;
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
Image addImages(DTL1 dtl1Im1, DTL1 dtl1Im2)
{
  // Takes the size of the biggest image
  Image im((dtl1Im1.domain().size() > dtl1Im2.domain().size() ? dtl1Im1.domain() : dtl1Im2.domain()));

  // Max values of DT for each DT
  DTL1::Value maxDT2 = (*std::max_element(dtl1Im1.constRange().begin(), dtl1Im1.constRange().end()));
  DTL1::Value maxDT1 = (*std::max_element(dtl1Im2.constRange().begin(), dtl1Im2.constRange().end()));

  // Compute grayscale coefficient depending on the max DT value
  int step = 255 / maxDT2;
  int step2 = 255 / maxDT1;

  int value, value2;

  for(int y = 0; y < im.domain().upperBound()[1]; ++y)
  {
    for(int x = 0; x < im.domain().upperBound()[0]; ++x)
    {
      value = dtl1Im1.operator()({x,y});
      value2 = dtl1Im2.operator()({x,y});
      
      // the point at x,y takes the maximum value between the two
      im.setValue({x,y}, max(value * step, value2 * step2));
    }
  }
  return im;
}

void imDTToGS(Image& imDT, int minValue, int maxValue)
{

  int step = 255 / (maxValue - minValue);
  float value;

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

void processDT(Image& imDT, bool isInterior)
{
  for(int y = 0; y < imDT.domain().upperBound()[1]; y++)
  {
    for(int x = 0; x < imDT.domain().upperBound()[0]; x++)
    {
      if(imDT.operator()({x,y}) == 1.){
        imDT.setValue({x,y}, isInterior ? -0.5 : 0.5);
      } else
      {
        if(isInterior && imDT.operator()({x,y}) != 0)
          imDT.setValue({x,y}, -imDT.operator()({x,y}));
      }
    }
  }
}

Image rotateBackwardNearestNeighbor(Image image, float angle)
{
  int maxX = image.domain().upperBound()[0];
  int maxY = image.domain().upperBound()[1];
  int minX, minY;
  PointVector<2, int> center(maxX / 2, maxY/ 2);

  int x1, y1;
  int x2, y2;
  int x3, y3;
  int x4, y4;

  angle = -angle;

  x1 = center[0] + (0 - center[0]) * cos(angle) - (0 - center[1]) * sin(angle);
  y1 = center[1] + (0 - center[0]) * sin(angle) + (0 - center[1]) * cos(angle);

  x2 = center[0] + (maxX - center[0]) * cos(angle) - (0 - center[1]) * sin(angle);
  y2 = center[1] + (maxX - center[0]) * sin(angle) + (0 - center[1]) * cos(angle);

  x3 = center[0] + (maxX - center[0]) * cos(angle) - (maxY - center[1]) * sin(angle);
  y3 = center[1] + (maxX - center[0]) * sin(angle) + (maxY - center[1]) * cos(angle);

  x4 = center[0] + (0 - center[0]) * cos(angle) - (maxY - center[1]) * sin(angle);
  y4 = center[1] + (0 - center[0]) * sin(angle) + (maxY - center[1]) * cos(angle);

  minX = min(min(x1,x2), min(x3,x4));
  minY = min(min(y1,y2), min(y3,y4));
  maxX = max(max(x1,x2), max(x3,x4));
  maxY = max(max(y1,y2), max(y3,y4));

  Z2i::Domain domain(Z2i::Point(minX,minY), Z2i::Point(maxX,maxY));

  Image rotIm(domain);

  angle = -angle;

  double backX, backY;

  for(int y = minY; y < rotIm.domain().upperBound()[1]; ++y)
  {
    for(int x = minX; x < rotIm.domain().upperBound()[0]; ++x)
    {
      backX = center[0] + (x - center[0]) * cos(angle) - (y - center[1]) * sin(angle);
      backY = center[1] + (x - center[0]) * sin(angle) + (y - center[1]) * cos(angle);

      backX = round(backX);
      backY = round(backY);


      if((backX >= image.domain().upperBound()[0]) || (backX < 0))
      {
        continue;
      }

      if((backY >= image.domain().upperBound()[1]) || (backY < 0)){
        continue;
      }
        

      rotIm.setValue({x,y}, image.operator()({(int)backX, (int)backY}));
    }
  }
  return rotIm;
}

float computeBilinearInterpolation(Image image, float x, float y)
{
  int x1 = floor(x);
  int x2 = x1 + 1;
  int y1 = floor(y);
  int y2 = y1 + 1;
  
  float Q11 = image.operator()({x1, y1});
  float Q12 = image.operator()({x1, y2});
  float Q21 = image.operator()({x2, y1});
  float Q22 = image.operator()({x2, y2});


  float res = (1/((x2 - x1) * (y2 - y1))) 
        * (Q11 * (x2 - x) * (y2 - y) 
        + Q21 * (x - x1) * (y2 - y) 
        + Q12 * (x2 - x) * (y- y1) 
        + Q22 * (x - x1) * (y - y1));

  return res;
}

Image rotateBackwardBilinearInterpolation(Image image, float angle)
{
  int maxX = image.domain().upperBound()[0];
  int maxY = image.domain().upperBound()[1];
  int minX, minY;
  PointVector<2, int> center(maxX / 2, maxY/ 2);

  int x1, y1;
  int x2, y2;
  int x3, y3;
  int x4, y4;

  angle = -angle;

  x1 = center[0] + (0 - center[0]) * cos(angle) - (0 - center[1]) * sin(angle);
  y1 = center[1] + (0 - center[0]) * sin(angle) + (0 - center[1]) * cos(angle);

  x2 = center[0] + (maxX - center[0]) * cos(angle) - (0 - center[1]) * sin(angle);
  y2 = center[1] + (maxX - center[0]) * sin(angle) + (0 - center[1]) * cos(angle);

  x3 = center[0] + (maxX - center[0]) * cos(angle) - (maxY - center[1]) * sin(angle);
  y3 = center[1] + (maxX - center[0]) * sin(angle) + (maxY - center[1]) * cos(angle);

  x4 = center[0] + (0 - center[0]) * cos(angle) - (maxY - center[1]) * sin(angle);
  y4 = center[1] + (0 - center[0]) * sin(angle) + (maxY - center[1]) * cos(angle);

  minX = min(min(x1,x2), min(x3,x4));
  minY = min(min(y1,y2), min(y3,y4));
  maxX = max(max(x1,x2), max(x3,x4));
  maxY = max(max(y1,y2), max(y3,y4));



  Z2i::Domain domain(Z2i::Point(minX,minY), Z2i::Point(maxX,maxY));
  Image rotIm(domain);

  angle = -angle;

  double backX, backY;

  for(int y = minY; y < rotIm.domain().upperBound()[1]; ++y)
  {
    for(int x = minX; x < rotIm.domain().upperBound()[0]; ++x)
    {
      backX = center[0] + (x - center[0]) * cos(angle) - (y - center[1]) * sin(angle);
      backY = center[1] + (x - center[0]) * sin(angle) + (y - center[1]) * cos(angle);

      if((backX > image.domain().upperBound()[0]) || (backX < 0))
        continue;
      if((backY > image.domain().upperBound()[1]) || (backY < 0))
        continue;

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

void processImage(Image& image, float angle, METHOD method)
{
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
  Binarizer b(image, 1, 135);
  Binarizer bInv(imInv, 1, 135);

  // DTL2 dtl2(&image.domain(), &b, &Z2i::l2Metric);
  DTL1 dtl1(&image.domain(), &b, &Z2i::l1Metric);
  DTL1 dtl1Inv(&imInv.domain(), &bInv, &Z2i::l1Metric);

  // Compute the max value of the DT
  DTL1::Value maxDT2 = (*std::max_element(dtl1.constRange().begin(), dtl1.constRange().end()));
  DTL1::Value maxDT1 = (*std::max_element(dtl1Inv.constRange().begin(), dtl1Inv.constRange().end()));

  // Create GS DT images
  Image imGS = createImageFromDT(dtl1, maxDT2, false);
  Image imInvGS = createImageFromDT(dtl1Inv, maxDT1, false);

  processDT(imGS, true);
  processDT(imInvGS, false);

  // Add both DT
  Image imAddDTL = addImages(imInvGS, imGS);

  if((method == NEAREST_NEIGHBOR) || (method == ALL))
  {
    cout << "Computing rotation using nearest neighbor -" << endl;
    Image rotIm = rotateBackwardNearestNeighbor(imAddDTL, angle);
    thresholdDTImage(rotIm, rotIm);
    saveImage(board, rotIm, 0, 255, "../output/rot_NN");
  }
  
  if((method == BILINEAR_INTERPOLATION) || (method == ALL))
  {
    cout << "Computing rotation using bilinear interpolation -" << endl;
    Image rotIm = rotateBackwardBilinearInterpolation(imAddDTL, angle);
    thresholdDTImage(rotIm, rotIm);
    saveImage(board, rotIm, 0, 255, "../output/rot_BLI");
  }

  imDTToGS(imAddDTL, -maxDT2, maxDT1);

  // Create a set for the image and its inverse
  Z2i::DigitalSet set(image.domain());
  Z2i::DigitalSet setInv(imInv.domain());
  DigitalSetInserter<Z2i::DigitalSet> inserter(set);
  DigitalSetInserter<Z2i::DigitalSet> inserterInv(setInv);
  DGtal::setFromImage(dtl1, inserter, 1, 135);
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
  // Check args
  if(argc != 4)
  {
    usage();
    return 0;
  }

  // Load image
  Image image = PGMReader<Image>::importPGM(argv[1]);
  
  if(strcmp(argv[3], "bli") == 0)
    processImage(image, -stof(argv[2]), BILINEAR_INTERPOLATION);
  else if(strcmp(argv[3],"nn") == 0)
    processImage(image, -stof(argv[2]), NEAREST_NEIGHBOR);
  else if(strcmp(argv[3], "all") == 0)
    processImage(image, -stof(argv[2]), ALL);
  else 
    usage();

  return 0;
}
