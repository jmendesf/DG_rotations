#include <iostream>
#include <fstream>
#include <algorithm>
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

using namespace std;
using namespace DGtal;

// 2D image definition
typedef ImageContainerBySTLVector<Z2i::Domain, int>  Image;
// Grayscale mapping
typedef GrayscaleColorMap<int> Gray;

// Create a binarizer
typedef functors::IntervalForegroundPredicate<Image> Binarizer;  

// Create and apply a distance transform to the binary image with the L2 norm
// Def L1 norm aswell
typedef DistanceTransformation<Z2i::Space, Binarizer, Z2i::L2Metric> DTL2;
typedef DistanceTransformation<Z2i::Space, Binarizer, Z2i::L1Metric> DTL1;

// Print if wrong args
void usage()
{
  cout << "usage: ./rotation <path_to_pgm_file>" << endl;
}

// Inverse the given image
void inverseImage(Image& image)
{
  for(int y = 0; y < image.domain().myUpperBound[1]; y++)
    for(int x = 0; x < image.domain().myUpperBound[0]; x++)
      image.setValue({x,y}, 255 - int(image.operator()({x,y})));
}

// Returns an image corresponding to the given DT
Image createImageFromDT(DTL1 dtl1, int maxValue)
{
  Image im(dtl1.domain());
  int step = maxValue != 0 ? 255 / maxValue : 0;
  int value;

  for(int y = 0; y < im.domain().upperBound()[1]; ++y)
    for(int x = 0; x < im.domain().upperBound()[0]; ++x)
    {
      value = dtl1.operator()({x,y});
      im.setValue({x,y}, value * step);
    }
  return im;
}

// Returns the addition of two images
Image addImages(Image im1, Image im2)
{
  Image im((im1.domain().size() > im2.domain().size() ? im1.domain() : im2.domain()));
  int value;

  for(int y = 0; y < im.domain().upperBound()[1]; ++y)
    for(int x = 0; x < im.domain().upperBound()[0]; ++x)
    {
      value = im1.operator()({x,y}) + im2.operator()({x,y});
      if(value > 255)
        value = 255;

      im.setValue({x, y}, value);
    }
  return im;
}

// Returns the addition of two DT in the form of an image
Image addImages(DTL1 dtl1Im1, DTL1 dtl1Im2)
{
  Image im((dtl1Im1.domain().size() > dtl1Im2.domain().size() ? dtl1Im1.domain() : dtl1Im2.domain()));

  DTL1::Value maxDT2 = (*std::max_element(dtl1Im1.constRange().begin(), dtl1Im1.constRange().end()));
  DTL1::Value maxDT1 = (*std::max_element(dtl1Im2.constRange().begin(), dtl1Im2.constRange().end()));

  int step = 255 / maxDT2;
  int step2 = 255 / maxDT1;

  for(int y = 0; y < im.domain().upperBound()[1]; ++y)
  {
    for(int x = 0; x < im.domain().upperBound()[0]; ++x)
    {
      int value = dtl1Im1.operator()({x,y});
      int value2 = dtl1Im2.operator()({x,y});
      
      im.setValue({x,y}, max(value * step, value2 * step2));
    }
  }
  return im;
}

void processImage(Image& image)
{
  // Create board and set its size
  Board2D board;
  board << image.domain();
  board.clear();

  Image imInv = image;

  // Import image to the board
  Display2DFactory::drawImage<Gray>(board, image, (unsigned char)0, (unsigned char)255);

  // Get the inverse
  inverseImage(imInv);

  // Threshold the image
  Binarizer b(image, 1, 135);
  Binarizer bInv(imInv, 1, 135);

  // DTL2 dtl2(&image.domain(), &b, &Z2i::l2Metric);
  DTL1 dtl1(&image.domain(), &b, &Z2i::l1Metric);
  DTL1 dtl1Inv(&imInv.domain(), &bInv, &Z2i::l1Metric);

  // Compute the max value of the DT
  DTL1::Value maxDT2 = (*std::max_element(dtl1.constRange().begin(), dtl1.constRange().end()));
  DTL1::Value maxDT1 = (*std::max_element(dtl1Inv.constRange().begin(), dtl1Inv.constRange().end()));

  // Hue color map
  // typedef DGtal::HueShadeColorMap<DTL1::Value,2> HueTwice;

  // Create GS DT images
  Image imGS = createImageFromDT(dtl1, maxDT2);
  Image imInvGS = createImageFromDT(dtl1Inv, maxDT1);
  // Add both DT
  Image imAddDTL = addImages(dtl1, dtl1Inv);

  // Save output 
  board.clear();
  Display2DFactory::drawImage<Gray>(board, imAddDTL, 0, 255);
  board.saveEPS("../output/im_add_DT");
  board.clear();
  Display2DFactory::drawImage<Gray>(board, imGS, 0, 255);
  board.saveEPS("../output/im_GS_DT.eps");
  board.clear();
  Display2DFactory::drawImage<Gray>(board, imInvGS, 0, 255);
  board.saveEPS("../output/im_inv_GS_DT");

  // Create a set to put the DT in
  // Set for the inverse of the image aswell
  Z2i::DigitalSet set(image.domain());
  Z2i::DigitalSet setInv(imInv.domain());
  DigitalSetInserter<Z2i::DigitalSet> inserter(set);
  DigitalSetInserter<Z2i::DigitalSet> inserterInv(setInv);
  DGtal::setFromImage(dtl1, inserter, 1, 135);
  DGtal::setFromImage(imInv, inserterInv, 1, 135);

  // Save both sets
  board.clear();
  board << set.domain() << set;
  board.saveEPS("../output/set.eps");
  board.clear();
  board << setInv.domain() << setInv;
  board.saveEPS("../output/set_Inv.eps");
}


int main(int argc, char** argv)
{
  // Check args
  if(argc != 2)
  {
    usage();
    return 0;
  }

  // Load image
  Image image = PGMReader<Image>::importPGM(argv[1]);
  processImage(image);

  return 0;
}
