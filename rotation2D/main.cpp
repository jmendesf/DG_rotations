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
typedef ImageContainerBySTLVector<Z2i::Domain, unsigned char>  Image;
// Grayscale mapping
typedef GrayscaleColorMap<unsigned char> Gray;

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

void processImage(Image& image, Image& imInv)
{
  // Create board and set its size
  Board2D board;
  board << image.domain();
  board.clear();

  // import image to the board
  Display2DFactory::drawImage<Gray>(board, image, (unsigned char)0, (unsigned char)255);

  // Create a binarizer
  typedef functors::IntervalForegroundPredicate<Image> Binarizer;

  inverseImage(imInv);

  // Threshold the image
  Binarizer b(image, 1, 135);
  Binarizer bInv(imInv, 1, 135);

  // Create and apply a distance transform to the binary image with the L2 norm
  // Def L1 norm aswell
  typedef DistanceTransformation<Z2i::Space, Binarizer, Z2i::L2Metric> DTL2;
  typedef DistanceTransformation<Z2i::Space, Binarizer, Z2i::L1Metric> DTL1;

  DTL2 dtl2(&image.domain(), &b, &Z2i::l2Metric);
  DTL1 dtl1(&imInv.domain(), &bInv, &Z2i::l1Metric);

  // Compute the max value of the DT
  DTL2::Value maxDT2 = (*std::max_element(dtl2.constRange().begin(), dtl2.constRange().end()));
  DTL1::Value maxDT1 = (*std::max_element(dtl1.constRange().begin(), dtl1.constRange().end()));

  // Hue color map
  typedef DGtal::HueShadeColorMap<DTL2::Value,2> HueTwice;

  // Clear board
  board.clear();

  // Import GS DT image to the board
  // Save as output
  Display2DFactory::drawImage<Gray>(board, dtl2, (DTL2::Value)0, (DTL2::Value)maxDT2);
  board.saveEPS("../output/gray_scale.eps");
  board.clear();
  Display2DFactory::drawImage<Gray>(board, dtl1, (DTL1::Value)0, (DTL1::Value)maxDT1);
  board.saveEPS("../output/gray_scale_inv.eps");

  // Create a set to put the DT in
  // Set for the inverse of the image aswell
  Z2i::DigitalSet set(image.domain());
  Z2i::DigitalSet setInv(imInv.domain());
  DigitalSetInserter<Z2i::DigitalSet> inserter(set);
  DigitalSetInserter<Z2i::DigitalSet> inserterInv(setInv);
  DGtal::setFromImage(dtl2, inserter, 1, 135);
  DGtal::setFromImage(imInv, inserterInv, 1, 135);

  // Save both sets
  board.clear();
  board << set.domain() << set;
  board.saveEPS("../output/set.eps");
  board.clear();
  board << setInv.domain() << setInv;
  board.saveEPS("../output/setInv.eps");

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
  Image imInv = PGMReader<Image>::importPGM(argv[1]);
  processImage(image, imInv);

  return 0;
}
