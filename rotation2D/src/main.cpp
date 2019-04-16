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
#include "../include/tools.h" 
#include "../include/images.h"
#include "../include/rotations.h"

using std::cout;
using std::endl;
using std::string;


// Print if wrong args
void usage(bool help)
{
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

  cout << endl;
  // Compute rotations
  // Nearest neighbor
  if((method == NEAREST_NEIGHBOR) || (method == ALL))
  {
    cout << "-- Computing rotation using nearest neighbor -";
    Image rotIm = rotateBackward(imAddDTL, angle, NEAREST_NEIGHBOR);
    thresholdDTImage(rotIm, rotIm);
    saveImage(board, rotIm, 0, 255, "../output/rot_NN.eps");
    cout << " done." << endl;
    cout << "   Output saved as ../output/rot_NN.eps" << endl << endl;
  }
  
  // Bilinear interpolation
  if((method == BILINEAR_INTERPOLATION) || (method == ALL))
  {
    cout << "-- Computing rotation using bilinear interpolation -";
    Image rotIm = rotateBackward(imAddDTL, angle, BILINEAR_INTERPOLATION);
    thresholdDTImage(rotIm, rotIm);
    saveImage(board, rotIm, 0, 255, "../output/rot_BLI.eps");
    cout << " done." << endl;
    cout << "   Output saved as ../output/rot_BLI.eps" << endl << endl;
  }

  if((method == BICUBIC_INTERPOLATION) || (method == ALL))
  {
    cout << "-- Computing rotation using bicubic interpolation -";
    Image rotIm = rotateBackward(imAddDTL, angle, BICUBIC_INTERPOLATION);
    thresholdDTImage(rotIm, rotIm);
    saveImage(board, rotIm, 0, 255, "../output/rot_BIC.eps");
    cout << " done." << endl;
    cout << "   Output saved as ../output/rot_BIC.eps" << endl;
  }
  cout << endl;

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
    processImage(image, -std::stof(argv[2]), BILINEAR_INTERPOLATION, std::stoi(argv[4]), std::stoi(argv[5]));
  else if(strcmp(argv[3],"nn") == 0)
    processImage(image, -std::stof(argv[2]), NEAREST_NEIGHBOR, std::stoi(argv[4]), std::stoi(argv[5]));
  else if(strcmp(argv[3], "bic") == 0)
    processImage(image, -std::stof(argv[2]), BICUBIC_INTERPOLATION, std::stoi(argv[4]), std::stoi(argv[5]));
  else if(strcmp(argv[3], "all") == 0)
    processImage(image, -std::stof(argv[2]), ALL, std::stoi(argv[4]), std::stoi(argv[5]));
  else 
    usage(false); 

  return 0;
}
