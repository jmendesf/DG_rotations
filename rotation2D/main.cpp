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

typedef ImageContainerBySTLVector<Z2i::Domain, unsigned char>  Image;
typedef GrayscaleColorMap<unsigned char> Gray;

int main(int argc, char** argv)
{
    
  if(argc != 2)
  {
    cout << "usage: ./rotation <path_to_pgm_file>" << endl;
    return 0;
  }

  Image image = PGMReader<Image>::importPGM(argv[1]);
  trace.info() << "Imported image: " << image << endl;

  return 0;
}
