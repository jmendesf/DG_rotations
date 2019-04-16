#pragma once

#include <iostream>
#include <fstream>
#include <algorithm>
#include <string>
#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/base/BasicFunctors.h"
#include "DGtal/kernel/BasicPointPredicates.h"
#include "DGtal/images/ImageContainerBySTLVector.h"
#include "DGtal/images/ImageHelper.h"
#include "DGtal/geometry/volumes/distance/DistanceTransformation.h"
#include "DGtal/images/IntervalForegroundPredicate.h"
#include "DGtal/io/boards/Board2D.h"
#include "DGtal/io/colormaps/HueShadeColorMap.h"
#include "DGtal/io/colormaps/GrayscaleColorMap.h"
#include "Board/Point.h"

using namespace DGtal;
using namespace functors;

// Contains all image type definitions and tools to process an input image

// 2D image definition
typedef ImageContainerBySTLVector<Z2i::Domain, float>  Image;
// Grayscale mapping
typedef GrayscaleColorMap<float> Gray;
// Binarizer
typedef functors::IntervalForegroundPredicate<Image> Binarizer;  

// Create and apply a distance transform to the binary image with the L2 norm
// Def L1 norm aswell
typedef DistanceTransformation<Z2i::Space, Binarizer, Z2i::L2Metric> DTL2;
typedef DistanceTransformation<Z2i::Space, Binarizer, Z2i::L1Metric> DTL1;


void inverseImage(Image& image);
Image createImageFromDT(DTL1 dtl1, int maxValue, bool toGS);
Image createImageFromDT(DTL2 dtl2, int maxValue, bool toGS);
Image addImages(Image im1, Image im2);
Image addImages(DTL2 dtl2Im1, DTL2 dtl2Im2);
void imDTToGS(Image& imDT, int minValue, int maxValue);
void thresholdDTImage(Image src, Image& dst);
void processDT(Image& imDT, bool isInterior);