#pragma once

#include "../include/images.h"
#include <DGtal/images/imagesSetsUtils/SetFromImage.h>
#include "DGtal/topology/KhalimskySpaceND.h"
#include <DGtal/topology/SurfelAdjacency.h>
#include <DGtal/topology/helpers/Surfaces.h>
#include <map>
#include <unordered_map>
#include "DGtal/topology/CubicalComplex.h"
#include "DGtal/topology/KhalimskyCellHashFunctions.h"

typedef int Integer;
typedef SpaceND<2, Integer> Z2;
typedef MetricAdjacency<Z2, 1> Adj4;
typedef MetricAdjacency<Z2, 2> Adj8;
typedef DigitalTopology<Adj4,Adj8> DT4_8;
typedef DigitalSetSelector<Z2i::Domain, BIG_DS+HIGH_BEL_DS>::Type DigitalSet;
typedef DGtal::Object<DT4_8, DigitalSet> ObjectType;

typedef KhalimskySpaceND<2,int> KSpace;
typedef std::map<Z2i::Cell, CubicalCellData> Map;
typedef CubicalComplex<KSpace, Map> CC;

// Returns a digital set from a thresholded image
Z2i::DigitalSet createDigitalSetFromImage(Image image);

// Returns a vector containing the connected components of the given object
std::vector<ObjectType> createObjectVector(ObjectType objT);

Z2i::Curve getBoundary(ObjectType& object);

std::vector<Z2i::SCell> getBoundaryVector(ObjectType obj);

KSpace initKSpace(Z2i::Point p1, Z2i::Point p2);

void getCCFromImage(Image im, CC& c, KSpace K);
