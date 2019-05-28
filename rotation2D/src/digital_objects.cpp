#include "../include/digital_objects.h"

using namespace Z2i;

DigitalSet createDigitalSetFromImage(Image image)
{
  Z2i::DigitalSet set2d(image.domain());
  SetFromImage<Z2i::DigitalSet>::append<Image>(set2d, image, 1, 255);
  return set2d;
}

std::vector<ObjectType> createObjectVector(ObjectType objT)
{
  std::vector<ObjectType> objects;
  std::back_insert_iterator<std::vector<ObjectType>> inserter(objects);
  objT.writeComponents(inserter);
  return objects;
}

Curve getBoundary(ObjectType& object)
{
  KSpace kSpace;
  kSpace.init(object.domain().lowerBound(),
              object.domain().upperBound(), 
              true);
  SCell cell = Surfaces<KSpace>::findABel(kSpace, object.pointSet(), 10000);
  std::vector<Point> boundaryPoints;
  SurfelAdjacency<2> SAdj(true);
  Surfaces<KSpace>::track2DBoundaryPoints(boundaryPoints, kSpace, SAdj, object.pointSet(), cell);

  Curve boundaryCurve;
  boundaryCurve.initFromPointsVector(boundaryPoints);
  return boundaryCurve;
}

std::vector<Z2i::SCell> getBoundaryVector(ObjectType obj)
{
  KSpace kSpace;
  kSpace.init(obj.domain().lowerBound(), obj.domain().upperBound(), true);

  SCell cell = Surfaces<KSpace>::findABel(kSpace, obj.pointSet(), 15000);
  std::vector<SCell> vectBdrySCell;
  SurfelAdjacency<2> sAdj(true);
  Surfaces<KSpace>::track2DBoundary(vectBdrySCell, kSpace, sAdj, obj.pointSet(), cell);
  return vectBdrySCell;
}

KSpace initKSpace(Point p1, Point p2)
{
  KSpace K;
  K.init(p1, p2, true);
  return K;
}