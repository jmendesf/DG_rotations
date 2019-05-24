#include "../include/digital_objects.h"

Z2i::DigitalSet createDigitalSetFromImage(Image image)
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