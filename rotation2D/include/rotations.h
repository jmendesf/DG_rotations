#pragma once

#include "../include/images.h"

// Enum on the interpolation method
enum INTERPOLATION_METHOD {
  NEAREST_NEIGHBOR,
  BILINEAR_INTERPOLATION,
  BICUBIC_INTERPOLATION,
  ALL
};

// Rotates the image using backward rotation technique
Image rotateBackward(Image image, float angle, INTERPOLATION_METHOD method);