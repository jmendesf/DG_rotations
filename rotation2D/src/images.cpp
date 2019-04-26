#include "../include/images.h"

// Inverse the given image
void inverseImage(Image &image) {
    for (int y = image.domain().lowerBound()[1]; y <= image.domain().upperBound()[1]; y++)
        for (int x = image.domain().lowerBound()[0]; x <= image.domain().upperBound()[0]; x++)
            image.setValue({x, y}, 255 - int(image.operator()({x, y})));
}

// Returns an image corresponding to the given DT
Image createImageFromDT(DTL1 dtl1, int maxValue, bool toGS) {
    Image im(dtl1.domain());
    int step = 1;
    if (toGS)
        step = maxValue != 0 ? 255 / maxValue : 0;

    float value;

    for (int y = im.domain().lowerBound()[1]; y <= im.domain().upperBound()[1]; ++y)
        for (int x = im.domain().lowerBound()[0]; x <= im.domain().upperBound()[0]; ++x) {
            value = dtl1.operator()({x, y});
            im.setValue({x, y}, toGS ? step * value : value);
        }
    return im;
}

Image createImageFromDT(DTL2 dtl2, int maxValue, bool toGS) {
    Image im(dtl2.domain());
    int step = 1;
    if (toGS)
        step = maxValue != 0 ? 255 / maxValue : 0;

    float value = 0;

    for (int y = im.domain().lowerBound()[1]; y <= im.domain().upperBound()[1]; ++y)
        for (int x = im.domain().lowerBound()[0]; x <= im.domain().upperBound()[0]; ++x) {
            value = dtl2.operator()({x, y});
            im.setValue({x, y}, toGS ? step * value : value);
        }
    return im;
}

// Returns the addition of two images
Image addImages(Image im1, Image im2) {
    // Takes the size of the biggest image
    Image im((im1.domain().size() > im2.domain().size() ? im1.domain() : im2.domain()));
    float value;

    for (int y = im.domain().lowerBound()[1]; y <= im.domain().upperBound()[1]; ++y)
        for (int x = im.domain().lowerBound()[0]; x <= im.domain().upperBound()[0]; ++x) {
            value = im1.operator()({x, y}) + im2.operator()({x, y});
            if (value == 0.) {
                value = -0.5;
            }

            // image value cannot be more than 255
            if (value > 255.)
                value = 255.;

            im.setValue({x, y}, value);
        }
    return im;
}

// Returns the addition of two DT in the form of an image
Image addImages(DTL2 dtl2Im1, DTL2 dtl2Im2) {
    // Takes the size of the biggest image
    Image im((dtl2Im1.domain().size() > dtl2Im2.domain().size() ? dtl2Im1.domain() : dtl2Im2.domain()));

    // Max values of DT for each DT
    DTL1::Value maxDT2 = (*std::max_element(dtl2Im1.constRange().begin(), dtl2Im1.constRange().end()));
    DTL1::Value maxDT1 = (*std::max_element(dtl2Im2.constRange().begin(), dtl2Im2.constRange().end()));

    // Compute grayscale coefficient depending on the max DT value
    int step = maxDT2 != 0 ? 255 / maxDT2 : 1;
    int step2 = maxDT1 != 0 ? 255 / maxDT1 : 1;

    float value, value2;

    for (int y = im.domain().lowerBound()[1]; y <= im.domain().upperBound()[1]; ++y) {
        for (int x = im.domain().lowerBound()[0]; x <= im.domain().upperBound()[0]; ++x) {
            value = dtl2Im1.operator()({x, y});
            value2 = dtl2Im2.operator()({x, y});

            // the point at x,y takes the maximum value between the two
            im.setValue({x, y}, std::max(value * step, value2 * step2));
        }
    }
    return im;
}

// Converts a DT image into its grayscale equivalent (values from 0 to 255)
void imDTToGS(Image &imDT, int minValue, int maxValue) {
    // avoid division by 0
    int step = (maxValue != minValue) ? 255 / (maxValue - minValue) : 1;
    float value;

    // compute equivalent for DT in GS
    for (int y = imDT.domain().lowerBound()[1]; y < imDT.domain().upperBound()[1]; ++y) {
        for (int x = imDT.domain().lowerBound()[0]; x < imDT.domain().upperBound()[0]; ++x) {
            value = imDT.operator()({x, y});

            if (imDT.operator()({x, y}) < 0)
                value = 255;

            imDT.setValue({x, y}, abs(step * value));
        }
    }
}

// if pixel is positive: 0 (background)
// if pixel is negative: 255 (foreground)
void thresholdDTImage(Image src, Image &dst) {
    for (int y = dst.domain().lowerBound()[1]; y <= dst.domain().upperBound()[1]; y++) {
        for (int x = dst.domain().lowerBound()[0]; x <= dst.domain().upperBound()[0]; x++) {
            if (src.domain().isInside(PointVector<2, int>(x, y)))
                dst.setValue({x, y}, src.operator()({x, y}) < 0 ? 255 : 0);
        }
    }
}

// Prepares the DT for rotation. 
// Contour pixels take the value of 0.5 if belonging to the BG, -0.5 if FG.
void processDT(Image &imDT, bool isInterior) {
    for (int y = imDT.domain().lowerBound()[1]; y <= imDT.domain().upperBound()[1]; y++) {
        for (int x = imDT.domain().lowerBound()[0]; x <= imDT.domain().upperBound()[0]; x++) {
            if(imDT.operator()({x,y}) == 1.){
               imDT.setValue({x,y}, isInterior ? -0.5 : 0.5);
             } else {
            if (isInterior && (imDT.operator()({x, y}) != 0))
                imDT.setValue({x, y}, -imDT.operator()({x, y}));
            }
        }
    }
}

// Supposed to resize a given image by modifying the domain
// Cant modify the domain yet
Z2i::Domain getResizedDomain(Image &image) {
    int xMin = 0, xMax = 0, yMin = 0, yMax = 0;
    bool breakFor = false;

    for (int y = image.domain().lowerBound()[1]; y <= image.domain().upperBound()[1]; ++y) {
        for (int x = image.domain().lowerBound()[0]; x <= image.domain().upperBound()[0]; ++x) {
            if (image.operator()({x, y}) < 0) {
                yMin = y;
                breakFor = true;
                break;
            }
        }
        if (breakFor)
            break;
    }
    breakFor = false;

    for (int x = image.domain().lowerBound()[0]; x < image.domain().upperBound()[0]; ++x) {
        for (int y = image.domain().lowerBound()[1]; y < image.domain().upperBound()[1]; ++y) {
            if (image.operator()({x, y}) < 0) {
                xMin = x;
                breakFor = true;
                break;
            }
        }
        if (breakFor)
            break;
    }
    breakFor = false;


    for (int x = image.domain().upperBound()[0]; x > image.domain().lowerBound()[0]; --x) {
        for (int y = image.domain().lowerBound()[1]; y < image.domain().upperBound()[1]; ++y) {
            if (image.operator()({x, y}) < 0) {
                xMax = x;
                breakFor = true;
                break;
            }
        }
        if (breakFor)
            break;
    }
    breakFor = false;

    for (int y = image.domain().upperBound()[1]; y > image.domain().lowerBound()[1]; --y) {
        for (int x = image.domain().lowerBound()[0]; x < image.domain().upperBound()[0]; ++x) {
            if (image.operator()({x, y}) < 0) {
                yMax = y;
                breakFor = true;
                break;
            }
        }
        if (breakFor)
            break;
    }

    return Z2i::Domain(Z2i::Point(xMin, yMin), Z2i::Point(xMax, yMax));
}

ImagePGM thresholdToPGM(Image image)
{
    ImagePGM pgm(getResizedDomain(image));

    for(int y = pgm.domain().lowerBound()[1]; y <= pgm.domain().upperBound()[1]; ++y)
    {
        for(int x = pgm.domain().lowerBound()[0]; x <= pgm.domain().upperBound()[0]; ++x)
        {
            pgm.setValue({x,y}, image.operator()({x,y}) >= 0 ? 0 : 255);
        }
    }
    return pgm;
}