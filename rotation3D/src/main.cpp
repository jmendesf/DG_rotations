#include <iostream>
#include "DGtal/base/Common.h"
#include "DGtal/kernel/SpaceND.h"
#include "DGtal/kernel/domains/HyperRectDomain.h"
#include "DGtal/images/ImageSelector.h"
#include "DGtal/io/readers/VolReader.h"
#include "boost/program_options.hpp"
#include <DGtal/io/DrawWithDisplay3DModifier.h>

#include "DGtal/io/Color.h"
#include "DGtal/io/colormaps/GradientColorMap.h"
#include "DGtal/io/DrawWithDisplay3DModifier.h"
#include "DGtal/io/viewers/Viewer3D.h"

#include "DGtal/kernel/SpaceND.h"
#include "DGtal/kernel/domains/HyperRectDomain.h"
#include "DGtal/images/ImageSelector.h"

#include "DGtal/geometry/volumes/distance/DistanceTransformation.h"
#include "DGtal/images/SimpleThresholdForegroundPredicate.h"
#include "DGtal/helpers/StdDefs.h"

#include "DGtal/helpers/StdDefs.h"

// #include "../include/tools.h"
// #include "../include/images.h"
// #include "../include/rotations.h"

using std::cout;
using std::endl;
using std::string;
using namespace DGtal;
using namespace DGtal::Z3i;
namespace po = boost::program_options;

typedef ImageSelector<Z3i::Domain, float>::Type Image;
typedef functors::SimpleThresholdForegroundPredicate<Image> Predicate;
typedef DistanceTransformation<Z3i::Space, Predicate, Z3i::L2Metric> DTL2;

void usage() {
    cout << "usage: ./rotation3D angle aX aY aZ visualisation" << endl;
    cout << "visualisation: shape, dt" << endl;
}

void findExtrema(DTL2 dtL2, float &min, float &max) {
    for (DTL2::ConstRange::ConstIterator it = dtL2.constRange().begin(),
                 itend = dtL2.constRange().end();
         it != itend;
         ++it) {
        if ((*it) < min)
            min = (*it);
        if ((*it) > max)
            max = (*it);
    }
}

void inverseImage(Image src, Image &dst) {
    for (int z = src.domain().lowerBound()[2]; z <= src.domain().upperBound()[2]; z++)
        for (int y = src.domain().lowerBound()[1]; y <= src.domain().upperBound()[1]; y++)
            for (int x = src.domain().lowerBound()[0]; x <= src.domain().upperBound()[0]; x++)
                dst.setValue({x, y, z}, 255 - src.operator()({x, y, z}));
}

void DTToImage(DTL2 dtL2, double maxValue, Image &dst) {
    int step = 1;
    float value = 0;

    for (Z3i::Domain::ConstIterator it = dtL2.domain().begin(), itend = dtL2.domain().end();
         it != itend;
         ++it) {
        dst.setValue(*it, dtL2(*it));
    }
}

void thresholdImage(Image src, Image &dst) {
    for (Z3i::Domain::ConstIterator it = src.domain().begin(), itend = src.domain().end();
         it != itend;
         ++it) {
        dst.setValue(*it, src(*it) > 0 ? 255 : 0);
    }
}

void processDT(Image &imDT, bool isInterior) {
    for (Z3i::Domain::ConstIterator it = imDT.domain().begin(), itend = imDT.domain().end();
         it != itend;
         ++it) {
        if (imDT(*it) == 1) {
            imDT.setValue(*it, isInterior ? -0.5 : 0.5);
        } else {
            if (isInterior && imDT(*it) != 0)
                imDT.setValue(*it, -imDT(*it));
        }
    }
}

void addDTImages(Image im1, Image im2, Image &dst) {
    float value;
    for (Z3i::Domain::ConstIterator it = dst.domain().begin(), itend = dst.domain().end();
         it != itend;
         ++it) {
        value = im1(*it) + im2(*it);
        if (value == 0.)
            value = -0.5;
        if (value > 255)
            value = 255;
        dst.setValue(*it, value);
    }
}

int main(int argc, char **argv) {
    float vecRotation[3];
    float angle;

    if (argc == 6) {
        angle = stof(argv[1]);
        vecRotation[0] = stof(argv[2]);
        vecRotation[1] = stof(argv[3]);
        vecRotation[2] = stof(argv[4]);
    } else {
        usage();
        return 0;
    }


    QApplication application(argc, argv);
    string inputFilename = "../samples/cat10.vol";
    cout << "- Rotation on shape " << inputFilename
         << " with " << angle << " rad angle and axis vector ("
         << vecRotation[0] << ","
         << vecRotation[1] << ","
         << vecRotation[2] << ")."
         << endl;


    Image image = VolReader<Image>::importVol(inputFilename);
    Z3i::Domain domain(image.domain().lowerBound(), image.domain().upperBound());

    cout << "-- Thresholding ";
    Image thresholdedIm(domain);
    thresholdImage(image, thresholdedIm);
    cout << "- done." << endl;

    cout << "-- Inversing ";
    Image inverse(domain);
    inverseImage(thresholdedIm, inverse);
    cout << "- done." << endl;

    Viewer3D<> viewer;

    cout << "-- Computing distance transforms " << endl;
    cout << "     - foreground ";
    Predicate pIm(thresholdedIm, 0);
    DTL2 dtL2(&domain, &pIm, &Z3i::l2Metric);
    cout << "- done" << endl;

    cout << "     - background ";
    Predicate pInv(inverse, 0);
    DTL2 dtL2Inv(&domain, &pInv, &Z3i::l2Metric);
    cout << "- done." << endl;

    float min = 0, max = 0, minInv = 0, maxInv = 0;
    Image dtL2Im(domain), dtL2ImInv(domain), DTAddIm(domain);

    cout << "-- Converting DTL2 to Image ";
    findExtrema(dtL2, min, max);
    DTToImage(dtL2, max, dtL2Im);
    findExtrema(dtL2Inv, minInv, maxInv);
    DTToImage(dtL2Inv, maxInv, dtL2ImInv);
    cout << "- done." << endl;

    cout << "-- processing image DT ";
    processDT(dtL2Im, true);
    processDT(dtL2ImInv, false);
    cout << "- done." << endl;

    cout << "-- Merging both DTs ";
    addDTImages(dtL2Im, dtL2ImInv, DTAddIm);
    cout << "- done." << endl;

    GradientColorMap<long> gradient(-max - 1, maxInv + 1);
    gradient.addColor(Color::Red);
    gradient.addColor(Color::Yellow);
    gradient.addColor(Color::Green);
    gradient.addColor(Color::Cyan);
    gradient.addColor(Color::Blue);
    gradient.addColor(Color::Magenta);
    gradient.addColor(Color::Red);
    float transp = 20.;

    viewer << SetMode3D((*(domain.begin())).className(), "Paving");

    if(strcmp(argv[5], "dt") == 0)
    {
        for (Z3i::Domain::ConstIterator it = domain.begin(), itend = domain.end();
             it != itend;
             ++it) {

            double valDist = DTAddIm((*it));
            if (valDist > -max) {
                if (valDist < 0)
                    valDist = maxInv - abs(valDist);
                Color c = gradient(valDist);
                viewer << CustomColors3D(Color((float) (c.red()),
                                               (float) (c.green()),
                                               (float) (c.blue(), transp)),
                                         Color((float) (c.red()),
                                               (float) (c.green()),
                                               (float) (c.blue()), transp));
                viewer << *it;
            }
        }
    } else if(strcmp(argv[5], "shape") == 0)
    {
        for (Z3i::Domain::ConstIterator it = domain.begin(), itend = domain.end();
             it != itend;
             ++it) {

            double valDist = image((*it));
            if (valDist > 0) {
                Color c = gradient(valDist);
                viewer << CustomColors3D(Color((float) (c.red()),
                                               (float) (c.green()),
                                               (float) (c.blue(), 255)),
                                         Color((float) (c.red()),
                                               (float) (c.green()),
                                               (float) (c.blue()), 255));
                viewer << *it;
            }
        }
    }

    else
    {
        cout << "Unknown argument " << argv[5] << ". Exiting." << endl;
        return 0;
    }



    viewer << Viewer3D<>::updateDisplay;
    viewer.show();

    return application.exec();
}