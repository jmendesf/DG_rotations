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
#include <DGtal/io/writers/VolWriter.h>

// #include "../include/tools.h"
// #include "../include/images.h"
// #include "../include/rotations.h"

using std::cout;
using std::endl;
using std::string;
using std::min;
using std::max;
using namespace DGtal;
using namespace DGtal::Z3i;
namespace po = boost::program_options;

typedef ImageSelector<Z3i::Domain, double>::Type Image;
typedef functors::SimpleThresholdForegroundPredicate<Image> Predicate;
typedef DistanceTransformation<Z3i::Space, Predicate, Z3i::L2Metric> DTL2;

void usage() {
    cout << "usage: ./rotation3D angle aX aY aZ visualisation interp" << endl;
    cout << "visualisation: shape, dt, rot" << endl;
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
            imDT.setValue(*it, isInterior ? -1 : 1);
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
        dst.setValue(*it, value);
    }
}

// for debugging
void printMatrix(double matrix[3][3]) {
    for (int i = 0; i < 3; i++) {
        for (int y = 0; y < 3; y++)
            cout << matrix[i][y] << " ";
        cout << endl;
    }
}

void computeRotationMatrix(double v1, double v2, double v3, double matrix[3][3], double angle) {

    // Using Rodrigues' formula:
    // R = I + sin(angle) * A + (1 - cos(angle)) * A^2
    // With R the rotation matrix,
    // I the identity
    // A the matrix determined by the unit vector such that:
    //     ( 0   -v3   v2)
    // A = ( v3   0   -v1)
    //     (-v2   v1   0 )

    double w = v1 * v1 + v2 * v2 + v3 * v3;

    // A matrix
    double A[3][3] = {{0,   -v3, v2},
                      {v3,  0,   -v1},
                      {-v2, v1,  0}};
    // A matrix squared
    double ASquared[3][3] =
            {
                    {
                            A[0][0] * A[0][0] + A[0][1] * A[1][0] + A[0][2] * A[2][0],
                            A[0][0] * A[0][1] + A[0][1] * A[1][1] + A[0][2] * A[2][1],
                            A[0][0] * A[0][2] + A[0][1] * A[1][2] + A[0][2] * A[2][2]
                    },
                    {
                            A[1][0] * A[0][0] + A[1][1] * A[1][0] + A[1][2] * A[2][0],
                            A[1][0] * A[0][1] + A[1][1] * A[1][1] + A[1][2] * A[2][1],
                            A[1][0] * A[0][2] + A[1][1] * A[1][2] + A[1][2] * A[2][2]
                    },
                    {
                            A[2][0] * A[0][0] + A[2][1] * A[1][0] + A[2][2] * A[2][0],
                            A[2][0] * A[0][1] + A[2][1] * A[1][1] + A[2][2] * A[2][1],
                            A[2][0] * A[0][2] + A[2][1] * A[1][2] + A[2][2] * A[2][2]
                    }
            };

    for (int i = 0; i < 3; ++i) {
        for (int y = 0; y < 3; ++y) {
            ASquared[i][y] /= w;
            A[i][y] /= std::sqrt(w);
        }
    }

    float cosCoeff = 1 - cos(angle);

    // first operations
    for (int i = 0; i < 3; i++) {
        for (int y = 0; y < 3; y++) {
            A[i][y] *= sin(angle);
            ASquared[i][y] *= cosCoeff;
        }
    }

    // add identity matrix
    A[0][0] += 1;
    A[1][1] += 1;
    A[2][2] += 1;

    // add both matrices to the final matrix
    for (int i = 0; i < 3; i++) {
        for (int y = 0; y < 3; y++) {
            matrix[i][y] = A[i][y] + ASquared[i][y];
        }
    }
}

PointVector<3, double> computeRotation(PointVector<3, int> p, double rotationMatrix[3][3]) {
    PointVector<3, double> result;
    int x = p[0], y = p[1], z = p[2];
    result[0] = x * rotationMatrix[0][0] + y * rotationMatrix[0][1] + z * rotationMatrix[0][2];
    result[1] = x * rotationMatrix[1][0] + y * rotationMatrix[1][1] + z * rotationMatrix[1][2];
    result[2] = x * rotationMatrix[2][0] + y * rotationMatrix[2][1] + z * rotationMatrix[2][2];

    return result;
}

void computeTrilinearRotation(double x, double y, double z, double &value, Image srcImage) {
    int maxX = srcImage.domain().upperBound()[0];
    int maxY = srcImage.domain().upperBound()[1];
    int maxZ = srcImage.domain().upperBound()[2];

    int x0 = floor(x);
    int y0 = floor(y);
    int z0 = floor(z);

    int x1 = (x0 + 1) > maxX ? x0 : x0 + 1;
    int y1 = (y0 + 1) > maxY ? y0 : y0 + 1;
    int z1 = (z0 + 1) > maxZ ? z0 : z0 + 1;

    double xD = x - x0;
    double yD = y - y0;
    double zD = z - z0;

    double c000 = srcImage.operator()({x0, y0, z0});
    double c001 = srcImage.operator()({x0, y0, z1});
    double c010 = srcImage.operator()({x0, y1, z0});
    double c011 = srcImage.operator()({x0, y1, z1});
    double c100 = srcImage.operator()({x1, y0, z0});
    double c101 = srcImage.operator()({x1, y0, z1});
    double c110 = srcImage.operator()({x1, y1, z0});
    double c111 = srcImage.operator()({x1, y1, z1});

    double c00 = c000 * (1 - xD) + c100 * xD;
    double c01 = c001 * (1 - xD) + c101 * xD;
    double c10 = c010 * (1 - xD) + c110 * xD;
    double c11 = c011 * (1 - xD) + c111 * xD;

    double c0 = c00 * (1 - yD) + c10 * yD;
    double c1 = c01 * (1 - yD) + c11 * yD;

    value = c0 * (1 - zD) + c1 * zD;
}

Image rotateBackward(Image image, double angle, double v1, double v2, double v3, string interp) {
    int maxX = image.domain().upperBound()[0];
    int maxY = image.domain().upperBound()[1];
    int maxZ = image.domain().upperBound()[2];
    int minX = image.domain().lowerBound()[0];
    int minY = image.domain().lowerBound()[1];
    int minZ = image.domain().lowerBound()[2];

    PointVector<3, int> center((maxX - minX) / 2,
                               (maxY - minY) / 2,
                               (maxZ - minZ) / 2);

    double matrixForward[3][3];
    double matrixBackward[3][3];

    cout << "-- Computing Rodrigues' rotation matrices " << endl;
    cout << "     - Forward";
    computeRotationMatrix(v1, v2, v3, matrixForward, angle);
    cout << " - done." << endl;
    cout << "     - Backward";
    computeRotationMatrix(v1, v2, v3, matrixBackward, -angle);
    cout << " - done." << endl;

    cout << "-- Determining rotated domain boundaries by forward rotation ";
    PointVector<3, double> p000 = computeRotation({minX, minY, minZ}, matrixForward);
    PointVector<3, double> p100 = computeRotation({maxX, minY, minZ}, matrixForward);
    PointVector<3, double> p010 = computeRotation({minX, maxY, minZ}, matrixForward);
    PointVector<3, double> p110 = computeRotation({maxX, maxY, minZ}, matrixForward);
    PointVector<3, double> p001 = computeRotation({minX, minY, maxZ}, matrixForward);
    PointVector<3, double> p101 = computeRotation({maxX, minY, maxZ}, matrixForward);
    PointVector<3, double> p111 = computeRotation({maxX, maxY, maxZ}, matrixForward);
    PointVector<3, double> p011 = computeRotation({minX, maxY, maxZ}, matrixForward);

    int minXR, minYR, minZR;
    int maxXR, maxYR, maxZR;

    minXR = min(min(min(p000[0], p100[0]), min(p010[0], p110[0])), min(min(p001[0], p101[0]), min(p111[0], p011[0])));
    minYR = min(min(min(p000[1], p100[1]), min(p010[1], p110[1])), min(min(p001[1], p101[1]), min(p111[1], p011[1])));
    minZR = min(min(min(p000[2], p100[2]), min(p010[2], p110[2])), min(min(p001[2], p101[2]), min(p111[2], p011[2])));

    maxXR = max(max(max(p000[0], p100[0]), max(p010[0], p110[0])), max(max(p001[0], p101[0]), max(p111[0], p011[0])));
    maxYR = max(max(max(p000[1], p100[1]), max(p010[1], p110[1])), max(max(p001[1], p101[1]), max(p111[1], p011[1])));
    maxZR = max(max(max(p000[2], p100[2]), max(p010[2], p110[2])), max(max(p001[2], p101[2]), max(p111[2], p011[2])));

    PointVector<3, int> lowerBound = {minXR, minYR, minZR};
    PointVector<3, int> upperBound = {maxXR, maxYR, maxZR};
    Z3i::Domain newDomain(lowerBound, upperBound);
    Image imRot(newDomain);
    cout << "- done." << endl;

    double backX = 0;
    double backY = 0;
    double backZ = 0;

    double value = 0;

    for (int z = newDomain.lowerBound()[2]; z <= newDomain.upperBound()[2]; ++z) {
        for (int y = newDomain.lowerBound()[1]; y <= newDomain.upperBound()[1]; ++y) {
            for (int x = newDomain.lowerBound()[0]; x <= newDomain.upperBound()[0]; ++x) {
                PointVector<3, double> backP = computeRotation({x, y, z}, matrixBackward);

                if (interp.compare("nn") == 0) {
                    backX = (int) round(backP[0]);
                    backY = (int) round(backP[1]);
                    backZ = (int) round(backP[2]);
                }

                if (backP[0] <= minX || backP[1] <= minY || backP[2] <= minZ) {
                    imRot.setValue({x, y, z}, 1);
                    continue;
                }
                if (backP[0] >= maxX || backP[1] >= maxY || backP[2] >= maxZ) {
                    imRot.setValue({x, y, z}, 1);
                    continue;
                }

                if (interp.compare("tril") == 0) {
                    computeTrilinearRotation(backP[0], backP[1], backP[2], value, image);
                }

                if (interp.compare("nn") == 0)
                    value = image.operator()({(int) backX, (int) backY, (int) backZ});

                imRot.setValue({x, y, z}, value);
            }
        }
    }
    cout << "- done." << endl;
    return imRot;
}

Image DTToGrayscale(Image src, float min, float max) {
    if (min < 0)
        min = abs(min);
    Image gs(src.domain());

    float step1, step2;

    if (min == 0)
        min = 1;

    step1 = 255 / min;

    if (max == 0)
        max = 1;
    step2 = 255 / max;

    for (Z3i::Domain::ConstIterator it = src.domain().begin(), itend = src.domain().end();
         it != itend;
         ++it) {
        if (src(*it) < 0) {
            gs.setValue(*it, src(*it) * step1);
        } else if (src(*it) > 0) {
            gs.setValue(*it, src(*it) * step2);
        } else
            gs.setValue(*it, 0);
    }
    return gs;
}

void thresholdDTImage(Image src, Image &dst) {
    for (Z3i::Domain::ConstIterator it = src.domain().begin(), itend = src.domain().end();
         it != itend;
         ++it) {
        dst.setValue(*it, src(*it) >= 0 ? 0 : 255);
    }
}

void initGrad(GradientColorMap<double> &gradient) {
    gradient.addColor(Color::Red);
    gradient.addColor(Color::Yellow);
    gradient.addColor(Color::Green);
    gradient.addColor(Color::Cyan);
    gradient.addColor(Color::Blue);
    gradient.addColor(Color::Magenta);
    gradient.addColor(Color::Red);
}

int main(int argc, char **argv) {
    float vecRotation[3];
    float angle;
    string interp;

    if (argc == 7) {
        angle = stod(argv[1]);
        interp = argv[6];
        vecRotation[0] = stod(argv[2]);
        vecRotation[1] = stod(argv[3]);
        vecRotation[2] = stod(argv[4]);
    } else {
        usage();
        return 0;
    }

    if ((interp.compare("all") != 0) &&
        interp.compare("nn") != 0 &&
        interp.compare("tril") != 0) {
        trace.error();
        cout << "invalid interpolation argument: " << interp << endl;
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

    cout << "- Interpolation method chosen: " << interp << "." << endl;

    //Image image = VolReader<Image>::importVol(inputFilename);
    //Z3i::Domain domain(image.domain().lowerBound(), image.domain().upperBound());

    PointVector<3, int> lowerBound = {-20, -20, -20};
    PointVector<3, int> upperBound = {19, 19, 19};

    Z3i::Domain domain(lowerBound, upperBound);
    Image image(domain);
    int cubeX = 10, cubeY = 10, cubeZ = 10;
    for (int z = domain.lowerBound()[2]; z <= domain.upperBound()[2]; ++z) {
        for (int y = domain.lowerBound()[1]; y <= domain.upperBound()[1]; ++y) {
            for (int x = domain.lowerBound()[0]; x <= domain.upperBound()[0]; ++x) {
                if (x <= cubeX && y <= cubeY && z <= cubeZ)
                    if (x >= -cubeX && y >= -cubeY && z >= -cubeZ)
                        image.setValue({x, y, z}, 150);
                    else
                        image.setValue({x, y, z}, 0);
            }
        }
    }

    cout << "-- Thresholding ";
    Image thresholdedIm(domain);
    thresholdImage(image, thresholdedIm);
    cout << "- done." << endl;

    cout << "-- Inversing ";
    Image inverse(domain);
    inverseImage(thresholdedIm, inverse);
    cout << "- done." << endl;

    Viewer3D<> viewer1;
    Viewer3D<> viewer2;
    viewer1.setWindowTitle("NN");
    viewer2.setWindowTitle("Trilinear");

    cout << "-- Computing distance transforms " << endl;
    cout << "     - Foreground ";
    Predicate pIm(thresholdedIm, 0);
    DTL2 dtL2(&domain, &pIm, &Z3i::l2Metric);
    cout << "- done" << endl;

    cout << "     - Background ";
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

    cout << "-- Processing image DT ";
    processDT(dtL2Im, true);
    processDT(dtL2ImInv, false);
    cout << "- done." << endl;

    cout << "-- Merging both DTs ";
    addDTImages(dtL2Im, dtL2ImInv, DTAddIm);
    cout << "- done." << endl;


    Image imRotNN = DTAddIm;
    Image imRotTril = DTAddIm;

    if (interp.compare("all") == 0 || interp.compare("nn") == 0) {
        cout << "-- Computing rotation using nearest neighbor interpolation" << endl;
        imRotNN = rotateBackward(DTAddIm, angle, vecRotation[0], vecRotation[1], vecRotation[2], "nn");
    }
    if (interp.compare("all") == 0 || interp.compare("tril") == 0) {
        cout << "-- Computing rotation using trilinear interpolation" << endl;
        imRotTril = rotateBackward(DTAddIm, angle, vecRotation[0], vecRotation[1], vecRotation[2], "tril");
    }

    GradientColorMap<double> gradient(0, 255);
    initGrad(gradient);
    float transp = 20.;

    int i = 0;
    cout << "-- Populating viewers" << endl;
    viewer1 << SetMode3D((*(domain.begin())).className(), "Paving");
    if (strcmp(argv[5], "dt") == 0) {
        for (Z3i::Domain::ConstIterator it = domain.begin(), itend = domain.end();
             it != itend;
             ++it) {

            double valDist = imRotNN((*it));
            if (valDist > -max) {
                if (valDist < 0)
                    valDist = maxInv - abs(valDist);
                Color c = gradient(valDist);
                viewer1 << CustomColors3D(Color((float) (c.red()),
                                                (float) (c.green()),
                                                (float) (c.blue(), transp)),
                                          Color((float) (c.red()),
                                                (float) (c.green()),
                                                (float) (c.blue()), transp));
                viewer1 << *it;
            }
        }
    } else if (strcmp(argv[5], "shape") == 0) {
        for (Z3i::Domain::ConstIterator it = domain.begin(), itend = domain.end();
             it != itend;
             ++it) {

            double valDist = image((*it));
            if (valDist > 0) {
                Color c = gradient(valDist);
                viewer1 << CustomColors3D(Color((float) (c.red()),
                                                (float) (c.green()),
                                                (float) (c.blue(), 255)),
                                          Color((float) (c.red()),
                                                (float) (c.green()),
                                                (float) (c.blue()), 255));
                viewer1 << *it;
            }
        }
    } else if (strcmp(argv[5], "rot") == 0) {
        Image threshImRotNN(imRotNN.domain());
        Image threshImRotTril(imRotTril.domain());
        thresholdDTImage(imRotNN, threshImRotNN);
        thresholdDTImage(imRotTril, threshImRotTril);

        if (interp.compare("all") == 0 || interp.compare("nn") == 0) {
            for (Z3i::Domain::ConstIterator it = imRotNN.domain().begin(), itend = imRotNN.domain().end();
                 it != itend;
                 ++it) {
                double valDist = threshImRotNN((*it));

                if (valDist > 254) {
                    Color c = gradient(valDist);
                    viewer1 << CustomColors3D(Color((float) (c.red()),
                                                    (float) (c.green()),
                                                    (float) (c.blue(), 255)),
                                              Color((float) (c.red()),
                                                    (float) (c.green()),
                                                    (float) (c.blue()), 255));
                    viewer1 << *it;
                }
            }
        }

        if (interp.compare("all") == 0 || interp.compare("tril") == 0) {
            for (Z3i::Domain::ConstIterator it = imRotNN.domain().begin(), itend = imRotNN.domain().end();
                 it != itend;
                 ++it) {
                double valDist = threshImRotTril((*it));
                if (valDist > 254) {
                    Color c = gradient(valDist);
                    viewer2 << CustomColors3D(Color((float) (c.red()),
                                                    (float) (c.green()),
                                                    (float) (c.blue(), 255)),
                                              Color((float) (c.red()),
                                                    (float) (c.green()),
                                                    (float) (c.blue()), 255));
                    viewer2 << *it;
                }
            }
        }

    } else {
        trace.error();
        cout << "Unknown argument " << argv[5] << ". Exiting." << endl;
        return 0;
    }

    if (interp.compare("all") == 0 || interp.compare("nn") == 0) {
        viewer1 << Viewer3D<>::updateDisplay;
        viewer1.show();
    }

    if (interp.compare("tril") == 0) {
        viewer2 << Viewer3D<>::updateDisplay;
        viewer2.show();
    }

    Image gsDT = DTToGrayscale(DTAddIm, -max, maxInv);
    Image gsDT2 = DTToGrayscale(imRotTril, -max, maxInv);
    Image gsDT3 = DTToGrayscale(imRotNN, -max, maxInv);
    VolWriter<Image, functors::Cast<unsigned char>>::exportVol("../output/starting_DT.vol", gsDT);
    VolWriter<Image, functors::Cast<unsigned char>>::exportVol("../output/tril_DT.vol", gsDT2);
    VolWriter<Image, functors::Cast<unsigned char>>::exportVol("../output/NN_DT.vol", gsDT3);

    cout << "-- Visualising with option " << argv[5] << "." << endl;
    int appRet = application.exec();
    return appRet;
}
