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

#include <DGtal/images/imagesSetsUtils/SetFromImage.h>
#include "DGtal/topology/KhalimskySpaceND.h"
#include <DGtal/topology/SurfelAdjacency.h>
#include <DGtal/topology/helpers/Surfaces.h>
#include <map>
#include <unordered_map>
#include "DGtal/topology/CubicalComplex.h"
#include "DGtal/topology/KhalimskyCellHashFunctions.h"

#include "DGtal/kernel/sets/DigitalSetInserter.h"
#include <DGtal/images/imagesSetsUtils/SetFromImage.h>

using std::cout;
using std::endl;
using std::string;
using std::min;
using std::max;
using namespace DGtal;
namespace po = boost::program_options;

typedef ImageSelector<Z3i::Domain, double>::Type Image;
typedef functors::SimpleThresholdForegroundPredicate<Image> Predicate;
typedef DistanceTransformation<Z3i::Space, Predicate, Z3i::L2Metric> DTL2;

typedef KhalimskySpaceND<3, int> KSpace;
typedef std::map<Z3i::Cell, CubicalCellData> Map;
typedef CubicalComplex<KSpace, Map> CC;

typedef SpaceND<3, int> Z3;
typedef MetricAdjacency<Z3, 1> Adj6;
typedef MetricAdjacency<Z3, 2> Adj18;
typedef MetricAdjacency<Z3, 3> Adj26;
typedef DigitalTopology<Adj6, Adj18> DT6_18;
typedef DigitalTopology<Adj6, Adj6> DT6_6;
typedef DigitalTopology<Adj26, Adj6> DT26_6;

typedef DGtal::Object<DT6_6, Z3i::DigitalSet> ObjectType;
typedef DGtal::Object<DT26_6, Z3i::DigitalSet> ObjectType26_6;

void usage() {
    cout << "usage: ./rotation3D angle aX aY aZ visualisation interp [shape] [shape param]" << endl;
    cout << "visualisation: shape, rot" << endl;
    cout << "shape: " << endl;
    cout << "cube: length" << endl;
    cout << "sph (sphere): radius" << endl;
    cout << "el (ellipsoid): a, b, c" << endl;
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

// for debuggingy
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

    cout << "       Computing Rodrigues' rotation matrices... " << endl;
    cout << "           - Forward";
    computeRotationMatrix(v1, v2, v3, matrixForward, angle);
    cout << " - done." << endl;
    cout << "           - Backward";
    computeRotationMatrix(v1, v2, v3, matrixBackward, -angle);
    cout << " - done." << endl;

    cout << "       Determining rotated domain boundaries by forward rotation ";
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

    int diffX = maxXR - minXR;
    int diffy = maxYR - minYR;
    int diffz = maxZR - minZR;

    int total = (diffX + 1) * (diffy + 1) * (diffz + 1);

    Z3i::Domain newDomain(lowerBound, upperBound);
    Image imRot(newDomain);
    cout << "- done." << endl;

    double backX = 0;
    double backY = 0;
    double backZ = 0;
    double value = 0;
    int count = 1;


    for (int z = newDomain.lowerBound()[2]; z <= newDomain.upperBound()[2]; ++z) {
        for (int y = newDomain.lowerBound()[1]; y <= newDomain.upperBound()[1]; ++y) {
            for (int x = newDomain.lowerBound()[0]; x <= newDomain.upperBound()[0]; ++x) {
                cout << "       Pointwise rotation... " << count << "/" << total << ".\r";
                cout.flush();
                count++;

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
    cout << endl;
    cout << endl;
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
    gradient.addColor(Color::Yellow);
    gradient.addColor(Color::Yellow);
    gradient.addColor(Color::Green);
    gradient.addColor(Color::Cyan);
    gradient.addColor(Color::Blue);
    gradient.addColor(Color::Magenta);
    gradient.addColor(Color::Yellow);
}

bool isInsideEllipsoid(double a, double b, double c, int x, int y, int z) {
    return (double) ((x / a) * (x / a) + (y / b) * (y / b) + (z / c) * (z / c)) <= 1.;
}

double distanceToPoint(int x1, int y1, int z1, int x2, int y2, int z2) {
    return sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1) + (z2 - z1) * (z2 - z1));
}

KSpace initKSpace(Z3i::Point p1, Z3i::Point p2) {
    KSpace K;
    K.init(p1, p2, true);
    return K;
}

void getCCFromImage(Image im, CC &c, KSpace K) {
    for (Z3i::Domain::ConstIterator it = im.domain().begin(), itend = im.domain().end(); it != itend; ++it)
        if (im(*it) > 0)
            c.insertCell(K.uSpel(*it));
    c.close();
}

Z3i::DigitalSet createDigitalSetFromImage(Image image) {
    Z3i::DigitalSet set3d(image.domain());
    SetFromImage<Z3i::DigitalSet>::append<Image>(set3d, image, 1, 255);
    return set3d;
}

std::vector<ObjectType> createObjectVector(ObjectType objT) {
    std::vector<ObjectType> objects;
    std::back_insert_iterator<std::vector<ObjectType>> inserter(objects);
    objT.writeComponents(inserter);
    return objects;
}

std::vector<ObjectType26_6> createObjectVector(ObjectType26_6 objT) {
    std::vector<ObjectType26_6> objects;
    std::back_insert_iterator<std::vector<ObjectType26_6>> inserter(objects);
    objT.writeComponents(inserter);
    return objects;
}

void objectTypeToCubicalComplex(ObjectType obj, CC &cc, KSpace k) {
    for (auto it = obj.begin(), itend = obj.end(); it != itend; ++it)
        cc.insertCell(k.uSpel(*it));
    cc.close();
}

int main(int argc, char **argv) {
    float vecRotation[3];
    float angle;
    string interp, shape;

    if (argc == 7 || argc == 8 || argc == 9 || argc == 11) {
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

    Adj6 adj6;
    Adj18 adj18;
    Adj26 adj26;

    DT6_6 dt6_6(adj6, adj6, JORDAN_DT);
    DT26_6 dt26_6(adj26, adj6, JORDAN_DT);

    cout << endl;
    QApplication application(argc, argv);
    string inputFilename = "../samples/bunny.vol";
    cout << "- Rotation on shape " << inputFilename
         << " with " << angle << " rad angle and axis vector ("
         << vecRotation[0] << ","
         << vecRotation[1] << ","
         << vecRotation[2] << ")."
         << endl;

    cout << "- Interpolation method chosen: " << interp << "." << endl << endl;

    Image image = VolReader<Image>::importVol(inputFilename);
    Z3i::Domain domain(image.domain().lowerBound(), image.domain().upperBound());

    if (argc == 8 || argc == 9 || argc == 11) {
        shape = argv[7];
        PointVector<3, int> lowerBound = {-20, -20, -20};
        PointVector<3, int> upperBound = {19, 19, 19};

        domain = Z3i::Domain(lowerBound, upperBound);
        image = Image(domain);

        if (shape == "cube" && (argc == 8 || argc == 9)) {
            double cubeD;
            if (argc == 8)
                cubeD = 1.5;
            else
                cubeD = stod(argv[8]);

            for (int z = domain.lowerBound()[2]; z <= domain.upperBound()[2]; ++z) {
                for (int y = domain.lowerBound()[1]; y <= domain.upperBound()[1]; ++y) {
                    for (int x = domain.lowerBound()[0]; x <= domain.upperBound()[0]; ++x) {
                        if (x <= cubeD && y <= cubeD && z <= cubeD)
                            if (x > -cubeD && y > -cubeD && z > -cubeD)
                                image.setValue({x, y, z}, 150);
                            else
                                image.setValue({x, y, z}, 0);
                    }
                }
            }
        } else if (shape == "sph" && (argc == 8 || argc == 9)) {
            double r;
            if (argc == 8)
                r = 10;
            else
                r = stod(argv[8]);

            for (int z = domain.lowerBound()[2]; z <= domain.upperBound()[2]; ++z) {
                for (int y = domain.lowerBound()[1]; y <= domain.upperBound()[1]; ++y) {
                    for (int x = domain.lowerBound()[0]; x <= domain.upperBound()[0]; ++x) {
                        if (distanceToPoint(x, y, z, 0, 0, 0) <= r)
                            image.setValue({x, y, z}, 150);
                        else
                            image.setValue({x, y, z}, 0);
                    }
                }
            }
        } else if (shape == "el" && (argc == 8 || argc == 11)) {
            double a, b, c;
            if (argc == 8)
                a = 10, b = 15, c = 7;
            else {
                a = stod(argv[8]);
                b = stod(argv[9]);
                c = stod(argv[10]);
            }
            for (int z = domain.lowerBound()[2]; z <= domain.upperBound()[2]; ++z) {
                for (int y = domain.lowerBound()[1]; y <= domain.upperBound()[1]; ++y) {
                    for (int x = domain.lowerBound()[0]; x <= domain.upperBound()[0]; ++x) {
                        if (isInsideEllipsoid(a, b, c, x, y, z))
                            image.setValue({x, y, z}, 150);
                        else
                            image.setValue({x, y, z}, 0);
                    }
                }
            }
        } else if (shape == "plane" && argc == 8) {
            for (int z = domain.lowerBound()[2] + 3; z <= domain.upperBound()[2] - 3; ++z) {
                for (int y = domain.lowerBound()[1] + 3; y <= domain.upperBound()[1] - 3; ++y) {
                    for (int x = domain.lowerBound()[0] + 3; x <= domain.upperBound()[0] - 3; ++x) {
                        if ((2*x + y + z) < 5  && (2*x + y + z) > -5 )
                            image.setValue({x, y, z}, 150);
                        else
                            image.setValue({x, y, z}, 0);
                    }
                }
            }
        } else {
            trace.error();
            cout << " invalid shape: " << argv[7] << endl;
            return 0;
        }
    }

    if (strcmp(argv[5], "shape") == 0) {
        GradientColorMap<double> gradient(0, 255);
        initGrad(gradient);

        Viewer3D<> viewer1;
        viewer1 << SetMode3D((*(domain.begin())).className(), "PavingWired");
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
        viewer1 << Viewer3D<>::updateDisplay;
        viewer1.show();
        int appRet = application.exec();
        return appRet;
    }

    cout << "-- Thresholding ";
    Image thresholdedIm(domain);
    thresholdImage(image, thresholdedIm);
    cout << "- done." << endl;

    cout << "-- Inversing ";
    Image inverse(domain);
    inverseImage(thresholdedIm, inverse);
    cout << "- done." << endl;

    cout << "-- Building digital sets..." << endl;
    Z3i::DigitalSet imSet = createDigitalSetFromImage(image);
    Z3i::DigitalSet imSetInverse = createDigitalSetFromImage(inverse);

    cout << "-- Building corresponding objects and computing connected components..." << endl;
    ObjectType imObject(dt6_6, imSet);
    ObjectType imObjectInv(dt6_6, imSetInverse);
    ObjectType26_6 imObjectInv26_6(dt26_6, imSetInverse);

    std::vector<ObjectType> imObjects = createObjectVector(imObject);
    std::vector<ObjectType> imObjectsInv = createObjectVector(imObjectInv);
    std::vector<ObjectType26_6> imObjectsInv26_6 = createObjectVector(imObjectInv26_6);
    
    cout << "-- Building cubical complexes ";
    KSpace K = initKSpace(image.domain().lowerBound(), image.domain().upperBound());

    CC ccImInv(K);
    getCCFromImage(inverse, ccImInv, K);
    cout << "- done." << endl;

    int imInvB0 = imObjectsInv26_6.size();
    int imInvB2 = imObjects.size();
    int imInvB1 = imInvB0 + imInvB2 - ccImInv.euler();

    int imB0 = imInvB2;
    int imB1 = imInvB1;
    int imB2 = imInvB0 - 1;
    
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

    std::vector<std::vector<ObjectType>> objComponents;
    std::vector<std::vector<ObjectType>> objInvComponents;

    int rotB0NN = -1, rotB1NN = -1, rotB2NN = -1;
    int rotB0Tril = -1, rotB1Tril = -1, rotB2Tril = -1;

    Image threshImRotNN(DTAddIm.domain());
    Image threshImRotTril(DTAddIm.domain());

    if (interp.compare("all") == 0 || interp.compare("nn") == 0) {
        cout << "-- Computing rotation using nearest neighbor interpolation" << endl;
        imRotNN = rotateBackward(DTAddIm, angle, vecRotation[0], vecRotation[1], vecRotation[2], "nn");
        threshImRotNN = Image(imRotNN.domain());
        thresholdDTImage(imRotNN, threshImRotNN);

        Image imRotInv(threshImRotNN.domain());
        inverseImage(threshImRotNN, imRotInv);

        cout << "-- Building corresponding digital objects..." << endl;
        Z3i::DigitalSet rotSet = createDigitalSetFromImage(threshImRotNN);
        Z3i::DigitalSet rotInvSet = createDigitalSetFromImage(imRotInv);

        ObjectType objTRot(dt6_6, rotSet);
        ObjectType objTRotInv(dt6_6, rotInvSet);
        ObjectType26_6 objTRotInv26_6(dt26_6, rotInvSet);

        std::vector<ObjectType> connectedComponents = createObjectVector(objTRot);
        std::vector<ObjectType26_6> rotObjectsInv26_6 = createObjectVector(objTRotInv26_6);

        KSpace kRot = initKSpace(imRotInv.domain().lowerBound(), imRotInv.domain().upperBound());
        CC ccRotInv(kRot);
        getCCFromImage(imRotInv, ccRotInv, kRot);

        int invB0 = rotObjectsInv26_6.size();
        int invB2 = connectedComponents.size();
        int invB1 = invB0 + invB2 - ccRotInv.euler();

        rotB0NN = invB2;
        rotB1NN = invB1;
        rotB2NN = invB0 - 1;

        objComponents.push_back(connectedComponents);
        objInvComponents.push_back(createObjectVector(objTRotInv));
    }

    if (interp.compare("all") == 0 || interp.compare("tril") == 0) {
        cout << "-- Computing rotation using trilinear interpolation" << endl;
        imRotTril = rotateBackward(DTAddIm, angle, vecRotation[0], vecRotation[1], vecRotation[2], "tril");
        threshImRotTril = Image(imRotTril.domain());
        thresholdDTImage(imRotTril, threshImRotTril);

        Image imRotInv(threshImRotTril.domain());
        inverseImage(threshImRotTril, imRotInv);

        cout << "-- Building corresponding digital object..." << endl;
        Z3i::DigitalSet rotSet = createDigitalSetFromImage(threshImRotTril);
        Z3i::DigitalSet rotInvSet = createDigitalSetFromImage(imRotInv);

        ObjectType objTRot(dt6_6, rotSet);
        ObjectType objTRotInv(dt6_6, rotInvSet);
        ObjectType26_6 objTRotInv26_6(dt26_6, rotInvSet);

        std::vector<ObjectType> connectedComponents = createObjectVector(objTRot);
        std::vector<ObjectType26_6> rotObjectsInv26_6 = createObjectVector(objTRotInv26_6);

        KSpace kRot = initKSpace(imRotInv.domain().lowerBound(), imRotInv.domain().upperBound());
        CC ccRotInv(kRot);
        getCCFromImage(imRotInv, ccRotInv, kRot);

        int invB0 = rotObjectsInv26_6.size();
        int invB2 = connectedComponents.size();
        int invB1 = invB0 + invB2 - ccRotInv.euler();

        rotB0Tril = invB2;
        rotB1Tril = invB1;
        rotB2Tril = invB0 - 1;

        objComponents.push_back(connectedComponents);
        objInvComponents.push_back(createObjectVector(objTRotInv));
    }

    GradientColorMap<double> gradient(0, 255);
    initGrad(gradient);

    viewer1 << SetMode3D((*(domain.begin())).className(), "PavingWired");
    viewer2 << SetMode3D((*(domain.begin())).className(), "PavingWired");

     if (strcmp(argv[5], "rot") == 0) {
        if (interp.compare("all") == 0 || interp.compare("nn") == 0) {
            for (int i = 1; i < objInvComponents[0].size(); i++) {
                for (auto it = objInvComponents[0][i].begin(), itend = objInvComponents[0][i].end();
                     it != itend;
                     ++it) {
                    Color c = Color::Blue;
                    viewer1 << CustomColors3D(Color((float) (c.red()),
                                                    (float) (c.green()),
                                                    (float) (c.blue(), 230)),
                                              Color((float) (c.red()),
                                                    (float) (c.green()),
                                                    (float) (c.blue()), 230));
                    viewer1 << *it;
                }
            }

            for (int i = 0; i < objComponents[0].size(); i++) {
                Color c = objComponents[0][i].size() < 10 ? Color::Red : Color::Yellow;
                for (auto it = objComponents[0][i].begin(), itend = objComponents[0][i].end();
                     it != itend;
                     ++it) {

                    viewer1 << CustomColors3D(Color((float) (c.red()),
                                                    (float) (c.green()),
                                                    (float) (c.blue(), 230)),
                                              Color((float) (c.red()),
                                                    (float) (c.green()),
                                                    (float) (c.blue()), 230));
                    viewer1 << *it;
                }
            }

        }

        if (interp.compare("all") == 0 || interp.compare("tril") == 0) {
            int trilIndex = (interp == "tril") ? 0 : 1;
            for (int i = 1; i < objInvComponents[trilIndex].size(); i++) {
                for (auto it = objInvComponents[trilIndex][i].begin(), itend = objInvComponents[trilIndex][i].end();
                     it != itend;
                     ++it) {
                    Color c = Color::Blue;
                    viewer2 << CustomColors3D(Color((float) (c.red()),
                                                    (float) (c.green()),
                                                    (float) (c.blue(), 230)),
                                              Color((float) (c.red()),
                                                    (float) (c.green()),
                                                    (float) (c.blue()), 230));
                    viewer2 << *it;
                }
            }

            for (int i = 0; i < objComponents[trilIndex].size(); i++) {
                Color c = objComponents[trilIndex][i].size() > 10 ? Color::Yellow : Color::Red;
                for (auto it = objComponents[trilIndex][i].begin(), itend = objComponents[trilIndex][i].end();
                     it != itend;
                     ++it) {
                    viewer2 << CustomColors3D(Color((float) (c.red()),
                                                    (float) (c.green()),
                                                    (float) (c.blue(), 230)),
                                              Color((float) (c.red()),
                                                    (float) (c.green()),
                                                    (float) (c.blue()), 230));
                    viewer2 << *it;
                }
            }
        }

    } else {
        trace.error();
        cout << "Unknown argument " << argv[5] << ". Exiting." << endl;
        return 0;
    }

    cout << "-- Viewers populated." << endl;

    if (interp.compare("all") == 0 || interp.compare("nn") == 0 || strcmp(argv[5], "shape") == 0) {
        viewer1 << Viewer3D<>::updateDisplay;
        viewer2 << Viewer3D<>::updateDisplay;
        viewer1.show();
    }

    if (interp.compare("tril") == 0 && strcmp(argv[5], "shape") != 0) {
        viewer2 << Viewer3D<>::updateDisplay;
        viewer2.show();
    }

    Image gsDT = DTToGrayscale(DTAddIm, -max, maxInv);
    Image gsDT2 = DTToGrayscale(imRotTril, -max, maxInv);
    Image gsDT3 = DTToGrayscale(imRotNN, -max, maxInv);

    VolWriter<Image, functors::Cast<unsigned char>>::exportVol("../output/starting_DT.vol", gsDT);
    VolWriter<Image, functors::Cast<unsigned char>>::exportVol("../output/tril_DT.vol", gsDT2);
    VolWriter<Image, functors::Cast<unsigned char>>::exportVol("../output/NN_DT.vol", gsDT3);

    cout << "-- All output saved in the output/ folder." << endl;

    cout << "-- Visualising with option " << argv[5] << "." << endl;

    cout << endl;
    cout << "==================================================" << endl << endl;
    cout << "-- Topological informations:" << endl;
    cout << "   Original image: " << endl;
    cout << "       B0: " << imB0 << endl;
    cout << "       B1: " << imB1 << endl;
    cout << "       B2: " << imB2 << endl << endl;

    /*
    int i = 0;
    for (auto connComp : imObjects) {
        cout << "               Volume (component #" << i++ << "): " << connComp.size() << endl;
    }
    
    cout << "       - Nb of cavities                        : " << imObjectsInv.size() - 1 << endl;
    cout << endl;
    */

    if (interp == "all") {
        for (int i = 0; i < 2; i++) {
            if (i == 0)
            {
              cout << "   Nearest neighbor rotation: " << endl;
              cout << "       B0: " << rotB0NN << endl;
              cout << "       B1: " << rotB1NN << endl;
              cout << "       B2: " << rotB2NN << endl << endl;
            }
              
            else
            {
              cout << "   Trilinear interpolation rotation: " << endl;
              cout << "       B0: " << rotB0Tril << endl;
              cout << "       B1: " << rotB1Tril << endl;
              cout << "       B2: " << rotB2Tril << endl;
            }
                
        }
    } else {
        cout << "   Rotated image: " << endl;
        cout << "       - Nb connected components               : " << objComponents[0].size() << endl;
        cout << "             B0: " << (interp == "nn" ? rotB0NN : rotB0Tril) << endl;  
        cout << "             B1: " << (interp == "nn" ? rotB1NN : rotB1Tril) << endl;
        cout << "             B2: " << (interp == "nn" ? rotB2NN : rotB2Tril) << endl;
    }

    int appRet = application.exec();
    cout << endl;
    return appRet;
}
