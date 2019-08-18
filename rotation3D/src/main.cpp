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
                dst.setValue({x, y, z}, 255 - src({x, y, z}));
}

void DTToImage(DTL2 dtL2, double maxValue, Image &dst) {
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

    double c000 = srcImage({x0, y0, z0});
    double c001 = srcImage({x0, y0, z1});
    double c010 = srcImage({x0, y1, z0});
    double c011 = srcImage({x0, y1, z1});
    double c100 = srcImage({x1, y0, z0});
    double c101 = srcImage({x1, y0, z1});
    double c110 = srcImage({x1, y1, z0});
    double c111 = srcImage({x1, y1, z1});

    double c00 = c000 * (1 - xD) + c100 * xD;
    double c01 = c001 * (1 - xD) + c101 * xD;
    double c10 = c010 * (1 - xD) + c110 * xD;
    double c11 = c011 * (1 - xD) + c111 * xD;

    double c0 = c00 * (1 - yD) + c10 * yD;
    double c1 = c01 * (1 - yD) + c11 * yD;

    value = c0 * (1 - zD) + c1 * zD;
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

PointVector<3, double> computeVector(PointVector<3, double> a, PointVector<3, double> b) {
    return PointVector<3, double>(b[0] - a[0], b[1] - a[1], b[2] - a[2]);
}

PointVector<3, double> planeCoefficients(PointVector<3, double> p1, PointVector<3, double> p2) {
    double r1 = p1[1] * p2[2] - p1[2] * p2[1];
    double r2 = p1[2] * p2[0] - p1[0] * p2[2];
    double r3 = p1[0] * p2[1] - p1[1] * p2[0];

    return PointVector<3, double>(r1, r2, r3);
}

double getDOnPlaneEquation(PointVector<3, double> p, double a, double b, double c) {
    return -(p[0] * a + p[1] * b + p[2] * c);
}

PointVector<4, double> getPlaneEquation(PointVector<3, double> a, PointVector<3, double> b, PointVector<3, double> c) {
    PointVector<3, double> ab = computeVector(a, b);
    PointVector<3, double> ac = computeVector(a, c);
    PointVector<3, double> coeffs = planeCoefficients(ab, ac);
    double d = getDOnPlaneEquation(a, coeffs[0], coeffs[1], coeffs[2]);

    return PointVector<4, double>(coeffs[0], coeffs[1], coeffs[2], d);
}

bool isInUpperSpace(double a, double b, double c, double d, PointVector<3, double> p) {
    return (p[0] * a + p[1] * b + p[2] * c + d > 0);
}

bool isForeground(unsigned char cube, PointVector<3, double> p) {
    switch (cube) {
        case 0:
            return false;
        case 0 ^ 255:
            return true;
        case 16:
            return !isInUpperSpace(0.25, 0.25, 0.25, -0.125, p);
        case 16 ^ 255:
            return isInUpperSpace(0.25, 0.25, 0.25, -0.125, p);
        case 48:
            return !isInUpperSpace(0, 0.5, 0.5, -0.25, p);
        case 48 ^ 255:
            return isInUpperSpace(0, 0.5, 0.5, -0.25, p);
        case 80:
            return (!isInUpperSpace(0.25, 0.25, 0.25, -0.125, p) || !isInUpperSpace(-0.25, 0.25, -0.25, 0.375, p));
        case 80 ^ 255:
            return (isInUpperSpace(0.25, 0.25, 0.25, -0.125, p) && !isInUpperSpace(-0.25, 0.25, -0.25, 0.375, p));
        case 35:
            return (!isInUpperSpace(-0.25, -0.25, 0.25, 0.125, p) && !isInUpperSpace(0, 0, 1, -0.5, p));
        case 35 ^ 255:
            return (isInUpperSpace(-0.25, -0.25, 0.25, 0.125, p) || isInUpperSpace(0, 0, 1, -0.5, p));
        case 51:
            return (p[2] <= .5);
        case 51 ^ 255:
            return (p[2] >= .5);
        case 163:
            return (!isInUpperSpace(.25, .25, -.25, .125, p) ||
                    (!isInUpperSpace(-0.25, -0.25, 0.25, 0.125, p) && !isInUpperSpace(0, 0, 1, -0.5, p)));
        case 163 ^ 255:
            return (isInUpperSpace(.25, .25, -.25, .125, p) &&
                    (isInUpperSpace(-0.25, -0.25, 0.25, 0.125, p) || isInUpperSpace(0, 0, 1, -0.5, p)));
        case 90:
            return (!isInUpperSpace(0.25, 0.25, 0.25, -0.125, p) || !isInUpperSpace(-0.25, -0.25, 0.25, 0.375, p) ||
                    !isInUpperSpace(-0.25, 0.25, -0.25, 0.375, p) || !isInUpperSpace(0.25, -0.25, -0.25, 0.375, p));
        case 90 ^ 255:
            return (isInUpperSpace(0.25, 0.25, 0.25, -0.125, p) && isInUpperSpace(-0.25, -0.25, 0.25, 0.375, p) &&
                    isInUpperSpace(-0.25, 0.25, -0.25, 0.375, p) && isInUpperSpace(0.25, -0.25, -0.25, 0.375, p));
        case 20:
            return (!isInUpperSpace(0.25, 0.25, 0.25, -0.125, p) || !isInUpperSpace(-0.25, -0.25, -0.25, 0.625, p));
        case 20 ^ 255:
            return (isInUpperSpace(0.25, 0.25, 0.25, -0.125, p) && isInUpperSpace(-0.25, -0.25, -0.25, 0.625, p));
        case 52:
            return (!isInUpperSpace(0, 0.5, 0.5, -0.25, p) || !isInUpperSpace(-0.25, -0.25, -0.25, 0.625, p));
        case 52 ^ 255:
            return (isInUpperSpace(0, 0.5, 0.5, -0.25, p) && isInUpperSpace(-0.25, -0.25, -0.25, 0.625, p));
        case 164:
            return (!isInUpperSpace(-0.25, 0.25, 0.25, 0.125, p) || !isInUpperSpace(-0.25, -0.25, -0.25, 0.625, p) ||
                    !isInUpperSpace(.25, .25, -.25, .125, p));
        case 164 ^ 255:
            return (isInUpperSpace(-0.25, 0.25, 0.25, 0.125, p) && isInUpperSpace(-0.25, -0.25, -0.25, 0.625, p) &&
                    isInUpperSpace(.25, .25, -.25, .125, p));
        case 150:
            return (!isInUpperSpace(0.5, 0.5, 0, -0.25, p) || !isInUpperSpace(-0.5, -0.5, 0, 0.75, p));
        case 150 ^ 255:
            return (isInUpperSpace(0.5, 0.5, 0, -0.25, p) && isInUpperSpace(-0.5, -0.5, 0, 0.75, p));
        case 27:
            return (!isInUpperSpace(0.25, -0.25, 0.25, -0.125, p));
        case 27 ^ 255:
            return (isInUpperSpace(0.25, -0.25, 0.25, -0.125, p));
        case 43:
            return (!isInUpperSpace(-.5, -.5, 0, .25, p) && !isInUpperSpace(-.25, -.75, .25, .125, p) &&
                    !isInUpperSpace(.25, -.25, .75, -.625, p) && !isInUpperSpace(.5, 0, .5, -.75, p));
        case 43 ^ 255:
            return (isInUpperSpace(-.5, -.5, 0, .25, p) || isInUpperSpace(-.25, -.75, .25, .125, p) ||
                    isInUpperSpace(.25, -.25, .75, -.625, p) || isInUpperSpace(.5, 0, .5, -.75, p));
        case 23:
            if (p[1] < (-2 * p[0] + 1))
                return (!isInUpperSpace(.5, 0, .5, -.25, p));
            if (p[1] < p[0] - .5)
                return false;
            if (p[1] < (-.5 * p[0] + 1))
                return (!isInUpperSpace(.25, -.25, .75, -.125, p));
            if (p[1] < (-p[0] + 1.5))
                return (!isInUpperSpace(-.25, -.75, .25, .625, p));
            return true;
        case 23 ^ 255:
            if (p[1] < (-2 * p[0] + 1))
                return (isInUpperSpace(.5, 0, .5, -.25, p));
            if (p[1] < p[0] - .5)
                return true;
            if (p[1] < (-.5 * p[0] + 1))
                return (isInUpperSpace(.25, -.25, .75, -.125, p));
            if (p[1] < (-p[0] + 1.5))
                return (isInUpperSpace(-.25, -.75, .25, .625, p));
            return false;
        default:
            cout << "not computed yet." << endl;
            return false;
    }
}

PointVector<3, double> computeRotation(int rotationMatrix[3][3], PointVector<3, double> p) {
    double a = rotationMatrix[0][0] * p[0] + rotationMatrix[0][1] * p[1] + rotationMatrix[0][2] * p[2];
    double b = rotationMatrix[1][0] * p[0] + rotationMatrix[1][1] * p[1] + rotationMatrix[1][2] * p[2];
    double c = rotationMatrix[2][0] * p[0] + rotationMatrix[2][1] * p[1] + rotationMatrix[2][2] * p[2];

    return PointVector<3, double>(a, b, c);
}

PointVector<3, double> computeRotation(int **rotationMatrix, PointVector<3, double> p) {
    double a = rotationMatrix[0][0] * p[0] + rotationMatrix[0][1] * p[1] + rotationMatrix[0][2] * p[2];
    double b = rotationMatrix[1][0] * p[0] + rotationMatrix[1][1] * p[1] + rotationMatrix[1][2] * p[2];
    double c = rotationMatrix[2][0] * p[0] + rotationMatrix[2][1] * p[1] + rotationMatrix[2][2] * p[2];

    return PointVector<3, double>(a, b, c);
}

int **transposeMatrix(int matrix[3][3]) {
    int **dstMatrix = new int *[3];
    for (int i = 0; i < 3; i++)
        dstMatrix[i] = new int[3];

    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            dstMatrix[j][i] = matrix[i][j];
    return dstMatrix;
}

unsigned char pointToBinaryValue(PointVector<3, double> p) {
    if (p[0] == -.5) {
        if (p[1] == -.5) {
            if (p[2] == -.5)
                return 16;
            else
                return 128;
        } else if (p[2] == -.5) {
            return 1;
        } else
            return 8;
    } else {
        if (p[1] == -.5) {
            if (p[2] == -.5)
                return 32;
            else
                return 64;
        } else if (p[2] == -.5) {
            return 2;
        } else
            return 4;
    }
}

unsigned char computeRotatedCube(unsigned char cube, int rotationMatrix[3][3]) {
    unsigned char resultingCube = 0;

    if ((cube & (unsigned char) 128) != 0)
        resultingCube += pointToBinaryValue(computeRotation(rotationMatrix, PointVector<3, double>(-.5, -.5, .5)));
    if ((cube & (unsigned char) 64) != 0)
        resultingCube += pointToBinaryValue(computeRotation(rotationMatrix, PointVector<3, double>(.5, -.5, .5)));
    if ((cube & (unsigned char) 32) != 0)
        resultingCube += pointToBinaryValue(computeRotation(rotationMatrix, PointVector<3, double>(.5, -.5, -.5)));
    if ((cube & (unsigned char) 16) != 0)
        resultingCube += pointToBinaryValue(computeRotation(rotationMatrix, PointVector<3, double>(-.5, -.5, -.5)));
    if ((cube & (unsigned char) 8) != 0)
        resultingCube += pointToBinaryValue(computeRotation(rotationMatrix, PointVector<3, double>(-.5, .5, .5)));
    if ((cube & (unsigned char) 4) != 0)
        resultingCube += pointToBinaryValue(computeRotation(rotationMatrix, PointVector<3, double>(.5, .5, .5)));
    if ((cube & (unsigned char) 2) != 0)
        resultingCube += pointToBinaryValue(computeRotation(rotationMatrix, PointVector<3, double>(.5, .5, -.5)));
    if ((cube & (unsigned char) 1) != 0)
        resultingCube += pointToBinaryValue(computeRotation(rotationMatrix, PointVector<3, double>(-.5, .5, -.5)));

    return resultingCube;
}

unsigned char computeRotatedCube(unsigned char cube, int **rotationMatrix) {
    unsigned char resultingCube = 0;

    if ((cube & (unsigned char) 128) != 0)
        resultingCube += pointToBinaryValue(computeRotation(rotationMatrix, PointVector<3, double>(-.5, -.5, .5)));
    if ((cube & (unsigned char) 64) != 0)
        resultingCube += pointToBinaryValue(computeRotation(rotationMatrix, PointVector<3, double>(.5, -.5, .5)));
    if ((cube & (unsigned char) 32) != 0)
        resultingCube += pointToBinaryValue(computeRotation(rotationMatrix, PointVector<3, double>(.5, -.5, -.5)));
    if ((cube & (unsigned char) 16) != 0)
        resultingCube += pointToBinaryValue(computeRotation(rotationMatrix, PointVector<3, double>(-.5, -.5, -.5)));
    if ((cube & (unsigned char) 8) != 0)
        resultingCube += pointToBinaryValue(computeRotation(rotationMatrix, PointVector<3, double>(-.5, .5, .5)));
    if ((cube & (unsigned char) 4) != 0)
        resultingCube += pointToBinaryValue(computeRotation(rotationMatrix, PointVector<3, double>(.5, .5, .5)));
    if ((cube & (unsigned char) 2) != 0)
        resultingCube += pointToBinaryValue(computeRotation(rotationMatrix, PointVector<3, double>(.5, .5, -.5)));
    if ((cube & (unsigned char) 1) != 0)
        resultingCube += pointToBinaryValue(computeRotation(rotationMatrix, PointVector<3, double>(-.5, .5, -.5)));

    return resultingCube;
}

std::map<unsigned char, int **> computeRotationMatrices() {
    std::map<unsigned char, int **> matrixFromBinary;
    std::vector<unsigned char> cubes = {16, 16 ^ 255, 48, 48 ^ 255, 80, 80 ^ 255, 35, 35 ^ 255, 51, 51 ^ 255,
                                        163, 163 ^ 255, 90, 90 ^ 255, 20, 20 ^ 255, 52, 52 ^ 255, 164, 164 ^ 255,
                                        164, 164 ^ 255, 150, 150 ^ 255, 27, 27 ^ 255, 43, 43 ^ 255, 23, 23 ^ 255};
    for (auto cube : cubes) {
        int rotationMatrix0[3][3] = {{1, 0, 0},
                                     {0, 1, 0},
                                     {0, 0, 1}};
        unsigned char res = computeRotatedCube(cube, rotationMatrix0);
        if (matrixFromBinary.find(res) == matrixFromBinary.end())
            matrixFromBinary[res] = transposeMatrix(rotationMatrix0);

        int rotationMatrix1[3][3] = {{1, 0, 0},
                                     {0, 0, -1},
                                     {0, 1, 0}};
        res = computeRotatedCube(cube, rotationMatrix1);
        if (matrixFromBinary.find(res) == matrixFromBinary.end())
            matrixFromBinary[res] = transposeMatrix(rotationMatrix1);

        int rotationMatrix2[3][3] = {{1, 0,  0},
                                     {0, -1, 0},
                                     {0, 0,  -1}};
        res = computeRotatedCube(cube, rotationMatrix2);
        if (matrixFromBinary.find(res) == matrixFromBinary.end())
            matrixFromBinary[res] = transposeMatrix(rotationMatrix2);

        int rotationMatrix3[3][3] = {{1, 0,  0},
                                     {0, 0,  1},
                                     {0, -1, 0}};
        res = computeRotatedCube(cube, rotationMatrix3);
        if (matrixFromBinary.find(res) == matrixFromBinary.end())
            matrixFromBinary[res] = transposeMatrix(rotationMatrix3);

        int rotationMatrix4[3][3] = {{0, -1, 0},
                                     {1, 0,  0},
                                     {0, 0,  1}};
        res = computeRotatedCube(cube, rotationMatrix4);
        if (matrixFromBinary.find(res) == matrixFromBinary.end())
            matrixFromBinary[res] = transposeMatrix(rotationMatrix4);

        int rotationMatrix5[3][3] = {{0, 0, 1},
                                     {1, 0, 0},
                                     {0, 1, 0}};
        res = computeRotatedCube(cube, rotationMatrix5);
        if (matrixFromBinary.find(res) == matrixFromBinary.end())
            matrixFromBinary[res] = transposeMatrix(rotationMatrix5);

        int rotationMatrix6[3][3] = {{0, 1, 0},
                                     {1, 0, 0},
                                     {0, 0, -1}};
        res = computeRotatedCube(cube, rotationMatrix6);
        if (matrixFromBinary.find(res) == matrixFromBinary.end())
            matrixFromBinary[res] = transposeMatrix(rotationMatrix6);

        int rotationMatrix7[3][3] = {{0, 0,  -1},
                                     {1, 0,  0},
                                     {0, -1, 0}};
        res = computeRotatedCube(cube, rotationMatrix7);
        if (matrixFromBinary.find(res) == matrixFromBinary.end())
            matrixFromBinary[res] = transposeMatrix(rotationMatrix7);

        int rotationMatrix8[3][3] = {{-1, 0,  0},
                                     {0,  -1, 0},
                                     {0,  0,  1}};
        res = computeRotatedCube(cube, rotationMatrix8);
        if (matrixFromBinary.find(res) == matrixFromBinary.end())
            matrixFromBinary[res] = transposeMatrix(rotationMatrix8);

        int rotationMatrix9[3][3] = {{-1, 0,  0},
                                     {0,  0,  -1},
                                     {0,  -1, 0}};
        res = computeRotatedCube(cube, rotationMatrix9);
        if (matrixFromBinary.find(res) == matrixFromBinary.end())
            matrixFromBinary[res] = transposeMatrix(rotationMatrix9);

        int rotationMatrix10[3][3] = {{-1, 0, 0},
                                      {0,  1, 0},
                                      {0,  0, -1}};
        res = computeRotatedCube(cube, rotationMatrix10);
        if (matrixFromBinary.find(res) == matrixFromBinary.end())
            matrixFromBinary[res] = transposeMatrix(rotationMatrix10);

        int rotationMatrix11[3][3] = {{-1, 0, 0},
                                      {0,  0, 1},
                                      {0,  1, 0}};
        res = computeRotatedCube(cube, rotationMatrix11);
        if (matrixFromBinary.find(res) == matrixFromBinary.end())
            matrixFromBinary[res] = transposeMatrix(rotationMatrix11);

        int rotationMatrix12[3][3] = {{0,  1, 0},
                                      {-1, 0, 0},
                                      {0,  0, 1}};
        res = computeRotatedCube(cube, rotationMatrix12);
        if (matrixFromBinary.find(res) == matrixFromBinary.end())
            matrixFromBinary[res] = transposeMatrix(rotationMatrix12);

        int rotationMatrix13[3][3] = {{0,  0,  1},
                                      {-1, 0,  0},
                                      {0,  -1, 0}};
        res = computeRotatedCube(cube, rotationMatrix13);
        if (matrixFromBinary.find(res) == matrixFromBinary.end())
            matrixFromBinary[res] = transposeMatrix(rotationMatrix13);

        int rotationMatrix14[3][3] = {{0,  -1, 0},
                                      {-1, 0,  0},
                                      {0,  0,  -1}};
        res = computeRotatedCube(cube, rotationMatrix14);
        if (matrixFromBinary.find(res) == matrixFromBinary.end())
            matrixFromBinary[res] = transposeMatrix(rotationMatrix14);

        int rotationMatrix15[3][3] = {{0,  0, -1},
                                      {-1, 0, 0},
                                      {0,  1, 0}};
        res = computeRotatedCube(cube, rotationMatrix15);
        if (matrixFromBinary.find(res) == matrixFromBinary.end())
            matrixFromBinary[res] = transposeMatrix(rotationMatrix15);

        int rotationMatrix16[3][3] = {{0, 0, -1},
                                      {0, 1, 0},
                                      {1, 0, 0}};
        res = computeRotatedCube(cube, rotationMatrix16);
        if (matrixFromBinary.find(res) == matrixFromBinary.end())
            matrixFromBinary[res] = transposeMatrix(rotationMatrix16);

        int rotationMatrix17[3][3] = {{0, 1, 0},
                                      {0, 0, 1},
                                      {1, 0, 0}};
        res = computeRotatedCube(cube, rotationMatrix17);
        if (matrixFromBinary.find(res) == matrixFromBinary.end())
            matrixFromBinary[res] = transposeMatrix(rotationMatrix17);

        int rotationMatrix18[3][3] = {{0, 0,  1},
                                      {0, -1, 0},
                                      {1, 0,  0}};
        res = computeRotatedCube(cube, rotationMatrix18);
        if (matrixFromBinary.find(res) == matrixFromBinary.end())
            matrixFromBinary[res] = transposeMatrix(rotationMatrix18);

        int rotationMatrix19[3][3] = {{0, -1, 0},
                                      {0, 0,  -1},
                                      {1, 0,  0}};
        res = computeRotatedCube(cube, rotationMatrix19);
        if (matrixFromBinary.find(res) == matrixFromBinary.end())
            matrixFromBinary[res] = transposeMatrix(rotationMatrix19);

        int rotationMatrix20[3][3] = {{0,  0,  -1},
                                      {0,  -1, 0},
                                      {-1, 0,  0}};
        res = computeRotatedCube(cube, rotationMatrix20);
        if (matrixFromBinary.find(res) == matrixFromBinary.end())
            matrixFromBinary[res] = transposeMatrix(rotationMatrix20);

        int rotationMatrix21[3][3] = {{0,  -1, 0},
                                      {0,  0,  1},
                                      {-1, 0,  0}};
        res = computeRotatedCube(cube, rotationMatrix21);
        if (matrixFromBinary.find(res) == matrixFromBinary.end())
            matrixFromBinary[res] = transposeMatrix(rotationMatrix21);
        int rotationMatrix22[3][3] = {{0,  0, 1},
                                      {0,  1, 0},
                                      {-1, 0, 0}};
        res = computeRotatedCube(cube, rotationMatrix22);
        if (matrixFromBinary.find(res) == matrixFromBinary.end())
            matrixFromBinary[res] = transposeMatrix(rotationMatrix22);

        int rotationMatrix23[3][3] = {{0,  1, 0},
                                      {0,  0, -1},
                                      {-1, 0, 0}};
        res = computeRotatedCube(cube, rotationMatrix23);
        if (matrixFromBinary.find(res) == matrixFromBinary.end())
            matrixFromBinary[res] = transposeMatrix(rotationMatrix23);

    }
    return matrixFromBinary;
}

bool pointFromImageIsForeground(unsigned char cube, int **matrix, PointVector<3, double> p) {
    if (cube == 0)
        return false;
    if (cube == 255)
        return true;

    unsigned char refCube = computeRotatedCube(cube, matrix);
    PointVector<3, double> repositionnedPoint(p[0] - .5, p[1] - .5, p[2] - .5);

    repositionnedPoint = computeRotation(matrix, repositionnedPoint);
    return isForeground(refCube, PointVector<3, double>(repositionnedPoint[0] + .5, repositionnedPoint[1] + .5,
                                                        repositionnedPoint[2] + .5));
}

bool computeRotationMC(PointVector<3, double> p, Image image, std::map<unsigned char, int **> matrixFromBinary) {
    int x = floor(p[0]);
    int y = floor(p[1]);
    int z = floor(p[2]);

    unsigned char cube = 0;
    if (image({x, y, z + 1}) > 0)
        cube += 128;
    if (image({x + 1, y, z + 1}) > 0)
        cube += 64;
    if (image({x + 1, y, z}) > 0)
        cube += 32;
    if (image({x, y, z}) > 0)
        cube += 16;
    if (image({x, y + 1, z + 1}) > 0)
        cube += 8;
    if (image({x + 1, y + 1, z + 1}) > 0)
        cube += 4;
    if (image({x + 1, y + 1, z}) > 0)
        cube += 2;
    if (image({x, y + 1, z}) > 0)
        cube += 1;

    PointVector<3, double> normP(p[0] - x, p[1] - y, p[2] - z);
    return pointFromImageIsForeground(cube, matrixFromBinary[cube], normP);
}

Image rotateBackward(Image image, double angle, double v1, double v2, double v3, string interp,
                     std::map<unsigned char, int **> matrixFromBinary) {
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
    int count = 0;

    for (int z = newDomain.lowerBound()[2]; z <= newDomain.upperBound()[2]; ++z) {
        for (int y = newDomain.lowerBound()[1]; y <= newDomain.upperBound()[1]; ++y) {
            for (int x = newDomain.lowerBound()[0]; x <= newDomain.upperBound()[0]; ++x) {
                cout << "       Pointwise rotation... " << count << "/" << total << ".\r";
                cout.flush();
                count++;

                PointVector<3, double> backP = computeRotation({x, y, z}, matrixBackward);

                if (backP[0] <= minX || backP[1] <= minY || backP[2] <= minZ) {
                    imRot.setValue({x, y, z}, 1);
                    continue;
                }
                if (backP[0] >= maxX || backP[1] >= maxY || backP[2] >= maxZ) {
                    imRot.setValue({x, y, z}, 1);
                    continue;
                }
                if (interp.compare("mc") == 0)
                    value = computeRotationMC(backP, image, matrixFromBinary) ? 255 : 0;

                imRot.setValue({x, y, z}, value);
            }
        }
    }
    cout << endl;
    cout << endl;
    return imRot;
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
                    value = image({(int) backX, (int) backY, (int) backZ});

                imRot.setValue({x, y, z}, value);

            }
        }
    }
    cout << endl;
    cout << endl;
    return imRot;
}

int main(int argc, char **argv) {
    float vecRotation[3];
    float angle;
    string interp, shape;

    std::map<unsigned char, int **> matrixFromBinary = computeRotationMatrices();

    if (argc == 7 || argc == 8 || argc == 9 || argc == 11) {
        angle = stod(argv[1]);
        interp = argv[6];
        vecRotation[0] = std::stod(argv[2]);
        vecRotation[1] = std::stod(argv[3]);
        vecRotation[2] = std::stod(argv[4]);
    } else {
        usage();
        return 0;
    }

    if (interp != "all" &&
        interp != "nn" &&
        interp != "tril" &&
        interp != "mc") {
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
                        if ((2 * x + y + z) < 5 && (2 * x + y + z) > -5)
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
        viewer1 << SetMode3D((*(domain.begin())).className(), "Paving");
        for (Z3i::Domain::ConstIterator it = domain.begin(), itend = domain.end();
             it != itend;
             ++it) {

            double valDist = image((*it));
            if (valDist > 0) {
                Color c = gradient(valDist);
                viewer1 << CustomColors3D(Color((float) (c.red()),
                                                (float) (c.green()),
                                                (float) (c.blue(), 50)),
                                          Color((float) (c.red()),
                                                (float) (c.green()),
                                                (float) (c.blue()), 50));
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
    Viewer3D<> viewer3;
    viewer1.setWindowTitle("NN");
    viewer2.setWindowTitle("Trilinear");
    viewer3.setWindowTitle("Marching Cubes");

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
    Image imRotMC = DTAddIm;

    std::vector<std::vector<ObjectType>> objComponents;
    std::vector<std::vector<ObjectType>> objInvComponents;

    int rotB0NN = -1, rotB1NN = -1, rotB2NN = -1;
    int rotB0Tril = -1, rotB1Tril = -1, rotB2Tril = -1;
    int rotB0MC = -1, rotB1MC = -1, rotB2MC = -1;

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

    if(interp == "all" || interp == "mc") {
        cout << "-- Computing rotation using Marching Cubes algorithm" << endl;
        imRotMC = rotateBackward(image, angle, vecRotation[0], vecRotation[1], vecRotation[2], "mc", matrixFromBinary);

        Image imRotInv(imRotMC.domain());
        inverseImage(imRotMC, imRotInv);

        cout << "-- Building corresponding digital objects..." << endl;
        Z3i::DigitalSet rotSet = createDigitalSetFromImage(imRotMC);
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

        rotB0MC = invB2;
        rotB1MC = invB1;
        rotB2MC = invB0 - 1;

        objComponents.push_back(connectedComponents);
        objInvComponents.push_back(createObjectVector(objTRotInv));
    }

    GradientColorMap<double> gradient(0, 255);
    initGrad(gradient);

    viewer1 << SetMode3D((*(domain.begin())).className(), "Paving");
    viewer2 << SetMode3D((*(domain.begin())).className(), "Paving");

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
                                                    (float) (c.blue(), 50)),
                                              Color((float) (c.red()),
                                                    (float) (c.green()),
                                                    (float) (c.blue()), 50));
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
                                                    (float) (c.blue(), 190)),
                                              Color((float) (c.red()),
                                                    (float) (c.green()),
                                                    (float) (c.blue()), 190));
                    viewer2 << *it;
                }
            }
        }

        if (interp.compare("all") == 0 || interp.compare("mc") == 0) {
            int mcIndex = (interp == "mc") ? 0 : 2;
            for (int i = 1; i < objInvComponents[mcIndex].size(); i++) {
                for (auto it = objInvComponents[mcIndex][i].begin(), itend = objInvComponents[mcIndex][i].end();
                     it != itend;
                     ++it) {
                    Color c = Color::Blue;
                    viewer3 << CustomColors3D(Color((float) (c.red()),
                                                    (float) (c.green()),
                                                    (float) (c.blue(), 230)),
                                              Color((float) (c.red()),
                                                    (float) (c.green()),
                                                    (float) (c.blue()), 230));
                    viewer3 << *it;
                }
            }

            for (int i = 0; i < objComponents[mcIndex].size(); i++) {
                Color c = objComponents[mcIndex][i].size() > 10 ? Color::Yellow : Color::Red;
                for (auto it = objComponents[mcIndex][i].begin(), itend = objComponents[mcIndex][i].end();
                     it != itend;
                     ++it) {
                    viewer3 << CustomColors3D(Color((float) (c.red()),
                                                    (float) (c.green()),
                                                    (float) (c.blue(), 190)),
                                              Color((float) (c.red()),
                                                    (float) (c.green()),
                                                    (float) (c.blue()), 190));
                    viewer3 << *it;
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
        //viewer1.show();
    }

    if (interp.compare("tril") == 0 && strcmp(argv[5], "shape") != 0) {
        viewer2 << Viewer3D<>::updateDisplay;
        //viewer2.show();
    }

    if (interp.compare("mc") == 0 && strcmp(argv[5], "shape") != 0) {
        viewer3 << Viewer3D<>::updateDisplay;
        viewer3.show();
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

    if (interp == "all") {
        for (int i = 0; i < 2; i++) {
            if (i == 0) {
                cout << "   Nearest neighbor rotation: " << endl;
                cout << "       B0: " << rotB0NN << endl;
                cout << "       B1: " << rotB1NN << endl;
                cout << "       B2: " << rotB2NN << endl << endl;
            } else if (i == 1){
                cout << "   Trilinear interpolation rotation: " << endl;
                cout << "       B0: " << rotB0Tril << endl;
                cout << "       B1: " << rotB1Tril << endl;
                cout << "       B2: " << rotB2Tril << endl;
            } else {
                cout << "   Marching Cubes method rotation: " << endl;
                cout << "       B0: " << rotB0MC << endl;
                cout << "       B1: " << rotB1MC << endl;
                cout << "       B2: " << rotB2MC << endl;
            }

        }
    } else {
        cout << "   Rotated image: " << endl;
        cout << "       - Nb connected components               : " << objComponents[0].size() << endl;
        cout << "             B0: " << (interp == "nn" ? rotB0NN : (interp == "tril" ? rotB0Tril : rotB0MC)) << endl;
        cout << "             B1: " << (interp == "nn" ? rotB1NN : (interp == "tril" ? rotB1Tril : rotB1MC)) << endl;
        cout << "             B2: " << (interp == "nn" ? rotB2NN : (interp == "tril" ? rotB2Tril : rotB2MC)) << endl;
    }

    int appRet = application.exec();
    cout << endl;
    return appRet;
}
