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

void findExtrema(DTL2 dtL2, float& min, float& max)
{
    for(DTL2::ConstRange::ConstIterator it = dtL2.constRange().begin(),
                itend=dtL2.constRange().end();
        it!=itend;
        ++it)
    {
        if(  (*it) < min )
            min=(*it);
        if( (*it) > max )
            max=(*it);
    }
}

void inverseImage(Image src, Image& dst) {
    for (int z = src.domain().lowerBound()[2]; z <= src.domain().upperBound()[2]; z++)
        for (int y = src.domain().lowerBound()[1]; y <= src.domain().upperBound()[1]; y++)
            for (int x = src.domain().lowerBound()[0]; x <= src.domain().upperBound()[0]; x++)
                dst.setValue({x, y, z}, 255 - src.operator()({x, y, z}));
}

void thresholdImage(Image src, Image& dst)
{
    for(Z3i::Domain::ConstIterator it = src.domain().begin(), itend=src.domain().end();
        it!=itend;
        ++it)
    {
        dst.setValue(*it, src(*it) > 0 ? 255 : 0);
    }

}

int main(int argc, char **argv) {
    QApplication application(argc, argv);
    string inputFilename = "../samples/cat10.vol";

    Image image = VolReader<Image>::importVol(inputFilename);
    Z3i::Domain domain(image.domain().lowerBound(), image.domain().upperBound());

    Image thresholdedIm(domain);
    thresholdImage(image, thresholdedIm);

    Image inverse(domain);
    inverseImage(thresholdedIm, inverse);

    Viewer3D<> viewer;

    Predicate aPredicate(thresholdedIm, 0);
    DTL2 dtL2(&domain, &aPredicate, &Z3i::l2Metric);

    float min = 0;
    float max = 0;
    findExtrema(dtL2, min, max);

    GradientColorMap<long> gradient( 0,30);
    gradient.addColor(Color::Red);
    gradient.addColor(Color::Yellow);
    gradient.addColor(Color::Green);
    gradient.addColor(Color::Cyan);
    gradient.addColor(Color::Blue);
    gradient.addColor(Color::Magenta);
    gradient.addColor(Color::Red);
    float transp = 30.;

    viewer << SetMode3D( (*(domain.begin())).className(), "Paving");
    for(Z3i::Domain::ConstIterator it = domain.begin(), itend=domain.end();
        it!=itend;
        ++it){

        double valDist= inverse( (*it) );
        cout << valDist << endl;
        Color c= gradient(valDist);

        if(valDist > 0)
        {
            viewer << CustomColors3D(Color((float)(c.red()),
                                           (float)(c.green()),
                                           (float)(c.blue(),transp)),
                                     Color((float)(c.red()),
                                           (float)(c.green()),
                                           (float)(c.blue()),transp));
            viewer << *it ;
        }

        // transp *= valDist;
        /*
        if(dtL2(*it)<=30 && image(*it)>0){
            viewer << CustomColors3D(Color((float)(c.red()),
                                           (float)(c.green()),
                                           (float)(c.blue(),transp)),
                                     Color((float)(c.red()),
                                           (float)(c.green()),
                                           (float)(c.blue()),transp));
            viewer << *it ;
        }*/

        // transp = 30.;
    }

    viewer << Viewer3D<>::updateDisplay;
    viewer.show();

    return application.exec();
}