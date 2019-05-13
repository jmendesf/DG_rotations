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

int main(int argc, char **argv) {
    QApplication application(argc, argv);
    string inputFilename = "../samples/lobster.vol";

    Image image = VolReader<Image>::importVol(inputFilename);
    Image GSim = image;
    Z3i::Domain domain(image.domain().lowerBound(), image.domain().upperBound());
    Viewer3D<> viewer;

    Predicate aPredicate(GSim, 0);
    DTL2 dtL2(&domain, &aPredicate, &Z3i::l2Metric);

    float min = 0;
    float max = 0;
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

    cout << "Min Value: " << min << endl;
    cout << "Max Value: " << max << endl;

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

        double valDist= dtL2( (*it) );
        Color c= gradient(valDist);
        transp *= valDist;
        if(dtL2(*it)<=30 && image(*it)>0){
            viewer << CustomColors3D(Color((float)(c.red()),
                                           (float)(c.green()),
                                           (float)(c.blue(),transp)),
                                     Color((float)(c.red()),
                                           (float)(c.green()),
                                           (float)(c.blue()),transp));
            viewer << *it ;
        }
        transp = 30.;
    }

    viewer << Viewer3D<>::updateDisplay;
    viewer.show();

    return application.exec();
}