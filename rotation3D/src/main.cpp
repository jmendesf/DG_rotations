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
namespace po = boost::program_options;
typedef ImageSelector<Z3i::Domain, unsigned char>::Type Image;


int main(int argc, char **argv) {
    QApplication application(argc,argv);

    string inputFilename = "../samples/Al.100.vol";


    Image image = VolReader<Image>::importVol(inputFilename);
    Z3i::Domain domain(image.domain().lowerBound(), image.domain().upperBound());

    Z3i::K3 ks;
    ks.init(image.domain().lowerBound(), image.domain().upperBound(), true);
    Viewer3D<Z3i::Space, Z3i::K3> viewer(ks);
    viewer.setCameraPosition(0, -800, 0);
    viewer.show();

    int thresholdMin = 30;
    int thresholdMax = 255;

    GradientColorMap<long> gradient( thresholdMin, thresholdMax);
    gradient.addColor(Color::Blue);
    gradient.addColor(Color::Green);
    gradient.addColor(Color::Yellow);
    gradient.addColor(Color::Red);

    for(Z3i::Domain::ConstIterator it = domain.begin(), itend=domain.end(); it!=itend; ++it){
        unsigned char  val= image( (*it) );
        Color c= gradient(val);
        if(val<=thresholdMax && val >=thresholdMin){
            viewer <<  CustomColors3D(DGtal::Color(255,0,0), DGtal::Color(255,0,0));
            viewer << *it;
        }
    }

    int nbNul = 0;
    int tot = 0;
    for (auto p : image.constRange())
    {
        if(int(p) != 0)
        {
            // cout << int(p) << endl;
        }
        else
            nbNul++;
        tot++;
    }
    cout << "Nb zero: " << nbNul << endl;
    cout << "Nb non null: " << tot - nbNul << endl;

    viewer.setAutoFillBackground(true);
    viewer << Viewer3D<>::updateDisplay;

    return application.exec();
}
