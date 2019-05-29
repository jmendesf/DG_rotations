#include <boost/algorithm/minmax_element.hpp>
#include "DGtal/kernel/sets/DigitalSetInserter.h"
#include "DGtal/io/readers/PGMReader.h"
#include "Board/Point.h"
#include "../include/tools.h"
#include "../include/images.h"
#include "../include/rotations.h"
#include "../include/digital_objects.h"
#include "DGtal/io/writers/PGMWriter.h"


using std::cout;
using std::endl;
using std::string;

// Print if wrong args
void usage(bool help) {
    cout << endl;
    cout << "-- Usage" << endl;
    cout << "   ./rotation2D <file> <angle> <method> [<minThresh> <maxThresh>]" << endl;
    cout << "   method: nn (nearest neighbor), bli (bilinear), bic (bicubic), all " << endl;

    if (help) {
        cout << endl;
        cout << "-- Help" << endl;
        cout << "   Threshold values are needed if the image is not binary" << endl;
        cout << "   example: ./rotation2D ../samples/contourS.pgm 2.5 all 1 135" << endl;
        cout << "   (if binary): ./rotation2D ../samples/sw.pgm 2.5 all" << endl;
    } else
        cout << " - ./rotation2D help for more" << endl;

    cout << endl;
}

// Save board from image
// Save to path
void saveImage(Board2D board, Image image, int minVal, int maxVal, string path) {
    board.clear();
    Display2DFactory::drawImage<Gray>(board, image, minVal, maxVal);
    board.saveEPS(path.c_str());
}

// Save board from set
// Save to path
void saveSet(Board2D board, Z2i::DigitalSet set, string path) {
    board.clear();
    board << set.domain() << set;
    board.saveEPS(path.c_str());
}

template<class T>
void sendToBoard( Board2D & board, T & p_Object, DGtal::Color p_Color) {
    board << CustomStyle( p_Object.className(), new DGtal::CustomFillColor(p_Color));
    board << p_Object;
}

void send1CellsVectorToBoard(std::vector<Z2i::SCell> bdryV, Board2D& board)
{
    for(auto cell : bdryV)
        sendToBoard(board, cell, Color::Red);

}

void processImage(Image &image, float angle, INTERPOLATION_METHOD method, int minThresh, int maxThresh) {
    Adj4 adj4;
    Adj8 adj8;
    DT4_8 dt4_8(adj4, adj8, JORDAN_DT);

    // Create board and set its size
    Board2D board;
    board << image.domain();
    board.clear();

    // copy the image and inverse it
    Image imInv = image;

    cout << endl;
    cout << "-- Inversing -";
    inverseImage(imInv);
    cout << " done." << endl;

    Z2i::DigitalSet imSet = createDigitalSetFromImage(image);
    Z2i::DigitalSet imInvSet = createDigitalSetFromImage(imInv);

    ObjectType imObject(dt4_8, imSet);
    ObjectType imInvObject(dt4_8, imInvSet);

    std::vector<ObjectType> imObjects = createObjectVector(imObject);
    std::vector<ObjectType> imInvObjects = createObjectVector(imInvObject);
    std::vector<Z2i::SCell> bdryVect = getBoundaryVector(imObjects[0]);
    std::vector<Z2i::SCell> bdryVectInv = getBoundaryVector(imInvObjects[0]);

    cout << "-- Building cubical complexes... "; 

    KSpace K = initKSpace(image.domain().lowerBound(), image.domain().upperBound());

    CC ccIm(K);
    CC ccImInv(K);
    getCCFromImage(image, ccIm, K);
    getCCFromImage(imInv, ccImInv, K);

    cout << "- done." << endl;

    board << CustomStyle( ccIm.className(), 
                            new CustomColors( Color(80,80,100), Color(180,180,200) ) )
          << ccIm;
    board.saveEPS( "../output/original_CC.eps" );

    board.clear(Color::White);
    send1CellsVectorToBoard(bdryVect, board);
    board.saveEPS("../output/1cells.eps");
    board.clear(Color::White);

    cout << "-- Thresholding -";
    // Threshold the image and its inverse
    Binarizer b(image, minThresh, maxThresh);
    Binarizer bInv(imInv, minThresh, maxThresh);
    cout << " done." << endl;

    cout << "-- Computing the distance transform -";
    // Declare distance transform operators
    DTL2 dtl2(&image.domain(), &b, &Z2i::l2Metric);
    DTL2 dtl2Inv(&imInv.domain(), &bInv, &Z2i::l2Metric);
    cout << " done." << endl;

    // Compute the max value of the DT
    DTL2::Value maxDT2 = (*std::max_element(dtl2.constRange().begin(), dtl2.constRange().end()));
    DTL2::Value maxDT1 = (*std::max_element(dtl2Inv.constRange().begin(), dtl2Inv.constRange().end()));

    // Create GS DT images
    Image imGS = createImageFromDT(dtl2, maxDT2, false);
    Image imInvGS = createImageFromDT(dtl2Inv, maxDT1, false);

    cout << "-- Exporting preprocessed data -";
    // Save the image and its inverse
    saveImage(board, image, 0, 255, "../output/pre_processing/image.eps");
    saveImage(board, imInv, 0, 255, "../output/pre_processing/inverse.eps");
    // Save both DTs
    // +1 is necessary for some float values are superior to maxDT1 (which is an integer)
    saveImage(board, imInvGS, 0, maxDT1 + 1, "../output/pre_processing/dm_inv.eps");
    saveImage(board, imGS, 0, maxDT2 + 1, "../output/pre_processing/dm_or.eps");
    saveImage(board, addImages(imGS, imInvGS), 0, std::max(maxDT1 + 1, maxDT2 + 1),
              "../output/pre_processing/dm_merging.eps");
    cout << " done." << endl;

    // Process the DT
    // imGS: negative values, 1 -> -0.5
    // imInvGS: 1 -> 0.5
    processDT(imGS, true);
    processDT(imInvGS, false);

    // Add both DT
    Image imAddDTL = addImages(imInvGS, imGS);

    // Built path string
    string path;

    Image imNN = imAddDTL;
    Image imBil = imAddDTL;
    Image imBic = imAddDTL;
    std::vector<CC> cCVector;
    // Compute rotations depending on the given argument
    // Nearest neighbor
    if ((method == NEAREST_NEIGHBOR) || (method == ALL)) {
        cout << "-- Computing rotation using nearest neighbor -";
        imNN = rotateBackward(imNN, angle, NEAREST_NEIGHBOR);
        Image rotIm(getResizedDomain(imNN));

        // Thresholding operations
        thresholdDTImage(imNN, rotIm);
        ImagePGM pgm = thresholdToPGM(imNN);

        // Saves
        PGMWriter<ImagePGM>::exportPGM("../output/rotation_NN/NN.pgm", pgm);
        path = "../output/rotation_NN/NN.eps";
        saveImage(board, rotIm, 0, 255, path);
        cout << " done." << endl;
        cout << "   Output saved as " << path << endl;

        cout << "   Building corresponding cubical complex... ";
        KSpace kNN = initKSpace(rotIm.domain().lowerBound(), rotIm.domain().upperBound());
        CC cNN(kNN);
        getCCFromImage(rotIm, cNN, kNN);
        cCVector.push_back(cNN);
        cout << "- done." << endl;
    }

    // Bilinear interpolation
    if ((method == BILINEAR_INTERPOLATION) || (method == ALL)) {
        cout << "-- Computing rotation using bilinear interpolation -";
        imBil = rotateBackward(imBil, angle, BILINEAR_INTERPOLATION);
        Image rotIm(getResizedDomain(imBil));

        // Thresholding operations
        thresholdDTImage(imBil, rotIm);
        ImagePGM pgm = thresholdToPGM(imBil);

        // Saves
        PGMWriter<ImagePGM>::exportPGM("../output/rotation_BIL/bil.pgm", pgm);
        path = "../output/rotation_BIL/bil.eps";
        saveImage(board, rotIm, 0, 255, path);
        cout << " done." << endl;
        cout << "   Output saved as " << path << endl;

        cout << "   Building corresponding cubical complex... ";
        KSpace kBil = initKSpace(rotIm.domain().lowerBound(), rotIm.domain().upperBound());
        CC cBil(kBil);
        getCCFromImage(rotIm, cBil, kBil);
        cCVector.push_back(cBil);
        cout << "- done." << endl;
    }

    // Bicubic interpolation
    if ((method == BICUBIC_INTERPOLATION) || (method == ALL)) {
        cout << "-- Computing rotation using bicubic interpolation -";
        imBic = rotateBackward(imBic, angle, BICUBIC_INTERPOLATION);
        Image rotIm(getResizedDomain(imBic));
        // Thresholding operations
        thresholdDTImage(imBic, rotIm);

        ImagePGM pgm = thresholdToPGM(imBic);

        // Saves
        PGMWriter<ImagePGM>::exportPGM("../output/rotation_BIC/bic.pgm", pgm);
        path = "../output/rotation_BIC/bic.eps";
        saveImage(board, rotIm, 0, 255, path);
        cout << " done." << endl;
        cout << "   Output saved as " << path << endl;

        cout << "   Building corresponding cubical complex... ";
        KSpace kBic = initKSpace(rotIm.domain().lowerBound(), rotIm.domain().upperBound());
        CC cBic(kBic);
        getCCFromImage(rotIm, cBic, kBic);
        cCVector.push_back(cBic);
        cout << "- done." << endl;
    }

    cout << endl;
    cout << "============================================================" << endl;
    cout << endl;
    cout << "-- Topological informations:" << endl;
    cout << "   Original image: " << endl;
    cout << "       - Nb connected components foreground (b0)  : " << imObjects.size()        << endl;
    cout << "       - Nb connected components background (b1)  : " << imInvObjects.size() - 1 << endl;
    cout << "       - Surface of the foreground(nb of pixels)  : " << imObject.size()         << endl;
    cout << "       - Surface of the background(nb of pixels)  : " << imInvObject.size()      << endl;
    cout << "       - Nb 1-cells of the boundary (foreground)  : " << bdryVect.size()         << endl;
    cout << "       - Nb 1-cells of the boundary (background)  : " << bdryVectInv.size()      << endl;
    cout << "       - Euler's number (foreground)              : " << ccIm.euler()            << endl;
    cout << "       - Euler's number (background)              : " << ccImInv.euler()         << endl;
    cout << endl;

    cout << "   Rotated image: " << endl;
    if(cCVector.size() == 3)
    {
        int count = 0;
        for(auto cc : cCVector)
        {
            switch(count++)
            {
                case 0: 
                    cout << "       - Euler's number (Nearest Neighbor)        : ";
                    break;
                case 1: 
                    cout << "       - Euler's number (Bil interpolation)       : ";
                    break;
                case 2: 
                    cout << "       - Euler's number (Bic interpolation)       : ";
                    break;
            }     
            cout << cc.euler() << endl;
        }
    } else 
    {
        cout << "       - Euler's number (rotated)                 : " << cCVector[0].euler() << endl;
    }

    // Convert original image to grayscale
    imDTToGS(imAddDTL, -maxDT2, maxDT1);

    // Create a set for the image and its inverse
    Z2i::DigitalSet set(image.domain());
    Z2i::DigitalSet setInv(imInv.domain());
    DigitalSetInserter <Z2i::DigitalSet> inserter(set);
    DigitalSetInserter <Z2i::DigitalSet> inserterInv(setInv);
    DGtal::setFromImage(dtl2, inserter, 1, 135);
    DGtal::setFromImage(imInv, inserterInv, 1, 135);
}

int main(int argc, char **argv) {
    if (argc == 2)
        if (strcmp(argv[1], "help") == 0) {
            usage(true);
            return 0;
        }

    // Check args
    if (argc != 6 && argc != 4) {
        usage(false);
        return 0;
    }

    // Load image
    Image image = PGMReader<Image>::importPGM(argv[1]);

    int minThresh = (argc == 4) ? 254 : std::stoi(argv[4]);
    int maxThresh = (argc == 4) ? 255 : std::stoi(argv[5]);


    // process depending on the user's choice
    if (strcmp(argv[3], "bli") == 0)
        processImage(image, -std::stof(argv[2]), BILINEAR_INTERPOLATION, minThresh, maxThresh);
    else if (strcmp(argv[3], "nn") == 0)
        processImage(image, -std::stof(argv[2]), NEAREST_NEIGHBOR, minThresh, maxThresh);
    else if (strcmp(argv[3], "bic") == 0)
        processImage(image, -std::stof(argv[2]), BICUBIC_INTERPOLATION, minThresh, maxThresh);
    else if (strcmp(argv[3], "all") == 0)
        processImage(image, -std::stof(argv[2]), ALL, minThresh, maxThresh);
    else
        usage(false);

    return 0;
}
