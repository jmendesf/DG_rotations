#include <boost/algorithm/minmax_element.hpp>
#include "DGtal/kernel/sets/DigitalSetInserter.h"
#include "DGtal/io/readers/PGMReader.h"
#include "Board/Point.h"
#include "../include/tools.h"
#include "../include/images.h"
#include "../include/rotations.h"
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

void processImage(Image &image, float angle, INTERPOLATION_METHOD method, int minThresh, int maxThresh) {
    // Create board and set its size
    Board2D board;
    board << image.domain();
    board.clear();
    // copy the image then inverse it
    Image imInv = image;
    inverseImage(imInv);

    saveImage(board, image, 0, 255, "../output/dt_set/im.eps");

    // Import image to the board
    Display2DFactory::drawImage<Gray>(board, image, (unsigned char) 0, static_cast<unsigned char>(255));

    // Threshold the image
    Binarizer b(image, minThresh, maxThresh);
    Binarizer bInv(imInv, minThresh, maxThresh);

    DTL2 dtl2(&image.domain(), &b, &Z2i::l2Metric);
    DTL2 dtl2Inv(&imInv.domain(), &bInv, &Z2i::l2Metric);

    // Compute the max value of the DT
    DTL2::Value maxDT2 = (*std::max_element(dtl2.constRange().begin(), dtl2.constRange().end()));
    DTL2::Value maxDT1 = (*std::max_element(dtl2Inv.constRange().begin(), dtl2Inv.constRange().end()));

    // Create GS DT images
    Image imGS = createImageFromDT(dtl2, maxDT2, false);
    Image imInvGS = createImageFromDT(dtl2Inv, maxDT1, false);

    // Process the DT
    // imGS: negative values, 1 -> -0.5
    // imInvGS: 1 -> 0.5
    processDT(imGS, true);
    processDT(imInvGS, false);

    board.clear();
    Display2DFactory::drawImage<Gray>(board, dtl2Inv, (DTL2::Value) 0,
                                      (DTL2::Value) maxDT1);

    // Add both DT
    Image imAddDTL = addImages(imInvGS, imGS);
    saveImage(board, imAddDTL, -maxDT2, maxDT1, "../DM_original.eps");


    cout << endl;

    // Built path string
    string path;

    Image imNN = imAddDTL;
    Image imBil = imAddDTL;
    Image imBic = imAddDTL;

    // Compute rotations depending on the given argument
    // Nearest neighbor
    if ((method == NEAREST_NEIGHBOR) || (method == ALL)) {
        cout << "-- Computing rotation using nearest neighbor -";
        imNN = rotateBackward(imNN, angle, NEAREST_NEIGHBOR);
        Image rotIm(getResizedDomain(imNN));

        thresholdDTImage(imNN, rotIm);
        ImagePGM pgm = thresholdToPGM(imNN);
        PGMWriter<ImagePGM>::exportPGM("../output/rotation_NN/NN.pgm", pgm);

        path = "../output/rotation_NN/NN.eps";
        saveImage(board, rotIm, 0, 255, path);
        cout << " done." << endl;
        cout << "   Output saved as " << path << endl;
    }

    // Bilinear interpolation
    if ((method == BILINEAR_INTERPOLATION) || (method == ALL)) {
        cout << "-- Computing rotation using bilinear interpolation -";
        imBil = rotateBackward(imBil, angle, BILINEAR_INTERPOLATION);

        Image rotIm(getResizedDomain(imBil));
        thresholdDTImage(imBil, rotIm);
        ImagePGM pgm = thresholdToPGM(imBil);

        PGMWriter<ImagePGM>::exportPGM("../output/rotation_BIL/bil.pgm", pgm);
        path = "../output/rotation_BIL/bil.eps";
        saveImage(board, rotIm, 0, 255, path);
        cout << " done." << endl;
        cout << "   Output saved as " << path << endl;
    }

    // Bicubic interpolation
    if ((method == BICUBIC_INTERPOLATION) || (method == ALL)) {
        cout << "-- Computing rotation using bicubic interpolation -";
        imBic = rotateBackward(imBic, angle, BICUBIC_INTERPOLATION);

        Image rotIm(getResizedDomain(imBic));
        thresholdDTImage(imBic, rotIm);
        ImagePGM pgm = thresholdToPGM(imBic);

        PGMWriter<ImagePGM>::exportPGM("../output/rotation_BIC/bic.pgm", pgm);
        path = "../output/rotation_BIC/bic.eps";
        saveImage(board, rotIm, 0, 255, path);
        cout << " done." << endl;
        cout << "   Output saved as " << path << endl;
    }

    cout << endl;

    // Convert original image to grayscale
    imDTToGS(imAddDTL, -maxDT2, maxDT1);

    // Create a set for the image and its inverse
    Z2i::DigitalSet set(image.domain());
    Z2i::DigitalSet setInv(imInv.domain());
    DigitalSetInserter<Z2i::DigitalSet> inserter(set);
    DigitalSetInserter<Z2i::DigitalSet> inserterInv(setInv);
    DGtal::setFromImage(dtl2, inserter, 1, 135);
    DGtal::setFromImage(imInv, inserterInv, 1, 135);

    // Save output
    saveImage(board, imAddDTL, 0, 255, "../output/dt_set/im_add_DT.eps");
    saveImage(board, imGS, 0, 255, "../output/dt_set/im_GS_DT.eps");
    saveImage(board, imInvGS, 0, 255, "../output/dt_set/im_inv_GS_DT");
    saveSet(board, set, "../output/dt_set/set.eps");
    saveSet(board, setInv, "../output/dt_set/set_inv.eps");
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
