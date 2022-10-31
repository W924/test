#include "Detector.h"
#include <strings.h>

using namespace std;

int main(int argc, char *argv[]) {
    // Clion 编译文件在 cmake-build-debug 中，所以是".."
    string infileFolder = "../input/Air20Quality20Data/Beijing/";
    string outfileFolder = "../output/Air20Quality20Data/Beijing/";

    int single_file_test = strcasecmp(argv[1], "-singleFile");
    int all_file_test = strcasecmp(argv[1], "-allFile");

    Detector detector;
    if(!single_file_test) {
        string infileName = argv[2];
        detector.singleFileTest(infileFolder, outfileFolder, infileName);
    }

    if(!all_file_test) {
        detector.allFileTest(infileFolder, outfileFolder);
    }

    return 0;
}