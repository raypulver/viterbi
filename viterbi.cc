#include "voxelizer.h"
#include <cmath>
#include <iostream>
#include <json-c/json.h>
#include "hmm.h"
using namespace std;
typedef enum _mode {
  GENERATE,
  SOLVE,
  ROTATE
} the_mode_t;

void woop() {}
static the_mode_t mode;
int main(int argc, char **argv) {
  if (argc > 1) {
    if (!strcmp(argv[argc - 1], "generate")) {
      mode = GENERATE;
    } else if (!strcmp(argv[argc - 1], "solve")) {
      mode = SOLVE;
    } else {
      mode = ROTATE;
    }
  } else mode = SOLVE;
  if (mode == SOLVE) {
  PNG<PNG_FORMAT_GA> *png = PNG<PNG_FORMAT_GA>::FromFile("littlecircle.png");
  size_t coords[2] = { (size_t) png->GetWidth(), (size_t) png->GetHeight() };
  HMM2D::Ptr hmm = Calculate2DHMM<PNG<PNG_FORMAT_GA>::Pixel>((PNG<PNG_FORMAT_GA>::Pixel *) png->GetPixelArray(), coords);
  GenProjections(png, hmm->xobs, hmm->yobs);
  cout << "Generated!" << endl;
  hmm2d_t *hmmc = HMM2DToC(hmm);
  woop();
  double start = clock();
  viterbi2d_result_t *result = viterbi2d_max(hmmc);
  cout << clock() - start << endl;
  } else if (mode == ROTATE) {
    cout << M_PI << endl;
    cout << sin(2*M_PI) << endl;
    cout << sin(4*M_PI) << endl;
  } else {
    WriteTriangle(256, 256, "littlecircle.png");
  }
}
