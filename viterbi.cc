#include "voxelizer.h"
#include <iostream>
using namespace std;
typedef enum _mode {
  GENERATE,
  SOLVE
} the_mode_t;

static the_mode_t mode;
int main(int argc, char **argv) {
  if (argc > 1) {
    if (!strcmp(argv[argc - 1], "generate")) {
      mode = GENERATE;
    } else if (!strcmp(argv[argc - 1], "solve")) {
      mode = SOLVE;
    }
  } else mode = SOLVE;
  if (mode == SOLVE) {
  PNG<PNG_FORMAT_GA> *png = PNG<PNG_FORMAT_GA>::FromFile("littlecircle.png");
  size_t coords[2] = { (size_t) png->GetWidth(), (size_t) png->GetHeight() };
  HMM2D::Ptr hmm = Calculate2DHMM<PNG<PNG_FORMAT_GA>::Pixel>((PNG<PNG_FORMAT_GA>::Pixel *) png->GetBuffer(), coords);
  GenProjections(png, hmm->xobs, hmm->yobs);
  Viterbi2DResult *result = Viterbi2DMax(hmm);
  PNG<PNG_FORMAT_GA> *reconstruction = Reconstruct(hmm, result);
  reconstruction->Write("result.png");
  } else {
    WriteCircle(8, 8, "littlecircle.png");
  }
}
