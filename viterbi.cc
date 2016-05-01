#include "voxelizer.h"
#include <cmath>
#include <iostream>
#include <json-c/json.h>
using namespace std;
typedef enum _mode {
  GENERATE,
  SOLVE,
  ROTATE
} the_mode_t;

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
  //GenProjections(png, hmm->xobs, hmm->yobs);
  //printf("%s\n\n", json_object_to_json_string(HMM2DToJsonObject(hmm)));
/*  hmm->Rotate(M_PI/4);
  printf("%s\n\n", json_object_to_json_string(HMM2DToJsonObject(hmm)));
  */
//  hmm->Rotate(M_PI/2);
  //printf("%s\n\n", json_object_to_json_string(HMM2DToJsonObject(hmm)));
  hmm->Rotate(M_PI/2);
  json_object *obj = json_object_new_object();
  json_object_object_add(obj, "positive", HMM2DToJsonObject(hmm));
  json_object_object_add(obj, "negative", HMM2DToJsonObject(Calculate2DHMMReverse(png->GetPixelArray(), coords)));
  printf("%s\n", json_object_to_json_string(obj));
  /*
  Viterbi2DResult *result = Viterbi2DMax(hmm);
  PNG<PNG_FORMAT_GA> *reconstruction = Reconstruct(hmm, result);
  reconstruction->Write("result.png");
  */
  } else if (mode == ROTATE) {
    cout << M_PI << endl;
    cout << sin(2*M_PI) << endl;
    cout << sin(4*M_PI) << endl;
  } else {
    WriteTriangle(500, 500, "littlecircle.png");
  }
}
