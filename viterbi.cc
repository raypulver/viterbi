#include "voxelizer.h"
#include <cmath>
#include <iostream>
#include <getopt.h>
#include <json-c/json.h>
#include "hmm.h"
using namespace std;
typedef enum _mode {
  GENERATE,
  SOLVE,
  ROTATE
} the_mode_t;

static void woop() {}
static the_mode_t mode;
static int dim = 0;
int main(int argc, char **argv) {
  if (argc > 1) {
    if (!strcmp(argv[argc - 2], "generate")) {
      mode = GENERATE;
			dim = atoi(argv[argc - 1]);
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
  //printf("%s\n\n", json_object_to_json_string(HMM2DToJsonObject(hmm)));
/*  hmm->Rotate(M_PI/4);
  printf("%s\n\n", json_object_to_json_string(HMM2DToJsonObject(hmm)));
  */
//  hmm->Rotate(M_PI/2);
  //printf("%s\n\n", json_object_to_json_string(HMM2DToJsonObject(hmm)));
  /*
  json_object *obj = json_object_new_object();
  json_object_object_add(obj, "original", HMM2DToJsonObject(hmm));
  hmm->Rotate(M_PI/2);
  json_object_object_add(obj, "rotated", HMM2DToJsonObject(hmm));
  */
  //printf("%s\n", json_object_to_json_string(obj));
  
  //Viterbi2DResult *result = Viterbi2DMax(hmm);
  hmm2d_t *hmmc = HMM2DToC(hmm);
	woop();
  double start = clock();
  viterbi2d_result_t *result = viterbi2d_max(hmmc);
  cout << clock() - start << endl;
  /*
  start = clock();
  viterbi2d_result_t *result_no_cache = viterbi2d_max_no_cache(hmmc);
  cout << clock() - start << endl;
   */
  //PNG<PNG_FORMAT_GA> *reconstruction = Reconstruct(hmm, result);
  //reconstruction->Write("result.png");
  } else if (mode == ROTATE) {
    cout << M_PI << endl;
    cout << sin(2*M_PI) << endl;
    cout << sin(4*M_PI) << endl;
  } else {
    WriteTriangle(dim, dim, "littlecircle.png");
  }
}
