#include "voxelizer.h"
#include <cmath>
#include <iostream>
#include <getopt.h>
#include <json-c/json.h>
#include "hmm.h"
#include "cache.h"
using namespace std;
typedef enum _mode {
  GENERATE,
  SOLVE,
  ROTATE
} the_mode_t;

void woop() {}
static the_mode_t mode;
static int dim = 0;
bool reconstructit = false;
int main(int argc, char **argv) {
  static struct option long_options[] = {
    {"reconstruct", no_argument, 0, 'r'},
    {0, 0, 0, 0}
  };
  int long_index = 0;
  int c;
  while ((c = getopt_long(argc, argv, "r", long_options, &long_index)) != -1) {
    switch (c) {
      case 'r':
        reconstructit = true;
        break;
    }
  }
  while (optind < argc) {
    if (!strcmp(argv[optind], "generate")) {
      mode = GENERATE;
      optind++;
      if (optind == argc) { cerr << "viterbi: must supply a dimension to generate" << endl; exit(1); }
      dim = atoi(argv[optind]);
    } else if (!strcmp(argv[optind], "solve")) {
      mode = SOLVE;
    } else if (!strcmp(argv[optind], "rotate")) {
      mode = ROTATE;
    }
    optind++;
  } 
  if (mode == SOLVE) {
  PNG<PNG_FORMAT_GA> *png = PNG<PNG_FORMAT_GA>::FromFile("out.png");
  size_t coords[2] = { (size_t) png->GetWidth(), (size_t) png->GetHeight() };
  HMM2D::Ptr hmm = Calculate2DHMM<PNG<PNG_FORMAT_GA>::Pixel>((PNG<PNG_FORMAT_GA>::Pixel *) png->GetPixelArray(), coords);
  GenProjectionsMapped(png, hmm);
  hmm2d_t *hmmc = HMM2DToC(hmm);
  double start = clock();
  cache_t *cache;
  viterbi2d_result_t *result = viterbi2d_max(hmmc, &cache);
  if (reconstructit) reconstruct(hmmc, cache, "reconstruction.png");
  check_cache(cache);
  cache_free(cache);
  woop();
  cout << clock() - start << endl;
  } else if (mode == ROTATE) {
  PNG<PNG_FORMAT_GA> *png = PNG<PNG_FORMAT_GA>::FromFile("out.png");
  size_t coords[2] = { (size_t) png->GetWidth(), (size_t) png->GetHeight() };
  HMM2D::Ptr hmm = Calculate2DHMM<PNG<PNG_FORMAT_GA>::Pixel>((PNG<PNG_FORMAT_GA>::Pixel *) png->GetPixelArray(), coords);
  hmm2d_t *hmmc = HMM2DToC(hmm);
  print_hmm(hmmc);
  matrix_t *rotated = rotate(hmmc, M_PI/4);
  for (size_t i = 0; i < rotated->x; ++i) {
    cout << "[ ";
    for (size_t j = 0; j < rotated->y; ++j) {
      cout << *matrix_el(rotated, i, j) << " ";
    }
    cout << "]" << endl;
  }
  } else {
    WriteTriangle(dim, dim, "out.png");
    PNG<PNG_FORMAT_GA> *png = PNG<PNG_FORMAT_GA>::FromFile("out.png");
    size_t coords[2] = { (size_t) png->GetWidth(), (size_t) png->GetHeight() };
    HMM2D::Ptr hmm = Calculate2DHMM<PNG<PNG_FORMAT_GA>::Pixel>((PNG<PNG_FORMAT_GA>::Pixel *) png->GetPixelArray(), coords);
    hmm2d_t *hmmc = HMM2DToC(hmm);
    print_hmm(hmmc);
    delete png;
  }
}
