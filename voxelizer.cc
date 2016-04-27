#include <limits>
#include <vector>
#include <cstddef>
#include <cstdarg>
#include <cstdlib>
#include <memory>
#include <algorithm>
#include <unistd.h>
#include <getopt.h>
#include <cmath>
#include <iostream>
#include <fstream>
#include <map>
#include <png.h>
#include <libgen.h>
#include <json-c/json.h>
#include "pdb.h"
#include "voxelizer.h"

using namespace std;

const char usage_format_string[] = "Usage: %s [-d dimensions] [-o output] input\nVoxelize a PDB file to voxel space of specified dimensions. Outputs 3D array of densities in JSON.";
char *output_filename = 0;
char *input_filename = 0;

void Usage() {
  Die(usage_format_string);
}

int main(int argc, char **argv) {
  base = basename(argv[0]);
  VMDPLUGIN_init();
  int x = 0, y = 0, z = 0, c;
  opterr = 0;
  vector<char *> filenames;
  vector<uint8_t> values;
  double radius = 0;

  static struct option long_options[] = {
    {"verbose", no_argument, 0, 'x'},
    {"dimensions", required_argument, 0, 'd'},
    {"output", required_argument, 0, 'o'},
    {"radius", required_argument, 0, 'r'},
    {"a-matrix", optional_argument, 0, 'a'},
    {"help", optional_argument, 0, 'h'},
    {0, 0, 0, 0}
  };

  int long_index = 0;

  char *a_matrix_filename = 0;
  bool output_a_matrix = false;

  while ((c = getopt_long(argc, argv, "vd:o:r:ha:", long_options, &long_index)) != -1) {
    switch (c) {
      case 'a':
        output_a_matrix = true;
        a_matrix_filename = optarg;
        break;
      case 'r':
        radius = atof(optarg);
        if (radius < 0) Die("Radius %f cannot be negative", radius);
        break;
      case 'd':
        strtodim(optarg, &x, &y, &z);
        if (x < 0 || y < 0 || z < 0) Die("Cannot supply negative dimensions");
        break;
      case 'o':
        output_filename = optarg;
        break;
      case 'h':
        Usage();
        break;
      case '?':
        if (optopt == 'a') {
          output_a_matrix = true;
        } else if (optopt == 'd' || optopt == 'r' || optopt == 'o') {
          Die("Option -%c requires an argument", optopt);
        } else if (isprint(optopt)) {
          Die("Unknown option '-%c'. Run '%s --help' to see options'", optopt, base);
        } else {
          Die("Unknown option character '\\x%x'. Run '%s --help' to see options'", optopt, base);
        }
        break;
      case ':':
        if (optopt == 'a') {
          output_a_matrix = true;
        }
        break;
    } 
  }
  while (optind < argc) {
    ParseFilename(argv[optind], filenames, values);
    optind++;
  }
  if (!filenames.size()) Die("Must supply input filename");
  MultiPDBVoxelizer mpv;
  for (unsigned i = 0; i < filenames.size(); ++i) {
    mpv.push_back(PDB::New(filenames[i], values[i]));
  }
  mpv.SetDimensions(x, y, z);
  mpv.CalculateSpan();
  mpv.SetRadius(radius);
  PNG<PNG_FORMAT_GA>::Pixel *voxels = mpv.Voxelize();
  PNG<PNG_FORMAT_GA>::Pixel *slice;
  string output_basename;
  if (output_filename) {
    output_basename = output_filename;
    for (int k = 0; k < z; ++k) {
      slice = ZSlice(voxels, x, y, k);
      PNG<PNG_FORMAT_GA> img (x, y, slice);
      if (!img.Write(output_basename + to_string(k) + ".png")) {
        Die("Failed to write %s", (output_basename + to_string(k) + ".png").c_str());
      }
    }
  }
  if (output_a_matrix) {
    size_t coords[] = { (size_t) x, (size_t) y, (size_t) z };
    shared_ptr<json_object> json_obj (CalculateHMMGroup(voxels, coords)->as_json_object(), &::json_object_put);
    if (a_matrix_filename) {
      const char *json = json_object_to_json_string(json_obj.get());
      ofstream out(a_matrix_filename);
      out << json;
      out.close();
    } else {
      printf("%s\n", json_object_to_json_string(json_obj.get()));
    }
  }
  delete[] voxels;
  return 0;
} 
