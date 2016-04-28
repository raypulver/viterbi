#ifndef VOXELIZER_H
#define VOXELIZER_H
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

using namespace std;

template <typename T> T log(T, string);

static inline void strtodim(char *str, int *x, int *y, int *z) {
  char *next;
  *x = atoi(strtok(str, "x"));
  if (!(next = strtok(NULL, "x"))) {
    *y = *x;
    *z = *x;
    return;
  }
  *y = atoi(next);
  *z = atoi(strtok(NULL, "x"));
}

static inline bool InSphere(double x, double h, double y, double k, double z, double l, double r) {
  return pow(x - h, 2) + pow(y - k, 2) + pow(z - l, 2) <= pow(r, 2);
}

template <typename T> T find_max(T *vals, size_t sz);

template <typename T> void triple_min_max(T *, size_t, T *, T *, T *, T *, T *, T *);

template <typename T> T sum(T *, size_t);

template <typename T> static inline T varmax(T a) { return a; }

template <typename T, typename... Types> static inline T varmax(T a, Types ...rest) {
  T last = varmax(rest...);
  if (last > a) return last;
  else return a;
}

template <int format> class PNG {
  void *buffer;
  void *colormap;
  png_image img;
  int width;
  int height;
  public:
    struct Pixel {
      uint8_t g;
      uint8_t a;
      uint8_t GetValue();
    };
    PNG(int, int, void *);
    ~PNG();
    int Write(string);
    int Write(const char *filename);
};
  
template <typename T> T *ZSlice(T *values, size_t x, size_t y, size_t z);

#define MAX_ERROR_FORMAT_STRING_SIZE (1 << 16)

extern char *base;

template <typename T, typename... Args> void Die(T, Args...);

template <typename T> void Die(T);

void Usage();

class PDB {
  friend class MultiPDBVoxelizer;
  molfile_timestep_t ts;
  molfile_atom_t *atoms;
  int natoms;
  int optflags;
  float xmin, xmax, ymin, ymax, zmin, zmax;
  void *handle;
  uint8_t density;
  public:
    typedef shared_ptr<PDB> Ptr;
    static Ptr New(char *, uint8_t);
    PDB(char *, uint8_t);
    ~PDB();
};

class MultiPDBVoxelizer {
  float xmin, xmax, ymin, ymax, zmin, zmax, xdiff, ydiff, zdiff, xadj, yadj, zadj, maxdim;
  double xratio, yratio, zratio, step, radius, vradius;
  int x, y, z, a, v, maxpxl;
  int xoffset, yoffset, zoffset;
  vector<PDB::Ptr> pdbs;
  public:
    void SetRadius(double r);
    void SetDimensions(int i, int j, int k);
    void push_back(PDB::Ptr);
    void CalculateSpan();
    PNG<PNG_FORMAT_GA>::Pixel *Voxelize();
};

template <typename T> void ParseFilename(char *fn, vector<char *> &filenames, vector<T> &values);

json_object *NewDoubleOrInt(double d);

struct HMM {
  typedef pair<uint8_t, size_t> State;
  typedef pair<State, State> Transition;
  typedef uint8_t Observation;
  typedef shared_ptr<HMM> Ptr;
  static Ptr New();
  vector<Observation> obs;
  vector<State> states;
  vector<double> matrix;
  vector<double> emit;
  vector<double> initial;
  json_object *as_json_object();
};

struct HMMGroup {
  typedef shared_ptr<HMMGroup> Ptr;
  static Ptr New();
  HMM::Ptr xpos, xneg, ypos, yneg, zpos, zneg;
  json_object *as_json_object();
};

void IncreaseOrDecrease(size_t &n, uint8_t sign);

template <typename T> HMM::Ptr CalculateHMM(T *items, size_t *coords, uint8_t fix, uint8_t sign);

template <typename T> HMMGroup::Ptr CalculateHMMGroup(T *items, size_t *dimensions);

#define INDEX(it) (distance(it.begin(), it))

struct ViterbiResult {
  ViterbiResult(ViterbiResult *plast, HMM::State *pptr, double pprobability);
  ~ViterbiResult();
  ViterbiResult *last;
  HMM::State *ptr;
  double probability;
};

struct HMM2D {
  enum class Direction {
    X, Y
  };
  typedef shared_ptr<HMM2D> Ptr;
  static Ptr New();
  vector<double> &GetTransition(Direction);
  vector<double> &GetInitial(Direction);
  typedef uint8_t PartialState;
  typedef vector<PartialState> State;
  vector<PartialState> states;
  vector<double> xtransition;
  vector<double> xinitial;
  vector<double> ytransition;
  vector<double> yinitial;
};

struct Permutation {
  vector<Permutation *> last;
  HMM2D::PartialState state;
  double probability;
  ~Permutation();
};

struct Viterbi2DResult {
  ~Viterbi2DResult();
  Viterbi2DResult *last;
  HMM2D::State x;
  HMM2D::State y;
  HMM2D::Direction direction;
  double probability;
};

ViterbiResult *Viterbi(HMM::Ptr m, const vector<HMM::Observation> &, size_t len, HMM::State *s, ViterbiResult *last);

ViterbiResult *Viterbi(HMM::Ptr m, const vector<HMM::Observation> &, size_t len, HMM::State *s, ViterbiResult *last);

ViterbiResult *ViterbiMax(HMM::Ptr m, vector<HMM::Observation> &);

ViterbiResult *ViterbiMax(HMM::Ptr m, const vector<HMM::Observation> &, size_t len, ViterbiResult *last);

json_object *StateToJsonObject(HMM::State s);

void PrintViterbiResult(ViterbiResult *vr);

Permutation *EM(HMM2D::Ptr a, HMM2D::Direction d, size_t len, HMM2D::PartialState s);
Permutation *EMMax(HMM2D::Ptr a, HMM2D::Direction d, size_t len);
#endif
