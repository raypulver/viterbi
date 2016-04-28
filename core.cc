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

template <typename T> T log(T arg, string str) {
  cout << str << ": " << arg << endl;
  return arg;
}

template <typename T> T find_max(T *vals, size_t sz) {
  size_t i;
  T max = numeric_limits<T>::min();
  for (i = 0; i < sz; ++i) {
    if (vals[i] > max) max = vals[i];
  }
  return max;
}

template <typename T> void triple_min_max(T *coords, size_t n, T *xmin, T *xmax, T *ymin, T *ymax, T *zmin, T *zmax) {
  size_t i;
  *xmin = numeric_limits<T>::max();
  *ymin = numeric_limits<T>::max();
  *zmin = numeric_limits<T>::max();
  *xmax = numeric_limits<T>::min();
  *ymax = numeric_limits<T>::min();
  *zmax = numeric_limits<T>::min(); 
  for (i = 0; i < n; ++i) {
    if (coords[i*3] < *xmin) *xmin = coords[i*3];
    if (coords[i*3] > *xmax) *xmax = coords[i*3];
    if (coords[i*3 + 1] < *ymin) *ymin = coords[i*3 + 1];
    if (coords[i*3 + 1] > *ymax) *ymax = coords[i*3 + 1];
    if (coords[i*3 + 2] < *zmin) *zmin = coords[i*3 + 2];
    if (coords[i*3 + 2] > *zmax) *zmax = coords[i*3 + 2];
  }
}

template <typename T> T sum(T *values, size_t sz) {
  size_t i;
  T sum = 0;
  for (i = 0; i < sz; ++i) {
    sum += values[i];
  }
  return sum;
}

template <int format> uint8_t PNG<format>::Pixel::GetValue() {
  return g;
}

template <int format> PNG<format>::PNG(int x, int y, void *buf) : buffer(buf) {
  img.width = x;
  img.height = y;
  img.version = PNG_IMAGE_VERSION;
  img.opaque = nullptr;
  img.format = format;
  img.flags = 0;
  img.colormap_entries = 0;
  colormap = nullptr;
}

template <int format> PNG<format>::~PNG() {
  if (buffer) free(buffer);
  if (colormap) free(colormap);
}

template <int format> int PNG<format>::Write(string filename) {
  return Write(filename.c_str());
}

template <int format> int PNG<format>::Write(const char *filename) {
  return png_image_write_to_file(&img, filename, 0, buffer, 0, colormap);
}
  
template <typename T> T *ZSlice(T *values, size_t x, size_t y, size_t z) {
  size_t i, j;
  T *retval = new T[x*y];
  for (i = 0; i < x; ++i) {
    for (j = 0; j < y; ++j) {
      retval[i*x + j] = values[i*x*y + j*y + z];
    }
  }
  return retval;
}
  
char *base;

template <typename T, typename... Args> void Die(T fmt, Args ...args) {
  char buf[MAX_ERROR_FORMAT_STRING_SIZE];
  memset(buf, 0, sizeof(buf));
  sprintf(buf, fmt, args...);
  string err = base;
  err += ": ";
  err += buf;
  cerr << err << endl;
  exit(1);
}

template <typename T> void Die(T arg) {
  Die("%s", arg);
}

PDB::Ptr PDB::New(char *filename, uint8_t density) { return PDB::Ptr(new PDB(filename, density)); }

PDB::PDB(char *filename, uint8_t dens) : density(dens) {
  handle = plugin.open_file_read(filename, "pdb", &natoms);
  if (!handle) Die("<VMDPLUGIN> open_file_read(\"%s\") failed.", filename);
  atoms = new molfile_atom_t[natoms];
  ts.coords = new float[natoms*3];
  plugin.read_structure(handle, &optflags, atoms);
  plugin.read_next_timestep(handle, natoms, &ts);
  triple_min_max(ts.coords, natoms, &xmin, &xmax, &ymin, &ymax, &zmin, &zmax);
}

PDB::~PDB() {
  plugin.close_file_read(handle);
  delete[] atoms;
  delete[] ts.coords; 
}

void MultiPDBVoxelizer::SetRadius(double r) {
  radius = r;
  vradius = step*r;
  step *= (double) (maxpxl + r*6)/maxpxl;
  vradius = r*step;
  xoffset += r*3;
  yoffset += r*3;
  zoffset += r*3;
}

void MultiPDBVoxelizer::SetDimensions(int i, int j, int k) { x = i, y = j, z = k, v = x*y*z, a = x*y; }
void MultiPDBVoxelizer::push_back(PDB::Ptr p) { pdbs.push_back(p); }
void MultiPDBVoxelizer::CalculateSpan() {
  xmin = numeric_limits<float>::max();
  ymin = numeric_limits<float>::max();
  zmin = numeric_limits<float>::max();
  xmax = numeric_limits<float>::min();
  ymax = numeric_limits<float>::min();
  zmax = numeric_limits<float>::min();
  for (auto it = pdbs.begin(); it != pdbs.end(); it++) {
    if (((*it)->xmin) < xmin) xmin = (*it)->xmin;
    if (((*it)->ymin) < ymin) ymin = (*it)->ymin;
    if (((*it)->zmin) < zmin) zmin = (*it)->zmin;
    if (((*it)->xmax) > xmax) xmax = (*it)->xmax;
    if (((*it)->ymax) > ymax) ymax = (*it)->ymax;
    if (((*it)->zmax) > zmax) zmax = (*it)->zmax;
  }
  xdiff = xmax - xmin;
  ydiff = ymax - ymin;
  zdiff = zmax - zmin;
  maxdim = varmax(xdiff, ydiff, zdiff);
  xratio = xdiff/maxdim;
  yratio = ydiff/maxdim;
  zratio = zdiff/maxdim;
  if (zratio == 1.0) maxpxl = z;
  if (yratio == 1.0) maxpxl = y;
  if (xratio == 1.0) maxpxl = x;
  xadj = xmin/xratio;
  yadj = ymin/yratio;
  zadj = zmin/zratio;
  step = maxdim/maxpxl;
  xoffset = xdiff*(1 - xratio)*step/2;
  yoffset = ydiff*(1 - yratio)*step/2;
  zoffset = xdiff*(1 - zratio)*step/2;
}

PNG<PNG_FORMAT_GA>::Pixel *MultiPDBVoxelizer::Voxelize() {
  int density = 0, maxdens = 0, i, j, k, l;
  double mincoords[3];
  double maxcoords[3];
  double center[3];
  double centercoords[3];
  PDB *winner = nullptr;
  PNG<PNG_FORMAT_GA>::Pixel *retval = new PNG<PNG_FORMAT_GA>::Pixel[v];
  for (i = 0; i < x; ++i) {
    for (j = 0; j < y; ++j) {
      for (k = 0; k < z; ++k) {
        winner = nullptr;
        maxdens = 0;
        for (auto it = pdbs.begin(); it != pdbs.end(); it++) {
          density = 0;
          for (l = 0; l < (*it)->natoms; l++) {
            center[0] = (*it)->ts.coords[l*3];
            center[1] = (*it)->ts.coords[l*3 + 1];
            center[2] = (*it)->ts.coords[l*3 + 2];
            mincoords[0] = xadj + (double) (i - xoffset)*step;
            centercoords[0] = mincoords[0] + step/2;
            maxcoords[0] = mincoords[0] + step;
            mincoords[1] = yadj + (double) (j - yoffset)*step;
            centercoords[1] = mincoords[1] + step/2;
            maxcoords[1] = mincoords[1] + step;
            mincoords[2] = zadj + (double) (k - zoffset)*step;
            centercoords[2] = mincoords[2] + step/2;
            maxcoords[2] = mincoords[2] + step;
            if ((center[0] >= mincoords[0] &&
              center[0] < maxcoords[0] &&
              center[1] >= mincoords[1] &&
              center[1] < maxcoords[1] &&
              center[2] >= mincoords[2] &&
              center[2] < maxcoords[2]) || 
              InSphere(centercoords[0], center[0], centercoords[1], center[1], centercoords[2], center[2], vradius*get_pte_vdw_radius((*it)->atoms[l].atomicnumber))) {
              ++density;
            }
          }
          if (density > maxdens) {
            maxdens = density;
            winner = it->get();
          }
        }
        if (!winner) { retval[i*a + j*y + k] = {0, 0}; }
        else { retval[i*a + j*y + k] = { winner->density, 0xff }; }
      }
    }
  }
  return retval;
}

template <typename T> void ParseFilename(char *fn, vector<char *> &filenames, vector<T> &values) {
  char *ptr = fn + strlen(fn);
  filenames.push_back(fn);
  for (; ptr >= fn; ptr--) {
    if (*ptr == ':') {
      *ptr = '\0';
      values.push_back(atoi(ptr + 1));
      return;
    }
  }
  values.push_back(100);
}

json_object *NewDoubleOrInt(double d) {
  if (!d) return json_object_new_int(0);
  return json_object_new_double(d);
}

HMM::Ptr HMM::New() { return HMM::Ptr(new HMM); }
json_object *HMM::as_json_object() {
  json_object *retval = json_object_new_object();
  json_object *state_index = json_object_new_array();
  for (auto it = states.begin(); it != states.end(); it++) {
    json_object *state = json_object_new_object();
    json_object_object_add(state, "density", json_object_new_int(it->first));
    json_object_object_add(state, "duration", json_object_new_int(it->second));
    json_object_array_add(state_index, state);
  }
  json_object_object_add(retval, "states", state_index);
  json_object *initial_states = json_object_new_array();
  json_object_object_add(retval, "initial", initial_states);
  for (auto it = initial.begin(); it != initial.end(); it++) {
    json_object_array_add(initial_states, NewDoubleOrInt(*it));
  }
  json_object *matrix_array = json_object_new_array();
  json_object_object_add(retval, "matrix", matrix_array); 
  for (size_t i = 0; i < states.size(); ++i) {
    json_object *row = json_object_new_array();
    json_object_array_add(matrix_array, row);
    for (size_t j = 0; j < states.size(); ++j) {
      json_object_array_add(row, NewDoubleOrInt(matrix[i*states.size() + j]));
    }
  }
  return retval;
}

HMMGroup::Ptr HMMGroup::New() { return HMMGroup::Ptr(new HMMGroup); }
json_object *HMMGroup::as_json_object() {
  json_object *retval = json_object_new_object();
  json_object_object_add(retval, "xpos", xpos->as_json_object());
  json_object_object_add(retval, "xneg", xneg->as_json_object());
  json_object_object_add(retval, "ypos", ypos->as_json_object());
  json_object_object_add(retval, "yneg", yneg->as_json_object());
  json_object_object_add(retval, "zpos", zpos->as_json_object());
  json_object_object_add(retval, "zneg", zneg->as_json_object());
  return retval;
}

void IncreaseOrDecrease(size_t &n, uint8_t sign) {
  if (sign) n--;
  else n++;
}

template <typename T> HMM::Ptr CalculateHMM(T *items, size_t *coords, uint8_t fix, uint8_t sign) {
  size_t permutecoords[3];
  size_t multiplier[3];
  size_t idx = 0;
  for (size_t m = 0; m < 3; ++m) {
    if (fix == m) continue;
    permutecoords[idx] = coords[m];
    if (m == 0) multiplier[idx] = coords[0]*coords[1];
    if (m == 1) multiplier[idx] = coords[1];
    if (m == 2) multiplier[idx] = 1;
    idx++; 
  }
  switch (fix) {
    case 0:
      multiplier[idx] = coords[0]*coords[1];
      break;
    case 1:
      multiplier[idx] = coords[1];
      break;
    case 2:
      multiplier[idx] = 1;
  }
  permutecoords[idx] = coords[fix];
  HMM::Ptr retval = HMM::New();
  uint8_t last_depth = 0;
  uint8_t depth = 0;
  size_t duration = 0;
  bool starting = true;
  bool initial = true;
  HMM::State last_state;
  vector<HMM::State> initial_states;
  vector<HMM::Transition> transitions;
  for (size_t i = 0; i < permutecoords[0]; ++i) {
    for (size_t j = 0; j < permutecoords[1]; ++j) {
      for (size_t k = (sign == 1 ? permutecoords[2] - 1: 0); (sign == 1 ? k != ((size_t) 0 - (size_t) 1) : k < permutecoords[2]); IncreaseOrDecrease(k, sign)) {
        depth = items[i*multiplier[0] + j*multiplier[1] + multiplier[2]*k].GetValue();
        if (!starting && (last_depth != depth)) {
          if (initial)
            initial_states.push_back(HMM::State(last_depth, duration));
          else transitions.push_back(HMM::Transition(last_state, HMM::State(last_depth, duration)));
          retval->states.push_back(HMM::State(last_depth, duration));
          last_state = HMM::State(last_depth, duration);
          duration = 0;
          initial = false;
        }
        duration++;
        last_depth = depth;
        starting = false;
      }
      if (initial) {
        initial_states.push_back(HMM::State(depth, duration));
      }
      else transitions.push_back(HMM::Transition(last_state, HMM::State(depth, duration)));
      retval->states.push_back(HMM::State(depth, duration));
      memset(&last_state, 0, sizeof(HMM::State)); 
      duration = 0;
      last_depth = 0;
      depth = 0;
      initial = true;
      starting = true;
    }
  }
  sort(retval->states.begin(), retval->states.end());
  auto it = unique(retval->states.begin(), retval->states.end(), [&] (HMM::State a, HMM::State b) -> bool {
    return ((a.first == b.first) && (b.second == a.second));
  });
  retval->states.resize(distance(retval->states.begin(), it));
  retval->matrix.resize(retval->states.size() * retval->states.size());
  map<HMM::State, size_t> state_map;
  for (auto it = retval->states.begin(); it != retval->states.end(); it++) {
    state_map[*it] = distance(retval->states.begin(), it);
  }
  for (auto it = transitions.begin(); it != transitions.end(); it++) {
    retval->matrix[state_map[it->first]*retval->states.size() + state_map[it->second]]++;
  }
  for (size_t i = 0; i < retval->states.size(); ++i) {
    double sum = 0;
    for (size_t j = 0; j < retval->states.size(); ++j) {
      sum += retval->matrix[i*retval->states.size() + j];
    }
    if (sum) for (size_t j = 0; j < retval->states.size(); ++j) {
      retval->matrix[i*retval->states.size() + j] /= sum;
    }
  }
  retval->initial.resize(retval->states.size());
  for (auto it = initial_states.begin(); it != initial_states.end(); it++) {
    retval->initial[state_map[*it]]++;
  }
  double sum = 0;
  for (size_t i = 0; i < retval->states.size(); ++i) {
    sum += retval->initial[i];
  }
  if (sum) for (size_t i = 0; i < retval->states.size(); ++i) {
    retval->initial[i] /= sum;
  }
  return retval;
}

template <typename T> HMMGroup::Ptr CalculateHMMGroup(T *items, size_t *dimensions) {
  HMMGroup::Ptr retval = HMMGroup::New();
  retval->xpos = CalculateHMM(items, dimensions, 0, 0);
  retval->xneg = CalculateHMM(items, dimensions, 0, 1);
  retval->ypos = CalculateHMM(items, dimensions, 1, 0);
  retval->yneg = CalculateHMM(items, dimensions, 1, 1);
  retval->zpos = CalculateHMM(items, dimensions, 2, 0);
  retval->zneg = CalculateHMM(items, dimensions, 2, 1);
  return retval;
}

Permutation *EMStartingWith(HMM2D::Ptr a, HMM2D::Direction d, size_t len, HMM2D::PartialState s) {
  vector<double> backup = a->GetInitial(d);
  fill(a->GetInitial(d).begin(), a->GetInitial(d).end(), 0);
  a->GetInitial(d)[s] = 1;
  Permutation* retval = EMMax(a, d, len + 1);
  Permutation* starting = retval->last[0];
  a->GetInitial(d) = backup;
  return starting;
}

void ForeachPermutation(Permutation *p, function<bool(HMM2D::State &, size_t, double)> fn, HMM2D::State &state, size_t &idx) {
  Permutation *tmp;
  if (p->last.size() == 0) {
    state.push_back(p->state);
    fn(state, idx, p->probability);
    idx++;
    state.resize(state.size() - 1);
  } else {
    if (p->state != 255) state.push_back(p->state);
    for (auto it = p->last.begin(); it != p->last.end(); it++) {
      ForeachPermutation(*it, fn, state, idx);
    }
    if (p->state != 255) state.resize(state.size() - 1);
    for (auto it = p->last.begin(); it != p->last.end(); it++) {
      ForeachPermutation(*it, fn, state, idx);
    }
    if (p->state != 255) state.resize(state.size() - 1);
  }
}

void ForeachPermutation(Permutation *p, function<bool(HMM2D::State &, size_t, double)> fn) {
  HMM2D::State permute;
  size_t i = 0;
  ForeachPermutation(p, fn, permute, i);
}

double SumThe2DState(HMM2D::State &state) {
  size_t sum = 0;
  for (auto s : state) {
    sum += s;
  }
  return sum;
}

double SumThe2DStateExceptFirst(HMM2D::State &state) {
  size_t sum = 0;
  HMM2D::State::const_iterator it = state.cbegin();
  it++;
  for (; it != state.cend(); it++) {
    sum += *it;
  }
  return sum;
}
    
double Prob(HMM2D::Ptr a, HMM2D::Direction d, size_t obs, HMM2D::State &state, size_t len) {
  double total = 0;
  size_t sum = SumThe2DState(state);
  int o = obs;
  o -= sum;
  if (o < 0) return 0;
  vector<double> initial = a->GetInitial(d);
  fill(a->GetInitial(d).begin(), a->GetInitial(d).end(), 0);
  a->GetInitial(d)[state.back()] = 1;
  Permutation *permute = EMMax(a, d, len - state.size());
  ForeachPermutation(permute, [&] (HMM2D::State &state, size_t i, double probability) -> bool {
    if (SumThe2DStateExceptFirst(state) == o) { total += probability; }
    return true;
  });
  delete permute;
  a->GetInitial(d) = initial;
  return total;
}
  
Viterbi2DResult::~Viterbi2DResult() {
  if (last) delete last;
}
   

Viterbi2DResult *Viterbi2D(HMM2D::Ptr a, vector<HMM::Observation> &xobs, vector<HMM::Observation> &yobs, size_t s, size_t t, HMM2D::State &i, HMM2D::State &j, Viterbi2DResult *last) {
  Viterbi2DResult *result = new Viterbi2DResult();
  result->probability = Prob(a, HMM2D::Direction::X,  xobs[s], i, yobs.size())*Prob(a, HMM2D::Direction::Y, yobs[t], j, xobs.size());
  if (!result->probability) return result;
  double max = 0;
  Permutation *xset = EMMax(a, HMM2D::Direction::X, i.size() - 2);
  Permutation *yset = EMMax(a, HMM2D::Direction::Y, j.size() - 2);
  ForeachPermutation(xset, [&] (HMM2D::State xpartial, size_t xidx, double) -> bool {
    ForeachPermutation(yset, [&] (HMM2D::State ypartial, size_t yidx, double) -> bool {
      if (xpartial.back() != ypartial.back()) { return true; }
      double xprob, yprob;
      double overall;
      Viterbi2DResult *xviterbi = Viterbi2D(a, xobs, yobs, s - 1, t, xpartial, j, result);
      Viterbi2DResult *yviterbi = Viterbi2D(a, xobs, yobs, s, t - 1, i, ypartial, result);
      xprob = a->xtransition[a->states.size()*xpartial.back() + i.back()]*xviterbi->probability;
      yprob = a->ytransition[a->states.size()*ypartial.back() + j.back()]*yviterbi->probability;
      overall = xprob*yprob;
      if (overall > max) {
        max = overall;
        result->x = xpartial;
        result->y = ypartial;
        if (yprob > xprob) {
          result->direction = HMM2D::Direction::Y;
          delete xviterbi;
          result->last = yviterbi;
        } else {
          result->direction = HMM2D::Direction::X;
          delete yviterbi;
          result->last = xviterbi;
        }
      } else {
        delete xviterbi;
        delete yviterbi;
      }
      return true;
    });
    return true;
  });
  result->probability *= max;
  return result;
}

HMM2D::Ptr HMM2D::New() {
  return HMM2D::Ptr(new HMM2D());
}

vector<double> &HMM2D::GetTransition(HMM2D::Direction d) {
  if (d == HMM2D::Direction::X) {
    return xtransition;
  } else return ytransition;
}

vector<double> &HMM2D::GetInitial(HMM2D::Direction d) {
  if (d == HMM2D::Direction::X) {
    return xinitial;
  } else return yinitial;
}

Permutation *EMMax(HMM2D::Ptr a, HMM2D::Direction d, size_t len) {
  Permutation *result = new Permutation();
  result->probability = 100;
  result->state = 255;
  for (size_t i = 0; i < a->states.size(); ++i) {
    Permutation *previous = EM(a, d, len, a->states[i]);
    if (!previous->probability) {}
    else result->last.push_back(previous);
  }
  return result;
}

Permutation *EM(HMM2D::Ptr a, HMM2D::Direction d, size_t len, HMM2D::PartialState s) {
  Permutation *result;
  if (!len) {
    result = new Permutation();
    result->state = s;
    result->probability = a->GetInitial(d)[s];
    return result;
  }
  double max = 0;
  for (size_t i = 0; i < a->states.size(); i++) {
    result = new Permutation();
    Permutation *previous = EM(a, d, len - 1, a->states[i]);
    result->probability = a->GetTransition(d)[i*a->states.size() + s] * previous->probability;
    if (result->probability > max) {
      max = result->probability;
      result->state = s;
    }
    if (!result->probability) {
      delete previous;
    }
    else result->last.push_back(previous);
  }
  result->probability = max;
  sort(result->last.begin(), result->last.end(), [] (Permutation *a, Permutation *b) -> bool {
    if (a->probability > b->probability) return true;
    else return false;
  });
  return result;
}


ViterbiResult::ViterbiResult(ViterbiResult *plast, HMM::State *pptr, double pprobability) : last(plast), ptr(pptr), probability(pprobability) {}
ViterbiResult::~ViterbiResult() { if (last) delete last; }

ViterbiResult *Viterbi(HMM::Ptr m, const vector<HMM::Observation> &obs, size_t len, HMM::State *s, ViterbiResult *last) {
  size_t idx = s - &m->states[0];
  if (len == 0) {
    return new ViterbiResult(nullptr, s, m->emit[idx*m->obs.size() + obs[0]]*m->initial[idx]);
  }
  double max = 0;
  HMM::State *state = nullptr;
  ViterbiResult *previous = nullptr;
  for (auto it = m->states.begin(); it != m->states.end(); it++) {
    auto idx = distance(m->states.begin(), it);
    ViterbiResult *result = Viterbi(m, obs, len - 1, &*it, last);
    double prob = m->emit[idx*m->obs.size() + obs[len]]*m->matrix[distance(m->states.begin(), it)*m->states.size() + distance(m->states.begin(), find(m->states.begin(), m->states.end(), *result->ptr))] * result->probability;
    if (!previous) {
      max = prob;
      state = &*it;
      previous = result;
    } else if (prob > max) {
      delete previous;
      max = prob;
      state = &*it;
      previous = result;
    } else {
      delete result;
    }
  }
  return new ViterbiResult(previous, state, max);
}

ViterbiResult *ViterbiMax(HMM::Ptr m, vector<HMM::Observation> &obs) {
  for (auto it = obs.begin(); it != obs.end(); it++) {
    *it = distance(m->obs.begin(), find(m->obs.begin(), m->obs.end(), *it));
  }
  return ViterbiMax(m, obs, obs.size() - 1, nullptr);
}

ViterbiResult *ViterbiMax(HMM::Ptr m, const vector<HMM::Observation> &obs, size_t len, ViterbiResult *last) {
  double max = 0;
  HMM::State *s = nullptr;
  ViterbiResult *previous = nullptr;
  ViterbiResult *retval = nullptr;
  for (auto it = m->states.begin(); it != m->states.end(); it++) {
    ViterbiResult *result = Viterbi(m, obs, len - 1, &*it, last);
    double prob = m->emit[distance(m->states.begin(), it)*m->obs.size() + obs[len]]*m->matrix[distance(m->states.begin(), it)*m->states.size() + distance(m->states.begin(), find(m->states.begin(), m->states.end(), *result->ptr))] * result->probability;
    if (!previous) {
      max = prob;
      s = &*it;
      previous = result;
    } else if (prob > max) {
      delete previous;
      max = prob;
      s = &*it;
      previous = result;
    } else {
      delete result;
    }
  }
  return new ViterbiResult(previous, s, max);
}

json_object *StateToJsonObject(HMM::State s) {
  json_object *retval = json_object_new_object();
  json_object_object_add(retval, "density", json_object_new_int(s.first));
  json_object_object_add(retval, "duration", json_object_new_int(s.second));
  return retval;
}

void PrintViterbiResult(ViterbiResult *vr) {
  json_object *result = json_object_new_array();
  while (vr) {
    json_object_array_add(result, StateToJsonObject(*vr->ptr));
    vr = vr->last;
  }
  cout << json_object_to_json_string(result) << endl;
  json_object_put(result);
}

template class PNG<PNG_FORMAT_GA>;

template void Die<char const*, char const*>(char const*, char const*);

template HMMGroup::Ptr CalculateHMMGroup<PNG<1>::Pixel>(PNG<1>::Pixel*, unsigned long*);
template void Die<char const*>(char const*);
template void Die<char const*, double>(char const*, double);
template void Die<char const*, int, char*>(char const*, int, char*);
template void Die<char const*, int>(char const*, int);
template void ParseFilename<unsigned char>(char*, std::vector<char*, std::allocator<char*> >&, std::vector<unsigned char, std::allocator<unsigned char> >&);
template PNG<1>::Pixel* ZSlice<PNG<1>::Pixel>(PNG<1>::Pixel*, unsigned long, unsigned long, unsigned long);

Permutation::~Permutation() {
  for (auto it = last.begin(); it != last.end(); it++) {
    delete *it;
  }
}
/*

Permutation *EM(HMM::Ptr a, size_t len, HMM2D::PartialState state) {
  Permutation *result;
  if (!len) {
    result = new Permutation();
    result->state = state;
    result->probability = a->initial[state];
    return result;
  }
  double max = 0;
  for (size_t i = 0; i < a->states.size(); i++) {
    result = new Permutation();
    Permutation *previous = EM(a, len - 1, i);
    result->probability = a->matrix[i * a->states.size() + state] * previous->probability;
    if (result->probability > max) {
      max = result->probability;
    }
    if (!result->probability) {
      delete previous;
    }
    else result->last.push_back(previous);
  }
  result->probability = max;
  sort(result->last.begin(), result->last.end(), [] (Permutation *a, Permutation *b) -> bool {
    if (a->probability > b->probability) return true;
    else return false;
  });
  return result;
} 

*/
