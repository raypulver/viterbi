#include "voxelizer.h"
#include <iostream>
using namespace std;
int main(int argc, char **argv) {
  HMM::State a (0, 0);
  HMM::State b (0, 1);
  HMM::Observation y = 0;
  HMM::Observation z = 1;
  HMM::Ptr h = HMM::New();
  h->states.push_back(a);
  h->states.push_back(b);
  h->initial.push_back(0.90);
  h->initial.push_back(0.10);
  h->obs.push_back(y);
  h->obs.push_back(z);
  h->emit.push_back(1);
  h->emit.push_back(0);
  h->emit.push_back(0);
  h->emit.push_back(1);
  h->matrix.push_back(0.90);
  h->matrix.push_back(0.10);
  h->matrix.push_back(0.10);
  h->matrix.push_back(0.90);
  vector<HMM::Observation> obs;
  obs.push_back(z);
  obs.push_back(z);
  obs.push_back(z);
  obs.push_back(z);
  obs.push_back(z);
  ViterbiResult *vr = ViterbiMax(h, obs);
  PrintViterbiResult(vr);
  cout << vr->probability << endl;
  return 0;
}
