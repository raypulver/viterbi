#ifndef VITERBI_H
#define VITERBI_H

#include <unistd.h>

typedef struct _viterbi2d_result_t {
  size_t x;
  size_t y;
  struct _viterbi2d_result_t *lastx;
  struct _viterbi2d_result_t *lasty;
  long double probability;
} viterbi2d_result_t;

#endif
