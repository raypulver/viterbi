#ifndef HMM_H
#define HMM_H

#include <unistd.h>
#include "viterbi.h"
#include "cache.h"
#include "matrix.h"
#include "vector.h"

typedef struct _hmm2d_t {
  size_t n;
  size_t *states;
  size_t t;
  size_t *obs;
  matrix_t *ax;
  matrix_t *ay;
  vector_t *pix;
  vector_t *piy;
  matrix_t *bx;
  matrix_t *by;
  obs_vector_t *xobs;
  obs_vector_t *yobs;
} hmm2d_t;

#ifdef __cplusplus
extern "C" {
#endif

viterbi2d_result_t *init_viterbi2d_result();
void viterbi2d_free(viterbi2d_result_t *res);
long state_to_idx(hmm2d_t *hmm, size_t k);
viterbi2d_result_t *viterbi2d(hmm2d_t *hmm, cache_t *cache, size_t t, size_t k);
viterbi2d_result_t *viterbi2d_max(hmm2d_t *hmm, cache_t **cache);
viterbi2d_result_t *viterbi2d_max_no_cache(hmm2d_t *hmm);
hmm2d_t *init_hmm2d();
void print_hmm(hmm2d_t *);
void reconstruct(hmm2d_t *hmm, cache_t *cache, const char *filename);

#ifdef __cplusplus
}
#endif

#endif
