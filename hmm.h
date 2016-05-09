#ifndef HMM_H
#define HMM_H
#define MPFR_PRECISION (1 << 7)
#define MPFR_RND MPFR_RNDN

#include <unistd.h>
#include <mpfr.h>

typedef struct _matrix_t {
  mpfr_t *data;
  size_t alloc;
  size_t x;
  size_t y;
} matrix_t;

typedef struct _vector_t {
  mpfr_t *data;
  size_t alloc;
  size_t len;
} vector_t;

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
  vector_t *xobs;
  vector_t *yobs;
} hmm2d_t;

typedef struct _viterbi2d_result_t {
  size_t x;
  size_t y;
  struct _viterbi2d_result_t *lastx;
  struct _viterbi2d_result_t *lasty;
  mpfr_t probability;
} viterbi2d_result_t;

#ifdef __cplusplus
extern "C" {
#endif

matrix_t *init_matrix(size_t x, size_t y);
mpfr_t *matrix_el(matrix_t *m, size_t x, size_t y);
vector_t *init_vector(size_t x);
void vector_free(vector_t *vec);
int vector_push(vector_t *vec, long double el);
mpfr_t *vector_el(vector_t *vec, size_t i);

viterbi2d_result_t *init_viterbi2d_result();
void viterbi2d_free(viterbi2d_result_t *res);
long state_to_idx(hmm2d_t *hmm, size_t k);
viterbi2d_result_t *viterbi2d(hmm2d_t *hmm, void *cache, size_t t, size_t k);
viterbi2d_result_t *viterbi2d_max(hmm2d_t *hmm);
viterbi2d_result_t *viterbi2d_max_no_cache(hmm2d_t *hmm);
hmm2d_t *init_hmm2d();

#ifdef __cplusplus
}
#endif

#endif
