#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "hmm.h"
#include "cache.h"

matrix_t *init_matrix(size_t x, size_t y) {
	size_t i;
  matrix_t *retval = (matrix_t *) malloc(sizeof(matrix_t));
  if (!retval) return retval;
  retval->alloc = x*y;
  retval->x = x;
  retval->y = y;
  retval->data = (mpfr_t *) calloc(retval->alloc, sizeof(mpfr_t));
	for (i = 0; i < retval->alloc; ++i) {
		mpfr_init2(retval->data[i], MPFR_PRECISION);
		mpfr_set_ui(retval->data[i], 0, MPFR_RND);
	}
  return retval;
}

vector_t *init_vector(size_t x) {
	size_t i, last_alloc;
  vector_t *retval = (vector_t *) malloc(sizeof(vector_t));
  if (!retval) return retval;
	last_alloc = retval->alloc;
  retval->alloc = x;
  retval->data = (mpfr_t *) calloc(retval->alloc, sizeof(mpfr_t));
	for (i = 0; i < retval->alloc; ++i) {
		mpfr_init2(retval->data[i], MPFR_PRECISION);
		mpfr_set_ui(retval->data[i], 0, MPFR_RND);
	}
  retval->len = 0;
  return retval;
}

viterbi2d_result_t *init_viterbi2d_result() {
  viterbi2d_result_t *retval = (viterbi2d_result_t *) malloc(sizeof(viterbi2d_result_t));
  memset(retval, 0, sizeof(viterbi2d_result_t));
	mpfr_init2(retval->probability, MPFR_PRECISION);
	mpfr_set_ui(retval->probability, 0, MPFR_RND);
  return retval;
}

int vector_push(vector_t *vec, long double d) {
  if (!vec) return 2;
  if (vec->len + 1 > vec->alloc) {
    vec->data = (mpfr_t *) realloc(vec->data, (vec->alloc << 1)*sizeof(mpfr_t));
    if (!vec->data) return 1;
    vec->alloc <<= 1;
  }
  mpfr_set_d(vec->data[vec->len], d, MPFR_RND);
  vec->len++;
  return 0;
}

mpfr_t *vector_el(vector_t *vec, size_t i) {
  if (i < vec->len) {
    return &vec->data[i];
  } else return NULL;
}

void vector_free(vector_t *vec) {
  free(vec->data);
  free(vec);
}

mpfr_t *matrix_el(matrix_t *m, size_t x, size_t y) {
  return &m->data[x*m->y + y];
}

void viterbi2d_free(viterbi2d_result_t *res) {
  if (res->lastx) viterbi2d_free(res->lastx);
  if (res->lasty) viterbi2d_free(res->lasty);
  free(res);
}

long state_to_idx(hmm2d_t *hmm, size_t k) {
  int i;
  for (i = 0; i < hmm->n; ++i) {
    if (hmm->states[i] == k) return i;
  }
  return -1;
}
viterbi2d_result_t *viterbi2d(hmm2d_t *hmm, void *cache, size_t t, size_t k) {
  size_t x, y;
  mpfr_t max, overall;
  long idx;
	int cmp;
  viterbi2d_result_t *xviterbi, *yviterbi, *retval;
  idx = state_to_idx(hmm, k);
  if (idx == -1) return NULL;
  if ((retval = cache_get(cache, t, idx))) {
		return retval;
	}
	mpfr_init2(max, MPFR_PRECISION);
	mpfr_set_ui(max, 0, MPFR_RND);
	mpfr_init2(overall, MPFR_PRECISION);
  retval = init_viterbi2d_result();
  if (t == 0) {
    retval->x = k;
    retval->y = k;
    mpfr_mul(retval->probability, *vector_el(hmm->pix, idx), *vector_el(hmm->piy, idx), MPFR_RND);
  } else if (t < hmm->xobs->len) {
    for (x = 0; x < hmm->n; ++x) {
      xviterbi = viterbi2d(hmm, cache, t - 1, hmm->states[x]);
      mpfr_mul(overall, xviterbi->probability, *matrix_el(hmm->ax, x, idx), MPFR_RND);
			mpfr_mul(overall, overall, *vector_el(hmm->piy, idx), MPFR_RND);
      if ((cmp = mpfr_cmp(overall, max)) > 0 || !(cmp || x)) {
        retval->lastx = xviterbi;
        retval->x = x;
        mpfr_set(retval->probability, overall, MPFR_RND);
        mpfr_swap(max, overall);
      }
    }
  } else if (!(t % hmm->xobs->len)) {
    for (y = 0; y < hmm->n; ++y) {
      yviterbi = viterbi2d(hmm, cache, t - hmm->xobs->len, hmm->states[y]);
      mpfr_mul(overall, yviterbi->probability, *matrix_el(hmm->ay, y, idx), MPFR_RND);
			mpfr_mul(overall, overall, *vector_el(hmm->pix, idx), MPFR_RND);
      if ((cmp = mpfr_cmp(overall, max)) > 0 || !(cmp || y)) {
        retval->lasty = yviterbi;
        retval->y = y;
        mpfr_set(retval->probability, overall, MPFR_RND);
        mpfr_swap(max, overall);
      }
    }
  } else {
    for (x = 0; x < hmm->n; ++x) {
      for (y = 0; y < hmm->n; ++y) {
        xviterbi = viterbi2d(hmm, cache, t - 1, hmm->states[x]);
        yviterbi = viterbi2d(hmm, cache, t - hmm->xobs->len, hmm->states[y]);
        mpfr_mul(overall, xviterbi->probability, *matrix_el(hmm->ax, x, idx), MPFR_RND);
			  mpfr_mul(overall, overall, yviterbi->probability, MPFR_RND);
				mpfr_mul(overall, overall, *matrix_el(hmm->ay, y, idx), MPFR_RND);
        if ((cmp = mpfr_cmp(overall, max)) > 0 || !(cmp || x || y)) {
          retval->lastx = xviterbi;
          retval->lasty = yviterbi;
          retval->x = x;
          retval->y = y;
          mpfr_set(retval->probability, overall, MPFR_RND);
					mpfr_swap(max, overall);
        }
      }
    }
  }
  cache_put(cache, t, idx, retval);
  return retval;
}

viterbi2d_result_t *viterbi2d_max(hmm2d_t *hmm) {
  size_t x;
  mpfr_t max;
  size_t len;
  cache_t *cache;
	int cmp;
  viterbi2d_result_t *result, *retval = NULL;
  len = hmm->xobs->len*hmm->yobs->len - 1;
  cache = init_cache(len, hmm->n);
	mpfr_init2(max, MPFR_PRECISION);
	mpfr_set_ui(max, 0, MPFR_RND);
  for (x = 0; x < hmm->n; ++x) {
    result = viterbi2d(hmm, cache, len, hmm->states[x]);
    if ((cmp = mpfr_cmp(result->probability, max)) > 0 || !(cmp || x)) {
      retval = result;
      mpfr_set(max, result->probability, MPFR_RND);
    }
  }
  check_cache(cache);
  cache_free(cache);
  return retval;
} 
hmm2d_t *init_hmm2d() {
  hmm2d_t *retval = (hmm2d_t *) malloc(sizeof(hmm2d_t));
  memset(retval, 0, sizeof(hmm2d_t));
  return retval;
}
