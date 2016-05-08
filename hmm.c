#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "hmm.h"
#include "cache.h"

matrix_t *init_matrix(size_t x, size_t y) {
  matrix_t *retval = (matrix_t *) malloc(sizeof(matrix_t));
  if (!retval) return retval;
  retval->alloc = x*y;
  retval->x = x;
  retval->y = y;
  retval->data = (long double *) calloc(x*y, sizeof(long double));
  return retval;
}

vector_t *init_vector(size_t x) {
  vector_t *retval = (vector_t *) malloc(sizeof(vector_t));
  if (!retval) return retval;
  retval->alloc = x;
  retval->data = (long double *) calloc(x, sizeof(long double));
  retval->len = 0;
  return retval;
}

viterbi2d_result_t *init_viterbi2d_result() {
  viterbi2d_result_t *retval = (viterbi2d_result_t *) malloc(sizeof(viterbi2d_result_t));
  memset(retval, 0, sizeof(viterbi2d_result_t));
  return retval;
}

int vector_push(vector_t *vec, long double el) {
  if (!vec) return 2;
  if (vec->len + 1 > vec->alloc) {
    vec->data = (long double *) realloc(vec->data, (vec->alloc << 1)*sizeof(long double));
    if (!vec->data) return 1;
    vec->alloc <<= 1;
  }
  vec->data[vec->len] = el;
  vec->len++;
  return 0;
}

long double *vector_el(vector_t *vec, size_t i) {
  if (i < vec->len) {
    return &vec->data[i];
  } else return NULL;
}

void vector_free(vector_t *vec) {
  free(vec->data);
  free(vec);
}

long double *matrix_el(matrix_t *m, size_t x, size_t y) {
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
  long double max, xprob, yprob, overall;
  long idx;
  viterbi2d_result_t *xviterbi, *yviterbi, *retval;
  idx = state_to_idx(hmm, k);
  if (idx == -1) return NULL;
  if ((retval = cache_get(cache, t, idx))) return retval;
  retval = init_viterbi2d_result();
  if (t == 0) {
    retval->x = k;
    retval->y = k;
    retval->probability = (*vector_el(hmm->pix, idx))*(*vector_el(hmm->piy, idx));
  } else if (t < hmm->xobs->len) {
    max = 0;
    for (x = 0; x < hmm->n; ++x) {
      xviterbi = viterbi2d(hmm, cache, t - 1, hmm->states[x]);
      overall = xviterbi->probability*(*matrix_el(hmm->ax, x, idx))*(*vector_el(hmm->piy, idx));
      if (overall > max) {
        retval->lastx = xviterbi;
        retval->x = x;
        retval->probability = overall;
        max = overall;
      }
    }
  } else if (!(t % hmm->xobs->len)) {
    max = 0;
    for (y = 0; y < hmm->n; ++y) {
      yviterbi = viterbi2d(hmm, cache, t - hmm->xobs->len, hmm->states[y]);
      overall = yviterbi->probability*(*matrix_el(hmm->ay, y, idx))*(*vector_el(hmm->pix, idx));
      if (overall > max) {
        retval->lasty = yviterbi;
        retval->y = y;
        retval->probability = overall;
        max = overall;
      }
    }
  } else {
    max = 0;
    for (x = 0; x < hmm->n; ++x) {
      for (y = 0; y < hmm->n; ++y) {
        xviterbi = viterbi2d(hmm, cache, t - 1, hmm->states[x]);
        yviterbi = viterbi2d(hmm, cache, t - hmm->xobs->len, hmm->states[y]);
        overall = xviterbi->probability*(*matrix_el(hmm->ax, x, idx)) * yviterbi->probability* (*matrix_el(hmm->ay, y, idx));
        if (overall > max) {
          retval->lastx = xviterbi;
          retval->lasty = yviterbi;
          retval->x = x;
          retval->y = y;
          retval->probability = overall;
        }
      }
    }
  }
  cache_put(cache, t, idx, retval);
  return retval;
}

viterbi2d_result_t *viterbi2d_no_cache(hmm2d_t *hmm, size_t t, size_t k) {
  size_t x, y;
  long double max, xprob, yprob, overall;
  long idx;
  viterbi2d_result_t *xviterbi, *yviterbi, *retval;
  idx = state_to_idx(hmm, k);
  if (idx == -1) return NULL;
  retval = init_viterbi2d_result();
  if (t == 0) {
    retval->x = k;
    retval->y = k;
    retval->probability = (*vector_el(hmm->pix, idx))*(*vector_el(hmm->piy, idx));
  } else if (t < hmm->xobs->len) {
    max = 0;
    for (x = 0; x < hmm->n; ++x) {
      xviterbi = viterbi2d_no_cache(hmm, t - 1, hmm->states[x]);
      overall = xviterbi->probability*(*matrix_el(hmm->ax, x, idx))*(*vector_el(hmm->piy, idx));
      if (overall > max) {
        if (retval->lastx) free(retval->lastx);
        retval->lastx = xviterbi;
        retval->x = x;
        retval->probability = overall;
        max = overall;
      } else viterbi2d_free(xviterbi);
    }
  } else if (!(t % hmm->xobs->len)) {
    max = 0;
    for (y = 0; y < hmm->n; ++y) {
      yviterbi = viterbi2d_no_cache(hmm, t - hmm->xobs->len, hmm->states[y]);
      overall = yviterbi->probability*(*matrix_el(hmm->ay, y, idx))*(*vector_el(hmm->pix, idx));
      if (overall > max) {
        if (retval->lasty) viterbi2d_free(retval->lasty);
        retval->lasty = yviterbi;
        retval->y = y;
        retval->probability = overall;
        max = overall;
      } else viterbi2d_free(yviterbi);
    }
  } else {
    max = 0;
    for (x = 0; x < hmm->n; ++x) {
      for (y = 0; y < hmm->n; ++y) {
        xviterbi = viterbi2d_no_cache(hmm, t - 1, hmm->states[x]);
        yviterbi = viterbi2d_no_cache(hmm, t - hmm->xobs->len, hmm->states[y]);
        overall = xviterbi->probability*(*matrix_el(hmm->ax, x, idx)) * yviterbi->probability* (*matrix_el(hmm->ay, y, idx));
        if (overall > max) {
          if (retval->lastx) viterbi2d_free(retval->lastx);
          if (retval->lasty) viterbi2d_free(retval->lasty);
          retval->lastx = xviterbi;
          retval->lasty = yviterbi;
          retval->x = x;
          retval->y = y;
          retval->probability = overall;
        }
      }
    }
  }
  return retval;
}

viterbi2d_result_t *viterbi2d_max_no_cache(hmm2d_t *hmm) {
  size_t x;
  long double max;
  size_t len;
  viterbi2d_result_t *result, *retval = NULL;
  len = hmm->xobs->len*hmm->yobs->len;
  max = 0;
  for (x = 0; x < hmm->n; ++x) {
    result = viterbi2d_no_cache(hmm, len, hmm->states[x]);
    if (result->probability > max) {
      if (retval) viterbi2d_free(retval);
      retval = result;
      max = result->probability;
    } else viterbi2d_free(result);
  }
  return retval;
} 

static void woop () {}

viterbi2d_result_t *viterbi2d_max(hmm2d_t *hmm) {
  size_t x;
  long double max;
  size_t len;
  cache_t *cache;
  viterbi2d_result_t *result, *retval = NULL;
  len = hmm->xobs->len*hmm->yobs->len;
  cache = init_cache(len, hmm->n);
  max = 0;
  for (x = 0; x < hmm->n; ++x) {
    result = viterbi2d(hmm, cache, len, hmm->states[x]);
    if (result->probability > max) {
      retval = result;
      max = result->probability;
    }
  }
  woop();
  check_cache(cache);
  cache_free(cache);
  return retval;
} 
hmm2d_t *init_hmm2d() {
  hmm2d_t *retval = (hmm2d_t *) malloc(sizeof(hmm2d_t));
  memset(retval, 0, sizeof(hmm2d_t));
  return retval;
}
