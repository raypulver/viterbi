#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <json-c/json.h>
#include "hmm.h"
#include "cache.h"

matrix_t *init_matrix(size_t x, size_t y) {
	size_t i;
  matrix_t *retval = (matrix_t *) malloc(sizeof(matrix_t));
  if (!retval) return retval;
  retval->alloc = x*y;
  retval->x = x;
  retval->y = y;
  retval->data = (long double *) calloc(retval->alloc, sizeof(long double));
	memset(retval->data, 0, retval->alloc);
  return retval;
}

void free_matrix(matrix_t *matrix) {
	free(matrix->data);
	free(matrix);
}

vector_t *init_vector(size_t x) {
  size_t i, last_alloc;
  vector_t *retval = (vector_t *) malloc(sizeof(vector_t));
  if (!retval) return retval;
  last_alloc = retval->alloc;
  retval->alloc = x;
  retval->data = (long double *) calloc(retval->alloc, sizeof(long double));
	memset(retval->data, 0, retval->alloc);
  retval->len = 0;
  return retval;
}

obs_vector_t *init_obs_vector(size_t x) {
	size_t i, last_alloc;
  obs_vector_t *retval = (obs_vector_t *) malloc(sizeof(obs_vector_t));
  if (!retval) return retval;
	last_alloc = retval->alloc;
  retval->alloc = x;
  retval->data = (size_t *) calloc(retval->alloc, sizeof(size_t));
	memset(retval->data, 0, sizeof(size_t)*retval->alloc);
  retval->len = 0;
  return retval;
}
viterbi2d_result_t *init_viterbi2d_result() {
  viterbi2d_result_t *retval = (viterbi2d_result_t *) malloc(sizeof(viterbi2d_result_t));
  memset(retval, 0, sizeof(viterbi2d_result_t));
  return retval;
}

int obs_vector_push(obs_vector_t *vec, size_t obs) {
  if (!vec) return 2;
  if (vec->len + 1 > vec->alloc) {
    vec->data = (size_t *) realloc(vec->data, (vec->alloc << 1)*sizeof(size_t));
    if (!vec->data) return 1;
    vec->alloc <<= 1;
  }
  vec->data[vec->len] = obs;
  vec->len++;
  return 0;
}
int vector_push(vector_t *vec, long double d) {
  if (!vec) return 2;
  if (vec->len + 1 > vec->alloc) {
    vec->data = (long double *) realloc(vec->data, (vec->alloc << 1)*sizeof(long double));
    if (!vec->data) return 1;
    vec->alloc <<= 1;
  }
  vec->data[vec->len] = d;
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
void obs_vector_free(obs_vector_t *vec) {
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

double prob(hmm2d_t *hmm, cache_t *cache, size_t t, size_t k, int isx) {
	size_t i, j, pos, remaining;
	double weightedsum = 0, sum = 0;
	if (isx) {
		pos = t/hmm->yobs->len;
		remaining = hmm->yobs->len - t % hmm->yobs->len;
		for (i = pos*hmm->yobs->len; i < (pos + 1)*hmm->yobs->len; ++i) {
			for (j = 0; j < hmm->n; ++j) {
				if (*cache_el(cache, i, j)) {
  			  weightedsum += hmm->states[j]*((*cache_el(cache, i, j))->probability);
			  }
			}
		}
	} else {
		pos = t % hmm->yobs->len;
		remaining = hmm->xobs->len - t/hmm->yobs->len;
		for (i = pos; i < hmm->yobs->len*hmm->xobs->len; i += hmm->yobs->len) {
			for (j = 0; j < hmm->n; ++j) {
				if (*cache_el(cache, i, j)) {
  				weightedsum += hmm->states[j]*(*cache_el(cache, i, j))->probability;
				}
			}
		}
	}
	if (isx) {
  	if (k <= hmm->xobs->data[pos] - weightedsum) return 1;
		else return 0;
	} else {
		if (k <= hmm->yobs->data[pos] - weightedsum) return 1;
		else return 0;
	}
}

  
viterbi2d_result_t *viterbi2d(hmm2d_t *hmm, void *cache, size_t t, size_t k) {
  size_t x, y;
  long double max, overall;
  long idx;
	int cmp;
  viterbi2d_result_t *xviterbi, *yviterbi, *retval;
  idx = state_to_idx(hmm, k);
  if (idx == -1) return NULL;
  if ((retval = cache_get(cache, t, idx))) {
		return retval;
	}
	max = 0;
  retval = init_viterbi2d_result();
	retval->probability = prob(hmm, cache, t, k, 1)*prob(hmm, cache, t, k, 0);
  if (t == 0) {
    retval->x = k;
    retval->y = k;
    retval->probability = sqrt(retval->probability*(*vector_el(hmm->pix, idx))*(*vector_el(hmm->piy, idx)));
  } else if (t < hmm->xobs->len) {
    for (x = 0; x < hmm->n; ++x) {
      xviterbi = viterbi2d(hmm, cache, t - 1, hmm->states[x]);
      overall = sqrt(retval->probability*xviterbi->probability*(*matrix_el(hmm->ax, x, idx))*(*vector_el(hmm->piy, idx)));
      if (overall > max || (overall == max && !x)) {
        retval->lastx = xviterbi;
        retval->x = x;
        retval->probability = overall;
        max = overall;
      }
    }
  } else if (!(t % hmm->xobs->len)) {
    for (y = 0; y < hmm->n; ++y) {
      yviterbi = viterbi2d(hmm, cache, t - hmm->xobs->len, hmm->states[y]);
      overall = sqrt(retval->probability*yviterbi->probability*(*matrix_el(hmm->ay, y, idx))*(*vector_el(hmm->pix, idx)));
			if (overall > max || (overall == max && !y)) {
        retval->lasty = yviterbi;
        retval->y = y;
        retval->probability = overall;
        max = overall;
      }
    }
  } else {
    for (x = 0; x < hmm->n; ++x) {
      for (y = 0; y < hmm->n; ++y) {
        xviterbi = viterbi2d(hmm, cache, t - 1, hmm->states[x]);
        yviterbi = viterbi2d(hmm, cache, t - hmm->xobs->len, hmm->states[y]);
				overall = sqrt(retval->probability*xviterbi->probability*(*matrix_el(hmm->ax, x, idx))*yviterbi->probability*(*matrix_el(hmm->ay, y, idx)));
        if (overall > max || (overall == max && !x && !y)) {
          retval->lastx = xviterbi;
          retval->lasty = yviterbi;
          retval->x = x;
          retval->y = y;
          retval->probability = overall;
					max = overall;
        }
      }
    }
  }
  cache_put(cache, t, idx, retval);
  return retval;
}

viterbi2d_result_t *viterbi2d_max(hmm2d_t *hmm) {
  size_t x;
  long double max;
  size_t len;
  cache_t *cache;
  viterbi2d_result_t *result, *retval = NULL;
	max = 0;
  len = hmm->xobs->len*hmm->yobs->len - 1;
  cache = init_cache(len, hmm->n);
  for (x = 0; x < hmm->n; ++x) {
    result = viterbi2d(hmm, cache, len, hmm->states[x]);
    if (result->probability > max || (result->probability == max && !x)) {
      retval = result;
			max = result->probability;
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
void print_hmm(hmm2d_t *hmm) {
	size_t i, j;
	json_object *retval, *xrow, *yrow, *pix = json_object_new_array(), *axmatrix = json_object_new_array(), *piy = json_object_new_array(), *aymatrix = json_object_new_array();
  retval = json_object_new_object();
	for (i = 0; i < hmm->n; ++i) {
		xrow = json_object_new_array();
		yrow = json_object_new_array();
		for (j = 0; j < hmm->n; ++j) {
			json_object_array_add(xrow, json_object_new_double(*matrix_el(hmm->ax, i, j)));
			json_object_array_add(yrow, json_object_new_double(*matrix_el(hmm->ay, i, j)));
		}
		json_object_array_add(pix, json_object_new_double(*vector_el(hmm->pix, i)));
		json_object_array_add(piy, json_object_new_double(*vector_el(hmm->piy, i)));
		json_object_array_add(axmatrix, xrow);
		json_object_array_add(aymatrix, yrow);
	}
	json_object_object_add(retval, "pix", pix);
	json_object_object_add(retval, "piy", piy);
	json_object_object_add(retval, "xtransition", axmatrix);
	json_object_object_add(retval, "ytransition", aymatrix);
	printf("%s", json_object_to_json_string(retval));
	json_object_put(retval);
}
