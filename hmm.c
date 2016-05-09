#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <json-c/json.h>
#include <gmp.h>
#include "hmm.h"
#include "cache.h"

matrix_t *init_matrix(size_t x, size_t y) {
	size_t i;
  matrix_t *retval = (matrix_t *) malloc(sizeof(matrix_t));
  if (!retval) return retval;
  retval->alloc = x*y;
  retval->x = x;
  retval->y = y;
  retval->data = (mpq_ptr) calloc(retval->alloc, sizeof(mpq_t));
	for (i = 0; i < retval->alloc; ++i) {
		mpq_init(&retval->data[i]);
	}
  return retval;
}

void free_matrix(matrix_t *matrix) {
  size_t i;
	for (i = 0; i < matrix->alloc; ++i) {
		mpq_clear(&matrix->data[i]);
	}
	free(matrix->data);
	free(matrix);
}

vector_t *init_vector(size_t x) {
  size_t i, last_alloc;
  vector_t *retval = (vector_t *) malloc(sizeof(vector_t));
  if (!retval) return retval;
  last_alloc = retval->alloc;
  retval->alloc = x;
  retval->data = (mpq_ptr) calloc(retval->alloc, sizeof(mpq_t));
  for (i = 0; i < retval->alloc; ++i) {
    mpq_init(&retval->data[i]);
  }
  retval->len = 0;
  return retval;
}

void vector_free(vector_t *vec) {
	size_t i;
	for (i = 0; i < vec->alloc; ++i) {
		mpq_clear(&vec->data[i]);
	}
	free(vec->data);
	free(vec);
}

viterbi2d_result_t *init_viterbi2d_result() {
  viterbi2d_result_t *retval = (viterbi2d_result_t *) malloc(sizeof(viterbi2d_result_t));
  memset(retval, 0, sizeof(viterbi2d_result_t));
	mpq_init(retval->probability);
	mpq_set_ui(retval->probability, 0, 1);
  return retval;
}

static void woop() {}

int vector_push(vector_t *vec, long double d) {
  if (!vec) return 2;
  if (vec->len + 1 > vec->alloc) {
    vec->data = (mpq_ptr) realloc(vec->data, (vec->alloc << 1)*sizeof(mpq_t));
    if (!vec->data) return 1;
    vec->alloc <<= 1;
  }
  mpq_set_d(&vec->data[vec->len], d);
  vec->len++;
  return 0;
}

mpq_ptr vector_el(vector_t *vec, size_t i) {
  if (i < vec->len) {
    return &vec->data[i];
  } else return NULL;
}

mpq_ptr matrix_el(matrix_t *m, size_t x, size_t y) {
  return &m->data[x*m->y + y];
}

void viterbi2d_free(viterbi2d_result_t *res) {
  if (res->lastx) viterbi2d_free(res->lastx);
  if (res->lasty) viterbi2d_free(res->lasty);
	mpq_clear(res->probability);
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
  mpq_t max, overall;
  long idx;
	int cmp;
  viterbi2d_result_t *xviterbi, *yviterbi, *retval;
  idx = state_to_idx(hmm, k);
  if (idx == -1) return NULL;
  if ((retval = cache_get(cache, t, idx))) {
		return retval;
	}
	mpq_init(max);
	mpq_set_ui(max, 0, 1);
	mpq_init(overall);
  retval = init_viterbi2d_result();
  if (t == 0) {
    retval->x = k;
    retval->y = k;
    mpq_mul(retval->probability, vector_el(hmm->pix, idx), vector_el(hmm->piy, idx));
		mpq_canonicalize(retval->probability);
  } else if (t < hmm->xobs->len) {
    for (x = 0; x < hmm->n; ++x) {
      xviterbi = viterbi2d(hmm, cache, t - 1, hmm->states[x]);
      mpq_mul(overall, xviterbi->probability, matrix_el(hmm->ax, x, idx));
			mpq_mul(overall, overall, vector_el(hmm->piy, idx));
			mpq_canonicalize(overall);
      if ((cmp = mpq_cmp(overall, max)) > 0) {
        retval->lastx = xviterbi;
        retval->x = x;
        mpq_set(retval->probability, overall);
        mpq_swap(max, overall);
      }
    }
  } else if (!(t % hmm->xobs->len)) {
    for (y = 0; y < hmm->n; ++y) {
      yviterbi = viterbi2d(hmm, cache, t - hmm->xobs->len, hmm->states[y]);
      mpq_mul(overall, yviterbi->probability, matrix_el(hmm->ay, y, idx));
			mpq_mul(overall, overall, vector_el(hmm->pix, idx));
			mpq_canonicalize(overall);
      if ((cmp = mpq_cmp(overall, max)) > 0) {
        retval->lasty = yviterbi;
        retval->y = y;
        mpq_set(retval->probability, overall);
        mpq_swap(max, overall);
      }
    }
  } else {
    for (x = 0; x < hmm->n; ++x) {
      for (y = 0; y < hmm->n; ++y) {
        xviterbi = viterbi2d(hmm, cache, t - 1, hmm->states[x]);
        yviterbi = viterbi2d(hmm, cache, t - hmm->xobs->len, hmm->states[y]);
        mpq_mul(overall, xviterbi->probability, matrix_el(hmm->ax, x, idx));
				mpq_canonicalize(overall);
			  mpq_mul(overall, overall, yviterbi->probability);
				mpq_canonicalize(overall);
				mpq_mul(overall, overall, matrix_el(hmm->ay, y, idx));
				mpq_canonicalize(overall);
        if ((cmp = mpq_cmp(overall, max)) > 0) {
          retval->lastx = xviterbi;
          retval->lasty = yviterbi;
          retval->x = x;
          retval->y = y;
          mpq_set(retval->probability, overall);
					mpq_swap(max, overall);
        }
      }
    }
  }
  cache_put(cache, t, idx, retval);
  return retval;
}

viterbi2d_result_t *viterbi2d_max(hmm2d_t *hmm) {
  size_t x;
  mpq_t max;
  size_t len;
  cache_t *cache;
	int cmp;
  viterbi2d_result_t *result, *retval = NULL;
  len = hmm->xobs->len*hmm->yobs->len - 1;
  cache = init_cache(len, hmm->n);
	mpq_init(max);
	mpq_set_ui(max, 0, 1);
  for (x = 0; x < hmm->n; ++x) {
    result = viterbi2d(hmm, cache, len, hmm->states[x]);
    if ((cmp = mpq_cmp(result->probability, max)) > 0) {
      retval = result;
      mpq_set(max, result->probability);
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
			json_object_array_add(xrow, json_object_new_double(mpq_get_d(matrix_el(hmm->ax, i, j))));
			json_object_array_add(yrow, json_object_new_double(mpq_get_d(matrix_el(hmm->ay, i, j))));
		}
		json_object_array_add(pix, json_object_new_double(mpq_get_d(vector_el(hmm->pix, i))));
		json_object_array_add(piy, json_object_new_double(mpq_get_d(vector_el(hmm->piy, i))));
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
