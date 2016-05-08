#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include "cache.h"

cache_t *init_cache(size_t sz, size_t k) {
  cache_t *retval = (cache_t *) malloc(sizeof(cache_t));
  retval->t = sz;
  retval->k = k;
  retval->alloc = sz*k;
  retval->data = (viterbi2d_result_t **) calloc(sz*k, sizeof(viterbi2d_result_t *));
  memset(retval->data, 0, sz*k*sizeof(viterbi2d_result_t *));
  return retval;
}

viterbi2d_result_t **cache_el(cache_t *cache, size_t t, size_t k) {
  if (t > cache->t || k > cache->k) return NULL;
  return &cache->data[t*cache->k + k];
}

void cache_free(cache_t *cache) {
  size_t t, k;
  for (t = 0; t < cache->t; ++t) {
    for (k = 0; k < cache->k; ++k) {
      free(*cache_el(cache, t, k));
    }
  }
  free(cache->data);
  free(cache);
}

viterbi2d_result_t *cache_get(cache_t *cache, size_t t, size_t k) {
  return *cache_el(cache, t, k);
}

void cache_put(cache_t *cache, size_t t, size_t k, viterbi2d_result_t *res) {
  *cache_el(cache, t, k) = res;
}

void cache_empty(cache_t *cache) {
  memset(cache->data, 0, cache->k*cache->t*sizeof(viterbi2d_result_t *));
}

int check_cache(cache_t *cache) {
  size_t x, y;
  for (x = 0; x < cache->t; ++x) {
    for (y = 0; y < cache->k; ++y) {
      if ((*cache_el(cache, x, y))->probability) printf("%zu %zu %Lf\n", x, y, (*cache_el(cache, x, y))->probability);
    }
  }
  return 1;
}
