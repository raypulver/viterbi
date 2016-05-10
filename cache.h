#ifndef CACHE_H
#define CACHE_H

#include "hmm.h"

typedef struct _cache_t {
  viterbi2d_result_t **data;
  size_t alloc;
  size_t t;
  size_t k;
} cache_t;

#ifdef __cplusplus
extern "C" {
#endif


cache_t *init_cache(size_t sz, size_t k);
void cache_free(cache_t *cache);
viterbi2d_result_t **cache_el(cache_t *cache, size_t t, size_t k);
viterbi2d_result_t *cache_get(cache_t *cache, size_t t, size_t k);
void cache_put(cache_t *cache, size_t t, size_t k, viterbi2d_result_t *res);
void cache_empty(cache_t *cache);
int check_cache(cache_t *cache);

#ifdef __cplusplus
}
#endif

#endif
