#ifndef CACHE_H
#define CACHE_H

#include <unistd.h>
#include <stdint.h>
#include "viterbi.h"

typedef struct _cache_entry_t {
  viterbi2d_result_t *ptr;
  uint8_t dirty;
} cache_entry_t;

typedef struct _cache_t {
  cache_entry_t *data;
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
void cache_el_mark(cache_t *cache, size_t t);
uint8_t cache_el_is_marked(cache_t *cache, size_t t);
viterbi2d_result_t *cache_get(cache_t *cache, size_t t, size_t k);
void cache_put(cache_t *cache, size_t t, size_t k, viterbi2d_result_t *res);
void cache_empty(cache_t *cache);
int check_cache(cache_t *cache);

#ifdef __cplusplus
}
#endif

#endif
