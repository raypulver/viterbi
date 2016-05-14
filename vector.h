#ifndef VECTOR_H
#define VECTOR_H

#include <unistd.h>

typedef struct _vector_t {
  long double *data;
  size_t alloc;
  size_t len;
} vector_t;

typedef struct _obs_vector_t {
	size_t *data;
	size_t alloc;
	size_t len;
} obs_vector_t;

#ifdef __cplusplus
extern "C" {
#endif

vector_t *init_vector(size_t x);
obs_vector_t *init_obs_vector(size_t x);
void vector_free(vector_t *vec);
void obs_vector_free(obs_vector_t *vec);
int vector_push(vector_t *vec, long double el);
int obs_vector_push(obs_vector_t *vec, size_t el);
long double *vector_el(vector_t *vec, size_t i);
size_t *obs_vector_el(obs_vector_t *vec, size_t i);

#ifdef __cplusplus
}
#endif

#endif
