#ifndef MATRIX_H
#define MATRIX_H

#include <unistd.h>

typedef struct _matrix_t {
  long double *data;
  size_t alloc;
  size_t x;
  size_t y;
} matrix_t;

#ifdef __cplusplus
extern "C" {
#endif

matrix_t *init_matrix(size_t x, size_t y);
long double *matrix_el(matrix_t *m, size_t x, size_t y);
void matrix_free(matrix_t *m);

#ifdef __cplusplus
}
#endif

#endif
