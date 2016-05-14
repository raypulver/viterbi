#include <stdio.h>
#include "expokit.h"
#include "matrix.h"

int main(int argc, char **argv) {
  size_t i, j;
  int flags;
  matrix_t *ret, *mymat = init_matrix(3, 3);
  *matrix_el(mymat, 0, 0) = 1;
  *matrix_el(mymat, 1, 1) = 2;
  *matrix_el(mymat, 2, 2) = 1;
  ret = matrix_exponentiate(mymat, 2, &flags);
  for (i = 0; i < mymat->x; ++i) {
    printf("[");
    for (j = 0; j < mymat->y; ++j) {
      printf("% lf ", (double) *matrix_el(ret, i, j));
    }
    printf("]\n");
  }
  printf("%i\n", flags);
}
