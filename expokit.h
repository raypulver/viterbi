#ifndef EXPOKIT_H
#define EXPOKIT_H

#include "matrix.h"

#ifdef __cplusplus

extern "C" void dgpadm_(int *ideg, int *m, double *t, double *H, int *ldh, double wsp[], int *lwsp, double ipiv[], int *iexph, int *ns, int *iflag);

#else

extern void dgpadm_(int *ideg, int *m, double *t, double *H, int *ldh, double wsp[], int *lwsp, double ipiv[], int *iexph, int *ns, int *iflag);

#endif

#ifdef __cplusplus
extern "C" {
#endif

matrix_t *matrix_exponentiate(matrix_t *m, double t, int *flags);

#ifdef __cplusplus
}
#endif
#endif
