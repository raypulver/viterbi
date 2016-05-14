#include "expokit.h"
int main(int argc, char **argv) {
  int order, ideg, lwsp, ldh, iexph, ns, iflag;
  double t;
  double *H, *wsp, *ipiv;
  DGPADM(&ideg, &order, &t, H, &ldh, wsp, &lwsp, ipiv, &iexph, &ns, &iflag);
}
  
