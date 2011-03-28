//#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "pmf.h"
#include "gamma.h"

struct pmf *gamma_generate(double gamma, unsigned int max, unsigned int step)
{
  unsigned int i;
  struct pmf *g;

  g = pmf_create(max, 0);
  for (i = 0; i < max; i += step) {
    pmf_set(g, i, 1.0 - pow(gamma, -(double)i));
    //printf("P{w > %d} <= %f\n", i, pow(gamma, -i));
  }

  return g;
}

#if 0
static double mypow(double a, int b)
{
  double res = 1.0;
  unsigned int i;

  if (b < 0) {
    b = -b;
    a = 1.0 / a;
  }
  for (i = 0; i < (unsigned int)b; i++) {
    res *= a;
  }

  return res;
}
#endif

double compute_sum(const struct pmf *y, double gamma)
{
  int i;
  double sum = 0.0;

  for (i = pmf_min(y); i <= pmf_max(y); i++) {
    sum += pow(gamma, i) * pmf_get(y, i);
#if 0
    if (pmf_get(y, i) > 0) {
      printf("%d: Sum = %f\n", sum);
    }
#endif
  }

  return sum;
}

double get_gamma(const struct pmf *y, unsigned int max)
{
  double gamma, mingamma, maxgamma;
  int cnt = 0;
  
  gamma = 2.0; mingamma = 1.0; maxgamma = 2.0; cnt = 0;
  while (1) {
    double eps;

    //printf("Trying gamma: %f %f\n", gamma, pmf_tail(y));
    if (((int)max > pmf_max(y)) && (pmf_tail(y) > 1e-10)) {
      eps = pow(gamma, max) * pmf_tail(y);
    } else {
      eps = 0;
    }
//printf("Epsilon: %f Sum: %f\n", eps, compute_sum(y, gamma));
    if (compute_sum(y, gamma) + eps <= 1.0) {
      mingamma = gamma;
      gamma = (gamma + maxgamma) / 2;
      cnt++;
    } else {
      maxgamma = gamma;
      gamma = (gamma + mingamma) / 2;
    }
    if (cnt == 5) {
      return mingamma;
    }
  }
}
