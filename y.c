#include <stdio.h>
#include <stdlib.h>

#include "pmf.h"
#include "y.h"

static double epsilon(const struct pmf *c, int cmax)
{
  int i;
  double sum = 0;

  for (i = 0; i <= cmax; i++) {
    sum += pmf_get(c, i);
  }

  return 1 - sum;
}

int y_max(unsigned int x_max, unsigned int z_max, unsigned int q)
{
  return x_max - z_max * q;
}

struct pmf *compute(const struct pmf *c, const struct pmf *z, int q, int n)
{
  struct pmf *y;
  int zi, yi;
  double sum = 0;
  int ymax = y_max(pmf_max(c), pmf_min(z), q);

  y = pmf_create(n, n / 2);

  for (yi = -n / 2; yi <= ymax; yi++) {
    double py;
    int zmin, zmax;

    py = 0;
    zmin = -yi / q - 1;
    zmax = (pmf_max(c) - yi) / q + 1;
    if (zmax > pmf_max(z)) zmax = pmf_max(z);
    if (zmin < 1) zmin = 1;
    for (zi = zmin; zi <= zmax; zi++) {
//      if ((yi >= 1 - q * pmf_max(z)) && (yi < pmf_max(c) -  q)) {
      if ((yi + zi * q >= 0) && (yi + zi * q <= pmf_max(c))) {
        py += pmf_get(c, yi + zi * q) * pmf_get(z, zi);
      }
#if 0
      if (yi + zi * q > C) {
        fprintf(stderr, "c = %d + %d * %d > %d\n", yi, zi, q, C);
      } else {
        fprintf(stderr, "OK\n");
      }
#endif
//      }
    }
    if (pmf_set(y, yi, py) < 0) {
      fprintf(stderr, "Error inserting %d, %f\n", yi, py);
      exit(-1);
    }
    sum += py;
  }
printf("Max{c} = %d\n",  pmf_max(c));
  printf("Sum: %f; Ymax = %d; espilon = %f < %f\n", sum, pmf_max(c) - pmf_max(z) * q, 1 - sum, 
  pmf_max(c) * epsilon(c, pmf_max(c)));
  if (1.0 - sum > 1e-10) {
    y->tail = 1.0 - sum;
  }
  if (pmf_check(y)) {
    printf("Bad PMF: %d | %f!!!\n", pmf_check(y), pmf_sum(y, 0));
  }

  printf("min{Y}: %d\tmax{Y}: %d\t Ymax: %d\n", pmf_min(y), pmf_max(y), ymax);

  return y;
}
