#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

#include "pmf.h"

#include "distr.h"
#include "compute.h"
#include "driver.h"

#define ITERATIONS 500

struct distribution *pmf2distr(struct pmf *p)
{
  struct distribution *d;

  d = malloc(sizeof(struct distribution));
  d->values = p->elems;

  return d;
}

double *matrix_generate_pseudo(struct pmf *exec, int period, int qs, int ts)
{
  double *matrix;
  int q1;
  struct distribution *e1;

  e1 = pmf2distr(exec);
  q1 = qs * (period / ts);
  matrix = pseudo_generate(e1, q1, PMF_SIZE);

  return matrix;
}

double *matrix_generate_generic(struct pmf *exec, struct pmf *interarrival, int qs, int ts)
{
  double *matrix;
  double *v;
  struct distribution *e1, *i1;

  e1 = pmf2distr(exec);
  i1 = pmf2distr(interarrival);
fprintf(stderr, "Transofrming\n");
  v = generic_transform(i1, ts, PMF_SIZE);
  if (v == NULL) {
    fprintf(stderr, "Failed to transform the interarrival vector...\n");
    return NULL;
  }
fprintf(stderr, "Generating\n");
  matrix = generic_generate(e1, v, qs, PMF_SIZE);

  if (matrix == NULL) {
    fprintf(stderr, "Matrix generation failed\n");
    free(v);

    return NULL;
  }
  free(v);

  return matrix;
}

double *dl_generate(double *matrix, int size, int qs, int ts)
{
  double *sol, *values1;

fprintf(stderr, "Solving\n");
  sol = solve(matrix, size, ITERATIONS);
  /*
  sol = create_compute_window(matrix, size, ITERATIONS);
  */
  if (sol == NULL) {
    fprintf(stderr, "Unable to solve the eigenvector problem!!!\n");
    free(matrix);

    return NULL;
  }

  values1 = malloc((ceil((double)size / qs) * ts + 1) * sizeof(double));
  if (values1 != NULL) {
    int dd;
    int i;

    memset(values1, 0, (ceil((double)size / qs) * ts + 1) * sizeof(double));
    for (i = 0; i < size; i++) {
      dd = ceil((double)i / qs) * ts;
      values1[dd] += sol[i];
    }
  }
  free(matrix);
  free(sol);

  return values1;
}
