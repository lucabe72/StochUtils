#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

#include "pmf.h"

#include "distr.h"
#include "compute.h"
#include "driver.h"

#define ITERATIONS 500

double *matrix_generate_pseudo(struct pmf *exec, int period, int qs, int ts)
{
  double *matrix;
  int q1;

  q1 = qs * (period / ts);
  matrix = pseudo_generate(exec, q1);

  return matrix;
}

double *matrix_generate_generic(struct pmf *exec, struct pmf *interarrival, int qs, int ts)
{
  double *matrix;
  double *v;

fprintf(stderr, "Transofrming\n");
  v = generic_transform(interarrival, ts);
  if (v == NULL) {
    fprintf(stderr, "Failed to transform the interarrival vector...\n");
    return NULL;
  }
fprintf(stderr, "Generating\n");
  matrix = generic_generate(exec, v, qs);

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
