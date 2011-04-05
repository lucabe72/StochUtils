#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

#include "distr.h"
#include "compute.h"

#define ITERATIONS 500

double *matrix_generate_pseudo(struct distribution *exec, int period, int qs, int ts)
{
  double *matrix;
  int q1;

  q1 = qs * (period / ts);
  matrix = pseudo_generate(exec, q1, 2000);

  return matrix;
}

double *matrix_generate_generic(struct distribution *exec, struct distribution *interarrival, int qs, int ts)
{
  double *matrix;
  double *v;

fprintf(stderr, "Transofrming\n");
  v = generic_transform(interarrival, ts, 2000);
  if (v == NULL) {
    fprintf(stderr, "Failed to transform the interarrival vector...\n");
    return NULL;
  }
fprintf(stderr, "Generating\n");
  matrix = generic_generate(exec, v, qs, 2000);

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
