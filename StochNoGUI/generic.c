#include <stdlib.h>
#include <stdio.h>

#include "pmf.h"

#include "distr.h"
#include "compute.h"

#define STEP 10
double *generic_transform(struct pmf *d_t, int period)
{
  double *v, *v1;
  int j, j1, th;
  int n = d_t->size;

  v = d_t->elems;
  v1 = malloc(sizeof(double) * n);
  if (v == NULL) {
    fprintf(stderr, "Error Allocating the temporary vector...\n");
    
    return NULL;
  }

  for (j = 0; j < n; j++) {
    v1[j] = 0;
  }
  j1 = 0;
  th = period;
/*  for (j = 0; j <= n; j++) {*/
  for (j = 0; j < n; j++) {
    if ((j < period) && (v[j] != 0)) {
      fprintf(stderr, "Error!!! Interarrival times MUST be greater than %d!!!\n",
		      period);
      return 0;
    }

    if (j >= th) {
      th += period;
      j1++;
    }
    v1[j1] += v[j];
  }

  return v1;
}

double *generic_generate(struct pmf *d_c, double *v, int q)
{
  int j, i, xx;
  double *mat;
  double *u = d_c->elems;
  int n = d_c->size;

  mat = malloc(sizeof(double) * n * n);
  if (mat == NULL) {
    fprintf(stderr, "Error Allocating the Matrix...\n");

    return NULL;
  }

  for (i = 0; i < n; i++) {
if (i % 100 == 0) printf("I: %d/%d\n", i, n);
    for (j = 0; j < n; j++) {
/*      mat[i][j] = 0; */
	*(mat + i * n + j) = 0;
      for (xx = 0; xx < n; xx++) {
        if (j <= xx * q) {
/*          mat[i][j] += u[i] * v[xx]; */
          *(mat + i * n + j) += u[i] * v[xx]; 
        } else {
          if ((i + xx * q - j >= 0) && (i + xx * q - j < n)) {
/*            mat[i][j] += u[i + xx * q - j] * v[xx];*/
            *(mat + i * n + j) += u[i + xx * q - j] * v[xx];
          } else {
/*            mat[i][j] = 0; */
            *(mat + i * n + j) = 0;
          }
        }
      }
    }
  }

  return mat;
}
