#include <stdlib.h>
#include <stdio.h>

#include "pmf.h"

#include "distr.h"

double *pseudo_generate(struct pmf *d, int q)
{
  int j, i;
  double *mat;
  int n = d->size;

  mat = malloc(sizeof(double) * n * n);
  if (mat == 0) {
    fprintf(stderr, "Error Allocating the Matrix...\n");
    return 0;
  }

  for(i = 0; i < n; i++) {
    for (j = 0; j < n; j++) {
      if (j <= q) {
/*        m[i][j] = d->values[i]; */
        *(mat + i * n + j) = d->elems[i];
      } else {
        if (i + q - j > 0) {
/*          m[i][j] = d->values[i + q - j]; */
          *(mat + i * n + j) = d->elems[i + q - j];
        } else {
/*          m[i][j] = 0; */
          *(mat + i * n + j) = 0;
        }
      }
    }
  }
  fprintf(stderr, "Matrix generated...\n");

  return mat;
}

double *solve(double *mat, int n, int iter)
{
  int k;
  int i, j;
  double *w, old, s, *sum;

  w = malloc(sizeof(double) * n);
  if (w == 0) {
    fprintf(stderr, "Error Allocating w...\n");
    return 0;
  }
  sum = malloc(sizeof(double) * n);
  if (sum == 0) {
    free(w);
    fprintf(stderr, "Error Allocating sum...\n");
    return 0;
  }

  for (j = 0; j < n; j++) {
    w[j] = 1;
    sum[j] = 0;
  }

  k = 0;
  old = 1000;
  s = 0;
  while(k < iter) {
    k++;
    old = w[0];
    s = 0;
    for (i = 0; i < n; i++) {
      sum[i] = 0;
      for (j = 0; j < n; j++) {
/*
        sum[i] += m[i][j] * w[j];
*/
        sum[i] += *(mat + i * n +j) * w[j];
      }
      s += sum[i];
    }
  
    for(i = 0; i < n; i++) {
      w[i] = sum[i] / s;
    }
  }
  free(sum);

  return w;
}
