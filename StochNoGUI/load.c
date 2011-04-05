#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "distr.h"
//#include "pdf.h"
#include "load.h"

struct distribution *pdf_load(char *name, int size)
{
  double *values;
  int n, done;
  FILE *f;
  float v;
  double sum;
  struct distribution *d;
  
  values = malloc(sizeof(double) * size);
  if (values == NULL) {
    fprintf(stderr, "Cannot (m)allocate %d = %d * %d bytes\n",
		    sizeof(double) * size,
		    sizeof(double),
		    size);
    return NULL;
  }
  memset(values, 0, size * sizeof(double));
  d = malloc(sizeof(struct distribution));
  if (d == NULL) {
    fprintf(stderr, "Cannot allocate new distribution\n");
    free(values);

    return NULL;
  }

  f = fopen(name, "r");
  done = 0;
  while (!done) {
    done = (fscanf(f, "%d %f\n", &n, &v) != 2);
printf("%d %f\n", n, v);
    if (!done) {
      if (n >= size) {
        fprintf(stderr, "Strange... %d > %d\n", n, size);
	n = size - 1;
      }
      if (n < 0) {
        fprintf(stderr, "%d < 0???\n", n);
        free(values);
	free(d);

	return NULL;
      }
      if (values[n] != 0) {
        fprintf(stderr, "Same value %d 2 times???\n", n);
        free(values);
	free(d);

	return NULL;
      }
      values[n] = v;
    }
  }

  sum = 0;
  for (n = 0; n < size; n++) {
    sum += values[n];
  }
  if (sum == 0) {
    fprintf(stderr, "Wrong file!!!\n");
    free(values);
    free(d);

    return NULL;
  }

  d->values = values;

  if (sum != 1.0) {
    printf("Sum: %f != 1!!!\n", sum);
    /*
    values[sd->max_distr_dim - 1] += (1.0 - sum);
    */
  }

/*
  for (n = 0; n < 20; n++) {
    printf("V[%d]: %f\n", n, values[n]);
  }
*/

  return d;
}
