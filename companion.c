#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <stdint.h>
#include <string.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_eigen.h>

#include "pmf.h"
#include "pmf-file.h"

#define Nc 100000

static struct pmf *load(const char *fname, int n)
{
  FILE *f;
  struct pmf *p;

  f = fopen(fname, "r");
  if (!f) {
    perror("FOpen");

    exit(-1);
  }

  p = pmf_create(n, 0);
  n = pmf_read(p, f);
  fclose(f);

  if (pmf_check(p)) {
    printf("Bad PMF: %d | %f!!!\n", pmf_check(p), pmf_sum(p, 0));

    exit(-1);
  }

  return p;
}

static const struct pmf *resample(const struct pmf *c, int q)
{
  struct pmf *c1;
  double sum;
  int i;

  c1 = pmf_create(Nc, 0);
  sum = 0;
  for (i = pmf_min(c); i <= (pmf_max(c) / q) * q + q; i++) {
    sum += pmf_get(c, i);
    if (i % q == 0) {
      pmf_set(c1, i, sum);
      sum = 0;
    }
  }

  return c1;
}


int main(int argc, char *argv[])
{
  const struct pmf *c, *c_orig;
  double *mydata;
  gsl_complex product;
  unsigned int q, size, step, i;
  gsl_matrix_view m;
  gsl_vector_complex *eval; 
  gsl_matrix *evec;
  gsl_eigen_nonsymm_workspace *w;

  if (argc < 3) {
    fprintf(stderr, "Usage: %s <c PMF> <step size> <q>\n", argv[0]);

    return -1;
  }
  c_orig = load(argv[1], Nc);
  step = atoi(argv[2]);
  q = atoi(argv[3]);

  //print(c, "c");
  c = resample(c_orig, step);
  if (q % step) {
    printf("Error: Q mod D = %d mod %d = %d != 0\n", q, step, q % step);

    return -1;
  }

  size = (pmf_max(c) - pmf_min(c)) / step;
  mydata = malloc(size * size * sizeof(double));
  memset(mydata, 0, size * size * sizeof(double));

  for (i = 0; i < size; i++) {
    if (i > 0) {
      int j = i - 1;

      mydata[i + j * size] = 1;
    }

    if ((size - i) * step + pmf_min(c) == q) {
      mydata[i + (size - 1) * size] = (1 - pmf_get(c, ((size - i) * step) + pmf_min(c))) / pmf_get(c, pmf_min(c));
    } else {
      mydata[i + (size - 1) * size] = -pmf_get(c, ((size - i) * step) + pmf_min(c)) / pmf_get(c, pmf_min(c));
    }
  }

  m = gsl_matrix_view_array(mydata, size, size);
  eval = gsl_vector_complex_alloc(size);
  evec = gsl_matrix_alloc(size, size);
  w = gsl_eigen_nonsymm_alloc(size);
  gsl_eigen_nonsymm(&m.matrix, eval, w);
  gsl_eigen_nonsymm_free (w);
     
  //gsl_eigen_nonsymm_sort (eval, evec, GSL_EIGEN_SORT_ABS_ASC);
       
  product = gsl_complex_rect(1, 0);
  printf("Size : %d\n", size);
  for (i = 0; i < size; i++) {
    gsl_complex eval_i = gsl_vector_complex_get (eval, i);
  
    if (gsl_complex_abs(eval_i) < 1 - 0.00000000001) {
      gsl_complex c = gsl_complex_sub(gsl_complex_rect(1, 0), eval_i);

      product = gsl_complex_mul(product, c);
    }
//    printf ("%d eigenvalue = %g %g %g\n", i, GSL_REAL(eval_i), GSL_IMAG(eval_i), gsl_complex_abs(eval_i));
//    printf("Product = %lf %lf\n", GSL_REAL(product), GSL_IMAG(product));
  }
  printf("P0 = %lf %lf\n", GSL_REAL(product), GSL_IMAG(product));
     
  gsl_vector_complex_free (eval);
  gsl_matrix_free (evec);

  return 0;
}

