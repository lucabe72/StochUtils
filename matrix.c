#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <stdint.h>
#include <math.h>

#include "pmf.h"
#include "pmf-file.h"
#include "gamma.h"
#include "y.h"
#include "cdf.h"
#include "qbdm.h"
#include "meschac/matrix.h"

#define Nc 30000
#define Nz 20

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

/*static void print(const struct pmf *p, char *v)
{
  int i;

  for (i = pmf_min(p); i <= pmf_max(p); i++) {
    if (pmf_get(p, i) > 1e-10) {
      printf("P{%s = %d }= %f\n", v, i, pmf_get(p, i));
    }
  }
  printf("P{%s > %d} = %f\n", v, pmf_max(p), pmf_tail(p));
}
*/

int main(int argc, char *argv[])
{
  struct pmf *c, *cdfc, *u;
  double avgc;
  int n,maxv,i,j;
  int qmin, qmax, q, qstep;
  MAT *mat,*B0,*A0,*A1,*A2,*R;
  VEC *X0;
  if (argc < 4) {
    fprintf(stderr, "Usage: %s <C PMF> <T PMF> Q\n", argv[0]);

    return -1;
  }
  c = load(argv[1], 8001);
  u = load(argv[2], 12);
  n=pmf_max(u);

  avgc = pmf_avg(c);

  qmin = avgc / pmf_max(u);
  qmax = pmf_max(c) / pmf_min(u);
  qstep = (qmax - qmin) / 10;
  printf("Q > %f / %d = %d\n", avgc, pmf_max(u), qmin - 1);
  printf("Q <= %d / %d = %d\n", pmf_max(c), pmf_min(u), qmax);
  //q= qmin - 1;
  q=atoi(argv[3]);
  printf("Q = %d\n", q);
  //1. compute cdf of U
  cdfc=cdf_pmf(c);

  //2. compute maxvalue
  //maxv=floor(pmf_max(c)/q)-pmf_min(u);
  //maxv=ceil(pmf_max(c)*1.0/q*1.0)-pmf_min(u);
  maxv=1;
  //maxv=ceil(pmf_max(c)/q);
  printf("%i %i %i MAX %i\n",pmf_max(c),pmf_min(u),q,maxv);
//  exit(0);
  mat=m_get(maxv*3,maxv*3);

  //3. compute matrix
  for (i=0; i<maxv*3; i++){
     for (j=0; j<maxv*3; j++){
	m_set_val(mat,i,j,matrix_prob2(i,j,q,cdfc,u));
	}
     }

  //m_output(mat);

  for (i=0; i<maxv*3; i++){
     for (j=0; j<maxv*3; j++)
       printf("%f ",m_get_val(mat, i,j) );
     printf("\n");
  }

  //4. Extract submatrix
  B0=m_get(maxv,maxv);
  A0=m_get(maxv,maxv);
  A1=m_get(maxv,maxv);
  A2=m_get(maxv,maxv);
  compute_matrixes(mat,maxv,B0,A0,A1,A2);

  m_output(B0);
  m_output(A0);
  m_output(A1);
  m_output(A2);
  //5. Compute matrix R
  R=m_get(maxv,maxv);
  R=computeR(R,A0,A1,A2,0.001);
  m_output(R);
  //6. Compute X0 from R
  X0=v_get(maxv);
  X0=computeX0(R,B0,A2,X0);
  v_output(X0);

  //7. Free the allocated memory
  m_free(B0);
  m_free(A0);
  m_free(A1);
  m_free(A2);
  m_free(mat);
  m_free(R);
  v_free(X0);
  return 0;
}
