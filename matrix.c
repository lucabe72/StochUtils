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

static int Q = 10000;
static int P = 20000;
static int T = 70000;
static int samples;


//#define DEBUG 1

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

static int opts_parse(int argc, char *argv[])
{
  int opt;

    while ((opt = getopt(argc, argv, "t:q:T:s:")) != -1) {
      switch (opt) {
        case 'q':
          Q = atoi(optarg);
	  break;
	case 't':
	  P = atoi(optarg);
	  break;
	case 'T':
	  T = atoi(optarg);
	  break;
	case 's':
	  samples = atoi(optarg);
	  break;
	default: /* ?~@~Y??~@~Y */
	  fprintf(stderr, "Usage: %s [-t nsecs] [-n] name\n", argv[0]);
	  exit(EXIT_FAILURE);
	}
    }

  return optind;
}



int main(int argc, char *argv[])
{
  struct pmf *c, *cdfc, *u;
  double avgc;
  int n,maxv,i,j;
  int qmin, qmax, qstep,t,forward,back;
  MAT *mat,*B0,*A0,*A1,*A2,*R;
  VEC *X0;
#ifdef DEBUG
  VEC *tmp;
#endif
  int opt;

  opt = opts_parse(argc, argv);
  c = load(argv[opt], Nc);
  u = load(argv[opt + 1], Nz);

  /*c = load(argv[1], 8001);
  u = load(argv[2], 12);*/
  n=pmf_max(u);

  avgc = pmf_avg(c);

  qmin = avgc / pmf_max(u);
  qmax = pmf_max(c) / pmf_min(u);
  qstep = (qmax - qmin) / 10;
#ifdef DEBUG
  printf("Q > %f / %d = %d\n", avgc, pmf_max(u), qmin - 1);
  printf("Q <= %d / %d = %d\n", pmf_max(c), pmf_min(u), qmax);
#endif
  //q= qmin - 1;
  //Q=atoi(argv[3]);
  t=atoi(argv[4]);

  //1. compute cdf of U
  cdfc=cdf_pmf(c);

  forward=ceil((1.0*pmf_max(c))/(1.0*Q))-pmf_min(u);
  back=-1*(ceil((1.0*pmf_min(c))/(1.0*Q))-pmf_max(u));

  //2. compute maxvalue
  maxv=max(forward,back);
  if (maxv<=0) maxv=1;

  mat=m_get(maxv*3,maxv*3);

  //3. compute matrix
  for (i=0; i<maxv*3; i++){
     for (j=0; j<maxv*3; j++){
	m_set_val(mat,i,j,matrix_prob3(i,j,Q,T,cdfc,u));
	}
     }

  //m_output(mat);
#ifdef DEBUG
  for (i=0; i<maxv*3; i++){
     for (j=0; j<maxv*3; j++)
       printf("%f ",m_get_val(mat, i,j) );
     printf("\n");
  }

  tmp=v_get(maxv);
  for (i=0; i<maxv; i++){
    get_row(mat,i,tmp);
    if (v_sum(tmp)!=1)
     {
     fprintf(stderr,"ROW SUM !=1\n");
     }
  }
#endif

  //4. Extract submatrix
  B0=m_get(maxv,maxv);
  A0=m_get(maxv,maxv);
  A1=m_get(maxv,maxv);
  A2=m_get(maxv,maxv);
  compute_matrixes(mat,maxv,B0,A0,A1,A2);
#ifdef DEBUG
  m_output(B0);
  m_output(A0);
  m_output(A1);
  m_output(A2);
#endif
  //5. Compute matrix R
  R=m_get(maxv,maxv);
  R=computeR(R,A0,A1,A2,0.001);
#ifdef DEBUG
  m_output(R);
#endif
  //6. Compute X0 from R
  X0=v_get(maxv);
  X0=computeX0(R,B0,A2,X0);
 #ifdef DEBUG
  v_output(X0);
  printf("PROB %f\n",v_sum(X0));
#endif

  for (i=1; i<=maxv; i++)
	printf("P{r = %i} %f\n",T*i,v_get_val(X0,i-1));
  if (v_sum(X0)!=1)
	printf("P{r >= %i} %f\n",T*i,1-v_sum(X0));


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
