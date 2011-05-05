#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <stdint.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>

#include "pmf.h"
#include "pmf-file.h"
#include "gamma.h"
#include "y.h"
#include "cdf.h"
#include "qbdm.h"
#include "meschac/matrix.h"

#define Nc 21000
#define Nz 12

static int Q = 10000;
static int T = 20000;
int ps=0;

//#define DEBUG 1

static inline uint64_t get_time(void)
{
  struct timeval tv;
  uint64_t res;

  gettimeofday(&tv, NULL);
  res = tv.tv_sec;
  res = res * 1000000 + tv.tv_usec;

  return res;
}



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

  /*if (pmf_check(p)) {
    printf("Bad PMF: %d | %f!!!\n", pmf_check(p), pmf_sum(p, 0));

    exit(-1);
  }*/

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

    while ((opt = getopt(argc, argv, "t:q:p")) != -1) {
      switch (opt) {
        case 'q':
          Q = atoi(optarg);
	  break;
	case 't':
	  T = atoi(optarg);
	  break;
	case 'p':
	  ps = 1;
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
  double avgc,sum;
  int n,maxv,i,j;
  int qmin, qmax, qstep,t,forward,back;
  uint64_t t1, t2, t3;
  MAT *mat,*B0,*A0,*A1,*A2,*R,*matp,*B0p,*A0p,*A1p,*A2p,*Rp;
  VEC *X0,*X0p,*X1;
#ifdef DEBUG
  VEC *tmp;
#endif
  int opt;
  int index;
  float tmp_prob;

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

  forward=ceil((1.0*pmf_max(c))/(1.0*Q));
  back=-1*ceil((1.0*pmf_min(c))/(1.0*Q));


  //2. compute maxvalue
  maxv=max(forward,back);
  if (maxv<=0) maxv=1;

  mat=m_get(maxv*3,maxv*3);
  t1=get_time();
  //3. compute matrix
  for (i=0; i<maxv*3; i++){
     for (j=0; j<maxv*3; j++){
	m_set_val(mat,i,j,matrix_prob_ts(i,j,Q,cdfc,u));
	}
     }
if (ps){
  matp=m_get(maxv*3,maxv*3);
  //3b. create a pessimistic matrix
  for (i=0; i<maxv*3; i++){
     sum=0;
     for (j=0; j<i; j++)
	sum+=m_get_val(mat,i,j);
     if (j>0)
     m_set_val(matp,i,j-1,sum);
     for (j=i; j<maxv*3; j++)
        m_set_val(matp,i,j,m_get_val(mat,i,j));
     }
}
#ifdef DEBUG
  if (ps){
  printf("PESSIMISTIC MATRIX: \n");
  for (i=0; i<maxv*3; i++){
     for (j=0; j<maxv*3; j++)
       printf("%f ",m_get_val(matp, i,j) );
     printf("\n");
  }
     printf("\n");
}
  printf("MATRIX: \n");
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
     //fprintf(stderr,"ROW SUM !=1\n");
     }
  }
#endif
  index=maxv;
  tmp_prob=m_get_val(mat,0,index);
  while (tmp_prob==0){
    index--;
    tmp_prob=m_get_val(mat,0,index);
  }
  forward=index+1;
  //printf("FORWARD!! %i\n",forward);
  index=1;
  tmp_prob=m_get_val(mat,index,0);
  while (tmp_prob>0){
    index++;
    tmp_prob=m_get_val(mat,index,0);
  }
  back=index-1;
  maxv=max(forward,back);
  //printf("%i \n",maxv);

  //4. Extract submatrix
  B0=m_get(maxv,maxv);
  A0=m_get(maxv,maxv);
  A1=m_get(maxv,maxv);
  A2=m_get(maxv,maxv);
  compute_matrixes(mat,maxv,B0,A0,A1,A2);
  if (ps){
  B0p=m_get(maxv,maxv);
  A0p=m_get(maxv,maxv);
  A1p=m_get(maxv,maxv);
  A2p=m_get(maxv,maxv);
  compute_matrixes(matp,maxv,B0p,A0p,A1p,A2p);
  }
#ifdef DEBUG
  m_output(B0);
  m_output(A0);
  m_output(A1);
  m_output(A2);
#endif
  //5. Compute matrix R
  R=m_get(maxv,maxv);
  R=computeR(R,A0,A1,A2,0.001);
  if (ps){
    Rp=m_get(maxv,maxv);
    Rp=computeR(Rp,A0p,A1p,A2p,0.001);
  }
  t2=get_time();
#ifdef DEBUG
  m_output(R);
#endif
  //6. Compute X0 from R
  X0=v_get(maxv);
  X1=v_get(maxv);
  X0=computeX0(R,B0,A2,X0);
  if (ps){
    X0p=v_get(maxv);
    X0p=computeX0(Rp,B0p,A2p,X0p);
  }
  t3=get_time();
#ifdef DEBUG
  v_output(X0);
  printf("PROB %f\n",v_sum(X0));
  if (ps){
  v_output(X0p);
  printf("PROBp %f\n",v_sum(X0p));
  }
#endif
  int cc=0;
  double prob=0;
  for (i=pmf_min(u); cc<maxv; i++){
        prob+=v_get_val(X0,cc);
	printf("P{d < %i} %f\n",T*(cc+1),prob);
	cc++;
  }
  if (v_sum(X0)!=1){
    int start=cc;
    vm_mlt(R,X0,X1);
    for (cc=0; cc<maxv; cc++){
      prob+=v_get_val(X1,cc);
      printf("P{d < %i} %f\n",T*(start+cc+1),prob);
    }
    printf("P{d >= %i} %f\n",T*(start+cc+1),1-v_sum(X0)-v_sum(X1));
}
  if (ps){
  printf("\nRESULTS FOR CONSERVATIVE APPROSSIMATION:\n");
  cc=0;
  prob=0;
  for (i=pmf_min(u); cc<=maxv; i++){
        prob+=v_get_val(X0p,cc);
	printf("P{d < %i} %f\n",T*(cc+1),prob);
	cc++;
  }
  if (v_sum(X0p)!=1)
	printf("P{d >= %i} %f\n",T*cc,1-v_sum(X0p));
  }
  //7. Free the allocated memory
  m_free(B0);
  m_free(A0);
  m_free(A1);
  m_free(A2);
  m_free(mat);
  if(ps) m_free(matp);
  m_free(R);
  v_free(X0);
  if (ps)v_free(X0p);

  printf("Ctime: %Lu\n", t3 - t1);
  printf("Ctime: %Lu %Lu\n", t2 - t1, t3 - t2);

  return 0;
}
