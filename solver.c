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
#include "meschac/matlab.h"

#define Nc 21000
#define Nz 12

static char ifile[50]="";
static int response=0;
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


static int opts_parse(int argc, char *argv[])
{
  int opt;

    while ((opt = getopt(argc, argv, "i:r")) != -1) {
      switch (opt) {
        case 'i':
          strcpy(ifile,optarg);
          break;
        case 'r':
          response=1;
          break;
	default: /* ?~@~Y??~@~Y */
	  fprintf(stderr, "Usage: %s [-t nsecs] [-n] name\n", argv[0]);
	  exit(EXIT_FAILURE);
	}
    }

  return optind;
}

void print(MAT *m)
{
unsigned int i,j;
unsigned int dim=m->n;
if (dim>25) dim=25;
for (i=0; i<dim; i++){
  for (j=0; j<dim; j++)
    printf("%f ",m_get_val(m,i,j));
  printf("\n");
}


}

int main(int argc, char *argv[])
{
  uint64_t t1, t2, t3;
  char *name;
  MAT *mat,*B0,*A0,*A1,*A2,*R;
  VEC *X0,*X1;
#ifdef DEBUG
  VEC *tmp;
#endif
  int maxv,back,forward;
  int opt;
  int index,i;
  float tmp_prob;

  opt = opts_parse(argc, argv);
  if (strcmp(ifile,"")==0){
     fprintf(stderr,"You must specify an input file\n");
     exit(0);  
  }
  FILE *fm=fopen(ifile,"r");
  mat=m_load(fm,&name);
  maxv=mat->n;
  
  //print(mat);
  if (m_get_val(mat,0,0)==1.0)
  {
  printf("Prob: %f\n",1.000000);
  printf("Ctime: %u\n", 0);
  printf("Ctime: %u %u\n", 0, 0);
  exit(0);
  } 
  index=maxv-1;
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
    if (index>=maxv) break;
    tmp_prob=m_get_val(mat,index,0);
  }
  back=index-1;
  maxv=max(forward,back);
  
  //printf("n=%i max(%i,%i)=%i\n",mat->n,back,forward,maxv);

/*  c = load(argv[opt], Nc);
  u = load(argv[opt + 1], Nz);

  mat=m_get(maxv*3,maxv*3);*/
  t1=get_time();
  //3. compute matrix

#ifdef DEBUG 
  tmp=v_get(maxv);
  for (i=0; i<maxv; i++){
    get_row(mat,i,tmp);
    if (v_sum(tmp)!=1)
     {
     //fprintf(stderr,"ROW SUM !=1\n");
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
  t2=get_time();
#ifdef DEBUG
  m_output(R);
#endif
  //6. Compute X0 from R
  X0=v_get(maxv);
  X1=v_get(maxv);
  X0=computeX0(R,B0,A2,X0);
  t3=get_time();
#ifdef DEBUG
  v_output(X0);
  printf("PROB %f\n",v_sum(X0));
#endif
  if (response){
  int cc=0;
  double prob=0;
  for (i=0; cc<maxv; i++){
        prob+=v_get_val(X0,cc);
	printf("P{d < %iT} %f\n",(cc+1),prob);
	cc++;
  }
  if (v_sum(X0)!=1){
    int start=cc;
    vm_mlt(R,X0,X1);
    for (cc=0; cc<maxv; cc++){
      prob+=v_get_val(X1,cc);
      printf("P{d < %iT} %f\n",(start+cc+1),prob);
    }
    printf("P{d >= %iT} %f\n",(start+cc+1),1-v_sum(X0)-v_sum(X1));
}
  }
  printf("Prob: %f\n",v_get_val(X0,0));
  //7. Free the allocated memory
  m_free(B0);
  m_free(A0);
  m_free(A1);
  m_free(A2);
  m_free(mat);
  m_free(R);
  v_free(X0);

  printf("Ctime: %Lu\n", t3 - t1);
  printf("Ctime: %Lu %Lu\n", t2 - t1, t3 - t2);

  return 0;
}
