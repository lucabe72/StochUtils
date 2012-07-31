#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <stdint.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>

#include "meschac/matrix.h"
#include "meschac/matlab.h"


static char ifile[50]="";

static int opts_parse(int argc, char *argv[])
{
  int opt;
    while ((opt = getopt(argc, argv, "i:")) != -1) {
      switch (opt) {
        case 'i':
          strcpy(ifile,optarg);
          break; 
        default: /* ?~@~Y??~@~Y */
          fprintf(stderr, "Usage: %s [-t nsecs] [-n] name\n", argv[0]);
          exit(EXIT_FAILURE);
        }
    }
  return optind;
}

static inline uint64_t get_time(void)
{
  struct timeval tv;
  uint64_t res;

  gettimeofday(&tv, NULL);
  res = tv.tv_sec;
  res = res * 1000000 + tv.tv_usec;

  return res;
}


int main(int argc, char *argv[])
{
  int h,n,opt;
  double ap,result;
  MAT *mat;
  VEC *vec=NULL;
  char *name;
  uint64_t t1, t2; 
 
  opt = opts_parse(argc, argv);
  
  if (strcmp(ifile,"")==0){
     fprintf(stderr,"You must specify an input file\n");
     exit(0);  
  }
  FILE *fm=fopen(ifile,"r");
  mat=m_load(fm,&name);
  t1=get_time();
  vec=get_row(mat,1,vec);
  
  ap=v_get_val(vec,0);
  n=vec->dim;
  result=1;
  for (h=2; h<n+1; h++)
    {
    result -= (h-1) * v_get_val(vec,h)/v_get_val(vec,0);
    }
  t2=get_time();
  printf("Prob: %.6f\n",result);
  printf("Ctime: %Lu\n", t2 - t1);
  return 0;
}
