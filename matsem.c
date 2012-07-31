#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <stdint.h>
#include <time.h>
#include <sys/time.h>

#include "meschac/matrix.h"
#include "meschac/matlab.h"

#define Nc 200
#define Nz 12

static int Q = 10000;
static int T = 20000;
static int z=4;
int ps=0;
char ifile[50]="";
char ofile[50]="";

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
    while ((opt = getopt(argc, argv, "t:q:z:o:i:p")) != -1) {
      switch (opt) {
        case 'q':
          Q = atoi(optarg);
          break;
        case 't':
          T = atoi(optarg);
          break;
        case 'z':
          z = atoi(optarg);
          break;
        case 'p':
          ps = 1;
          break;
        case 'i':
          strcpy(ifile,optarg);
          break;
        case 'o':
          strcpy(ofile,optarg);
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
int opt;
uint64_t time_start,time_end;
int n,i,j;
double sum;
char *name;
MAT *matrix,*m_red;

opt = opts_parse(argc, argv);
if (strcmp(ifile,"")==0)
  {
  fprintf(stderr,"You must specify an input file\n");
  exit(0);
  }
if (strcmp(ofile,"")==0)
  {
  fprintf(stderr,"You must specify an output file\n");
  exit(0);
  }
FILE *rm=fopen(ifile,"r");

FILE *fm=fopen(ofile,"w");
matrix=m_load(rm,&name);
fclose(rm);


n=matrix->n;

     

m_red=m_get(n,n);

time_start=get_time();
for (i=0; i<n; i++){
     sum=0;
     for (j=0; j<i; j++){
        sum+=m_get_val(matrix,i,j);
        m_set_val(m_red,i,j,0.0);
     }
     if (j>0)
       m_set_val(m_red,i,j-1,sum);
     for (j=i; j<n; j++){
        m_set_val(m_red,i,j,m_get_val(matrix,i,j));
     } 
}
time_end=get_time();

m_save(fm,m_red,"matrice");
fclose(fm);
m_free(matrix);
m_free(m_red);
return (0);
}
