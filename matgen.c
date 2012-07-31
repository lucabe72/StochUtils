#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <stdint.h>
#include <time.h>
#include <sys/time.h>

#include "pmf.h"
#include "cdf.h"
#include "pmf-file.h"
#include "models.h"
#include "meschac/matrix.h"
#include "meschac/matlab.h"

#define Nc 100000
#define Nz 12

static int Q = 10000;
static int T = 20000;
static int z = 0;
static int d = 0;
static char ofile[50] = "";
static int model = 0;


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
     *     printf("Bad PMF: %d | %f!!!\n", pmf_check(p), pmf_sum(p, 0));
     *         exit(-1);
     *           }*/
    return p;
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


static int opts_parse(int argc, char *argv[])
{
    int opt;
    while ((opt = getopt(argc, argv, "t:q:z:d:o:erl")) != -1) {
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
	case 'd':
	    d = atoi(optarg);
	    break;
	case 'o':
	    strcpy(ofile, optarg);
	    break;
	case 'r':
	    model = RTSSMODEL;
	    break;
	case 'e':
	    model = ETFAMODEL;
	    break;
	case 'l':
	    model = LASTMODEL;
	    break;
	default:		/* ?~@~Y??~@~Y */
	    fprintf(stderr, "Example periodic: %s -o ofile -z n cfile\n",
		    argv[0]);
	    fprintf(stderr, "Example generic: %s -o ofile cfile ufile\n",
		    argv[0]);
	    fprintf(stderr, "Parameters:\n");
	    fprintf(stderr, "-q n:\t Budget value\n");
	    fprintf(stderr, "-t n:\t Server period value\n");
	    fprintf(stderr, "-d n:\t Delta\n");

	    exit(EXIT_FAILURE);
	}
    }
    return optind;
}

void print(MAT *m)
{
unsigned int i,j;
for (i=0; i<m->n; i++){
  for (j=0; j<m->m; j++)
    printf("%f ",m_get_val(m,i,j));
  printf("\n");
}
}

int main(int argc, char *argv[])
{
    int opt;
    struct pmf *c;
    struct pmf *u;
    struct pmf *c1;
    uint64_t time_start, time_end;
    int n, i, j;
    Real prob = 0.0;
    MAT *matrix;

    opt = opts_parse(argc, argv);

    if (strcmp(ofile, "") == 0) {
	fprintf(stderr, "You must specify an output file.mat\n");
	exit(0);
    }
    c = load(argv[opt], Nc);
    c1 = pmf2cdf(c);

//print(c1,"c1");

    time_start = get_time();
    switch (model) {
    case ETFAMODEL:{
	    if (z == 0)
		u = load(argv[opt + 1], Nc);
	    else {
		u = pmf_create(Nc, 0);
		pmf_set(u, z, 1.0);

	    }
	    n = 100;
	    matrix = m_get(n, n);
	    for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
		    prob = prob_efta(i, j, Q, c, u);
		    m_set_val(matrix, i, j, prob);
		}
	    }
	    break;
	}
    case LASTMODEL:{
	    if (d == 0)
		d = Q;
	    n = pmf_max(c) / d * 3;
	    matrix = m_get(n, n);
	    for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
		    prob = prob_last(i, j, Q, z, d, c1);
		    m_set_val(matrix, i, j, prob);
		}
	    }
	    break;
	}
    case RTSSMODEL:{
	    n = 10;
	    matrix = m_get(n, n);
	    for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
		    prob = prob_rtss(i, j, Q, z, c);
		    m_set_val(matrix, i, j, prob);
		}
	    }
	    break;
	}
    default:
	fprintf(stderr, "Choose a model please\n");
	exit(0);
    }
    time_end = get_time();


/*for (i=0; i<n; i++){
     for (j=0; j<n; j++){
       printf("%f ",m_get_val(matrix,i,j));
     }
       printf("\n");
    
}*/
    //print(matrix);
    FILE *fm = fopen(ofile, "w");
    m_save(fm, matrix, "matrix");
    fclose(fm);
    return (0);
}
