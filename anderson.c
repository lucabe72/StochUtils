#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <stdint.h>
#include <time.h>
#include <sys/time.h>

#include "pmf.h"
#include "cdf.h"
#include "pmf-file.h"

#define Nc 100000
#define Nz 12

static int Q = 10000;
static int T = 20000;
static int z = 0;


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
    while ((opt = getopt(argc, argv, "t:q:z:")) != -1) {
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


int main(int argc, char *argv[])
{
    int opt;
    struct pmf *c;
    uint64_t time_start, time_end;
    double var = 0;
    double avg = 0;
    opt = opts_parse(argc, argv);

    time_start = get_time();
    c = load(argv[opt], Nc);
    var = pmf_var(c);
    avg = pmf_avg(c);

    printf("VAR %f\n",var); 
    printf("((%f/(2*%d*(%d-%f)*(1-0.5))+3)*%d\n",var,Q,Q,avg,T);
    printf("%f\n",var/((2*Q*(Q-avg))+3)*T); 
    time_end = get_time();
    
    return (0);
}
