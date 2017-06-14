#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "pmf.h"
#include "pmf-file.h"
#include "pmf-sample.h"

#define N 500000

static int max = 1000;
static unsigned long long int eps_c1;

double rnd(void)
{
  double v;

  v = rand();

  return v / (double)RAND_MAX;
}

void rnd_init(unsigned long int seed)
{
  srand(seed);
}

int pmf_rnd(const struct pmf *c)
{
  double v;
  int res;

  v = rnd();

  res = 0;
  while ((v > 0) && res <= pmf_max(c)) {
    v -= pmf_get(c, res++);
  }

  if (v > 0) {
    return N - 1;
  }

  return --res;
}

static long long unsigned int e(int n)
{
  long long unsigned int res = 1;
  int i;

  for (i = 0; i < n; i++) {
    res *= 10ULL;
  }

  return res;
}

static int opts_parse(int argc, char *argv[])
{
  int opt;

  while ((opt = getopt(argc, argv, "s:m:e:")) != -1) {
    switch (opt) {
      case 'e':
	eps_c1 = e(atoi(optarg));
        break;
      case 'm':
        max = atoi(optarg);
        break;
      case 's':
        rnd_init(atoi(optarg));
        break;
      default: /* ’?’ */
        fprintf(stderr, "Wrong parameter %c\n", opt);
        exit(EXIT_FAILURE);
     }
  }

  return optind;
}

int main(int argc, char *argv[])
{
  FILE *f;
  struct pmf *p, *res;
  int n, i, optind;
  const double epsilon =  1e-10;
  double sum;

  optind = opts_parse(argc, argv);
  if (argc - optind < 1) {
    fprintf(stderr, "Usage: %s <PMF>\n", argv[0]);

    return -1;
  }
  f = fopen(argv[optind], "r");
  if (!f) {
    perror("FOpen");

    return -1;
  }

  p = pmf_create(N, 0);
  n = pmf_read(p, f);

  if (pmf_check(p)) {
    printf("Bad PMF: %d | %f!!!\n", pmf_check(p), pmf_sum(p, 0));

    return -1;
  }
  printf("Read %d samples\n", n);
  for (i = 0; i < pmf_max(p) + 1; i++) {
    if (pmf_get(p, i) > epsilon) {
      if (eps_c1) {
        /* fprintf(stderr, "%f -> %f\n",  pmf_get(p, i), pmf_get(p, i) * ((double)(eps_c1 - 1) / (double)eps_c1)); */
        pmf_set(p, i, pmf_get(p, i) * ((double)(eps_c1 - 1) / (double)eps_c1));
      }
      fprintf(stderr, "P{x = %d} = %f\n", i, pmf_get(p, i));
    }
  }
  fprintf(stderr, "PMF Check[1]: %d\n", pmf_check(p));
  if (eps_c1) {
    p->tail = 1.0 / (double)eps_c1;
  }
  fprintf(stderr, "PMF Check[2]: %d\n", pmf_check(p));

  res = pmf_create(/*pmf_max(p) + 1*/ N, 0);
  for (i = 0; i < max; i++) {
    int val;

    val = pmf_rnd(p);
    //printf("%d\n", val);
    pmf_collect(res, val);
  }
  pmf_normalise(res);
  for (i = 0; i < pmf_max(res) + 1; i++) {
    if (pmf_get(res, i) > epsilon) {
      printf("%d\t %f\n", i, pmf_get(res, i));
    }
  }
  if (pmf_max(res) <= pmf_max(p)) {
    sum = p->tail;
  } else {
    sum = 0;
  }
  for (i = pmf_max(res) + 1; i < pmf_max(p) + 1; i++) {
    sum += pmf_get(p, i);
  }
  fprintf(stderr, "#Max: %d\n", pmf_max(res));
  fprintf(stderr, "P{x > %d} <= %1.20f   (= %1.20f)\n", pmf_max(res), 1.0 / max, sum);

  double eps=(pmf_avg(res)+(pmf_std(res) / sqrt(max))*3)/(pmf_max(res));  
  fprintf(stderr, "P{x > %d} <= %1.20f \n",pmf_max(p) ,eps);
  int rep=pow(((pmf_std(res)*3) / (eps*pmf_max(res)-pmf_avg(res))),2);
  fprintf(stderr, "n = %i \n",rep);


  pmf_free(p);
  pmf_free(res);

  return 0;
}

