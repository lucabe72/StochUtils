#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>

#include "pmf.h"
#include "pmf-file.h"
#include "pmf-sample.h"

#define N 500000

static int max = 1000;

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
  while (v > 0) {
    v -= pmf_get(c, res++);
  }

  return --res;
}

static int opts_parse(int argc, char *argv[])
{
  int opt;

  while ((opt = getopt(argc, argv, "s:m:")) != -1) {
    switch (opt) {
      case 'm':
        max = atoi(optarg);
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
  double sum = 0;

  optind = opts_parse(argc, argv);
printf("OI: %d\n", optind);
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
  for (i = 0; i < pmf_max(p) + 1; i++) {
    fprintf(stderr, "P{x = %d} = %f\n", i, pmf_get(p, i));
  }

  res = pmf_create(pmf_max(p) + 1, 0);
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

  return 0;
}

