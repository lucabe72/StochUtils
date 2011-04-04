#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>

#include "pmf.h"

static int offs = 100;
static int size = 500;
static int samples = 10000;
static int verbose;

static void pmf_rnd(struct pmf *p, int offs, int size, int samples)
{
  int i;
  double sum = 0;

  for (i = offs; i < offs + size; i++) {
    unsigned long long int val;

    val = rand();
    val = (val* samples) / RAND_MAX;

    pmf_set(p, i, val);
    sum += val;
  }

  for (i = offs; i < offs + size; i++) {
    pmf_set(p, i, pmf_get(p, i) / sum);
  }
}

static int opts_parse(int argc, char *argv[])
{
  int opt;

  while ((opt = getopt(argc, argv, "o:s:S:v")) != -1) {
    switch (opt) {
      case 'o':
        offs = atoi(optarg);
        break;
      case 's':
        size = atoi(optarg);
        break;
      case 'S':
        samples = atoi(optarg);
        break;
      case 'v':
        verbose = 1;
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
  struct pmf *p;
  int i;

  opts_parse(argc, argv);
  p = pmf_create(size + offs, 0);
  pmf_rnd(p, offs, size, samples);

  if (pmf_check(p)) {
    printf("Bad PMF: %d | %f!!!\n", pmf_check(p), pmf_sum(p, 0));

    return -1;
  }
  for (i = 0; i <= pmf_max(p); i++) {
    if (verbose) {
      printf("P{x = %d} = %f\n", i, pmf_get(p, i));
    } else {
      printf("%d\t%f\n", i, pmf_get(p, i));
    }
  }
  if (verbose) {
    printf("P{x >= %d} = %f\n", pmf_max(p), pmf_tail(p));
  }

  return 0;
}
