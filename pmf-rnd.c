#include <stdio.h>
#include <stdlib.h>

#include "pmf.h"
#include "pmf-file.h"
#include "pmf-sample.h"

#define N 100000

static int max = 1000;

double rnd(void)
{
  double v;

  v = rand();

  return v / (double)RAND_MAX;
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

int main(int argc, char *argv[])
{
  FILE *f;
  struct pmf *p, *res;
  int n, i;
  const double epsilon =  1e-10;

  if (argc < 2) {
    fprintf(stderr, "Usage: %s <PMF>\n", argv[0]);

    return -1;
  }
  if (argc > 2) {
    max = atoi(argv[2]);
    if (argc > 3) {
      srand(atoi(argv[3]));
    }
  }
  f = fopen(argv[1], "r");
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

