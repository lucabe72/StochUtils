#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <stdint.h>

#include "pmf.h"
#include "pmf-file.h"

#define Nc 100000
#define N 2

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

  if (pmf_check(p)) {
    printf("Bad PMF: %d | %f!!!\n", pmf_check(p), pmf_sum(p, 0));

    exit(-1);
  }

  return p;
}

static void print(const struct pmf *p, char *v)
{
  int i;

  for (i = pmf_min(p); i <= pmf_max(p); i++) {
    if (pmf_get(p, i) > 1e-10) {
      printf("P{%s = %d }= %f\n", v, i, pmf_get(p, i));
    }
  }
  printf("P{%s > %d} = %f\n", v, pmf_max(p), pmf_tail(p));
}

static const struct pmf *resample(const struct pmf *c, int q)
{
  struct pmf *c1;
  double sum;
  int i;

  c1 = pmf_create(Nc, 0);
  sum = 0;
  for (i = pmf_min(c); i <= (pmf_max(c) / q) * q + q; i++) {
    sum += pmf_get(c, i);
    if (i % q == 0) {
      pmf_set(c1, i, sum);
      sum = 0;
    }
  }

  return c1;
}

int main(int argc, char *argv[])
{
  struct pmf *c;
  const struct pmf *c1;
  int q;

  if (argc < 3) {
    fprintf(stderr, "Usage: %s <q> <c PMF>\n", argv[0]);

    return -1;
  }
  q = atoi(argv[1]);
  c = load(argv[2], Nc);

  print(c, "c");
  c1 = resample(c, q);
  print(c1, "c1");

  return 0;
}
