#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <stdint.h>

#include "pmf.h"
#include "pmf-file.h"
#include "cdf.h"

#define Nc 30000
#define Nz 20

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

//    exit(-1);
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

int main(int argc, char *argv[])
{
  struct pmf *c, *fc;

  if (argc < 2) {
    fprintf(stderr, "Usage: %s <PMF>\n", argv[0]);

    return -1;
  }
  c = load(argv[1], Nc);

  print(c, "c");
  
  fc = pmf2cdf(c);
  print(fc, "C");

  return 0;
}

