#include <stdio.h>

#include "pmf.h"
#include "pmf-file.h"

#define N 100

int main(int argc, char *argv[])
{
  FILE *f;
  struct pmf *p;
  int n, i;

  f = fopen(argv[1], "r");
  if (!f) {
    perror("FOpen");

    return -1;
  }

  p = pmf_create(N, 0);
  n = pmf_read(p, f);

  pmf_set_samples(p, 1000000);

  if (pmf_check(p)) {
    printf("Bad PMF: %d | %f!!!\n", pmf_check(p), pmf_sum(p, 0));

    return -1;
  }
  for (i = 0; i < N; i++) {
    printf("P{x = %d} = %f\n", i, pmf_get(p, i));
  }
  printf("P{x >= %d} = %f\n", N, pmf_tail(p));

  return 0;
}
