#include <stdlib.h>
#include <stdio.h>

#include "pmf.h"
#include "pmf-file.h"
#include "generic.h"

#define N 100

int main(int argc, char *argv[])
{
  FILE *f;
  struct pmf *t;
  double *z;
  int n, i, ts;

  if (argc < 3) {
    fprintf(stderr, "Usage: %s <i PMF> ts\n", argv[0]);

    return -1;
  }
  f = fopen(argv[1], "r");
  if (!f) {
    perror("FOpen");

    return -1;
  }
  ts = atoi(argv[2]);

  t = pmf_create(N * ts, 0);
  n = pmf_read(t, f);


  if (pmf_check(t)) {
    printf("Bad PMF: %d | %f!!!\n", pmf_check(t), pmf_sum(t, 0));

    return -1;
  }
  z = generic_transform(t, ts);
  for (i = 0; i < N; i++) {
    printf("%d %14.15f\n", i, z[i]);
  }

  return 0;
}
