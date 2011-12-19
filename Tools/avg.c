#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <stdint.h>

#include "pmf.h"
#include "pmf-file.h"
#include "gamma.h"
#include "y.h"

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

int main(int argc, char *argv[])
{
  struct pmf *c, *z;
  double avgc, avgz;
  int qmin, qmax, q, qstep;

  if (argc < 3) {
    fprintf(stderr, "Usage: %s <C PMF> <Z PMF>\n", argv[0]);

    return -1;
  }
  c = load(argv[1], Nc);
  z = load(argv[2], Nz);

  print(z, "z");
  print(c, "c");

  avgc = pmf_avg(c);
  avgz = pmf_avg(z);

  qmin = avgc / avgz;
  qmax = pmf_max(c) / pmf_min(z);
  qstep = (qmax - qmin) / 10;
  printf("Q > %f / %f = %d\n", avgc, avgz, qmin - 1);
  printf("Q <= %d / %d = %d\n", pmf_max(c), pmf_min(z), qmax);


  for (q = qmin; q <= qmax; q += qstep) {
    char name[16];
    struct pmf *y;
    double gamma;
    const double eps = 5e-5;

    y = compute(c, z, q);
    for (gamma = 1 + 1e-5; gamma < 1 + 50 * eps; gamma += eps) {
      printf("Q: %d\t Gamma: %f \t sum: %f\n", q, gamma, compute_sum(y, gamma));
    }
    sprintf(name, "Y%d", q);
    print(y, name);
  }

  return 0;
}
