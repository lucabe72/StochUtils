#include <sys/time.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <stdint.h>

#include "pmf.h"
#include "y.h"
#include "z.h"
#include "v.h"
#include "dl.h"
#include "gamma.h"
#include "pmf-file.h"
#include "pmf-modify.h"

#define C 20000
#define Nc 30000
#define Nz 20
static int Q = 10000;
static int P = 20000;
static int T = 70000;
static int samples;

static inline uint64_t get_time(void)
{
  struct timeval tv;
  uint64_t res;

  gettimeofday(&tv, NULL);
  res = tv.tv_sec;
  res = res * 1000000 + tv.tv_usec;

  return res;
}

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

static void print(const struct pmf *p, int n, int offs, char v)
{
  int i;

  for (i = offs; i < n + offs; i++) {
    if (pmf_get(p, i) > 1e-10) {
      printf("P{%c = %d} = %f\n", v, i, pmf_get(p, i));
    }
  }
  printf("P{%c >= %d} = %f\n", v, n, pmf_tail(p));
}

static int opts_parse(int argc, char *argv[])
{
  int opt;

  while ((opt = getopt(argc, argv, "t:q:T:s:")) != -1) {
    switch (opt) {
      case 'q':
        Q = atoi(optarg);
        break;
      case 't':
        P = atoi(optarg);
        break;
      case 'T':
        T = atoi(optarg);
        break;
      case 's':
        samples = atoi(optarg);
        break;
      default: /* ’?’ */
        fprintf(stderr, "Usage: %s [-t nsecs] [-n] name\n", argv[0]);
        exit(EXIT_FAILURE);
     }
  }

  return optind;
}

int main(int argc, char *argv[])
{
  struct pmf *t, *c, *z, *y, *gamma_pmf, *v, *dl;
  double gamma;
  uint64_t t1, t2;
  int opt;

  opt = opts_parse(argc, argv);
  c = load(argv[opt], Nc);
  t = load(argv[opt + 1], Nz);
  z = z_generate(t, P);
 
  print(z, Nz, 0, 'z');

  //c1 = pmf_restrict(c, 20000);
  if (samples) pmf_set_samples(c, samples);
  print(c, Nc, 0, 'c');

  t1 = get_time();
  y = compute(c, z, Q, 20 * Nc);
  gamma = get_gamma(y, T);
  gamma_pmf = gamma_generate(gamma, T, /*T / 1000*/1);
#if 0
  printf("C = %d\n", C);
  y = compute(c, z, Q, 20 * Nc, C, Nz);
  print(y, 20 * Nc,  -20 * Nc / 2, 'y');
  //pmf_set_samples(p, 1000000);
#endif
  v = v_compute(gamma_pmf, c, T, Q / 10);

  dl = stochdl_compute(v, Q, P);
  t2 = get_time();
  print(y, 20 * Nc,  -20 * Nc / 2, 'y');
  printf("Gamma: %f\n", gamma);
//  print(gamma_pmf, T, 0, 'g');
  print(v, T, 0, 'v');
  print(dl, T / Q * P, 0, 'd');
  printf("Ctime: %Lu\n", t2 - t1);

  return 0;
}
