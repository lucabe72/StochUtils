#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>

#include "pmf.h"
#include "pmf-file.h"

#include "distr.h"
#include "driver.h"

static int qs = 10000;
static int ts = 20000;
static int p = -1;


struct pmf *pdf_load(const char *name, int size)
{
  struct pmf *p;
  FILE *f;

  f = fopen(name, "r");
  if (!f) {
    perror("FOpen");

    return NULL;
  }

  p = pmf_create(size, 0);
  pmf_read(p, f);

  return p;
}

static inline uint64_t get_time(void)
{
  struct timeval tv;
  uint64_t res;

  gettimeofday(&tv, NULL);
  res = tv.tv_sec;
  res = res * 1000000 + tv.tv_usec;

  return res;
}

static int opts_parse(int argc, char *argv[])
{
  int opt;

  while ((opt = getopt(argc, argv, "t:q:p:")) != -1) {
    switch (opt) {
      case 'q':
        qs = atoi(optarg);
        break;
      case 't':
        ts = atoi(optarg);
        break;
      case 'p':
        p = atoi(optarg);
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
  struct pmf *c, *t;
  int i;
  double *m, *dl, sum;
  int opt;
  uint64_t t1, t2;

  opt = opts_parse(argc, argv);
  if (argc - opt < 1) {
    fprintf(stderr, "Usage: %s <c PMF> <t PMF> Q T\n", argv[0]);

    return -1;
  }
printf("Loading C %s\n", argv[opt]);
  c = pdf_load(argv[opt], PMF_SIZE);
  if (p <= 0) {
printf("Loading T %s\n", argv[opt + 1]);
    t = pdf_load(argv[opt + 1], PMF_SIZE);
  }
  t1 = get_time();
  if (p > 0) {
//    printf("P: %d\n", p);
    m = matrix_generate_pseudo(c, p, qs, ts);
  } else {
    m = matrix_generate_generic(c, t, qs, ts);
  }
  if (m == NULL) {
    fprintf(stderr, "Cannot generate matrix!\n");

    return -1;
  }

  dl = dl_generate(m, PMF_SIZE, qs, ts);
  t2 = get_time();
  sum = 0;
  for (i = 0; i <= (PMF_SIZE / qs + 1) * ts; i++) {
    sum += dl[i];
    if (dl[i] > 0) {
      printf("P{d = %d} = %f | %f\n", i, dl[i], sum);
    }
  }
  printf("Ctime: %Lu\n", t2 - t1);

  free(dl);

  return 0;
}
