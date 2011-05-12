#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#ifdef GSL
#include <stdarg.h>
#include <gsl/gsl_randist.h>
#endif

#include "pmf.h"
#include "pmf-sample.h"
#include "pmf-modify.h"

#define RANDOM		0
#define UNIFORM		1
#define GAUSSIAN	2
#define GAMMA		3
#define BETA		4
#define EXPONENTIAL	5
#define LAPLACE		6
#define EXPPOWER	7
#define CAUCHY		8
#define RAYLEIGH	9
#define LANDAU		10
#define LOGNORMAL	11

#define MIN_DIS (-100)
#define MAX_DIS 100

static int offs = 100;
static int size = 500;
static int samples = 10000;
static int verbose;
static int seed;
static int distro = 0;
static double par1 = 0;
static double par2 = 0;

void help(void)
{
  printf("Usage:\n");
  printf("./pmfgen <opts>\n\n");
  printf("Common options: \n");
  printf("-o offset\t Minimum value of the pmf \n");
  printf("-s size  \t Size the pmf \n");
  printf("-S samples\t Number of point of pmf \n");

  printf("\nDistro specific: \n");
  printf("-d 0\t Random pmf \n");
  printf("-d 1\t Uniform pmf \n");
  printf("-d 2 -p sigma\t Gaussian pmf \n");
  printf("-d 3 -p alpha -P beta\t Gamma pmf \n");
  printf("-d 4 -p alpha -P beta\t Beta pmf \n");
  printf("-d 5 -p nu\t Exponential pmf \n");
  printf("-d 6 -p alpha\t Laplace pmf \n");
  printf("-d 7 -p alpha -P beta\t ExpPower pmf \n");
  printf("-d 8 -p alpha\t Cauchy pmf \n");
  printf("-d 9 -p sigma\t Rayleigh pmf \n");
  printf("-d 10 \t Landau pmf \n");
  printf("-d 11 -p zeta -P sigma\t LogNormal pmf \n");


}

#ifdef GSL
static double gen_number(int value, int distro, va_list ap)
{
  double val = 0.0;
  switch (distro) {
  case RANDOM:{
      val = rand() * 1.0 / RAND_MAX;
      break;
    }
  case UNIFORM:{
      double max = va_arg(ap, double);
      double min = va_arg(ap, double);
      double x = (value - MIN_DIS) * 1.0;
      val = gsl_ran_flat_pdf(x, 0, max - min);
      break;
    }
  case GAUSSIAN:{
      double sigma = va_arg(ap, double);
      val = gsl_ran_gaussian_pdf(value, sigma);
      break;
    }
  case GAMMA:{
      double alpha = va_arg(ap, double);
      double beta = va_arg(ap, double);
      double x = (value - MIN_DIS) * 1.0;
      val = gsl_ran_gamma_pdf(x, alpha, beta);
      fprintf(stderr, "%f %f %f %f\n", x, alpha, beta, val);
      break;
    }
  case BETA:{
      double alpha = va_arg(ap, double);
      double beta = va_arg(ap, double);
      double x = (value - MIN_DIS) * 1.0 / (MAX_DIS - MIN_DIS);
      val = gsl_ran_beta_pdf(x, alpha, beta);
      break;
    }
  case EXPONENTIAL:{
      double mu = va_arg(ap, double);
      double x = (value - MIN_DIS) * 1.0;
      val = gsl_ran_exponential_pdf(x, mu);
    }
  case LAPLACE:{
      double alpha = va_arg(ap, double);
      val = gsl_ran_laplace_pdf(value, alpha);
      break;
    }
  case EXPPOWER:{
      double alpha = va_arg(ap, double);
      double beta = va_arg(ap, double);
      double x = (value - MIN_DIS) * 1.0;
      val = gsl_ran_exppow_pdf(x, alpha, beta);
      break;
    }
  case CAUCHY:{
      double alpha = va_arg(ap, double);
      val = gsl_ran_cauchy_pdf(value, alpha);
      break;
    }
  case RAYLEIGH:{
      double sigma = va_arg(ap, double);
      double x = (value - MIN_DIS + 1) * 1.0;
      val = gsl_ran_rayleigh_pdf(x, sigma);
      break;
    }
  case LANDAU:{
      val = gsl_ran_landau_pdf(value);
      break;
    }
  case LOGNORMAL:{
      double zeta = va_arg(ap, double);
      double sigma = va_arg(ap, double);
      double x = (value - MIN_DIS + 1) * 1.0;
      val = gsl_ran_lognormal_pdf(x, zeta, sigma);
      break;
    }
  default:
    {
//fprintf(stderr,"Distro not yet implemented\n");
      help();
    }

  }
  va_end(ap);
  return val;
}

static void pmf_gen(struct pmf *p, int offs, int distro, ...)
{
  int i;
  double val;
  va_list ap;
  va_start(ap, distro);

  for (i = MIN_DIS; i < MAX_DIS; i++) {
    val = gen_number(i, distro, ap);
    pmf_set(p, i + MAX_DIS, val);
  }
  pmf_shift(p, offs);
  pmf_normalise(p);
  va_end(ap);
}

#else

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

#endif

static int opts_parse(int argc, char *argv[])
{
  int opt;

  while ((opt = getopt(argc, argv, "d:p:P:o:s:S:r:vh")) != -1) {
    switch (opt) {
      case 'd':
        distro = atoi(optarg);
        break;
      case 'o':
        offs = atoi(optarg);
        break;
      case 's':
        size = atoi(optarg);
        break;
      case 'S':
        samples = atoi(optarg);
        break;
      case 'r':
        seed = atoi(optarg);
        break;
      case 'p':
        par1 = atof(optarg);
        break;
      case 'P':
        par2 = atof(optarg);
        break;
      case 'v':
        verbose = 1;
        break;
      case 'h':
        help();
        exit(EXIT_FAILURE);
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
  if (seed) {
    srand(seed);
  }
  p = pmf_create(size + offs, 0);
#ifdef GSL
  pmf_gen(p, offs, distro, par1, par2, NULL);
#else
  pmf_rnd(p, offs, size, samples);
#endif

  if (pmf_check(p)) {
    printf("Bad PMF: %d | %f!!!\n", pmf_check(p), pmf_sum(p, 0));

    return -1;
  }
  for (i = 0; i <= pmf_max(p); i++) {
    if (verbose) {
      printf("P{x = %d} = %f\n", i, pmf_get(p, i));
    } else {
      printf("%d\t%14.15f\n", i, pmf_get(p, i));
    }
  }
  if (verbose) {
    printf("P{x >= %d} = %f\n", pmf_max(p), pmf_tail(p));
  }

  return 0;
}
