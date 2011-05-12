#include <stdlib.h>

#include "pmf.h"

const double epsilon = 1e-10;

double pmf_avg(const struct pmf *p)
{
  unsigned int i;
  double avg = 0;

  for (i = 0; i < p->size; i++) {
    avg += ((i - p->offset) * p->elems[i]);
  }

  return avg;
}

void pmf_shift(struct pmf *p, int start)
{
  int i;
  if (pmf_min(p) > start) {

  } else {
    for (i = pmf_max(p); i > 0; i--)
      pmf_set(p, i + start - pmf_min(p), pmf_get(p, i));
    for (i = 0; i <= start; i++)
      pmf_set(p, i, 0.0);
  }
}


struct pmf *pmf_create(unsigned int size, unsigned int offs)
{
  struct pmf *res;
  unsigned int i;

  res = malloc(sizeof(struct pmf));

  if (res == NULL) {
    return NULL;
  }
  res->tail = 0.0;
  res->min = size - offs;
  res->max = -offs;
  res->size = size;
  res->offset = offs;
  res->elems = malloc(sizeof(double) * size);
  if (res->elems == NULL) {
    free(res);

    return NULL;
  }
  for (i = 0; i < size; i++) {
    res->elems[i] = 0;
  }

  return res;
}

int pmf_set(struct pmf *d, int val, double p)
{
  if (val > (int)d->size - (int)d->offset) {
    return -1;
  }
  if (p > 1e-10) {
    if (val > d->max) {
      d->max = val;
    }
    if (val < d->min) {
      d->min = val;
    }
  }
  d->elems[val + d->offset] = p;

  return 1;
}

int pmf_set_samples(struct pmf *d, int samples)
{
  unsigned int i;

  if (d->tail > 0) {
    return -1;
  }

  d->tail = 1.0 / samples;		/* epsilon = 1 / samples */
  for (i = 0; i < d->size; i++) {
    if (d->elems[i] > epsilon) {
      d->elems[i] = d->elems[i] * samples / (samples + 1);
    }
  }

  return 1;
}

double pmf_sum(const struct pmf *d, unsigned int i)
{
  double p;

  p = pmf_tail(d);
  for (; i < d->size; i++) {
    p += d->elems[i];
  }

  return p;
}

int pmf_check(const struct pmf *d)
{
  double p;

  p = pmf_sum(d, 0);

  if (p < 1.0 - epsilon) return 1;
  if (p > 1.0 + epsilon) return 2;

  return 0;
}
