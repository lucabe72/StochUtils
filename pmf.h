#ifndef PMF_H
#define PMF_H

struct pmf {
  double *elems;
  unsigned int size;
  unsigned int offset;
  int min;
  int max;
  double tail;
  //unsigned long long int samples;
};

double pmf_avg(const struct pmf *p);
double pmf_std(const struct pmf *p);
double pmf_var(const struct pmf *p);
struct pmf *pmf_create(unsigned int size, unsigned int offs);
void pmf_free(struct pmf *p);
int pmf_set(struct pmf *d, int val, double p);
int pmf_set_samples(struct pmf *d, int samples);
double pmf_sum(const struct pmf *d, unsigned int i);
int pmf_check(const struct pmf *d);

static inline double pmf_get(const struct pmf *d, int i)
{
  return d->elems[i + d->offset];
}

static inline double pmf_tail(const struct pmf *d)
{
  return d->tail;
}

static inline int pmf_max(const struct pmf *d)
{
  return d->max;
}

static inline int pmf_min(const struct pmf *d)
{
  return d->min;
}
#endif
