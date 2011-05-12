#ifndef CDF_H
#define CDF_H

struct pmf *pmf2cdf(struct pmf *p);
struct pmf *cdf2pmf(struct pmf *c);

static inline double cdf_get(const struct pmf *d, int i)
{
  if (i<0) return 0.0;
  if (i>=d->max) return 1.0;
  return d->elems[i + d->offset];
}
#endif
