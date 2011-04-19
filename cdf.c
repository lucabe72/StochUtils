#include "pmf.h"
#include "cdf.h"

struct pmf *pmf2cdf(struct pmf *p)
{
  struct pmf *res;
  double sum = 0;
  int i;
  
  res = pmf_create(p->size, 0);
  for (i = pmf_min(p); i <= pmf_max(p); i++) {
    sum += pmf_get(p, i);
    pmf_set(res, i, sum);
  }

  return res;
}

struct pmf *cdf2pmf(struct pmf *c)
{
  struct pmf *res;
  double sum = 0;
  int i;
  
  res = pmf_create(c->size, 0);
  for (i = pmf_min(c); i <= pmf_max(c); i++) {
    pmf_set(res, i, pmf_get(c, i) - sum);
    sum += pmf_get(c, i);
  }

  return res;
}
