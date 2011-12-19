#include "pmf.h"
#include "dl.h"

struct pmf *stochdl_compute(const struct pmf *v, unsigned int q, unsigned int t)
{
  struct pmf *dl;
  unsigned int maxdl;
  unsigned int i;

  maxdl = (pmf_max(v) + 1) / q * t + t;
  dl = pmf_create(maxdl, 0);
  for (i = 0; i < maxdl; i += t) {
    pmf_set(dl, i, pmf_get(v, (i / t) * q));
  }
//  dl->tail = pmf_get(v, pmf_max(v) / q * t);
  return dl;
}
