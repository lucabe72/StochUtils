#include "pmf.h"
#include "pmf-modify.h"

struct pmf *pmf_restrict(const struct pmf *p, int max)
{
  unsigned int i;
  struct pmf *new;
  double sum;

  new = pmf_create(max, p->offset);
  sum= 0;

  for (i = 0; i < max + p->offset; i++) {
    pmf_set(new, i, pmf_get(p, i));
    sum += pmf_get(p, i);
  }
  if (1.0 - sum > 1e-10) {
    new->tail = 1.0 - sum;
  }

  return new;
}
