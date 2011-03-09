#include "pmf.h"
#include "pmf-sample.h"

int pmf_collect(struct pmf *s, int i)
{
  return pmf_set(s, i, pmf_get(s, i) + 1);
}

void pmf_normalise(struct pmf *s)
{
  int i;
  double sum;

  sum = 0;
  for (i = pmf_min(s); i <= pmf_max(s); i++) {
    sum += pmf_get(s, i);
  }
  for (i = pmf_min(s); i <= pmf_max(s); i++) {
    pmf_set(s, i, pmf_get(s, i) / sum);
  }
}
