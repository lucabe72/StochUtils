#include <stdio.h>
#include <stdlib.h>

#include "pmf.h"
#include "z.h"

struct pmf *z_generate(const struct pmf *t, int ts)
{
  double sum;
  unsigned int j, th;
  struct pmf *z;

  z = pmf_create(pmf_max(t) / ts * ts, 0);

  th = ts;
  sum = 0;
  for (j = pmf_min(t); j < pmf_max(t); j++) {
    if (j < ts) {
      fprintf(stderr, "Error!!! Interarrival times MUST be greater than %d!!!\n",
		      ts);
      return NULL;
    }

    if (j >= th) {
      pmf_set(z, th - ts, sum);
      th += ts;
      sum = 0;
    }
    sum += pmf_get(t, j);
  }

  return z;
}

