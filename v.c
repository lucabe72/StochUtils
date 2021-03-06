#include <math.h>
#include "pmf.h"
#include "v.h"

struct pmf *v_compute(const struct pmf *w, const struct pmf *c, unsigned int max, unsigned int step)
{
  unsigned int i, vi;
  struct pmf *v;

  v = pmf_create(max, 0);

  for (vi = 0; vi < max; vi += step) {
    double p;
    unsigned int cmax;
    unsigned int min;

    p = 0;
    cmax = pmf_max(c);
    min = cmax < vi ? cmax : vi;
//    if (cmax > vi) {
//      cmax = vi;
//    }
#if 0
    for (i = 0; i <= cmax; i++) {
      if (i <= vi) {
        p += pmf_get(c, i) * (1.0 - pmf_get(w, (int)vi - (int)i));
      } else {
        p += pmf_get(c, i);
      }
    }
#else
    for (i = 0; i <= min; i++) {
      p += pmf_get(c, i) * (1.0 - pmf_get(w, (int)vi - (int)i));
    }
    for (i = vi + 1; i <= cmax; i++) {
      p += pmf_get(c, i);
    }
#endif
    p += pmf_tail(c);
    pmf_set(v, vi, 1.0 - p);
  }
  if (1.0 - pmf_get(v, pmf_max(v)) > 1e-10) {
    v->tail = 1.0 - pmf_get(v, pmf_max(v));
  }

  return v;
}
