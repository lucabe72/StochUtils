#include <stdio.h>

#include "pmf.h"
#include "pmf-file.h"

int pmf_read(struct pmf *d, FILE *f)
{
  int cnt = 0;

  while (!feof(f)) {
    int res;
    unsigned int i;
    double p;

    res = fscanf(f, "%u\t%lf\n", &i, &p);
    if (res == 2) {
      if (pmf_set(d, i, p) < 0) {
        fprintf(stderr, "Failed to insert %u: PMF(%u)=%f\n", cnt, i, p);

        return -1;
      }
    } else if (res != 0) {
      perror("fscanf");

      return -1;
    }
    cnt++;
  }

  return 1;
}
