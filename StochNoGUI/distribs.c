#include <stdlib.h>

#include "distr.h"

#define MAX 2000

void init_distribs(struct sys_distr *sd)
{
  sd->first = NULL;
  sd->n_distr = 0;
  sd->max_distr_dim = MAX;
}
