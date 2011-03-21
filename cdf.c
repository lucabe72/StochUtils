#include "cdf.h"
#include "pmf.h"

struct pmf *cdf_pmf(struct pmf *p)
{
struct pmf *res;
double sum=0;
int i;
res = pmf_create(p->size, 0);
for (i = pmf_min(p); i <= pmf_max(p); i++) {
	sum=sum+pmf_get(p, i);
	pmf_set(res,i,sum);
}
return res;
};
