#include "models.h"
#include "cdf.h"

double prob_efta(int i, int j, int q, struct pmf *p, struct pmf *u)
{
    int z = 0;
    unsigned int h = 0;
    //int k = pmf_min(u);
    double prob = 0, pv;
    for (h = 0; h < u->size; h++) {
	z = (h - u->offset);
	pv = u->elems[h];
	if (pv > 0) {
	    if (i + 1 - z > 0)
		prob +=
		    pv * (cdf_get(p, (j - i + z) * q) -
			  cdf_get(p, (j - i + z - 1) * q));
	    else
		prob +=
		    pv * (cdf_get(p, (j + 1) * q) - cdf_get(p, (j) * q));
	}
    }
    return prob;
}

double prob_rtss(int h, int k, int q, int z, struct pmf *u)
{
    double prob = 0.0;
    int maxc = pmf_max(u);
    int minc = pmf_min(u);

    if (h - z * q <= 0) {
	prob = pmf_get(u, k);
	return prob;
    }
    if ((k - h + z * q >= minc) && (k - h + z * q <= maxc))
	prob = pmf_get(u, k - h + z * q);
    else
	prob = 0.0;
    return prob;
}

double prob_last(int i, int j, int q, int n, int d, struct pmf *c)
{
    double prob = 0, first, second;
    if (n * q + (j - i) * d > 0)
	first = cdf_get(c, n * q + (j - i) * d);
    else
	first = 0.0;
    if (j == 0)
	prob = first;
    else {
	if (n * q + (j - i - 1) >= 0)
	    second = cdf_get(c, n * q + (j - i - 1) * d);
	else
	    second = 0;
	if (first != 0)
	    prob = first - second;
	else
	    prob = 0;
    }
    return prob;

}
