#include "pmf.h"

#define ETFAMODEL 1
#define RTSSMODEL 2
#define LASTMODEL 3

double prob_efta(int i, int j, int q, struct pmf *p, struct pmf *u);
double prob_rtss(int h, int k, int q, int z, struct pmf *u);
double prob_last(int h, int k, int q, int z, int d, struct pmf *c);
