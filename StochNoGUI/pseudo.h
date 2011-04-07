#ifndef PSEUDO_H
#define PSEUDO_H

double *pseudo_generate(struct pmf *d, int q);
double *solve(double *mat, int n, int iter);

#endif	/* PSEUDO_H */
