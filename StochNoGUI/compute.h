#ifndef __COMPUTE_H__
#define __COMPUTE_H__

double *pseudo_generate(struct pmf *d, int q);
double *generic_generate(struct pmf *d_c, double *t, int q);
double *generic_transform(struct pmf *d_t, int t);
double *create_compute_window(double *mat, int n, int iterations);
double *solve(double *mat, int n, int iter);

#endif
