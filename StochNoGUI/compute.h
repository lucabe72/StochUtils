#ifndef __COMPUTE_H__
#define __COMPUTE_H__

double *pseudo_generate(struct distribution *d, int q, int n);
double *generic_generate(struct distribution *d_c, double *t, int q, int n);
double *generic_transform(struct distribution *d_t, int t, int n);
double *create_compute_window(double *mat, int n, int iterations);
double *solve(double *mat, int n, int iter);

#endif
