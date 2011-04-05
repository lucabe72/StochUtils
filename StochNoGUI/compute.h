#ifndef __COMPUTE_H__
#define __COMPUTE_H__

double *create_compute_window(double *mat, int n, int iterations);
double *pseudo_generate(struct distribution *d, int q, int n);
double *generic_generate(struct distribution *d_c, double *t, int q, int n);
double *generic_transform(struct distribution *d_t, int t, int n);
double *solve(double *mat, int n, int iter);

double *matrix_generate_generic(struct distribution *exec, struct distribution *interarrival, int qs, int ts);
double *matrix_generate_pseudo(struct distribution *exec, int period, int qs, int ts);
double *dl_generate(double *matrix, int size, int qs, int ts);

#endif
