double *matrix_generate_generic(struct pmf *exec, struct pmf *interarrival, int qs, int ts);
double *matrix_generate_pseudo(struct pmf *exec, int period, int qs, int ts);
double *dl_generate(double *matrix, int size, int qs, int ts);
