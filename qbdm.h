#include "pmf.h"
#include "cdf.h"
#include "meschac/matrix.h"

double matrix_prob(int i, int j, int n, int q, struct pmf *p);
double matrix_prob2(int i, int j, int q, struct pmf *p, struct pmf *u);
double matrix_prob3(int i, int j, int q, int t, struct pmf *p, struct pmf *u);
void compute_matrixes(MAT *mat,int dim,MAT *B,MAT *A0,MAT *A1,MAT *A2);
MAT *computeR(MAT *mat, MAT *A0, MAT *A1,MAT *A2, double epsilon);
VEC *computeX0(MAT *R,MAT *B0, MAT *A2, VEC *X0);
