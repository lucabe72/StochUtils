#include "qbdm.h"
#include "cdf.h"
#include "pmf.h"
#include <malloc.h>
#include <math.h>
#include "meschac/matrix2.h"


double matrix_prob(int i, int j, int n, int q, struct pmf *p)
{
    if (j == 0) {
	if (pmf_max(p) < ((n - i) * q))
	    return 1;
	if (n - i > 0)
	    return cdf_get(p, (n - i) * q);
	else
	    return 0;
    } else if (i == j)
	return cdf_get(p, n * q) - cdf_get(p, (n - 1) * q);

    else
	return cdf_get(p, (n + j - i) * q) - cdf_get(p,
						     (n + j - i - 1) * q);
}

double matrix_prob2(int i, int j, int q, struct pmf *p, struct pmf *u)
{
    int n=0,k=0;
    double prob=0;
//    printf("(%i,%i),%i,%i=",i,j,q,pmf_max(u));
    if (j == 0) {
	/*if (pmf_max(p) < ((pmf_max(u) - i) * q))
	    return 1;*/
	//if (n - i > 0){
	for (k = 0; k < u->size; k++) {
	  n=(k - u->offset);
//	printf("size: %i, k=%i,off=%i n=%i\n",u->size,k,u->offset,n);
	  prob +=  u->elems[k]*cdf_get(p, (n - i) * q);
	  }

	    /*for (n=i; n<pmf_max(u); n++)
	      prob+=pmf_get(u,n)*cdf_get(p, (n - i) * q);*/
//	    printf("%f\n",prob);
	    return prob;
	/*    }
	else
	    return 0;*/
    } else if (i == j){
	for (k = 0; k < u->size; k++) {
	  n=(k - u->offset);
//	printf("size: %i, k=%i,off=%i n=%i\n",u->size,k,u->offset,n);
	  prob +=  u->elems[k]*(cdf_get(p, n * q) - cdf_get(p, (n - 1) * q));
	}

        /*for (n=i; n<pmf_max(u); n++)
	   prob+=pmf_get(u,n)*(cdf_get(p, n * q) - cdf_get(p, (n - 1) * q));*/
//	    printf("%f\n",prob);
	return prob;
	}
    else
    {
    for (k = 0; k < u->size; k++) {
    	n=(k - u->offset);
//	printf("size: %i, k=%i,off=%i n=%i\n",u->size,k,u->offset,n);
        prob +=  u->elems[k]*(cdf_get(p, (n + j - i) * q) - cdf_get(p,(n + j - i - 1) * q)) ;
	  }

    	/*for (n=i; n<pmf_max(u); n++)
	   prob+=pmf_get(u,n)*(cdf_get(p, (n + j - i) * q) - cdf_get(p,(n + j - i - 1) * q));
    	for (n=i; n<pmf_max(u); n++)
	   prob+=pmf_get(u,n)*(cdf_get(p, (n + j - i) * q) - cdf_get(p,(n + j - i - 1) * q));*/
//	    printf("%f\n",prob);
    	return prob;
    }
}




void compute_matrixes(MAT * mat, int dim, MAT * B, MAT * A0,
		      MAT * A1, MAT * A2)
{
    int i, j, h, k;
    for (i = 0; i < dim; i++) {
	for (j = 0; j < dim; j++)
	    m_set_val(B, i, j, m_get_val(mat, i, j));
	for (j = dim, h = 0; j < dim * 2; j++, h++)
	    m_set_val(A0, i, h, m_get_val(mat, i, j));
    }
    for (i = dim, k = 0; i < dim * 2; i++, k++) {
	for (j = 0; j < dim; j++)
	    m_set_val(A2, k, j, m_get_val(mat, i, j));
	for (j = dim, h = 0; j < dim * 2; j++, h++)
	    m_set_val(A1, k, h, m_get_val(mat, i, j));
    }
}



MAT *computeR(MAT * mat, MAT * A0, MAT * A1, MAT * A2, double epsilon)
{
    int m = mat->m;
    int f = 1;
    double d;
    MAT *RT, *Rnew, *Rdiff, *res1, *res2, *res3, *res4;
    RT = m_get(m, m);
    RT = m_zero(RT);
    res1 = m_get(m, m);
    res2 = m_get(m, m);
    res3 = m_get(m, m);
    res4 = m_get(m, m);
    Rnew = m_get(m, m);
    Rdiff = m_get(m, m);
    while (f == 1) {
	res1 = m_mlt(RT, A1, res1);
	res2 = m_mlt(RT, RT, res2);
	res3 = m_add(A0, res1, res3);
	res4 = m_mlt(res2, A2, res4);
	Rnew = m_add(res3, res4, Rnew);
	Rdiff = m_sub(RT, Rnew, Rdiff);
	d = m_norm_inf(Rdiff);
	if (d < epsilon)
	    f = 0;
	else
	    m_move(Rnew, 0, 0, Rnew->m, Rnew->n, RT, 0, 0);
    }
    m_free(res1);
    m_free(res2);
    m_free(res3);
    m_free(res4);
    m_free(Rnew);
    m_free(Rdiff);
    return RT;
}

MAT *pinv(MAT * M, MAT * out)
{
    MAT *U, *Ut, *V, *S, *res1, *res2;
    VEC *b;
    unsigned int i = 0;
    U = m_get(M->m, M->m);
    Ut = m_get(M->m, M->m);
    V = m_get(M->n, M->n);
    S = m_get(M->m, M->n);
    S = m_zero(S);
    b = v_get(M->n);
    res1 = m_get(M->m, M->n);
    res2 = m_get(M->m, M->n);
    svd(M, U, V, b);
    for (i = 0; i < b->dim; i++)
	m_set_val(S, i, i, 1 / v_get_val(b, i));
    Ut = m_transp(U, Ut);
    res1 = m_mlt(Ut, S, res1);
    res2 = m_mlt(res1, V, res2);
    out = m_transp(res2, out);
    m_free(U);
    m_free(Ut);
    m_free(V);
    m_free(S);
    m_free(res1);
    m_free(res2);
    v_free(b);
    return out;
}


VEC *computeX0(MAT * R, MAT * B0, MAT * A2, VEC * v)
{
    int m = R->m;
    MAT *eye, *res, *res1, *res2, *res3, *res4, *res5, *M, *pi, *ones,
	*zero;
    res = m_get(m, m);
    res1 = m_get(m, m);
    res2 = m_get(m, m);
    res3 = m_get(m, m);
    ones = m_get(m, 1);
    ones = m_ones(ones);
    eye = m_get(m, m);
    eye = m_ident(eye);
    zero = m_get(1, m);
    zero = m_zero(zero);
    res4 = m_get(1, m + 1);
    res5 = m_get(1, m);
    pi = m_get(m, m);
    M = m_get(m, m + 1);
    // first part of M
    res = m_add(B0, m_sub(m_mlt(R, A2, res), eye, res), res);
    m_move(res, 0, 0, res->m, res->m, M, 0, 0);

    //second part of M
    res1 = m_sub(eye, R, res1);
    res2 = m_inverse(res1, res2);
    res3 = m_mlt(res2, ones, res3);

    m_move(res3, 0, 0, res3->m, res3->n, M, 0, m);

    res4 = m_zero(res4);
    m_set_val(res4, 0, m, 1);
    pi = pinv(M, pi);
    //togheter
    res5 = m_mlt(res4, pi, res5);
    //compute x0
    v = mv_move(res5, 0, 0, res5->m, res5->n, v, 0);
    m_free(eye);
    m_free(res);
    m_free(res1);
    m_free(res2);
    m_free(res3);
    m_free(res4);
    m_free(res5);
    m_free(M);
    m_free(pi);
    m_free(ones);
    m_free(zero);
    return v;
}
