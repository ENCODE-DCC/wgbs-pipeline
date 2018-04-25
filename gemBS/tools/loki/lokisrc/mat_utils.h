#ifndef _MAT_UTILS_H_
#define _MAT_UTILS_H_

#define BB(b,i,j) ((b)[i*(i+1)/2+j])
#define BB1(b,i,j) ((i<j)?BB(b,j,i):BB(b,i,j))
#define IDX(i,j) (i*(i+1)/2+j)
int chol_fact(double *, double *, int , double *);
void chol_solve(double *, double *, double *, int);
void chol_inv(double *, double *, int);

#endif
