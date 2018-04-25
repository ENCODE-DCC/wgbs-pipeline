#include <config.h>
#include <math.h>
#include "mat_utils.h"

/*
 * Compute Cholesky factorization of symmetric matrix. 
 * Inputs:  a - Lower triangle of matrix 
 *          n - Size of matrix 
 *
 * Outputs: b - Cholesky factor of a
 *        det - log determinant of a 
 *
 * Returns 1 if a is not positive definite, otherwise returns 0 
 *
 * Note that b can point to the same location as a, otherwise a is unmodified 
 */
int chol_fact(double *a, double *b, int n, double *ldet)
{
	int i, j, k, kk;
	double zz, piv;

	kk=n*(n+1)/2;
	*ldet=0.0;
	for(i=n-1;i>=0;i--) {
		zz=a[--kk];
		for(k=i+1;k<n;k++) zz-=BB(b,k,i)*BB(b,k,i);
		if(zz<=0.0) return -1;
		b[kk]=piv=sqrt(zz);
		*ldet+=log(zz);
		for(j=i-1;j>=0;j--) {
			zz=a[--kk];
			for(k=i+1;k<n;k++) zz-=BB(b,k,i)*BB(b,k,j);
			b[kk]=zz/piv;
		}
	}
	return 0;
}

/*
 * Solve A*y=z 
 *
 * Inputs:  a - Cholesky factor of A computed by chol_fact()
 *          y - Right hand side vector 
 *          n - Size of A and y 
 *
 * Outputs: z - Solution of A*y=z 
 *
 * a is not modified.  The same is true for y unless y and z point to the same
 * location, which is OK 
 */
void chol_solve(double *a, double *y, double *z, int n)
{
	int i, j, kk;
	double zz;

	for(i=n-1;i>=0;i--) {
		zz=y[i];
		for(j=i+1;j<n;j++) zz-=BB(a,j,i)*z[j];
		z[i]=zz/BB(a,i,i);
	}
	kk=0;
	for(i=0;i<n;i++) {
		zz=z[i];
		for(j=0;j<i;j++) zz-=a[kk++]*z[j];
		z[i]=zz/a[kk++];
	}
}

/*
 * Compute inverse of symmetrix pd matrix A 
 *
 * Inputs:  a - Cholesky factor of A computed by chol_fact() n - Size of A
 *
 * Outputs: b - (half stored) inverse of A 
 *
 * Note that b can point to the same location as a, otherwise a is unmodified 
 */
void chol_inv(double *a, double *b, int n)
{
    int i, j, k, kk;
    double zz;

    kk=n*(n+1)/2;
    for(i=n-1;i>=0;i--) {
		 zz=BB(a,i,i);
		 b[--kk]=1.0/zz;
		 for(j=i-1;j>=0;j--)	{
			 zz=0.0;
			 for(k=i;k>j;k--) zz-=BB(a,k,j)*BB(b,i,k);
			 b[--kk]=zz/BB(a,j,j);
		 }
    }
    kk=n*(n+1)/2;
    for(i=n-1;i>=0;i--) {
		 for(j=i;j>=0;j--) {
	    zz=0.0;
	    for(k=j;k>=0;k--) zz+=BB(b,i,k)*BB(b,j,k);
	    b[--kk]=zz;
		 }
    }
}

