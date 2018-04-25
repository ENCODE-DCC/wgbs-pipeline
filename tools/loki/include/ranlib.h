#ifndef _RANLIB_H_
#define _RANLIB_H_

/* Prototypes for all user accessible RANLIB routines */

/* #define ranf genrand*/
#define ranf genrand
#define safe_ranf safe_genrand
#define init_ranf sgenrand

extern double genrand(void);
extern unsigned int genint(void);
extern double safe_genrand(void);
extern double test_genrand(void);
extern int set_mt_idx(int);
extern void sgenrand(unsigned int sd);
double ppnd(double,int *);
double trunc_normal(const double,const double,const int,int *);
void dirichlet(const int,const double *,double *);
extern void advnst(int k);
extern double genbet(const double,const double);
extern double genchi(const double);
extern double genexp(double const);
extern double genf(const double,const double);
extern double gengam(const double,const double);
extern void genmul(const int,const double *,const int,int *);
extern double gennch(const double,const double);
extern double gennf(const double,const double,const double);
extern double gennor(const double,const double);
extern void genprm(int *,const int);
extern double genunf(const double,const double);
extern int ignbin(const int,const double);
extern int ignnbn(const int,const double);
extern int ignpoi(const double);
extern int ignuin(const int,const int);
extern int mltmod(const int,const int,const int);
extern double sgamma(const double);
extern double sexpo(void);
extern double snorm(void);

#endif
