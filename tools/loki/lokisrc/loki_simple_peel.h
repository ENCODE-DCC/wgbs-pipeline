#ifndef _LOKI_SIMPLE_PEEL_H_
#define _LOKI_SIMPLE_PEEL_H_

double loki_simple_sample(const struct Simple_Element *,const int,pen_func,lk_ulong **,double **,struct R_Func *,struct loki *);
double loki_simple_peelop(const struct Simple_Element *,const int,const int,pen_func,lk_ulong **,double **,struct R_Func *,struct loki *);
double loki_simple_peelop_x(const struct Simple_Element *,const int,const int,pen_func,lk_ulong **,double **,struct R_Func *,struct loki *);
double peel_to_par(const struct Simple_Element *,const int,pen_func,lk_ulong **,struct R_Func *,struct loki *);

#endif
