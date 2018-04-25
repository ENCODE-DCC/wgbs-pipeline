#ifndef _GET_PAR_PROBS_H_
#define _GET_PAR_PROBS_H_

double get_par_probs(double *,const int,struct Marker *,pen_func,lk_ulong **,double **,
							struct R_Func *,struct loki *);
double get_par_probs_x(double *,const int,struct Marker *,pen_func,lk_ulong **,double **,
							struct R_Func *,struct loki *);
double get_trait_par_probs(double *,const int,const int,trait_pen_func *,double **,
									struct R_Func *,struct loki *);



#endif
