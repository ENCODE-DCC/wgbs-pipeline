#ifndef _HANDLE_RES_H_
#define _HANDLE_RES_H_

#define RES_PRIOR_V0 1.0
#define RES_PRIOR_S0 1.0

extern double Calc_Res_Ratio(double,double,struct loki *);
extern double Calc_Resprop(struct loki *);
extern double Calc_CensResLike(struct loki *);
extern double Calc_ResLike(struct loki *);
extern double Sample_ResVar(struct loki *);
extern double Calc_Var_Prior(double,const struct loki *);
extern double Recalc_Res(int,struct loki *);


#endif
