#ifndef _KINSHIP_H_
#define _KINSHIP_H_

struct kin {
  double p[9];
  double expected,sd;
  int n;
};

double kinship(const int, const int,const struct Id_Record *);
void calc_kin(int *,struct loki *);
double calc_gen(int,int *,int *,int);

#endif
