#ifndef _PED_UTILS_H_
#define _PED_UTILS_H_

#ifndef _LOKI_STRUCT_H_
#include <loki_struct.h>
#endif

typedef struct nucfam {
  struct Id_Record *father;
  struct Id_Record *mother;
  struct Id_Record **kids;
  struct f_gtype *gtypes;
  struct f_gtype **err_gtypes;
  struct gen_err_rec *gen_errs;
  int *status;
  int n_gen_errs;
  int flag;
  int n_gt;
  int comp;
  int n_err;
  int n_err1[N_LINK_TYPES];
} nuc_fam;

struct lk_fam {
  int *cp_start;
  nuc_fam *families;	
};

void prune_ped_for_locus(struct Locus *,struct loki *);
void prune_ped(struct loki *);
void build_families(struct loki *);
int check_ped(struct loki *);
void get_loops(struct loki *);

#endif
