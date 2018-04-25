#ifndef _GEN_ELIM_H_
#define _GEN_ELIM_H_

#include "loki.h"

typedef struct f_gtype {
	gtype par[2];
} fam_gtype;

typedef struct mf_gtype {
	ord_gtype gt[4];
} fam_mgtype;

typedef struct g_err {
	struct g_err *next;
	struct Id_Record *id;
	nuc_fam *fam;
	int err_type;
} gen_elim_err;

typedef struct rf_err {
	struct rf_err *next;
	nuc_fam *fam;
	int *ngens;
	gtype **gtypes;
} gen_elim_rf;

struct gen_err_rec {
	int n_errors;	
	struct Id_Record *id_list;
};

#define GEN_ELIM_X_HET_MALE 1 /* Heterozygous male for X-linked data */
#define GEN_ELIM_Y_HET_MALE 2 /* Heterozygous male for Y-linked data */
#define GEN_ELIM_MIT_HET_FEMALE 3 /* Heterozygous female for mitochondrial data */
#define GEN_ELIM_Y_OBS_FEMALE 4 /* Observed female for Y-linked data */
#define GEN_ELIM_MIT_OBS_MALE 5 /* Observed male for mitochondrial data */
#define GEN_ELIM_HALF_OBS 6 /* Half observed genotype */
#define GEN_ELIM_PASS1 7 /* Error in pass 1 (id is a child) */
#define GEN_ELIM_PASS2_KID 8 /* Error in pass 2 (id is a child) */
#define GEN_ELIM_PASS2_PAR 9 /* Error in pass 2 (id is a parent) */
#define GEN_ELIM_PASS2_YM 10 /* Error in pass 2 (Y or Mitochondrial locus) */

int Gen_Elim(struct loki *);
void del_gtypes(struct Id_Record *,gtype **,int *);
int add_ind_gtype(gtype **,int,int *,const int *);
void recode_alleles(int,struct loki *);
gen_elim_err *process_mark(struct Marker *,struct loki *,int);
void process_mark_errors(struct Marker *,struct loki *,int);
gen_elim_err *process_mark_x(struct Marker *,struct loki *,int,int);
gen_elim_err *process_mark_y(struct Marker *,struct loki *,int,int);
gen_elim_err *add_gen_elim_err(gen_elim_err *,nuc_fam *,struct Id_Record *,int);
gen_elim_err *invert_err_list(gen_elim_err *);
void handle_errors(struct loki *,gen_elim_err **);
void sort_fam_error(int *,int,struct loki *);
void sort_ind_error(int *,int,int *);
void set_recode(struct Marker *,struct loki *,int);
gen_elim_err *gen_elim_marker(int,struct loki *); 

#endif

