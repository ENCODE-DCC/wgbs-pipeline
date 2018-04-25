#ifndef _LOKI_H_
#define _LOKI_H_

/****************************************************************************
 *                                                                          *
 *     Loki - Programs for genetic analysis of complex traits using MCMC    *
 *                                                                          *
 *             Simon Heath - University of Washington                       *
 *                                                                          *
 *                       July 1997                                          *
 *                                                                          *
 * loki.h:                                                                  *
 *                                                                          *
 ****************************************************************************/

#include "loki_struct.h"
#include "lk_long.h"

#define LUMPED_ALLELE "__LUMPED__"

struct gtype_data {
  struct gtype_data *next;
  int *gdata;
};

struct Id_Record {
  struct Id_Record **kids;
  double **res; /* Residuals */
  double **pseudo_qt; /* Pseudo quantitative trait (for censored or categorical models) */
  double **vv; /* Weights for student t-model */
  double *bv; /* Breeding values */
  double *bvsum,*bvsum2; /* Stores for average and variance of BV estimates */
  struct id_data *data,**data1;
  struct nucfam *family;
  struct gtype_data gdata;
  double tp[4],tpp[2][2];
  int sire,dam;
  int idx;
  int proband;
  int n_fam;
  int comp; /* Component */
  int fam_code; /* Original family designation */
  int population; /* Original population if specified */
  int rfp;
  int sex,affected,group;
  int nkids;
  int flag;
  int n_rec;
  int *gt_idx; /* Index for original genotypes for this individual in marker->orig_gt */
  int n_gt_sets; /* Number of sets of genotypes for this individual */  
  int allele[2];
};

struct population {
  struct population *next;
  char *name;
  char *ethnicity;
  int idx;
};

struct output_gen {
  struct output_gen *next;
  char *file;
  int link_group;
};

union arg_type {
  char *string;
  int value;
};

struct Id_Recode {
  union arg_type *recode;
  int flag;
};

struct IBD_List {
  double *pos;
  int idx,size;
};

struct Locus {
  double pos[2];
  double phys_pos;
  int *gt; /* genotypes */
  int *seg[2],*seg_bk[2]; /* segregation pattern */
  int *genes[2]; /* Founder genes */
  int *pruned_flag;
  int *founder_flag;
  char *name;
  char *chrom_seg; /* Chromosome segment */
  unsigned long model_flag; /* Which models are affected by this locus? */
  double **eff;
  double *lk_store;
  double *variance;
  double **freq;
  double *aff_freq,*diff_freq;
  int n_alleles; /* In entire pedigree */
  int link_group;
  int flag;
  int index;
  int type;
};

typedef struct {
  int mat;
  int pat;
} gtype;

typedef struct {
  int *gt;
  int *order;
} ord_gtype;

typedef struct {
  ord_gtype mat;
  ord_gtype pat;
} mgtype;

#define MARKER_MICROSAT 1
#define MARKER_SNP 2
#define MARKER_SUPER 3

struct Marker {
  struct Locus locus;
  int marker_type;
  char *name; /* 'Common' name */
  char *rs_id; /* rs number */
  char *cng_id; /* CNG id */ 
  char *mid; /* CNG database id */ 
  int **group;
  double **counts;
  struct Model_Term **mterm;
  int *haplo;
  int *count_flag;
  char **recode;
  signed char **freq_set;
  int *orig_gt; /* original genotypes */
  int *m_flag;
  int *n_all1;   /* In each component */
  int *nhaps[2];
  int *ngens;
  int **allele_trans;
 /* For super-locus */  lk_ulong **all_set;
  lk_ulong *req_set[2];
  lk_ulong *temp[2];
  gtype **gtypes;
  int pos_set;
  
  /* For super-locus */
  struct Marker *parent; /* For members of super-loci */
  int n_children;
  int *child_flag; /* Whether child locus should be included in analysis */
  struct Marker **children;
};

struct Link {
  char *name;
  int *mk_index;
  struct IBD_List *ibd_list;
  struct Link *real_chr; /* In case this is a pseudochromosome */
  double r1[2],r2[2]; /* Map range */
  int ibd_est_type;
  int n_markers;
  int type;
  int sample_pos;
  int range_set[2];
};

struct Variable {
  char *name;
  union arg_type *recode;
  int type;
  int n_levels,index,rec_flag;
};

struct Model_Var {
  int type;
  int var_index;
};

struct Model_Term {
  double *eff;
  struct Model_Var *vars;
  int n_vars;
  int df;
  int out_flag;
};

struct Model {
  struct Model_Term *term;
  struct Model_Var var;
  int n_terms;
  int polygenic_flag;
};

struct lk_ped {
  int ped_size;
  int n_comp;
  int n_genetic_groups;
  int family_id;
  int n_orig_families;
  int n_pops;
  int *singleton_flag;
  int *comp_size;
  int *comp_ngenes;
  int *comp_start;
  struct population *populations;
  struct Id_Record *id_array;
  struct Id_Recode id_recode;
  struct Id_Recode fam_recode;
  struct Id_Recode group_recode;
};

struct lk_markers {
  int n_markers;
  int n_links;
  int sex_map;
  int n_id_recs;
  double total_maplength[2];
  struct Marker *marker;
  struct Link *linkage;
};

struct lk_model {
  int n_models;
  int no_overdominant;
  int censored_flag;
  int censor_mode;
  int use_student_t;
  int tau_mode;
  int tloci_mean_set;
  int n_random;
  int *founder_flag;
  int *pruned_flag;
  char *qtl_name;
  char *tl_name;
  unsigned long *rand_flag;
  struct Variable **rand_list;
  double res_nu; /* Student t parameter */
  double **c_var;
  double tloci_mean;
  double *residual_var;
  double *additive_var;
  double *residual_var_limit;
  double *additive_var_limit;
  int *res_var_set;
  int *add_var_set;
  int *grand_mean_set;
  double mjrgene_var;
  double *grand_mean;
  double *tau;
  double *tau_beta;
  struct Model *models;
  struct SparseMatRec **AIMatrix;
  struct Locus *tlocus;
};

struct lk_data {
  int n_id_records;
  int n_nonid_records;
  struct Variable *id_variable;
  struct Variable *nonid_variable;
};

#define DEFAULT_MAX_TLOCI 16
#define DEFAULT_N_TLOCI 8
#define PRUNE_LOCUS_SPECIFIC 1
#define NO_PRUNE 0
#define DEFAULT_PRUNE_OPTION PRUNE_LOCUS_SPECIFIC

void ReadBinFiles(char **,int,struct loki *);
int ReadXMLFile(struct loki *);
int ReadParam(FILE *,char *,struct loki *);
void Calculate_NRM(struct loki *);
void AllocMarkerStruct(struct loki *);
void AllocEffects(struct loki *);
void LokiSetup(struct loki *);
void FreeStuff(void);
void InitValues(struct loki *);
void SampleLoop(int,int,struct loki *);
void delete_traitlocus(struct Locus *);
int get_new_traitlocus(const int,struct loki *);
void delete_all_traitloci(struct loki *loki);
void loki_identity(double *,int,int,const struct Id_Record *);
char *fget_line(FILE *);
void init_traitlocus(struct loki *,int);
void init_trait_loci(struct loki *);

#endif
