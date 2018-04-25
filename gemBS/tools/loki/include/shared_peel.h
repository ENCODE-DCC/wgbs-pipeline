#ifndef _SHARED_PEEL_H_
#define _SHARED_PEEL_H_

#define HAD_M 0x1000
#define HAD_P 0x2000
#define HAP_DAT 0x4000
#define HAP_JNT 0x8000
#define IN_RF 0x10000
#define FIXED_FLAG  0x20000
#define HAP_FND  0x40000

#define X_MAT 0
#define X_PAT 1

#define X_MM_PM 0
#define X_MM_PP 1
#define X_MP_PM 2
#define X_MP_PP 3

#define PEEL_SIMPLE 1
#define PEEL_COMPLEX 2
#define PEEL_INITIAL 3
#define FENRIS_PEEL_SIMPLE 8

#define TL_NAME "#traitlocus#"

union Peelseq_Pointer {
  struct Complex_Element *complex;
  struct Simple_Element *simple;
  struct Initial_Element *initial;
  struct Fenris_Simple_Element *fsimple;
};

struct Peelseq_Head {
  union Peelseq_Pointer ptr;
  int type;
};

#define INITIAL_TRIPLET 1
#define INITIAL_DUPLET 2
#define INITIAL_DATA 3
#define INITIAL_FOUNDER 4

struct Initial_Element {
  struct Peelseq_Head next;
  int involved[3];           /* ids of above (i>0 - maternal allele of i, i<0 - paternal allele of -i) */
  int flags[3];              /* Flags for involved alleles */
  int rf_idx;               /* Included R-function (from simple peeling) */
  int type;
  int n_involved;          /* No. alleles involved in the operation */
  int out_index;           /* Index of output R-Function */
};

struct Simple_Element {
  struct Peelseq_Head next;
  int *off;
  int sire;
  int dam;
  int n_off;
  int pivot;
  int out_index;
};

struct Fenris_Simple_Element {
  struct Peelseq_Head next;
  int *off;
  int *rf;
  int sire;
  int dam;
  int n_off;
  int pivot;
  int out_index;
};

struct Complex_Element {
  struct Peelseq_Head next;
  int *involved;           /* ids of above (i>0 - maternal allele of i, i<0 - paternal allele of -i) */
  int *flags;              /* Flags for involved alleles */
  int *index;         /* Indices for R-Functions */
  int n_peel;              /* No. alleles to be peeled */
  int n_out;               /* No. output alleles */
  int n_involved;          /* No. alleles involved in the operation */
  int out_index;           /* Index of output R-Function */
  int n_rfuncs;            /* How many R-Functions to be combined */
};

void free_peelseq(struct Peelseq_Head *pp);
extern int num_bits(const int n_all);

#endif
