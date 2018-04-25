#ifndef _LOKI_PEEL_H_
#define _LOKI_PEEL_H_

#include "shared_peel.h"
#include "bin_tree.h"

#define OP_SAMPLING 128

#define SAMPLED_MAT 128
#define SAMPLED_PAT 256

#define LOCUS_SAMPLED 128
#define RFMASK_OK 256
#define TL_UNLINKED 1
#define TL_LINKED 2
#define TL_UPDATING 4

struct R_Func
{
  int n_ind;
  int n_terms;
  lk_ulong mask[2];
  int mask1[2];
  int *id_list;
  lk_ulong *index;
  double *p;
  int flag;
};

struct hash_data
{
  lk_ulong index;
  double p;
};

struct hash_block
{
  struct hash_block *next;
  struct bin_node *elements;
  struct hash_data *hd;
  int size,ptr;
};

struct peel_mem_block
{
  struct peel_mem_block *next;
  lk_ulong *index;
  double *val;
  size_t size,ptr;
};

struct fset
{
  int pat_gene[2];
  int mat_gene[2];
  double p;
};

struct peel_mem
{
  struct fset *s0;
  int *s1;
  double *s2;
  lk_ulong *s3;
  int *s4;
  void **s5;
  void **s6;
  double *s7;
  int *s8;
  struct R_Func *r_funcs;
};

struct lk_peel {
  struct peel_mem workspace;
  struct peel_mem_block *first_mem_block[2];
  struct peel_mem_block *mem_block[2];
  struct Peelseq_Head **peelseq_head;
  struct Peelseq_Head *peelseq_head_gen;
  struct Peelseq_Head *peelseq_head_gen_x;
  int max_peel_off; /* Max offspring in a single family peel op */
  int max_peel_ops; /* Max peel operations (single component) */
  int max_rfuncs; /* Max R-Functions (single component) */
  int max_out; /* Max total output individuals (single component) */
  int **tl_group;
  double **seg_count;
  double *aff_freq;
};

typedef void pen_func(double *,int,struct Locus *,int,int,const struct loki *);
typedef void trait_pen_func(double *,int,struct Locus *,const struct loki *);

extern pen_func penetrance;
extern trait_pen_func s_penetrance,s_penetrance1;
double q_penetrance(int,int,struct Locus *,const struct loki *);

/* Peeling output level - controlled by lower 3 bits in peel_trace */
#define TRACE_LEVEL_0 0
#define TRACE_LEVEL_1 1
#define TRACE_LEVEL_2 2
#define TRACE_LEVEL_3 3
#define TRACE_LEVEL_4 4
#define TRACE_MASK 7

#define CHK_PEEL(x) (((*peel_trace)&TRACE_MASK)>=(x))
#define MRK_MBLOCK 0
#define TRT_MBLOCK 1
#define MB_SIZE 4096

void prt_peel_trace(const int, const char *, ...);
double peel_locus(struct Locus **,int,int,int,struct loki *);
void free_complex_mem(void);
void set_sort_sex(const int);
void peel_alloc(struct loki *);
double loki_complex_peelop(const struct Complex_Element *,const int,const int,pen_func,const int,struct R_Func *,double **,struct loki *);
double loki_trait_complex_peelop(const struct Complex_Element *,const int,const int,struct R_Func *,trait_pen_func *,double **,struct loki *);
void get_rf_memory(struct R_Func *,size_t,int,struct loki *);
lk_ulong get_index1(int,int *,const int);
struct peel_mem_block *get_new_memblock(size_t,int);
#endif
