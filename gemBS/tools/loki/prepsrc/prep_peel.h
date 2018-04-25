#ifndef _PREP_PEEL_H_
#define _PREP_PEEL_H_

#include "lk_long.h"
#include "shared_peel.h"
#include "bin_tree.h"

#define HAS_DATA 1
#define IS_PRUNED 2
#define IS_FIXED 4
#define RF_INCLUDED 8
#define WAS_PIVOT 16
#define PEELED 32
#define IS_SINGLETON 64
#define HAS_GDATA 128
#define STABLE_FLAGS 255
#define HAP_MAT 256
#define HAP_PAT 512
#define HAP_M 1024
#define HAP_P 2048

#define CHK_PEEL(x) ((trace_peel&TRACE_MASK)>=x)
#define TRACE_LEVEL0 0
#define TRACE_LEVEL1 1
#define TRACE_LEVEL2 2
#define TRACE_LEVEL3 3
#define TRACE_LEVEL4 4
#define TRACE_MASK 7

struct R_Func
{
	int n_ind;
	int n_terms;
	int flag;
	lk_ulong mask;
	lk_ulong mask1;
	int *id_list;
	lk_ulong *index;
	struct Peelseq_Head *peel_elem;
};

struct hash_block
{
	struct hash_block *next;
	struct bin_node *elements;
	lk_ulong *idx;
	int size,ptr;
};

extern int do_peel_op(const struct Complex_Element *element,struct R_Func *r_func,const int n_all,const int id,lk_ulong **all_set,lk_ulong *req_set[]);
extern void print_orig_allele_id(FILE *fptr,const int i);
extern void free_hash_blocks(void);
extern int total_terms,total_comb,trace_restrict;

#endif
