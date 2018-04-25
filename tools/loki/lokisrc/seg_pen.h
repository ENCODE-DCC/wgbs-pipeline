#ifndef _SEG_PEN_H_
#define _SEG_PEN_H_

extern int pass_founder_genes1(const int,int,int,const struct loki *);
extern int pass_founder_genes1a(const int,const int *,const int,int,const struct loki *);
extern int pass_founder_genes1b(const int,const int *,const int,const struct loki *);
extern int pass_founder_genes2(const int locus,const int comp,int **seg,const struct loki *);
extern void pass_founder_genes(struct Locus *,const struct loki *);
extern double seg_pen(struct Locus *,int,int *,int,const struct loki *);
extern void seg_alloc(struct loki *);
extern void setup_obslist(struct loki *);
extern void seg_dealloc(struct loki *);
extern void pass_founder_genes_alloc(void);
extern void pass_founder_genes_dealloc(void);
extern void seg_init_freq(const struct Locus *,const struct loki *);
extern void seg_sample_freq(const struct Locus *,const struct loki *);
extern void seg_update_aff_freq(const struct Locus *,const struct loki *);

struct nuc_family {
	int *kids;
	int nkids;
};

struct cg_stack {
	int id;
	int kid_ptr;
	int par_flag;
};

#endif
