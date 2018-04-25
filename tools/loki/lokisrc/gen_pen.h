#ifndef _GEN_PEN_H_
#define _GEN_PEN_H_

struct gp_rfunc {
	struct gp_rfunc *next;
	int *inv;
	double *p;
	int n_inv,nc,flag;
};

struct gp_rfunc_ptr {
	struct gp_rfunc_ptr *next;
	struct gp_rfunc *rf;
};

struct gp_rfnode {
	struct gp_rfnode *next;
	struct gp_rfunc *rf;
	int y;
};

struct deg_list {
	struct deg_list *next,*prev,*abs_list;
	int deg,gene;
};

struct gp_peel_op {
	struct gp_peel_op *next;
	struct gp_rfunc *rfl;
	int *inv;
	int *pflag;
	int n_inv;
	int n_peel;
};

/* Maximum number of genes in an RFunction for gen_pen()
 * Note that memory and patience will run out long before this does...*/
#define MAX_GP_RF_GENES 64

extern double gen_pen(struct Locus *,int,int *,int,const struct loki *);
void alloc_gen_pen(int);
void free_gen_pen(void);

#endif
