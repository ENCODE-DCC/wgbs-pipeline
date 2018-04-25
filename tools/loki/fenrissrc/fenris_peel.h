/* A cut down version of the Loki R-Functions */
struct FR_Func
{
	int *id_list;
	double *p;
	int n_ind;
	int flag;
};

/* Penetrance parameter structure */
struct pen_par
{
	double **freq;
	double *e;
};

struct triplet_peel
{
	int *kids;
	int ids,idd;
	int nkids;
	int nc;
	int g[4];
	double z[4];
};

typedef double fenris_pen_func(double *,int,int,struct pen_par *); /* p,id,locus,pen_params */

#define RESCALE_LIMIT 0.01

#define PEN_MODEL_EQUAL 0
#define PEN_MODEL_PROP 1
#define PEN_MODEL_EMP 2

fenris_pen_func fpen_emodel_equal; /* Equal mis-match probabilities */
fenris_pen_func fpen_emodel_prop; /* Mis-match probabilities proportional to genotype frequencies */
fenris_pen_func fpen_emodel_emp; /* Empirical model of Sobel et al. */

void estimate_freq(int,double *);
void peel_freq(int,struct Peelseq_Head *,int *,int,double *);
struct Peelseq_Head *peel_init(int *);
fenris_pen_func *get_pen_model(int);
double fenris_simple_peel(struct Fenris_Simple_Element *,int,double **,fenris_pen_func *,struct pen_par *);
double fenris_simple_distribute(struct Fenris_Simple_Element *,int,double **,double **,fenris_pen_func *,struct pen_par *);
double fenris_complex_peel(struct Complex_Element *,int,double **,fenris_pen_func *,struct pen_par *);
double fenris_complex_distribute(struct Complex_Element *,int,double **,double **,fenris_pen_func *,struct pen_par *);
double get_par_probs(double *,int,int,fenris_pen_func *,struct pen_par *,double **);
void fenris_peel_alloc(struct Peelseq_Head *);
void fenris_peel_free(void);
void free_peel_sequence(struct Peelseq_Head *);
void print_peelseq_element(FILE *,struct Peelseq_Head *);
void free_fenris_complex_peel(void);
void setup_fenris_complex_peel(void);

extern double **kval;
extern int *fc_state;
