#define ST_BLOCK_SIZE 256
#define ST_MIN_BITS 1
#define ST_MAX_BITS 10
#define INITIAL_NODE_BITS (ST_MIN_BITS)
#define INITIAL_NODE_SIZE (1<<(INITIAL_NODE_BITS))
#define MAX_FSP_RF_GENES 64
#define ST_FREE_STACK_SIZE 32

struct int_store {
	struct int_store *next;
	int *p;
	int ptr;
};

struct state_node1 {
	int *state; /* List of states of daughter nodes */
	int n; /* No. daughter nodes */
	int size;
};

struct state_node {
	union {
		struct state_node *next; /* List of pointers of daughter nodes */
		struct state_node1 *last; /* List of pointers of daughter nodes */
	} data;
	int *state; /* List of states of daughter nodes */
	int n; /* No. daughter nodes */
	int size;
};

struct deg_list {
	struct deg_list *next,*prev,*abs_list;
	int deg,gene;
};

struct gp_rfunc {
	struct gp_rfunc *next;
	int *inv;
	double *p;
	int n_inv,flag;
};

struct gp_rfunc_ptr {
	struct gp_rfunc_ptr *next;
	struct gp_rfunc *rf;
};

struct gp_rfnode {
	struct gp_rfnode *next;
	struct gp_rfunc *rf;
	int x;
};

double fseg_pen(int **,int **,double **,int *,double **,int,int,int);
void setup_fseg_pen(int,int);
