#ifndef _MIN_DEG_H_
#define _MIN_DEG_H_

struct mat_elem {
	struct mat_elem *next;
	int x;
};

struct deg_list {
	struct deg_list *next,*last,*abs_list;
	int node;
};

struct pair_wt {
	int pair_node;
	double wt;
};

int *min_deg(int,int *,int *,int *,int);
int *greedy(int,int *,int *,int *,double *,struct pair_wt *,int,double *);

#endif
