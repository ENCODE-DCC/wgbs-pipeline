struct relate_atom {
	int type;
	int n;
	int deg[2];
};

struct relate_off {
	struct relate_off *next;
	struct relate_atom *atoms;
	int x;
	int n;
	double p;
};

#define REL_HALF 1
#define REL_FULL 2


