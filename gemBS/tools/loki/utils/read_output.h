struct ro_mark {
	struct ro_mark *next;
	char *name;
	double pos[2];
};

struct ro_link {
	struct ro_link *next;
	char *name;
	double map_r1[2];
	double map_r2[2];
	struct ro_mark *markers;
};

struct ro_cov {
	struct ro_cov *next;
	char *name;
};

struct ro_data {
	struct ro_link *link;
	struct ro_cov *cov;
	double total_map[2];
	char *model;
	int format;
	int sex_map;
	int n_genetic_groups;
	int n_cov_col;
	int n_qtl[2];
};

int read_header(FILE *,struct ro_data *,string **);

