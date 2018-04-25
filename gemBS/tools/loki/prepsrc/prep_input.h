#ifndef _PREP_FILE_COM_H_
#define _PREP_FILE_COM_H_

struct prep_file {
	struct prep_file *next;
	char *name;
	struct parse_var_list *varlist;
	int *fixed_start,*fixed_end;
	string *fs,*rs,*gs;
	int skip;
	int ncol;
	struct DataBlock *data;
};

#define N_SPECIAL_VAR 4

struct special_vars {
	struct parse_var *var[N_SPECIAL_VAR];
	string *string[N_SPECIAL_VAR];
	int type[N_SPECIAL_VAR];
};

#define WAIT_LINK 1
#define WAIT_FILE 2

struct typeI_wait {
	int type;
	union {
		struct lk_link *lk;
		struct prep_file *file;
	} ptr;
};

struct prep_file *prep_file_com(struct parse_clause *,struct parse_term *,struct parse_var_list *,struct special_vars *,struct loki *);

#endif
