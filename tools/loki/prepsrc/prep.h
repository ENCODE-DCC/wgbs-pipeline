#ifndef _PREP_H_
#define _PREP_H_

#include "loki_struct.h"
#include "lk_long.h"
#include "string_utils.h"

struct sex_define {
	struct sex_define *next;
	struct parse_var_list *sex_exp;
	struct parse_var *sex_elem;
};

struct lk_ped {
	struct Id_Record *id_array;
	struct Family *family;
	struct label_data **ped_recode;
	struct label_data **family_recode;
	struct sex_define *sex_def;
	struct parse_var_list *ped_vars;
	int ped_size;
	int fam_flag;
};

struct lk_marker {
	struct lk_marker *next;
	struct parse_term *term;
	struct parse_term *haps[2];
	struct lk_link *link;
	int type;
};

struct lk_link {
	struct lk_link *next;
	string *name;
	int type;
};

struct Missing {
	struct Missing *next;
	struct parse_term *miss;
	struct parse_var_list *vl;
	char *scope;
};
	
struct lk_data {
	struct prep_file *infiles;
	struct Missing *missing;
};

struct lk_assign {
	struct lk_assign *next;
	struct parse_term *var;
	struct parse_term *arg;
};

struct lk_depend_node {
	union {
		struct lk_file *file;
		struct parse_var *vv;
		struct parse_var_elem *vv1;
	} elem;
	int type;
};

struct lk_depend {
	struct lk_depend *next;
	struct lk_depend_node *node;
};

struct lk_prep {
	struct lk_ped *pedigree;
	struct lk_marker *markers;
	struct lk_link *links;
	struct lk_data *data;
	struct lk_assign *models;
	struct lk_assign *assignments;
};

#define LK_MARKER 0
#define LK_QTL 1

#define SAFE_PTR(x) ((x)?(x):NULL_STR)
#endif
