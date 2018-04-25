#ifndef _LOKI_SCAN_H_
#define _LOKI_SCAN_H_

union var_ptr
{
	struct Marker *marker;
	struct Link *link;
	struct Variable *var;
};
	
struct lk_variable
{
	union var_ptr var;
	int type;
};

struct num_array {
	struct id_data *x;
	int size,ptr;
};

struct clause_atom {
	struct clause_atom *next;
	int i;
	char *v;
};

struct string_list {
	struct string_list *next;
	char *v;
};

#define LK_TYPE_MARKER 1
#define LK_TYPE_LINK 2
#define LK_TYPE_IDVAR 3
#define LK_TYPE_NONIDVAR 4

#define MAX_INCLUDE 16

void include_param_file(char *);
extern int iflag,list_ptr;
extern char *fname_list[MAX_INCLUDE+1];

#endif
