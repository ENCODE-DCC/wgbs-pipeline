#ifndef _PARSER_H_
#define _PARSER_H_

#include "string_utils.h"

typedef union {
	int i;
	double x;
	string *str;
	struct parse_var *vv;
	struct parse_var_elem *vv1;
	struct expr_op *eop;
	struct func_op *fop;
} expr_elem;
	
typedef struct {
	char *fname;
	int line;
	int col;
} filepos;

struct deferred_assign {
	struct deferred_assign *next;
	struct parse_term *var;
	struct parse_term *expr;
};
	
struct parse_term {
	expr_elem elem;
	struct deferred_assign *defer;
	struct deferred_assign *clause_assign;
	int type;
};

struct format_clause {
	struct format_clause *next;
	struct format_clause *elem;
	int type;
	int multiplier;
};

struct parse_clause {
	struct parse_clause *next;
	union {
		struct parse_term *term;
		struct format_clause *format;
	} clause;
	int type;
};

struct parse_var {
	expr_elem arg;
	string *name;
	int type;
	int ltype;
	int size;
};

struct parse_var_elem {
	struct parse_var *head;
	struct parse_term *ex;
	int type;
};

struct parse_var_list {
	struct parse_var_list *next;
	struct parse_term term;
	struct parse_clause *clause;
	filepos pos;
};

struct expr_op {
	struct parse_term *exp1,*exp2;
	int op;
};

typedef struct parse_f_def {
	struct parse_f_def *next;
	char *com;
	union {
		double (*f_real)(double);
		int (*f_int)(int);
		string *(*f_string)(string *);
	} f;
	int type;
} func_def;

struct func_op {
	struct parse_term *ex;
	func_def *f;
};

typedef struct {
	void (*assign)(struct parse_term *,struct parse_term *);
	void (*ctypeI)(char *,struct parse_clause *,struct parse_term *,struct parse_var_list *);
	void (*ctypeII)(int,char *,struct parse_term *,struct parse_var_list *);
	void (*ctypeIII)(char *,struct parse_term *,struct parse_term *);
	void (*ctypeV)(char *,int,struct parse_var_list *);
	int *line;
	char *fname;
	/* Starting positions of last command, and beginning of last token */
	filepos start_pos;
	filepos begin_tok;
} parse_handler;

int Parse_File(FILE *,char **,char **,func_def *,parse_handler *);
func_def *register_string_function(func_def *,string *(*)(string *),char *);
func_def *register_int_function(func_def *,int (*)(int),char *);
func_def *register_real_function(func_def *,double (*)(double),char *);
void set_typeI_kludge(void);
void free_functions(func_def *);
void init_parse_utils(parse_handler *);
void do_deferred(struct parse_term *);
void print_parse_clause(FILE *,struct parse_clause *);
string *parse_clause_str(struct parse_clause *);
string *convert_to_string(struct parse_term *);
string *format_clause_str(struct format_clause *);
struct parse_term *check_immediate(struct parse_term *);
struct parse_term *do_func(func_def *,struct parse_term *);
struct parse_term *get_term_list(struct parse_term *,int *);
struct parse_term *try_resolve_expr(struct parse_term *);
struct parse_term *parse_make_texp(int,double,string *,int);
struct parse_term *do_op(struct parse_term *,struct parse_term *,int);
struct parse_term *do_op1(struct parse_term *,struct parse_term *,int);
struct parse_term *do_op2(struct parse_term *,struct parse_term *,int,int);
struct parse_term *get_var_term(struct parse_var *);
struct parse_term *get_new_term(void);
struct parse_var *find_parse_var(struct parse_term *);
struct format_clause *make_format_clause(int,struct format_clause *);
struct parse_clause *convert_format_clause(struct format_clause *);
struct parse_clause *make_clause(struct parse_term *);
int convert_numeric(struct parse_term *);
void assign_var(struct parse_term *,struct parse_term *);
void free_parse_term(struct parse_term *);
void free_parse_var(struct parse_var *);
struct parse_term *copy_parse_term(struct parse_term *);
void print_all_vars(void);
void parerror_fp(const filepos *);
void parerror1(const char *, ...);
void parerror2(const char *, ...);
void parerror3(filepos *,const char *, ...);
void print_parse_term_name(FILE *,struct parse_term *);
struct parse_var *insert_prep_var(string *);
string *parse_term_name(struct parse_term *);
string *var_list_string(struct parse_var_list *);
struct parse_term *get_array_term(struct parse_var *,struct parse_term *);
struct parse_var *get_array_var(struct parse_var *,struct parse_term **);
struct parse_var_list *get_array_var_list(struct parse_var *,struct parse_term *);
struct parse_var_list *make_var_list(struct parse_term *,struct parse_clause *);
struct parse_var_list *addto_var_list(struct parse_var_list *,struct parse_var_list *);
struct parse_var_list *make_loop_clause(struct parse_var_list *,struct parse_var *,struct parse_term *);
string *check_string(struct parse_term *);
int convert_string(struct parse_term *);
int convert_to_int(struct parse_term *,int *);
void free_parse_var_list(struct parse_var_list *);
void print_var_list(FILE *,struct parse_var_list *);
void print_parse_exp(FILE *,struct parse_term *);
string *parse_exp_str(struct parse_term *);
char *parse_exp_cstr(struct parse_term *);
void print_parse_op(FILE *,int);
void free_parse_clause(struct parse_clause *);
void com_ctypeI(struct parse_var *,struct parse_clause *,struct parse_term *,struct parse_var_list *);
void com_ctypeII(int, struct parse_var *,struct parse_term *,struct parse_var_list *);
void com_ctypeIII(struct parse_var *,struct parse_term *,struct parse_term *);
void com_ctypeV(struct parse_var *,int,struct parse_var_list *);
void com_assign(struct parse_term *,struct parse_term *);
void handle_clause(struct parse_clause *);
void set_var_list_flag(struct parse_var_list *,int);
void set_var_list_vector_flag(struct parse_var_list *);
struct parse_term *parse_incr_decr(struct parse_term *,int);
struct parse_var_list *merge_var_varlist(struct parse_var *,struct parse_clause *,struct parse_var_list *);
string *parse_op_str(int);

#define PREFLAG 0x4000
#define POSTFLAG 0x8000

#define PRE_INCR (PREFLAG|'+')
#define PRE_DECR (PREFLAG|'-')
#define POST_INCR (POSTFLAG|'+')
#define POST_DECR (POSTFLAG|'-')

#define VT_TYPES 0x0fff
#define VT_VECTOR 0x10000
#define VT_IMMED 0x20000
#define VT_SCALAR 0x40000
#define VT_ARRAY 0x80000
#define VT_ARRAY_ELEMENT 0x100000
#define VT_SIZE_SET 0x200000

#define var_name(x) (get_cstring(parse_term_name(x)))
#define IS_VAR(x) ((x->type&VT_TYPES)==VARIABLE)

#define FC_SKIP 0
#define FC_READ 1
#define FC_SUBCLAUSE 2

#define PC_TERM 1
#define PC_FORMAT_CLAUSE 2

#define WORKING_ZERO (1.0e-16) /* For deciding if a float can be converted to a real */
#define MAX_ARRAY_IDX 0x100000 /* For sanity checking */

#endif

