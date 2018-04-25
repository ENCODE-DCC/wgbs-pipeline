#ifndef _SCAN_H_
#define _SCAN_H_

/****************************************************************************
 *                                                                          *
 *     Loki - Programs for genetic analysis of complex traits using MCMC    *
 *                                                                          *
 *             Simon Heath - University of Washington                       *
 *                                                                          *
 *                       March 1997                                         *
 *                                                                          *
 * scan.h:                                                                  *
 *                                                                          *
 ****************************************************************************/

#include "loki_struct.h"
#include "bin_tree.h"

extern int symbol_lookup(char *,int);
extern void print_scan_warn(char *, ...);
extern void print_scan_err(char *, ...);
extern void print_version_and_exit(void);
extern void ReadData(char *);
extern void init_stuff(char **);
extern char *get_marker_name(const int);

union arg_type {
	long value;
	double rvalue;
	struct bin_node *var;
	struct var_element *element;
	char *string;
};

struct express {
	union arg_type arg;
	int type;
};

struct var_element {
	union arg_type arg;
	int n_levels,index,oindex,type;
};

struct scan_data {
	struct var_element *element;
	char *name;
	int n_elements; /* 1 for scalars, >=1 for arrays */
	int vtype;
};

struct var_list {
	struct var_list *next;
	struct bin_node *var;
	int index;
};

struct label_data {
	union	{
		char *string;
		long value;
	} data;
	int type,index,balance;
};

struct InFile {
	struct InFile *next;
	char *name;
	int shell_flag;
	int nvar,ncol,n_records,id_col,family_col;
	struct var_element **element;
	struct format *format;
	struct fformat *fformat;
	struct DataBlock *data;
};

struct Id_Record {
	int *kids,*famlist;
	int *haplo[2];
	int *mkflag;
	struct id_data *data,**data1;
	int sire,dam;
	int nfam,nkids;
	int fam_code;
	int rf_idx;
	int sex,affected,proband;
	int family;
	int component;
	int group;
	int flag,nrec,order;
	int fc_flag;
	int nhaps[2],ngens,ngens1,sg[2];
	int allele[2];
};

struct Family {
	int *kids;
	int sire,dam;
	int nkids;
};

struct format {
	int n_atoms,line;
	struct format_atom *f_atoms;
};

struct Marker {
	int index;
	int link;
	int **order,*o_size;
	int **allele_trans;
	int n_sub_elements;
	int pos_set[3];
	double pos[3];
	struct scan_data *var;
	struct var_element *element,*hap_element[2];
	struct var_element **sub_element;
};

struct format_clause {
	int n_atoms,fc_size;
	struct format_atom **f_atoms;
};

struct fformat {
	char *fs;
	char *rs;
	char *gs;
	int skip;
};

struct format_atom {
	int size,pos;
};

struct model {
	struct model *next;
	struct scan_data *trait;
	struct model_list *model_list;
	int index;
};

struct model_list {
	struct model_list *next;
	struct var_element **element;
	int nvar;
};
		
struct Link {
	struct Link *next;
	char *name;
	struct var_element **element;
	int n_loci;
	int type;
};

struct Miss {
	struct Miss *next;
	struct express Missing;
	struct var_element **element;
	char *scope;
	int nvar;
};

struct op_stack {
	union arg_type arg;
	int type;
};

struct operation {
	struct operation *next;
	union arg_type arg;
	int type,op;
};

struct Restrict {
	struct Restrict *next;
	struct operation *Op_List;
	struct var_element **element;
	int nvar;
	int flag;
};

struct Censor {
	struct Censor *next;
	struct operation *Op_List;
	struct var_element *element;
};

struct gt_data {
	struct label_data *node1,*node2;
};

union DataRec {
	struct label_data *node;
	struct gt_data *gt_data;
	long value;
	double rvalue;
};

struct DataBlock {
	struct DataBlock *next;
	unsigned char *type;
	union DataRec *records;
	int blocksize,record_ptr;
};

struct recode_table_tag {
	struct label_data *node;
	int index;
};

struct sex_def {
	struct sex_def *next;
	struct express *sex_exp[2];
	struct var_element *sex_elem;
};

extern void free_nodes(void);
extern void setup_pedigree(int,struct recode_table_tag *,char *);
extern void recode_factors(int,struct recode_table_tag *);
extern void recode_marker_data(int,struct recode_table_tag *);
extern int check_missing(int,int,int,struct DataBlock *);
extern void Check_Inbr(char *);
extern void check_vars(struct bin_node *,int *,void (*func)(struct bin_node *,int *));
extern void free_infile(struct InFile *);
extern void prune_pedigree(char *);
extern void count_components(char *);
extern void free_restrict(struct Restrict *);
extern void affected_data(void);
extern void proband_data(void);
extern void restrict_data(void);
extern void print_orig_alleles(FILE *,const int,const int,const int);
extern void censored_data(void);
extern void cleanup_unused(void);
extern void count_relationships(void);
extern void Output_Data(void);
extern void Output_Raw_Data(void);
extern void print_orig_id(FILE *,const int,const int);
extern void print_orig_id1(FILE *,const int,const int);
extern int print_orig_family(FILE *,const int,const int);
extern void print_orig_triple(FILE *,const int);
extern void match_records(void);
extern void free_var(struct bin_node *);
extern void WriteReport(char *);
extern int Genotype_Elimination(int,char *,int);
extern void InitFamilies(char *);
extern void check_ymark(void);
extern void free_op(struct operation *);
extern void yyerror(char *);
extern void yyerror1(char *);
extern void check_sex(void),check_sex2(void),count_loops(char *),handle_groups(char *),check_inbreeding(char *);
extern struct label_data *find_node(const void *,int,int);
extern int ReadControl(FILE *,char *,char **);

extern struct InFile *Infiles;
extern struct bin_node *root_var;
extern struct Link *links;
extern struct Miss *Miss;
extern struct Restrict *Restrictions;
extern struct Censor *Censored;
extern struct operation *Affected,*Unaffected,*Proband;
extern struct model *Models;
extern struct Marker *markers,*traitlocus;
extern struct Id_Record *id_array;
extern struct var_element **var_factors;
extern struct label_data **ped_recode,**family_recode,**group_recode,***factor_recode;
extern struct Family *family;
extern struct remember *RemBlock;

extern int *rec_tab,*rec_tab1;
extern char linebuf[],linebuf1[],*Filter,*gsformat,*rsformat,*fsformat,*OutputFile,*OutputLaurFile,*OutputRawFile,*ErrorDir;
extern int nfiles,scan_error,scan_error_n,scan_warn_n,file_skip;
extern int lineno,lineno1,at_file,at_model,at_use,begin_comm;
extern int max_scan_errors,max_scan_warnings,nrm_flag,n_genetic_groups;
extern int ped_size,pruned_ped_size,*ped_recode1,n_factors,n_markers,n_id_records,n_nonid_records;
extern int n_comp,*comp_size,n_orig_families,n_families,verbose,family_id,sig_caught,catch_sigs,sig_quiet;
extern struct sex_def *sex_def;
extern struct var_element *group_elem;

#define MAX_LOOP 16
#define MAX_INCLUDE 16
#define LINEBUFSIZE 4096

extern int loop_level,loop_ptr[MAX_LOOP],loop_stat[MAX_LOOP],loop_record,loop_stack_size,comp_sflag;
extern int loop_main_ptr,in_loopclause,loop_clause_end,loop_clause_step,loop_clause_ptr,trace_peel,strip_vars;
extern struct var_element *loop_clause_element,**id_elements,**nonid_elements;

extern FILE *yyin;

#define SCAN_ERR 1
#define PED_ERR 8
#define FILE_ERR 16
#define FORMAT_ERR 32
#define LINK_ERR 64

#define NUM_PREP_SYSTEM_VAR 16

extern int syst_var[NUM_PREP_SYSTEM_VAR];

#define PRUNE_OPTION 0
#define RECODE_OPTION 1
#define NO_EXTRA_ALLELE 2
#define PEEL_OPTION 3
#define TRACE_RESTRICT 4
#define TRACE_CENSORED 5
#define TRACE_AFFECTED 6
#define CORRECT_ERRORS 7
#define PEEL_TRACE 8
#define MULTIPLE_RECORDS 9
#define MULTIVARIATE_TEST 10
#define ERROR_CHECK 11
#define NO_DEFAULT_MISSING 12
#define SKIP_BAD_REALS 13
#define SKIP_BAD_INTS 14
#define IGNORE_CASE 15

#endif
