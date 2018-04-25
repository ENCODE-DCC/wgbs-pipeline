%{
/****************************************************************************
 *                                                                          *
 *     Loki - Programs for genetic analysis of complex traits using MCMC    *
 *                                                                          *
 *             Simon Heath - University of Washington                       *
 *                                                                          *
 *                       July 1997                                          *
 *                                                                          *
 * param_parse.y:                                                           *
 *                                                                          *
 * yacc source for parameter file parser.                                   *
 *                                                                          *
 ****************************************************************************/

#include <config.h>
#include <stdlib.h>
#include <string.h>
#ifdef USE_DMALLOC
#include <dmalloc.h>
#endif
#include <stdarg.h>
#include <stdio.h>
#include <math.h>
#include <sys/time.h>
	
#include "utils.h"
#include "loki.h"
#include "loki_scan.h"
#include "loki_ibd.h"
#include "loki_utils.h"
#include "lk_malloc.h"
#include "shared_peel.h"
#include "mat_utils.h"
#include "meiosis_scan.h"
#include "ranlib.h"
#include "snprintf.h"

static struct loki *loki;
static struct Marker *freq_marker;
static int check_variance(double,int);
static void print_scan_warn(char *, ...);
static void set_position(struct lk_variable *, double, double);
static struct lk_variable *find_var(char *, int, int);
static struct Marker *check_marker(struct lk_variable *);
static int find_allele(char *, struct Marker *);
static int find_group(char *,int, int);
static int find_trait(struct lk_variable *);
static void set_output_gen(char *,char *);
static struct IBD_List *add_ibd_list(double,struct IBD_List *);
static void set_freq(struct Marker *, double,int);
static void set_map_range(char *,double,double, int);
static void set_tloci(int,int);	
static void set_ibd_list(char *,struct IBD_List *,int);
static void set_ibd_markers(char *);
static void set_output(struct lk_variable *);
static void set_group_order(int);
static void set_analyze(char *);
static void set_ibd_mode(char *);
static void set_pseudo_chrome(struct string_list *,struct clause_atom *);
static struct num_array *make_num_array(int);
static void free_num_array(struct num_array *);
static void add_real_to_num_array(struct num_array *,double);
static void add_int_to_num_array(struct num_array *,int);
static void set_syst_var(int,struct num_array *);	
static int group_ptr,*group_order,group_counter,freq_allele,c_flag;
static struct clause_atom *new_clause_atom(char *,int);
static struct string_list *new_string_list(char *);

extern void yyerror(char *s),print_scan_err(char *fmt, ...);
extern int yyparse(void),yylex(void),lineno,lineno1,tokenpos;
extern char *yytext,linebuf[];

extern FILE *yyin;
int scan_error_n,iflag;
static int max_scan_errors=30,n_gen_grps,start_flag;
static int scan_warn_n,max_scan_warnings=30;
%}

%union  {
	char *string;
	int value;
	double rvalue;
	struct IBD_List *rlist;
	struct lk_variable *lk_var;
	struct num_array *num_array;
	struct clause_atom *clause;
	struct string_list *string_list;
}

%token RESIDUAL GENETIC VARIANCE POSITION FREQUENCY VIRTUAL
%token START MEAN ITERATIONS SAMPLE FROM OUTPUT MAP TOTAL
%token SEED SFILE SEEDFILE TRAIT LOCI SET SYSTEM_VAR TIMECOM
%token ESTIMATE IBD GROUP ORDER MALE FEMALE LIMIT AFFECTED
%token PHENO GENO COUNTS DUMP TYPE ANALYZE NORMAL STUDENT_T
%token HAPLO INCLUDE FUNCTION HALDANE KOSAMBI POLYGENIC
%token MARKERS GRID COMPRESS DIR PSEUDO
  
%token <string> STRING
%token <value> INTEGER SYSTEM_VAR
%token <rvalue> REAL
%type <rvalue> rnum
%type <num_array> array
%type <rlist> ibdlist
%type <lk_var> lkvar
%type <value> allele group trait_var
%type <clause> clause clause_atom
%type <string_list> pchrom_list pchrom_list_atom
  
%%

parmfile: {lineno1=lineno;} command1 {iflag=0;}
       | error {iflag=0;}
       | parmfile {lineno1=lineno;} command1 {iflag=0;}
       | parmfile error {iflag=0;}
       ;

command1: command
       | command_a
       | START {start_flag=2;} command {start_flag=1;}
       ;

command: resvarcommand
       | addvarcommand
       | positioncommand
       | frequencycommand
       | meancommand
       ;

command_a: itercommand
       | samplecommand
       | outputcommand
       | mapcommand
       | seedcommand
       | tlocicommand
       | setcommand
       | ibdcommand
       | groupcommand
       | limitresvarcommand
       | limitaddvarcommand
       | analyzecommand
		 | includecommand
       | aff_freqcommand
       | limit_timecommand
		 | compresscommand
       | pseudocommand
       ;

samplecommand: SAMPLE FROM INTEGER {loki->params.sample_from[0]=loki->params.sample_from[1]=$3;}
       | SAMPLE FROM INTEGER opt_comma INTEGER {loki->params.sample_from[0]=$3; loki->params.sample_from[1]=$5;}
       | START OUTPUT INTEGER {loki->params.sample_from[0]=loki->params.sample_from[1]=$3;}
       | START OUTPUT INTEGER opt_comma INTEGER {loki->params.sample_from[0]=$3; loki->params.sample_from[1]=$5;}
       ;

setcommand: SET SYSTEM_VAR array {set_syst_var($2,$3); free_num_array($3);}
       | SET STRING array { free($2); free_num_array($3); yyerror("Unrecognized system variable"); }
       ;

includecommand: INCLUDE {iflag=1;} STRING {include_param_file($3);}
       ;

opt_comma:
       | ','
       ;

pchrom_list_atom: STRING {$$=new_string_list($1);}
       ;

pchrom_list: pchrom_list_atom
       | pchrom_list ',' pchrom_list_atom { $3->next=$1; $$=$3; }
       ;

clause_atom: STRING '=' INTEGER {$$=new_clause_atom($1,$3);}
       | START '=' INTEGER {$$=new_clause_atom(strdup("INIT"),$3);}
       | FREQUENCY '=' INTEGER {$$=new_clause_atom(strdup("FREQ"),$3);}
       ;

clause: clause_atom
       | clause ';' clause_atom { $3->next=$1; $$=$3; }
       ;

pseudocommand: PSEUDO {set_pseudo_chrome(0,0);}
       | PSEUDO pchrom_list {set_pseudo_chrome($2,0);}
       | PSEUDO '[' clause ']' {set_pseudo_chrome(0,$3);}
       | PSEUDO '[' clause ']' pchrom_list {set_pseudo_chrome($5,$3);}
       ;

ibdcommand: ESTIMATE IBD STRING opt_comma ibdlist { set_ibd_list($3,$5,IBD_EST_DISCRETE); free($3);}
       | ESTIMATE IBD ibdlist { set_ibd_list(0,$3,IBD_EST_DISCRETE); }
       | ESTIMATE IBD MARKERS STRING { set_ibd_markers($4); free($4);}
       | ESTIMATE IBD MARKERS { set_ibd_markers(0); }
       | ESTIMATE IBD GRID STRING opt_comma ibdlist {set_ibd_list($4,$6,IBD_EST_GRID); free($4);}
       | ESTIMATE IBD GRID ibdlist { set_ibd_list(0,$4,IBD_EST_GRID); }
       ;

compresscommand: COMPRESS IBD OUTPUT {loki->params.compress_ibd=1;}
		 | COMPRESS OUTPUT IBD {loki->params.compress_ibd=1;}
		 ;
		 
aff_freqcommand: ESTIMATE AFFECTED FREQUENCY {loki->params.est_aff_freq=1;}
       ;

analyzecom: STRING {set_analyze($1); free($1);}
       | AFFECTED {set_analyze("AFFECTED");}
       | IBD {set_analyze("IBD");}
       ;

analyzelist: analyzecom
       | analyzelist ',' analyzecom
       ;

analyzecommand: ANALYZE {loki->params.analysis=0;} analyzelist
       ;

limit_timecommand: TIMECOM LIMIT rnum {loki->params.limit_time=$3,loki->params.limit_timer_type=ITIMER_REAL;}
       | LIMIT TIMECOM rnum {loki->params.limit_time=$3,loki->params.limit_timer_type=ITIMER_REAL;}
       | LIMIT VIRTUAL TIMECOM rnum {loki->params.limit_time=$4,loki->params.limit_timer_type=ITIMER_VIRTUAL;}
       | VIRTUAL TIMECOM LIMIT rnum {loki->params.limit_time=$4,loki->params.limit_timer_type=ITIMER_VIRTUAL;}
       ;

seedcommand: SEEDFILE STRING {if(loki->names[LK_SEEDFILE]) free(loki->names[LK_SEEDFILE]); loki->names[LK_SEEDFILE]=$2;}
       | SEEDFILE STRING ',' INTEGER {if(loki->names[LK_SEEDFILE]) free(loki->names[LK_SEEDFILE]); loki->names[LK_SEEDFILE]=$2; if($4) loki->params.ranseed_set|=1; else loki->params.ranseed_set&=~1;}
       | SEED SFILE STRING {if(loki->names[LK_SEEDFILE]) free(loki->names[LK_SEEDFILE]); loki->names[LK_SEEDFILE]=$3;}
       | SEED SFILE STRING ',' INTEGER {if(loki->names[LK_SEEDFILE]) free(loki->names[LK_SEEDFILE]); loki->names[LK_SEEDFILE]=$3; if($5) loki->params.ranseed_set|=1; else loki->params.ranseed_set&=~1;}
       | SEED INTEGER {
			 if($2<0) yyerror("Seedvalue out of range");
			 else {
				 init_ranf($2);
				 loki->params.ranseed_set|=2;
			 }
		 }
       ;
mapcommand: MAP STRING rnum ',' rnum {set_map_range($2,$3,$5,-1); free($2);}
       | MALE MAP STRING rnum ',' rnum {set_map_range($3,$4,$6,X_PAT); free($3); }
       | FEMALE MAP STRING rnum ',' rnum {set_map_range($3,$4,$6,X_MAT); free($3); }
       | TOTAL MAP rnum {set_map_range(0,$3,$3,-1);}
       | TOTAL MAP rnum ',' rnum {set_map_range(0,$3,$5,-2);}
       | TOTAL MALE MAP rnum {set_map_range(0,$4,$4,X_PAT);}
       | TOTAL FEMALE MAP rnum {set_map_range(0,$4,$4,X_MAT);}
       | MAP TOTAL rnum {set_map_range(0,$3,$3,-1);}
       | MAP TOTAL rnum ',' rnum {set_map_range(0,$3,$5,-2);}
       | MALE MAP TOTAL rnum {set_map_range(0,$4,$4,X_PAT);}
       | FEMALE MAP TOTAL rnum {set_map_range(0,$4,$4,X_MAT);}
       | MAP FUNCTION HALDANE {loki->params.map_function=MAP_HALDANE;}
       | MAP FUNCTION KOSAMBI {loki->params.map_function=MAP_KOSAMBI;}
       ;

tlocicommand: TRAIT LOCI INTEGER {set_tloci(-1,$3);}
       | START TRAIT LOCI INTEGER {set_tloci(-2,$4);}
       | TRAIT LOCI INTEGER ',' INTEGER {set_tloci($3,$5);}
       | TRAIT LOCI MEAN rnum {loki->models->tloci_mean=$4; loki->models->tloci_mean_set=1;}
       ;

itercommand: ITERATIONS INTEGER { loki->params.num_iter=$2; };

outputcommand: OUTPUT FREQUENCY INTEGER {loki->params.sample_freq[0]=loki->params.sample_freq[1]=$3;}
       | OUTPUT FREQUENCY INTEGER ',' INTEGER {loki->params.sample_freq[0]=$3; loki->params.sample_freq[1]=$5;}
       | OUTPUT FREQUENCY STRING {if(loki->names[LK_FREQFILE]) free(loki->names[LK_FREQFILE]); loki->names[LK_FREQFILE]=$3;}
       | SAMPLE FREQUENCY INTEGER {loki->params.sample_freq[0]=loki->params.sample_freq[1]=$3;}
       | SAMPLE FREQUENCY INTEGER ',' INTEGER {loki->params.sample_freq[0]=$3; loki->params.sample_freq[1]=$5;}
       | OUTPUT PHENO STRING {loki->names[LK_PHENFILE]=$3; }
       | OUTPUT GENO STRING {set_output_gen($3,0);} 
       | OUTPUT GENO STRING ',' STRING {set_output_gen($3,$5); free($5); }
       | OUTPUT lkvarlist
       | OUTPUT TYPE INTEGER {loki->params.output_type=$3;}
       | OUTPUT SFILE STRING {if(loki->names[LK_OUTPUTFILE]) free(loki->names[LK_OUTPUTFILE]); loki->names[LK_OUTPUTFILE]=$3;}
       | OUTPUT POSITION SFILE STRING {if(loki->names[LK_POSFILE]) free(loki->names[LK_POSFILE]); loki->names[LK_POSFILE]=$4;}
       | OUTPUT IBD SFILE STRING {if(loki->names[LK_IBDFILE]) free(loki->names[LK_IBDFILE]); loki->names[LK_IBDFILE]=$4;}
       | OUTPUT IBD DIR STRING {if(loki->names[LK_IBDDIR]) free(loki->names[LK_IBDDIR]); loki->names[LK_IBDDIR]=$4;}
       | DUMP SFILE STRING {if(loki->names[LK_DUMPFILE]) free(loki->names[LK_DUMPFILE]); loki->names[LK_DUMPFILE]=$3;}
       | DUMP FREQUENCY INTEGER {loki->params.dump_freq=$3;}
		 | OUTPUT HAPLO STRING {loki->params.output_haplo=1;if(loki->names[LK_HAPLOFILE]) free(loki->names[LK_HAPLOFILE]); loki->names[LK_HAPLOFILE]=$3;}
		 | OUTPUT HAPLO {loki->params.output_haplo=1;}
		 | OUTPUT POLYGENIC STRING {if(loki->names[LK_POLYFILE]) free(loki->names[LK_POLYFILE]); loki->names[LK_POLYFILE]=$3;}
		 | OUTPUT IBD STRING {set_ibd_mode($3); free($3);}
		 | OUTPUT IBD STRING ',' STRING {set_ibd_mode($3); free($3); set_ibd_mode($5); free($5);}
		 ;

resvarcommand: RESIDUAL VARIANCE rnum { if(!check_variance($3,0)) {loki->models->res_var_set[0]=start_flag; loki->models->residual_var[0]=$3;} }
       | RESIDUAL VARIANCE trait_var rnum { if($3>=0 && !check_variance($4,1)) {loki->models->res_var_set[$3]=start_flag; BB(loki->models->residual_var,$3,$3)=$4;} }
       ;

limitresvarcommand: RESIDUAL VARIANCE LIMIT rnum { if(!check_variance($4,0)) loki->models->residual_var_limit[0]=$4; }
          | RESIDUAL VARIANCE LIMIT trait_var rnum { if($4>=0 && !check_variance($5,0)) loki->models->residual_var_limit[$4]=$5; }
          | LIMIT RESIDUAL VARIANCE rnum { if(!check_variance($4,0)) loki->models->residual_var_limit[0]=$4; }
          | LIMIT RESIDUAL VARIANCE trait_var rnum { if($4>=0 && !check_variance($5,0)) loki->models->residual_var_limit[$4]=$5; }
          ;

addvarcommand: GENETIC VARIANCE rnum { if(!check_variance($3,0)) {loki->models->add_var_set[0]=start_flag; loki->models->additive_var[0]=$3;} }
          | GENETIC VARIANCE trait_var rnum { if($3>=0 && !check_variance($4,0)) {loki->models->add_var_set[$3]=start_flag; BB(loki->models->additive_var,$3,$3)=$4;} }
          ;

limitaddvarcommand: GENETIC VARIANCE LIMIT rnum { if(!check_variance($4,0)) loki->models->additive_var_limit[0]=$4; }
          | GENETIC VARIANCE LIMIT trait_var rnum { if($4>=0 && !check_variance($5,0)) loki->models->additive_var_limit[$4]=$5; }
          | LIMIT GENETIC VARIANCE rnum { if(!check_variance($4,0)) loki->models->additive_var_limit[0]=$4; }
          | LIMIT GENETIC VARIANCE trait_var rnum { if($4>=0 && !check_variance($5,0)) loki->models->additive_var_limit[$4]=$5; }
          ;

meancommand: MEAN rnum { if(loki->models->n_models>1) yyerror("Model must be specified when multiple models are present");
	                      else if(!loki->models->n_models) print_scan_warn("No model present - MEAN command ignored\n");
	                      else {loki->models->grand_mean_set[0]=start_flag; loki->models->grand_mean[0]=$2; } }
          | MEAN trait_var rnum { if(!loki->models->n_models) print_scan_warn("No model present - MEAN command ignored\n");
				                      else if($2>=0) {loki->models->grand_mean_set[$2]=start_flag; loki->models->grand_mean[$2]=$3;} }
          ;

positioncommand: POSITION lkvar rnum { set_position($2,$3,$3); }
        | POSITION lkvar rnum ',' rnum { loki->markers->sex_map=1; set_position($2,$3,$5); }
        ;

frequencycommand: FREQUENCY {c_flag=0;} lkmarker freqlist
		  | FREQUENCY {c_flag=1;} COUNTS lkmarker freqlist
        | COUNTS {c_flag=1;} lkmarker freqlist
        ;

trait_var: lkvar {$$=find_trait($1);} ;

lkvar: STRING { $$=find_var($1,0,0); }
       | STRING '(' INTEGER ')' { $$=find_var($1,$3,1); }
       ;

lkvarlist: lkvar { if($1) {set_output($1); free($1);} }
       | lkvarlist ',' lkvar { if($3) {set_output($3); free($3);} }
       ;

lkmarker: lkvar { freq_marker=check_marker($1); } ;

groupcommand: GROUP ORDER {group_ptr=0;} grouplist {if(group_ptr<n_gen_grps) print_scan_err("Line %d: Too few groups in order statement\n",lineno1);}
       ;

group: INTEGER { $$=find_group(yytext,$1,1); }
       | REAL { $$=find_group(yytext,0,0); }
       | STRING { $$=find_group($1,0,0); free($1);}
       ;

grouplist: group {set_group_order($1);}
       | grouplist ',' group {set_group_order($3);}
       ;

allele: INTEGER { $$=find_allele(yytext,freq_marker); }
       | REAL { $$=find_allele(yytext,freq_marker); }
       | STRING { $$=find_allele($1,freq_marker); free($1); }
       ;

freqlist: allele ',' {group_counter=0; freq_allele=$1;} freqlist1 {if(freq_marker && group_counter<n_gen_grps) print_scan_err("Line %d: Too few frequencies specified\n",lineno);}
        | freqlist allele ',' {group_counter=0; freq_allele=$2;} freqlist1 {if(freq_marker && group_counter<n_gen_grps) print_scan_err("Line %d: Too few frequencies specified\n",lineno);}
        ;

freqlist1: rnum {set_freq(freq_marker,$1,freq_allele);}
        | '*' {group_counter++;}
        | freqlist1 ',' rnum {set_freq(freq_marker,$3,freq_allele);}
        | freqlist1 ',' '*' {group_counter++;}
        ;

ibdlist: rnum { $$=add_ibd_list($1,0); }
        | ibdlist ',' rnum { $$=add_ibd_list($3,$1); }
        ;

array: REAL {add_real_to_num_array(($$=make_num_array(8)),$1);}
        | INTEGER {add_int_to_num_array(($$=make_num_array(8)),$1);}
        | array ',' REAL { add_real_to_num_array($1,$3); }
        | array ',' INTEGER { add_int_to_num_array($1,$3); }
        ;

rnum: REAL
    | INTEGER { $$=(double)$1; }
    ;

%%

static void set_syst_var(int com,struct num_array *r)
{
	int i;
	double x[4];
	
	if(!r) ABT_FUNC("passed zero pointer\n");
	if(com==SYST_MSCAN_PROBS) {
		if(r->ptr!=4) yyerror("System variable MSCAN_PROBS requires 4 parameters\n");
		else {
			for(i=0;i<4;i++) x[i]=(r->x[i].flag==ST_INTEGER?(double)r->x[i].data.value:r->x[i].data.rvalue);
			i=set_mscan_probs(x);
			if(i) yyerror("Illegal parameters specified for MSCAN_PROBS\n");
		}
	} else {
		if(r->ptr!=1) yyerror("Multiple values not supported for this system variable\n");
		else loki->sys.syst_var[com]=r->x[0];
	}
}

static struct num_array *make_num_array(int n)
{
	struct num_array *r;
	
	if(n<1) n=8;
	if(!(r=malloc(sizeof(struct num_array)))) ABT_FUNC(MMsg);
	r->ptr=0;
	r->size=n;
	if(!(r->x=malloc(sizeof(struct id_data)*n))) ABT_FUNC(MMsg);
	return r;
}

static void free_num_array(struct num_array *r)
{
	if(!r) ABT_FUNC("passed zero pointer\n");
	if(r->x) free(r->x);
	free(r);
}

static void add_real_to_num_array(struct num_array *r,double x)
{
	if(r->ptr>=r->size) {
		if(r->ptr>r->size) ABT_FUNC("num_array corrupted?\n");
		r->size<<=1;
		if(!(r->x=realloc(r->x,sizeof(struct id_data)*r->size))) ABT_FUNC(MMsg);
	}
	r->x[r->ptr].data.rvalue=x;
	r->x[r->ptr++].flag=ST_REAL;
}

static void add_int_to_num_array(struct num_array *r,int i)
{
	if(r->ptr>=r->size) {
		if(r->ptr>r->size) ABT_FUNC("num_array corrupted?\n");
		r->size<<=1;
		if(!(r->x=realloc(r->x,sizeof(struct id_data)*r->size))) ABT_FUNC(MMsg);
	}
	r->x[r->ptr].data.value=i;
	r->x[r->ptr++].flag=ST_INTEGER;
}

static int find_trait(struct lk_variable *lkv)
{
	int i,type,mod;
	struct Variable *var;
	
	if(!loki->models->n_models || (lkv->type!=LK_TYPE_IDVAR && lkv->type!=LK_TYPE_NONIDVAR)) {
		yyerror("Not a trait variable");
		return -1;
	}
	for(mod=0;mod<loki->models->n_models;mod++) {
		type=loki->models->models[mod].var.type;
		i=loki->models->models[mod].var.var_index;
		var=(type&ST_CONSTANT)?loki->data->id_variable+i:loki->data->nonid_variable+i;
		if(var==lkv->var.var) break;
	}
	if(mod==loki->models->n_models) {
		yyerror("Not a trait variable");
		mod= -1;
	}
	return mod;
}

static void set_analyze(char *p)
{
	int i;
	char *com[]={"AFFECTED","NULL","IBD",0};
		
	if(p) {
		i=0;
		while(com[i]) {
			if(!strcasecmp(com[i],p)) break;
			i++;
		}
		if(com[i]) loki->params.analysis|=(1<<i);
		else yyerror("Invalid parameter to analyze statement");
	}
}

static void set_ibd_mode(char *p)
{
	int i;
	char *com[]={"LOKI","MERLIN","SOLAR","SINGLEPOINT","SINGLE","SINGLE_POINT",0};
		
	if(p) {
		i=0;
		while(com[i]) {
			if(!strcasecmp(com[i],p)) break;
			i++;
		}
		if(com[i]) {
			if(i<3) loki->params.ibd_mode=i;
			else loki->params.ibd_mode |=4;
		} else yyerror("Unknown ibd mode");
	}
}

static void set_group_order(int gp)
{
	int i;
	
	if(gp<0) return;
	for(i=0;i<group_ptr;i++) if(group_order[i]==gp)	{
		yyerror("Group repeated in order statement");
		return;
	}
	if(group_ptr>=n_gen_grps) yyerror("Too many groups - internal error?");
	else group_order[group_ptr++]=gp;
}

static struct IBD_List *add_ibd_list(double x,struct IBD_List *p)
{
	if(!p) {
		if(!(p=malloc(sizeof(struct IBD_List)))) ABT_FUNC(MMsg);
		p->idx=0;
		p->size=32;
		if(!(p->pos=malloc(sizeof(double)*p->size))) ABT_FUNC(MMsg);
	}
	if(p->size==p->idx) {
		p->size*=2;
		if(!(p->pos=realloc(p->pos,sizeof(double)*p->size))) ABT_FUNC(MMsg);
	}
	p->pos[p->idx++]=x;
	return p;
}

static int check_link_name(char *name)
{
	int i=0;
	
	if(name)	{
		for(i=0;i<loki->markers->n_links;i++) if(loki->markers->linkage[i].name && !strcasecmp(name,loki->markers->linkage[i].name)) break;
	} else for(i=0;i<loki->markers->n_links;i++) if(!loki->markers->linkage[i].name) break;
	if(i==loki->markers->n_links) {
		(void)fprintf(stderr,"Warning: linkage group %s not found; command ignored\n",name);
		i= -1;
	}
	return i;
}

static struct clause_atom *new_clause_atom(char *v,int i)
{
	struct clause_atom *c;
	
	c=lk_malloc(sizeof(struct clause_atom));
	c->next=0;
	c->v=v;
	c->i=i;
	return c;
}

static struct string_list *new_string_list(char *v)
{
	struct string_list *s;
	
	s=lk_malloc(sizeof(struct string_list));
	s->next=0;
	s->v=v;
	return s;
}

static void set_pseudo_chrome(struct string_list *s,struct clause_atom *c)
{
	int i;
	static char *opt[]={"INIT","FREQ",0};

	if(c) {
		while(c) {
			i=0;
			while(opt[i]) {
				if(!strcasecmp(c->v,opt[i])) break;
				i++;
			}
			switch(i) {
			 case 0:
				loki->params.pseudo_start=c->i;
				break;
			 case 1:
				loki->params.pseudo_freq=c->i;
				break;
			 default:
				(void)fprintf(stderr,"Warning: Unknown option %s in PSEUDOCHROMOSOME command; option ignored\n",c->v);
				break;
			}
			free(c->v);
			c=c->next;
		}
		free_list(c,0);
	}
	if(s) {
		while(s) {
			i=check_link_name(s->v);
			if(i>=0) {
				loki->markers->linkage[i].type|=LINK_MIRRORED;
				loki->params.pseudo_flag=1;
			}
			free(s->v);
			s=s->next;
		}
		free_list(s,0);
	} else {
		for(i=0;i<loki->markers->n_links;i++) loki->markers->linkage[i].type|=LINK_MIRRORED;
		loki->params.pseudo_flag=1;
	}
}
  
static void set_output_gen(char *file,char *link)
{
	int i;
	struct output_gen *p;

	if(link)	{
		i=check_link_name(link);
		if(i<0) return;
		i++;
	} else i=0;
	p=loki->params.Output_Gen;
	if(!(loki->params.Output_Gen=malloc(sizeof(struct output_gen)))) ABT_FUNC(MMsg);
	loki->params.Output_Gen->next=p;
	loki->params.Output_Gen->file=file;
	loki->params.Output_Gen->link_group=i;
}

static void check_previous_list(int i)
{
	if(loki->markers->linkage[i].ibd_est_type) {
		(void)fprintf(stderr,"Warning: overwriting previous IBD settings for linkage group %s\n",loki->markers->linkage[i].name);
		if(loki->markers->linkage[i].ibd_list) {
			free(loki->markers->linkage[i].ibd_list->pos);
			free(loki->markers->linkage[i].ibd_list);
		}
		loki->markers->linkage[i].ibd_list=0;
		loki->markers->linkage[i].ibd_est_type=0;
	}
}

static void set_ibd_list(char *name,struct IBD_List *p,int type)
{
	int i=0,k;
	
	if(type==IBD_EST_GRID) {
		i=-1;
		if(p->idx<3) (void)fprintf(stderr,"Warning: too few parameters (%d) for IBD Grid (3 required); IBD request ignored\n",p->idx);
		else if(p->idx>3) (void)fprintf(stderr,"Warning: too many parameters (%d) for IBD Grid (3 required); IBD request ignored\n",p->idx);
		else if(fabs(p->pos[2])<IBD_MIN_GRID_STEP) (void)fprintf(stderr,"Warning: step size (%g) for IBD Grid < IBD_MIN_GRID_STEP (%g) in loki_ibd.h ; IBD request ignored\n",p->pos[2],IBD_MIN_GRID_STEP);
		else {
			k=1+(int)(.5+(p->pos[1]-p->pos[0])/p->pos[2]);
			if(k>IBD_MAX_GRID) (void)fprintf(stderr,"Warning: grid evaluations requested (%d) for IBD Grid > IBD_MAX_GRID (%d) in loki_ibd.h ; IBD request ignored\n",k,IBD_MAX_GRID);
			else i=0;
		}
	}
	if(!i) i=check_link_name(name);
	if(i<0) {
		free(p->pos);
		free(p);
	} else {
		check_previous_list(i);
		loki->markers->linkage[i].ibd_list=p;
		loki->markers->linkage[i].ibd_est_type=type;
	}
}

static void set_ibd_markers(char *name) 
{
	int i;
	
	i=check_link_name(name);
	if(i>=0) {
		for(;i<loki->markers->n_links;i++) {
			check_previous_list(i);
			loki->markers->linkage[i].ibd_est_type=IBD_EST_MARKERS;
			if(name) break;
		}
	}
}

static void set_tloci(int a,int b)
{
	if(a<0) {
		if(a== -1) loki->params.min_tloci=loki->params.max_tloci=b;
		else loki->params.start_tloci=b;
	} else if(a<b) {
		loki->params.min_tloci=a;
		loki->params.max_tloci=b;
	} else {
		loki->params.min_tloci=b;
		loki->params.max_tloci=a;
	}
}

static struct Marker *check_marker(struct lk_variable *lkvar)
{
	struct Marker *mk=0;
	
	if(!lkvar) return 0;
	if(lkvar->type!=LK_TYPE_MARKER)
		yyerror("Attempting to set frequency of a non-marker");
	else mk=lkvar->var.marker;
	free(lkvar);
	return mk;
}

static void set_output(struct lk_variable *lkvar)
{
	int i,j,type,mod;
	struct Variable *var;
	
	if(!lkvar) return;
	for(mod=0;mod<loki->models->n_models;mod++) {
		for(i=0;i<loki->models->models[mod].n_terms;i++) {
			type=loki->models->models[mod].term[i].vars[0].type;
			j=loki->models->models[mod].term[i].vars[0].var_index;
			if(type&ST_MARKER) {
				if(lkvar->type==LK_TYPE_MARKER && lkvar->var.marker==loki->markers->marker+j) {
					loki->models->models[mod].term[i].out_flag=1;
				}
			} else if(type&(ST_TRAITLOCUS|ST_ID|ST_SIRE|ST_DAM)) continue;
			else {
				if(type&ST_CONSTANT) var=loki->data->id_variable+j;
				else var=loki->data->nonid_variable+j;
				if((lkvar->type==LK_TYPE_IDVAR || lkvar->type==LK_TYPE_NONIDVAR) && lkvar->var.var==var) {
					loki->models->models[mod].term[i].out_flag=1;
				}
			}
		}
	}
}

static void set_map_range(char *name,double r1,double r2,int flag)
{
	int i;
	double t;
	static char *sexstr[2]={"female","male"};
	
	if(flag!= -1) loki->markers->sex_map=1;
	if(name)	{
		i=check_link_name(name);
		if(i>=0) {
			if(r2<r1) {
				t=r2;
				r2=r1;
				r1=t;
			}
			if(flag== -1) {
				loki->markers->linkage[i].r1[0]=loki->markers->linkage[i].r1[1]=r1;
				loki->markers->linkage[i].r2[0]=loki->markers->linkage[i].r2[1]=r2;
				loki->markers->linkage[i].range_set[0]=loki->markers->linkage[i].range_set[1]=1;
				message(INFO_MSG,"Map range for linkage group '%s' set to %g-%gcM\n",name,r1,r2);
			} else {
				loki->markers->linkage[i].r1[flag]=r1;
				loki->markers->linkage[i].r2[flag]=r2;
				loki->markers->linkage[i].range_set[flag]=1;
				message(INFO_MSG,"Map range (%s) for linkage group '%s' set to %g-%gcM\n",sexstr[flag],name,r1,r2);
			}
		}
	} else {
		if(flag<0) {
			loki->markers->total_maplength[X_PAT]=r1;
			loki->markers->total_maplength[X_MAT]=r2;
			message(INFO_MSG,"Total (genome) map length set to (%g,%g)cM\n",r1,r2);
		} else {
			loki->markers->total_maplength[flag]=r1;
			message(INFO_MSG,"Total (genome) %s map length set to %gcM\n",sexstr[flag],r1);
		}
	}
}

static int find_group(char *p,int gp, int flag)
{
	int i;
	char *s;
	struct Id_Recode *rec;
	
	rec=&loki->pedigree->group_recode;
	if(!rec->recode) return -1;
	if(rec->flag==ST_STRING) {
		for(i=0;i<n_gen_grps;i++) if(!(strcasecmp(p,rec->recode[i].string))) return i;
	} else {
		if(!flag) {
			gp=strtol(p,&s,10);
			if(!(*s)) for(i=0;i<n_gen_grps;i++) if(gp==rec->recode[i].value) return i;
		} else for(i=0;i<n_gen_grps;i++) if(gp==rec->recode[i].value) return i;
	}
	yyerror("Group not found\n");
	return -1;
}

static int find_allele( char *p, struct Marker *mk)
{
	int i,j;
	
	if(!mk) return -1;
	j=mk->locus.n_alleles-1;
	for(i=0;i<j;i++) if(!(strcmp(p,mk->recode[i]))) return i;
	return -1;
}

static struct lk_variable *find_var(char *p, int idx, int flag)
{
	int i,j=0;
	struct lk_variable *lkv=0;
	char *p1;
	
	if(flag) {
		if(idx<1) {
			yyerror("Illegal array index");
			free(p);
			return 0;
		}
		i=(int)strlen(p)+(int)(log((double)idx)/log(10.0)+4.00001);
		if(!(p1=malloc((size_t)i))) ABT_FUNC(MMsg);
		snprintf(p1,i,"%s(%d)",p,idx);
		free(p);
		p=p1;
	}
	for(i=0;i<loki->markers->n_markers;i++) if(!strcasecmp(loki->markers->marker[i].name,p)) {
		j=LK_TYPE_MARKER;
		break;
	}
	if(!j) for(i=0;i<loki->data->n_id_records;i++) if(!strcasecmp(loki->data->id_variable[i].name,p)) {
		j=LK_TYPE_IDVAR;
		break;
	}
	if(!j) for(i=0;i<loki->data->n_nonid_records;i++) if(!strcasecmp(loki->data->nonid_variable[i].name,p)) {
		j=LK_TYPE_NONIDVAR;
		break;
	}
	if(!j) for(i=0;i<loki->markers->n_links;i++) if(loki->markers->linkage[i].name && !strcasecmp(loki->markers->linkage[i].name,p)) {
		j=LK_TYPE_LINK;
		break;
	}
	if(j) {
		if(!(lkv=malloc(sizeof(struct lk_variable)))) ABT_FUNC(MMsg);
		lkv->type=j;
		switch(j) {
		 case LK_TYPE_MARKER:
			lkv->var.marker=loki->markers->marker+i;
			break;
		 case LK_TYPE_LINK:
			lkv->var.link=loki->markers->linkage+i;
			break;
		 case LK_TYPE_IDVAR:
			lkv->var.var=loki->data->id_variable+i;
			break;
		 case LK_TYPE_NONIDVAR:
			lkv->var.var=loki->data->nonid_variable+i;
			break;
		}
	}
	free(p);
	return lkv;
}

static void set_freq(struct Marker *mk, double freq,int allele)
{
	static int fg,fg1;
	int i=0;
	
	group_counter++;
	if(mk) {
		if(n_gen_grps>1) {
			if(!group_ptr)	{
				if(!fg) yyerror("Genetic group order not set");
				fg=1;
				return;
			}
			if(group_counter>n_gen_grps) {
				if(!fg1) yyerror("Too many frequencies specified (only 1 per genetic group)");
				fg1=1;
				return;
			}
			i=group_order[group_counter-1];
			if(i<0) return;
		}
		if(freq<0.0) {
			yyerror("Invalid (negative) frequency");
			return;
		}
		if(allele>=0 && freq==0.0)	{
			yyerror("Can not set frequency of observed allele to zero\n");
			return;
		}
		if(allele>=0) {
			mk->locus.freq[i][allele]=freq;
		} else {
			allele=mk->locus.n_alleles-1;
			mk->locus.freq[i][allele]+=freq;
		}
		mk->freq_set[i][allele]=start_flag;
		mk->count_flag[i]=c_flag;
	}
	fg1=0;
}

static void set_position(struct lk_variable *lkvar, double pos1, double pos2)
{
	if(lkvar) {
		if(lkvar->type!=LK_TYPE_MARKER) {
			yyerror("Attempting to set position of a non-marker");
			return;
		}
		lkvar->var.marker->locus.pos[X_PAT]=pos1;
		lkvar->var.marker->locus.pos[X_MAT]=pos2;
		lkvar->var.marker->pos_set=start_flag;
		free(lkvar);
	}
}

static int check_variance(const double v,const int fg)
{
	if(v<=0.0) {
		yyerror("Variance must be positive");
		return 1;
	}
	if(loki->models->n_models>1 && !fg) {
		yyerror("Must specify which model when multiple models are present");
		return 1;
	}
	if(!loki->models->n_models) {
		print_scan_warn("No model present - VARIANCE command ignored\n");
		return 1;
	}
	return 0;
}

void print_scan_err(char *fmt, ...)
{
	va_list args;
	
	va_start(args,fmt);
	(void)vfprintf(stderr,fmt,args);
	va_end(args);
	if((++scan_error_n)>=max_scan_errors) abt(__FILE__,__LINE__,"Too many errors - aborting\n");
}

void print_scan_warn(char *fmt, ...)
{
	va_list args;
	
	if(scan_warn_n<max_scan_warnings)
	{
		va_start(args,fmt);
		(void)vfprintf(stderr,fmt,args);
		va_end(args);
	}
	scan_warn_n++;
}

void yyerror(char *s)
{
     int i;
	
	print_scan_err("Line %d: %s\n%s\n",lineno,s,linebuf);
	if(scan_error_n<=max_scan_errors)
	{
		for(i=1;i<tokenpos;i++) (void)putc('-',stderr);
		(void)fputs("^\n",stderr);
	}
}

int yywrap(void)
{
	return 1;
}

int ReadParam(FILE *fptr,char *cname,struct loki *lk_struct)
{
	int i,j;
	void yy_cleanup(void);

#if YYDEBUG
	yydebug=1;
#endif
	loki=lk_struct;
	loki->params.start_tloci=-1;
	start_flag=1;
	loki->params.max_tloci=DEFAULT_MAX_TLOCI;
	loki->params.prune_option=PRUNE_LOCUS_SPECIFIC;
	loki->params.analysis=0;
	loki->params.output_type=DEFAULT_OUTPUT_TYPE;
	for(i=0;i<2;i++) loki->params.sample_freq[i]=1;
	fname_list[0]=cname;
	list_ptr=0;
	for(i=0;i<NUM_SYSTEM_VAR;i++) loki->sys.syst_var[i].flag=0;
	n_gen_grps=loki->pedigree->n_genetic_groups;
	if(!(group_order=malloc(sizeof(int)*n_gen_grps))) ABT_FUNC(MMsg);
	yyin=fptr;
	if((i=yyparse())) {
	  (void)fprintf(stderr,"Error: yyparse returned error %d\n",i);
	  scan_error_n++;
	}
	yy_cleanup();
	if(group_order) free(group_order);
	if(loki->params.start_tloci<0) loki->params.start_tloci=loki->params.min_tloci;
	else if(loki->params.start_tloci<loki->params.min_tloci || loki->params.start_tloci>loki->params.max_tloci) {
		(void)fprintf(stderr,"ReadParam(): Starting no. trait loci (%d) is outside set range (%d-%d)\n",loki->params.start_tloci,loki->params.min_tloci,loki->params.max_tloci);
		scan_error_n++;
	}
	if(loki->models->n_models>1 && loki->params.output_type<2) {
		(void)fprintf(stderr,"ReadParam(): Ouput type %d not supported with multilpe trait loci\n",loki->params.output_type);
		scan_error_n++;
	}
	if(!loki->sys.syst_var[SYST_IBD_OUTPUT].flag) {
		j=loki->params.ibd_mode;
		if(loki->params.compress_ibd) j|=COMPRESS_IBD;
		if(j) {
			loki->sys.syst_var[SYST_IBD_OUTPUT].flag=ST_INTEGER;
			loki->sys.syst_var[SYST_IBD_OUTPUT].data.value=j;
		}
	}
	if(scan_error_n) return 1;
	return 0;
}
