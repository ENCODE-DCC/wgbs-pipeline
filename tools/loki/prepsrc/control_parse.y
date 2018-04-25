%{
/****************************************************************************
 *                                                                          *
 *     Loki - Programs for genetic analysis of complex traits using MCMC    *
 *                                                                          *
 *             Simon Heath - University of Washington                       *
 *                                                                          *
 *                       March 1997                                         *
 *                                                                          *
 * control_parse.y:                                                         *
 *                                                                          *
 * yacc source for control file parser.                                     *
 *                                                                          *
 * Copyright (C) Simon C. Heath 1997, 2000, 2002                            *
 * This is free software.  You can distribute it and/or modify it           *
 * under the terms of the Modified BSD license, see the file COPYING        *
 *                                                                          *
 ****************************************************************************/

#include <stdarg.h>
#include <stdlib.h>
#include <string.h>
#ifdef USE_DMALLOC
#include <dmalloc.h>
#endif
#include <stdio.h>
#include <math.h>
#include <ctype.h>
#include <assert.h>

#include "config.h"
#include "utils.h"
#include "lk_malloc.h"
	
#ifndef YYDEBUG
#define YYDEBUG 0
#endif
%}

%union 
{
	char *string;
	struct bin_node *var;
	int value;
	double rvalue;
	struct format_clause *format_clause;
	struct fformat *fformat;
	struct format_atom *f_atom;
	struct model_list *model_list;
	struct var_list *var_list;
	struct var_element *element;
	struct express *express;
}

%token FILEC MARKER LOCUS TRAIT RANDOM PEDIGREE LOG SUPER
%token MODEL FILTER LINK MISSING FACTOR BREAK DOLOOP WHILE
%token USE WHERE ORSYMBOL ANDSYMBOL NEQSYMBOL LEQSYMBOL GEQSYMBOL
%token NOTSYMBOL LOGICAL SHELL ARRAY PRINTEXP INCLUDE RAWOUTPUT
%token LOOP_CLAUSE_START LOOP_CLAUSE_END CONSTANT MULTIPLE RSFORMAT FSFORMAT
%token SKIPFORMAT GSFORMAT CENSORED GROUP SET GENDER AFFECTED OUTPUT 
%token ERRORDIR LAUROUTPUT UNAFFECTED POSITION FREQUENCY PROBAND

%token <string> STRING
%token <var> VARIABLE ASSIGN ARRAY_VAR
%token <value> INTEGER SYSTEM_VAR
%token <rvalue> REAL

%type <string> complex_string complex_string1 filename_string
%type <format_clause> formatlist
%type <fformat> fformatlist fformat
%type <f_atom> format
%type <model_list> modellist
%type <var_list> fsingle_vlist filevarlist varlist interactionlist single_vlist simple_varlist loop_clause
%type <element> single_element assignment
%type <express> expression condition
%type <value> linkcom1
%type <var> variable

%left ORSYMBOL
%left ANDSYMBOL
%left '=' NEQSYMBOL
%left '<' '>' LEQSYMBOL GEQSYMBOL
%left '+' '-'
%left '*' '.' '/'
%left NOTSYMBOL UMINUS

%{
#include "scan.h"
#include "scanner.h"

static struct format_atom *make_f_atom(int,int);
static struct format_clause *add_f_atom(struct format_clause *,struct format_atom *);
static struct format_clause *add_f_list(struct format_clause *,struct format_clause *,int);
static struct format *setup_format(struct format_clause *);
static struct bin_node *create_var(char *);
static struct model_list *add_to_model(struct model_list *,struct var_list *);
static struct var_list *add_to_var_list(struct var_list *,struct bin_node *,struct express *);
static struct var_list *add_var_lists(struct var_list *,struct var_list *);
static struct var_element *get_element(struct bin_node *,struct express *);
static struct var_element *assign_var(struct bin_node *,struct express *,struct express *);
static struct express *alloc_express(void);
static struct express *do_express_op(struct express *,struct express *,int);
static struct express *do_logical_op(struct express *,struct express *,int);
static struct fformat *add_fformat(struct fformat *,struct fformat *);
static struct fformat *create_fformat(void *,int);
static void begin_looping(struct var_element *, struct express *, struct express *);
static void free_vlist(struct var_list *);
static void do_ped_com(struct var_list *);
static void set_position(struct var_element *,struct express *,struct express *,struct express *);
static void add_restriction(struct var_list *);
static void add_censored(struct var_element *,const int);
static void do_file_com(char *,struct format_clause *,struct fformat *,struct var_list *);
static void set_locus_array(struct bin_node *);
static void set_locus_element(struct var_element *);
static void set_slocus_element(struct var_element *,struct var_list *);
static void set_haplo_element(struct var_element *,struct var_element *);
static void do_link_com(char *,int, struct var_list *);
static void do_missing_com(struct express *,struct var_list *,char *);
static void change_type(int,struct var_list *);
static void do_model_com(struct model_list *,struct bin_node *,struct express *);
static void add_operation(void *,int,int);
static void set_array_var(struct scan_data *,struct express *);
static void check_element_add_op(struct var_element *);
static void enter_loop(void);
static void start_loopclause(void);
static void do_while_com(struct express *);
static void print_exp(struct express *);
static void new_command(void);
static void set_sex(struct var_element *,struct express *,struct express *);
static void set_group(struct var_element *);
static int count_var_list(struct var_list *),shell_flag;
static struct format_atom *f_atom_list;
static struct var_element *pedlist[4];
static int f_atom_n,f_atom_size=32,pedflag;
struct operation *Affected,*Unaffected,*Proband;
struct bin_node *root_var;
struct InFile *Infiles;
struct Link *links;
struct Miss *Miss;
struct Restrict *Restrictions;
struct Censor *Censored;
struct model *Models;
static struct operation *Op_List;
struct Marker *markers,*traitlocus;
static struct var_element *hap_list[2];
struct express *sex_exp[2];
struct sex_def *sex_def;
struct var_element *group_elem;
static char *string_copy(char *s1,char *s2);
static char *LogFile;
static struct marker_info *m_info;
	
int scan_error,scan_error_n,scan_warn_n;
int max_scan_errors=30,max_scan_warnings=30,n_markers,iflag,file_skip;
char *Filter,*ErrorDir,*rsformat,*fsformat,*gsformat,*OutputFile,*OutputRawFile,*OutputLaurFile;
int loop_level,loop_ptr[MAX_LOOP],loop_stat[MAX_LOOP],loop_record,loop_stack_size=256;
int loop_main_ptr,in_loopclause,loop_clause_end,loop_clause_step,loop_clause_ptr;
int syst_var[NUM_PREP_SYSTEM_VAR];
int family_id;
struct var_element *loop_clause_element;
struct token_store *loop_stack;
	  
%}

%%

comfile: command { new_command(); }
       | error { new_command(); }
		 | comfile BREAK
       | comfile BREAK command { new_command(); }
       | comfile BREAK error { new_command(); }
       ;

command: filecommand
       | filtercommand
       | includecommand
	    | logcommand
	    | pedcommand
	    | locicommand
	    | changetypecommand
	    | linkcommand
	    | modelcommand
	    | missingcommand
       | usecommand
       | wherecommand
       | assigncommand
       | arraycommand
       | printcommand
       | docommand
       | whilecommand
       | censorcommand
       | affectedcommand
       | groupcommand
       | setcommand
       | defformatcommand
       | sexcommand
       | outputcommand
       | errordircommand
       | positioncommand
	    ;

assigncommand: assignment {}
       ;

assignment: ASSIGN '=' expression { $$=assign_var($1,0,$3); if($3) free($3);}
       | ASSIGN '(' expression ')' '=' expression { $$=assign_var($1,$3,$6); if($3) free($3); if($6) free($6);}
       ;

expression: REAL { $$=alloc_express(); $$->type=ST_REAL; $$->arg.rvalue=$1; }
       | INTEGER { $$=alloc_express(); $$->type=ST_INTEGER; $$->arg.value=$1; }
       | STRING { $$=alloc_express(); $$->type=ST_STRING; $$->arg.string=$1; }
       | single_element { $$=alloc_express(); $$->type=$1->type; if($1->type==ST_STRING) $$->arg.string=string_copy(0,$1->arg.string);
		  else $$->arg=$1->arg; }
       | expression '+' expression { $$=do_express_op($1,$3,'+'); }
       | expression '-' expression { $$=do_express_op($1,$3,'-'); }
       | expression '*' expression { $$=do_express_op($1,$3,'*'); }
       | expression '/' expression { $$=do_express_op($1,$3,'/'); }
       | '-' expression %prec UMINUS { $$=do_express_op($2,0,'-'); }
       | '(' expression ')' { $$=$2; }
       ;

setcommand: SET SYSTEM_VAR INTEGER { syst_var[$2]=$3; }
       | SET variable INTEGER { yyerror("Unrecognized system variable"); }
       ;

docommand: DOLOOP { enter_loop(); }
       ;

includecommand: INCLUDE {iflag=1;} complex_string {include_control_file($3);}
       ;

censorcommand: CENSORED single_element WHERE '(' res_condition ')' { add_censored($2,1); at_use=0;}
       ;

affectedcommand: AFFECTED WHERE '(' res_condition ')' { add_censored(0,0); at_use=0;}
		|  UNAFFECTED WHERE '(' res_condition ')' { add_censored(0,2); at_use=0;}
		|  PROBAND WHERE '(' res_condition ')' { add_censored(0,3); at_use=0;}
       ;

sexcommand: GENDER single_element expression ',' expression { set_sex($2,$3,$5); }
       ;

whilecommand: WHILE '(' condition ')' { do_while_com($3); if($3) free($3);}
       ;

positioncommand: POSITION single_element expression { set_position($2,$3,0,0); }
       | POSITION single_element expression ',' expression { set_position($2,0,$3,$5); }
       | POSITION single_element expression ',' expression ',' expression { set_position($2,$3,$5,$7); }
       | POSITION single_element ',' expression { set_position($2,0,0,$4); }
       | POSITION single_element ',' ',' expression { set_position($2,0,0,$5); }
       ;

condition: expression
          | condition '=' condition {$$=do_logical_op($1,$3,'=');}
          | condition NEQSYMBOL condition {$$=do_logical_op($1,$3,NEQSYMBOL);}
          | condition LEQSYMBOL condition {$$=do_logical_op($1,$3,LEQSYMBOL);}
          | condition GEQSYMBOL condition {$$=do_logical_op($1,$3,GEQSYMBOL);}
          | condition '<' condition {$$=do_logical_op($1,$3,'<');}
          | condition '>' condition {$$=do_logical_op($1,$3,'>');}
          | condition ORSYMBOL condition {$$=do_logical_op($1,$3,ORSYMBOL);}
          | condition ANDSYMBOL condition {$$=do_logical_op($1,$3,ANDSYMBOL);}
          | NOTSYMBOL condition {$$=do_logical_op($2,0,NOTSYMBOL);}
          ;

printcommand: PRINTEXP printlist { (void)fputc('\n',stdout); }
       ;

printlist: expression { print_exp($1); if($1) free($1); }
       | printlist ',' expression { print_exp($3); if($3) free($3); }
       ;

defformatcommand: RSFORMAT complex_string {if(rsformat) free(rsformat); rsformat=$2;}
       | FSFORMAT complex_string {if(fsformat) free(fsformat); fsformat=$2;}
       | GSFORMAT complex_string {if(gsformat) free(gsformat); gsformat=$2;}
       | SKIPFORMAT INTEGER {file_skip=$2;}
       ;

arraycommand: ARRAY arraylist;

arraylist: VARIABLE '(' expression ')' {set_array_var($1->data,$3); if($3) free($3); }
       | arraylist ',' VARIABLE '(' expression ')' {set_array_var($3->data,$5); if($5) free($5); }
       ;

usecommand: USE varlist WHERE '(' res_condition ')' {add_restriction($2);at_use=0;}
          | USE WHERE '(' res_condition ')' {add_restriction(0);at_use=0;}
       ;
		 
wherecommand: WHERE '(' res_condition ')' USE varlist {add_restriction($6);at_use=0;}
            | WHERE '(' res_condition ')' {add_restriction(0);at_use=0;}
        ;
		  
res_condition: INTEGER {add_operation(&($1),INTEGER,0);}
          | REAL {add_operation(&($1),REAL,0);}
          | STRING {add_operation($1,STRING,0);}
          | single_element {if($1) check_element_add_op($1);}
          | res_condition '+' res_condition {add_operation(0,0,'+');}
          | res_condition '-' res_condition {add_operation(0,0,'-');}
          | res_condition '*' res_condition {add_operation(0,0,'*');}
          | res_condition '/' res_condition {add_operation(0,0,'/');}
          | '(' res_condition ')'
          | '-' res_condition %prec UMINUS {add_operation(0,0,UMINUS);}
          | res_condition '=' res_condition {add_operation(0,0,'=');}
          | res_condition NEQSYMBOL res_condition {add_operation(0,0,NEQSYMBOL);}
          | res_condition LEQSYMBOL res_condition {add_operation(0,0,LEQSYMBOL);}
          | res_condition GEQSYMBOL res_condition {add_operation(0,0,GEQSYMBOL);}
          | res_condition '<' res_condition {add_operation(0,0,'<');}
          | res_condition '>' res_condition {add_operation(0,0,'>');}
          | res_condition ORSYMBOL res_condition {add_operation(0,0,ORSYMBOL);}
          | res_condition ANDSYMBOL res_condition {add_operation(0,0,ANDSYMBOL);}
          | NOTSYMBOL res_condition {add_operation(0,0,NOTSYMBOL);}
          ;

filecommand: FILEC filename_string ',' filevarlist { do_file_com($2,0,0,$4); }
       | FILEC '[' formatlist ']' filename_string ',' filevarlist { do_file_com($5,$3,0,$7); }
       | FILEC '[' fformatlist ']' filename_string ',' filevarlist { do_file_com($5,0,$3,$7); }
	  ;

filename_string: complex_string { $$=$1; }
       | SHELL '(' complex_string ')' { $$=$3; shell_flag=1; }
       ;

formatlist: format {$$=add_f_atom(0,$1); }
       | formatlist ',' format {$$=add_f_atom($1,$3); }
	  | INTEGER '(' formatlist ')' {$$=add_f_list(0,$3,$1); }
	  | formatlist ',' INTEGER '(' formatlist ')' {$$=add_f_list($1,$5,$3); }
       ;

fformatlist: fformat
	  | fformatlist ',' fformat {$$=add_fformat($1,$3); }
	  | fformatlist ';' fformat {$$=add_fformat($1,$3); }
	  ;

fformat: FSFORMAT complex_string {$$=create_fformat($2,2); }
       | RSFORMAT complex_string {$$=create_fformat($2,1); }
       | GSFORMAT complex_string {$$=create_fformat($2,4); }
       | SKIPFORMAT INTEGER {$$=create_fformat(&$2,3); }
       ;

format: INTEGER {$$=make_f_atom($1,0);}
       | INTEGER 'x' {$$=make_f_atom($1,1);}
       | 'x' {$$=make_f_atom(1,1);}
       | INTEGER error {$$=make_f_atom(0,1); scan_error|=FORMAT_ERR; }
       | error {$$=make_f_atom(0,1); scan_error|=FORMAT_ERR; }
	  ;

logcommand: LOG complex_string
	  { if(LogFile) free(LogFile); LogFile=$2; }
     ;
	  
outputcommand: OUTPUT complex_string { if(OutputFile) free(OutputFile); OutputFile=$2; }
     | LAUROUTPUT complex_string { if(OutputLaurFile) free(OutputLaurFile); OutputLaurFile=$2; }
     | RAWOUTPUT complex_string { if(OutputRawFile) free(OutputRawFile); OutputRawFile=$2; }
;

errordircommand: ERRORDIR complex_string 
	  { if(ErrorDir) free(ErrorDir); ErrorDir=$2;}
     ;
	  
missingcommand: MISSING expression { do_missing_com($2,0,0); free($2); }
              | MISSING expression varlist { do_missing_com($2,$3,0); free($2); }
              | MISSING expression ',' varlist { do_missing_com($2,$4,0); free($2); }
              | MISSING '[' complex_string ']' expression { do_missing_com($5,0,$3); free($5); }
       ;
		 
pedcommand: PEDIGREE varlist { do_ped_com($2); }
          ;

locicommand: MARKER LOCUS locuslist
           | TRAIT LOCUS varlist { change_type(ST_TRAITLOCUS,$3); free_vlist($3);}
           | SUPER LOCUS slocuslist
          ;
			 
changetypecommand: CONSTANT varlist { change_type(ST_CONSTANT,$2); free_vlist($2);}
       | MULTIPLE varlist { change_type(ST_MULTIPLE,$2); free_vlist($2);}
       | RANDOM varlist { change_type(ST_RANDOM|ST_FACTOR,$2); free_vlist($2);}
       | FACTOR varlist { change_type(ST_FACTOR,$2); free_vlist($2);}
       | REAL varlist {change_type(ST_REALTYPE,$2); free_vlist($2);}
       | INTEGER varlist {change_type(ST_INTTYPE,$2); free_vlist($2);}
       ;
		 
linkcom1: LINK { $$=LINK_AUTO; }
       | LINK '[' 'x' ']' { $$=LINK_X; }
       | LINK '[' 'y' ']' { $$=LINK_Y; }
	    ;

linkcommand: linkcom1 varlist { do_link_com(0,$1,$2); }
         | linkcom1 complex_string1 ',' varlist { do_link_com($2,$1,$4); }
         | linkcom1 complex_string1 { do_link_com($2,$1,0); }
         ;

filtercommand: FILTER complex_string  {
	if(Filter) {
	  print_scan_warn("Line %d: Warning - Filter defined twice\n",lineno);
	  free(Filter);
   }
	Filter=$2; }
       ;
	  
modelcommand: MODEL variable '=' modellist {do_model_com($4,$2,0);}
       | MODEL ARRAY_VAR '(' expression ')' '=' modellist {do_model_com($7,$2,$4); if($4) free($4); }
       ;

groupcommand: GROUP single_element {set_group($2);}
       ;

modellist: modellist '+' single_vlist { $$=add_to_model($1,$3); }
       | modellist '+' interactionlist { $$=add_to_model($1,$3); }
	  | interactionlist { $$=add_to_model(0,$1); }
	  | single_vlist { $$=add_to_model(0,$1); }
       ;

interactionlist: single_vlist '*' single_vlist { $$=add_var_lists($1,$3); }
               | single_vlist '.' single_vlist { $$=add_var_lists($1,$3); }
               | interactionlist '*' single_vlist { $$=add_var_lists($1,$3); }
               | interactionlist '.' single_vlist { $$=add_var_lists($1,$3); }
               ;

variable: VARIABLE
       | GENDER { $$=create_var("SEX"); }
       ;

single_vlist: variable { $$=add_to_var_list(0,$1,0); }
       | ARRAY_VAR { $$=add_to_var_list(0,$1,0); }
       | ARRAY_VAR '(' expression ')' { $$=add_to_var_list(0,$1,$3); if($3) free($3); }
       ;

single_element: variable { $$=get_element($1,0); }
       | ARRAY_VAR '(' expression ')' { $$=get_element($1,$3); if($3) free($3); }
       ;

locuslist: locus | locuslist ',' locus;

locus: single_element lociclause { if($1) set_locus_element($1); }
		 | single_element { if($1) set_locus_element($1); }
		 | ARRAY_VAR { if($1) set_locus_array($1); }
       ;

lociclause: '[' single_element ']' { if($2) set_haplo_element($2,0); }
         | '[' single_element ',' single_element ']' { if($2) set_haplo_element($2,$4); }
         | '[' error ']'
	    ;

slocuslist: slocus | slocuslist ',' locus;

slocus: single_element '[' varlist ']' { if($1) set_slocus_element($1,$3); }
        | single_element { yyerror ("No marker list for Super Locus\n"); }
        | single_element '[' error ']'
        ;

simple_varlist: single_vlist
       | simple_varlist ',' single_vlist { $$=add_var_lists($1,$3); }
	  ;

open_bracket: '(' { start_loopclause(); }
            ;
	
loop_clause1: open_bracket simple_varlist ',' BREAK assignment ',' expression ')' { free_vlist($2); begin_looping($5,$7,0); }
           | open_bracket simple_varlist ',' BREAK assignment ',' expression ',' expression ')' { free_vlist($2); begin_looping($5,$7,$9); }
           ;

loop_clause: loop_clause1 LOOP_CLAUSE_START simple_varlist LOOP_CLAUSE_END { $$=$3; in_loopclause=0; }
           | loop_clause1 { $$=0; in_loopclause=0; }
           ;

fsingle_vlist: { $$=add_to_var_list(0,0,0); }
       | variable { $$=add_to_var_list(0,$1,0); }
       | ARRAY_VAR { $$=add_to_var_list(0,$1,0); }
       | ARRAY_VAR '(' expression ')' { $$=add_to_var_list(0,$1,$3); if($3) free($3); }
       ;

filevarlist: fsingle_vlist
       | loop_clause
       | filevarlist ',' fsingle_vlist { $$=add_var_lists($1,$3); }
       | filevarlist ',' loop_clause { $$=add_var_lists($1,$3); }
	  ;

varlist: single_vlist
       | loop_clause
       | varlist ',' single_vlist { $$=add_var_lists($1,$3); }
       | varlist ',' loop_clause { $$=add_var_lists($1,$3); }
	  ;

complex_string1: STRING { $$ = $1; }
       | complex_string1 '+' STRING { $$ = string_copy($1,$3); free($3); }
       | single_element '+' complex_string1 { if($1 && ($1->type&ST_STRING)) $$ = string_copy($3,$1->arg.string); else $$=$3; }
       ;

complex_string: STRING { $$ = $1; }
       | single_element { if($1 && ($1->type&ST_STRING)) $$ = string_copy(0,$1->arg.string); else $$=0; }
       | complex_string '+' STRING { $$ = string_copy($1,$3); free($3); }
       | complex_string '+' single_element { if($3 && ($3->type&ST_STRING)) $$ = string_copy($1,$3->arg.string); else $$=$1; }
       ;
%%

static void enter_loop(void)
{	
	if(loop_level<MAX_LOOP)	{
		loop_stat[loop_level]=loop_record;
		loop_ptr[loop_level]=loop_main_ptr;
		loop_level++;
		if(!loop_record) loop_record=1;
	} else yyerror("Too many nested loops\n");
}

static void start_loopclause(void)
{	
	if(loop_level<MAX_LOOP)	{
		loop_stat[loop_level]=loop_record;
		loop_ptr[loop_level]=loop_main_ptr;
		in_loopclause=1;
		loop_level++;
		if(!loop_record) loop_record=1;
	} else yyerror("Too many nested loops\n");
}

static void begin_looping(struct var_element *element, struct express *exp1, struct express *exp_2)
{	
	int er=0,i;
	
	if(element && exp1) {
		if(exp1->type==ST_INTEGER) loop_clause_end=(int)exp1->arg.value;
		else er=1;
		if(!er && exp_2) {
			if(exp_2->type==ST_INTEGER) loop_clause_step=(int)exp_2->arg.value;
			else er=1;
		} else loop_clause_step=1;
		if(element->type&ST_INTEGER) {
			i=(int)element->arg.value;
			if(loop_clause_step<0) {
				if(i<loop_clause_end) er= -1;
			} else if(i>loop_clause_end) er= -1;
		} else er=1;
	} else er=2;
	if(er) {
		switch(er) {
		 case 1:
			yyerror("Loop variable not integer type\n");
			break;
		 case 2:
			yyerror("Syntax error\n");
			break;
		}
		loop_record=loop_stat[--loop_level];
		if(!loop_record) loop_main_ptr=loop_level?loop_ptr[loop_level-1]:0;
		in_loopclause=0;
	} else {
		in_loopclause= -1;
		loop_record= -1;
		loop_clause_ptr=loop_main_ptr;
		loop_clause_element=element;
		loop_main_ptr=loop_ptr[loop_level-1];
	}
	if(exp1) free(exp1);
	if(exp_2) free(exp_2);
}

static int check_pos(struct express *e,double *x)
{
	int er=0;
	
	if(e->type==ST_INTEGER) *x=(double)e->arg.value;
	else if(e->type==ST_REAL) *x=e->arg.rvalue;
	else {
		yyerror1("Non-numeric arguments to position command\n");
		er=1;
	}
	return er;
}

static void set_position(struct var_element *elem,struct express *sex_avg,struct express *male,struct express *female)
{
	int i,set[3],er=0;
	double x[3];
	struct marker_info *mi;
	
	/* First check if positions are valid */
	if(sex_avg) {
		if(check_pos(sex_avg,x+2)) er=1;
		else set[2]=1;
	} else set[2]=0;
	if(!er && female) {
		if(check_pos(male,x)) er=1;
		else set[0]=1;
	} else set[0]=0;
	if(!er && male) {
		if(check_pos(female,x+1)) er=1;
		else set[1]=1;
	} else set[1]=0;
	/* If OK, then check if this element has already got an info entry */
	if(!er && (set[0] || set[1] || set[2])) {
		mi=m_info;
		while(mi) {
			if(mi->element==elem) break;
			mi=mi->next;
		}
		if(mi) {
			if((mi->pos_set[0] && set[0]) || (mi->pos_set[1] && set[1]) || (mi->pos_set[2] && set[2])) {
				yyerror1("Error: duplicate positions set for marker\n");
				er=1;
			} else {
				for(i=0;i<3;i++) mi->pos_set[i]|=set[i];
			}
		} else {
			mi=lk_malloc(sizeof(struct marker_info));
			mi->next=m_info;
			m_info=mi;
			for(i=0;i<3;i++) mi->pos_set[i]=set[i];
			mi->element=elem;
		}
		if(!er) for(i=0;i<3;i++) if(set[i]) mi->pos[i]=x[i];
	}
}

static int if_true(struct express *express)
{
	int l=0;
	
	if(express)	{
		switch(express->type) {
		 case ST_INTEGER:
			l=(express->arg.value!=0);
			break;
		 case ST_REAL:
			l=(express->arg.rvalue!=0.0);
			break;
		 case ST_STRING:
			if(express->arg.string && express->arg.string[0]) l=1;
		}
	}
	return l;
}

static void do_while_com(struct express *express)
{
	int i;
	
	if(loop_level)	{
		if(!scan_error_n && if_true(express)) {
			loop_record= -1;
			loop_main_ptr=loop_ptr[loop_level-1];
		} else {
			loop_record=loop_stat[--loop_level];
			if(!loop_record) {
				i=loop_main_ptr-1;
				loop_main_ptr=loop_level?loop_ptr[loop_level-1]:0;
				for(;i>=loop_main_ptr;i--)	{
					if(loop_stack[i].token==STRING) free(loop_stack[i].yylval.string);
				}
			}
		}
	} else yyerror("WHILE outside of do loop\n");
}

static void print_exp(struct express *express)
{
	if(!express) return;
	switch(express->type) {
	 case ST_STRING:
		(void)fputs(express->arg.string,stdout);
		free(express->arg.string);
		break;
	 case ST_INTEGER:
		(void)printf("%ld",express->arg.value);
		break;
	 case ST_REAL:
		(void)printf("%g",express->arg.rvalue);
		break;
	}
}

static void set_sex(struct var_element *elem,struct express *exp_1,struct express *exp_2)
{
	struct sex_def *se;
	
	if(!exp_1 || !exp_2) yyerror1("Null arguments to sex command\n");
	else {
		if(exp_1->type != exp_2->type) yyerror1("Arguments to sex command of different type\n");
		else if(exp_1->type!=ST_INTEGER && exp_1->type!=ST_STRING) yyerror1("Arguments to sex command of invalid type\n");
		else {
			se=lk_malloc(sizeof(struct sex_def));
			se->sex_exp[0]=exp_1;
			se->sex_exp[1]=exp_2;
			se->sex_elem=elem;
			elem->type|=(ST_SEX|ST_FACTOR|ST_CONSTANT|ST_DATA);
			if(exp_1->type==ST_INTEGER) elem->type|=ST_INTTYPE;
			se->next=sex_def;
			sex_def=se;
		}
	}
}

static void set_group(struct var_element *elem)
{
	if(group_elem) yyerror1("Error: Multiple group commands");
	else {
		group_elem=elem;
		group_elem->type|=(ST_GROUP|ST_FACTOR|ST_CONSTANT);
	}
}

static void do_ped_com(struct var_list *vlist)
{
	int i,j,n;
	struct var_list *vlist1;
	struct scan_data *sd;
	
	if(pedflag) {
		yyerror1("Error: Multiple pedigree commands");
		scan_error|=PED_ERR;
	}
	pedflag=1;
	n=0;
	while(vlist) {
		sd=vlist->var->data;
		if(sd->vtype&ST_ARRAY) {
			if(vlist->index) {
				if(n<4) pedlist[n]=sd->element+vlist->index-1;
				n++;
			} else for(i=0;i<sd->n_elements;i++)	{
				if(n<4) pedlist[n]=sd->element+i;
				n++;
			}
		} else {
			if(n<4) pedlist[n]=sd->element;
			n++;
		}
		vlist1=vlist->next;
		free(vlist);
		vlist=vlist1;
	}
	if(n!=3 && n!=4) {
		yyerror1("Error: Wrong no. variables for pedigree command (3 or 4 required)");
		scan_error|=PED_ERR;
	} else {
		for(i=1;i<n;i++) {
			for(j=0;j<i;j++) if(pedlist[i]==pedlist[j]) {
				yyerror1("Error: Repeated variables in pedigree command");
				scan_error|=PED_ERR;
			}
		}
		if(n==4) {
			pedlist[0]->type|=(ST_FAMILY|ST_FACTOR|ST_CONSTANT);
			family_id=1;
			i=1;
		} else i=0;
		pedlist[i++]->type|=(ST_ID|ST_FACTOR|ST_CONSTANT);
		pedlist[i++]->type|=(ST_SIRE|ST_FACTOR|ST_CONSTANT);
		pedlist[i++]->type|=(ST_DAM|ST_FACTOR|ST_CONSTANT);
	}
}

static struct express *alloc_express(void)
{
	struct express *e;
	
	e=lk_malloc(sizeof(struct express));
	return e;
}

static struct express *do_logical_op(struct express *ex1,struct express *ex2,int op)
{
	int i=0,l=0;
	double rv1,rv2;
	char *s1,*s2;
	
	s1=s2=0;
	if(ex1->type&ST_STRING) {
		s1=ex1->arg.string;
		i++;
	}
	if(ex2 && ex2->type&ST_STRING) {
		s2=ex2->arg.string;
		i++;
	}
	if(i==2)	{
		switch(op) {
		 case '=':
			l=mystrcmp(s1,s2)?0:1;
			break;
		 case NEQSYMBOL:
			l=mystrcmp(s1,s2)?1:0;
			break;
		 case '<':
			l=mystrcmp(s1,s2)<0?1:0;
			break;
		 case '>':
			l=mystrcmp(s1,s2)>0?1:0;
			break;
		 case LEQSYMBOL:
			l=mystrcmp(s1,s2)>=0?1:0;
			break;
		 case GEQSYMBOL:
			l=mystrcmp(s1,s2)>=0?1:0;
			break;
		 case ORSYMBOL:
			l=(s1 || s2);
			break;
		 case ANDSYMBOL:
			l=(s1 && s2);
			break;
		 default:
			ABT_FUNC("Internal error - invalid string op\n");
		}
	} else if(i && !ex2)	{
		assert(op==NOTSYMBOL);
		l=s1?0:1;
	} else if(i) {
		switch(op) {
		 case ORSYMBOL:
			if(ex1->type&ST_STRING) l=(s1 || ex2->arg.value);
			else l=(s2 || ex1->arg.value);
			break;
		 case ANDSYMBOL:
			if(ex1->type&ST_STRING) l=(s1 && ex2->arg.value);
			else l=(s2 && ex1->arg.value);
			break;
		 default:
			ABT_FUNC("Internal error - invalid string op\n");
		}
	} else {
		rv1=rv2=0.0;
		if(ex1->type&ST_INTEGER) rv1=(double)ex1->arg.value;
		else if(ex1->type&ST_REAL) rv1=ex1->arg.rvalue;
		if(ex2) {
			if(ex2->type&ST_INTEGER) rv2=(double)ex2->arg.value;
			else if(ex2->type&ST_REAL) rv2=ex2->arg.rvalue;
		}
		switch(op) {
		 case '=':
			l=(rv1==rv2);
			break;
		 case NEQSYMBOL:
			l=(rv1!=rv2);
			break;
		 case '<':
			l=(rv1<rv2);
			break;
		 case '>':
			l=(rv1>rv2);
			break;
		 case LEQSYMBOL:
			l=(rv1<=rv2);
			break;
		 case GEQSYMBOL:
			l=(rv1>=rv2);
			break;
		 case ORSYMBOL:
			l=(rv1 || rv2);
			break;
		 case ANDSYMBOL:
			l=(rv1 && rv2);
			break;
		 case NOTSYMBOL:
			l=(rv1==0.0);
			break;
		 default:
			ABT_FUNC("Internal error - invalid op\n");
		}
	}
	if(ex2) free(ex2);
	ex1->type=ST_INTEGER;
	ex1->arg.value=l;
	return ex1;
}

static struct express *do_express_op(struct express *ex1,struct express *ex2,int op)
{
	double rv1,rv2;
	int i;
	
	if(ex1->type&ST_STRING)	{
		if(ex2 && ex2->type&ST_STRING) {
			if(op!='+') yyerror("Illegal string operation\n");
			else {
				ex1->arg.string=string_copy(ex1->arg.string,ex2->arg.string);
				free(ex2->arg.string);
			}
		} else yyerror("Can't mix numeric and string expressions\n");
	} else if(ex2 && ex2->type&ST_STRING) yyerror("Can't mix numeric and string expressions\n");
	else {
		rv1=rv2=0.0;
		if(ex1->type&ST_INTEGER) rv1=(double)ex1->arg.value;
		else if(ex1->type&ST_REAL) rv1=ex1->arg.rvalue;
		if(ex2) {
			if(ex2->type&ST_INTEGER) rv2=(double)ex2->arg.value;
			else if(ex2->type&ST_REAL) rv2=ex2->arg.rvalue;
		}
		switch(op) {
		 case '+': 
			rv1+=rv2;
			break;
		 case '-':
			if(ex2) rv1-=rv2;
			else rv1= -rv1;
			break;
		 case '*':
			rv1*=rv2;
			break;
		 case '/':
			if(rv2==0.0) {
				yyerror("Divide by zero error\n");
				rv1=0.0;
			} else rv1/=rv2;
			break;
		}
		i=(int)rv1;
		if((double)i==rv1) {
			ex1->type=ST_INTEGER;
			ex1->arg.value=i;
		} else {
			ex1->type=ST_REAL;
			ex1->arg.rvalue=rv1;
		}
	}
	if(ex2) free(ex2);
	return ex1;
}

static void check_element_add_op(struct var_element *element)
{
	switch(element->type&(ST_REAL|ST_INTEGER|ST_STRING)) {
	 case ST_STRING:
		add_operation(string_copy(0,element->arg.string),STRING,0);
		break;
	 case ST_INTEGER:
		add_operation(&element->arg,INTEGER,0);
		break;
	 case ST_REAL:
		add_operation(&element->arg,REAL,0);
		break;
	 case 0:
		add_operation(element,VARIABLE,0);
		break;
	 default:
		ABT_FUNC("Internal error - illegal element type\n");
	}
}

static int check_index(struct scan_data *sd,struct express *express)
{
	int i;

	if(express->type!=ST_INTEGER) {
		if(in_loopclause<=0) yyerror("Non-integral expression for array index");
	} else if(sd->vtype&ST_ARRAY) {
		i=(int)express->arg.value;
		if(i<1 || i>sd->n_elements) {
			if(in_loopclause<=0) yyerror("Array index out of bounds");
		} else return i;
	} else yyerror("Not an array");
	return 0;
}

static struct var_element *get_element(struct bin_node *node,struct express *express)
{
	int i;
	struct scan_data *sd;
	
	sd=node->data;
	if(express) {
		if(!(i=check_index(sd,express))) return 0;
		return sd->element+i-1;
	} else {
		if(sd->vtype&ST_ARRAY)	{
			yyerror("Illegal reference to array");
			return 0;
		}
	}
	return sd->element;
}

static void set_array_var(struct scan_data *sd,struct express *express)
{
	int i;

	if(express->type!=ST_INTEGER) yyerror("Non-integral expression for array size");
	else if((i=(int)express->arg.value)<1) yyerror("Illegal array size");
	else if(sd->vtype) yyerror("Can't redefine variable");
	else {
		sd->vtype|=ST_ARRAY;
		sd->n_elements=i;
		free(sd->element);
		sd->element=lk_calloc((size_t)sd->n_elements,sizeof(struct var_element));
	}
}

static int count_var_list(struct var_list *vlist)
{
	int i=0;
	struct scan_data *sd;
	
	while(vlist) {
		sd=vlist->var?vlist->var->data:0;
		if(sd && (sd->vtype&ST_ARRAY) && !vlist->index) i+=sd->n_elements;
		else i++;
		vlist=vlist->next;
	}
	return i;
}

static struct var_element *assign_var(struct bin_node *node,struct express *ix,struct express *express)
{
	struct var_element *element;
	struct scan_data *sd;
	
	if(!express) return 0;
	if(!(element=get_element(node,ix))) return 0;
	switch(express->type) {
	 case ST_STRING:
	     element->arg.string=express->arg.string;
		RemBlock=AddRemem(element->arg.string,RemBlock);
		break;
	 case ST_REAL:
	 case ST_INTEGER:
	     element->arg=express->arg;
		break;
	 case 0:
	     yyerror1("Undefined assignment\n");
		element->type=0;
		element->arg.string=0;
		break;
	 default:
		ABT_FUNC(IntErr);
	}
	if(!ix) {
		sd=node->data;
		sd->vtype|=ST_SCALAR;
	}
	element->type=express->type;
	return element;
}

void check_vars(struct bin_node *node,int *i,void check_func(struct bin_node *,int *))
{
	if(node->left) {
		check_vars(node->left,i,check_func);
	}
	check_func(node,i);
	if(node->right) {
		check_vars(node->right,i,check_func);
	}
}

static void check_vars_1(struct bin_node *node,void check_func(struct bin_node *))
{
	if(node->left) {
		check_vars_1(node->left,check_func);
	}
	check_func(node);
	if(node->right) {
		check_vars_1(node->right,check_func);
	}
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
	
	if(scan_warn_n<max_scan_warnings) {
		va_start(args,fmt);
		(void)vfprintf(stderr,fmt,args);
		va_end(args);
	}
	scan_warn_n++;
}

static void add_operation(void *arg,int type, int op)
{
	struct operation *o;
	
	assert(arg || !type);
	o=lk_malloc(sizeof(struct operation));
	o->next=Op_List;
	Op_List=o;
	o->type=type;
	o->op=op;
	switch(type) {
	 case VARIABLE:
		o->arg.element= (struct var_element *)arg;
		break;
	 case INTEGER:
		o->arg.value= *(int *)arg;
		break;
	 case REAL:
		o->arg.rvalue= *(double *)arg;
		break;
	 case STRING:
		o->arg.string= (char *)arg;
		break;
	}
}

static void new_command(void)
{
	shell_flag=in_loopclause=0;
	Op_List=0;
	iflag=0;
}

static struct model_list *add_to_model(struct model_list *model,struct var_list *vlist)
{
	struct var_list *vlist1;
	struct scan_data *sd;
	struct model_list *m1;
	int i;
	
	m1=lk_malloc(sizeof(struct model_list));
	if(vlist) {
		i=count_var_list(vlist);
		m1->element=lk_malloc(sizeof(void *)*i);
		i=0;
		while(vlist) {
			sd=vlist->var->data;
			if(sd->vtype&ST_ARRAY)	{
				if(vlist->index) {
					sd->element[vlist->index-1].type|=ST_MODEL;
					sd->element[vlist->index-1].index=vlist->index;
					sd->element[vlist->index-1].oindex=vlist->index;
					m1->element[i++]=sd->element+vlist->index-1;
				} else yyerror("Error - Can't use whole arrays as model parameters");
			} else {
				sd->element[0].type|=ST_MODEL;
				sd->element[0].index=0;
				m1->element[i++]=sd->element;
			}
			vlist1=vlist->next;
			free(vlist);
			vlist=vlist1;
		}
		m1->nvar=i;
	} else ABT_FUNC("Nothing to add...\n");
	m1->next=model;
	return m1;
}

static char *string_copy(char *s1,char *s2)
{
   size_t sz;	
	if(s1) { 
	   sz=strlen(s2);
		s1=lk_realloc(s1,strlen(s1)+sz+1);
		(void)strncat(s1,s2,sz);
	} else {
	  sz=strlen(s2)+1;
		s1=lk_malloc(sz);
		(void)strncpy(s1,s2,sz);
	}
	return s1;
}

static struct format *setup_format(struct format_clause *fc)
{
	int i,n=0,pp=0;
	struct format_atom **fa;
	struct format *format;
	
	fa=fc->f_atoms;
	for(i=0;i<fc->n_atoms;i++) if(!fa[i]->pos) n++;
	if(!n) {
		if(!(scan_error&FORMAT_ERR)) yyerror("Error - Empty format clause");
		free(fa);
		free(fc);
		scan_error|=FORMAT_ERR;
		scan_error_n++;
		return 0;
	}
	format=lk_malloc(sizeof(struct format));
	format->line=lineno;
	format->f_atoms=lk_malloc(sizeof(struct format_atom)*n);
	for(i=n=0;i<fc->n_atoms;i++) {
		if(!fa[i]->pos) {
			format->f_atoms[n].size=fa[i]->size;
			format->f_atoms[n++].pos=pp;
		}
		pp+=fa[i]->size;
	}
	free(fa);
	format->n_atoms=n;
	f_atom_n=0;
	free(fc);
	return format;
}

static struct format_atom *make_f_atom(int n,int flag)
{
	if(f_atom_n>=f_atom_size) {
		f_atom_size*=2;
		f_atom_list=lk_realloc(f_atom_list,sizeof(struct format_atom)*f_atom_size);
	}
	f_atom_list[f_atom_n].size=n;
	f_atom_list[f_atom_n].pos=flag;
	return &f_atom_list[f_atom_n++];
}

static struct fformat *add_fformat(struct fformat *f1,struct fformat *f2)
{
	if(f2->rs) {
		if(f1->rs) free(f1->rs);
		f1->rs=f2->rs;
	}
	if(f2->fs) {
		if(f1->fs) free(f1->fs);
		f1->fs=f2->fs;
	}
	if(f2->gs) {
		if(f1->gs) free(f1->gs);
		f1->gs=f2->gs;
	}	
	if(f2->skip) f1->skip=f2->skip;
	free(f2);
	return f1;
}

static struct fformat *create_fformat(void *p,int fg)
{
	struct fformat *ff;
	int *i;
	
	ff=lk_malloc(sizeof(struct fformat));
	ff->rs=ff->fs=ff->gs=0;
	ff->skip=0;
	switch(fg) {
	 case 1:
		ff->rs=p;
		break;
	 case 2:
		ff->fs=p;
		break;
	 case 3:
		i=p;
		ff->skip=*i;
		break;
	 case 4:
		ff->gs=p;
		break;
	 default:
		ABT_FUNC("Internal error - incorrect flag\n");
	}
	return ff;
}

static struct format_clause *add_f_atom(struct format_clause *fc,struct format_atom *fa)
{
	if(!fc) {
		fc=lk_malloc(sizeof(struct format_clause));
		fc->fc_size=16;
		fc->n_atoms=0;
		fc->f_atoms=lk_malloc(sizeof(struct format_atom *)*fc->fc_size);
	}
	if(fc->n_atoms>=fc->fc_size) {
		fc->fc_size*=2;
		fc->f_atoms=lk_realloc(fc->f_atoms,sizeof(struct format_atom *)*fc->fc_size);
	}
	fc->f_atoms[fc->n_atoms++]=fa;
	return fc;
}

static struct format_clause *add_f_list(struct format_clause *fc,struct format_clause *fc1,int n)
{
	int sz,i,j;
	
	sz=fc1->n_atoms*n;
	if(!fc) {
		fc=lk_malloc(sizeof(struct format_clause));
		fc->fc_size=16;
		if(sz>16) fc->fc_size=sz;
		fc->n_atoms=0;
		fc->f_atoms=lk_malloc(sizeof(struct format_atom *)*fc->fc_size);
	} else {
		if(sz>(fc->fc_size-fc->n_atoms)) {
			fc->fc_size=sz+fc->n_atoms;
			fc->f_atoms=lk_realloc(fc->f_atoms,sizeof(struct format_atom *)*fc->fc_size);
		}
	}
	for(i=0;i<n;i++) for(j=0;j<fc1->n_atoms;j++)
	  fc->f_atoms[fc->n_atoms++]=fc1->f_atoms[j];
	free(fc1->f_atoms);
	free(fc1);
	return fc;
}

static struct bin_node *alloc_var(char *p)
{
	struct bin_node *node;
	struct scan_data *sd;
	int i;
	
	node=lk_malloc(sizeof(struct bin_node));
	node->left=node->right=0;
	node->balance=0;
	sd=lk_malloc(sizeof(struct scan_data));
	node->data=sd;
	sd->vtype=0;
	i=(int)strlen(p);
	sd->name=lk_malloc((size_t)i+1);
	sd->name[i--]=0;
	for(;i>=0;i--) sd->name[i]=toupper((int)p[i]);
	sd->n_elements=1;
	sd->element=lk_calloc(1,sizeof(struct var_element));
	sd->element->arg.element=0;
	return node;
}

static struct bin_node *find_var(char *p,struct bin_node *node,struct bin_node **node1,int *balanced)
{
	int i;
	struct scan_data *sd;
	
	sd=node->data;
	if((i=strcasecmp(p,sd->name))) {
		if(i<0) {
			if(node->left) {
				node->left=find_var(p,node->left,node1,balanced);
			} else {
				*node1=node->left=alloc_var(p);
				*balanced=0;
			}
			if(!(*balanced)) {
				switch(node->balance) {
				 case -1:
					node=rotate_left(node);
					*balanced=1;
					break;
				 case 0:
					node->balance=-1;
					break;
				 case 1:
					node->balance=0;
					*balanced=1;
				}
			}
		} else {
			if(node->right) {
				node->right=find_var(p,node->right,node1,balanced);
			} else {
				*node1=node->right=alloc_var(p);
				*balanced=0;
			}
			if(!(*balanced)) {
				switch(node->balance) {
				 case -1:
					node->balance=0;
					*balanced=1;
					break;
				 case 0:
					node->balance=1;
					break;
				 case 1:
					node=rotate_right(node);
					*balanced=1;
				}
			}
		}
	} else {
		*node1=node;
		*balanced=1;
	}
	return node;
}

static void Check_var(struct bin_node *node)
{
	int i;
	struct var_element *element;
	struct scan_data *sd;
	char *nbuf;
	
	if(node->left) Check_var(node->left);
	sd=node->data;
	i=strlen(sd->name)+4+log((double)(sd->n_elements+1))/log(10.0);
	nbuf=lk_malloc((size_t)i);
	for(i=0;i<sd->n_elements;i++) {
		if(sd->vtype&ST_ARRAY) (void)snprintf(nbuf,i,"%s(%d)",sd->name,i+1);
		else (void)strncpy(nbuf,sd->name,(size_t)i);
		element=sd->element+i;
		if(!(element->type&(ST_DATA|ST_TRAITLOCUS|ST_LINKED|ST_LUMPED))) {
			if(element->type&(ST_MODEL|ST_SEX|ST_ID|ST_FAMILY|ST_SIRE|ST_DAM|ST_TRAIT|ST_GROUP))
			  print_scan_err("Error: No data for variable %s\n",nbuf);
		}
		if((element->type&ST_DATA) && (element->type&ST_TRAITLOCUS))
		  print_scan_err("Error: Variable %s can not have data\n",nbuf);
		else if((element->type&ST_LINKED) && !(element->type&(ST_TRAITLOCUS|ST_MARKER)))
		  print_scan_err("Error: Variable %s is not a locus and so can not be linked\n",nbuf);
		else if((element->type&ST_TRAIT) && (element->type&(ST_GROUP|ST_MARKER|ST_TRAITLOCUS|ST_ID|ST_FAMILY|ST_SIRE|ST_DAM|ST_HAPLO|ST_LUMPED|ST_LINKED|ST_STRING|ST_REAL|ST_INTEGER)))
		  print_scan_err("Error: Variable %s inappropriate type for trait\n",nbuf);
		else if((element->type&ST_TRAITLOCUS) && (element->type&(ST_SEX|ST_GROUP|ST_CENSORED|ST_RANDOM|ST_MARKER|ST_ID|ST_FAMILY|ST_SIRE|ST_DAM|ST_HAPLO|ST_STRING|ST_REAL|ST_INTEGER|ST_REALTYPE|ST_INTTYPE)))
		  print_scan_err("Error: Variable %s inappropriate type for trait locus\n",nbuf);
		else if((element->type&ST_MARKER) && (element->type&(ST_SEX|ST_GROUP|ST_CENSORED|ST_REAL|ST_INTEGER|ST_RANDOM|ST_ID|ST_FAMILY|ST_SIRE|ST_DAM|ST_HAPLO|ST_STRING|ST_INTEGER)))
		  print_scan_err("Error: Variable %s inappropriate type for marker\n",nbuf);
		else if((element->type&ST_HAPLO) && (element->type&(ST_SEX|ST_GROUP|ST_CENSORED|ST_REAL|ST_RANDOM|ST_ID|ST_FAMILY|ST_SIRE|ST_DAM|ST_STRING|ST_REAL|ST_INTEGER)))
		  print_scan_err("Error: Variable %s inappropriate type for haplotype\n",nbuf);
		else if((element->type&ST_RANDOM) && (element->type&(ST_SEX|ST_GROUP|ST_STRING|ST_REAL|ST_INTEGER|ST_REAL)))
		  print_scan_err("Error: Variable %s inappropriate type to be random\n",nbuf);
		else if((element->type&(ST_INTTYPE|ST_REALTYPE)) && (element->type&(ST_STRING|ST_REAL|ST_INTEGER)))
		  print_scan_err("Error: Type collision for variable %s\n",nbuf);
		else if((element->type&ST_INTTYPE) && (element->type&ST_REALTYPE))
		  print_scan_err("Error: Real variable %s can not also be integer type\n",nbuf);
		else if((element->type&(ST_STRING|ST_REAL)) && (element->type&(ST_SEX|ST_ID|ST_FAMILY|ST_SIRE|ST_DAM)))
		  print_scan_err("Error: Variable %s can not be a pedigree or sex variable\n",nbuf);
		else if((element->type&ST_REAL) && (element->type&(ST_ID|ST_FAMILY|ST_SIRE|ST_DAM|ST_SEX)))
		  print_scan_err("Error: Real variable %s can not be a pedigree or sex variable\n",nbuf);
		else if((element->type&ST_FACTOR) && (element->type&ST_REAL))
		  print_scan_err("Error: Real variable %s can not be a factor\n",nbuf);
		else if((element->type&ST_CONSTANT)&&(element->type&ST_MULTIPLE))
		  print_scan_err("Error: Variable %s can not be in multiple records and be constant\n",nbuf);
		else if(element->type&(ST_SEX|ST_ID|ST_FAMILY|ST_SIRE|ST_DAM|ST_TRAITLOCUS|ST_GROUP|ST_LINKED|ST_MODEL|ST_TRAIT))
		  element->type|=ST_REQUIRED;
		else if(element->type&ST_HAPLO) {
			if(element->arg.element && element->arg.element->type&ST_LINKED) {
				element->type|=ST_REQUIRED;
				if(!(element->type&ST_DATA)) print_scan_err("Error: No data for variable %s\n",nbuf);
			}
		}
		if(element->type&ST_MARKER) n_markers++;
		if(!(element->type&(ST_CONSTANT|ST_MULTIPLE))) 
		  element->type|=syst_var[MULTIPLE_RECORDS]?ST_MULTIPLE:ST_CONSTANT;
		if(element->type&(ST_MARKER|ST_REQUIRED|ST_RESTRICT))	{
			if(!(element->type&(ST_HAPLO|ST_LUMPED))) element->arg.var=node;
		} else element->type=0;
	}
	free(nbuf);
	if(node->right) Check_var(node->right);
}

static struct bin_node *create_var(char *p)
{
	int k;
	struct bin_node *node;
	
	if(!root_var) node=root_var=alloc_var(p);
	else {
		root_var=find_var(p,root_var,&node,&k);
	}
	return node;
}

int symbol_lookup(char *p,int fg)
{
	static char *Coms[] = {"FILE","LOCUS","LOCI","MARKER","DISCRETE","MODEL","PEDIGREE","LOG",
		"FILTER","MISSING","MODEL","LINK","RANDOM","TRAIT","WHERE","USE",
		"REAL","INTEGER","SHELL","ARRAY","PRINT","DO","WHILE","CONSTANT",
		"MULTIPLE","CENSORED","GROUP","SET","SEX","AFFECTED","UNAFFECTED","PROBAND","OUTPUT","INCLUDE","ERRORDIR",
		"LAUROUTPUT","RAWOUTPUT","POSITION","FREQUENCY","SUPER",(char *)0};
	static int Com_token[] = {FILEC,LOCUS,LOCUS,MARKER,FACTOR,MODEL,PEDIGREE,LOG,
		FILTER,MISSING,MODEL,LINK,RANDOM,TRAIT,WHERE,USE,
		REAL,INTEGER,SHELL,ARRAY,PRINTEXP,DOLOOP,WHILE,CONSTANT,
		MULTIPLE,CENSORED,GROUP,SET,GENDER,AFFECTED,UNAFFECTED,PROBAND,OUTPUT,INCLUDE,ERRORDIR,
		LAUROUTPUT,RAWOUTPUT,POSITION,FREQUENCY,SUPER,SYSTEM_VAR,VARIABLE,ARRAY_VAR};
	static char *Syst[] = {"PRUNE_OPTION","RECODE_OPTION","NO_EXTRA_ALLELE",
		"PEEL_OPTION","TRACE_RESTRICT","TRACE_CENSORED","TRACE_AFFECTED",
		"CORRECT_ERRORS","TRACE_PEEL","MULTIPLE_RECORDS","MULTIVARIATE_TEST",
		"ERROR_CHECK","NO_DEFAULT_MISSING","SKIP_BAD_REALS","SKIP_BAD_INTS","IGNORE_CASE",(char *)0};
	int i=0,j=0;
	static struct scan_data *sd;
	
	while(Coms[i])	{
		if(!strcasecmp(Coms[i],p)) break;
		i++;
	}
	at_file=0;
	if(Com_token[i]==FILEC || Com_token[i]==LINK) at_file=1;
	if(Com_token[i]==SYSTEM_VAR) {
		i++;
		while(Syst[j])	{
			if(!strcasecmp(Syst[j],p))	{
				yylval.value=j;
				i--;
				break;
			}
			j++;
		}
	}
	if(Com_token[i]==VARIABLE) {
		if(fg==1 && begin_comm) {
			begin_comm=0;
			return BREAK;
		}
		yylval.var=create_var(p);
		sd=yylval.var->data;
		if(sd->vtype&ST_ARRAY) i++;
		if(fg==1) {
			begin_comm=1;
			(void)strncpy(linebuf1,linebuf,LINEBUFSIZE);
			lineno1=lineno;
		}
	} else if(begin_comm && Com_token[i]!=SYSTEM_VAR && Com_token[i]!=LOCUS && Com_token[i]!=SHELL
				 && !(at_use==1 && Com_token[i]==WHERE) && !(at_use==2 && Com_token[i]==USE))	{
		begin_comm=0;
		at_use=0;
		return BREAK;
	} else {
		begin_comm=1;
		(void)strncpy(linebuf1,linebuf,LINEBUFSIZE);
		lineno1=lineno;
		if(Com_token[i]==MODEL) at_model=1;
		else at_model=0;
		if(Com_token[i]==USE || Com_token[i]==CENSORED || Com_token[i]==AFFECTED || Com_token[i]==UNAFFECTED || Com_token[i]==PROBAND) at_use|=1;
		else if(Com_token[i]==WHERE) at_use|=2;
		else at_use=0;
	}
	return Com_token[i];
}

static struct var_list *add_to_var_list(struct var_list *vlist,struct bin_node *node,struct express *express)
{
	struct var_list *vlist1,*vlist2;
	struct scan_data *sd=0;
	int i;

	if(node)	sd=node->data;
	if(express) i=check_index(sd,express);
	else {
		i=0;
		if(sd && !(sd->vtype&ST_ARRAY)) sd->vtype|=ST_SCALAR;
	}
	vlist1=lk_malloc(sizeof(struct var_list));
	vlist1->next=0;
	vlist1->var=node;
	vlist1->index=i;
	vlist2=vlist;
	if(vlist2) {
		while(vlist2->next) vlist2=vlist2->next;
		vlist2->next=vlist1;
	} else vlist=vlist1;
	return vlist;
}

struct var_list *add_var_lists(struct var_list *vlist,struct var_list *vlist1)
{
	struct var_list *vlist2;
	
	vlist2=vlist;
	if(vlist2) {
		while(vlist2->next) vlist2=vlist2->next;
		vlist2->next=vlist1;
	} else vlist=vlist1;
	return vlist;
}

static void set_locus_array(struct bin_node *node)
{
	struct scan_data *sd;
	int i;
	
	sd=node->data;
	if(sd->vtype&ST_ARRAY) {
		for(i=0;i<sd->n_elements;i++) {
			set_locus_element(sd->element+i);
		}
	} else yyerror("Not an array");
}

static void set_locus_element(struct var_element *element)
{
	element->type|=(ST_MARKER|ST_FACTOR|ST_CONSTANT);
	if(hap_list[0]) {
		if(hap_list[0]->arg.element && hap_list[0]->arg.element!=element)	{
			yyerror1("Haplotype vector (left) used twice");
			hap_list[0]->arg.element=0;
		} else hap_list[0]->arg.element=element;
	}
	if(hap_list[1]) {
		if(hap_list[1]->arg.element && hap_list[1]->arg.element!=element)	{
			yyerror1("Haplotype vector (right) used twice");
			hap_list[1]->arg.element=0;
		} else hap_list[1]->arg.element=element;
	}
	hap_list[0]=hap_list[1]=0;
}

static void set_slocus_element(struct var_element *element,struct var_list *vlist)
{
	int j;
	struct scan_data *sd;
	
	element->type|=(ST_MARKER|ST_FACTOR|ST_CONSTANT);
	while(vlist) {
		sd=vlist->var->data;
		assert(sd);
		if(sd->vtype&ST_ARRAY) {
			if(vlist->index) {
				sd->element[vlist->index-1].type|=ST_LUMPED;
				sd->element[vlist->index-1].arg.element=element;
			} else for(j=0;j<sd->n_elements;j++) {
				sd->element[j].type|=ST_LUMPED;
				sd->element[j].arg.element=element;
			}
		} else {
			sd->element[0].type|=ST_LUMPED;
			sd->element[0].arg.element=element;
		}
		vlist=vlist->next;
	}
}

static void set_haplo_element(struct var_element *element,struct var_element *element1)
{
	element->type|=(ST_HAPLO|ST_FACTOR|ST_CONSTANT);
	if(element1) element1->type|=(ST_HAPLO|ST_FACTOR|ST_CONSTANT);
	hap_list[0]=element;
	hap_list[1]=element1;
}

static void do_file_com(char *fname,struct format_clause *fc,struct fformat *ff,struct var_list *vlist)
{
	int i,j;
	struct InFile *file;
	struct format *format;
	struct var_list *vlist1;
	struct var_element *element;
	struct scan_data *sd;
	
	if(!vlist) {
		yyerror1("No variables listed for FILE command\n");
		return;
	} else if(!fname) {
		free_vlist(vlist);
		return;
	} else if(!fname[0]) {
		yyerror1("Zero length filename for FILE command\n");
		free_vlist(vlist);
		return;
	}
	file=Infiles;
	Infiles=lk_calloc(1,sizeof(struct InFile));
	Infiles->next=file;
	Infiles->nvar=count_var_list(vlist);
	Infiles->element=lk_malloc(sizeof(void *)*Infiles->nvar);
	i=0;
	while(vlist) {
		if(vlist->var) {
			sd=vlist->var->data;
			if(sd->vtype&ST_ARRAY) {
				if(vlist->index) {
					element=sd->element+vlist->index-1;
					element->type|=ST_DATA;
					Infiles->element[i++]=element;
				} else {
					for(j=0;j<sd->n_elements;j++) {
						element=sd->element+j;
						element->type|=ST_DATA;
						Infiles->element[i++]=element;
					}
				}
			} else {
				element=sd->element;
				element->type|=ST_DATA;
				Infiles->element[i++]=element;
			}
		} else Infiles->element[i++]=0;
		vlist1=vlist->next;
		free(vlist);
		vlist=vlist1;
	}
	if(fc) {
		format=setup_format(fc);
		Infiles->format=format;
		if(!(scan_error&FORMAT_ERR)) {
			if(format->n_atoms<i) {
				(void)printf("format->n_atoms = %d\n",format->n_atoms);
				(void)printf("i = %d\n",i);
				print_scan_err("Line %d: Error - Too many variables for format clause\n",format->line);
				scan_error|=FORMAT_ERR;
			} else if(format->n_atoms>i)
				  print_scan_warn("Line %d: Warning - Too few variables for format clause\n",format->line);
		}
	} else if(ff) Infiles->fformat=ff;
	Infiles->name=fname;
	Infiles->shell_flag=shell_flag;
}

static void change_type(int type,struct var_list *vlist)
{
	int j;
	struct scan_data *sd;
	
	
	while(vlist) {
		sd=vlist->var->data;
		if(sd->vtype&ST_ARRAY) {
			if(vlist->index) sd->element[vlist->index-1].type|=type;
			else for(j=0;j<sd->n_elements;j++)
			  sd->element[j].type|=type;
		} else sd->element[0].type|=type;
		vlist=vlist->next;
	}
}

static void free_vlist(struct var_list *vlist)
{
	struct var_list *vlist1;
	
	while(vlist) {
		vlist1=vlist->next;
		free(vlist);
		vlist=vlist1;
	}
}

static void do_link_com(char *s,int type,struct var_list *vlist)
{
	struct Link *l,*l1,**ll;
	struct var_list *vlist1;
	struct var_element *element;
	struct scan_data *sd=0;
	int i,j,k;
	
	if(vlist) sd=vlist->var->data;
	if(!s && sd) {
		if(sd->vtype&ST_ARRAY && vlist->index) {
			element=sd->element+vlist->index-1;
			if(element->type&ST_STRING) {
				s=element->arg.string;
				vlist1=vlist->next;
				free(vlist);
				vlist=vlist1;
			}
		} else {
			element=sd->element;
			if(element->type&ST_STRING) {
				s=element->arg.string;
				vlist1=vlist->next;
				free(vlist);
				vlist=vlist1;
			}
		}
	}
	ll=&links;
	while(*ll) {
		l=*ll;
		if(s) {
			if(l->name) {
				if(!strcasecmp(s,l->name)) break;
			}
		} else if(!l->name) break;
		ll=&l->next;
	}
	if(*ll) l1=*ll;
	else {
		l1=lk_malloc(sizeof(struct Link));
		l1->next=0;
		l1->name=s;
		l1->n_loci=0;
		l1->element=0;
		l1->type=-1;
		*ll=l1;
	}
	i=count_var_list(vlist);
	if(l1->type>=0 && l1->type!=type) print_scan_err("Error: Linkage group has inconsistent linkage type\n");
	l1->type=type;
	if(i) {
		k=i+l1->n_loci;
		if(l1->element) {
			l1->element=lk_realloc(l1->element,sizeof(void *)*k);
		} else l1->element=lk_malloc(sizeof(void *)*k);
		while(vlist) {
			sd=vlist->var->data;
			if(sd->vtype&ST_ARRAY) {
				if(vlist->index) {
					element=sd->element+vlist->index-1;
					for(i=0;i<l1->n_loci;i++) {
						if(l1->element[i]==element) break;
					}
					if(i==l1->n_loci) {
						if(element->type&ST_LINKED) {
							print_scan_err("Error: %s(%d) appears in multiple linkage groups\n",sd->name,vlist->index);
							scan_error|=LINK_ERR;
						} else {
							element->type|=ST_LINKED;
							l1->element[l1->n_loci++]=element;
						}
					}
				} else {
					for(j=0;j<sd->n_elements;j++) {
						element=sd->element+j;
						for(i=0;i<l1->n_loci;i++) {
							if(l1->element[i]==element) break;
						}
						if(i==l1->n_loci) {
							if(element->type&ST_LINKED) {
								print_scan_err("Error: %s(%d) appears in multiple linkage groups\n",sd->name,vlist->index);
								scan_error|=LINK_ERR;
							} else {
								element->type|=ST_LINKED;
								l1->element[l1->n_loci++]=element;
							}
						}
						sd->vtype|=ST_LINKED;
					}
				}
			} else {
				element=sd->element;
				for(i=0;i<l1->n_loci;i++) {
					if(l1->element[i]==element) break;
				}
				if(i==l1->n_loci) {
					if(element->type&ST_LINKED) {
						print_scan_err("Error: %s appears in multiple linkage groups\n",sd->name);
						scan_error|=LINK_ERR;
					} else {
						element->type|=ST_LINKED;
						l1->element[l1->n_loci++]=element;
					}
				}
			}
			vlist1=vlist->next;
			free(vlist);
			vlist=vlist1;
		}
		if(l1->n_loci<k) l1->element=lk_realloc(l1->element,sizeof(void *)*l1->n_loci);
	}
}

static void do_missing_com(struct express *expr,struct var_list *vlist,char *s1)
{
	struct var_list *vlist1;
	struct scan_data *sd;
	struct var_element **elem;
	struct Miss *m;
	int i,j;
	char *p;
	
	if(s1) {
		assert(!vlist);
		if(s1[0]==0) {
			print_scan_err("Empty scope - MISSING directive ignored\n");
			if(expr->type==ST_STRING) free(expr->arg.string);
			free(s1);
			return;
		}
		qstrip(s1);
		p=s1;
		i=j=0;
		while(*p) {
			switch(toupper((int)*p)) {
			 case '!':
			 case 'F':
			 case 'G':
			 case 'P':
			 case 'C':
			 case 'R':
			 case 'I':
				break;
			 default: i=1;
			}
			if(i) break;
			p++;
		}
		if(*p) {
			j=1;
			print_scan_err("Illegal character '%c' in MISSING scope\n",*p);
		} else if(*(--p)=='!') {
			j=1;
			print_scan_err("MISSING scope can not end with a '!'\n",*p);
		}
		if(j) {
			free(s1);
			if(expr->type==ST_STRING) free(expr->arg.string);
			return;
		}
	}
	m=Miss;
	Miss=lk_malloc(sizeof(struct Miss));
	Miss->Missing.arg=expr->arg;
	Miss->Missing.type=expr->type;
	Miss->next=m;
	Miss->element=0;
	Miss->scope=0;
	if((i=count_var_list(vlist))) {
		elem=lk_malloc(sizeof(void *)*i);
		i=0;
		while(vlist) {
			sd=vlist->var->data;
			if(sd->vtype&ST_ARRAY) {
				if(vlist->index) elem[i++]=sd->element+vlist->index-1;
				else for(j=0;j<sd->n_elements;j++) elem[i++]=sd->element+j;
			} else elem[i++]=sd->element;
			Miss->element=elem;
			vlist1=vlist->next;
			free(vlist);
			vlist=vlist1;
		}
	} else if(s1) Miss->scope=s1;
	Miss->nvar=i;
}

static void do_model_com(struct model_list *mlist,struct bin_node *node,struct express *express)
{
	struct model *model,*model1;
	struct var_element *element;
	struct scan_data *sd;
	
	sd=node->data;
	model=lk_malloc(sizeof(struct model));
	model->next=0;
   if(Models) {
		model1=Models;
		while(model1->next) model1=model1->next;
		model1->next=model;
	} else Models=model;
	model->trait=sd;
	if(!express) {
		model->index=0;
		sd->element[0].type|=ST_TRAIT;
	} else {
		element=get_element(node,express);
		if(element) {
			model->index=(int)express->arg.value;
			element->index=element->oindex=model->index;
			element->type|=ST_TRAIT;
		} else model->trait=0;
	}
	model->model_list=mlist;
}

static void add_censored(struct var_element *element,const int fg)
{
	struct operation *ops;
	struct Censor *cen;
	
	if(fg==1 && !element) {
		print_scan_err("Error: Nothing to censor!\n");
		return;
	}
	/* Reverse list order (really return list to original order! */
	ops=Op_List=reverse_list(Op_List);
	switch(fg) {
	 case 1:
		cen=lk_malloc(sizeof(struct Censor));
		cen->next=Censored;
		Censored=cen;
		cen->Op_List=ops;
		cen->element=element;
		element->type|=ST_CENSORED;
		break;
	 case 0:
		if(Affected) {
			print_scan_warn("Warning - new affected statement overrules previous statement\n");
			free_op(Affected);
		}
		Affected=ops;
		break;
	 case 2:
		if(Unaffected) {
			print_scan_warn("Warning - new unaffected statement overrules previous statement\n");
			free_op(Unaffected);
		}
		Unaffected=ops;
		break;
	 case 3:
		if(Proband) {
			print_scan_warn("Warning - new proband statement overrules previous statement\n");
			free_op(Proband);
		}
		Proband=ops;
		break;
	}
	ops=Op_List;
	while(ops) {
		if(ops->type==VARIABLE) ops->arg.element->type|=ST_RESTRICT;
		ops=ops->next;
	}
}

static void add_restriction(struct var_list *vlist)
{
	struct operation *ops;
	struct Restrict *res;
	struct var_list *vlist1;
	struct scan_data *sd;
	int i,j;
	
	/* Reverse list order (really return list to original order! */
	Op_List=ops=reverse_list(Op_List);
	res=lk_malloc(sizeof(struct Restrict));
	res->next=Restrictions;
	Restrictions=res;
	res->Op_List=ops;
	if((res->nvar=count_var_list(vlist))) {
		res->element=lk_malloc(sizeof(void *)*res->nvar);
		i=0;
		while(vlist) {
			sd=vlist->var->data;
			if(sd->vtype&ST_ARRAY) {
				if(vlist->index) res->element[i++]=sd->element+vlist->index-1;
				else for(j=0;j<sd->n_elements;j++) res->element[i++]=sd->element+j;
			} else res->element[i++]=sd->element;
			vlist1=vlist->next;
			free(vlist);
			vlist=vlist1;
		}	
	} else res->element=0;
	while(ops) {
		if(ops->type==VARIABLE) ops->arg.element->type|=ST_RESTRICT;
		ops=ops->next;
	}
}

static void find_markers(struct bin_node *node,int *i)
{
	int j;
	struct scan_data *sd;
	
	sd=node->data;
	for(j=0;j<sd->n_elements;j++) if(sd->element[j].type&ST_MARKER) {
		if(sd->element[j].type&ST_REQUIRED) {
			markers[*i].element=sd->element+j;
			markers[*i].var=sd;
			markers[(*i)++].index=sd->n_elements>1?j+1:0;
		}
	}
}

static void find_trait_loci(struct bin_node *node,int *i)
{
	int j;
	struct scan_data *sd;
	
	sd=node->data;
	for(j=0;j<sd->n_elements;j++) if(sd->element[j].type&ST_TRAITLOCUS) {
		if(traitlocus) {
			traitlocus[*i].element=sd->element+j;
			traitlocus[*i].var=sd;
			traitlocus[(*i)].index=j+1;
		}
		(*i)++;
	}
}

static void handle_super_loci(struct bin_node *node)
{
	int k;
	struct scan_data *sd;
	
	sd=node->data;
	for(k=0;k<sd->n_elements;k++) if(sd->element[k].type&ST_LUMPED) {
		if(sd->element[k].arg.element->type&ST_LINKED) {
			if(!(sd->element[k].type&ST_MARKER)) {
				if(sd->n_elements>1) print_scan_err("Error: Variable '%s(%d)' is not a marker and so can not be part of a super locus\n",sd->name,k+1);
				else print_scan_err("Error: Variable '%s' is not a marker and so can not be part of a super locus\n",sd->name);
			} else sd->element[k].type|=ST_LINKED;
		}
	}
}

static void find_haplo(struct bin_node *node)
{
	int j,k;
	struct scan_data *sd;
	
	sd=node->data;
	for(k=0;k<sd->n_elements;k++) if(sd->element[k].type&ST_HAPLO) {
		for(j=0;j<n_markers;j++) if(sd->element[k].arg.element==markers[j].element) {
			if(!markers[j].hap_element[0]) markers[j].hap_element[0]=sd->element+k;
			else if(!markers[j].hap_element[1])	markers[j].hap_element[1]=sd->element+k;
			else {
				if(markers[j].index) print_scan_err("Error: marker %s(%d) has >2 haplotype vectors associated with it\n",markers[j].var->name,markers[j].index);
				else print_scan_err("Error: marker %s has >2 haplotype vectors associated with it\n",markers[j].var->name);
			}
			break;
		}
		assert(j<n_markers);
	}
}

static void find_lumped1(struct bin_node *node)
{
	int j,k;
	struct scan_data *sd;
	
	sd=node->data;
	for(k=0;k<sd->n_elements;k++) if(sd->element[k].type&ST_LUMPED) {
		if(sd->element[k].type&ST_LINKED) {
			for(j=0;j<n_markers;j++) if(sd->element[k].arg.element==markers[j].element) {
				if(markers[j].hap_element[0]) {
					if(markers[j].index) print_scan_err("Error: super locus %s(%d) has haplotype vectors associated with it\n",markers[j].var->name,markers[j].index);
					else print_scan_err("Error: super locus %s has haplotype vectors associated with it\n",markers[j].var->name);
				}
				markers[j].n_sub_elements++;
				break;
			}
			assert(j<n_markers);
		}
	}
}

static void find_lumped2(struct bin_node *node)
{
	int j,k;
	struct scan_data *sd;
	
	sd=node->data;
	for(k=0;k<sd->n_elements;k++) if(sd->element[k].type&ST_LUMPED) {
		if(sd->element[k].type&ST_LINKED) {
			for(j=0;j<n_markers;j++) if(sd->element[k].arg.element==markers[j].element) {
				markers[j].sub_element[markers[j].n_sub_elements++]=sd->element+k;
				break;
			}
		}
	}
}

static void correct_linkage(struct bin_node *node)
{
	int j,k,k1;
	struct scan_data *sd;
	char *p1,*p2;
	
	sd=node->data;
	for(k=0;k<sd->n_elements;k++) if(sd->element[k].type&ST_LUMPED) {
		if(sd->element[k].type&ST_LINKED) {
			for(j=0;j<n_markers;j++) if(sd->element+k==markers[j].element) break;
			assert(j<n_markers);
			for(k1=0;k1<n_markers;k1++) if(sd->element[k].arg.element==markers[k1].element) {
				break;
			}
			assert(k1<n_markers);
			if(markers[j].link) {
				p1=get_marker_name(j);
				p2=get_marker_name(k1);
				print_scan_err("Error: Marker '%s' is a part of superlocus '%s' so can not appear directly in a LINK statement\n",p1,p2);
				free(p1);
				free(p2);
			}
			markers[j].link=markers[k1].link;
		}
	}
}

static void strip_names(struct bin_node *node)
{
	char *p;
	int i;
	struct scan_data *sd;
	
	sd=node->data;
	if((p=sd->name)) {
		i=strlen(p);
		if(i>2) {
			if(p[i-1]=='_' && p[0]=='_') {
				p[i-1]=0;
				memmove(p,p+1,i-1);
			}
		}
	}
}

int ReadControl(FILE *fptr,char *cname,char **lfile)
{
	int i,j,k;
	void yy_cleanup(void);
	struct InFile *infile,**infile_p;
	struct Restrict *res,*res1,**res_p;
	struct Censor *cen,**cen_p;
	struct var_element *elem;
	struct Link *linkp;
	struct operation *ops;
	struct express tmp_expr;
	struct marker_info *mi;
	
	yyin=fptr;
	fname_list[0]=cname;
	list_ptr=0;
	f_atom_list=lk_malloc(sizeof(struct format_atom)*f_atom_size);
	for(i=0;i<NUM_SYSTEM_VAR;i++) syst_var[i]=0;
	syst_var[PRUNE_OPTION]=syst_var[RECODE_OPTION]=2;
	syst_var[ERROR_CHECK]=1;
	if((i=yyparse())) print_scan_err("Error: yyparse returned error %d\n",i);
	yy_cleanup();
	if(strip_vars) check_vars_1(root_var,strip_names);
	/* Sanity check! */
	if(!scan_error_n)	{
		if(!Infiles) print_scan_err("Error: No input files specified\n");
		if(!pedflag) print_scan_err("Error: No pedigree variables specified\n");
		else {
			for(i=0;i<3;i++) if(pedlist[i+family_id]->type&ST_INTTYPE) break;
			if(i<3) for(i=0;i<3;i++) pedlist[i+family_id]->type|=ST_INTTYPE;
		}
		if(root_var) {
			check_vars_1(root_var,handle_super_loci);
			Check_var(root_var);
		}
		/* Flag variables used as the operands to a restriction statement *whose result is used* as ST_REQUIRED */
		res=0;
		while(res!=Restrictions) {	
			res1=Restrictions;
			while(res1->next!=res) res1=res1->next;
			for(i=j=0;i<res1->nvar;i++) if(res1->element[i]->type&ST_REQUIRED) {
				j=1;
				break;
			}
			if(!res1->nvar || j) {
				ops=res1->Op_List;
				while(ops) {
					if(ops->type==VARIABLE) ops->arg.element->type|=ST_REQUIRED;
					ops=ops->next;
				}
			}
			res=res1;
		}
		/* Delete restrict structures that are not used */
		res=Restrictions;
		res_p= &Restrictions;
		while(res) {
			for(i=j=0;i<res->nvar;i++) if(res->element[i]->type&ST_REQUIRED) {
				j=1;
				break;
			}
			if(res->nvar && !j) {
				*res_p=res->next;
				free_restrict(res);
				res= *res_p;
			} else {
				res_p= &res->next;
				res=res->next;
			}
		}
		if(Unaffected && !Affected) print_scan_err("Error: Unaffected definition without affected definition\n");
		/* Flag variables used in censored statements as required.  Delete unused censored statements */
		cen=Censored;
		cen_p= &Censored;
		while(cen) {
			if(cen->element->type&ST_TRAIT) {
				ops=cen->Op_List;
				while(ops) {
					if(ops->type==VARIABLE) ops->arg.element->type|=ST_REQUIRED;
					ops=ops->next;
				}
				cen_p= &cen->next;
				cen=cen->next;
			} else {
				*cen_p=cen->next;
				free_op(cen->Op_List);
				free(cen);
				cen= *cen_p;
			}
		}
		if((ops=Affected)) {
			while(ops) {
				if(ops->type==VARIABLE) ops->arg.element->type|=ST_REQUIRED;
				ops=ops->next;
			}
		}
		if((ops=Unaffected)) {
			while(ops) {
				if(ops->type==VARIABLE) ops->arg.element->type|=ST_REQUIRED;
				ops=ops->next;
			}
		}
		if((ops=Proband)) {
			while(ops) {
				if(ops->type==VARIABLE) ops->arg.element->type|=ST_REQUIRED;
				ops=ops->next;
			}
		}
		/* Check file structures - remove ones that aren't needed */
		infile=Infiles;
		infile_p= &Infiles;
		while(infile) {
			infile->ncol=0;
			for(i=0;i<infile->nvar;i++) {
				elem=infile->element[i];
			   if(elem) elem->type&=~ST_FLAG;
			}
			for(k=infile->nvar-1;k>=0;k--) {
				elem=infile->element[k];
				if(elem && (elem->type&ST_REQUIRED)) break;
			}
			for(i=j=0;i<infile->nvar;i++)	{
				elem=infile->element[i];
				if(elem) {
					if((elem->type&ST_MARKER) && !(elem->type&ST_REQUIRED)) {
						if(j<k) elem->type|=(ST_REQUIRED|ST_NOT_REALLY_REQUIRED);
						else elem->type=0;
					}
					if(elem->type&ST_REQUIRED) {
						if(elem->type&ST_FLAG) {
							print_scan_err("Error: Duplicate variables for file %s\n",infile->name);
							break;
						}
						elem->type|=ST_FLAG;
						if(elem->type&ST_ID)	{
							j|=1;
							infile->id_col=infile->ncol;
						} else if(elem->type&ST_FAMILY) {
							j|=2;
							infile->family_col=infile->ncol;
						}
						infile->ncol++;
					} else infile->element[i]=0;
				}
			}
			for(i=0;i<infile->nvar;i++) if(infile->element[i]) infile->element[i]->type&=~ST_FLAG;
			if(!(j&1)) print_scan_err("Error: No id column for file %s\n",infile->name);
			else if(family_id && j!=3) print_scan_err("Error: No family column for file %s\n",infile->name);
			if(infile->ncol==1) {
				*infile_p=infile->next;
				free_infile(infile);
				infile= *infile_p;
			} else {
				infile_p= &infile->next;
				infile=infile->next;
			}
		}
		if(!Infiles) print_scan_err("Error: No input files with data\n");
		free(f_atom_list);
		/* Count markers and link up with haplotype vectors */
		if(n_markers) {
			markers=lk_calloc((size_t)n_markers,sizeof(struct Marker));
			for(i=0;i<n_markers;i++) {
				markers[i].allele_trans=0;
			   markers[i].order=0;
				markers[i].o_size=0;
				markers[i].pos_set[0]=markers[i].pos_set[1]=0;
				markers[i].hap_element[0]=markers[i].hap_element[1];
				markers[i].sub_element=0;
				markers[i].n_sub_elements=0;
			}
			i=0;
			if(root_var) {
				check_vars(root_var,&i,find_markers);
				check_vars_1(root_var,find_haplo);
				check_vars_1(root_var,find_lumped1);
				for(i=0;i<n_markers;i++) if(markers[i].n_sub_elements) {
					markers[i].sub_element=lk_malloc(sizeof(void *)*markers[i].n_sub_elements);
					markers[i].n_sub_elements=0;
				}
				check_vars_1(root_var,find_lumped2);
			}
			n_markers=i;
			while(m_info) {
				mi=m_info->next;
				for(i=0;i<n_markers;i++) if(markers[i].element==m_info->element) {
					for(k=0;k<3;k++) {
						markers[i].pos_set[k]=m_info->pos_set[k];
						markers[i].pos[k]=m_info->pos[k];
					}
					break;
				}
				free(m_info);
				m_info=mi;
			}
			for(i=0;i<n_markers;i++) {
				if(!markers[i].element || markers[i].element->type&ST_NOT_REALLY_REQUIRED) continue;
				linkp=0;
				j=0;
				linkp=links;
				while(linkp) {
					j++;
					for(k=0;k<linkp->n_loci;k++) {
						if(linkp->element[k]==markers[i].element) {
							markers[i].link=j;
							break;
						}
					}
					if(k<linkp->n_loci) break;
					linkp=linkp->next;
				}
				if(!linkp && !(markers[i].element->type&ST_LUMPED)) {
					if(markers[i].var->vtype&ST_ARRAY)
					  abt(__FILE__,__LINE__,"%s(): No linkage group specified for candidate gene %s(%d)\n",__func__,markers[i].var->name,markers[i].index);
					else abt(__FILE__,__LINE__,"%s(): No linkage group specified for candidate gene %s\n",__func__,markers[i].var->name);
				}
				if(markers[i].n_sub_elements) {
					assert(!markers[i].hap_element[0]);
					if(markers[i].element->type&ST_DATA) {
						if(markers[i].var->vtype&ST_ARRAY)
						  print_scan_err("Error: super locus %s(%d) can not have its own data\n",markers[i].var->name,markers[i].index);
						else
						  print_scan_err("Error: super locus %s can not have its own data\n",markers[i].var->name);
					}
				} else if(markers[i].hap_element[0]) {
					if(markers[i].element->type&ST_DATA) {
						if(markers[i].var->vtype&ST_ARRAY)
						  print_scan_err("Error: marker variable %s(%d) can not have both genotype and haplotype data\n",markers[i].var->name,markers[i].index);
						else
						  print_scan_err("Error: marker variable %s can not have both genotype and haplotype data\n",markers[i].var->name);
					}
					if(markers[i].hap_element[0]->type&ST_INTTYPE) markers[i].hap_element[1]->type|=ST_INTTYPE;
					if(markers[i].hap_element[1] && markers[i].hap_element[1]->type&ST_INTTYPE) markers[i].hap_element[0]->type|=ST_INTTYPE;
				} else {
					if(!(markers[i].element->type&ST_DATA)) {
						if(markers[i].var->vtype&ST_ARRAY)
						  print_scan_err("Error: marker variable %s(%d) has no data\n",markers[i].var->name,markers[i].index);
						else
						  print_scan_err("Error: marker variable %s has no data\n",markers[i].var->name);
					}
				}
			}
		}
		i=0;
		if(root_var) {
			check_vars_1(root_var,correct_linkage);
			check_vars(root_var,&i,find_trait_loci);
		}
		if(i) {
			if(i>1) print_scan_err("Error: multiple trait loci indicated\n");
			else {
				traitlocus=lk_calloc(1,sizeof(struct Marker));
				traitlocus->sub_element=0;
				traitlocus->n_sub_elements=0;
				traitlocus->order=0;
				traitlocus->o_size=0;
				i=0;
				check_vars(root_var,&i,find_trait_loci);
			}
		}
		if(Models && Models->next && !syst_var[MULTIVARIATE_TEST]) {
			print_scan_err("Error: Multiple models not currently supported\n");
		}
	}
	*lfile=LogFile;
	if(!scan_error_n && !Miss && !syst_var[NO_DEFAULT_MISSING]) {
		tmp_expr.arg.string=strdup("0");
		tmp_expr.type=ST_STRING;
		do_missing_com(&tmp_expr,0,strdup("PF"));
	}
	return scan_error_n;
}
