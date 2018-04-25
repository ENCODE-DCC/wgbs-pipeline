%pure_parser
%{
/****************************************************************************
 *                                                                          *
 *     Loki - Programs for genetic analysis of complex traits using MCMC    *
 *                                                                          *
 *                     Simon Heath - CNG, Evry                              *
 *                                                                          *
 *                           June 2003                                      *
 *                                                                          *
 * parse.y:                                                                 *
 *                                                                          *
 * Basic control file parser                                                *
 *                                                                          *
 * Copyright (C) Simon C. Heath 2003                                        *
 * This is free software.  You can distribute it and/or modify it           *
 * under the terms of the Modified BSD license, see the file COPYING        *
 *                                                                          *
 ****************************************************************************/
#include <stdlib.h>
#include <string.h>
#ifdef USE_DMALLOC
#include <dmalloc.h>
#endif
#include <stdio.h>
#include <stdarg.h>
#include <ctype.h>
#include <math.h>
#include <sys/stat.h>
	
#include "string_utils.h"
#include "bin_tree.h"
#include "parser.h"
#include "parse.tab.h"

#undef YYLSP_NEEDED
#undef YYLEX_PARAM

#define BUFSIZE 640
#define LOOKAHEAD 128

static int parerror(char *);
static int parlex(yystype *);
static char *key_err="Can not use keyword '%s' as a variable\n";
static parse_handler *handle;
	
#define ALPHA_BIT 1
#define DIGIT_BIT 2
#define E_BIT 4
#define USCORE_BIT 8
#define PERIOD_BIT 16
#define PM_BIT 32
#define MATH_BIT 64
#define REL_BIT 128
#define LBRACK_BIT 256
#define RBRACK_BIT 512
#define COMMA_BIT 1024
#define X_BIT 2048
#define TOKEN_BITS (PM_BIT-1)
#define FORMAT_BITS (DIGIT_BIT|X_BIT|LBRACK_BIT|RBRACK_BIT|COMMA_BIT)

struct lex_info {
	struct lex_info *next,*last;
	int line_no;
	int col_no,col_no1;
	int eof;
	char *name;
	int ptr,ptr1;
	int last_token;
	int flag;
	yystype last_lval;
	FILE *fptr;
	char buf[BUFSIZE];
};

static void include_file(string *);
static void copy_filepos(filepos *,filepos *);

static struct lex_info *lex_info;
static struct bin_node *root_var;
static int lex_tab[256];
static char **key_tab;
static func_def *func_list;
static struct parse_term *vtmp;
	
%}

%union {
	int i;
	double x;
	string *str;
	struct parse_term *term;
	struct parse_var *vv;
	struct parse_var_list *vl;
	struct format_clause *fc;
	struct parse_clause *pc;
	func_def *f;
}

%token KEYWORD KEYWORD_B ALIAS FOR IF ELSE DO WHILE INCLUDE NE_TOKEN LEQ_TOKEN GEQ_TOKEN CMP_TOKEN
%token VARIABLE EXPR_OP INF XNAN DUMMY DOUBLE_DOT EMPTY_VAR EXPO CLAUSE_X EQ_TOKEN SUB_EXPR
  
%token <str> STRING
%token <i> INTEGER CLAUSE_INT PREFIX POSTFIX MATH_SHORT KEYWORD KEYWORD_B
%token <x> REAL
%token <vv> WORD WORD_B
%token <f> FUNCTION
  
%type <term> gen_expr gen_expr_list gen_expr_list1 gen_element var1 array_index array_index_elem 
%type <term> imm_expr imm_element imm_element1
%type <vl> var var2 var_list array_var_list loop_clause var_list1 poly_commas tok tok_list
%type <str> str_elem str_elem1 str_expr 
%type <fc> format_atom format_clause
%type <pc> clause_elem clause_content clause
  
%nonassoc LOWER_THAN_ELSE
%nonassoc ELSE

%left ','
%right MATH_SHORT '='
%left '|'
%left '&'
%nonassoc EQ_TOKEN NE_TOKEN
%nonassoc '<' '>' LEQ_TOKEN GEQ_TOKEN
%left '-' '+' '.'
%left '*' '/'
%right '!' UMINUS
%right EXPO
%nonassoc PREFIX POSTFIX 
  
%%

control_file: statement
     | directive
     | control_file statement
     ;

directive: include_com
     ;

statement: ind_statement {copy_filepos(&handle->start_pos,&handle->begin_tok);}
     | '{' statement_block '}'
     ;

statement_block: ind_statement {copy_filepos(&handle->start_pos,&handle->begin_tok);}
     | statement_block ind_statement {copy_filepos(&handle->start_pos,&handle->begin_tok);}
     ;

ind_statement: ctypeI
     | ctypeII
     | ctypeIII
     | ctypeIV
     | ctypeV
     | ctype_bad
     | for_statement
     | if_statement
     | do_statement
     | while_statement
     | ';'
     ;

/* General commands */
ctypeI: WORD imm_expr {com_ctypeI($1,0,$2,0);}
     | WORD imm_expr ',' var_list {com_ctypeI($1,0,$2,$4);}
     | WORD var_list {com_ctypeI($1,0,0,$2);}
     | WORD clause imm_expr {com_ctypeI($1,$2,$3,0);}
     | WORD clause imm_expr ',' var_list {com_ctypeI($1,$2,$3,$5);}
     | WORD clause var_list {com_ctypeI($1,$2,0,$3);}
     ;

/* Position/Frequency commands */
ctypeII: WORD WORD tok_list
     | KEYWORD WORD WORD tok_list
     ;

assignment: var1 '=' gen_expr {if($1 && $3) com_assign($1,$3);}
     ;

/* Assignment or model commands */
ctypeIII: WORD WORD '=' gen_expr {}
     | WORD WORD_B array_index ')' '=' gen_expr {}
     | assignment
     | assignment ',' gen_expr_list {free_parse_term($3);}
     | var1 POSTFIX {free_parse_term(parse_incr_decr($1,PREFLAG|($2)));}
     | PREFIX var1 {free_parse_term(parse_incr_decr($2,PREFLAG|($1)));}
     | var1 MATH_SHORT gen_expr {vtmp=do_op(copy_parse_term($1),$3,$2);
		                           if($1 && vtmp) com_assign($1,vtmp);} 
     ;

/* Restrict/affected type commands */
ctypeIV: WORD KEYWORD_B gen_expr ')'
     | KEYWORD var_list KEYWORD '(' gen_expr ')'
     | KEYWORD KEYWORD_B gen_expr ')'
	  | KEYWORD_B gen_expr ')' KEYWORD var_list
     ;

/* Locus type commands */
ctypeV: WORD KEYWORD var_list
     ;

/* Some common errors */
ctype_bad: KEYWORD '=' {parerror1(key_err,key_tab[$1]);} gen_expr
     | KEYWORD_B '=' {parerror1(key_err,key_tab[$1]);} gen_expr
     | WORD KEYWORD '=' {parerror1(key_err,key_tab[$2]);} gen_expr
     | WORD KEYWORD_B '=' {parerror1(key_err,key_tab[$2]);} gen_expr
     ;

format_atom: CLAUSE_INT {$$=make_format_clause($1,0);}
     | CLAUSE_X {$$=make_format_clause(-1,0);}
     | CLAUSE_INT CLAUSE_X {} {$$=make_format_clause(-$1,0);}
     | CLAUSE_INT '(' format_clause ')' {$$=make_format_clause($1,reverse_list($3));}
     | '(' format_clause ')' {$$=$2;}
     ;

format_clause: format_atom
     | format_clause ',' format_atom {$3->next=$1; $$=$3;}
     ;

clause_elem: gen_expr_list1 {$$=make_clause($1);}
     | format_clause {$$=convert_format_clause(reverse_list($1));}
     ;

clause_content: clause_elem 
     | clause_content ';' clause_elem {$3->next=$1; $$=$3;}
     ;

clause: '[' clause_content ']' {$$=reverse_list($2);}
     | '[' ']' {$$=0;}
     ;

for_statement: FOR '(' opt_expr ';' opt_expr ';' opt_expr ')' statement
     ;

if_statement: IF '(' gen_expr ')' statement %prec LOWER_THAN_ELSE
     | IF '(' gen_expr ')' statement ELSE statement
     ;

do_statement: DO statement WHILE '(' gen_expr ')'
     ;

while_statement: WHILE '(' gen_expr ')' statement
     ;

opt_expr: /* Empty */
     | gen_expr
     ;

include_com: INCLUDE str_expr DUMMY {include_file($2);}
     ;

tok_list: tok
     | ',' tok {$$=addto_var_list(make_var_list(0,0),$2);}
     | poly_commas tok {$$=addto_var_list(addto_var_list($1,make_var_list(0,0)),$2);}
     | tok_list ',' tok {$$=addto_var_list($1,$3);}
     | tok_list poly_commas tok {$$=addto_var_list(addto_var_list($1,$2),$3);}
     ;

tok: imm_element {$$=make_var_list($1,0); copy_filepos(&($$->pos),&handle->begin_tok);}
     ;

array_index_elem: gen_expr
     | gen_expr DOUBLE_DOT gen_expr
     | DOUBLE_DOT gen_expr {$$=$2;}
     | gen_expr DOUBLE_DOT
     ;

array_index: array_index_elem
     | array_index ',' array_index_elem 
     ;

var1: WORD {$$=get_var_term($1);}
     | WORD_B array_index ')' {$$=get_array_term($1,$2);}
     ;

var2: WORD {$$=make_var_list(get_var_term($1),0);}
     | WORD_B array_index ')' {$$=get_array_var_list($1,$2);}
     | /* Empty */ {$$=make_var_list(0,0);}
     ;

loop_clause: '(' array_var_list ',' WORD '=' gen_expr ')' {$$=make_loop_clause($2,$4,$6);}
     | '(' array_var_list ')' {$$=make_loop_clause($2,0,0);}
     ; 

array_var_list: var2
     | array_var_list ',' var2 {$$=addto_var_list($1,$3);}
     | KEYWORD_B gen_expr ')' {parerror1(key_err,key_tab[$1]); $$=0;}
     | KEYWORD {parerror1(key_err,key_tab[$1]); $$=0;}
     | array_var_list ',' KEYWORD_B gen_expr ')' {parerror1(key_err,key_tab[$3]); $$=0;}
     | array_var_list ',' KEYWORD {parerror1(key_err,key_tab[$3]); $$=0;}
     ;

var: var1 {$$=make_var_list($1,0); copy_filepos(&($$->pos),&handle->begin_tok);}
     | var1 clause {$$=make_var_list($1,$2);}
     ;

var_list1: var
     | loop_clause
     ;

poly_commas: ',' ',' {$$=make_var_list(0,0);}
     | poly_commas ',' {$$=addto_var_list($1,make_var_list(0,0));}
     ;

var_list: var_list1
     | ',' var_list1 {$$=addto_var_list(make_var_list(0,0),$2);}
     | poly_commas var_list1 {$$=addto_var_list(addto_var_list($1,make_var_list(0,0)),$2);}
     | var_list ',' var_list1 {$$=addto_var_list($1,$3);}
     | var_list poly_commas var_list1 {$$=addto_var_list(addto_var_list($1,$2),$3);}
     ;

imm_expr: imm_element
     | imm_element1 '+' gen_expr {$$=do_op1($1,$3,'+');}
     | imm_element1 '.' gen_expr {$$=do_op1($1,$3,'.');}
     | imm_element1 '-' gen_expr {$$=do_op1($1,$3,'-');}
     | imm_element1 '/' gen_expr {$$=do_op1($1,$3,'/');}
     | imm_element1 '*' gen_expr {$$=do_op1($1,$3,'*');}
     | imm_element1 '<' gen_expr {$$=do_op1($1,$3,'<');}
     | imm_element1 '>' gen_expr {$$=do_op1($1,$3,'>');}
     | imm_element1 '&' gen_expr {$$=do_op1($1,$3,'&');}
     | imm_element1 '|' gen_expr {$$=do_op1($1,$3,'|');}
     | imm_element1 EQ_TOKEN gen_expr {$$=do_op1($1,$3,EQ_TOKEN);}
     | imm_element1 NE_TOKEN gen_expr {$$=do_op1($1,$3,NE_TOKEN);}
     | imm_element1 LEQ_TOKEN gen_expr {$$=do_op1($1,$3,LEQ_TOKEN);}
     | imm_element1 GEQ_TOKEN gen_expr {$$=do_op1($1,$3,GEQ_TOKEN);}
     | imm_element1 EXPO gen_expr {$$=do_op1($1,$3,EXPO);}
     | '-' gen_expr %prec UMINUS {$$=do_op1($2,0,'-');}
     | '!' gen_expr %prec UMINUS {$$=do_op1($2,0,'!');}
     | '(' imm_expr ')' {$$=check_immediate($2);}
     | FUNCTION '(' imm_expr ')' {$$=do_func($1,$3);}
     ;

str_op: '+'
     | '.'
     ;

str_expr: str_elem
     | str_expr str_op str_elem1 {$$=add_strings($1,$3);}
     | var1 str_op str_elem1 {$$=add_strings(check_string($1),$3);}
     ;

str_elem: STRING
     | FUNCTION '(' str_expr ')' {$$=$1->f.f_string($3);}
     ;
  
str_elem1: str_elem
     | var1 {$$=check_string($1);} 
     ;
  
gen_expr_list: gen_expr {if($1->defer) do_deferred($1); $$=$1; }
     | gen_expr_list ',' gen_expr {if($3->defer) do_deferred($3); $$=do_op($1,$3,',');}
     ;

gen_expr_list1: gen_expr
     | gen_expr_list1 ',' gen_expr {$$=do_op($1,$3,SUB_EXPR);}
     ;

gen_expr: gen_expr '+' gen_expr {$$=do_op($1,$3,'+');}
     | gen_expr '-' gen_expr {$$=do_op($1,$3,'-');}
     | gen_expr '.' gen_expr {$$=do_op($1,$3,'.');}
     | gen_expr '*' gen_expr {$$=do_op($1,$3,'*');}
     | gen_expr '/' gen_expr {$$=do_op($1,$3,'/');}
     | gen_expr '=' gen_expr {$$=do_op($1,$3,'=');}
     | gen_expr EQ_TOKEN gen_expr {$$=do_op($1,$3,EQ_TOKEN);}
     | gen_expr '<' gen_expr {$$=do_op($1,$3,'<');}
     | gen_expr '>' gen_expr {$$=do_op($1,$3,'>');}
     | gen_expr '|' gen_expr {$$=do_op($1,$3,'|');}
     | gen_expr '&' gen_expr {$$=do_op($1,$3,'&');}
     | gen_expr EXPO gen_expr {$$=do_op($1,$3,EXPO);}
     | gen_expr LEQ_TOKEN gen_expr {$$=do_op($1,$3,LEQ_TOKEN);}
     | gen_expr GEQ_TOKEN gen_expr {$$=do_op($1,$3,GEQ_TOKEN);}
     | gen_expr NE_TOKEN gen_expr {$$=do_op($1,$3,NE_TOKEN);}
     | gen_expr MATH_SHORT gen_expr {
		  if(IS_VAR($1)) {
			  vtmp=do_op(copy_parse_term($1),$3,$2); if($1 && vtmp) com_assign(copy_parse_term($1),vtmp); $$=$1;
		  } else {
			  parerror1("Can't assign to constant: '=' used instead of '==' perhaps?\n");
			  free_parse_term($3);
			  $$=$1;
		  } }
     | '-' gen_expr %prec UMINUS {$$=do_op($2,0,'-');}
     | '!' gen_expr %prec UMINUS {$$=do_op($2,0,'!');}
     | '(' gen_expr_list ')' {$$=$2;}
     | gen_element
     ;

gen_element: imm_element
     | KEYWORD {parerror1(key_err,key_tab[$1]);$$=0;}
     | KEYWORD_B {parerror1(key_err,key_tab[$1]);$$=0;}
     | var1
     | var1 POSTFIX {$$=parse_incr_decr($1,POSTFLAG|($2));}
     | PREFIX var1 {$$=parse_incr_decr($2,PREFLAG|($1));}
     | FUNCTION '(' gen_expr ')' {$$=do_func($1,$3);}
     ;

imm_element: INTEGER {$$=parse_make_texp($1,0,0,INTEGER);}
     | REAL {$$=parse_make_texp(0,$1,0,REAL);}
     | STRING {$$=parse_make_texp(0,0,$1,STRING);}
     ;

imm_element1: imm_element
     | var1 {$$=check_immediate($1);}
     ;
%%

static void set_filepos(filepos *fp)
{
	fp->fname=lex_info->name;
	fp->line=lex_info->line_no;
	fp->col=lex_info->col_no;
	if(!fp->col) fp->col++;
}

static void copy_filepos(filepos *fp1,filepos *fp2)
{
	memcpy(fp1,fp2,sizeof(filepos));
}

#ifdef FUNC_NAME
#undef FUNC_NAME
#endif
#define FUNC_NAME "include_file"
static void include_file(string *s)
{
	struct lex_info *li;
	struct stat st1,st2;
	
	if(s) {
		/* Check that we will not end up in a loop... */
		if(stat(s->s,&st1)) {
			fprintf(stderr,"File Error.  Couldn't stat include file '%s' for input: ",s->s);
			perror(0);
			ABT_FUNC(AbMsg);
		}
		li=lex_info;
		while(li) {
			if(stat(li->name,&st2)) {
				fprintf(stderr,"File Error.  Couldn't stat include file '%s': ",li->name);
				perror(0);
				ABT_FUNC(AbMsg);
			}
			if(st1.st_ino==st2.st_ino && st1.st_dev==st2.st_dev) break;
			li=li->last;
		}
		if(li) {
			fprintf(stderr,"Recursive include files - '%s'\n",s->s);
			li=lex_info;
			while(li) {
				if(!stat(li->name,&st2) && (st1.st_ino==st2.st_ino && st1.st_dev==st2.st_dev))
				  fprintf(stderr,"         -> included from '%s' <-\n",li->name);
				else fprintf(stderr,"            included from '%s'\n",li->name);
				li=li->last;
			}
			ABT_FUNC(AbMsg);
		}
		li=lex_info->next;
		if(!li) {
			if(!(li=malloc(sizeof(struct lex_info)))) ABT_FUNC(MMsg);
			li->next=0;
			li->last=lex_info;
			lex_info->next=li;
		}
		li->name=copy_cstring(s);
		free_string(s);
		li->ptr=li->ptr1=li->col_no=li->col_no1=li->line_no=li->eof=li->last_token=li->flag=0;
		if(!(li->fptr=fopen(li->name,"r"))) abt(__FILE__,__LINE__,"%s(): File Error.  Couldn't open '%s' for input\n",FUNC_NAME,li->name);
		printf("-->Moving to include file '%s'\n",li->name);
		lex_info=li;
		handle->fname=lex_info->name;
		handle->line=&lex_info->line_no;
	}
}
  
static int read_buf(void) 
{
	int s;
	
	if(lex_info->eof) return 0;
	s=(int)fread(lex_info->buf+lex_info->ptr1,1,BUFSIZE-lex_info->ptr1,lex_info->fptr);
	lex_info->ptr1+=s;
	lex_info->ptr=0;
	if(feof(lex_info->fptr)) lex_info->eof=1;
	return s;
}

static int get_char(int flag) 
{
	int i;
	
	if(lex_info->flag) {
		lex_info->flag=0;
		while(lex_info) {
			fclose(lex_info->fptr); /* Yes, close file and move back */
			free(lex_info->name);
			lex_info=lex_info->last;
			printf("--> Moving back to include file '%s'\n",lex_info->name);
			handle->fname=lex_info->name;
			handle->line=&lex_info->line_no;
			if(lex_info->last_token) break;
		}
	}
	/* Check if we have enough characters in the buffer */
	while(((i=lex_info->ptr1-lex_info->ptr))<(flag?1:LOOKAHEAD)) {
#ifdef DEBUG
		if(i<0) ABT_FUNC("Buffer overrun...\n");
#endif
		if(flag) return -1;
		if(!(lex_info->eof)) {
			/* Shift unread part of buffer to beginning */
			if(lex_info->ptr && i) {
				memmove(lex_info->buf,lex_info->buf+lex_info->ptr,i);
				lex_info->ptr1=i;
				lex_info->ptr=0;
			}
			/* Read in more characters (if possible) */
			if(read_buf()) break; /* OK, refilled buffer */
		}
		/* Can't read in more, are we at the end of the buffer ? */
		if(!i) {
			/* Yes, is there a parent file to return to? */
			if(lex_info->last) {
				i=-1;
				lex_info->flag=1;
				break;
			} else {
				lex_info->buf[lex_info->ptr1]=0; 
				break;
			}
		} else break;
	}
	if(i>=0) {
		/* Get next character */
		if(lex_info->ptr<lex_info->ptr1) {
			i=lex_info->buf[lex_info->ptr++];
			lex_info->col_no++;
		} else i=EOF;
		/* Handle DOS or Mac line ending sequences */
		if(i=='\r') {
			if(lex_info->ptr<lex_info->ptr1 && lex_info->buf[lex_info->ptr]=='\n') 
			  lex_info->ptr++;
			i='\n';
		}
		if(i=='\n') {
			lex_info->line_no++;
			lex_info->col_no1=lex_info->col_no;
			lex_info->col_no=0;
		}
	} else i=0;
	return i;
}

#ifdef FUNC_NAME
#undef FUNC_NAME
#endif
#define FUNC_NAME "alloc_var"
static struct bin_node *alloc_var(string *s)
{
	struct bin_node *node;
	struct parse_var *vv;
	
	if(!(node=malloc(sizeof(struct bin_node)))) ABT_FUNC(MMsg);
	node->left=node->right=0;
	node->balance=0;
	if(!(vv=malloc(sizeof(struct parse_var)))) ABT_FUNC(MMsg);
	node->data=vv;
	vv->name=s;
	vv->type=VARIABLE;
	vv->size=0;
	return node;
}

static struct bin_node *find_var(string *s,struct bin_node *node,struct bin_node **node1,int *balanced)
{
	int i;
	struct parse_var *vv;
	
	vv=node->data;
	if((i=strcasecmp(s->s,vv->name->s))) {
		if(i<0) {
			if(node->left) {
				node->left=find_var(s,node->left,node1,balanced);
			} else {
				*node1=node->left=alloc_var(s);
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
				node->right=find_var(s,node->right,node1,balanced);
			} else {
				*node1=node->right=alloc_var(s);
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
		free_string(s);
	}
	return node;
}

static int check_token(string *s,yystype *lval,char c)
{
	int i,j;
	static char *builtin[]={"FOR","IF","ELSE","WHILE","DO","INCLUDE","ALIAS","AND","OR","NOT","EQ","NE",0};
	int builtin_tok[]={FOR,IF,ELSE,WHILE,DO,INCLUDE,ALIAS,'&','|','!',EQ_TOKEN,NE_TOKEN,0};
	struct bin_node *node;
	func_def *fl;
	
	i=-1;
	while(builtin[++i]) if(!strcasecmp(s->s,builtin[i])) break;
	j=builtin_tok[i];;
	if(!j && key_tab) {
		i=-1;
		while(key_tab[++i]) if(!strcasecmp(s->s,key_tab[i])) break;
		if(key_tab[i]) {
			j=KEYWORD;
			lval->i=i;
			free_string(s);
		}
	}
	if(!j && func_list && c=='(') {
		fl=func_list;
		i=0;
		while(fl) {
			if(!strcasecmp(s->s,fl->com)) break;
			fl=fl->next;
			i++;
		}
		if(fl) {
			j=FUNCTION;
			lval->f=fl;
			free_string(s);
		}
	}
	if(!j) {
		if(!root_var) node=root_var=alloc_var(s);
		else root_var=find_var(s,root_var,&node,&i);
		j=WORD;
		lval->vv=(struct parse_var *)node->data;
	}
	return j;
}

static void copy_lval(yystype *s1,yystype *s2,int tok)
{
	switch(tok) {
	 case INTEGER:
		s1->i=s2->i;
		break;
	 case REAL:
		s1->x=s2->x;
		break;
	 case WORD:
		s1->vv=s2->vv;
		break;
	 case STRING:
		s1->str=copy_string(s2->str);
		break;
	}
}

#ifdef FUNC_NAME
#undef FUNC_NAME
#endif
#define FUNC_NAME "parlex"
static int parlex(yystype *lval)
{
	string *s=0;
	int c,c1,c2,state,res,pt1;
	static int inc_state,postfix_flag,in_clause;
	struct lex_info *li;
	filepos fp;
	
	for(;;) {
		c=get_char(0);
		if(c==EOF) {
			c=0;
			break;
		}
		if(isspace((int)c)) {
			postfix_flag=0;
			continue;
		}
		/* Remove comments */
		c1=lex_info->buf[lex_info->ptr];
		if(c=='#' || (c=='/' && c1=='/')) {
			do c=get_char(0); while(c && c!='\n');
		} else if(c=='/' && c1=='*') {
			do c=get_char(0); while (c!='*' || lex_info->buf[lex_info->ptr]!='/');
			(void)get_char(0);
		} else break;
	}
	res=0;
	set_filepos(&fp);
	/* Strings */
	if(c=='\"' || c=='\'') {
		for(;;) {
			c1=get_char(0);
			if(c1==EOF) c1=0;
			if(c1==c || !c1 || c1=='\n') break;
			if(c1=='\\') { /* Check for escape sequences */
				c1=get_char(0);
				switch(c1) {
				 case 'a':
					c1='\a';
					break;
				 case 'b':
					c1='\b';
					break;
				 case 'f':
					c1='\f';
					break;
				 case 'n':
					c1='\n';
					break;
				 case 'r':
					c1='\r';
					break;
				 case 't':
					c1='\t';
					break;
				 case 'v':
					c1='\v';
					break;
				}
			}
			s=add_to_string(s,c1);
		}
		if(c1!=c) parerror1("Unterminated string\n");
		c=res=STRING;
		lval->str=s;
	} else if(c=='[' && !in_clause) { /* Format Clauses */
		in_clause=1;
	} else if(c==']' && in_clause) {
		in_clause=0;
	}
	c1=lex_tab[c];
	if(!res) {
		if(in_clause && (c1&FORMAT_BITS)) {
			if(c1&DIGIT_BIT) {
				pt1=lex_info->ptr-1;
				do {
					c=get_char(1);
					if(c<0) ABT_FUNC("Parse error - number too long\n");
					c1=lex_tab[c];
				} while(c1&DIGIT_BIT);
				res=CLAUSE_INT;
				lval->i=atoi(lex_info->buf+pt1);
			} else if(c1&X_BIT) {
				c2=lex_info->buf[lex_info->ptr];
				if(!(lex_tab[c2]&(TOKEN_BITS|REL_BIT|MATH_BIT|LBRACK_BIT))) {
					res=CLAUSE_X;
					c=get_char(0);
				}
			}
		}
		if(!res && (c1&TOKEN_BITS)) {
			if(c1&(ALPHA_BIT|USCORE_BIT)) state=1;
			else if(c1&DIGIT_BIT) state=2;
			else if((c1&PERIOD_BIT)&&!postfix_flag) state=7;
			else state=0;
			/* Starting position of token */
			pt1=lex_info->ptr-1;
			while(state>0) {
				c=get_char(1);
				if(c<0) ABT_FUNC("Parse error - token too long\n");
				c1=lex_tab[c];
				switch(state) {
				 case 1: /* WORD */
					if(!(c1&(ALPHA_BIT|USCORE_BIT|DIGIT_BIT))) {
						res=WORD;
						state=0;
					}
					break;
				 case 2: /* INTEGER */
					if(c1&E_BIT) {
						res=INTEGER;
						state=4;
					} else if(c1&PERIOD_BIT) state=3;
					else if(!(c1&DIGIT_BIT)) {
						res=INTEGER;
						state=0;
					}
					break;
				 case 3: /* REAL */
					if(c1&E_BIT) {
						res=REAL;
						state=4;
					} else if(!(c1&DIGIT_BIT)) {
						res=REAL;
						state=0;
					}
					break;
				 case 4: /* REAL (exponential) */
					if(c1&DIGIT_BIT) state=5;
					else if(c1&PM_BIT) state=6;
					else {
						lex_info->ptr--;
						if(lex_info->col_no) lex_info->col_no--;
						else {
							lex_info->line_no--;
							lex_info->col_no=lex_info->col_no1;
						}
						state=0;
					}
					break;
				 case 5: /* REAL (exp. after digit) */
					if(!(c1&DIGIT_BIT)) {
						res=REAL;
						state=0;
					}
					break;
				 case 6: /* REAL (exp. after +-) */
					if(c1&DIGIT_BIT) {
						res=REAL;
						state=4;
					} else {
						lex_info->ptr-=2;
						state=0;
					}
					break;
				 case 7: /* After initial . */
					if(c1&DIGIT_BIT) {
						res=REAL;
						state=3;
					} else {
						state=0;
						res='.';
					}
				}
			}
		}
	}
	if(res) {
		copy_filepos(&handle->begin_tok,&fp);
	}
	if(res && res!=STRING) {
		lex_info->ptr--;
		if(lex_info->col_no) lex_info->col_no--;
		else {
			lex_info->line_no--;
			lex_info->col_no=lex_info->col_no1;
		}
		switch(res) {
		 case WORD:
			s=addn_to_string(s,lex_info->buf+pt1,lex_info->ptr-pt1);
			for(;;) {
				c1=lex_info->buf[lex_info->ptr];
				if(c1 && isspace(c1)) (void)get_char(0);
				else break;
			}
			res=check_token(s,lval,c1);
			if(res==WORD || res==KEYWORD) {
				if(c1=='(') {
					res=(res==WORD?WORD_B:KEYWORD_B);
					lex_info->col_no++;
					lex_info->ptr++;
				}
			}
			break;
		 case INTEGER:
			lval->i=atoi(lex_info->buf+pt1);
			break;
		 case REAL:
			lval->x=atof(lex_info->buf+pt1);
			break;
		}
		c=res;
	} else if(c1&REL_BIT) {
		switch(c) {
		 case '<':
			c1=lex_info->buf[lex_info->ptr];
			if(c1=='=') {
				c=LEQ_TOKEN;
				lex_info->ptr++;
				lex_info->col_no++;
			} else if(c1=='>') {
				c=NE_TOKEN;
				lex_info->ptr++;
				lex_info->col_no++;
			}
			break;
		 case '>':
			if(lex_info->buf[lex_info->ptr]=='=') {
				c=GEQ_TOKEN;
				lex_info->ptr++;
				lex_info->col_no++;
			}
			break;
		 case '=':
			if(lex_info->buf[lex_info->ptr]=='=') {
				lex_info->ptr++;
				lex_info->col_no++;
				c=EQ_TOKEN;
			}
			break;
			/* Handle doubled symbols (no special meaning...) */
		 case '|':
			if(lex_info->buf[lex_info->ptr]=='|') {
				lex_info->ptr++;
				lex_info->col_no++;
			}
			break;
		 case '&':
			if(lex_info->buf[lex_info->ptr]=='&') {
				lex_info->ptr++;
				lex_info->col_no++;
			}
			break;
		 case '!':
			if(lex_info->buf[lex_info->ptr]=='=') {
				c=NE_TOKEN;
				lex_info->ptr++;
				lex_info->col_no++;
			}
		}
	} else if(c1&MATH_BIT) {
		c1=lex_info->buf[lex_info->ptr];
		if(c1==c) {
			if(c=='+' || c=='-') {
				lex_info->ptr++;
				lex_info->col_no++;
				lval->i=c;
				c=(postfix_flag?POSTFIX:PREFIX);
			} else if(c=='.') {
				c=DOUBLE_DOT;
				lex_info->ptr++;
				lex_info->col_no++;
			} else if(c=='*') {
				c=EXPO;
				lex_info->ptr++;
				lex_info->col_no++;
			}
		} else if(c1=='=') {
			lex_info->ptr++;
			lex_info->col_no++;
			lval->i=c;
			c=MATH_SHORT;
		}
	}
	if(!c && lex_info->flag) {
		li=lex_info->last;
		while(li) {
			c=li->last_token;
			if(c) {
				copy_lval(lval,&li->last_lval,c);
				break;
			}
			li=li->last;
		}
	}
	/* Kludge to avoid sending next token if we've finished an include statement 
	 * As this is a pure parse we can not use yyclearin to discard lookup token... */
	if(inc_state==1) {
		if(c==STRING || c==WORD) inc_state=2;
		else inc_state=3;
	} else if(inc_state==2) {
		if(c=='+') inc_state=1;
		else inc_state=3;
	} 
	if(inc_state==3) {
		if(!lex_info->flag) {
			/* Store last token, to be used when we return to this file */
			printf("storing token '%d' for file '%s'\n",c,lex_info->name);
			lex_info->last_token=c;
			copy_lval(&lex_info->last_lval,lval,c);
		} else lex_info->last_token=0;
		c=DUMMY;
		inc_state=0;
	}
	postfix_flag=(c==WORD || c==')')?1:0;
	return c;
}

static int parerror(char *s)
{
	int c;
	
	c=lex_info->col_no;
	if(!c) c++;
	if(s) fprintf(stderr,"%s, line %d column %d: %s\n",lex_info->name,lex_info->line_no+1,c,s);
	else fprintf(stderr,"%s, line %d column %d: ",lex_info->name,lex_info->line_no+1,c);
	return 0;
}

/* Report file position given by fp */
void parerror_fp(const filepos *fp)
{
	if(stdout) (void)fflush(stdout);
	if(fp->fname) fprintf(stderr,"%s, line %d column %d: ",fp->fname,fp->line+1,fp->col);
	else fputs("<Unknown position>: ",stderr);
}

/* Report file position at beginning of last token*/
void parerror1(const char *fmt, ...)
{
	va_list args;

	parerror_fp(&handle->begin_tok);
	va_start(args,fmt);
	(void)vfprintf(stderr,fmt,args);
	va_end(args);
}

/* Report file position at beginning of last command */
void parerror2(const char *fmt, ...)
{
	va_list args;

	parerror_fp(&handle->start_pos);
	va_start(args,fmt);
	(void)vfprintf(stderr,fmt,args);
	va_end(args);
}

void parerror3(filepos *fp,const char *fmt, ...)
{
	va_list args;

	parerror_fp(fp);
	va_start(args,fmt);
	(void)vfprintf(stderr,fmt,args);
	va_end(args);
}

/* Set up lexer table */
static void init_lex_tab(void) {
	int i,j;
	
	for(i=0;i<256;i++) {
		j=0;
		if(isalpha(i)) {
			j=ALPHA_BIT;
			if(i=='e' || i=='E') j|=E_BIT;
			else if(i=='x' || i=='X') j|=X_BIT;
		} else if(isdigit(i)) j=DIGIT_BIT;
		else switch(i) {
		 case '+':
		 case '-':
			j=MATH_BIT|PM_BIT;
			break;
		 case '_':
			j=USCORE_BIT;
			break;
		 case '.':
			j=MATH_BIT|PERIOD_BIT;
			break;
		 case ',':
			j=COMMA_BIT;
			break;
		 case '(':
			j=LBRACK_BIT;
			break;
		 case ')':
			j=RBRACK_BIT;
			break;
		 case '<':
		 case '>':
		 case '=':
		 case '!':
		 case '|':
		 case '&':
			j=REL_BIT;
			break;
		 case '/':
		 case '*':
			j=MATH_BIT;
			break;
		}
		lex_tab[i]=j;
	}
}

static void init_func_list(func_def *funcs)
{
	funcs=register_real_function(funcs,log,"log");
	funcs=register_real_function(funcs,log10,"log10");
	funcs=register_real_function(funcs,exp,"exp");
	funcs=register_real_function(funcs,sqrt,"sqrt");
	funcs=register_real_function(funcs,fabs,"abs");
	funcs=register_real_function(funcs,cbrt,"cbrt");
	funcs=register_real_function(funcs,sin,"sin");
	funcs=register_real_function(funcs,cos,"cos");
	funcs=register_real_function(funcs,tan,"tan");
	funcs=register_real_function(funcs,sinh,"sinh");
	funcs=register_real_function(funcs,cosh,"cosh");
	funcs=register_real_function(funcs,tanh,"tanh");
	funcs=register_real_function(funcs,asin,"asin");
	funcs=register_real_function(funcs,acos,"acos");
	funcs=register_real_function(funcs,atan,"atan");
	func_list=funcs;
}

static void free_vars(struct bin_node *node)
{
	struct parse_var *v;
	
	if(node->left) free_vars(node->left);
	if(node->right) free_vars(node->right);
	v=node->data;
	free_parse_var(v);
	free(node);
}

static void free_parse_stuff(void)
{
	if(root_var) free_vars(root_var);
	free(lex_info);
}

#ifdef FUNC_NAME
#undef FUNC_NAME
#endif
#define FUNC_NAME "Parse_File"
int Parse_File(FILE *fptr,char **keywords,func_def *funcs,parse_handler *hand)
{
	int err;
	
#if YYDEBUG
	pardebug=1;
#endif
	handle=hand;
	init_parse_utils(hand);
	if(!(lex_info=calloc((size_t)1,sizeof(struct lex_info)))) ABT_FUNC(MMsg);
	lex_info->next=lex_info->last=0;
	lex_info->fptr=fptr;
	handle->line=&lex_info->line_no;
	lex_info->name=handle->fname;
	handle->start_pos.fname=handle->fname;
	handle->start_pos.col=1;
	handle->start_pos.line=1;
	init_lex_tab();
	key_tab=keywords;
	init_func_list(funcs);
	atexit(free_parse_stuff);
	err=parparse();
	return err;
}
