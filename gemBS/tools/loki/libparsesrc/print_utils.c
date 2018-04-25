/****************************************************************************
 *                                                                          *
 *     Loki - Programs for genetic analysis of complex traits using MCMC    *
 *                                                                          *
 *                     Simon Heath - CNG, Evry                              *
 *                                                                          *
 *                           July 2003                                      *
 *                                                                          *
 * print_utils.c:                                                           *
 *                                                                          *
 * Copyright (C) Simon C. Heath 2003                                        *
 * This is free software.  You can distribute it and/or modify it           *
 * under the terms of the Modified BSD license, see the file COPYING        *
 *                                                                          *
 ****************************************************************************/
#include <config.h>
#include <stdlib.h>
#include <string.h>
#ifdef USE_DMALLOC
#include <dmalloc.h>
#endif
#include <stdio.h>
#include <math.h>
	
#include "parser.h"
#include "parse.tab.h"

void print_parse_op(FILE *fptr,int op)
{
	string *s;

	s=parse_op_str(op);
	fputs(get_cstring(s),fptr);
	free_string(s);
}

string *parse_op_str(int op)
{
	string *s;
	
	switch(op) {
	 case '.':
	 case '+':
	 case '-':
	 case '*':
	 case '/':
	 case '<':
	 case '>':
	 case '=':
	 case ',':
	 case '!':
	 case '&':
	 case '|':
		s=add_to_string(0,op);
		break;
	 case LEQ_TOKEN:
		s=addn_to_string(0,"<=",2);
		break;
	 case EXPO:
		s=addn_to_string(0,"**",2);
		break;
	 case GEQ_TOKEN:
		s=addn_to_string(0,">=",2);
		break;
	 case NE_TOKEN:
		s=addn_to_string(0,"!=",2);
		break;
	 case SUB_EXPR:
		s=addn_to_string(0,",",1);
		break;
	 default:
		s=add_to_string(0,'?');
	}
	return s;
}

void print_parse_exp(FILE *fptr,struct parse_term *ex)
{
	string *s;
	
	s=parse_exp_str(ex);
	if(s) {
		fputs(get_cstring(s),fptr);
		free_string(s);
	}
}

char *parse_exp_cstr(struct parse_term *ex)
{
	string *s;
	char *p=0;
	
	s=parse_exp_str(ex);
	if(s) p=extract_cstring(s);
	return p;
}

string *convert_to_string(struct parse_term *ex)
{
	int type;
	
	string *s=0;
	
	type=ex->type&VT_TYPES;
	switch(type) {
	 case INTEGER:
		s=int_to_string(ex->elem.i);
		break;
	 case STRING:
		s=copy_string(ex->elem.str);
		break;
	 case REAL:
		s=real_to_string(ex->elem.x);
		break;
	 case INF:
		s=addn_to_string(0,"<INF>",5);
		break;
	 case XNAN:
		s=addn_to_string(0,"<NAN>",5);
		break;
	}
	return s;
}

static string *parse_exp_str1(struct parse_term *ex,int outside)
{
	struct expr_op *eop;
	struct func_op *fop;
	int type;
	string *s=0;
	
	if(!ex) return s;
	type=ex->type&VT_TYPES;
	switch(type) {
	 case INTEGER:
		s=int_to_string(ex->elem.i);
		break;
	 case STRING:
		s=add_to_string(0,'\"');
		s=add_strings(s,copy_string(ex->elem.str));
		s=add_to_string(s,'\"');
		break;
	 case REAL:
		s=real_to_string(ex->elem.x);
		break;
	 case VARIABLE:
		s=parse_term_name(ex);
		break;
	 case INF:
		s=addn_to_string(0,"<INF>",5);
		break;
	 case XNAN:
		s=addn_to_string(0,"<NAN>",5);
		break;
	 case FUNCTION:
		fop=ex->elem.fop;
		s=addn_to_string(0,fop->f->com,strlen(fop->f->com));
		add_to_string(s,'(');
		s=add_strings(s,parse_exp_str(fop->ex));
		add_to_string(s,')');
		break;
	 case EXPR_OP:
		eop=ex->elem.eop;
		if(eop->op&(PREFLAG|POSTFLAG)) {
			s=parse_term_name(eop->exp1);
			if(eop->op&PREFLAG) {
				s=paddn_to_string(s,eop->op==PRE_INCR?"++":"--",2);
			} else if(eop->op&POSTFLAG) {
				s=addn_to_string(s,eop->op==POST_INCR?"++":"--",2);
			}
		} else {
			if(!eop->exp2) { /* Unary operation */
				s=parse_op_str(eop->op);
				s=add_strings(s,parse_exp_str(eop->exp1));
			} else { /* Binary operation */
				if(eop->op==SUB_EXPR) {
					s=add_strings(s,parse_exp_str1(eop->exp1,1));
					s=add_strings(s,parse_op_str(eop->op));
					s=add_strings(s,parse_exp_str1(eop->exp2,1));
				} else {
					if(!outside) s=add_to_string(0,'(');
					s=add_strings(s,parse_exp_str1(eop->exp1,0));
					s=add_strings(s,parse_op_str(eop->op));
					s=add_strings(s,parse_exp_str1(eop->exp2,0));
					if(!outside) s=add_to_string(s,')');
				}
			}
			break;
		}
	}
	return s;
}

string *parse_exp_str(struct parse_term *ex)
{
	return parse_exp_str1(ex,1);
}

void print_var_list(FILE *fptr,struct parse_var_list *vl)
{
	string *s;
	
	s=var_list_string(vl);
	fputs(get_cstring(s),fptr);
	free_string(s);
}

string *var_list_string(struct parse_var_list *vl)
{
	int i=0;
	string *s=0;
	
	while(vl) {
		if(i++) {
			s=add_to_string(s,',');
			s=add_strings(s,parse_exp_str(&vl->term));
		} else s=parse_exp_str(&vl->term);
		if(vl->clause) {
			s=add_to_string(s,'[');
			s=add_strings(s,parse_clause_str(vl->clause));
			s=add_to_string(s,']');
		}
		vl=vl->next;
	}
	return s;
}

void print_parse_term_name(FILE *fptr,struct parse_term *v)
{
	string *s;
	
	s=parse_term_name(v);
	if(s) {
		fputs(get_cstring(s),fptr);
		free_string(s);
	}
}

string *parse_term_name(struct parse_term *v)
{
	string *s;
	int type;
	
	type=v->type&VT_TYPES;
	if(type==EMPTY_VAR) return 0;
	if(type==VARIABLE) {
		if(!v->elem.vv) s=addn_to_string(0,"<NULL>",6);
		else {
			if(v->type&VT_ARRAY_ELEMENT) {
				s=copy_string1(0,v->elem.vv1->head->name);
				s=add_to_string(s,'(');
				s=add_strings(s,parse_exp_str(v->elem.vv1->ex));
				s=add_to_string(s,')');
			} else s=copy_string(v->elem.vv->name);
		}
	} else s=addn_to_string(0,"<OOOK!>",7);
	return s;
}

string *format_clause_str(struct format_clause *fc)
{
	int flag=0;
	string *s=0;
	
	while(fc) {
		if(flag++) add_to_string(s,',');
		switch(fc->type) {
		 case FC_SKIP:
			if(fc->multiplier>1) s=add_strings(s,int_to_string(fc->multiplier));
			s=add_to_string(s,'X');
			break;
		 case FC_READ:
			s=add_strings(s,int_to_string(fc->multiplier));
			break;
		 case FC_SUBCLAUSE:
			s=add_strings(s,int_to_string(fc->multiplier));
			s=add_to_string(s,'(');
			s=add_strings(s,format_clause_str(fc->elem));
			s=add_to_string(s,')');
			break;
		}
		fc=fc->next;
	}
	return s;
}

void print_format_clause(FILE *fptr,struct format_clause *fc)
{
	string *s;
	
	s=format_clause_str(fc);
	if(s) {
		fputs(get_cstring(s),fptr);
		free_string(s);
	}
}

string *parse_clause_str(struct parse_clause *pc)
{
	int flag=0;
	string *s=0;
	
	while(pc) {
		if(flag++) add_to_string(s,';');
		switch(pc->type) {
		 case PC_TERM:
			s=add_strings(s,parse_exp_str(pc->clause.term));
			break;
		 case PC_FORMAT_CLAUSE:
			s=add_strings(s,format_clause_str(pc->clause.format));
			break;
		}
		pc=pc->next;
	}
	return s;
}

void print_parse_clause(FILE *fptr,struct parse_clause *pc)
{
	string *s;
	
	s=parse_clause_str(pc);
	if(s) {
		fputs(get_cstring(s),fptr);
		free_string(s);
	}
}
