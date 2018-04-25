/****************************************************************************
 *                                                                          *
 *     Loki - Programs for genetic analysis of complex traits using MCMC    *
 *                                                                          *
 *                     Simon Heath - CNG, Evry                              *
 *                                                                          *
 *                           June 2003                                      *
 *                                                                          *
 * parse_utils.c:                                                           *
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
#include <ctype.h>
#include <math.h>
#ifdef HAVE_IEEEFP_H
#include <ieeefp.h>
#endif
#include <assert.h>

#include "parser.h"
#include "parse.tab.h"
#include "lk_malloc.h"

static parse_handler *handle;
static int handling_clause;

static void assign_call(struct parse_term *v,struct parse_term *v1)
{
	struct deferred_assign *def;
	
	if(!handling_clause) {
		if(handle->assign) handle->assign(v,v1);
	} else {
		def=lk_malloc(sizeof(struct deferred_assign));
		def->next=v->clause_assign;
		v->clause_assign=0;
		def->var=copy_parse_term(v);
		v->clause_assign=def;
		def->expr=copy_parse_term(v1);
	}
}

void init_parse_utils(parse_handler *h)
{
	handle=h;
}

void do_deferred(struct parse_term *term)
{
	struct deferred_assign *def,*def1;
	struct parse_term *v,*v1;
	
	def=term->defer;
	while(def) {
		def1=def->next;
		v=def->var;
		v1=def->expr;
		assign_var(v,v1);
		assign_call(v,v1);
		if(v1->defer) do_deferred(v1);
		free_parse_term(v1);
		free_parse_term(v);
		free(def);
		def=def1;
	}
	term->defer=0;
}

void com_assign(struct parse_term *v,struct parse_term *t)
{
	t=try_resolve_expr(t);
	if(t) {
		assign_var(v,t);
		assign_call(v,t);
		if(t->defer) do_deferred(t);
		free_parse_term(t);
	}
	free_parse_term(v);
}

struct parse_term *parse_incr_decr(struct parse_term *v,int type)
{
	int j;
	struct parse_term *v1,*v2;
	struct deferred_assign *def;
	
	j=(v->type&VT_TYPES);
	assert(j==VARIABLE);
	v2=parse_make_texp(1,0,0,INTEGER);
	v1=do_op(copy_parse_term(v),v2,(type==POST_INCR || type==PRE_INCR)?'+':'-');
	v1=try_resolve_expr(v1);
	if(v1) {
		/* If postfix, defer assignment until later */
		if(type==POST_INCR || type==POST_DECR) {
			def=lk_malloc(sizeof(struct deferred_assign));
			def->next=v->defer;
			v->defer=def;
			def->var=copy_parse_term(v);
			def->expr=v1;
		} else { /* Whereas prefix we do immediately */
			assign_var(v,v1);
			assign_call(v,v1);
			if(v1->defer) do_deferred(v1);
			free_parse_term(v1);
		}
	}
	return v;
}

void handle_clause(struct parse_clause *c)
{
	handling_clause=1;
	while(c) {
		if(c->type==PC_TERM) {
			c->clause.term=try_resolve_expr(c->clause.term);
		}
		c=c->next;
	}
	handling_clause=0;
}

struct parse_var_list *merge_var_varlist(struct parse_var *v,struct parse_clause *c,struct parse_var_list *vl)
{
	struct parse_var_list *vl1;
	
	vl1=vl;
	while(vl1->next) vl1=vl1->next;
	vl1->next=make_var_list(get_var_term(v),c);
	return vl;
}

struct parse_term *check_immediate(struct parse_term *t)
{
	int type=0;
	
	if(t) {
		type=(t->type&VT_TYPES);
		if(type==VARIABLE || type==EXPR_OP) {
			t=try_resolve_expr(t);
			if(t) type=(t->type&VT_TYPES);
		}
	}
	if(t) {
		if(!(type==INTEGER || type==STRING || type==REAL || type==INF || type==XNAN)) {
			free_parse_term(t);
			t=0;
		}
	}
	return t;
}

void com_ctypeI(struct parse_var *v,struct parse_clause *c,struct parse_term *s,struct parse_var_list *vl)
{
	struct parse_var_list *vl1;
	int type;
	
	if(c) handle_clause(c);
	if(handle->ctypeI) {
		if(vl) {
			vl=reverse_list(vl);
			/* If arg not set, see if first item in list is an immediate var */
			if(!s && !vl->clause) {
 				type=vl->term.type&VT_TYPES;
				if(type==STRING || type==INTEGER || type==REAL) {
					s=copy_parse_term(&vl->term);
					vl1=vl->next;
					if(type==STRING) free_string(vl->term.elem.str);
					free(vl);
					vl=vl1;
				}
			}
			/* Check if remaining var_list contains just vectors */
			vl1=vl;
			while(vl1) {
 				type=vl1->term.type&VT_TYPES;
				if(type!=VARIABLE && type!=EMPTY_VAR) {
					parerror3(&vl1->pos,"Non-vector element in variable list\n");
					if(type==STRING) free_string(vl1->term.elem.str);
					vl1->term.type=EMPTY_VAR;
				}
				vl1=vl1->next;
			}
		}
		handle->ctypeI(v?get_cstring(v->name):0,c,s,vl);
	} else free_parse_var_list(vl);
	if(s) free_parse_term(s);
	if(c) free_parse_clause(c);
}

void com_ctypeII(int keyword,struct parse_var *com,struct parse_term *arg,struct parse_var_list *vl)
{
	if(handle->ctypeII) {
		if(vl) vl=reverse_list(vl);
		handle->ctypeII(keyword,get_cstring(com->name),arg,vl);
	} else free_parse_var_list(vl);
	if(arg) free_parse_term(arg);
}

void com_ctypeIII(struct parse_var *com,struct parse_term *v,struct parse_term *arg)
{
	if(handle->ctypeIII) {
		arg=try_resolve_expr(arg);
		handle->ctypeIII(get_cstring(com->name),v,arg);
	}
	if(v) free_parse_term(v);
	if(arg) free_parse_term(arg);
}

void com_ctypeV(struct parse_var *com,int keyword,struct parse_var_list *vl)
{
	if(handle->ctypeV) {
		if(vl) vl=reverse_list(vl);
		handle->ctypeV(get_cstring(com->name),keyword,vl);
	} else free_parse_var_list(vl);
}

struct parse_term *get_new_term(void) 
{
	struct parse_term *term;
	
	term=lk_malloc(sizeof(struct parse_term));
	term->type=0;
	term->defer=term->clause_assign=0;
	return term;
}

struct parse_var *find_parse_var(struct parse_term *tm)
{
	struct parse_var *v=0;
	
	if(tm) {
		if(tm->type&VT_ARRAY_ELEMENT) v=get_array_var(tm->elem.vv1->head,&tm->elem.vv1->ex);
		else v=tm->elem.vv;
	}
	return v;
}

int convert_numeric(struct parse_term *v)
{
	char *p;
	double x,i,fr;
	
	assert(v->type==STRING);
	x=strtod(get_cstring(v->elem.str),&p);
	free_string(v->elem.str);
	fr=modf(x,&i);
	if(fabs(fr)<WORKING_ZERO) {
		v->elem.i=(int)i;
		v->type=INTEGER;
	} else {
		v->elem.x=x;
		v->type=REAL;
	}
	return v->type;
}

struct parse_clause *convert_format_clause(struct format_clause *fc)
{
	struct parse_clause *pc;
	
	pc=lk_malloc(sizeof(struct parse_clause));
	pc->next=0;
	pc->type=PC_FORMAT_CLAUSE;
	pc->clause.format=fc;
	return pc;
}

struct parse_clause *make_clause(struct parse_term *ex)
{
	struct parse_clause *pc;
	
	pc=lk_malloc(sizeof(struct parse_clause));
	pc->next=0;
	pc->type=PC_TERM;
	pc->clause.term=ex;
	return pc;
}

struct format_clause *make_format_clause(int i,struct format_clause *fc)
{
	struct format_clause *fc1;
	
	fc1=lk_malloc(sizeof(struct format_clause));
	fc1->next=0;
	fc1->multiplier=abs(i);
	if((fc1->elem=fc)) fc1->type=FC_SUBCLAUSE;
	else fc1->type=i>0?FC_READ:FC_SKIP;
	return fc1;
}

static int convert_real(struct parse_term *v)
{
	assert(v->type==INTEGER);
	v->elem.x=(double)v->elem.i;
	return (v->type=REAL);
}

struct parse_term *parse_make_texp(int i,double x,string *s,int type)
{
	struct parse_term *expr;

	expr=get_new_term();
	switch(type) {
	 case INTEGER:
		expr->elem.i=i;
		break;
	 case REAL:
		expr->elem.x=x;
		break;
	 case STRING:
		expr->elem.str=s;
		break;
	}
	expr->type=type;
	return expr;
}

int convert_string(struct parse_term *t)
{
	switch(t->type) {
	 case INTEGER:
		t->elem.str=int_to_string(t->elem.i);
		break;
	 case REAL:
		t->elem.str=real_to_string(t->elem.x);
		break;
	 case INF:
		t->elem.str=addn_to_string(0,"<INF>",5);
		break;
	 case XNAN:
		t->elem.str=addn_to_string(0,"<XNAN>",5);
		break;
	}
	t->type=STRING;
	return STRING;
}

static void do_immediate_op(struct parse_term *v1,struct parse_term *v2,int op)
{
	int t1,t2;
	double x1,x2;
	string *s;
	
	t1=v1->type;
	if(v2) {
		t2=v2->type;
		/* Handle type promotion */
		if(op=='.') {
			if(t1!=STRING) t1=convert_string(v1);
			if(t2!=STRING) t2=convert_string(v2);
		} else if(t1!=t2 || (t1==STRING && !(op=='+' || op=='<' || op=='>' || op==EQ_TOKEN || op==NE_TOKEN || op==GEQ_TOKEN || op==LEQ_TOKEN))) {
			if(t1==INF || t1==XNAN) {
				if(t2==STRING) free_string(v2->elem.str);
				v2->type=t2=t1;
			} else if(t2==INF || t2==XNAN) {
				if(t1==STRING) free_string(v1->elem.str);
				v1->type=t1=t2;
			}
			if(t1!=INF && t1!=XNAN) {
				if(t1==STRING) t1=convert_numeric(v1);
				if(t2==STRING) t2=convert_numeric(v2);
				if(t1!=t2) {
					if(t1==INTEGER) t1=convert_real(v1);
					else t2=convert_real(v2);
				}
			}
		}
		switch(t1) {
		 case INF:
		 case XNAN:
			break;
		 case STRING:
			s=v1->elem.str;
			switch(op) {
			 case '+':
			 case '.':
				v1->elem.str=add_strings(s,copy_string(v2->elem.str)); 
				break;
			 case EQ_TOKEN:
				v1->elem.i=(!string_cmp(s,v2->elem.str))?1:0;
				v1->type=INTEGER;
				break;
			 case NE_TOKEN:
				v1->elem.i=string_cmp(s,v2->elem.str)?1:0;
				v1->type=INTEGER;
				break;
			 case '<':
				v1->elem.i=(string_cmp(s,v2->elem.str)<0)?1:0;
				v1->type=INTEGER;
				break;
			 case '>':
				v1->elem.i=(string_cmp(s,v2->elem.str)>0)?1:0;
				v1->type=INTEGER;
				break;
			 case LEQ_TOKEN:
				v1->elem.i=(string_cmp(s,v2->elem.str)<=0)?1:0;
				v1->type=INTEGER;
				break;
			 case GEQ_TOKEN:
				v1->elem.i=(string_cmp(s,v2->elem.str)>=0)?1:0;
				v1->type=INTEGER;
				break;
			}
			if(v1->type==INTEGER) free_string(s);
			break;
		 case INTEGER:
			switch(op) {
			 case '+':
				v1->elem.i+=v2->elem.i;
				break;
			 case '-':
				v1->elem.i-=v2->elem.i;
				break;
			 case '/':
				if(!v2->elem.i) {
					v1->elem.i=0;
					v1->type=INF;
				} else {
					x1=(double)v1->elem.i;
					x2=(double)v2->elem.i;
					x1/=x2;
					x2=modf(x1,&x2);
					if(fabs(x2)<WORKING_ZERO) v1->elem.i/=v2->elem.i;
					else {
						v1->elem.x=x1;
						v1->type=REAL;
					}
				}
				break;
			 case '*':
				v1->elem.i*=v2->elem.i;
				break;
			 case EXPO:
				if(v1->elem.i<0 || (!v1->elem.i && v2->elem.i)) {
					v1->type=XNAN;
				} else {
					x1=(double)v1->elem.i;
					x2=(double)v2->elem.i;
					x1=pow(x1,x2);
					if(isinf(x1)) v1->type=INF;
					else if(isnan(x1)) v1->type=XNAN;
					else {
						v1->type=REAL;
						v1->elem.x=x1;
					}
				}
				break;
			 case '<':
				v1->elem.i=v1->elem.i<v2->elem.i?1:0;
				break;
			 case '>':
				v1->elem.i=v1->elem.i>v2->elem.i?1:0;
				break;
			 case EQ_TOKEN:
				v1->elem.i=v1->elem.i==v2->elem.i?1:0;
				break;
			 case LEQ_TOKEN:
				v1->elem.i=v1->elem.i<=v2->elem.i?1:0;
				break;
			 case GEQ_TOKEN:
				v1->elem.i=v1->elem.i<=v2->elem.i?1:0;
				break;
			 case NE_TOKEN:
				v1->elem.i=v1->elem.i!=v2->elem.i?1:0;
				break;
			 case '|':
				v1->elem.i=v1->elem.i||v2->elem.i?1:0;
				break;
			 case '&':
				v1->elem.i=v1->elem.i&&v2->elem.i?1:0;
				break;
			}
			break;
		 case REAL:
			switch(op) {
			 case '+':
				v1->elem.x+=v2->elem.x;
				break;
			 case '-':
				v1->elem.x-=v2->elem.x;
				break;
			 case '/':
				if(fabs(v2->elem.x)<WORKING_ZERO) {
					v1->elem.x=0.0;
					v1->type=INF;
				} else v1->elem.x/=v2->elem.x;
				break;
			 case '*':
				v1->elem.x*=v2->elem.x;
				break;
			 case EXPO:
				x1=v1->elem.x;
				x2=v2->elem.x;
				if(x1<0.0 || (x1==0.0 && x2!=0.0)) {
					v1->type=XNAN;
				} else {
					x1=pow(x1,x2);
					if(isinf(x1)) v1->type=INF;
					else if(isnan(x1)) v1->type=XNAN;
					else v1->elem.x=x1;
				}
				break;
			 case '<':
				v1->elem.i=v1->elem.x<v2->elem.x?1:0;
				v1->type=INTEGER;
				break;
			 case '>':
				v1->elem.i=v1->elem.x>v2->elem.x?1:0;
				v1->type=INTEGER;
				break;
			 case EQ_TOKEN:
				v1->elem.i=fabs(v1->elem.x-v2->elem.x)<WORKING_ZERO?1:0;
				v1->type=INTEGER;
				break;
			 case LEQ_TOKEN:
				v1->elem.i=v1->elem.x<=v2->elem.x?1:0;
				v1->type=INTEGER;
				break;
			 case GEQ_TOKEN:
				v1->elem.i=v1->elem.x<=v2->elem.x?1:0;
				v1->type=INTEGER;
				break;
			 case NE_TOKEN:
				v1->elem.i=fabs(v1->elem.x-v2->elem.x)>WORKING_ZERO?1:0;
				v1->type=INTEGER;
				break;
			 case '|':
				v1->elem.i=v1->elem.x||v2->elem.x?1:0;
				v1->type=INTEGER;
				break;
			 case '&':
				v1->elem.i=v1->elem.x&&v2->elem.x?1:0;
				v1->type=INTEGER;
				break;
			}
		}
	} else {
		assert(op=='-' || op=='!');
		if(t1==STRING) t1=convert_numeric(v1);
		if(op=='-') {
			if(t1==INTEGER) v1->elem.i=-v1->elem.i;
			else v1->elem.x=-v1->elem.x;
		} else if(op=='!') {
			if(t1==INTEGER ) v1->elem.i=v1->elem.i?0:1;
			else {
				v1->elem.i=fabs(v1->elem.x)>WORKING_ZERO?0:1;
				v1->type=INTEGER;
			}
		}
	}
}

struct parse_term *do_op1(struct parse_term *arg1,struct parse_term *arg2,int op)
{
	int fg=0;
	
	if(arg2) {
		if((arg2=check_immediate(arg2))) fg=1;
	} else if((arg1=check_immediate(arg1))) fg=1;
	if(!fg) {
		parerror1("Got vector when expecting immediate argument\n");
		if(arg1) free_parse_term(arg1);
		if(arg2) free_parse_term(arg2);
	}
	return fg?do_op(arg1,arg2,op):0;
}

struct parse_term *do_op2(struct parse_term *arg1,struct parse_term *arg2,int op,int cl_flag)
{
	int tmp;
	struct parse_term *tm;
	
	tmp=handling_clause;
	handling_clause=cl_flag;
	tm=do_op(arg1,arg2,op);
	handling_clause=tmp;
	return tm;
}

struct parse_term *do_op(struct parse_term *arg1,struct parse_term *arg2,int op)
{
	int c1=0,c2=0,fg;
	struct expr_op *eop;
	struct deferred_assign *def;
	struct parse_term *tmp;

	if(!arg1) {
		if(arg2) return arg2;
		return 0;
	}
	if(!arg2 && op!='-' && op!='!') return arg1;
	if(op!=SUB_EXPR) {
		if(arg1->type==INTEGER || arg1->type==STRING || arg1->type==REAL || arg1->type==INF || arg1->type==XNAN) c1=1;
		if(!arg2 || arg2->type==INTEGER || arg2->type==STRING || arg2->type==REAL || arg2->type==INF || arg2->type==XNAN) c2=1;
	}
	if(op==',' && arg1->type!=EXPR_OP) {
		free_parse_term(arg1);
		arg1=arg2;
	} else if(op=='=') {
		if(c1) {
			parerror1("Can't assign to constant: '=' used instead of '==' perhaps?\n");
			free_parse_term(arg2);
		} else {
			if(!c2) arg2=try_resolve_expr(arg2);
			if(arg2) {
				assign_var(arg1,arg2);
				assign_call(arg1,arg2);
				if(arg2->defer) do_deferred(arg2);
				free_parse_term(arg2);
			}
		}
	} else {
		if(c1 && c2) { /* Are both immediate variables? */
			do_immediate_op(arg1,arg2,op);
			if(arg2) {
				if(arg2->defer) {
					def=arg2->defer;
					while(def->next) def=def->next;
					def->next=arg1->defer;
					arg1->defer=arg2->defer;
					arg2->defer=0;
				}
				free_parse_term(arg2);
			}
		} else {
			fg=0;
			if(op=='/') {
				tmp=copy_parse_term(arg2);
				if(!c2) tmp=try_resolve_expr(tmp);
				if(tmp) {
					if((tmp->type&VT_TYPES)==STRING) (void)convert_numeric(tmp);
					switch(tmp->type&VT_TYPES) {
					 case INTEGER:
						if(!tmp->elem.i) fg=1;
						else if(tmp->elem.i==1) fg=2;
						break;
					 case REAL:
						if(fabs(tmp->elem.x)<WORKING_ZERO) fg=1;
						else if(fabs(tmp->elem.x-1.0)<WORKING_ZERO) fg=2;
						break;
					}
					if(fg) {
						if(fg==1) arg1->type=INF|VT_VECTOR;
						free_parse_term(arg2);
					}
					free_parse_term(tmp);
				}
			}
			if(!fg) {
				/* Make new op */
				eop=lk_malloc(sizeof(struct expr_op));
				eop->exp1=arg1;
				eop->exp2=arg2;
				eop->op=op;
				arg1=get_new_term();
				arg1->type=EXPR_OP;
				arg1->elem.eop=eop;
			}
		}
	}	
	return arg1;
}

void set_flag(struct parse_term *ex,int flag)
{
	int type;
	struct expr_op *eop;
	struct parse_var *vv;
	
	if(!ex) return;
	type=ex->type&VT_TYPES;
	if(type==VARIABLE) {
		ex->type|=flag;
		vv=find_parse_var(ex);
		vv->type|=flag;
	} else if(type==EXPR_OP) {
		eop=ex->elem.eop;
		if(eop->exp1) set_flag(eop->exp1,flag);
		if(eop->exp2) set_flag(eop->exp2,flag);
	}
}

void set_var_list_flag(struct parse_var_list *vl,int flag)
{
	while(vl) {
		set_flag(&vl->term,flag);
		vl=vl->next;
	}
}

void set_vector_flag(struct parse_term *ex)
{
	set_flag(ex,VT_VECTOR);
}

void set_var_list_vector_flag(struct parse_var_list *vl)
{
	set_var_list_flag(vl,VT_VECTOR);
}

void free_parse_var(struct parse_var *vv)
{
	if(!vv) return;
	switch(vv->type&VT_TYPES) {
	 case STRING:
		free_string(vv->arg.str);
		break;
	 case VARIABLE:
		if((vv->type&VT_ARRAY)&&vv->size) free(vv->arg.vv);
		break;
	}
	free_string(vv->name);
	free(vv);
}

struct parse_term *copy_parse_term(struct parse_term *t)
{
	struct parse_term *t1;
	struct parse_var_elem *elem;
	struct expr_op *eop;
	struct func_op *fop;
	t1=get_new_term();
	t1->type=t->type;
	t1->elem=t->elem;
	t1->defer=t->defer;
	t1->clause_assign=t->clause_assign;
	switch(t->type&VT_TYPES) {
	 case VARIABLE:
		if(t->type&VT_ARRAY_ELEMENT) {
			elem=lk_malloc(sizeof(struct parse_var_elem));
			elem->head=t->elem.vv1->head;
			elem->type=t->elem.vv1->type;
			elem->ex=copy_parse_term(t->elem.vv1->ex);
			t1->elem.vv1=elem;
		}
		break;
	 case EXPR_OP:
		eop=lk_malloc(sizeof(struct expr_op));
		eop->exp1=copy_parse_term(t->elem.eop->exp1);
		eop->exp2=copy_parse_term(t->elem.eop->exp2);
		eop->op=t->elem.eop->op;
		t1->elem.eop=eop;
		break;
	 case FUNCTION:
		fop=lk_malloc(sizeof(struct func_op));
		fop->ex=copy_parse_term(t->elem.fop->ex);
		fop->f=t->elem.fop->f;
		t->elem.fop=fop;
		break;
	 case STRING:
		t1->elem.str=copy_string(t->elem.str);
		break;
	}
	return t1;
}

void free_defer(struct deferred_assign *d)
{
	struct deferred_assign *d1;
	
	while(d) {
		d1=d->next;
		free_parse_term(d->var);
		free_parse_term(d->expr);
		free(d);
		d=d1;
	}
}

void free_parse_term(struct parse_term *vv)
{
	struct expr_op *eop;
	
	if(vv) {
		switch(vv->type&VT_TYPES) {
		 case VARIABLE:
			if(vv->type&VT_ARRAY_ELEMENT) {
				free_parse_term(vv->elem.vv1->ex);
				free(vv->elem.vv1);
			}
			break;
		 case EXPR_OP:
			eop=vv->elem.eop;
			free_parse_term(eop->exp1);
			free_parse_term(eop->exp2);
			free(eop);
			break;
		 case FUNCTION:
			free_parse_term(vv->elem.fop->ex);
			free(vv->elem.fop);
			break;
		 case STRING:
			free_string(vv->elem.str);
			break;
		}
		if(vv->clause_assign) free_defer(vv->clause_assign);
		free(vv);
	}
	return;
}

void free_parse_var_list(struct parse_var_list *vl)
{
	struct parse_var_list *vl1;

	while(vl) {
		vl1=vl->next;
		if(vl->clause) free_parse_clause(vl->clause);
		if(vl->term.type&VT_ARRAY_ELEMENT) {
			free_parse_term(vl->term.elem.vv1->ex);
			free(vl->term.elem.vv1);
		}
		if((vl->term.type&VT_TYPES)==STRING) free_string(vl->term.elem.str);
		free(vl);
		vl=vl1;
	}
}

/* Split a term into sub-expressions */
struct parse_term *get_term_list(struct parse_term *term,int *i)
{
	int fg=0;
	static struct parse_term *list;
	static int list_size;
	struct expr_op *eop;
	
	if(term) {
		if(!list) {
			list_size=4;
			if(list_size<=*i) list_size=(*i)*1.5;
			list=lk_malloc(sizeof(struct parse_term)*list_size);
		} else if(list_size<=*i) {
			list_size=(*i)*1.5;
			list=lk_realloc(list,sizeof(struct parse_term)*list_size);
		}
		if((term->type&VT_TYPES)==EXPR_OP) {
			eop=term->elem.eop;
			if(eop->op==SUB_EXPR) {
				assert(eop->exp1 && eop->exp2);
				get_term_list(eop->exp1,i);
				get_term_list(eop->exp2,i);
				fg=1;
			}
		}
		if(!fg) {
			list[*i].elem=term->elem;
			list[*i].defer=term->defer;
			list[*i].clause_assign=term->clause_assign;
			list[(*i)++].type=term->type;
		}
	} else {
		if(list) {
			free(list);
			list=0;
		}
	}
	return list;
}

void assign_var(struct parse_term *term,struct parse_term *ex)
{
	int type;
	struct parse_var *vv,*base_vv=0;
	
	if(handling_clause) return;
	if(term->type&VT_ARRAY_ELEMENT) {
		base_vv=term->elem.vv1->head;
		vv=get_array_var(term->elem.vv1->head,&term->elem.vv1->ex);
	} else {
		vv=term->elem.vv;
		if(vv->type&VT_ARRAY) {
			parerror2("Array variable '%s' used in scalar context\n",var_name(term));
			return;
		}
		vv->type|=VT_SCALAR;
	}
	if(!vv) return;
	if((vv->type&VT_TYPES)==STRING) free_string(vv->arg.str);
	type=(ex->type&VT_TYPES);
	if(type==EXPR_OP || type==VARIABLE) {
		/* Check for type clash */
		if((vv->type&VT_IMMED) || (base_vv && (base_vv->type&VT_IMMED))) {
			parerror2("Variable '%s' used as both immediate and vector\n",var_name(term));
		}
		vv->type&=~(VT_TYPES|VT_IMMED);
		vv->type|=VARIABLE|VT_VECTOR;
		if(base_vv) base_vv->type|=VT_VECTOR;
		if(type==VARIABLE) ex->type|=VT_VECTOR;
		if(type==EXPR_OP) set_vector_flag(ex);
	} else {
		if((vv->type&VT_VECTOR) || (base_vv && (base_vv->type&VT_VECTOR))) {
			parerror2("Variable '%s' used as both immediate and vector\n",var_name(term));
		}
		vv->type&=~(VT_TYPES|VT_VECTOR);
		vv->type|=(ex->type|VT_IMMED);
		vv->arg=ex->elem;
		if(base_vv) base_vv->type|=VT_IMMED;
		if(type==STRING) vv->arg.str=copy_string(ex->elem.str);
	}
}

string *check_string(struct parse_term *vv)
{
	string *s=0;
	int type;
	
	if(vv) vv=try_resolve_expr(vv);
	if(vv) {
		type=vv->type&VT_TYPES;
		if(type==STRING) s=copy_string(vv->elem.str); 
		else if(type==VARIABLE) parerror1("Unknown variable '%s'\n",var_name(vv));
		else parerror1("Can't resolve string expression\n");
		free_parse_term(vv);
	}
	return s;
}

int convert_to_int(struct parse_term *ex,int *err)
{
	char *p;
	double x=0.0,x1,fr;
	int i=0,type;
	struct parse_var *v;
	expr_elem *ee;
	
	*err=0;
	if(!ex) *err=-1;
	else {
		type=ex->type&VT_TYPES;
		if(type==VARIABLE) {
			v=ex->elem.vv;
			ee=&v->arg;
			type=v->type&VT_TYPES;
		} else ee=&ex->elem;
		if(type==INTEGER) i=ee->i;
		else {
			/* Try and convert into integer form */
			if(type==EXPR_OP || type==VARIABLE) *err=1;
			else if(type==STRING) {
				x=strtod(get_cstring(ee->str),&p);
				if(*p) *err=2;
			} else x=ee->x;
			fr=modf(x,&x1);
			if(fabs(fr)<WORKING_ZERO) i=(int)x1;
			else *err=3;
		}
	}
	return i;
}

static int check_array_size(struct parse_var *vv,int i)
{
	int er=0,j;
	
	if(handling_clause) er=1;
	else {
		if(vv->type&VT_SCALAR) {
			parerror1("Scalar variable '%s' used in array context\n",get_cstring(vv->name));
			er=1;
		} else {
			if(vv->size && (vv->type&VT_ARRAY)) {
				if(i>vv->size) {
					if(vv->type&VT_SIZE_SET) {
						parerror1("The size of array '%s' has been set to %d, and can not be increased\n",get_cstring(vv->name),vv->size);
						er=1;
					} else {
						vv->arg.vv=lk_realloc(vv->arg.vv,sizeof(struct parse_var)*i);
						memset(vv->arg.vv+vv->size,0,sizeof(struct parse_var)*(i-vv->size));
						for(j=vv->size;j<i;j++) vv->arg.vv[i].name=vv->name;
						vv->size=i;
					}
				}
			} else {
				vv->arg.vv=lk_calloc((size_t)i,sizeof(struct parse_var));
				vv->size=i;
				vv->type|=VT_ARRAY;
				for(j=0;j<i;j++) vv->arg.vv[i].name=vv->name;
			}
		}
	}
	return er;
}

struct parse_var *get_array_var(struct parse_var *vv,struct parse_term **ex)
{
	int i,er;
	struct parse_var *elem=0;
	
	if(!(vv&& *ex)) return 0;
	*ex=try_resolve_expr(*ex);
	if(!*ex) return 0;
	i=convert_to_int(*ex,&er);
	if(er) parerror1("Array index for '%s' is not an integer expression\n",get_cstring(vv->name));
	if(!er && (i<1 || i>MAX_ARRAY_IDX)) {
		parerror1("Invalid array index (%d) for '%s'\n",i,get_cstring(vv->name));
		er=1;
	}
	if(!er) er=check_array_size(vv,i);
	if(!er) {
		elem=vv->arg.vv+i-1;
		if(!elem->type) {
			elem->type=VARIABLE|VT_ARRAY_ELEMENT;
			elem->arg.vv=0;
		}
	}
	return er?0:elem;
}

struct parse_term *get_array_term(struct parse_var *vv,struct parse_term *ex)
{
	struct parse_term *term=0;
	struct parse_var_elem *vv1;
	
	if(!(vv&&ex)) return 0;
	vv1=lk_malloc(sizeof(struct parse_var_elem));
	vv1->head=vv;
	vv1->ex=ex;
	term=get_new_term();
	term->type=VARIABLE|VT_ARRAY_ELEMENT;
	term->elem.vv1=vv1;
	return term;
}

static struct parse_var_list *get_new_varlist(void)
{
	struct parse_var_list *vl;
	
	vl=lk_malloc(sizeof(struct parse_var_list));
	vl->next=0;
	vl->term.elem.vv=0;
	vl->term.type=0;
	vl->term.defer=vl->term.clause_assign=0;
	vl->clause=0;
	vl->pos.fname=0;
	vl->pos.line=0;
	vl->pos.col=0;
	return vl;
}

struct parse_var_list *get_array_var_list(struct parse_var *vv,struct parse_term *ex)
{
	struct parse_var_list *vl=0;
	struct parse_var_elem *em;
	
	if(ex && vv) {
		if(vv->type&VT_SCALAR) {
			parerror1("Scalar variable '%s' used in array context\n",get_cstring(vv->name));
		} else {
			vv->type|=VT_VECTOR|VT_ARRAY;
			em=lk_malloc(sizeof(struct parse_var_elem));
			em->head=vv;
			em->ex=ex;
			vl=get_new_varlist();
			vl->term.type=VARIABLE|VT_ARRAY_ELEMENT;
			vl->term.elem.vv1=em;
		}
	}
	return vl;
}

struct parse_var_list *make_var_list(struct parse_term *v,struct parse_clause *clause)
{
	struct parse_var_list *vl=0;
	
	if(v) {
		v=try_resolve_expr(v);
		if(v) {
			if((v->type&VT_TYPES)==VARIABLE) v->type|=VT_VECTOR;
			vl=get_new_varlist();
			vl->clause=clause;
			vl->term.elem=v->elem;
			vl->term.type=v->type;
			free(v);
		}
	}
	if(!v) {
		vl=get_new_varlist();
		vl->clause=clause;
		vl->term.type=EMPTY_VAR;
	}
	if(clause) handle_clause(clause);
	return vl;
}

static void free_format_clause(struct format_clause *fc)
{
	struct format_clause *fc1;
	
	while(fc) {
		fc1=fc->next;
		if(fc->type==FC_SUBCLAUSE) free_format_clause(fc->elem);
		free(fc);
		fc=fc1;
	}
}

void free_parse_clause(struct parse_clause *pc)
{
	struct parse_clause *pc1;
	
	while(pc) {
		pc1=pc->next;
		switch(pc->type) {
		 case PC_TERM:
			free_parse_term(pc->clause.term);
			break;
		 case PC_FORMAT_CLAUSE:
			free_format_clause(pc->clause.format);
			break;
		}
		free(pc);
		pc=pc1;
	}
}

struct parse_var_list *addto_var_list(struct parse_var_list *vl1,struct parse_var_list *vl2) 
{
	struct parse_var_list *vl;
	
	if(vl2) {
		vl=vl2;
		while(vl->next) vl=vl->next;
		vl->next=vl1;
	} else vl2=vl1;
	return vl2;
}

struct parse_var_list *make_loop_clause(struct parse_var_list *vl,struct parse_var *v,struct parse_term *ex)
{
	struct parse_var_list *vl1=0,*vl2,*vl3,*vl_bk;
	struct parse_term *v1,*term,*expr;
	struct parse_var_elem *vv1=0;
	int i,j,ix[3],er=0;
	
	if(vl) {
		/* Make up list to interpolate into existing list */
		if(v) {
			if(v->type&VT_VECTOR) {
				parerror1("Illegal use of vector variable '%s' as loop index\n",get_cstring(v->name));
				er=1;
			} else if(v->type&VT_ARRAY) {
				parerror1("Illegal use of array variable '%s' in scalar context\n",get_cstring(v->name));
				er=1;
			}
		}
		if(!er) {
			if(v && ex) {
				i=0;
				expr=get_term_list(ex,&i);
				if(i<2) {
					parerror1("Too few expressions found in implied loop definition\n");
					er=1;
				} else if(i>3) {
					parerror1("Too many expressions found in implied loop definition\n");
					er=1;
				}
				if(!er) {
					for(j=0;j<i;j++) {
						ix[j]=convert_to_int(expr+j,&er);
						if(er) {
							parerror1("Non-integer expression in implied loop definition\n");
							break;
						}
					}
				}
				if(!er) {
					if(i<3) ix[2]=1; /* Set default step-size */
					else if(!ix[2]) parerror1("Zero step-size in implied loop definition\n");
					else {
						if(ix[2]<0) {
							i=ix[0];
							ix[0]=ix[1];
							ix[1]=i;
							ix[2]=-ix[2];
						}
						if(ix[0]<1 || ix[0]>MAX_ARRAY_IDX || ix[1]<1 || ix[1]>MAX_ARRAY_IDX) {
							parerror1("Illegal index values in implied loop definition\n");
							ix[2]=0;
						}
					}
				}
			} else ix[0]=ix[1]=ix[2]=1;
			if(ix[2]>0 && ix[1]>=ix[0]) {
				if(v) {
					v->type&=~VT_TYPES;
					v->type|=INTEGER;
				}
				vl_bk=vl;
				for(i=ix[0];i<=ix[1];i+=ix[2]) {
					if(v) v->arg.i=i;
					vl=vl_bk;
					vl3=0;
					do {
						v1=&vl->term;
						if(v1->type&VT_ARRAY_ELEMENT) {
							vv1=lk_malloc(sizeof(struct parse_var_elem));
							vv1->head=v1->elem.vv1->head;
							vv1->ex=v1->elem.vv1->ex;
							term=copy_parse_term(vv1->ex);
							/* Try and evaluate index */
							term=try_resolve_expr(term);
							vv1->ex=term;
							if(!term) v1=0;
						}
						if(v1) {
							vl2=get_new_varlist();
							vl2->next=vl3;
							vl3=vl2;
							if(v1->type&VT_ARRAY_ELEMENT) vl2->term.elem.vv1=vv1;
							else vl2->term.elem=v1->elem;
							vl2->term.type=v1->type;
							vl2->clause=vl->clause;
						}
						vl=vl->next;
					} while(vl);
					while(vl3) {
						vl2=vl3->next;
						vl3->next=vl1;
						vl1=vl3;
						vl3=vl2;
					}
				}
				vl=vl_bk;
				do {
					v1=&vl->term;
					if(v1->type&VT_ARRAY_ELEMENT) check_array_size(v1->elem.vv1->head,ix[1]);
					vl=vl->next;
				} while(vl);
				vl=vl_bk;
			}
		}
		free_parse_var_list(vl);
		if(ex) free_parse_term(ex);
	}
	return vl1;
}

struct parse_term *get_var_term(struct parse_var *v)
{
	struct parse_term *term=0;
	
	if(v) {
		term=get_new_term();
		term->type=VARIABLE;
		term->elem.vv=v;
	}
	return term;
}

struct parse_term *try_resolve_expr(struct parse_term *term)
{
	struct parse_term *t1,*t2;
	struct expr_op *eop;
	struct func_op *fop;
	struct parse_var *vv;
	int type;

	if(term) {
		type=term->type&VT_TYPES;
		if(type==VARIABLE) {
			if(term->type&VT_ARRAY_ELEMENT) vv=get_array_var(term->elem.vv1->head,&term->elem.vv1->ex);
			else vv=term->elem.vv;
			if(vv) {
				type=vv->type&VT_TYPES;
				if(type==INTEGER || type==STRING || type==REAL || type==INF || type==XNAN) {
					t1=get_new_term();
					t1->type=type;
					t1->defer=term->defer;
					t1->clause_assign=term->clause_assign;
					if(type==STRING) t1->elem.str=copy_string(vv->arg.str);
					else t1->elem=vv->arg;
					term->defer=term->clause_assign=0;
					free_parse_term(term);
					term=t1;
				}
			} else {
				free_parse_term(term);
				term=0;
			}
		} else if(type==EXPR_OP) {
			eop=term->elem.eop;
			t1=try_resolve_expr(eop->exp1);
			if(eop->exp2) t2=try_resolve_expr(eop->exp2);
			else t2=0;
			t1=do_op(t1,t2,eop->op);
			free(eop);
			free(term);
			term=t1;
		} else if(type==FUNCTION) {
			fop=term->elem.fop;
			t1=try_resolve_expr(fop->ex);
			t1=do_func(fop->f,t1);
			free(fop);
			free(term);
			term=t1;
		}
	}
	return term;
}
