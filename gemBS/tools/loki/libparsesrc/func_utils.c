/****************************************************************************
 *                                                                          *
 *     Loki - Programs for genetic analysis of complex traits using MCMC    *
 *                                                                          *
 *                     Simon Heath - CNG, Evry                              *
 *                                                                          *
 *                           July 2003                                      *
 *                                                                          *
 * func_utils.c:                                                            *
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
#ifdef HAVE_IEEEFP_H
#include <ieeefp.h>
#endif
#include <assert.h>

#include "config.h"
#include "parser.h"
#include "parse.tab.h"
#include "lk_malloc.h"

static func_def *add_new_func(func_def **pp)
{
	func_def *p;
	
	if(*pp) {
		p=*pp;
		while(p->next) p=p->next;
		p->next=lk_malloc(sizeof(func_def));
		p=p->next;
	} else {
		p=lk_malloc(sizeof(func_def));
		*pp=p;
	}
	p->next=0;
	return p;
}

func_def *register_string_function(func_def *p,string *(*f)(string *),char *com)
{
	func_def *p1;

	p1=add_new_func(&p);
	p1->com=com;
	p1->type=STRING;
	p1->f.f_string=f;
	return p;
}

func_def *register_int_function(func_def *p,int (*f)(int),char *com)
{
	func_def *p1;

	p1=add_new_func(&p);
	p1->com=com;
	p1->type=INTEGER;
	p1->f.f_int=f;
	return p;
}

func_def *register_real_function(func_def *p,double (*f)(double),char *com)
{
	func_def *p1;

	p1=add_new_func(&p);
	p1->com=com;
	p1->type=REAL;
	p1->f.f_real=f;
	return p;
}

static void do_immediate_func(func_def *f,struct parse_term *ex)
{
	int type,type1;
	double z;
	
	type=ex->type;
	if(type!=INF && type!=XNAN) {
		type1=f->type;
		if(type!=type1) {
			switch(type1) {
			 case INTEGER:
				if(type==STRING) type=convert_numeric(ex);
				if(type==REAL)	{
					ex->elem.i=(int)ex->elem.x;
					ex->type=INTEGER;
				}
				ex->elem.i=f->f.f_int(ex->elem.i);
				break;
			 case REAL:	
				if(type==STRING) type=convert_numeric(ex);
				if(type==INTEGER) {
					ex->elem.x=(double)ex->elem.i;
					ex->type=REAL;
				}
				z=f->f.f_real(ex->elem.x);
				if(isnan(z)) ex->type=XNAN;
				else if(isinf(z)) ex->type=INF;
				else ex->elem.x=z;
				break;
			 case STRING:
				ex->type=STRING;
				if(type==REAL) ex->elem.str=real_to_string(ex->elem.x);
				else if(type==INTEGER) ex->elem.str=int_to_string(ex->elem.i);
				ex->elem.str=f->f.f_string(ex->elem.str);
				break;
			}
		} else {
			switch(type1) {
			 case INTEGER:
				ex->elem.i=f->f.f_int(ex->elem.i);
				break;
 			 case REAL:
				z=f->f.f_real(ex->elem.x);
				if(isnan(z)) ex->type=XNAN;
				else if(isinf(z)) ex->type=INF;
				else ex->elem.x=z;
				break;
 			 case STRING:
				ex->elem.str=f->f.f_string(ex->elem.str);
				break;
			}
		}
	}
}

struct parse_term *do_func(func_def *f,struct parse_term *ex)
{
	int type;
	struct func_op *fop;

	type=ex->type;
	if(type==INTEGER || type==STRING || type==REAL || type==INF || type==XNAN) {
		do_immediate_func(f,ex);
	} else { /* Make new op */
		fop=lk_malloc(sizeof(struct func_op));
		fop->ex=ex;
		fop->f=f;
		ex=get_new_term();
		ex->type=FUNCTION;
		ex->elem.fop=fop;
	}
	return ex;
}

void free_functions(func_def *f)
{
	func_def *f1;
	
	while(f) {
		f1=f->next;
		free(f);
		f=f1;
	}
}
