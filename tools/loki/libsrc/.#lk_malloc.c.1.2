/****************************************************************************
 *                                                                          *
 *     Loki - Programs for genetic analysis of complex traits using MCMC    *
 *                                                                          *
 *                   Simon Heath - CNG, Paris                               *
 *                                                                          *
 *                        March 2004                                        *
 *                                                                          *
 * lk_malloc.c:                                                             *
 *                                                                          *
 * Wrappers for malloc calls with configurable default responses            *
 *                                                                          *
 * Copyright (C) Simon C. Heath 2004                                        *
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
#include <assert.h>

#include "config.h"
#include "utils.h"
#include "lk_malloc.h"

void *_lk_malloc(size_t size,const char *func,char *file,int line)
{
	void *p;
	
#ifndef NDEBUG
	if(!size) abt(file,line,"%s(): lk_malloc called with zero size argument\n",func);
#endif
#ifdef __DMALLOC_H__
	p=dmalloc_malloc(file,line,size,DMALLOC_FUNC_MALLOC,0,0);
#else
	p=malloc(size);
#endif
	if(!p) abt(file,line,"%s(): out of memory - %d %s requested\n",func,size,size==1?"byte":"bytes");
	return p;
}

void *_lk_realloc(void *p,size_t size,const char *func,char *file,int line)
{
	void *p1;
	
#ifndef NDEBUG
	if(!size) abt(file,line,"%s(): lk_realloc called with zero size argument\n",func);
	if(!p) abt(file,line,"%s(): lk_realloc called with zero pointer argument\n",func);
#endif
#ifdef __DMALLOC_H__
	p1=dmalloc_realloc(file,line,p,size,DMALLOC_FUNC_REALLOC,0);
#else
	p1=realloc(p,size);
#endif
	if(!p1) abt(file,line,"%s(): out of memory - %d %s requested\n",func,size,size==1?"byte":"bytes");
	return p1;
}

void *_lk_calloc(size_t number, size_t size,const char *func,char *file,int line)
{
	void *p;
	
#ifndef NDEBUG
	if(!(size && number)) abt(file,line,"%s(): lk_calloc called with arguments %d,%d\n",func,number,size);
#endif
#ifdef __DMALLOC_H__
	p=dmalloc_malloc(file,line,number*size,DMALLOC_FUNC_CALLOC,0,0);
#else
	p=calloc(number,size);
#endif
	if(!p) abt(file,line,"%s(): out of memory - %d x %d %s requested\n",func,number,size,number*size==1?"byte":"bytes");
	return p;
}
