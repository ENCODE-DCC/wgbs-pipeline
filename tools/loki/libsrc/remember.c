/****************************************************************************
 *                                                                          *
 *     Loki - Programs for genetic analysis of complex traits using MCMC    *
 *                                                                          *
 *             Simon Heath - University of Washington                       *
 *                                                                          *
 *                       March 1997                                         *
 *                                                                          *
 * remember.c:                                                              *
 *                                                                          *
 * Routines for keeping track of allocated memory that is used in weird     *
 * ways so that it is not easy to keep track of it.                         *
 *                                                                          *
 * Copyright (C) Simon C. Heath 1997, 2000, 2002                            *
 * This is free software.  You can distribute it and/or modify it           *
 * under the terms of the Modified BSD license, see the file COPYING        *
 *                                                                          *
 ****************************************************************************/

#include <config.h>
#include <stdlib.h>
#ifdef USE_DMALLOC
#include <dmalloc.h>
#endif
#include <stdio.h>
#include <assert.h>

#include "utils.h"
#include "lk_malloc.h"
#include "loki_struct.h"

/* Adds a memory pointer (from malloc()) to list.  Used by FreeStuff()
 * So all can be free()'d with one call */
struct remember *AddRemem(void *p,struct remember *rblock)
{
	struct remember *pr;
	
	assert(p && rblock);
	if(rblock->pos==REMSIZE) {
		/* If allocation of new block fails, we'll just return old block pointer so new
		 * location is not remembered */
		if(!(pr=lk_malloc(sizeof(struct remember)))) return rblock;
		rblock->next=pr;
		rblock=pr;
		rblock->pos=0;
		rblock->next=0;
	}
	rblock->mem[rblock->pos++]=p;
	return rblock;
}

void FreeRemem(struct remember *rblock)
{
	int i;
	struct remember *pr;
	
	while(rblock) {
		pr=rblock->next;
		for(i=0;i<rblock->pos;i++)	{
			free(rblock->mem[i]);
		}
		free(rblock);
		rblock=pr;
	}
}
