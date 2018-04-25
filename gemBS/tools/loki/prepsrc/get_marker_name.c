/****************************************************************************
 *                                                                          *
 *     Loki - Programs for genetic analysis of complex traits using MCMC    *
 *                                                                          *
 *             Simon Heath - CNG, Paris                                     *
 *                                                                          *
 *                        August 2002                                       *
 *                                                                          *
 * get_marker_name.c:                                                       *
 *                                                                          *
 * Copyright (C) Simon C. Heath 1997, 2000, 2002                            *
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

#include "utils.h"
#include "scan.h"
#include "snprintf.h"

#define TBUF_SIZE 32
char *get_marker_name(const int locus)
{
	char *s,tbuf[TBUF_SIZE];
	size_t i;
	struct Marker *mark;
	if(locus==n_markers) mark=traitlocus;
	else mark=markers+locus;
	
	if(mark->var->vtype&ST_ARRAY)	{
		(void)snprintf(tbuf,TBUF_SIZE,"%d",mark->index);
		i=strlen(tbuf)+strlen(mark->var->name)+3;
		if(!(s=malloc(i))) ABT_FUNC(MMsg);
		(void)snprintf(s,i,"%s(%s)",mark->var->name,tbuf);
	} else {
		i=strlen(mark->var->name)+1;
		if(!(s=malloc(i))) ABT_FUNC(MMsg);
		(void)strncpy(s,mark->var->name,i);
	}
	return s;
}

