/****************************************************************************
 *                                                                          *
 *     Loki - Programs for genetic analysis of complex traits using MCMC    *
 *                                                                          *
 *                     Simon Heath - CNG, Evry                              *
 *                                                                          *
 *                         January 2003                                     *
 *                                                                          *
 * compat.c:                                                                *
 *                                                                          *
 * Compatability utilities                                                  *
 *                                                                          *
 * Copyright (C) Simon C. Heath 2003                                        *
 * This is free software.  You can distribute it and/or modify it           *
 * under the terms of the Modified BSD license, see the file COPYING        *
 *                                                                          *
 ****************************************************************************/

#include <config.h>
#include <stdlib.h>
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif
#include <string.h>
#include <stdio.h>
#include <ctype.h>

#include "utils.h"
#include "compat.h"
#include "scan.h"

int check_format(char *name)
{
	char *format_names[]={"LOKI","QTDT","MERLIN","SOLAR","LINKAGE",0};
	int format_types[]={LOKI_FORMAT,QTDT_FORMAT,QTDT_FORMAT,SOLAR_FORMAT,LINKAGE_FORMAT};
	
	char *p,*p1;
	int i=0;
	
	p1=name+2;
	while(*p1 && isspace((int)*p1)) p1++;
	p=format_names[0];
	while(p) {
		if(!strcasecmp(p,p1)) break;
		p=format_names[++i];
	}
	if(!p) {
		fprintf(stderr,"Unknown format type (%s)\n",p1);
		exit(EXIT_FAILURE);
	}
	return format_types[i];
}

