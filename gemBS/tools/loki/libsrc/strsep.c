/****************************************************************************
 *                                                                          *
 *     Loki - Programs for genetic analysis of complex traits using MCMC    *
 *                                                                          *
 *                   Simon Heath - CNG, Evry                                *
 *                                                                          *
 *                        March 2003                                        *
 *                                                                          *
 * strsep.c:                                                                *
 *                                                                          *
 * Find first occurrence of a token in a string                             *
 *                                                                          *
 * Copyright (C) Simon C. Heath 2003                                        *
 * This is free software.  You can distribute it and/or modify it           *
 * under the terms of the Modified BSD license, see the file COPYING        *
 *                                                                          *
 ****************************************************************************/

#include <config.h>
#include <stdlib.h>
#include <stdio.h>
#include "libhdr.h"

char *my_strsep(char **s,const char *d)
{
	char *p,*s1,c,c1;
	const char *dp;
	
	p=*s;
	if(p) {
		s1=p;
		for(;;) {
			c=*p++;
			dp=d;
			do {
				c1=*dp++;
				if(c1==c) {
					if(!c) p=0;
					else p[-1]=0;
					*s=p;
					return s1;
				}
			} while(c1);
		}
	}
	return p;
}
