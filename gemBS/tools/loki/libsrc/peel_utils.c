/****************************************************************************
 *                                                                          *
 *     Loki - Programs for genetic analysis of complex traits using MCMC    *
 *                                                                          *
 *             Simon Heath - University of Washington                       *
 *                                                                          *
 *                        August 1997                                       *
 *                                                                          *
 * peel_util.c:                                                             *
 *                                                                          *
 * Copyright (C) Simon C. Heath 1997, 2000, 2002                            *
 * This is free software.  You can distribute it and/or modify it           *
 * under the terms of the Modified BSD license, see the file COPYING        *
 *                                                                          *
 ****************************************************************************/

#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "shared_peel.h"

int num_bits(const int n_all)
{
	int nb;	

	nb=(int)(0.99999+log((double)n_all)/log(2.0));
	return nb;
}

void free_peelseq(struct Peelseq_Head *pp)
{
	struct Simple_Element *simple_em;
	struct Fenris_Simple_Element *fsimple_em;
	struct Complex_Element *complex_em;
	void *old_ptr=0;
	
	while(pp->type) {
		switch(pp->type) {
		 case PEEL_SIMPLE:
			simple_em=pp->ptr.simple;
			if(old_ptr) free(old_ptr);
			pp= &simple_em->next;
			if(simple_em->off) free(simple_em->off);
			old_ptr=(void *)simple_em;
			break;
		 case FENRIS_PEEL_SIMPLE:
			fsimple_em=pp->ptr.fsimple;
			if(old_ptr) free(old_ptr);
			pp= &fsimple_em->next;
			if(fsimple_em->off) free(fsimple_em->off);
			old_ptr=(void *)fsimple_em;
			break;
		 case PEEL_COMPLEX:
			complex_em=pp->ptr.complex;
			if(old_ptr) free(old_ptr);
			pp= &complex_em->next;
			free(complex_em->involved);
			old_ptr=(void *)complex_em;
			break;
		}
	}
	if(old_ptr) free(old_ptr);
}
