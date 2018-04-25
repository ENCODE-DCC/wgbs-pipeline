/****************************************************************************
 *                                                                          *
 *     Loki - Programs for genetic analysis of complex traits using MCMC    *
 *                                                                          *
 *                    Simon Heath - CNG, France                             *
 *                                                                          *
 *                         August 2002                                      *
 *                                                                          *
 * fenris_utils.c:                                                          *
 *                                                                          *
 * Copyright (C) Simon C. Heath 2002                                        *
 * This is free software.  You can distribute it and/or modify it           *
 * under the terms of the Modified BSD license, see the file COPYING        *
 *                                                                          *
 ****************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#ifdef USE_DMALLOC
#include <dmalloc.h>
#endif

#include "utils.h"
#include "loki.h"
#include "loki_peel.h"
#include "min_deg.h"
#include "fenris_peel.h"
#include "fenris.h"

void print_peelseq_element(FILE *fptr,struct Peelseq_Head *p)
{
	int i,j,n_peel,n_inv,n_off,*inv,*flag;
	struct Fenris_Simple_Element *simple_em;
   struct Complex_Element *complex_em; 
	
	if(p->type==FENRIS_PEEL_SIMPLE) {
		simple_em=p->ptr.fsimple;
		print_orig_id(fptr,simple_em->sire);
		fputc(',',fptr);
		print_orig_id(fptr,simple_em->dam);
		fputs("->",fptr);
		if(simple_em->pivot) {
			print_orig_id(fptr,simple_em->pivot);
			fprintf(fptr," [%d]",simple_em->out_index);
		} else fputc('*',fptr);
		fputc(' ',fptr);
		n_off=simple_em->n_off;
		for(i=0;i<n_off;i++) {
			fputc(i?',':'(',fptr);
			print_orig_id(fptr,simple_em->off[i]+1);
		}
		fputs(") ",fptr);
		for(j=i=0;i<2+n_off;i++) if(simple_em->rf[i]>=0) {
			fputc(j++?',':'[',fptr);
			fprintf(fptr,"%d",simple_em->rf[i]);
		}
		if(j) fputc(']',fptr);
		fputc('\n',fptr);
	} else if(p->type==PEEL_COMPLEX) {
		complex_em=p->ptr.complex;
		n_inv=complex_em->n_involved;
		n_peel=complex_em->n_peel;
		inv=complex_em->involved;
		flag=complex_em->flags;
		for(i=0;i<n_peel;i++) {
			if(i) fputc(',',fptr);
			if(flag[i]&1) fputc('<',fptr);
			print_orig_id(fptr,inv[i]+1);
			if(flag[i]&1) fputc('>',fptr);
		}
		fputs("->",fptr);
		if(i<n_inv) {
			for(;i<n_inv;i++) {
				if(i>n_peel) fputc(',',fptr);
				if(flag[i]&1) fputc('<',fptr);
				print_orig_id(fptr,inv[i]+1);
				if(flag[i]&1) fputc('>',fptr);
			}
			fprintf(fptr," [%d]",complex_em->out_index);
		} else fputc('*',fptr);
		fputc(' ',fptr);
		if(complex_em->n_rfuncs) {
			fputc('C',fptr);
			for(i=0;i<complex_em->n_rfuncs;i++) {
				fputc(i?',':'[',fptr);
				fprintf(fptr,"%d",complex_em->index[i]);
			}
			fputc(']',fptr);
		}
		fputc('\n',fptr);
	}
}
