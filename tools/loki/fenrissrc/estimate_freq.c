/****************************************************************************
 *                                                                          *
 *     Loki - Programs for genetic analysis of complex traits using MCMC    *
 *                                                                          *
 *                    Simon Heath - CNG, France                             *
 *                                                                          *
 *                         August 2002                                      *
 *                                                                          *
 * estimate_freq.c:                                                         *
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
#include "fenris.h"
#include "fenris_peel.h"

/* Get initial estimate of allele frequencies by simply counting alleles */
static void count_alleles(double **counts,int locus)
{
	struct Marker *mark;
	int i,j,k,k1,k2,n_all,comp,cs,*a_trans,grp;
	double *fq,n;
	
	mark=marker+locus;
	n_all=mark->locus.n_alleles;
	for(grp=0;grp<n_genetic_groups;grp++) {
		if(mark->count_flag[grp]) for(i=0;i<n_all;i++) counts[grp][i]=mark->counts[grp][i];
		else for(i=0;i<n_all;i++) counts[grp][i]=1.0;
	}
	for(i=comp=0;comp<n_comp;comp++) {
		cs=comp_size[comp];
		a_trans=allele_trans[locus][comp];
		for(j=0;j<cs;j++,i++) {
			k=mark->haplo[i];
			if(k) {
				grp=id_array[i].group-1; /* Individual's genetic group */
				k1=a_trans[(k>>16)-1]; /* Translate back to original (before downcoding) allele ids */
				k2=a_trans[(k&65535)-1];
#ifdef DEBUG
				if(grp<0 || grp>n_genetic_groups) ABT_FUNC("Illegal group id\n");
				if(k1<0 || k1>=n_all || k2<0 || k2>=n_all) ABT_FUNC("Illegal allele values\n");
#endif
				counts[grp][k1]++;
				counts[grp][k2]++;
			}
		}
	}
	/* Get estimates */
	for(grp=0;grp<n_genetic_groups;grp++) {
		for(n=0.0,j=0;j<n_all;j++) n+=counts[grp][j];
		n=1.0/(double)n;
		fq=mark->locus.freq[grp];
		for(j=0;j<n_all;j++) fq[j]=counts[grp][j]*n;
	}
}

void estimate_freq(int pen_type,double *err_probs)
{
	int locus,n_all,i,j,k,interactive,*nops;
	struct Marker *mark;
	struct Peelseq_Head *pp;
	double **counts;
	
	/* Check if any estimation required */
	for(locus=0;locus<n_markers;locus++) {
		mark=marker+locus;
		n_all=mark->locus.n_alleles;
		for(j=k=0;k<n_genetic_groups;k++) {
			for(i=0;i<n_all;i++) if(mark->freq_set[k][i]!=1) j++;
			if(j) break;
		}
		if(j) break;
	}
	if(locus==n_markers) return; /* No, every frequency has been specified */
	fputs("Estimating allele frequencies\n",stdout);
	for(locus=i=0;locus<n_markers;locus++) {
		n_all=marker[locus].locus.n_alleles;
		if(n_all>i) i=n_all;
	}
	if(!(nops=malloc(sizeof(int)*n_comp))) ABT_FUNC(MMsg);
	pp=peel_init(nops);
	fenris_peel_alloc(pp);
	/* Allocate space for counts */
	if(!(counts=malloc(sizeof(void *)*n_genetic_groups))) ABT_FUNC(MMsg);
	if(!(counts[0]=malloc(sizeof(double)*n_genetic_groups*i))) ABT_FUNC(MMsg);
	for(i=1;i<n_genetic_groups;i++) counts[i]=counts[i-1]+i;
	interactive=isatty(fileno(stdout));
	for(locus=0;locus<n_markers;locus++) {
		mark=marker+locus;
		n_all=mark->locus.n_alleles;
		for(j=k=0;k<n_genetic_groups;k++) {
			for(i=0;i<n_all;i++) if(mark->freq_set[k][i]!=1) j++;
			if(j) break;
		}
		if(!j) continue;
		if(interactive) {
			print_marker_name(stdout,locus);
			fputs(": Getting initial estimate (allele counting)\r",stdout);
			fflush(stdout);
		}
		count_alleles(counts,locus);
		printf("\n");
		for(i=0;i<n_all;i++) {
			printf("%d %g\n",i,marker[locus].locus.freq[0][i]);
		}
		if(interactive) {
			print_marker_name(stdout,locus);
			fputs(": Getting ML estimate (peeling)             \r",stdout);
			fflush(stdout);
		}
		peel_freq(locus,pp,nops,pen_type,err_probs);
		if(interactive) {
			print_marker_name(stdout,locus);
			fputs(": Done                                      \r",stdout);
			fflush(stdout);
		}
	}
 	fenris_peel_free();
	free_peel_sequence(pp);
	free(nops);
	free(counts[0]);
	free(counts);
	if(interactive) fputc('\n',stdout);
}
