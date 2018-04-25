/****************************************************************************
 *                                                                          *
 *     Loki - Programs for genetic analysis of complex traits using MCMC    *
 *                                                                          *
 *                    Simon Heath - CNG, France                             *
 *                                                                          *
 *                         August 2002                                      *
 *                                                                          *
 * peel_freq.c:                                                             *
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
#include <math.h>
#include <float.h>

#ifndef DBL_MAX
#define DBL_MAX MAXDOUBLE
#endif

#include "utils.h"
#include "loki.h"
#include "loki_peel.h"
#include "fenris_peel.h"
#include "fenris.h"
#include "min_deg.h"

void peel_freq(int locus,struct Peelseq_Head *pp,int *n_ops,int pen_type,double *err_probs)
{
	int i,j,k,grp,n_all,comp,nop,*a_trans,*aflag,nn_all;
	struct Peelseq_Head *pp1,**peel_list=0;
	struct Fenris_Simple_Element *simple_em;
   struct Complex_Element *complex_em; 
	struct pen_par ppar;
	fenris_pen_func *pen;
	double **freq,**nfreq,**nfreq1,**prior,*fq,*fq1,z,z1,like,z2,oldlike=-DBL_MAX;
	
	/* Find maximum number of peel_ops required, and allocate enough storage to keep them
	 * This allows us to 'replay' a peeling sequence to calculate expectations etc. */
	for(i=comp=0;comp<n_comp-singleton_flag;comp++) if(n_ops[comp]>i) i=n_ops[comp];
	if(i) if(!(peel_list=malloc(sizeof(void *)*i))) ABT_FUNC(MMsg);
	/* Number of alleles at locus */
	nn_all=marker[locus].locus.n_alleles;
 	if(!(aflag=malloc(sizeof(int)*nn_all))) ABT_FUNC(MMsg);
	/* Storage for current frequency estimates (required because of allele lumping) */
	if(!(freq=malloc(sizeof(void *)*4*n_genetic_groups))) ABT_FUNC(MMsg);
	nfreq=freq+n_genetic_groups;
	nfreq1=nfreq+n_genetic_groups;
	prior=nfreq1+n_genetic_groups;
	if(!(freq[0]=malloc(sizeof(double)*4*n_genetic_groups*nn_all))) ABT_FUNC(MMsg);
	for(i=1;i<n_genetic_groups*4;i++) freq[i]=freq[i-1]+nn_all;
	/* Initialize penetrance model */
	pen=get_pen_model(pen_type);
	ppar.e=err_probs;
	ppar.freq=freq;
	for(;;) {
		for(grp=0;grp<n_genetic_groups;grp++) {
			if(marker[locus].count_flag[grp]) {
				for(j=0;j<nn_all;j++) nfreq[grp][j]=marker[locus].counts[grp][j];
			} else for(j=0;j<nn_all;j++) nfreq[grp][j]=1.0;
			for(j=0;j<nn_all;j++) prior[grp][j]=nfreq[grp][j];
		}
		/* Clear R-Function pointers */
		for(i=0;i<ped_size;i++) id_array[i].rfp=-1;
		/* Loop through each component */
		like=0.0;
		for(comp=0;comp<n_comp-singleton_flag;comp++) {
			nop=0;
			pp1=pp+comp;
			/* Get number of alleles in this component */
			n_all=marker[locus].n_all1[comp];
			/* Calculate current frequency estimates, taking account of allele lumping 
			 * in this component */
			z=1.0;
			a_trans=allele_trans[locus][comp];
			for(i=0;i<nn_all;i++) aflag[i]=0;
			for(j=0;j<n_all-1;j++) aflag[a_trans[j]]=1;
			for(grp=0;grp<n_genetic_groups;grp++) {
				fq=marker[locus].locus.freq[grp];
				for(j=0;j<n_all-1;j++) {
					k=a_trans[j];
					freq[grp][j]=fq[k];
					z-=fq[k];
				}
				freq[grp][j]=z;
				for(j=0;j<n_all;j++) nfreq1[grp][j]=0.0;
			}
			/* Peel! */
			z=0.0;
			while(pp1->type) {
				print_peelseq_element(stdout,pp1);
				if(pp1->type==FENRIS_PEEL_SIMPLE) {
					simple_em=pp1->ptr.fsimple;
					if(simple_em->pivot) {
						j=simple_em->pivot-1;
						i=(id_array[j].rfp>=0)?1:0;
						peel_list[nop++]=pp1;
						z+=fenris_simple_peel(simple_em,locus,freq,pen,&ppar);
					} else z+=fenris_simple_distribute(simple_em,locus,freq,nfreq1,pen,&ppar);
					pp1=&simple_em->next;
				}  else if(pp1->type==PEEL_COMPLEX) {
					complex_em=pp1->ptr.complex;
					if(complex_em->n_involved!=complex_em->n_peel) {
						peel_list[nop++]=pp1;
						z+=fenris_complex_peel(complex_em,locus,freq,pen,&ppar);
					} else z+=fenris_complex_distribute(complex_em,locus,freq,nfreq1,pen,&ppar);
					pp1=&complex_em->next;
				} else ABT_FUNC("Unknown peel type\n");
			}
			printf("Comp %d, like %.10f\n",comp,z);
			like+=z;
			/* Distribute evidence */
			for(i=nop-1;i>=0;i--) {
				pp1=peel_list[i];
				if(pp1->type==FENRIS_PEEL_SIMPLE) {
					simple_em=pp1->ptr.fsimple;
					(void)fenris_simple_distribute(simple_em,locus,freq,nfreq1,pen,&ppar);
				}  else {
					complex_em=pp1->ptr.complex;
					(void)fenris_complex_distribute(complex_em,locus,freq,nfreq1,pen,&ppar);
				}
			}
			/* Add frequency estimates */
			for(grp=0;grp<n_genetic_groups;grp++) {
				z=1.0;
				fq=marker[locus].locus.freq[grp];
				for(j=0;j<n_all-1;j++) {
					k=a_trans[j];
					nfreq[grp][k]+=nfreq1[grp][j];
					z-=fq[k];
				}
				z1=nfreq1[grp][j];
				for(j=0;j<nn_all;j++) if(!aflag[j]) nfreq[grp][j]+=z1*fq[j]/z;
			}
		}
		z1=0.0;
		for(grp=0;grp<n_genetic_groups;grp++) {
			fq=marker[locus].locus.freq[grp];
			fq1=nfreq[grp];
			z=z2=0.0;
			for(z=0.0,i=0;i<nn_all;i++) {
				z+=fq1[i];
				z2+=prior[grp][i];
				z1-=lgamma(prior[grp][i]);
			} 
			z1+=lgamma(z2);
			z=1.0/(z-(double)nn_all);
			for(i=0;i<nn_all;i++) {
				fq[i]=(fq1[i]-1.0)*z;
				z1+=log(fq[i])*(prior[grp][i]-1.0);
/*				printf("BB: %d %.10f %.10f\n",i,fq[i],fq1[i]); */
			}
		}
		printf("like1=%.10f, z1=%g, like=%.10f\n",like,z1,like+z1);
		if((like-oldlike)<1.0e-12) {
			if(oldlike-like>1.0e-5) fprintf(stderr,"%g OOOK!\n",oldlike-like);
			break;
		}
		oldlike=like;
		break;
	}
	for(i=0;i<nn_all;i++) {
		printf("BB: %d %g\n",i,marker[locus].locus.freq[0][i]);
	}
	if(peel_list) free(peel_list);
	min_deg(0,0,0,0);
	free(aflag);
	free(freq[0]);
	free(freq);
}
