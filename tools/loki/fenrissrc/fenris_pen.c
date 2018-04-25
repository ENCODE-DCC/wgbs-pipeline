/****************************************************************************
 *                                                                          *
 *     Loki - Programs for genetic analysis of complex traits using MCMC    *
 *                                                                          *
 *                    Simon Heath - CNG, France                             *
 *                                                                          *
 *                         August 2002                                      *
 *                                                                          *
 * fenris_pen.c:                                                            *
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

/* Penetrance for error model with equal error probs:
 * p(g|y) = 1-e + e/m   if g==y
 *          e/m if g!=y
 * where m is the number of possible genotypes
 */
double fpen_emodel_equal(double *p,int id,int locus,struct pen_par *par)
{
	int i,nn_all,n_all,comp,ngen,ng1,g1,g2,k;
	double e,z,pp=0.0;
	struct Marker *mark;

	mark=marker+locus;
	k=mark->haplo[id];
	if(!k) return 1.0;
	g1=(k>>16)-1;
	g2=(k&65535)-1;
#ifdef DEBUG
	/* This should be enforced by read_binfiles.c */
	if(g1<g2) ABT_FUNC("Bad order for genotype\n"); 
#endif
	k=g1*(g1+1)/2+g2;
	e=par->e[0];
	comp=id_array[id].comp;
	/* No. markers and genotypes for locus */
	nn_all=marker[locus].locus.n_alleles; 
	ngen=nn_all*(nn_all+1)/2;
	/* No. alleles in this component */
	n_all=marker[locus].n_all1[comp];
	ng1=n_all*(n_all+1)/2;
	z=e/(e+(1.0-e)*(double)ngen);
	for(i=0;i<k;i++) pp+=(p[i]*=z); /* genotypes < observed */
	i++; /* observed genotype */
	pp++;
	for(;i<ng1;i++) pp+=(p[i]*=z); /* genotypes > observed */
	return pp;
}

/* Penetrance for error model with error probs proportional to genotype frequencies:
 * p(g|y) = 1-e   if g==y
 *          e*Pr(y)/(1-Pr(g)) if g!=y
 *
 * where Pr(x) is the genotype frequency of x
 */
double fpen_emodel_prop(double *p,int id,int locus,struct pen_par *par)
{
	int i,j,k,k1,n_all,comp,grp,g1,g2;
	double e,*freq,z,z1,pp=0.0;
	struct Marker *mark;
	
	mark=marker+locus;
	k=mark->haplo[id];
	if(!k) return 1.0;
	g1=(k>>16)-1;
	g2=(k&65535)-1;
#ifdef DEBUG
	/* This should be enforced by read_binfiles.c */
	if(g1<g2) ABT_FUNC("Bad order for genotype\n"); 
#endif
	k=g1*(g1+1)/2+g2;
	e=par->e[0];
	comp=id_array[id].comp;
	grp=id_array[id].group-1;
	/* No. alleles in this component */
	n_all=marker[locus].n_all1[comp];
	/* No. genotypes */
	freq=par->freq[grp]; 
	z=(g1==g2)?freq[g1]*freq[g1]:freq[g1]*freq[g2]*2.0;
	z1=e/(1.0-z);
	for(k1=i=0;i<g1;i++) {
		z=freq[i];
		for(j=0;j<i;j++) p[k1++]*=z1*2.0*z*freq[j];
		pp+=(p[k1++]*=z1*z*z);
	}
	z=freq[i];
	for(j=0;j<g2;j++) pp+=(p[k1++]*=z1*2.0*z*freq[j]);
	pp+=(p[k1++]*=1.0-e);
	for(j=g2+1;j<i;j++) pp+=(p[k1++]*=z1*2.0*z*freq[j]);
	pp+=(p[k1++]*=z1*z*z);
	for(i=g1+1;i<n_all;i++) {
		z=freq[i];
		for(j=0;j<i;j++) pp+=(p[k1++]*=z1*2.0*z*freq[j]);
		pp+=(p[k1++]*=z1*z*z);
	}
	return pp;
}

/* Penetrance for error model with empirical error probs (Sobel, Papp & Lange 2002):
 */
double fpen_emodel_emp(double *p,int id,int locus,struct pen_par *par)
{
	int i,j,k,k1,nn_all,n_all,comp,ngen,g1,g2;
	double *e,n1,n2,n12,n23=0.0,z,pp=0.0;
	struct Marker *mark;
	
	mark=marker+locus;
	k=mark->haplo[id];
	if(!k) return 1.0;
	g1=(k>>16)-1;
	g2=(k&65535)-1;
#ifdef DEBUG
	/* This should be enforced by read_binfiles.c */
	if(g1<g2) ABT_FUNC("Bad order for genotype\n"); 
#endif
	e=par->e;
	comp=id_array[id].comp;
	/* No. markers and genotypes for locus */
	nn_all=marker[locus].locus.n_alleles; 
	ngen=nn_all*(nn_all+1)/2;
	n1=1.0/(double)(nn_all-1);
	if(nn_all>2) {
		n2=1.0/(double)(nn_all-2);
		n12=2.0/(double)((nn_all-1)*(nn_all-2));
		if(nn_all>3) n23=2.0/(double)((nn_all-3)*(nn_all-2));
	}
	/* No. alleles in this component */
	n_all=marker[locus].n_all1[comp];
	if(n_all<nn_all) n_all--;
	if(g1==g2) { /* Homozygous observed genotype */
		z=(nn_all==2)?1.0-e[0]-e[2]:1.0-e[0]-e[2]-e[4];
		for(k1=i=0;i<n_all;i++) {
			for(j=0;j<i;j++) {
				if(g1==i || g1==j) pp+=(p[k1++]*=e[0]*n1);
				else pp+=(p[k1++]*=n12*e[4]);
			}
			if(g1==i) pp+=(p[k1++]*=z);
			else pp+=(p[k1++]*=e[2]*n1);
		}
	} else { /* Heterozygous observed genotype */
		if(nn_all==2) z=1.0-e[3];
		else if(nn_all==3) z=1.0-e[1]-e[3]-e[4];
		else z=1.0-e[1]-e[2]-e[3]-e[4];
		for(k1=i=0;i<n_all;i++) {
			for(j=0;j<i;j++) {
				if((g1==i && g2==j)||(g1==j && g2==i)) pp+=(p[k1++]*=z);
				else {
					if(g1==i || g1==j || g2==i || g2==j) pp+=(p[k1++]*=.5*e[1]*n2);
					else pp+=(p[k1++]*=e[2]*n23);
				}
			}
			if(g1==i || g2==i) pp+=(p[k1++]*=0.5*e[3]);
			else pp+=(p[k1++]*=e[4]*n2);
		}
	}
	if(n_all<nn_all) { /* Lumping */
		z=(double)(nn_all-n_all);
		if(g1==g2) { /* Homozygous observed genotype */
			for(j=0;j<i;j++) {
				if(g1==j) pp+=(p[k1++]*=z*e[0]*n1);
				else pp+=(p[k1++]*=z*n12*e[4]);
			}
			pp+=(p[k1++]*=(z*e[2]*n1+.5*z*(z-1.0)*n12*e[4]));
		} else { /* Heterozygous observed genotype */
			for(j=0;j<i;j++) {
				if(g1==j || g2==j) pp+=(p[k1++]*=z*.5*e[1]*n2);
				else pp+=(p[k1++]*=z*e[2]*n23);
			}
			pp+=(p[k1++]*=(z*e[4]*n2+.5*z*(z-1.0)*e[2]*n23));
		}
	}
	return pp;
}
