/****************************************************************************
 *                                                                          *
 *     Loki - Programs for genetic analysis of complex traits using MCMC    *
 *                                                                          *
 *                    Simon Heath - CNG, France                             *
 *                                                                          *
 *                         September 2002                                   *
 *                                                                          *
 * fenris_complex_peel.c:                                                   *
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

#include "utils.h"
#include "loki.h"
#include "loki_peel.h"
#include "fenris_peel.h"
#include "fenris.h"

#define GT(x,y) ((x)>(y)?(x)*((x)+1)/2+(y):(y)*((y)+1)/2+(x))

struct triplet_peel **fams;
static int *kids,max_kids,max_fam,*trans_all1,*trans_all2;

void setup_fenris_complex_peel()
{
	int i,j,k,n_all,n_gen;
	
	n_all=0;
	for(i=0;i<n_markers;i++) {
		j=marker[i].locus.n_alleles;
		if(j>n_all) n_all=j;
	}
	n_gen=n_all*(n_all+1)/2;
	if(!(trans_all1=malloc(sizeof(int)*n_gen*2))) ABT_FUNC(MMsg);
	trans_all2=trans_all1+n_gen;
	for(i=k=0;i<n_all;i++) {
		for(j=0;j<=i;j++) {
			trans_all1[k]=i;
			trans_all2[k++]=j;
		}
	}
}

void free_fenris_complex_peel(void) 
{
	int i;
	
	if(kids) {
		free(kids);
		kids=0;
	}
	if(fams) {
		for(i=0;i<max_fam;i++) free(fams[i]);
		free(fams);
		fams=0;
	}
	if(trans_all1) {
		free(trans_all1);
		trans_all1=0;
	}
}

static int set_up_families(int n,int *inv,int *flags)
{
	int i,j,k,ids,idd,nfam=0,nkids=0,*ti;
	struct triplet_peel *t;
	
	/* Set up family triplets */
	for(i=0;i<n;i++) {
		flags[i]&=1;
		if(flags[i]) nkids++;
	}
	if(nkids) {
		/* Allocate memory for kids */
		if(nkids>max_kids) {
			if(kids) free(kids);
			if(!(kids=malloc(sizeof(int)*nkids))) ABT_FUNC(MMsg);
			max_kids=nkids;
		}
		/* Count number of families */
		for(i=0;i<n;i++) if((flags[i]&3)==1) {
			ids=id_array[inv[i]].sire-1;
			idd=id_array[inv[i]].dam-1;
			nfam++;
			for(j=i+1;j<n;j++) if((flags[j]&3)==1) {
				k=inv[j];
				if(id_array[k].sire-1==ids && id_array[k].dam-1==idd) flags[j]|=2;
			}
		}
		/* Allocate space if required */
		if(nfam>max_fam) {
			if(fams) {
				if(!(fams=realloc(fams,sizeof(void *)*nfam))) ABT_FUNC(MMsg);
			} else if(!(fams=malloc(sizeof(void *)*nfam))) ABT_FUNC(MMsg);
			for(;max_fam<nfam;max_fam++) {
				if(!(fams[max_fam]=malloc(sizeof(struct triplet_peel)))) ABT_FUNC(MMsg);
			}
		}
		ti=kids;
		/* Populate family structures */
		for(nfam=i=0;i<n;i++) if((flags[i]&5)==1) {
			ids=id_array[inv[i]].sire-1;
			idd=id_array[inv[i]].dam-1;
			fams[nfam]->nkids=1;
			fams[nfam]->kids=ti;
			*ti++=i;
			for(j=0;j<n;j++) if(i!=j) {
				if(inv[j]==idd) fams[nfam]->idd=j;
				else if(inv[j]==ids) fams[nfam]->ids=j;
				else if(j>i && (flags[j]&1)) {
					if(id_array[inv[j]].sire-1==ids && id_array[inv[j]].dam-1==idd) {
						*ti++=j;
						fams[nfam]->nkids++;
						flags[j]|=4;
					}
				}
			}
			flags[fams[nfam]->ids]|=16;
			flags[fams[nfam]->idd]|=16;
			nfam++;
		}
		if(nfam>1) {
			/* Sort on sire ids - bubble sort, but since normally nfam<=1, this should not be a problem! 
			 * Only required in the case of 'overlapping families' - where a kid in one family is a parent
			 * in another family.  In this case we want the family where the individual is a kid to come
			 * before the other one.  We can effectively do this by sorting on sire id */
			for(i=nfam-1;i>=0;i--) {
				for(k=j=0;j<i;j++) {
					if(fams[j]->ids>fams[j+1]->ids) {
						t=fams[j];
						fams[j]=fams[j+1];
						fams[j+1]=t;
						k++;
					}
				}
				if(!k) break;
			}
		}
		for(i=0;i<nfam;i++) {
			fputs("Fam: ",stdout);
			print_orig_id(stdout,inv[fams[i]->ids]+1);
			fputc(' ',stdout);
			print_orig_id(stdout,inv[fams[i]->idd]+1);
			fputc(' ',stdout);
			for(j=0;j<fams[i]->nkids;j++) {
				print_orig_id(stdout,inv[fams[i]->kids[j]]+1);
				fputc(' ',stdout);
			}
			fputc('\n',stdout);
		}
	}
	return nfam;
}

static void update_fam(int j,int *c,int *nn,int *flags,int *fc_perm,int nfam,int *fam_idx)
{
	int k,k1,ids,idd,g1,g2,g3,g4,*gg,kid;
	double *zz;
	
	ids=fams[j]->ids;
	idd=fams[j]->idd;
	k=fc_state[ids];
	if(flags[fc_perm[ids]]&1) k=fams[fam_idx[ids]]->g[k];
	g1=trans_all1[k];
	g2=trans_all2[k];
	k=fc_state[idd];
	if(flags[fc_perm[idd]]&1) k=fams[fam_idx[idd]]->g[k];
	g3=trans_all1[k];
	g4=trans_all2[k];
	zz=fams[j]->z;
	gg=fams[j]->g;
	if(g1==g2) {
		if(g3==g4) {
			fams[j]->nc=1;
			gg[0]=GT(g1,g3);
			zz[0]=1.0;
		} else {
			fams[j]->nc=2;
			gg[0]=GT(g1,g3);
			gg[1]=GT(g1,g4);
			zz[0]=zz[1]=0.5;
		}
	} else if(g3==g4) {
		fams[j]->nc=2;
		gg[0]=GT(g1,g3);
		gg[1]=GT(g2,g3);
		zz[0]=zz[1]=0.5;
	} else if(g1==g3 && g2==g4) {
		fams[j]->nc=3;
		gg[0]=GT(g1,g3);
		gg[1]=GT(g1,g4);
		gg[2]=GT(g2,g4);
		zz[0]=zz[2]=0.25;
		zz[1]=0.5;
	} else {
		fams[j]->nc=4;
		gg[0]=GT(g1,g3);
		gg[1]=GT(g1,g4);
		gg[2]=GT(g2,g3);
		gg[3]=GT(g2,g4);
		zz[0]=zz[1]=zz[2]=zz[3]=0.25;
	}
	for(k=0;k<fams[j]->nkids;k++) {
		kid=fams[j]->kids[k];
		fc_state[kid]=0;
		if(flags[fc_perm[kid]]&16) {
			for(k1=j+1;k1<nfam;k1++) {
				if(fams[k1]->ids==kid || fams[k1]->idd==kid) update_fam(k1,c,nn,flags,fc_perm,nfam,fam_idx);
			}
		}
		nn[kid]=fams[j]->nc;
		c[kid]=nn[kid]-1;
	}
}

double fenris_complex_peel(struct Complex_Element *element,int locus,double **freq,fenris_pen_func *pen,struct pen_par *ppar)
{
	int i,j,k,id,*inv,*flags,n_inv,n_peel,n_out,n_rf,*index,nfam,start;
	int nc,ids,idd;
	int n_all,n_gen,comp,kid,*fc_perm,*c,*t,*nn,*fam_idx;
	double z,prob=0.0;

	inv=element->involved;
	n_inv=element->n_involved;
	n_peel=element->n_peel;
	n_out=n_inv-n_peel;
	flags=element->flags;
	n_rf=element->n_rfuncs;
	comp=id_array[inv[0]].comp;
	n_all=marker[locus].n_all1[comp];
	n_gen=n_all*(n_all+1)/2;
	index=element->index;
	fc_perm=fc_state+n_inv;
	c=fc_perm+n_inv;
	t=c+n_inv;
	nn=t+n_inv;
	fam_idx=nn+n_inv;
	/* First set up (single person) functions for individuals to be peeled */
	for(i=0;i<n_peel;i++) {
		id=inv[i];
		z=get_par_probs(kval[i],id,locus,pen,ppar,freq);
		prob+=log(z);
	}
	/* Find family structures that need peeling */
	nfam=set_up_families(n_inv,inv,flags);
	/* Build permutation of inv so that parents come first, and kids come just after parents */
	j=n_inv;
	for(i=0;i<n_inv;i++) fam_idx[i]=-1;
	for(i=0;i<nfam;i++) {
		ids=fams[i]->ids;
		idd=fams[i]->idd;
		if(!(flags[ids]&8)) {
			fc_perm[--j]=ids;
			fams[i]->ids=j;
			flags[ids]|=8;
		} else {
			for(k=j;k<n_inv;k++) if(fc_perm[k]==ids) break;
			fams[i]->ids=k;
		}
		if(!(flags[idd]&8)) {
			fc_perm[--j]=idd;
			fams[i]->idd=j;
			flags[idd]|=8;
		} else {
			for(k=j;k<n_inv;k++) if(fc_perm[k]==idd) break;
			fams[i]->idd=k;
		}
		for(k=0;k<fams[i]->nkids;k++) {
			kid=fams[i]->kids[k];
			if(flags[kid]&8) ABT_FUNC("Shouldn't happen\n");
			fc_perm[--j]=kid;
			fam_idx[j]=i;
			fams[i]->kids[k]=j;
			flags[kid]|=8;
		}
	}
	for(i=0;i<n_inv;i++) if(!(flags[i]&8)) fc_perm[--j]=i;
	fputs("fc_perm: ",stdout);
	for(i=0;i<n_inv;i++) {
		if(i) fputc(' ',stdout);
		print_orig_id(stdout,inv[fc_perm[i]]+1);
	}
	fputc('\n',stdout);
	/* Initialize states */
	for(i=0;i<n_inv;i++) {
		j=fc_perm[i];
		fc_state[i]=0;
		t[i]=i;
		nn[i]=(flags[j]&1)?1:n_gen; /* No. states is ngen apart from kids who have one possible state initially */
		c[i]=nn[i]-1;
	}
	for(i=0;i<nfam;i++) {
		fams[i]->nc=1;		  
		fams[i]->g[0]=0;
		fams[i]->z[0]=1.0;
	}
	start=0;
	while(start<n_inv && !c[start]) start++;
#ifdef DEBUG
	if(i==n_inv) ABT_FUNC("No variable sites\n");
#endif
	t[0]=start;
	nc=0;
	return;
	/* Main loop */
	for(;;) {
		nc++;
		for(j=0;j<n_inv;j++) {
			k=fc_state[j];
			if(flags[fc_perm[j]]&1) k=fams[fam_idx[j]]->g[k];
			printf("%d ",k);
		}
		printf("\n");
		k=t[start];
		t[start]=start;
		if(!c[k]) break;
		c[k]--;
		if(!c[k]) {
			j=k+1;
			while(j<n_inv && nn[j]==1) j++;
			if(j<n_inv) {
				t[k]=t[j];
				t[j]=j;
				c[k]=nn[k]-1;
			}
		}
		fc_state[k]=(fc_state[k]+1)%nn[k];
		if(flags[fc_perm[k]]&16) { /* This is a parent, so we must update kids */
			for(j=0;j<nfam;j++) {
				if(fams[j]->ids==k || fams[j]->idd==k) update_fam(j,c,nn,flags,fc_perm,nfam,fam_idx);
			}
			k=0;
			while(k<n_inv && nn[k]==1) k++;
			for(j=0;j<=k && j<start;j++) t[j]=j;
			if(k!=start) {
				if(k<start) {
					if(t[start]==start) t[k]=k;
					t[start]=start;
				} 
				start=k;
			}
		}
	}
	printf("nc=%d\n",nc);
	return prob;
}

double fenris_complex_distribute(struct Complex_Element *element,int locus,double **freq,double **nfreq,fenris_pen_func *pen,struct pen_par *ppar)
{
	return 0.0;
}
