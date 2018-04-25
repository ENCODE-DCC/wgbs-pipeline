/****************************************************************************
*                                                                          *
*     Loki - Programs for genetic analysis of complex traits using MCMC    *
*                                                                          *
*             Simon Heath - University of Washington                       *
*                                                                          *
*                       July 1997                                          *
*                                                                          *
*         Massively updated by Simon Heath at MSKCC, June 2001             *
* seg_pen.c:                                                               *
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
#include <math.h>
#include <stdio.h>
#include <assert.h>

#include "utils.h"
#include "loki.h"
#include "ranlib.h"
#include "loki_peel.h"
#include "handle_res.h"
#include "seg_pen.h"
#include "lk_malloc.h"
#include "gen_pen.h"

static double **seg_freq;
static int *group,*first,*next,*alleles[2],*gpflag,*gpsize,**obslist,**cpstart,*fg_type;
static int *gplist,*gplist1,ngroups,*seg_alls,*lookup,*lookup1,*lookup_a,*lookup1_a;
static double *pp[2];

static void setup_groups(struct loki *loki)
{
	int i,i1,j,k,comp,locus,nfd=0,ntl,n_comp,n_markers;
	int **temp_p,*temp_p1;
	struct Id_Record *id_array;
	struct Locus *loc;
	
	id_array=loki->pedigree->id_array;
	ntl=loki->params.n_tloci?1:0;
	n_comp=loki->pedigree->n_comp;
	n_markers=loki->markers->n_markers;
	for(locus=0;locus<n_markers+ntl;locus++) {
		loc=(locus<n_markers)?&loki->markers->marker[locus].locus:loki->models->tlocus;
		for(i=comp=0;comp<n_comp;comp++) {
			for(j=0;j<loki->pedigree->comp_size[comp];j++) {
				i1=i+j;
				k=id_array[i1].sire;
				if(!k || loc->pruned_flag[k-1]) nfd+=2;
			}
			i+=j;
		}
	}
	if(!nfd) return;
	temp_p=lk_malloc(sizeof(void *)*(n_markers+ntl)*n_comp);
	loki->sys.RemBlock=AddRemem(temp_p,loki->sys.RemBlock);
	temp_p1=lk_malloc(sizeof(int)*nfd);
	loki->sys.RemBlock=AddRemem(temp_p1,loki->sys.RemBlock);
	for(locus=0;locus<n_markers+ntl;locus++) {
		if(locus<n_markers) {
			loki->markers->marker[locus].group=temp_p;
			loc=&loki->markers->marker[locus].locus;
		} else {
			loki->peel->tl_group=temp_p;
			loc=loki->models->tlocus;
		}
		for(i=comp=0;comp<n_comp;comp++)	{
			temp_p[comp]=temp_p1;
			for(nfd=j=0;j<loki->pedigree->comp_size[comp];j++) {
				i1=i+j;
				k=id_array[i1].sire;
				if(!k || loc->pruned_flag[k-1]) {
					temp_p1[nfd++]=id_array[i1].group-1;
					temp_p1[nfd++]=id_array[i1].group-1;
					if(!id_array[i1].group) {
						ABT_FUNC("Internal error: null group\n");
					}
				}
			}
			i+=j;
			temp_p1+=nfd;
		}
		temp_p+=n_comp;
	}
}

void seg_dealloc(struct loki *lk)
{
	int i;
	
	if(cpstart && cpstart[0]) {
		free(cpstart[0]);
		free(cpstart);
	}
	if(obslist) {
		for(i=0;i<lk->markers->n_markers+(lk->params.n_tloci?1:0);i++) if(obslist[i]) free(obslist[i]);
		free(obslist);
	}
	if(first) free(first);
	if(seg_freq) {
		if(seg_freq[0]) free(seg_freq[0]);
		free(seg_freq);
	}
	if(lk->peel->aff_freq) free(lk->peel->aff_freq);
	if(fg_type) free(fg_type);
	if(pp[0]) free(pp[0]);
	if(seg_alls) free(seg_alls);
	if(lookup) free(lookup);
	free_gen_pen();
	pass_founder_genes_dealloc();
}

void setup_obslist(struct loki *loki)
{
	int i,j,k,x,comp,nl,n_markers,n_comp,*comp_size;
	struct Marker *marker;
	struct Id_Record *id_array;
	
	marker=loki->markers->marker;
	comp_size=loki->pedigree->comp_size;
	id_array=loki->pedigree->id_array;
	n_markers=loki->markers->n_markers;
	n_comp=loki->pedigree->n_comp;
	if(loki->params.n_tloci) k=2;
	nl=n_markers+(loki->params.n_tloci?1:0);
	for(x=0;x<nl;x++) {
		for(i=j=comp=0;comp<n_comp;comp++) {
			cpstart[x][comp]=j;
			if(x<n_markers) {
				if(marker[x].mterm && marker[x].mterm[0]) {
					for(k=0;k<comp_size[comp];k++,i++) if(marker[x].haplo[i] || id_array[i].res[0]) j++;
				} else {
					for(k=0;k<comp_size[comp];k++,i++) if(marker[x].haplo[i]) j++;
				}
			} else {
				if(id_array[i].res) for(k=0;k<comp_size[comp];k++,i++) if(id_array[i].res[0]) j++;
			}
			j++;
		}
		cpstart[x][comp]=j;
		if(j) {
			if(obslist[x]) free(obslist[x]);
			obslist[x]=lk_malloc(sizeof(int)*j);
			for(i=j=comp=0;comp<n_comp;comp++) {
				if(x<n_markers) {
					if(marker[x].mterm && marker[x].mterm[0]) {
						for(k=0;k<comp_size[comp];k++,i++) if(marker[x].haplo[i] || id_array[i].res[0]) {
							obslist[x][j++]=i;
						}
					} else {
						for(k=0;k<comp_size[comp];k++,i++) if(marker[x].haplo[i]) {
							obslist[x][j++]=i;
						}
					}
				} else {
					for(k=0;k<comp_size[comp];k++,i++) if(id_array[i].res && id_array[i].res[0]) obslist[x][j++]=i;
				}
						obslist[x][j++]=-1;
			}
		} else obslist[x]=0;
	}
}

void seg_alloc(struct loki *loki)
{
	int i,j=0,k=0,k1,k2,nl,n_markers,n_genetic_groups,n_comp,*comp_size;
	struct Marker *marker;
	struct Id_Record *id_array;
	
	n_genetic_groups=loki->pedigree->n_genetic_groups;
	if(n_genetic_groups>1) setup_groups(loki);
	marker=loki->markers->marker;
	comp_size=loki->pedigree->comp_size;
	id_array=loki->pedigree->id_array;
	n_markers=loki->markers->n_markers;
	n_comp=loki->pedigree->n_comp;
	if(loki->params.n_tloci) k=2;
	for(i=0;i<n_markers;i++) {
		j=loki->markers->marker[i].locus.n_alleles;
		if(j>k) k=j;
	}
	if(k) {
		seg_alls=lk_malloc(sizeof(int)*k);
		seg_freq=lk_malloc(sizeof(void *)*2*n_genetic_groups);
		loki->peel->seg_count=seg_freq+n_genetic_groups;
		seg_freq[0]=lk_malloc(sizeof(double)*k*2*n_genetic_groups);
		for(i=1;i<n_genetic_groups*2;i++) seg_freq[i]=seg_freq[i-1]+k;
		if(loki->params.est_aff_freq) {
			loki->peel->aff_freq=lk_malloc(sizeof(double)*k);
		}
	}
	for(k=k1=i=0;i<n_comp;i++) {
		j=comp_size[i];
		if(j>k) k=j;
		j=loki->pedigree->comp_ngenes[i];
		if(j>k1) k1=j;
	}
	if(loki->params.est_aff_freq) fg_type=lk_malloc(sizeof(int)*k1);
	if(loki->params.n_tloci) alloc_gen_pen(k1);
	first=lk_malloc(sizeof(int)*k1*9);
	next=first+k1; /* First and next keep track of group members */
	group=next+k1; /* Group membership for each gene */
	gpflag=group+k1; /* Whether a group has 2 or 1 possible states */
	gpsize=gpflag+k1; /* Size of group */
	alleles[0]=gpsize+k1; /* Allele allocations for each gene */
	alleles[1]=alleles[0]+k1;
	gplist=alleles[1]+k1; /* List of active groups */
	gplist1=gplist+k1; /* Inverse of gplist */
	nl=n_markers+(loki->params.n_tloci?1:0);
	if(nl) {
		cpstart=lk_malloc(sizeof(void *)*nl);
		cpstart[0]=lk_malloc(sizeof(int)*nl*(n_comp+1));
		for(i=1;i<nl;i++) cpstart[i]=cpstart[i-1]+n_comp+1;
		obslist=lk_malloc(sizeof(void *)*nl);
		for(i=0;i<nl;i++) obslist[i]=0;
		setup_obslist(loki);
	}
	pp[0]=lk_malloc(sizeof(double)*k1*2);
	pp[1]=pp[0]+k1;
	pass_founder_genes_alloc();
	k=2;
	for(i=0;i<n_markers;i++) {
		if(marker[i].locus.n_alleles>k) k=marker[i].locus.n_alleles;
	}
	i=(k+1)*(k+2)/2;
	j=(k+1)*k/2;
 	lookup=lk_malloc(sizeof(int)*(i+j)*2);
	lookup1=lookup+i;
	lookup_a=lookup1+i;
	lookup1_a=lookup_a+j;
	for(k1=k2=i=0;i<=k;i++) {
		for(j=0;j<=i;j++) {
			lookup[k1]=i;
			lookup1[k1++]=j;
			if(i && j) {
				lookup_a[k2]=i;
				lookup1_a[k2++]=j;
			}
		}
	}
}

/* Change group g1 to group g2, link groups together, and switch allele allocations */
static void change_and_switch_group(int g1,int g2)
{
	int k,k1,k2;
	
	/* Change group and allele allocations for members of g1 */
	k1=gpflag[g2]-1;
	k2=k1^1;
	k=first[g1];
	while(next[k]>=0) {
		group[k]=g2;
		alleles[k1][k]=alleles[k2][k];
		k=next[k];
	}
	group[k]=g2;
	alleles[k1][k]=alleles[k2][k];
	/* Add members of g1 to g2 list */
	next[k]=first[g2];
	first[g2]=first[g1];
	/* Update log probabilities for group g2 */
	pp[k1][g2]+=pp[k2][g1];
	/* Update size of g2 */
	gpsize[g2]+=gpsize[g1];
	/* Remove g1 from active list */
	k=gplist1[g1];
	k1=gplist[k]=gplist[--ngroups];
	gplist1[k1]=k;
}

/* Change group g1 to group g2, link groups together, and swap allele allocations */
static void change_and_swap_group(int g1,int g2)
{
	int k,k1;
	
	/* Change group and allele allocations for members of g1 */
	k=first[g1];
	while(next[k]>=0) {
		group[k]=g2;
		k1=alleles[0][k];
		alleles[0][k]=alleles[1][k];
		alleles[1][k]=k1;
		k=next[k];
	}
	group[k]=g2;
	k1=alleles[0][k];
	alleles[0][k]=alleles[1][k];
	alleles[1][k]=k1;
	/* Add members of g1 to g2 list */
	next[k]=first[g2];
	first[g2]=first[g1];
	/* Update log probabilities for group g2 */
	pp[0][g2]+=pp[0][g1];
	pp[1][g2]+=pp[1][g1];
	/* Update size of g2 */
	gpsize[g2]+=gpsize[g1];
	/* Remove g1 from active list */
	k=gplist1[g1];
	k1=gplist[k]=gplist[--ngroups];
	gplist1[k1]=k;
}

/* Change group g1 to group g2 and link groups together */
static void change_group(int g1,int g2)
{
	int k,k1;
	
	/* Change group for members of g1 */
	k=first[g1];
	while(next[k]>=0) {
		group[k]=g2;
		k=next[k];
	}
	group[k]=g2;
	/* Add members of g1 to g2 list */
	next[k]=first[g2];
	first[g2]=first[g1];
	/* Update log probabilities for group g2 */
	pp[0][g2]+=pp[0][g1];
	pp[1][g2]+=pp[1][g1];
	/* Update size of g2 */
	gpsize[g2]+=gpsize[g1];
	/* Remove g1 from active list */
	k=gplist1[g1];
	k1=gplist[k]=gplist[--ngroups];
	gplist1[k1]=k;
}

void seg_init_freq(const struct Locus *loc,const struct loki *loki)
{
	int i,j,n_all;
	double z;
	struct Marker *mark;
	
	mark=(loc->type&ST_MARKER)?loki->markers->marker+loc->index:0;
	n_all=loc->n_alleles;
	for(j=0;j<loki->pedigree->n_genetic_groups;j++) {
		z=1.0/(double)n_all;;
		if(mark && mark->count_flag[j]) for(i=0;i<n_all;i++) loki->peel->seg_count[j][i]=mark->counts[j][i]+z;
		else for(i=0;i<n_all;i++) loki->peel->seg_count[j][i]=z;
	}
	if(loki->peel->aff_freq && mark) for(i=0;i<n_all;i++) loki->peel->aff_freq[i]=0.0;
}

void seg_sample_freq(const struct Locus *loc,const struct loki *loki)
{
	int i,i1,grp,n_all,n_genetic_groups;
	signed char *freq_set;
	double z,z1,z2,*count1,*freq;
	struct Marker *mark;
	
	mark=(loc->type&ST_MARKER)?loki->markers->marker+loc->index:0;	
	n_genetic_groups=loki->pedigree->n_genetic_groups;
	n_all=loc->n_alleles;
	for(grp=0;grp<n_genetic_groups;grp++) {
		z1=0.0;
		z=1.0;
		count1=loki->peel->seg_count[grp];
		freq=loc->freq[grp];
		if(!mark) {
			for(i=0;i<n_all;i++) z1+=count1[i];
			for(i=0;i<n_all-1;i++) {
				z1-=count1[i];
				freq[i]=z*genbet(count1[i],z1);
				z-=freq[i];
			}
			freq[i]=z; 
		} else {
			freq_set=mark->freq_set[grp];
			for(i1=-1,i=0;i<n_all;i++)	{
				if(freq_set[i]!=1) {
					z1+=count1[i];
					i1=i;
				} else z-=freq[i];
			}
			if(i1>=0) {
				for(i=0;i<i1;i++) if(freq_set[i]!=1) {
					z1-=count1[i];
					z2=genbet(count1[i],z1);
					freq[i]=z*z2;
					z-=freq[i];
				}
					freq[i]=z;
			}
		}
	}
}

void seg_update_aff_freq(const struct Locus *loc,const struct loki *loki)
{
	int i,n_all;
	double z;
	
	if(!loki->peel->aff_freq) return;
	n_all=loc->n_alleles;
	z=0.0;
	for(i=0;i<n_all;i++) {
		z+=loki->peel->aff_freq[i];
	}
	for(i=0;i<n_all;i++) {
		loc->aff_freq[i]=z>0.0?loki->peel->aff_freq[i]/z:0.0;
		loc->diff_freq[i]+=(loc->freq[0][i]<loc->aff_freq[i])?1.0:0.0;
	}
}

/* Calculate probability of segregation pattern conditional on 
* (a) genotype data (marker loci)
* (b) genotype data + trait data + sampled genotypes for individuals with trait data but no genotype data (candidate genes)
* (c) trait data + sampled genotypes for individuals with trait data (trait loci)
*/
double seg_pen(struct Locus *loc,int comp,int *err,int flag,const struct loki *loki)
{ 
	int i,i1,j,k,k1,k2,s,s1,a1,a2,g1,g2,o1=0,o2=0,o1a,o2a,o2b,*fflag,cs,cst;
	int k_sw[]={-1,1,1,-1,2,-1,3,-1,2,0,-1,-1,-1,-1,-1,-1};
	int k1_sw[]={-1,1,2,-1,1,-1,3,-1,2,0,-1,-1,-1,-1,-1,-1};
	int ngenes,*genesm,*genesp,*obs,*trans=0,*grp=0,*gt,*hap,si;
	int nn_all,n_all,idd,ids,**seg,sample_freq,locus_type,n_markers,n_genetic_groups;
	double p,z,z1,*ff_1=0,*ff_2=0,**ff1=0,*count1,ppen=0.0;
	struct Marker *mark;
	struct Id_Record *id_array;
	lk_ulong lump=0L,a;
	
	n_markers=loki->markers->n_markers;
	n_genetic_groups=loki->pedigree->n_genetic_groups;
	id_array=loki->pedigree->id_array;
	si=loki->params.si_mode;
	cs=loki->pedigree->comp_size[comp];
	cst=loki->pedigree->comp_start[comp];
	sample_freq=0;
	/* Marker locus */
	if(loc->type&ST_MARKER) {
		mark=loki->markers->marker+loc->index;
		hap=mark->haplo;
		n_all=mark->n_all1[comp];
		nn_all=loc->n_alleles;
		if(loc->eff) locus_type=1;
		else locus_type=0;
		s1=cpstart[loc->index][comp];
		s=cpstart[loc->index][comp+1]-s1;
		obs=obslist[loc->index]+s1;
		if(n_all>=2 && s) {
			/* Set up frequency arrays */
			trans=mark->allele_trans[comp];
			ff1=loc->freq;
			if(n_genetic_groups>1) grp=mark->group[comp];
			lump=0L;
			if(n_all==nn_all) for(i=0;i<nn_all;i++) seg_alls[i]=i;
			else {
				for(i=0;i<nn_all;i++)  seg_alls[i]=-1;
				for(i=0;i<nn_all;i++)  {
					k1=trans[i];
					if(k1>=0) seg_alls[k1]=i;
					else lump|=(1L<<i);
				}
			}
			for(k=0;k<n_genetic_groups;k++) {
				z=0.0;
				for(i=0;i<nn_all;i++)  {
					k1=seg_alls[i];
					if(k1>=0) seg_freq[k][k1]=log(ff1[k][i]);
					else z+=ff1[k][i];
				}
				if(n_all<nn_all) seg_freq[k][n_all-1]=log(z);
			}
			if(flag&4) for(k=0;k<n_genetic_groups;k++) {
				for(i=j=0;i<nn_all;i++) if(mark->freq_set[k][i]!=1) j++;
				if(j) sample_freq=4;
			}
				sample_freq|=(flag&8);
		}
		/* Trait locus */
	} else {
		lump=0L;
		hap=0;
		n_all=nn_all=loc->n_alleles;
		locus_type=2;
		s1=cpstart[n_markers][comp];
		s=cpstart[n_markers][comp+1]-s1;
		obs=obslist[n_markers]+s1;
		if(s>0) {
			/* Set up frequency arrays */
			if(n_genetic_groups>1) grp=loki->peel->tl_group[comp];
			ff1=loc->freq;
			for(k=0;k<n_genetic_groups;k++)
				for(i=0;i<n_all;i++) seg_freq[k][i]=log(ff1[k][i]);
			sample_freq=flag&12;
		}
	}
	if(n_all<2 || !s) {
		*err=0;
		return 0.0;
	}
	if(n_genetic_groups==1) ff_1=ff_2=seg_freq[0];
	*err=-1;
	genesm=loc->genes[X_MAT];
	genesp=loc->genes[X_PAT];
	gt=loc->gt;
	ngenes=loki->pedigree->comp_ngenes[comp];
	seg=loc->seg;
	j=cst;
	if(flag&8) {
		memset(fg_type,0,sizeof(int)*ngenes);
		j=cst;
		for(i=0;i<cs;i++,j++) {
			if(id_array[j].affected==2) {
				fg_type[genesm[j]-1]|=2;
				fg_type[genesp[j]-1]|=2;
			} else if(id_array[j].affected==1) {
				fg_type[genesm[j]-1]|=1;
				fg_type[genesp[j]-1]|=1;
			}
		}
	}
	/* Blank group array */
	(void)memset(group,-1,sizeof(int)*ngenes);
	/* Clear group counts */
	j=ngroups=0;
	if(!locus_type) {
		while((i=*(obs++))>=0) {
			k=hap[i];
			o1=lookup[k];
			o2=lookup1[k];
			a1=genesm[i]-1;
			a2=genesp[i]-1;
			if(n_genetic_groups>1) {
				ff_1=seg_freq[grp[a1]];
				ff_2=seg_freq[grp[a2]];
			}
			/* Is individual autozygous ? */
			if(a1!=a2) {
				/* Not autozygous, check if in existing group */
				g1=group[a1];
				g2=group[a2];
				if(g1>=0) {
					k=gpflag[g1];
					if(g2>=0) {
						/* Both genes in groups, same group ? */
						if(g1==g2) {
							/* Yes, group fixed ? */
							if(k--) {
								o1a=alleles[k][a1];
								o2a=alleles[k][a2];
								/* Yes, check match */
								if((o1a!=o1 || o2a!=o2) &&
									(o1a!=o2 || o2a!=o1)) return 4.0;
								pp[k][g1]+=ppen;
							} else {
								/* No, find which alternate matches */
								k1=3;
								o1a=alleles[0][a1];
								o2a=alleles[0][a2];
								if((o1a==o1 && o2a==o2)||(o1a==o2 && o2a==o1)) k1&=1;
								o1a=alleles[1][a1];
								o2a=alleles[1][a2];
								if((o1a==o1 && o2a==o2)||(o1a==o2 && o2a==o1)) k1&=2;
								if(k1==3) return 5.0; /* No match */
								if(!k1) {	
									pp[0][g1]+=ppen;
									pp[1][g1]+=ppen;
								} else {
									gpflag[g1]=k1--;
									pp[k1][g1]+=ppen;
								}
							}
						} else {
							/* Genes are in different groups. Groups fixed ? */
							k1=gpflag[g2];
							if(k) {
								o1a=alleles[k-1][a1];
								if(k1) {
									o2a=alleles[k1-1][a2];
									/* Both groups fixed, check allocations match */
									if((o1a!=o1 || o2a!=o2) &&
										(o1a!=o2 || o2a!=o1)) {
										return 6.0;
									}
								} else {
									/* Only first group fixed.  Check which possible
									* allocations for second group match */
									k1=0;
									if(o1a==o1) {
										if(alleles[0][a2]==o2) k1=1;
										else if(alleles[1][a2]==o2) k1=2;
									} else if(o1a==o2) {
										if(alleles[0][a2]==o1) k1=1;
										else if(alleles[1][a2]==o1) k1=2;
									}
									if(!k1) {
										/* No match */
										return 7.0;
									}
								}
							} else if(k1) {
								/* Only second group fixed.  Check which possible
								* allocations for first group match */
								o2a=alleles[k1-1][a2];
								k=0;
								if(o2a==o2) {
									if(alleles[0][a1]==o1) k=1;
									else if(alleles[1][a1]==o1) k=2;
								} else if(o2a==o1) {
									if(alleles[0][a1]==o2) k=1;
									else if(alleles[1][a1]==o2) k=2;
								}
								/* No match */
								if(!k) return 8.0;
							} else {
								/* No groups fixed.  See which allocations match up */
								k1=0;
								o1a=alleles[0][a1];
								o2a=alleles[0][a2];
								o2b=alleles[1][a2];
								if(o1a==o1) {
									if(o2a==o2) k1|=1;
									if(o2b==o2) k1|=2;
								}
								if(o1a==o2) {
									if(o2a==o1) k1|=1;
									if(o2b==o1) k1|=2;
								}
								o1a=alleles[1][a1];
								if(o1a==o1) {
									if(o2a==o2) k1|=4;
									if(o2b==o2) k1|=8;
								}
								if(o1a==o2) {
									if(o2a==o1) k1|=4;
									if(o2b==o1) k1|=8;
								}
								if(!k1) return 9.0; /* No match */
								/* Convert from k1 to appropriate gpflag entries */
								k=k_sw[k1];
								k1=k1_sw[k1];
								assert(k>=0 && k1>=0);
							}
							/* Are allocation flags the same ? */
							if(k!=k1) {
								/* No, so we must flip the allocations for all members of
								* one of the groups (the smallest) */
								if(gpsize[g1]<gpsize[g2]) {
									gpflag[g2]=k1;
									change_and_switch_group(g1,g2);
								} else {
									gpflag[g1]=k;
									change_and_switch_group(g2,g1);
								}
							} else {
								/* Yes, change group membership for smallest group */ 
								if(gpsize[g1]<gpsize[g2]) {
									if(k<3) {
										gpflag[g2]=k;
										change_group(g1,g2);
									} else {
										gpflag[g2]=0;
										change_and_swap_group(g1,g2);
									}
								} else {
									if(k<3) {
										gpflag[g1]=k;
										change_group(g2,g1);
									} else {
										gpflag[g1]=0;
										change_and_swap_group(g2,g1);
									}
								}
							}
							g1=group[a1];
							if((k=gpflag[g1])) pp[k-1][g1]+=ppen;
							else {
								pp[0][g1]+=ppen;
								pp[1][g1]+=ppen;
							} 
						}
					} else {
						/* Only a1 in group, a2 is new */
						/* Add a2 to group */
						next[a2]=first[g1];
						first[g1]=a2;
						group[a2]=g1;
						gpsize[g1]++;
						if(k--) {
							/* Group for a1 is fixed, check matches */
							o1a=alleles[k][a1];
							if(o1a==o1) alleles[k][a2]=o2;
							else if(o1a==o2) alleles[k][a2]=o1;
							else return 10.0; /* No match */
							pp[k][g1]+=ff_2[alleles[k][a2]-1]+ppen;
						} else {
							/* Group for a1 not fixed, check both allocations for matches */
							for(k1=k=0;k1<2;k1++) {
								o1a=alleles[k1][a1];
								if(o1a==o1) alleles[k1][a2]=o2;
								else if(o1a==o2) alleles[k1][a2]=o1;
								else k|=2-k1;
							}
							if(k==3) return 11.0; /* No match */
							gpflag[g1]=k;
							/* Add in contrib. to pp from gene a2 */
							if(k--) {
								pp[k][g1]+=ff_2[alleles[k][a2]-1]+ppen;
							} else {
								pp[0][g1]+=ff_2[alleles[0][a2]-1]+ppen;
								pp[1][g1]+=ff_2[alleles[1][a2]-1]+ppen;
							}
						}
					}
				} else if(g2>=0) {
					/* Only a2 in group, a1 is new */
					/* Add a1 to group */
					next[a1]=first[g2];
					first[g2]=a1;
					group[a1]=g2;
					gpsize[g2]++;
					k1=gpflag[g2];
					if(k1--) {
						/* Group for a2 is fixed, check matches */
						o1a=alleles[k1][a2];
						if(o1a==o1) alleles[k1][a1]=o2;
						else if(o1a==o2) alleles[k1][a1]=o1;
						else return 12.0; /* No match */
						pp[k1][g2]+=ff_1[alleles[k1][a1]-1]+ppen;
					} else {
						/* Group for a2 not fixed, check both allocations for matches */
						for(k1=k=0;k<2;k++) {
							o2a=alleles[k][a2];
							if(o2a==o1) alleles[k][a1]=o2;
							else if(o2a==o2) alleles[k][a1]=o1;
							else k1|=2-k;
						}
						if(k1==3) return 13.0; /* No match */
						gpflag[g2]=k1;
						/* Add in contrib. to pp from gene a1 */
						if(k1--) {
							pp[k1][g2]+=ff_1[alleles[k1][a1]-1]+ppen;
						} else {
							pp[0][g2]+=ff_1[alleles[0][a1]-1]+ppen;
							pp[1][g2]+=ff_1[alleles[1][a1]-1]+ppen;
						}
					}
					/* A first for both genes */
				} else { 
					if(o1!=o2) {
						/* Heterozygote */
						alleles[0][a1]=alleles[1][a2]=o1;
						alleles[1][a1]=alleles[0][a2]=o2;
						gpflag[j]=0;
						pp[0][j]=ff_1[o1-1]+ff_2[o2-1]+ppen;
						pp[1][j]=ff_2[o1-1]+ff_2[o2-1]+ppen;
					} else {
						/* Homozygote */
						alleles[0][a1]=alleles[0][a2]=o1;
						gpflag[j]=1;
						pp[0][j]=ff_1[o1-1]+ff_2[o1-1]+ppen;
					}
					/* Create new group */
					first[j]=a1;
					next[a1]=a2;
					next[a2]=-1;
					gpsize[j]=2;
					/* Add to active list */
					gplist1[j]=ngroups;
					group[a1]=group[a2]=j;
					gplist[ngroups++]=j++;
				}
			} else {
				/* Yes, Check homozygosity */
				if(o1!=o2) {
					return 1.0;
				}
				g1=group[a1];
				/* Already in group ? */
				if(g1>=0) {
					/* Yes, is group is fixed ? */
					k=gpflag[g1];
					if(k--) {
						/* Yes, Check allele matches ? */
						if(alleles[k][a1]!=o1) return 2.0;
						pp[k][g1]+=ppen;
					} else {
						/* No, find which alternate matches current allele */
						k=0;
						if(alleles[0][a1]==o1) k=1;
						else if(alleles[1][a1]==o1) k=2;
						if(!k) return 3.0; /* No match */
						/* Fix group appropriately */
						gpflag[g1]=k;
						pp[k-1][g1]+=ppen;
					}
				} else {
					/* No, create new group */
					first[j]=a1;
					gpflag[j]=gpsize[j]=1;
					next[a1]=-1;
					/* Add in contrib. to pp from this gene */
					pp[0][j]=ff_1[o1-1]+ppen;
					alleles[0][a1]=o1;
					/* Add group to active list */
					gplist[ngroups]=j;
					gplist1[j]=ngroups++;
					/* Set group membership for a1 and increment group counter */
					group[a1]=j++;
				}
			}
		}
		/* All allocations have been made.  Go through active group list
		* and accumulate the probability */
		p=0.0;
		/* First check if we are sampling */
		if(!flag) {
			/* No, just calculate the probability */
			for(j=0;j<ngroups;j++) {
				k=gplist[j];
				k1=gpflag[k];
				if(k1--) {
					/* Group is fixed (only 1 possible allocation) */
					p+=pp[k1][k];
					assert(!isnan(p));
				} else {
					/* 2 possible allocations for group */
					p+=addlog(pp[0][k],pp[1][k]);
					assert(!isnan(p));
				}
			}
		} else {
			/* Yes, sample founder alleles, pass down pedigree and update segregation pattern */
			for(j=0;j<ngroups;j++) {
				k=gplist[j];
				k1=gpflag[k];
				/* Group is fixed (only 1 possible allocation) */
				if(k1--) { 
					p+=pp[k1][k];
					assert(!isnan(p));
					/* 2 possible allocations for group */
				} else {
					z=addlog(pp[0][k],pp[1][k]);
					p+=z;
					assert(!isnan(p));
					/* Sample one of the possible allocations */
					z=exp(pp[0][k]-z);
					gpflag[k]=(ranf()<z)?1:2;
				}
			}
			/* Sample non-observed founder genes */
			k=gplist[0];
			k1=gpflag[k]-1;
			for(a1=0;a1<ngenes;a1++) {
				if(group[a1]<0) {
					if(n_genetic_groups>1) ff_1=seg_freq[grp[a1]];
					if(ff_1[0]<0.0) for(i=0;i<n_all;i++) ff_1[i]=exp(ff_1[i]);
					z=ranf();
					z1=0.0;
					for(i=0;i<n_all;i++) {
						z1+=ff_1[i];
						if(z<=z1) break;
					}
					group[a1]=k;
					alleles[k1][a1]=i+1;
				}
			}
			/* Pass founder alleles down through pedigree */
			i=cst;
			seg=loc->seg;
			fflag=loc->founder_flag;
			for(i1=0;i1<cs;i1++,i++) {
				a1=genesm[i]-1;
				a2=genesp[i]-1;
				g1=group[a1];
				k=gpflag[g1]-1;
				id_array[i].allele[X_MAT]=alleles[k][a1];
				g1=group[a2];
				k1=gpflag[g1]-1;
				id_array[i].allele[X_PAT]=alleles[k1][a2];
				if((flag&2) && !si) {
					/* Remove ambiguous segregations */
					if(fflag[i]) seg[X_MAT][i]=seg[X_PAT][i]=-1;
					else {
						idd=id_array[i].dam;
						k=id_array[i].allele[X_MAT];
						assert(k==id_array[idd-1].allele[seg[X_MAT][i]]);
						if(id_array[idd-1].allele[X_MAT]==id_array[idd-1].allele[X_PAT]) seg[X_MAT][i]=-2;
						ids=id_array[i].sire;
						k=id_array[i].allele[X_PAT];
						assert(k==id_array[ids-1].allele[seg[X_MAT][i]]);
						if(id_array[ids-1].allele[X_MAT]==id_array[ids-1].allele[X_PAT]) seg[X_PAT][i]=-2;
					}
				}
				/* Update stored genotypes */
				if(loc->pruned_flag[i]) k=0;
				else {
					k=id_array[i].allele[X_MAT];
					k1=id_array[i].allele[X_PAT];
					if(k>k1) k=k*(k-1)/2+k1;
					else k=k1*(k1-1)/2+k;
				}
				loc->gt[i]=k;
				/* Note that we don't have to update residuals because the genotypes of individuals with trait data
					* are not changed */
			}
			/* Get allele counts for frequency update */
			if(sample_freq&4) {
				if(n_genetic_groups==1) {
					count1=loki->peel->seg_count[0];
					ff_1=ff1[0];
					for(i=0;i<ngenes;i++) {
						g1=group[i];
						k=gpflag[g1]-1;
						k1=alleles[k][i]-1;
						k2=trans[k1];
						if(k2<0) {
							k2=0;
							a=lump;
							z=0.0;
							while(a) {
								if(a&1) z+=ff_1[k2];
								a>>=1;
								k2++;
							}
							z=1.0/z;
							k2=0;
							a=lump;
							while(a) {
								if(a&1) count1[k2]+=z*ff_1[k2];
								a>>=1;
								k2++;
							}
						} else count1[k2]+=1.0;
					}
				} else {
					for(i=0;i<ngenes;i++) {
						count1=loki->peel->seg_count[grp[i]];
						ff_1=ff1[grp[i]];
						g1=group[i];
						k=gpflag[g1]-1;
						k1=alleles[k][i]-1;
						k2=trans[k1];
						if(k2<0) {
							k2=0;
							a=lump;
							z=0.0;
							while(a) {
								if(a&1) z+=ff_1[k2];
								a>>=1;
								k2++;
							}
							z=1.0/z;
							k2=0;
							a=lump;
							while(a) {
								if(a&1) count1[k2]+=z*ff_1[k2];
								a>>=1;
								k2++;
							}
						} else count1[k2]+=1.0;
					}
				}
			}
			if(sample_freq&8) {
				for(i=0;i<ngenes;i++) {
					g1=group[i];
					k=gpflag[g1]-1;
					k1=alleles[k][i]-1;
					if(fg_type[i]&2) loki->peel->aff_freq[k1]+=1.0;
				}
			}
		}
	} else {
		while((i=*(obs++))>=0) {
			/* Candidate gene */
			if(locus_type==1) { 
				if(!hap[i]) {
					k=gt[i]-1;
					assert(k>=0);
					o1=lookup_a[k];
					o2=lookup1_a[k];
					assert(o1>=1 && o2>=1 && o1<=nn_all && o2<=nn_all);
				} else {
					k=hap[i];
					o1=lookup[k];
					o2=lookup1[k];
				}
				ppen=q_penetrance(i,k,loc,loki);
			} else {
				/* Trait locus */
				k=gt[i];
				assert(k>0);
				o1=lookup_a[k];
				o2=lookup1_a[k];
				assert(o1>=1 && o2>=1 && o1<=nn_all && o2<=nn_all);
				ppen=q_penetrance(i,k,loc,loki);
			}
			a1=genesm[i]-1;
			a2=genesp[i]-1;
			if(n_genetic_groups>1) {
				ff_1=seg_freq[grp[a1]];
				ff_2=seg_freq[grp[a2]];
			}
			/* Is individual autozygous ? */
			if(a1!=a2) {
				/* Not autozygous, check if in existing group */
				g1=group[a1];
				g2=group[a2];
				if(g1>=0) {
					k=gpflag[g1];
					if(g2>=0) {
						/* Both genes in groups, same group ? */
						if(g1==g2) {
							/* Yes, group fixed ? */
							if(k--) {
								o1a=alleles[k][a1];
								o2a=alleles[k][a2];
								/* Yes, check match */
								if((o1a!=o1 || o2a!=o2) &&
									(o1a!=o2 || o2a!=o1)) {
									return 4.0;
								}
								pp[k][g1]+=ppen;
							} else {
								/* No, find which alternate matches */
								for(k1=k=0;k<2;k++) {
									o1a=alleles[k][a1];
									o2a=alleles[k][a2];
									if((o1a!=o1 || o2a!=o2) &&
										(o1a!=o2 || o2a!=o1)) k1|=2-k;
								}
								if(k1==3) return 5.0; /* No match */
								gpflag[g1]=k1; /* Fix group */
								if(k1--) pp[k1][g1]+=ppen;
								else for(k1=0;k1<2;k1++) pp[k1][g1]+=ppen;
							}
						} else {
							/* Genes are in different groups. Groups fixed ? */
							k1=gpflag[g2];
							if(k) {
								o1a=alleles[k-1][a1];
								if(k1) {
									o2a=alleles[k1-1][a2];
									/* Both groups fixed, check allocations match */
									if((o1a!=o1 || o2a!=o2) &&
										(o1a!=o2 || o2a!=o1)) {
										return 6.0;
									}
								} else {
									/* Only first group fixed.  Check which possible
									* allocations for second group match */
									for(k1=0;k1<2;k1++) {
										o2a=alleles[k1][a2];
										if((o1a==o1 && o2a==o2) || 
											(o1a==o2 && o2a==o1)) break;
									}
									if(k1++==2) {
										return 7.0; /* No match */
									}
								}
							} else if(k1) {
								/* Only second group fixed.  Check which possible
								* allocations for first group match */
								o2a=alleles[k1-1][a2];
								for(k=0;k<2;k++) {
									o1a=alleles[k][a1];
									if((o1a==o1 && o2a==o2) || 
										(o1a==o2 && o2a==o1)) break;
								}
								if(k++==2) {
									return 8.0; /* No match */
								}
							} else {
								/* No groups fixed.  See which allocations match up */
								for(k1=k=0;k<4;k++) {
									o1a=alleles[(k&2)>>1][a1];
									o2a=alleles[k&1][a2];
									if((o1a==o1 && o2a==o2) ||
										(o1a==o2 && o2a==o1)) k1|=1<<k;
								}
								if(!k1) return 9.0; /* No match */
								/* Convert from k to appropriate gpflag entries */
								k=k_sw[k1];
								k1=k1_sw[k1];
								assert(k>=0 && k1>=0);
							}
							/* Are allocation flags the same ? */
							if(k!=k1) {
								/* No, so we must flip the allocations for all members of
								* one of the groups (the smallest) */
								if(gpsize[g1]<gpsize[g2]) {
									gpflag[g2]=k1;
									change_and_switch_group(g1,g2);
								} else {
									gpflag[g1]=k;
									change_and_switch_group(g2,g1);
								}
							} else {
								/* Yes, change group membership for smallest group */ 
								if(gpsize[g1]<gpsize[g2]) {
									if(k<3) {
										gpflag[g2]=k;
										change_group(g1,g2);
									} else {
										gpflag[g2]=0;
										change_and_swap_group(g1,g2);
									}
								} else {
									if(k<3) {
										gpflag[g1]=k;
										change_group(g2,g1);
									} else {
										gpflag[g1]=0;
										change_and_swap_group(g2,g1);
									}
								}
							}
							g1=group[a1];
							k=gpflag[g1];
							if(k) pp[k-1][g1]+=ppen;
							else for(k=0;k<2;k++) pp[k][g1]+=ppen;
						}
					} else {
						/* Only a1 in group, a2 is new */
						if(k) {
							/* Group for a1 is fixed, check matches */
							o1a=alleles[k-1][a1];
							if(o1a==o1) alleles[k-1][a2]=o2;
							else if(o1a==o2) alleles[k-1][a2]=o1;
							else return 10.0; /* No match */
						} else {
							/* Group for a1 not fixed, check both allocations for matches */
							for(k1=k=0;k1<2;k1++) {
								o1a=alleles[k1][a1];
								if(o1a==o1) alleles[k1][a2]=o2;
								else if(o1a==o2) alleles[k1][a2]=o1;
								else k|=2-k1;
							}
							if(k==3) return 11.0; /* No match */
							gpflag[g1]=k;
						}
						/* Add a2 to group */
						next[a2]=first[g1];
						first[g1]=a2;
						group[a2]=g1;
						gpsize[g1]++;
						/* Add in contrib. to pp from gene a2 */
						if(k--) pp[k][g1]+=ff_2[alleles[k][a2]-1]+ppen;
						else for(k=0;k<2;k++) pp[k][g1]+=ff_2[alleles[k][a2]-1]+ppen;
					}
				} else if(g2>=0) {
					/* Only a2 in group, a1 is new */
					k1=gpflag[g2];
					if(k1) {
						/* Group for a2 is fixed, check matches */
						o1a=alleles[k1-1][a2];
						if(o1a==o1) alleles[k1-1][a1]=o2;
						else if(o1a==o2) alleles[k1-1][a1]=o1;
						else return 12.0; /* No match */
					} else {
						/* Group for a2 not fixed, check both allocations for matches */
						for(k1=k=0;k<2;k++) {
							o2a=alleles[k][a2];
							if(o2a==o1) alleles[k][a1]=o2;
							else if(o2a==o2) alleles[k][a1]=o1;
							else k1|=2-k;
						}
						if(k1==3) return 13.0; /* No match */
						gpflag[g2]=k1;
					}
					/* Add a1 to group */
					next[a1]=first[g2];
					first[g2]=a1;
					group[a1]=g2;
					gpsize[g2]++;
					/* Add in contrib. to pp from gene a1 */
					if(k1--) pp[k1][g2]+=ff_1[alleles[k1][a1]-1]+ppen; 
					else for(k1=0;k1<2;k1++) pp[k1][g2]+=ff_1[alleles[k1][a1]-1]+ppen;
					/* A first for both genes */
				} else { 
					if(o1!=o2) {
						/* Heterozygote */
						alleles[0][a1]=alleles[1][a2]=o1;
						alleles[1][a1]=alleles[0][a2]=o2;
						gpflag[j]=0;
						pp[0][j]=ff_1[o1-1]+ff_2[o2-1]+ppen;
						pp[1][j]=ff_2[o1-1]+ff_1[o2-1]+ppen;
					} else {
						/* Homozygote */
						alleles[0][a1]=alleles[0][a2]=o1;
						gpflag[j]=1;
						pp[0][j]=ff_1[o1-1]+ff_2[o1-1]+ppen;
					}
					/* Create new group */
					first[j]=a1;
					next[a1]=a2;
					next[a2]=-1;
					gpsize[j]=2;
					/* Add to active list */
					gplist[ngroups]=j;
					gplist1[j]=ngroups++;
					group[a1]=group[a2]=j++;
				}
			} else {
				/* Yes, Check homozygosity */
				if(o1!=o2) {
					return 1.0;
				}
				g1=group[a1];
				/* Already in group ? */
				if(g1>=0) {
					/* Yes, is group is fixed ? */
					k=gpflag[g1];
					if(k--) {
						/* Yes, Check allele matches ? */
						if(alleles[k][a1]!=o1) return 2.0;
						pp[k][g1]+=ppen;
					} else {
						/* No, find which alternate matches current allele */
						for(k=0;k<2;k++) {
							o1a=alleles[k][a1];
							if(o1a==o1) break;
						}
						if(k==2) return 3.0; /* No match */
						/* Fix group appropriately */
						gpflag[g1]=k+1;
						pp[k][g1]+=ppen;
					}
				} else {
					/* No, create new group */
					first[j]=a1;
					gpflag[j]=gpsize[j]=1;
					next[a1]=-1;
					/* Add in contrib. to pp from this gene */
					pp[0][j]=ff_1[o1-1]+ppen;
					alleles[0][a1]=o1;
					/* Add group to active list */
					gplist[ngroups]=j;
					gplist1[j]=ngroups++;
					/* Set group membership for a1 and increment group counter */
					group[a1]=j++;
				}
			}
		}
		/* All allocations have been made.  Go through active group list
		* and accumulate the probability */
		p=0.0;
		/* First check if we are sampling */
		if(!flag) {
			/* No, just calculate the probability */
			for(j=0;j<ngroups;j++) {
				k=gplist[j];
				k1=gpflag[k];
				if(k1--) {
					/* Group is fixed (only 1 possible allocation) */
					p+=pp[k1][k];
					assert(!isnan(p));
				} else {
					/* 2 possible allocations for group */
					p+=addlog(pp[0][k],pp[1][k]);
					assert(!isnan(p));
				}
			}
		} else {
			/* Yes, sample founder alleles, pass down pedigree and update segregation pattern */
			for(j=0;j<ngroups;j++) {
				k=gplist[j];
				k1=gpflag[k];
				if(k1--) {
					/* Group is fixed (only 1 possible allocation) */
					p+=pp[k1][k];
					assert(!isnan(p));
				} else {
					/* 2 possible allocations for group */
					z=addlog(pp[0][k],pp[1][k]);
					p+=z;
					assert(!isnan(p));
					/* Sample one of the possible allocations */
					z=exp(pp[0][k]-z);
					gpflag[k]=(ranf()<z)?1:2;
				}
			}
			/* Sample non-observed founder genes */
			k=gplist[0];
			k1=gpflag[k]-1;
			for(a1=0;a1<ngenes;a1++) {
				if(group[a1]<0) {
					if(n_genetic_groups>1) ff_1=seg_freq[grp[a1]];
					if(ff_1[0]<0.0) for(i=0;i<n_all;i++) ff_1[i]=exp(ff_1[i]);
					z=ranf();
					z1=0.0;
					for(i=0;i<n_all;i++) {
						z1+=ff_1[i];
						if(z<=z1) break;
					}
					group[a1]=k;
					alleles[k1][a1]=i+1;
				}
			}
			/* Pass founder alleles down through pedigree */
			i=cst;
			seg=loc->seg;
			fflag=loc->founder_flag;
			for(i1=0;i1<cs;i1++,i++) {
				a1=genesm[i]-1;
				a2=genesp[i]-1;
				g1=group[a1];
				k=gpflag[g1]-1;
				id_array[i].allele[X_MAT]=alleles[k][a1];
				g1=group[a2];
				k1=gpflag[g1]-1;
				id_array[i].allele[X_PAT]=alleles[k1][a2];
				if((flag&2) && !si) {
					/* Remove ambiguous segregations */
					if(fflag[i]) seg[X_MAT][i]=seg[X_PAT][i]=-1;
					else {
						idd=id_array[i].dam;
						k=id_array[i].allele[X_MAT];
						assert(k==id_array[idd-1].allele[seg[X_MAT][i]]);
						if(id_array[idd-1].allele[X_MAT]==id_array[idd-1].allele[X_PAT]) seg[X_MAT][i]=-2;
						ids=id_array[i].sire;
						k=id_array[i].allele[X_PAT];
						assert(k==id_array[ids-1].allele[seg[X_MAT][i]]);
						if(id_array[ids-1].allele[X_MAT]==id_array[ids-1].allele[X_PAT]) seg[X_PAT][i]=-2;
					}
				}
				/* Update stored genotypes */
				if(loc->pruned_flag[i]) k=0;
				else {
					k=id_array[i].allele[X_MAT];
					k1=id_array[i].allele[X_PAT];
					if(k>k1) k=k*(k-1)/2+k1;
					else k=k1*(k1-1)/2+k;
				}
				loc->gt[i]=k;
				/* Note that we don't have to update residuals becuase the genotypes of individuals with trait data
					* are not changed */
			}
			/* Get allele counts for frequency update */
			if(sample_freq) {
				if(n_genetic_groups==1) {
					count1=loki->peel->seg_count[0];
					ff_1=ff1[0];
					for(i=0;i<ngenes;i++) {
						g1=group[i];
						k=gpflag[g1]-1;
						k1=alleles[k][i]-1;
						k2=trans[k1];
						if(k2<0) {
							k2=0;
							a=lump;
							z=0.0;
							while(a) {
								if(a&1) z+=ff_1[k2];
								a>>=1;
								k2++;
							}
							z=1.0/z;
							k2=0;
							a=lump;
							while(a) {
								if(a&1) count1[k2]+=z*ff_1[k2];
								a>>=1;
								k2++;
							}
						} else count1[k2]+=1.0;
					}
				} else {
					for(i=0;i<ngenes;i++) {
						count1=loki->peel->seg_count[grp[i]];
						ff_1=ff1[grp[i]];
						g1=group[i];
						k=gpflag[g1]-1;
						k1=alleles[k][i]-1;
						k2=trans[k1];
						if(k2<0) {
							k2=0;
							a=lump;
							z=0.0;
							while(a) {
								if(a&1) z+=ff_1[k2];
								a>>=1;
								k2++;
							}
							z=1.0/z;
							k2=0;
							a=lump;
							while(a) {
								if(a&1) count1[k2]+=z*ff_1[k2];
								a>>=1;
								k2++;
							}
						} else count1[k2]+=1.0;
					}
				}
			}
		}
	}
	*err=0;
	assert(!isnan(p));
	return p;
}
