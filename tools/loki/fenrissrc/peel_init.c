/****************************************************************************
 *                                                                          *
 *     Loki - Programs for genetic analysis of complex traits using MCMC    *
 *                                                                          *
 *                    Simon Heath - CNG, France                             *
 *                                                                          *
 *                         August 2002                                      *
 *                                                                          *
 * peel_init.c:                                                             *
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

static struct family *get_families(int comp,int *nfam)
{
	int i,j,k,kid,nf=0,cs,id,idd,ids,nkids,*kids;
	struct family *fam;
	
	cs=comp_size[comp];
	id=comp_start[comp];
	/* Find nuclear families */
	for(i=0;i<cs;i++) {
		id_array[id+i].flag=0;
		id_array[id+i].rfp=-1;
	}
	for(nf=nkids=i=0;i<cs;i++) if(!id_array[id+i].flag) {
		ids=id_array[id+i].sire;
		idd=id_array[id+i].dam;
		if(ids||idd) {
			j=ids?ids-1:idd-1;
			for(k=0;k<id_array[j].nkids;k++) {
				kid=id_array[j].kids[k];
				if(id_array[kid].sire==ids && id_array[kid].dam==idd) {
					id_array[kid].flag=1;
					nkids++;
				}
			}
			nf++;
		}
	}
	if(!(kids=malloc(sizeof(int)*nkids))) ABT_FUNC(MMsg);
	if(!(fam=malloc(sizeof(struct family)*nf))) ABT_FUNC(MMsg);
	for(nf=i=0;i<cs;i++) id_array[id+i].flag=0;
	for(i=0;i<cs;i++) if(id_array[id+i].flag!=2) {
		ids=id_array[id+i].sire;
		idd=id_array[id+i].dam;
		if(ids||idd) {
			j=ids?ids-1:idd-1;
			fam[nf].nkids=0;
			fam[nf].kids=kids;
			fam[nf].ids=ids;
			fam[nf].idd=idd;
			for(k=0;k<id_array[j].nkids;k++) {
				kid=id_array[j].kids[k];
				if(id_array[kid].sire==ids && id_array[kid].dam==idd) {
					id_array[kid].flag=2;
					fam[nf].nkids++;
					*kids++=kid;
				}
			}
			nf++;
		}
	}
	*nfam=nf;
	return fam;
}

static void calc_order(int comp,struct family *fam,int nfam,int *order)
{
	int i,j,nk,cs,id;
	
	cs=comp_size[comp];
	id=comp_start[comp];
	for(i=0;i<cs;i++) order[i]=0;
	for(i=0;i<nfam;i++) {
		order[fam[i].ids-1-id]++;
		order[fam[i].idd-1-id]++;
		nk=fam[i].nkids;
		for(j=0;j<nk;j++) order[fam[i].kids[j]-id]++;
	}
}

static int add_to_involved(const int i,int n,int *inv)
{
	int j;
	
	for(j=0;j<n;j++) if(i==inv[j]) break;
	if(j==n) inv[n++]=i;
	return n;
}

static int check_involved(const int i,int n,int *inv)
{
	int j;
	
	for(j=0;j<n;j++) if(i==inv[j]) break;
	return (j==n)?0:1;
}

static int get_involved(const int i,int *inv,int n_inv,struct FR_Func *rfunc,int n_rf)
{
	int j,j1,k,kid,nk,par,flag,sex;
	
	flag=n_inv;
	if(flag) {
		if(!check_involved(i,n_inv,inv)) return -1;
	} else inv[n_inv++]=i;
	if(id_array[i].flag&2) { /* Put in parents */
		par=id_array[i].sire;
		if(flag) {
			if(!check_involved(par-1,n_inv,inv)) return -1;
		} else inv[n_inv++]=par-1;
		par=id_array[i].dam;
		if(flag) {
			if(!check_involved(par-1,n_inv,inv)) return -1;
		} else inv[n_inv++]=par-1;
		id_array[i].flag|=4;
	}
	/* Check kids */
	nk=id_array[i].nkids;
	if(nk) {
		sex=id_array[i].sex;
		for(k=0;k<nk;k++) {
			kid=id_array[i].kids[k];
			if(id_array[kid].flag&2) {
				par=(sex==1)?id_array[kid].dam:id_array[kid].sire;
				if(flag) {
					if(!check_involved(par-1,n_inv,inv)) return -1;
				} else n_inv=add_to_involved(par-1,n_inv,inv);
				if(flag) {
					if(!check_involved(kid,n_inv,inv)) return -1;
				} else n_inv=add_to_involved(kid,n_inv,inv);
				id_array[kid].flag|=4;
			}
		}
	}
	/* Check R-Functions */
	for(k=0;k<n_rf;k++) if(!(rfunc[k].flag&6)) {
		rfunc[k].flag=0;
		for(j=0;j<rfunc[k].n_ind;j++) if(rfunc[k].id_list[j]==i) break;
		if(j<rfunc[k].n_ind) {
			for(j=0;j<rfunc[k].n_ind;j++) {
				j1=rfunc[k].id_list[j];
				if(flag) {
					if(!check_involved(j1,n_inv,inv)) return -1;
				} else n_inv=add_to_involved(j1,n_inv,inv);
			}
			rfunc[k].flag=1;
		}
	}
	for(k=0;k<n_rf;k++) if(rfunc[k].flag==1) rfunc[k].flag=2;
	for(k=0;k<n_inv;k++) {
		j=inv[k];
		if(id_array[j].flag&0x4) id_array[j].flag^=0x16;
	}
	return n_inv;
}

static int get_complex_peel_sequence(struct Peelseq_Head *pp,int comp,struct family *fam,int *famlist,int nf)
{
	int i,x,j,k,k1,kid,cs,id,n_rf,r_func_size,n_ind,ids,idd,n_rf1;
	int *mat,ptr,n_nfd,*perm,*trans,par,n_inv,*inv,n_peel,sex,nk,n_ops;
	struct FR_Func *rfuncs;
	struct family *ff;
	struct Complex_Element *element;
	
	/* Count no. of unpeeled individuals */
	cs=comp_size[comp];
	id=comp_start[comp];
	for(i=0;i<cs;i++) id_array[id+i].flag=0;
	for(n_ind=i=0;i<nf;i++) {
		ff=fam+famlist[i];
		ids=ff->ids;
		idd=ff->idd;
		if(!id_array[ids-1].flag) {
			id_array[ids-1].flag=1;
			n_ind++;
		}
		if(!id_array[idd-1].flag) {
			id_array[idd-1].flag=1;
			n_ind++;
		}
		for(k=0;k<ff->nkids;k++) {
			kid=ff->kids[k];
			if(!id_array[kid].flag) {
				id_array[kid].flag=1;
				n_ind++;
			}
		}
	}
	if(!n_ind) ABT_FUNC("No unpeeled individuals left\n");
	/* Assemble matrix of connections between individuals */
	if(!(perm=malloc(sizeof(int)*(n_ind+cs)))) ABT_FUNC(MMsg);
	trans=perm+n_ind;
	for(n_nfd=i=k=0;i<cs;i++) if(id_array[id+i].flag) {
		par=id_array[id+i].sire;
		if(par && id_array[par-1].flag) n_nfd++;
		perm[k]=id+i;
		trans[i]=k++;
	}
	/* Allocate matrix */
	k=n_ind+1+n_nfd*2+nf;
	if(!(mat=malloc(sizeof(int)*k))) ABT_FUNC(MMsg);
	ptr=n_ind+1;
	for(x=0;x<n_ind;x++) {
		/* Individual to deal with */
		i=perm[x];
		mat[x]=ptr;
		/* Put in parents (if not already peeled) */
		par=id_array[i].sire;
		if(par && id_array[par-1].flag) {
			mat[ptr++]=trans[par-1-id];
			par=id_array[i].dam;
			mat[ptr++]=trans[par-1-id];
			id_array[i].flag|=2;
		}
		/* Put in spouse(s) (if spouse(s) < self) */
		nk=id_array[i].nkids;
		if(nk) {
			sex=id_array[i].sex;
			if(!sex || sex>2) ABT_FUNC("Illegal sex!\n");
			for(k=0;k<nk;k++) {
				kid=id_array[i].kids[k];
				par=(sex==1)?id_array[kid].dam:id_array[kid].sire;
				if(par-1<i) {
					if(!(id_array[par-1].flag&64)) {
						mat[ptr++]=trans[par-1-id];
						id_array[par-1].flag|=64;
					}
				}
			}
		}
	}
	mat[x]=ptr;
	/* Get minimum degree ordering (into trans) */
	min_deg(n_ind,mat,trans,0);
	i=trans[1];
	trans[1]=trans[4];
	trans[4]=i;
	free(mat);
	/* Allocate space for some R-Functions (we'll increase this later if necessary) */
	r_func_size=nf;
	if(!(rfuncs=malloc(sizeof(struct FR_Func)*r_func_size))) ABT_FUNC(MMsg);
	n_rf=n_ops=0;
	if(!(inv=malloc(sizeof(int)*n_ind))) ABT_FUNC(MMsg);
	/* Assemble peel ops */
	for(x=0;x<n_ind;x++) {
		/* Get pivot */
		i=perm[trans[x]];
		/* Get nodes connected to pivot */
		n_inv=get_involved(i,inv,0,rfuncs,n_rf);
		/* Check if subsequent pivots can be peeled at the same time (i.e., without adding extra people) */
		for(k=x+1;k<n_ind;k++) if(get_involved(perm[trans[k]],inv,n_inv,rfuncs,n_rf)<0) {
			for(k1=0;k1<n_inv;k1++) {
				j=inv[k1];
				if(id_array[j].flag&0x4) id_array[j].flag&=~0x4;
			}
			break;
		}
		n_peel=k-x;
		x=k-1;
		if(!(element=malloc(sizeof(struct Complex_Element)))) ABT_FUNC(MMsg);
		pp->type=PEEL_COMPLEX;
		pp->ptr.complex=element;
		element->n_peel=n_peel;
		element->n_involved=n_inv;
		for(n_rf1=k=0;k<n_rf;k++) if(rfuncs[k].flag==2) n_rf1++;
		element->n_rfuncs=n_rf1;
		if(!(element->involved=malloc(sizeof(int)*(n_inv*2+n_rf1)))) ABT_FUNC(MMsg);
		element->flags=element->involved+n_inv;
		element->index=element->flags+n_inv;
		for(k1=k=0;k<n_inv;k++) {
			j=inv[k];
			element->flags[k]=0;
			if(id_array[j].flag&0x10) {
				element->flags[k]=1;
				id_array[j].flag&=~0x10;
			}
			element->involved[k]=j;
			if(id_array[j].rfp>=0) 	id_array[j].rfp=-1;
		}
		for(k=0;k<n_rf;k++) if(rfuncs[k].flag==2) {
			element->index[k1++]=k;
			rfuncs[k].flag|=4;
		}
		if(n_peel==n_inv) element->out_index=-1;
		else {
			n_ops++;
			element->out_index=n_rf;
			if(n_rf>=r_func_size) {
				r_func_size<<=1;
				if(!(rfuncs=realloc(rfuncs,sizeof(struct FR_Func)*r_func_size))) ABT_FUNC(MMsg);
			}
			if(!(rfuncs[n_rf].id_list=malloc(sizeof(int)*(n_inv-n_peel)))) ABT_FUNC(MMsg);
			rfuncs[n_rf].n_ind=n_inv-n_peel;
			rfuncs[n_rf].flag=0;
			for(k=n_peel;k<n_inv;k++) rfuncs[n_rf].id_list[k-n_peel]=inv[k];
			element->out_index=n_rf++;
		}
		pp= &element->next;
		pp->type=0;
	}
	for(i=0;i<n_rf;i++) free(rfuncs[i].id_list);
	free(rfuncs);
	free(inv);
	free(perm);
	return n_ops;
}

static void get_peel_sequence(struct Peelseq_Head *pp,int comp,struct family *fam,int nfam,int *n_ops)
{
	int i,j,k,id,*order,*famlist,ids,idd,nk,kid,pivot,nf,n_rf;
	struct family *ff;
	struct Fenris_Simple_Element *element;
	
	if(!(order=malloc(sizeof(int)*(nfam+comp_size[comp])))) ABT_FUNC(MMsg);
	famlist=order+comp_size[comp];
	/* Calculate order of individuals */
	calc_order(comp,fam,nfam,order);
	/* Look for simple peeling sequence */
	n_rf=0;
	for(i=0;i<nfam;i++) famlist[i]=i;
	id=comp_start[comp];
	nf=nfam;
	for(i=0;i<nf;i++) {
		ff=fam+famlist[i];
		ids=ff->ids;
		idd=ff->idd;
		pivot=0;
		if(order[ids-1-id]>1) pivot=ids;
		if(order[idd-1-id]>1) pivot=pivot?-1:idd;
		if(pivot>=0) {
			nk=ff->nkids;
			for(j=0;j<nk;j++) {
				kid=ff->kids[j];
				if(order[kid-id]>1) pivot=pivot?-1:kid+1;
				if(pivot<0) break;
			}
		}
		/* Does family have a single pivot (or no pivot)? */
		if(pivot>=0) {
			famlist[i]=famlist[--nf];
			if(!(element=malloc(sizeof(struct Fenris_Simple_Element)))) ABT_FUNC(MMsg);
			pp->type=FENRIS_PEEL_SIMPLE;
			pp->ptr.fsimple=element;
			element->sire=ids;
			element->dam=idd;
			element->pivot=pivot;
			nk=ff->nkids;
			/* Store offspring list in element.
			 * First we put terminal individuals (apart from pivot).
			 * Second we put kids with R-Functions (apart from pivot).
			 * Last we put the pivot, if a kid. */
			if(!(element->off=malloc(sizeof(int)*(2+nk*2)))) ABT_FUNC(MMsg);
			element->rf=element->off+nk;
			element->n_off=nk;
			element->rf[0]=id_array[ids-1].rfp;
			element->rf[1]=id_array[idd-1].rfp;
			for(k=j=0;j<nk;j++) {
				kid=ff->kids[j];
				if(id_array[kid].rfp<0 && kid+1!=pivot) {
					element->rf[2+k]=id_array[kid].rfp;
					element->off[k++]=kid;
				}
			}
			for(j=0;j<nk;j++) {
				kid=ff->kids[j];
				if(id_array[kid].rfp>=0 && kid+1!=pivot) {
					element->rf[2+k]=id_array[kid].rfp;
					element->off[k++]=kid;
				}
			}
			if(k<nk) {
				element->rf[2+k]=id_array[pivot-1].rfp;
				element->off[k]=pivot-1;
			}
			if(pivot) {
				order[pivot-id-1]--;
				element->out_index=n_rf;
				id_array[pivot-1].rfp=n_rf++;
			} else element->out_index=-1;
			print_peelseq_element(stdout,pp);
			pp=&element->next;
			pp->type=0;
			i=-1;
		}
	}
	*n_ops=n_rf;
	if(nf) *n_ops+=get_complex_peel_sequence(pp,comp,fam,famlist,nf);
	free(order);
}

struct Peelseq_Head *peel_init(int *n_ops)
{
	int i,comp,*nfams;
	struct Peelseq_Head *pp;
	struct family **families;
	
	/* Initialize family structures */
	if(!(families=malloc(sizeof(void *)*n_comp))) ABT_FUNC(MMsg);
 	if(!(nfams=malloc(sizeof(int)*n_comp))) ABT_FUNC(MMsg);
	for(i=0;i<n_comp-singleton_flag;i++) families[i]=get_families(i,nfams+i);
	/* Determine peeling sequence for each component */
	if(!(pp=malloc(sizeof(struct Peelseq_Head)*(n_comp-singleton_flag)))) ABT_FUNC(MMsg);
	for(comp=0;comp<n_comp-singleton_flag;comp++) {
		get_peel_sequence(pp+comp,comp,families[comp],nfams[comp],n_ops+comp);
	}
	for(i=0;i<n_comp-singleton_flag;i++) {
		free(families[i][0].kids);
		free(families[i]);
	}
	free(families);
	free(nfams);
	return pp;
}

void free_peel_sequence(struct Peelseq_Head *pp)
{
	int i;
	struct Peelseq_Head *p;
	
	for(i=0;i<n_comp-singleton_flag;i++) {
		p=pp+i;
		free_peelseq(p);
	}
	free(pp);
}
