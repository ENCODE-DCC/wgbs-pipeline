/****************************************************************************
*                                                                          *
*     Loki - Programs for genetic analysis of complex traits using MCMC    *
*                                                                          *
*                     Simon Heath - CNG, Evry                              *
*                                                                          *
*                           May 2003                                       *
*                                                                          *
* get_peelseq.c:                                                           *
*                                                                          *
* Get peeling sequence                                                     *  
*                                                                          *
* Copyright (C) Simon C. Heath 2003                                        *
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
#ifdef HAVE_LIMITS_H
#include <limits.h>
#endif
#include <assert.h>

#include "utils.h"
#include "loki.h"
#include "loki_peel.h"
#include "loki_utils.h"
#include "get_peelseq.h"
#include "lk_malloc.h"
#include "min_deg.h"

static struct loki *lk;

static void *get_new_elem(struct Peelseq_Head **p,int type)
{
	struct Simple_Element *s_elem;
	struct Complex_Element *c_elem;
	
	(*p)->type=type;
	if(type==PEEL_SIMPLE) {
		s_elem=lk_malloc(sizeof(struct Simple_Element));
		(*p)->ptr.simple=s_elem;
		*p=&s_elem->next;
		(*p)->type=0;
		return (void *)s_elem;
	} else if(type==PEEL_COMPLEX) {
		c_elem=lk_malloc(sizeof(struct Complex_Element));
		(*p)->ptr.complex=c_elem;
		*p=&c_elem->next;
		(*p)->type=0;
		return (void *)c_elem;
	} else ABT_FUNC("Invalid element type\n");
	return 0;
}

static int *assemble_matrix(int size,int *inv,int *g_perm,int *trans,int k2,const int linktype,struct Id_Record *id_array)
{
	int i,x,i1,j,k,n_inv,ptr,sex,*mat,mat_size;
	struct Id_Record *kid,**kids;
	
	if(!size) return 0;
	mat_size=size*2+1;
	mat=lk_malloc(sizeof(int)*mat_size);
	ptr=size+1;
	for(x=0;x<size;x++) {
		i1=g_perm[x];
		i=abs(i1)-1;
		n_inv=0;
		sex=id_array[i].sex;
		assert(sex>=0 && sex<=2);
		assert(linktype==LINK_AUTO || sex);
		/* Add first allele */
		inv[n_inv++]=i1;
		/* And second allele if required */
		k=i1<0?FIXED_MAT:FIXED_PAT;
		if((id_array[i].flag&(JNT_LINK|k))==JNT_LINK) inv[n_inv++]=-i1;
		/* Put in parents */
		if(id_array[i].flag&PAR_LINK) {
			if(i1>0) { 
				/* Mother */
				j=id_array[i].dam;
				if(!(id_array[j-1].flag&FIXED_MAT))	inv[n_inv++]=j;
				if(!(id_array[j-1].flag&FIXED_PAT))	inv[n_inv++]=-j;
			} else {
				/* Father (if required) */
				j=id_array[i].sire;
				if(linktype==LINK_AUTO) {
					if(!(id_array[j-1].flag&FIXED_MAT))	inv[n_inv++]=j;
					if(!(id_array[j-1].flag&FIXED_PAT))	inv[n_inv++]=-j;
				} else if(sex==2 && !(id_array[j-1].flag&FIXED_MAT)) inv[n_inv++]=j;
			}
		}
		/* Kids */
		if(id_array[i].flag&KID_LINK) {
			kids=id_array[i].kids;
			assert(sex);
			for(j=0;j<id_array[i].nkids;j++) {
				kid=kids[j];
				if(kid->flag&PAR_LINK) {
					if(sex==2 && !(kid->flag&FIXED_MAT)) inv[n_inv++]=kid->idx+1;
					else if((linktype==LINK_AUTO || kid->sex==2) && !(kid->flag&FIXED_PAT)) inv[n_inv++]=-(kid->idx+1);
				}
			}
		}
		/* Add to matrix */
		/* Count no. entries for this row */
		for(j=0,k=1;k<n_inv;k++) {
			inv[k]=trans[k2+inv[k]];
			if(inv[k]<x) j++;
		}
		/* Check space */
		if(mat_size<ptr+j) {
			do mat_size*=1.5; while(mat_size<ptr+j);
			mat=lk_realloc(mat,sizeof(int)*mat_size);
		}
		mat[x]=ptr;
		for(k=1;k<n_inv;k++)	if(inv[k]<x) mat[ptr++]=inv[k];
	}
	mat[x]=ptr;
	return mat;
}

static int *min_degree(int *g_perm,int *trans,int k2,int *inv,int size,const int linktype,int *nhaps[2],int *ngens,struct Id_Record *id_array)
{
	int i,j,k,*mat,*order,*order_bk,*t,*group,best_fg=0;
	double *wt,cost,best=-1.0,z;
	struct pair_wt *wt1=0;
	
	wt=malloc(sizeof(double)*size);
	if(ngens) {
		wt1=malloc(sizeof(struct pair_wt)*size);
		for(i=0;i<size;i++) wt1[i].pair_node=-1;
	}
	if(nhaps[0]) {
		for(i=0;i<size;i++) {
			j=g_perm[i];
			if(j>0) {
				wt[i]=log((double)nhaps[X_MAT][j-1]);
				if((id_array[j-1].flag&(JNT_LINK|FIXED_PAT))==JNT_LINK) {
					k=trans[k2-j];
					wt1[i].pair_node=k;
					wt1[k].pair_node=i;
					wt1[i].wt=wt1[k].wt=log((double)ngens[j-1]);
				}
			} else wt[i]=log(nhaps[X_PAT][-1-j]);
		}
	} else {
		z=log(2.0);
		for(i=0;i<size;i++) wt[i]=0;
	}
	order=lk_malloc(sizeof(int)*size*2);
	order_bk=lk_malloc(sizeof(int)*size*2);
	group=order+size;
	for(i=0;i<4;i++) {
		mat=assemble_matrix(size,inv,g_perm,trans,k2,linktype,id_array);
		assert(mat);
		greedy(size,mat,order,group,wt,wt1,i,&cost);
		if(best<0 || cost<best) {
			best=cost;
			best_fg=i;
			t=order;
			order=order_bk;
			group=order+size;
			order_bk=t;
		}
		free(mat);
	}
	free(order);
	if(wt1) free(wt1);
	free(wt);
	return order_bk;
}

static void make_ind_rf(struct R_Func *rf,struct Id_Record *id,int ltype)
{
	int k,ng,fg=0,*p;
	
	if((ltype==LINK_AUTO || id->sex==2) && !(id->flag&FIXED_PAT)) fg|=1;
	if(!(id->flag&FIXED_MAT)) fg|=2;
	ng=(fg==3)?2:1;
	p=lk_malloc(sizeof(int)*ng);
	rf->id_list=p;
	k=id->idx+1;
	ng=0;
	if(fg&2) p[ng++]=k;
	if(fg&1) p[ng++]=-k;
	rf->n_ind=ng;
	rf->flag=1;
}

static rf_node *get_new_rfn(rf_node **f_list)
{
	rf_node *ff;
	
	ff=*f_list;
	if(ff) *f_list=ff->next;
	else ff=lk_malloc(sizeof(rf_node));
	ff->next=0;
	return ff;
}

static int add_to_list(int *list,int *n,int c)
{
	int i;
	  
	for(i=0;i<(*n);i++) if(list[i]==c) break;
	if(i==(*n)) list[(*n)++]=c;
	return i;
}


struct Peelseq_Head *get_peelseq(struct Locus *loc,struct loki *loki,int ltype)
{
	int i,i1,j,k,k1,k2,comp,n_comp,*prune,fam1,fam2,*g_perm,*order,*inv,n_out;
	int nf,cs,cst,*ngens=0,p_idx[2],sx,pivot,n_rf,n_genes,n_ops;
	int ped_size,*group,k3,k4,l,sex,n_peel,*flags,kk,tot_out,*nhaps[2]={0,0};
	double wt,tot_wt=0.0,max_wt=0.0;
	struct Peelseq_Head *peel,*pp;
	nuc_fam *fam;
	struct Id_Record *id_array,*par[2],**kids,*kid,*id;
	struct Simple_Element *elem,*single_elem;
	struct Complex_Element *c_elem;
	struct Marker *mark=0;
	static struct R_Func *rf;
	static rf_node *free_rfn_list,**rfn;
	static int rf_size,*trans,*perm;
	rf_node *rfp,*rfp1,*rfp_list,**rfpp;
	int peel_option=3,n_all,n_bits;
	
	if(!loc) {
		if(trans) {
			free(trans);
			free(perm);
			free(rfn);
			free(rf);
			rfp=free_rfn_list;
			while(rfp) {
				rfp1=rfp->next;
				free(rfp);
				rfp=rfp1;
			}
			rf_size=0;
			rf=0;
			rfn=0;
			free_rfn_list=0;
			trans=perm=0;
		}
		return 0;
	}
	message(DEBUG_MSG,"Generating %s peel sequence for locus %s\n",ltype==LINK_AUTO?"autosomal":"X-linked",loc->name);
	prune=loc->pruned_flag;
	fam=loki->family->families;
	fam1=0;
	n_comp=loki->pedigree->n_comp;
	ped_size=loki->pedigree->ped_size;
	peel=lk_malloc(sizeof(struct Peelseq_Head)*n_comp);
	id_array=loki->pedigree->id_array;
	if(loc->type==ST_MARKER) {
		mark=loki->markers->marker+loc->index;
		ngens=mark->ngens;
		if(peel_option&1) for(i=0;i<2;i++) nhaps[i]=mark->nhaps[i];
	}
	if(!trans) {
		trans=lk_malloc(sizeof(int)*(ped_size+1)*2);
		rf_size=loki->family->cp_start[1];
		if(rf_size<16) rf_size=16;
		rf=lk_malloc(sizeof(struct R_Func)*rf_size);
		k=loki->pedigree->comp_size[0];
		k1=k3=loki->family->cp_start[1];
		for(comp=1;comp<n_comp;comp++) {
			k2=loki->pedigree->comp_size[comp];
			if(k2>k) k=k2;
			k2=loki->family->cp_start[comp+1];
			if(k2-k3>k1) k1=k2-k3;
			k3=k2;
		}
		rfn=lk_malloc(sizeof(void *)*2*k);
		perm=lk_malloc(sizeof(int)*k1);
	}
	for(comp=0;comp<n_comp;comp++) {
		single_elem=0;
		n_ops=n_rf=tot_out=0;
		pp=peel+comp;
		pp->type=0;
		fam2=loki->family->cp_start[comp+1];
		for(nf=0,i=fam1;i<fam2;i++) {
			/* Check if pruned */
			par[X_PAT]=fam[i].father;
			par[X_MAT]=fam[i].mother;
			kids=fam[i].kids;
			k=k1=0;
			while((kid=kids[k1++])) if(!prune[kid->idx]) k++;
			if(k) {
				/* Family not pruned */
				if(!par[X_PAT] || prune[par[X_PAT]->idx]) {
					/* Singleton family ? */
					/* If parents are pruned, there should be only 1 unpruned kid.
					* If kid has unpruned kids then it will be handled in the that family, otherwise
					* set up a singleton op */
					if(par[X_PAT]) {
						assert(prune[par[X_MAT]->idx] && k==1);
						k1=0;
						while((kid=kids[k1++])) {
							if(!prune[kid->idx]) break;
						}
						for(k1=0;k1<kid->nkids;k1++) {
							if(!prune[kid->kids[k1]->idx]) break;
						}
						/* Child has unpruned kids, so don't set up singleton op */
						if(k1<kid->nkids) k=0;
					}
					if(k) {
						if(single_elem) {
							elem=single_elem;
							elem->n_off+=k;
							elem->off=lk_realloc(elem->off,sizeof(int)*elem->n_off);
						} else {
							elem=get_new_elem(&pp,PEEL_SIMPLE);
							n_ops++;
							elem->sire=elem->dam=elem->pivot=0;
							elem->out_index=-1;
							elem->off=malloc(sizeof(int)*k);
							elem->n_off=k;
							single_elem=elem;
						}
						k=elem->n_off;
						if(k>loki->peel->max_peel_off) loki->peel->max_peel_off=k;
						k1=0;
						while((kid=kids[k1++])) if(!prune[kid->idx]) elem->off[--k]=kid->idx+1;
					}
				} else perm[nf++]=i;
			} else if(par[X_PAT]) { 
				/* Check if parents need to be mopped up in a singleton op */
				for(k1=0;k1<2;k1++) {
					if(!prune[par[k1]->idx] && !par[k1]->family) {
						for(k2=0;k2<par[k1]->nkids;k2++) {
							if(!prune[par[k1]->kids[k2]->idx]) break;
						}
						if(k2==par[k1]->nkids) k++;
					}
				}
				if(k) {
					if(single_elem) {
						elem=single_elem;
						elem->n_off+=k;
						elem->off=lk_realloc(elem->off,sizeof(int)*elem->n_off);
					} else {
						elem=get_new_elem(&pp,PEEL_SIMPLE);
						n_ops++;
						elem->sire=elem->dam=elem->pivot=0;
						elem->out_index=-1;
						elem->off=lk_malloc(sizeof(int)*k);
						elem->n_off=k;
						single_elem=elem;
					}
					k=elem->n_off;
					if(k>loki->peel->max_peel_off) loki->peel->max_peel_off=k;
					for(k1=0;k1<2;k1++) {
						if(!prune[par[k1]->idx] && !par[k1]->family) {
							for(k2=0;k2<par[k1]->nkids;k2++) {
								if(!prune[par[k1]->kids[k2]->idx]) break;
							}
							if(k2==par[k1]->nkids) elem->off[--k]=par[k1]->idx+1;
						}
					}
				}
			}
		}
		/* Get order for each pedigree member */
		cs=loki->pedigree->comp_size[comp];
		cst=loki->pedigree->comp_start[comp];
		for(j=0;j<cs;j++) {
			id_array[cst+j].flag=0;
			id_array[cst+j].rfp=-1;
		}
		for(i=0;i<nf;i++) {
			fam1=perm[i];
			p_idx[X_MAT]=fam[fam1].mother->idx;
			p_idx[X_PAT]=fam[fam1].father->idx;
			kids=fam[fam1].kids;
			k=-1;
			while((kid=kids[++k])) if(!prune[kid->idx]) {
				id_array[kid->idx].flag++;
			}
				/* Add parent links */
				id_array[p_idx[X_MAT]].flag++; 
			/* Don't add father for X_linked data *unless* he has female (unpruned, unfixed) offspring */
			id_array[p_idx[X_PAT]].flag++;
		}
		/* Look for simple peeling sequence */
		for(i=nf-1;i>=0;i--) {
			if(!(peel_option&2)) continue;
			fam1=perm[i];
			par[X_MAT]=fam[fam1].mother;
			par[X_PAT]=fam[fam1].father;
			pivot=0;
			for(sx=0;sx<2;sx++) {
				if((!ngens || ngens[par[sx]->idx]>1) && par[sx]->flag>1) {
					if(pivot) pivot=-2;
					else pivot=par[sx]->idx+1;
				}
			}
			kids=fam[fam1].kids;
			k=-1;
			k1=0;
			while((kid=kids[++k])) {
				if(!prune[kid->idx]) {
					k1++;
					if(pivot<0) {
						if(!ngens || ngens[kid->idx]>1) break;
					} else if((!ngens || ngens[kid->idx]>1) && kid->flag>1) {
						if(pivot) break;
						pivot=kid->idx+1;
					}
				}
			}
			if(!kid) {
					elem=get_new_elem(&pp,PEEL_SIMPLE);
					n_ops++;
					elem->sire=par[X_PAT]->idx+1;
					elem->dam=par[X_MAT]->idx+1;
					elem->pivot=pivot;
					for(k=0;k<2;k++) {
						if(par[k]->rfp>=0) par[k]->rfp=-1;
						if(par[k]->idx+1==pivot || pivot<0) par[k]->flag--;
						else if(par[k]->flag>1) {
							par[k]->flag--;
							if(k) elem->sire=-elem->sire;
							else elem->dam=-elem->dam;
						} else par[k]->flag=DONE_LINK;
					}
					assert(k1);
					elem->off=lk_malloc(sizeof(int)*k1);
					k=-1;
					k1=0;
					while((kid=kids[++k])) if(!prune[kid->idx]) {
						elem->off[k1]=kid->idx+1;
						if(kid->rfp>=0) kid->rfp=-1;
						if(pivot==kid->idx+1) kid->flag--;
						else if(kid->flag>1) {
							kid->flag--;
							elem->off[k1]=-elem->off[k1];
						} else kid->flag=DONE_LINK;
						k1++;
					}
						elem->n_off=k1;
					if(k1>loki->peel->max_peel_off) loki->peel->max_peel_off=k1;
					if(pivot) {
						if(pivot>0) {
							id_array[pivot-1].rfp=n_rf;
							elem->out_index=n_rf++;
							tot_out+=2;
						} else {
							elem->out_index=n_rf;
							par[X_PAT]->rfp=n_rf++;
							par[X_MAT]->rfp=n_rf++;
							tot_out+=4;
						}
					} else elem->out_index=-1;
					perm[i]=perm[--nf];
					i=nf;
				}
		}
			if(nf) {
				message(DEBUG_MSG,"Generating complex peel sequence for %d famil%s in component %d\n",nf,nf==1?"y":"ies",comp+1);
				if(loc->type==ST_MARKER) n_all=mark->n_all1[comp];
				else n_all=loc->n_alleles;
				n_bits=num_bits(n_all);
				/* Count how many genes are left (rough overestimate at this stage) */
				n_genes=0;
				for(i=0;i<nf;i++) {
					fam1=perm[i];
					/*							print_orig_triple(stdout,fam[fam1].father->idx+1);
					printf(" %d\n",mark->haplo[fam[fam1].father->idx]);
					print_orig_triple(stdout,fam[fam1].mother->idx+1);
					printf(" %d\n",mark->haplo[fam[fam1].mother->idx]); */
					kids=fam[fam1].kids;
					k=k1=k2=0;
					while((kid=kids[k++])) if(!prune[kid->idx] && (kid->flag&~DONE_LINK)<=1) {
						/*								print_orig_triple(stdout,kid->idx+1); 
						printf(" %d\n",mark->haplo[kid->idx]); */
						if(kid->sex==2) k1++;
						else k2++;
					}
						n_genes+=(ltype==LINK_AUTO?4+2*(k1+k2):3+2*k1+k2);
				}
				if(n_genes>cs*2) n_genes=cs*2;
				assert(n_genes);
				/* Allocate space for R-Functions */
				if(rf_size<n_rf) {
					rf_size=n_rf*1.5;
					free(rf);
					rf=lk_malloc(sizeof(struct R_Func)*rf_size);
				}
				for(i=0;i<n_genes;i++) rfn[i]=0;
				/* Allocate space and store gene list
					Also make required R-Functions */
				g_perm=lk_malloc(sizeof(int)*3*n_genes);
				inv=g_perm+n_genes;
				flags=inv+n_genes;
				memset(flags,0,sizeof(int)*n_genes);
				n_genes=0;
#ifndef NDEBUG
				for(j=0;j<cs;j++) {
					k=cst+j+1;
					trans[k+ped_size]=-1;
					trans[ped_size-k]=-1;
				}
#endif
				for(i=0;i<nf;i++) {
					fam1=perm[i];
					par[X_MAT]=fam[fam1].mother;
					par[X_PAT]=fam[fam1].father;
					if(!(par[X_MAT]->flag&(KID_LINK|JNT_LINK|DONE_LINK))) {
						j=par[X_MAT]->idx+1;
						k=par[X_MAT]->rfp;
						if(!nhaps[0] || nhaps[X_MAT][j-1]>1) {
							trans[j+ped_size]=n_genes;
							g_perm[n_genes++]=j;
						} else {
							assert(k<0);
							par[X_MAT]->flag|=FIXED_MAT;
						}
						if(!nhaps[0] || nhaps[X_PAT][j-1]>1) {
							trans[ped_size-j]=n_genes;
							g_perm[n_genes++]=-j;
						} else {
							par[X_MAT]->flag|=FIXED_PAT;
							assert(k<0);
						}
						par[X_MAT]->flag|=(KID_LINK|JNT_LINK|DONE_LINK);
						if(k>=0) {
							assert(!ngens || ngens[j-1]!=1); /* Make sure fixed individual is not chosen as pivot */
							make_ind_rf(rf+k,par[X_MAT],ltype);
							k1=n_genes;
							if(!(par[X_MAT]->flag&FIXED_PAT)) {
								rfn[--k1]=get_new_rfn(&free_rfn_list);
								rfn[k1]->rf_idx=k;
							}
							if(!(par[X_MAT]->flag&FIXED_MAT)) {
								rfn[--k1]=get_new_rfn(&free_rfn_list);
								rfn[k1]->rf_idx=k;
							}
						}
					}
					if(!(par[X_PAT]->flag&(KID_LINK|JNT_LINK|DONE_LINK))) {
						j=par[X_PAT]->idx+1;
						k=par[X_PAT]->rfp;
						if(!nhaps[0] || nhaps[X_MAT][j-1]>1) {
							trans[j+ped_size]=n_genes;
							g_perm[n_genes++]=j;
						} else {
							assert(k<0);
							par[X_PAT]->flag|=FIXED_MAT;
						}
						if(ltype==LINK_AUTO) {
							if(!nhaps[0] || nhaps[X_PAT][j-1]>1) {
								trans[ped_size-j]=n_genes;
								g_perm[n_genes++]=-j;
							} else {
								assert(k<0);
								par[X_PAT]->flag|=FIXED_PAT;
							}
							par[X_PAT]->flag|=JNT_LINK;
						}
						par[X_PAT]->flag|=DONE_LINK;
						if(k>=0) {
							assert(!ngens || ngens[j-1]!=1); /* Make sure fixed individual is not chosen as pivot */
							make_ind_rf(rf+k,par[X_PAT],ltype);
							k1=n_genes;
							if(ltype==LINK_AUTO && !(par[X_PAT]->flag&FIXED_PAT)) {
								rfn[--k1]=get_new_rfn(&free_rfn_list);
								rfn[k1]->rf_idx=k;
							}
							if(!(par[X_PAT]->flag&FIXED_MAT)) {
								rfn[--k1]=get_new_rfn(&free_rfn_list);
								rfn[k1]->rf_idx=k;
							}
						}
					}
					kids=fam[fam1].kids;
					k=k1=0;
					while((kid=kids[k++])) if(!prune[kid->idx]) {
						if((kid->flag&~LINK_FLAGS)<=1) {
							j=kid->idx+1;
							k2=kid->rfp;
							if(!nhaps[0] || nhaps[X_MAT][j-1]>1) {
								trans[j+ped_size]=n_genes;
								g_perm[n_genes++]=j;
							} else kid->flag|=FIXED_MAT;
							if(ltype==LINK_AUTO || kid->sex==2) {
								if(!nhaps[0] || nhaps[X_PAT][j-1]>1) {
									trans[ped_size-j]=n_genes;
									g_perm[n_genes++]=-j;
								} else kid->flag|=FIXED_PAT;
							}
							if(k2>=0) {
								assert(!ngens || ngens[j-1]!=1); /* Make sure fixed individual is not chosen as pivot */
								make_ind_rf(rf+k2,kid,ltype);
								j=n_genes;
								if((ltype==LINK_AUTO || kid->sex==2) && !(kid->flag&FIXED_PAT)) {
									rfn[--j]=get_new_rfn(&free_rfn_list);
									rfn[j]->rf_idx=k2;
								}
								if(!(kid->flag&FIXED_MAT)) {
									rfn[--j]=get_new_rfn(&free_rfn_list);
									rfn[j]->rf_idx=k2;
								}
							}
							if(ltype==LINK_AUTO || kid->sex==2) {
								kid->flag|=(PAR_LINK|JNT_LINK|DONE_LINK);
								k1++;
							} else kid->flag|=(MAT_LINK|DONE_LINK);
						} else {
							if(ltype==LINK_AUTO || kid->sex==2) {
								k1++;
								kid->flag|=PAR_LINK;
							} else kid->flag|=MAT_LINK;
						}
					}
						if(k1) par[X_PAT]->flag|=KID_LINK;
				}
				/* Get order for peeling genes */
				order=min_degree(g_perm,trans,ped_size,inv,n_genes,ltype,nhaps,ngens,id_array);
				group=order+n_genes;
				for(i=0;i<n_genes;i++) {
					j=order[i];
					k=g_perm[j];
					j=1;
					inv[0]=k;
					for(;i<n_genes-1;i++) {
						if(group[i+1]) break;
						if(n_bits*(j+1)>(int)LK_LONG_BIT) break;
						inv[j++]=g_perm[order[i+1]];
					}
					n_peel=j;
					rfp_list=0;
					k4=0;
					for(k1=0;k1<n_peel;k1++) {
						k=inv[k1];
						id=id_array+abs(k)-1;
						if(!id->sire || prune[id->sire-1]) flags[k1]|=HAP_FND;
						sex=id->sex;
						/* Add in other allele if required */
						if(id->flag&JNT_LINK) {
							flags[k1]|=HAP_JNT;
							flags[add_to_list(inv,&j,-k)]|=HAP_JNT;
							id->flag&=~JNT_LINK;
						}
						/* Put in parents */
						if(k>0) {
							if(id->flag&MAT_LINK) {
								add_to_list(inv,&j,id->dam);
								add_to_list(inv,&j,-id->dam);
								flags[k1]|=HAD_M;
								id->flag&=~MAT_LINK;
							}
						} else {
							if(id->flag&PAT_LINK) {
								if(ltype==LINK_AUTO) {
									add_to_list(inv,&j,id->sire);
									add_to_list(inv,&j,-id->sire);
									flags[k1]|=HAD_P;
								} else if(sex==2) {
									add_to_list(inv,&j,id->sire);
									flags[k1]|=HAD_P;
								}
								id->flag&=~PAT_LINK;
							}
						}
						/* and kids... */
						if(id->flag&KID_LINK) {
							kids=id->kids;
							if(sex==2) {
								for(i1=0;i1<id->nkids;i1++) {
									kid=kids[i1];
									if(kid->flag&MAT_LINK) {
										k2=add_to_list(inv,&j,kid->idx+1);
										flags[k2]|=HAD_M;
										kid->flag&=~MAT_LINK;
									}
								}
							} else {
								for(i1=0;i1<id->nkids;i1++) {
									kid=kids[i1];
									if(kid->flag&PAT_LINK) {
										k2=add_to_list(inv,&j,-(kid->idx+1));
										flags[k2]|=HAD_P;
										kid->flag&=~PAT_LINK;
									}
								}
							}
						}
						/* And R-Functions for the peeled genes */
						kk=trans[k+ped_size];
						assert(kk>=0);
						rfp=rfn[kk];
						if(rfp) flags[k1]|=IN_RF;
						while(rfp) {
							k2=rfp->rf_idx;
							rfp1=rfp->next;
							if(rf[k2].flag) {
								for(k3=0;k3<rf[k2].n_ind;k3++) {
									l=rf[k2].id_list[k3];
									if(l!=k) {
										kk=add_to_list(inv,&j,l);
										flags[kk]|=IN_RF;
									}
								}
								kk=trans[k+ped_size];
								rfp->next=rfp_list;
								rfp_list=rfp;
								k4++;
								rf[k2].flag=0;
							} else {
								rfp->next=free_rfn_list;
								free_rfn_list=rfp;
							}
							rfp=rfp1;
						}
					}
					k3=n_peel;
					for(k=k3;k<j;k++) {
						k1=inv[k];
						k2=0;
						if(k1<0) {
							k2=(id_array[-1-k1].flag&FIXED_PAT);
						} else k2=(id_array[k1-1].flag&FIXED_MAT);
						if(k2&FIXED_FLAGS) {
							flags[k]|=FIXED_FLAG;
							k2=flags[k];
							inv[k]=inv[k3];
							flags[k]=flags[k3];
							inv[k3]=k1;
							flags[k3++]=k2;
						}
					}
					/* Calculate weight */
					if(ngens) {
						wt=1.0;
						for(k=0;k<j;k++) {
							k1=inv[k];
							for(k2=0;k2<j;k2++) if(inv[k2]==-k1) break;
							if(k2<j) {
								if(k1>0) wt*=(double)ngens[k1-1];
							} else {
								if(k1<0) wt*=(double)nhaps[X_PAT][-1-k1];
								else wt*=(double)nhaps[X_MAT][k1-1];
							}
						}
					} else {
						k2=2<<j;
						wt=(double)k2;
					}
					tot_wt+=wt;
					if(wt>max_wt) max_wt=wt;
					/* Construct peel op */
					c_elem=get_new_elem(&pp,PEEL_COMPLEX);
					n_ops++;
					c_elem->involved=lk_malloc(sizeof(int)*(j*2+k4));
					c_elem->flags=c_elem->involved+j;
					c_elem->index=c_elem->flags+j;
					for(k=n_peel;k<j;k++) {
						k2=inv[k];
						if(flags[k]&FIXED_FLAG) continue;
						k1=trans[k2+ped_size];
						assert(k1>=0);
						if(flags[k]&IN_RF) {
							rfp=rfn[k1];
							rfpp=rfn+k1;
							while(rfp) {
								rfp1=rfp->next;
								k2=rfp->rf_idx;
								if(!rf[k2].flag) {
									rfp->next=free_rfn_list;
									free_rfn_list=rfp;
									*rfpp=rfp1;
								} else rfpp=&rfp->next;
								rfp=rfp1;
							}
						}
						if(n_peel<j) {
							rfp=get_new_rfn(&free_rfn_list);
							rfp->rf_idx=n_rf;
							rfp->next=rfn[k1];
							rfn[k1]=rfp;
						}
					}	
					for(k=0;k<j;k++) {
						c_elem->flags[k]=flags[k];
						c_elem->involved[k]=inv[k];
						flags[k]=0;
					}
					k4=0;
					while(rfp_list) {
						rfp=rfp_list->next;
						k2=rfp_list->rf_idx;
						c_elem->index[k4++]=k2;
						free(rf[k2].id_list);
						rf[k2].id_list=0;
						rfp_list->next=free_rfn_list;
						free_rfn_list=rfp_list;
						rfp_list=rfp;
					}
					c_elem->n_involved=j;
					c_elem->n_peel=n_peel;
					c_elem->n_rfuncs=k4;
					/* Make output R-Function */
					n_out=0;
					for(k=n_peel;k<j;k++) {
						if(!(c_elem->flags[k]&FIXED_FLAG)) n_out++;
					}
					c_elem->n_out=n_out;
					if(n_out) {
						if(n_rf==rf_size) {
							rf_size*=1.5;
							rf=lk_realloc(rf,sizeof(struct R_Func)*rf_size);
						}
						rf[n_rf].id_list=lk_malloc(sizeof(int)*n_out);
						for(k=0,k1=j-n_out;k1<j;k1++) rf[n_rf].id_list[k++]=inv[k1];
						rf[n_rf].flag=1;
						rf[n_rf].n_ind=k;
						c_elem->out_index=n_rf++;
						tot_out+=k;
					} else c_elem->out_index=-1;
				}
				for(i1=cs-1;i1>=0;i1--) {
					k=cst+i1;
					id=id_array+k;
					k1=id->flag;
					if(k1&(JNT_LINK|PAR_LINK)) {
						assert(k1&FIXED_FLAGS); /* Is individual peeled properly? */
						sex=id->sex;
						j=0;
						if(k1&JNT_LINK) {
							inv[j]=k+1;
							flags[j++]|=HAP_JNT;
							inv[j]=-(k+1);
							flags[j++]|=HAP_JNT;
						} 
						if(k1&MAT_LINK) {
							flags[add_to_list(inv,&j,k+1)]|=HAD_M;
							id->flag&=~MAT_LINK;
							if(id_array[id->dam-1].flag&JNT_LINK) {
								flags[add_to_list(inv,&j,id->dam)]|=HAP_JNT;
								flags[add_to_list(inv,&j,-id->dam)]|=HAP_JNT;
								id_array[id->dam-1].flag&=~JNT_LINK;
							} else {
								add_to_list(inv,&j,id->dam);
								add_to_list(inv,&j,-id->dam);
							}
						}
						if(k1&PAT_LINK) {
							flags[add_to_list(inv,&j,-(k+1))]|=HAD_P;
							id->flag&=~PAT_LINK;
							if(ltype==LINK_AUTO) {
								if(id_array[id->sire-1].flag&JNT_LINK) {
									flags[add_to_list(inv,&j,id->sire)]|=HAP_JNT;
									flags[add_to_list(inv,&j,-id->sire)]|=HAP_JNT;
									id_array[id->sire-1].flag&=~JNT_LINK;
								} else {
									add_to_list(inv,&j,id->sire);
									add_to_list(inv,&j,-id->sire);
								}
							} else if(sex==2) add_to_list(inv,&j,id->sire);
						}
						for(k=0;k<j;k++) flags[k]|=FIXED_FLAG;
						c_elem=get_new_elem(&pp,PEEL_COMPLEX);
						n_ops++;
						c_elem->involved=lk_malloc(sizeof(int)*j*2);
						c_elem->flags=c_elem->involved+j;
						c_elem->index=0;
						for(k=0;k<j;k++) {
							k1=inv[k];
							id=id_array+abs(k1)-1;
							if(k1<0) {
								if(!(id->flag&DONE_FND_PAT) && (!id->sire || prune[id->sire-1])) {
									flags[k]|=HAP_FND;
									id->flag|=DONE_FND_PAT;
								}
							} else {
								if(!(id->flag&DONE_FND_MAT) && (!id->dam || prune[id->dam-1])) {
									flags[k]|=HAP_FND;
									id->flag|=DONE_FND_MAT;
								}
							}
							c_elem->flags[k]=flags[k];
							c_elem->involved[k]=inv[k];
							flags[k]=0;
						}
						c_elem->n_involved=j;
						c_elem->n_peel=c_elem->n_rfuncs=c_elem->n_out=0;
						c_elem->out_index=-1;
					}
				}
				free(order);
				free(g_perm);
			}
			if(n_ops>loki->peel->max_peel_ops) loki->peel->max_peel_ops=n_ops;
			if(n_rf>loki->peel->max_rfuncs) loki->peel->max_rfuncs=n_rf;
			if(tot_out>loki->peel->max_out) loki->peel->max_out=tot_out;
			fam1=fam2;
	}
		if(tot_wt) message(DEBUG_MSG,"%s - total cost = %g, max op = %g\n",loc->name,tot_wt,max_wt);
		return peel;
}

static void Free_Peel_Seq(void)
{
	int i,comp;
	
	if(lk->peel) {
		if(lk->peel->peelseq_head) {
			for(i=0;i<lk->markers->n_markers;i++) if(lk->peel->peelseq_head[i]) {
				for(comp=0;comp<lk->pedigree->n_comp;comp++)
					free_peelseq(lk->peel->peelseq_head[i]+comp);
				free(lk->peel->peelseq_head[i]);
			}
				free(lk->peel->peelseq_head);
		}
		if(lk->peel->peelseq_head_gen) {
			for(comp=0;comp<lk->pedigree->n_comp;comp++)
				free_peelseq(lk->peel->peelseq_head_gen+comp);
			free(lk->peel->peelseq_head_gen);
		}
		if(lk->peel->peelseq_head_gen_x) {
			for(comp=0;comp<lk->pedigree->n_comp;comp++)
				free_peelseq(lk->peel->peelseq_head_gen_x+comp);
			free(lk->peel->peelseq_head_gen_x);
		}
		free(lk->peel);
		lk->peel=0;
	}
}

void Get_Peel_Seq(struct loki *loki)
{
	int i,n_markers,ltype,flag=0;
	struct Locus *loc;
	
	message(INFO_MSG,"Generating marker peeling sequences\n");
	loki->peel=lk_calloc((size_t)1,sizeof(struct lk_peel));
	loki->peel->seg_count=0;
	loki->peel->aff_freq=0;
	loki->peel->tl_group=0;
	loki->peel->workspace.r_funcs=0;
	loki->peel->workspace.s0=0;
	loki->peel->workspace.s2=loki->peel->workspace.s7=0;
	loki->peel->workspace.s1=loki->peel->workspace.s4=loki->peel->workspace.s8=0;
	loki->peel->workspace.s3=0;
	loki->peel->workspace.s5=loki->peel->workspace.s6=0;
	for(i=0;i<2;i++) {
		loki->peel->first_mem_block[i]=0;
		loki->peel->mem_block[i]=0;
	}
	n_markers=loki->markers->n_markers;
	if(n_markers) {
		loki->peel->peelseq_head=lk_malloc(sizeof(void *)*n_markers);
		for(i=0;i<n_markers;i++) {
			loc=&loki->markers->marker[i].locus;
			ltype=loki->markers->linkage[loc->link_group].type;
			if(ltype&LINK_PSEUDO) {
				loki->peel->peelseq_head[i]=0;
				continue;
			}
			ltype&=LINK_TYPES_MASK;
			if(ltype==LINK_AUTO || ltype==LINK_X || ltype==LINK_Z) {
				loki->peel->peelseq_head[i]=get_peelseq(loc,loki,ltype);
				/*	print_peelseq(stdout,loki->peel->peelseq_head[i]); */
			}
		}
	} else loki->peel->peelseq_head=0;
	if(loki->models->tlocus) {
		for(i=0;i<loki->markers->n_links;i++) {
			ltype=(loki->markers->linkage[i].type&LINK_TYPES_MASK);
			if(ltype==LINK_X) flag=1;
		}
		loki->peel->peelseq_head_gen=get_peelseq(loki->models->tlocus,loki,LINK_AUTO);
		loki->peel->peelseq_head_gen_x=flag?get_peelseq(loki->models->tlocus,loki,LINK_X):0;
	} else loki->peel->peelseq_head_gen=loki->peel->peelseq_head_gen_x=0;
	/* Call cleanup routines */
	get_peelseq(0,0,0);
	min_deg(0,0,0,0,0);
	lk=loki;
	if(atexit(Free_Peel_Seq)) message(WARN_MSG,"Unable to register exit function Free_Peel_Seq()\n");
}
