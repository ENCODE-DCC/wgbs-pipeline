/****************************************************************************
 *                                                                          *
 *     Loki - Programs for genetic analysis of complex traits using MCMC    *
 *                                                                          *
 *                    Simon Heath - CNG, France                             *
 *                                                                          *
 *                         September 2002                                   *
 *                                                                          *
 * fseg_pen.c:                                                              *
 *                                                                          *
 * Calculate penetrance for a (squashed) inheritance vector                 *
 *                                                                          *
 * Copyright (C) Simon C. Heath 2002                                        *
 * This is free software.  You can distribute it and/or modify it           *
 *                                                                          *
 * under the terms of the Modified BSD license, see the file COPYING        *
 *                                                                          *
 ****************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#ifdef USE_DMALLOC
#include <dmalloc.h>
#endif
#include <math.h>

#include "utils.h"
#include "loki.h"
#include "loki_peel.h"
#include "fenris.h"
#include "fenris_peel.h"
#include "calc_pen.h"

static double **gpen,**gpen1,**rf_ptr,*tmp_xx;
static int *order,*olist,*gflag,*gp_inv,rf_idx_n,**fact,**rfx,
  *tt,*cc,*cc1,*gstate,*rfxn,*ofact,*n_all1,*trans1,*trans2;
static struct deg_list *deglist,**first;
static struct gp_rfnode **adj_mat,*freelist;
static struct gp_rfunc *free_rf_list[MAX_FSP_RF_GENES];
static struct gp_rfunc_ptr **rfuncs,*free_rfp_list;
static lk_ulong *required;
  
void setup_fseg_pen(int ng,int n_all)
{
	int i;
	
	if(!(adj_mat=malloc(sizeof(void *)*ng))) ABT_FUNC(MMsg);
	if(!(gpen=malloc(sizeof(void *)*ng*2))) ABT_FUNC(MMsg);
	gpen1=gpen+ng;
	if(!(gpen[0]=malloc(sizeof(double)*ng*n_all))) ABT_FUNC(MMsg);
	for(i=1;i<ng;i++) gpen[i]=gpen[i-1]+n_all;
	if(!(first=malloc(sizeof(void *)*ng))) ABT_FUNC(MMsg);
	if(!(required=malloc(sizeof(lk_ulong)*ng))) ABT_FUNC(MMsg);
	if(!(rfuncs=malloc(sizeof(void *)*ng))) ABT_FUNC(MMsg);
	if(!(fact=malloc(sizeof(void *)*2*ng))) ABT_FUNC(MMsg);
	rfx=fact+ng;
	for(i=0;i<ng;i++) fact[i]=0;
	if(!(order=malloc(sizeof(int)*(ng*11+n_all*2)))) ABT_FUNC(MMsg);
	olist=order+ng;
	gflag=olist+ng;
	gp_inv=gflag+ng;
	tt=gp_inv+ng;
	cc=tt+ng;
	cc1=cc+ng;
	gstate=cc1+ng;
	rfxn=gstate+ng;
	ofact=rfxn+ng;
	n_all1=ofact+ng;
	trans1=n_all1+ng;
	trans2=trans1+n_all;
	if(!(deglist=malloc(sizeof(struct deg_list)*ng))) ABT_FUNC(MMsg);
	for(i=0;i<ng;i++) deglist[i].gene=i;
	if(!(tmp_xx=malloc(sizeof(double)*n_all*n_all))) ABT_FUNC(MMsg);
}

static struct gp_rfunc *get_new_rfunc(int n,int n_all)
{
	int nc;
	struct gp_rfunc *rf;

	if(!(rf=malloc(sizeof(struct gp_rfunc)))) ABT_FUNC(MMsg);
#ifdef TRACE
	if(CHK_TRACE(TRACE_LEVEL_3)) printf("Allocating rf %p\n",rf);
#endif
	rf->n_inv=n;
	nc=(int)(.5+exp(log((double)n_all)*(double)n));
	if(!(rf->p=malloc(sizeof(double)*nc))) ABT_FUNC(MMsg);
	if(!(rf->inv=malloc(sizeof(int)*n))) ABT_FUNC(MMsg);
	return rf;
}

static struct gp_rfnode *add_node(int x,int y)
{
	struct gp_rfnode *p,**pp;
	
	pp=adj_mat+y;
	p=*pp;
	while(p) {
		if(p->x>=x) break;
		pp=&p->next;
		p=*pp;
	}
	if(!p || p->x!=x) {
		if(freelist) {
			p=freelist;
			freelist=freelist->next;
		} else if(!(p=malloc(sizeof(struct gp_rfnode)))) ABT_FUNC(MMsg);
		p->next=*pp;
		p->x=x;
		p->rf=0;
		*pp=p;
	}
	return p;
}

static struct gp_rfnode *find_node(int x,int y)
{
	struct gp_rfnode *p;

	p=adj_mat[y];
	while(p) {
		if(p->x>=x) break;
		p=p->next;
	}
	if(p && p->x!=x) p=0;
	return p;
}

static void del_node(int x,int y)
{
	struct gp_rfnode *p,**pp;

	pp=adj_mat+y;
	p=*pp;
	while(p->x!=x) {
		pp=&p->next;
		p=*pp;
	}
	*pp=p->next;
	p->next=freelist;
	freelist=p;
}

static void del_row(int y)
{
	struct gp_rfnode *p;
	
	p=adj_mat[y];
	if(p) {
		while(p->next) p=p->next;
		p->next=freelist;
		freelist=adj_mat[y];
		adj_mat[y]=0;
	}
}

static void del_row1(int y)
{
	struct gp_rfnode *p;
	
	p=adj_mat[y];
	if(p) {
		while(p->next) {
			del_node(y,p->x);
			p=p->next;
		}
		del_node(y,p->x);
		p->next=freelist;
		freelist=adj_mat[y];
		adj_mat[y]=0;
	}
}

static int calc_degree(int i)
{
	int i1,j,j1,j2,k=0,k1,k2;
	struct deg_list *p;
	struct gp_rfnode *rfn;

	rfn=adj_mat[i];
	gp_inv[k++]=i;
	while(rfn) {
		gp_inv[k++]=rfn->x;
		rfn=rfn->next;
	}
	k2=k;
	gp_inv[k]=-1;
	for(j=1;j<k;j++) {
		i1=gp_inv[j];
		if(!gflag[i1]) {
			rfn=adj_mat[i1];
			if(rfn->next) {
				j2=1;
				while(rfn) {
					k1=rfn->x;
					if(k1!=i) {
						if(gp_inv[j2]==i1) j2++;
						if(gp_inv[j2++]!=k1) break;
						rfn=rfn->next;
					} else {
						rfn=rfn->next;
						if(!rfn && gp_inv[j2]==i1) j2++;
					}
				}
				if(!rfn && j2==k-1 && gp_inv[j2]==i1) k1=1;
				else k1=(!rfn&&j2==k)?1:0;
			} else k1=k==2?1:0;
			if(k1) {
#ifdef TRACE
				if(CHK_TRACE(TRACE_LEVEL_2)) (void)printf("Absorbing gene %d -> %d\n",i1,i);
#endif
				gflag[i1]|=1;
				p=deglist[i1].abs_list;
				if(p) {
					while(p->abs_list) p=p->abs_list;
					p->abs_list=deglist[i].abs_list;
				} else deglist[i1].abs_list=deglist[i].abs_list;
				deglist[i].abs_list=deglist+i1;
				if(deglist[i1].deg>=0) {
					j1=deglist[i1].deg;
#ifdef TRACE
					if(CHK_TRACE(TRACE_LEVEL_2)) (void)printf("Removing %d from degree list %d\n",i1,j1);
#endif
					if(deglist[i1].prev) deglist[i1].prev->next=deglist[i1].next;
					else first[j1]=deglist[i1].next;
					if(deglist[i1].next) deglist[i1].next->prev=deglist[i1].prev;
					deglist[i1].deg=-1;
				}
				rfn=adj_mat[i1];
				while(rfn) {
					k1=rfn->x;
					if(k1!=i && gflag[k1]) {
						j2=deglist[k1].deg;
						if(j2>=0) {
#ifdef TRACE
							if(CHK_TRACE(TRACE_LEVEL_2)) (void)printf("Moving %d from degree list %d -> %d\n",k1,j2,j2-1);
							
#endif
#ifdef DEBUG
							if(!j2) ABT_FUNC("Internal error - degree shouldn't be zero\n");
#endif
							if(deglist[k1].prev) deglist[k1].prev->next=deglist[k1].next;
							else first[j2]=deglist[k1].next;
							if(deglist[k1].next) deglist[k1].next->prev=deglist[k1].prev;
							deglist[k1].next=first[--j2];
							deglist[k1].prev=0;
							if(first[j2]) first[j2]->prev=deglist+k1;
							first[j2]=deglist+k1;
							deglist[k1].deg--;
						}
					}
					rfn=rfn->next;
				}
			}
		}
	}
	return k-1;
}

static void free_rfunc_list(struct gp_rfunc *rfl)
{
	int k;
	struct gp_rfunc *rf;
	
	while(rfl) {
		rf=rfl->next;
		k=rfl->n_inv-2;
#ifdef TRACE
		if(CHK_TRACE(TRACE_LEVEL_3)) (void)printf("[%s:%d] %s() - Adding %p to free_rf_list[%d]\n",__FILE__,__LINE__,__func__,(void *)rfl,k);
#endif
		rfl->next=free_rf_list[k];
		free_rf_list[k]=rfl;
		rfl=rf;
	}
}

/* Only called if we exit with an incompatibility */
static void clean_up(int ngenes,struct gp_rfunc *rfl)
{
	int i;
	struct gp_rfunc_ptr *rfp1,*rfp2;
	struct gp_rfunc *rf;
	
	for(i=0;i<ngenes;i++) {
		rfp1=rfuncs[i];
		while(rfp1) {
			rfp2=rfp1->next;
			rf=rfp1->rf;
			if(!rf->flag) {
				rf->flag=1;
				rf->next=rfl;
				rfl=rf;
			}
#ifdef TRACE
			if(CHK_TRACE(TRACE_LEVEL_3)) (void)printf("[%s:%d] %s() - Adding %p to free_rfp_list\n",__FILE__,__LINE__,__func__,(void *)rfp1);
#endif
			rfp1->next=free_rfp_list;
			free_rfp_list=rfp1;
			rfp1=rfp2;
		}
  		del_row(i);
	}
	while(rfl) {
		rf=rfl->next;
		i=rfl->n_inv-1;
#ifdef TRACE
		if(CHK_TRACE(TRACE_LEVEL_3)) (void)printf("[%s:%d] %s() - Adding %p to free_rf_list[%d]\n",__FILE__,__LINE__,__func__,(void *)rfl,i);
#endif
		rfl->next=free_rf_list[i];
		free_rf_list[i]=rfl;
		rfl=rf;
	}
}

double fseg_pen(int *state[2],int **alls,double **obspen,int *grps,double **freq,int nobs,int ng,int n_all)
{
	int i,j,k,k1,k2=0,o1,o2,o3,n_gen,no,tri_piv,dual_piv;
	int i1,i2,j1,n_inv,n_peel,n_rf,n_out,nc=0,*inv,*inv1,nl1,nl2,nl3;
	struct gp_rfnode *p,*p1=0,*p2;
	struct gp_rfunc *rf,*rf1,*rfl;
	struct gp_rfunc_ptr *rfp,*rfp1,*rfp2,**rfp3;
	double *pd1,*pd2,*pd3,*pd12,*pd13,*pd23,z,z1,z2,zz,l=0.0,*xx=0;
	double **tp1;
	struct deg_list *degp;
	lk_ulong a;
	static int bcount[]={0,1,1,2,1,2,2,3,1,2,2,3,2,3,3,4};

#ifdef TRACE
	struct deg_list *degp1;
#endif
	  
	n_gen=n_all*n_all;
	for(i=0;i<ng;i++) {
		required[i]=0;
		adj_mat[i]=0;
		order[i]=0;
		gflag[i]=0;
	}
	/* Make set of observed alleles for each founder gene */
	for(i=0;i<nobs;i++) {
		o1=state[0][i];
		o2=state[1][i];
		k1=alls[0][i];
		k2=alls[1][i];
		a=(1<<k1)|(1<<k2);
		required[o1]|=a;
		required[o2]|=a;
	}
	for(i=0;i<ng;i++) {
		a=required[i];
		for(k=j=0;j<n_all>>2;j++) {
			k+=bcount[(int)(a&0xf)];
			a>>=4;
		}
		for(j=0;j<(n_all&3);j++) {
			k+=(a&1);
			a>>=1;
		}
		if(n_all-k>1) n_all1[i]=k+1;
		else {
			n_all1[i]=n_all;
			required[i]=(1<<ng)-1;
		}
	}
	/* Set up adjacancy matrix */
	for(i=0;i<nobs;i++) {
		o1=state[0][i];
		o2=state[1][i];
		nl1=n_all1[o1];
		nl2=n_all1[o2];
		if(!gflag[o1]) {
			gflag[o1]=1;
			pd1=freq[grps[o1]];
			pd2=gpen[o1];
			if(nl1<n_all) {
				a=required[o1];
				k1=nl1-1;
				pd2[k1]=0.0;
				for(k=j=0;j<n_all;j++,a>>=1) {
					if(a&1) pd2[k++]=pd1[j];
					else pd2[k1]+=pd1[j];
				}
			} else for(j=0;j<n_all;j++) pd2[j]=pd1[j];
		}
		if(!gflag[o2]) {
			gflag[o2]=1;
			pd1=freq[grps[o2]];
			pd2=gpen[o2];
			if(nl2<n_all) {
				a=required[o2];
				k1=nl2-1;
				pd2[k1]=0.0;
				for(k=j=0;j<n_all;j++,a>>=1) {
					if(a&1) pd2[k++]=pd1[j];
					else pd2[k1]+=pd1[j];
				}
			} else for(j=0;j<n_all;j++) pd2[j]=pd1[j];
		}
		if(o1!=o2) {
			p=add_node(o1,o2);
			p1=add_node(o2,o1);
			if(!p->rf) {
				order[o1]++;
				order[o2]++;
				if((rf=free_rf_list[0])) free_rf_list[0]=rf->next;
				else rf=get_new_rfunc(2,n_all);
				rf->inv[0]=o1;
				rf->inv[1]=o2;
				rf->flag=0;
				p1->rf=p->rf=rf;
				pd2=rf->p;
				pd1=obspen[i];
				for(k1=j=0;j<n_all;j++) {
					for(k=0;k<j;k++,k1++) pd2[j*n_all+k]=pd2[k*n_all+j]=pd1[k1];
					pd2[j*(n_all+1)]=pd1[k1++];
				}
			} else {
				pd1=obspen[i];
				pd2=p->rf->p;
				z1=0.0;
				for(k1=j=0;j<n_all;j++) {
					for(k=0;k<j;k++,k1++) {
						z1+=z=pd1[k1];
						pd2[j*n_all+k]*=z;
						pd2[k*n_all+j]*=z;
					}
					z1+=pd2[j*(n_all+1)]*=pd1[k1++];
				}
				if(z1<RESCALE_LIMIT) {
					l+=log(z1);
					z1=1.0/z1;
					for(k=0;k<n_gen;k++) pd2[k]*=z1;
				}
			}
		} else {
			pd1=obspen[i];
			pd2=gpen[o1];
			z1=0.0;
			if(nl1<n_all) {
				a=required[o1];
				for(k1=k=j=0;j<n_all;j++,a>>=1) {
					if(a&1) z1+=(pd2[k++]*=pd1[k1]);
					else z1+=(z=pd1[k1]);
					k1+=j+2;
				}
				pd2[k]*=z;
			} else {
				for(k1=j=0;j<n_all;j++) {
					z1+=(pd2[j]*=pd1[k1]);
					k1+=j+2;
				}
			}
			if(z1<RESCALE_LIMIT) {
				l+=log(z1);
				z1=1.0/z1;
				for(k=0;k<nl1;k++) pd2[k]*=z1;
			}
		}
	}
	/* Have to downcode R-Functions now */
	for(i=0;i<ng;i++) {
		p=adj_mat[i];
		while(p) {
			nl1=n_all1[i];
			if(nl1<n_all) {
				a=required[i];
				for(k1=k=0;k<n_all;k++,a>>=1) {
					if(a&1) trans1[k1++]=k;
					else k2=k;
				}
				trans1[k1]=k2;
			} else for(k=0;k<n_all;k++) trans1[k]=k;
			j=p->x;
			if(j<i) {
				nl2=n_all1[j];
				if(nl1<n_all || nl2<n_all) {
#ifdef TRACE
					if(CHK_TRACE(TRACE_LEVEL_2)) printf("Downcoding RF(%d,%d)->%dx%d\n",j,i,nl2,nl1);
#endif
					a=required[j];
					if(nl2<n_all) {
						for(k1=k=0;k<n_all;k++,a>>=1) {
							if(a&1) trans2[k1++]=k;
							else k2=k;
						}
						trans2[k1]=k2;
					} else for(k=0;k<n_all;k++) trans2[k]=k;
					pd2=p->rf->p;
					memcpy(tmp_xx,pd2,n_gen*sizeof(double));
					for(j=0;j<nl2;j++) {
						j1=trans2[j];
						pd3=tmp_xx+j1*n_all;
						for(k=0;k<nl1;k++) {
							k1=trans1[k];
							*pd2++=pd3[k1];
						}
					}
				}
			} else break;
			p=p->next;
		}
	}
/*	for(i=0;i<ng;i++) {
		p=adj_mat[i];
		while(p) {
			j=p->x;
			if(j>i) {
				pd1=gpen[i];
				pd2=gpen[j];
				pd12=p->rf->p;
				nl1=n_all1[i];
				nl2=n_all1[j];
				z=0.0;
				k=0;
				for(i1=0;i1<nl1;i1++) for(j1=0;j1<nl2;j1++) z+=pd1[i1]*pd2[j1]*pd12[k++];
				printf("%d %d %.15f %.15f\n",i,j,z,log(z));
			}
			p=p->next;
		}
	} */
	for(no=i=0;i<ng;i++) if(gflag[i]) {
		if(!order[i]) {
#ifdef TRACE
		if(CHK_TRACE(TRACE_LEVEL_2)) printf("Peeling %d->*\n",i);
#endif
			pd1=gpen[i];
			nl1=n_all1[i];
			for(z=j=0;j<nl1;j++) z+=pd1[j];
			l+=log(z);
		} else olist[no++]=i;
	}
	while(no) {
/*		break; */
#ifdef TRACE
		if(CHK_TRACE(TRACE_LEVEL_3)) {
			for(i=0;i<no;i++) {
				j=olist[i];
				printf("%d: ",j);
				p=adj_mat[j];
				while(p) {
					printf("%d ",p->x);
					p=p->next;
				}
				printf(" (%d)\n",order[j]);
			}
		}
#endif
		tri_piv=dual_piv=-1;
		for(i=0;i<no;i++) {
			j=olist[i];
			if(order[j]==1) break;
			if(tri_piv<0 && order[j]==2) {
				p=adj_mat[j];
				o1=p->x;
				o2=p->next->x;
				p1=find_node(o1,o2);
				if(p1) tri_piv=j;
				else dual_piv=j;
			}
		}
		if(i<no) {
			/* Peel singly connected pivot */
			o1=j;
			p=adj_mat[o1];
			o2=p->x;
			if(order[o2]>1) {
#ifdef TRACE
				if(CHK_TRACE(TRACE_LEVEL_2)) printf("Peeling %d->%d\n",o1,o2);
#endif
				pd1=gpen[o1];
				pd2=gpen[o2];
				rf=p->rf;
				pd12=rf->p;
				zz=0.0;
				nl1=n_all1[o1];
				nl2=n_all1[o2];
				if(o2<o1) {
					for(k1=j=0;j<nl2;j++) {
						z=0.0;
						for(k=0;k<nl1;k++) z+=pd1[k]*pd12[k1++];
						zz+=pd2[j]*=z;
					}
				} else {
					for(j=0;j<nl2;j++) {
						z=0.0;
						k1=j;
						for(k=0;k<nl1;k++,k1+=nl2) z+=pd1[k]*pd12[k1];
						zz+=pd2[j]*=z;
					}
				}
				if(zz<RESCALE_LIMIT) {
					l+=log(zz);
					zz=1.0/zz;
					for(j=0;j<nl2;j++) pd2[j]*=zz;
				}
				rf->next=free_rf_list[0];
				free_rf_list[0]=rf;
				p->next=freelist;
				freelist=p;
				adj_mat[o1]=0;
				no--;
				for(;i<no;i++) olist[i]=olist[i+1];
				order[o2]--;
				del_node(o1,o2);
			} else {
#ifdef TRACE
				if(CHK_TRACE(TRACE_LEVEL_2)) 	printf("Peeling %d,%d->*\n",o1,o2);
#endif
				pd1=gpen[o1];
				pd2=gpen[o2];
				rf=p->rf;
				pd12=rf->p;
				z1=0.0;
				nl1=n_all1[o1];
				nl2=n_all1[o2];
				if(o1<o2) {
					for(k1=j=0;j<nl1;j++) {
						z=0.0;
						for(k=0;k<nl2;k++) z+=pd2[k]*pd12[k1++];
						z1+=pd1[j]*z;
					}
				} else {
					for(j=0;j<nl2;j++) {
						k1=j;
						z=0.0;
						for(k=0;k<nl1;k++,k1+=n_all) z+=pd1[k]*pd12[k1];
						z1+=pd2[j]*z;
					}
				}
				l+=log(z1);
				del_row(o1);
				del_row(o2);
				rf->next=free_rf_list[0];
				free_rf_list[0]=rf;
				j=i++;
				for(;i<no;i++) {
					k=olist[i];
					if(k!=o2) olist[j++]=k;
				}
				no=j;
			}
		} else if(tri_piv>=0) {
			/* Peel triangular clique */
			o1=tri_piv;
			p=adj_mat[o1];
			o2=p->x;
			o3=p->next->x;
			pd1=gpen[o1];
			pd12=p->rf->p;
			pd13=p->next->rf->p;
			pd23=p1->rf->p;
			nl1=n_all1[o1];
			nl2=n_all1[o2];
			nl3=n_all1[o3];
			if(order[o2]==2 && order[o3]==2) {
#ifdef TRACE
				if(CHK_TRACE(TRACE_LEVEL_2)) printf("Peeling %d,%d,%d->*\n",o1,o2,o3);
#endif
				pd2=gpen[o2];
				pd3=gpen[o3];
				z=0.0;
				for(i=i1=0;i<nl1;i++) {
					z1=pd1[i];
					for(j1=j=0;j<nl2;j++,i1++) {
						z2=pd2[j]*pd12[i1];
						i2=i*nl3;
						for(k=0;k<nl3;k++,j1++,i2++) {
							z+=z1*z2*pd3[k]*pd23[j1]*pd13[i2];
						}
					}
				}
				l+=log(z);
				del_row(o1);
				del_row(o2);
				del_row(o3);
				p->rf->next=p->next->rf;
				p->next->rf->next=p1->rf;
				p1->rf->next=free_rf_list[0];
				free_rf_list[0]=p->rf;
				for(i=j=0;i<no;i++) {
					k=olist[i];
					if(k!=o1 && k!=o2 && k!=o3) olist[j++]=k;
				}
				no=j;
			} else if(order[o2]==2 || order[o3]==2) {
				if(order[o3]==2) {
					k=o3;
					o3=o2;
					o2=k;
					pd3=pd12;
					pd12=pd13;
					pd13=pd3;
					k=nl3;
					nl3=nl2;
					nl2=k;
				}
				pd2=gpen[o2];
				pd3=gpen[o3];
#ifdef TRACE
				if(CHK_TRACE(TRACE_LEVEL_2)) printf("Peeling %d,%d->%d\n",o1,o2,o3); 
#endif
				zz=0.0;
				if(o3<o1) {
					for(k=k1=0;k<nl3;k++) {
						z2=0.0;
						for(i=i1=0;i<nl1;i++,k1++) {
							z=pd1[i]*pd13[k1];
							j1=k*nl2;
							for(j=0;j<nl2;j++,i1++,j1++) z2+=z*pd2[j]*pd12[i1]*pd23[j1];
						}
						zz+=pd3[k]*=z2;
					}
				} else if(o3<o2) {
					for(k=0;k<nl3;k++) {
						z2=0.0;
						k1=k;
						for(i=i1=0;i<nl1;i++,k1+=nl3) {
							z=pd1[i]*pd13[k1];
							j1=k*nl2;
							for(j=0;j<nl2;j++,i1++,j1++) z2+=z*pd2[j]*pd12[i1]*pd23[j1];
						}
						zz+=pd3[k]*=z2;
					}
				} else {
					for(k=0;k<nl3;k++) {
						z2=0.0;
						k1=k;
						for(i=i1=0;i<nl1;i++,k1+=nl3) {
							z=pd1[i]*pd13[k1];
							j1=k;
							for(j=0;j<nl2;j++,i1++,j1+=nl3) z2+=z*pd2[j]*pd12[i1]*pd23[j1];
						}
						zz+=pd3[k]*=z2;
					}
				}
				if(zz<RESCALE_LIMIT) {
					l+=log(zz);
					zz=1.0/zz;
					for(i=0;i<nl3;i++) pd3[i]*=zz;
				}
				del_row(o1);
				del_row(o2);
				p->rf->next=p->next->rf;
				p->next->rf->next=p1->rf;
				p1->rf->next=free_rf_list[0];
				free_rf_list[0]=p->rf;
				order[o3]-=2;
				del_node(o1,o3);
				del_node(o2,o3);
				for(j=i=0;i<no;i++) {
					k=olist[i];
					if(k!=o1 && k!=o2) olist[j++]=k;
				}
				no=j;
			} else {
#ifdef TRACE
				if(CHK_TRACE(TRACE_LEVEL_2)) printf("Peeling %d->%d,%d\n",o1,o2,o3);
#endif
				zz=0.0;
				if(o1<o2) {
					for(j=j1=0;j<nl2;j++) {
						for(k=0;k<nl3;k++,j1++) {
							z2=0.0;
							i1=j;
							k1=k;
							for(i=0;i<nl1;i++,i1+=nl2,k1+=nl3) z2+=pd1[i]*pd12[i1]*pd13[k1];
							zz+=pd23[j1]*=z2;
						}
					}
				} else if(o1<o3) {
					for(j=j1=0;j<nl2;j++) {
						for(k=0;k<nl3;k++,j1++) {
							z2=0.0;
							i1=j*nl1;
							k1=k;
							for(i=0;i<nl1;i++,i1++,k1+=nl3) z2+=pd1[i]*pd12[i1]*pd13[k1];
							zz+=pd23[j1]*=z2;
						}
					}
				} else {
					for(j=j1=0;j<nl2;j++) {
						for(k=k1=0;k<nl3;k++,j1++) {
							z2=0.0;
							i1=j*nl1;
							for(i=0;i<nl1;i++,i1++,k1++) z2+=pd1[i]*pd12[i1]*pd13[k1];
							zz+=pd23[j1]*=z2;
						}
					}
				}
				if(zz<RESCALE_LIMIT) {
					l+=log(zz);
					zz=1.0/zz;
					k=nl2*nl3;
					for(i=0;i<k;i++) pd23[i]*=zz;
				}
				del_row(o1);
				order[o2]--;
				order[o3]--;
				del_node(o1,o3);
				del_node(o1,o2);
				p->rf->next=p->next->rf;
				p->next->rf->next=free_rf_list[0];
				free_rf_list[0]=p->rf;
				for(j=i=0;i<no;i++) {
					k=olist[i];
					if(k!=o1) olist[j++]=k;
				}
				no--;
			}
		} else if(dual_piv>=0) {
			o1=dual_piv;
			p=adj_mat[o1];
			o2=p->x;
			o3=p->next->x;
			pd1=gpen[o1];
			pd12=p->rf->p;
			pd13=p->next->rf->p;
			nl1=n_all1[o1];
			nl2=n_all1[o2];
			nl3=n_all1[o3];
			p1=add_node(o2,o3);
			p2=add_node(o3,o2);
			if((rf=free_rf_list[0])) free_rf_list[0]=rf->next;
			else rf=get_new_rfunc(2,n_all);
			rf->inv[0]=o2;
			rf->inv[1]=o3;
			rf->flag=0;
			p1->rf=p2->rf=rf;
			pd23=p1->rf->p;
#ifdef TRACE
				if(CHK_TRACE(TRACE_LEVEL_2)) printf("Peeling %d->%d,%d\n",o1,o2,o3); 
#endif
			zz=0.0;
			if(o1<o2) {
				for(j=j1=0;j<nl2;j++) {
					for(k=0;k<nl3;k++,j1++) {
						z2=0.0;
						i1=j;
						k1=k;
						for(i=0;i<nl1;i++,i1+=nl2,k1+=nl3) z2+=pd1[i]*pd12[i1]*pd13[k1];
						zz+=pd23[j1]=z2;
					}
				}
			} else if(o1<o3) {
				for(j=j1=0;j<nl2;j++) {
					for(k=0;k<nl3;k++,j1++) {
						z2=0.0;
						i1=j*nl1;
						k1=k;
						for(i=0;i<nl1;i++,i1++,k1+=nl3) z2+=pd1[i]*pd12[i1]*pd13[k1];
						zz+=pd23[j1]=z2;
					}
				}
			} else {
				for(j=j1=0;j<nl2;j++) {
					for(k=k1=0;k<nl3;k++,j1++) {
						z2=0.0;
						i1=j*nl1;
						for(i=0;i<nl1;i++,i1++,k1++) z2+=pd1[i]*pd12[i1]*pd13[k1];
						zz+=pd23[j1]=z2;
					}
				}
			}
			if(zz<RESCALE_LIMIT) {
				l+=log(zz);
				zz=1.0/zz;
				k=nl2*nl3;
				for(i=0;i<k;i++) pd23[i]*=zz;
			}
			del_row(o1);
			p->rf->next=p->next->rf;
			p->next->rf->next=free_rf_list[0];
			free_rf_list[0]=p->rf;
			del_node(o1,o3);
			del_node(o1,o2);
			for(j=i=0;i<no;i++) {
				k=olist[i];
				if(k!=o1) olist[j++]=k;
			}
			no--;
		} else break;
	}
	/* Any genes left ? */
	if(no) {
		/* Yes - we have to do 'complex' gene peeling */
		for(i=0;i<no;i++) {
			first[i]=0;
			j=olist[i];
			deglist[j].abs_list=0;
			deglist[j].deg=-1;
			gflag[j]=0;
			rfp1=0;
			p=adj_mat[j];
			while(p) {
				if((rfp=free_rfp_list)) free_rfp_list=rfp->next;
				else {
					if(!(rfp=malloc(sizeof(struct gp_rfunc_ptr)))) ABT_FUNC(MMsg);
#ifdef TRACE
					if(CHK_TRACE(TRACE_LEVEL_3)) printf("Allocating rfp %p\n",rfp);
#endif
				}
				rfp->rf=p->rf;
				rfp->next=rfp1;
				rfp1=rfp;
				p=p->next;
			}
			rfuncs[j]=rfp1;
		}
		for(i=0;i<no;i++) {
			i1=olist[i];
			if(!gflag[i1]) {
				k=calc_degree(i1);
				if(k>=0) {
					j=0;
					for(k1=1;k1<=k;k1++) if(!gflag[gp_inv[k1]]) j++;
					deglist[i1].deg=j;
					deglist[i1].next=first[j];
					deglist[i1].prev=0;
					if(first[j]) first[j]->prev=deglist+i1;
					first[j]=deglist+i1;
				}
			}
		}
#ifdef TRACE
		if(CHK_TRACE(TRACE_LEVEL_3)) {
			(void)fputs("Degree lists\n",stdout);
			for(i=0;i<no;i++) {
				degp=first[i];
				if(degp) {
					(void)printf("%d:",i);
					while(degp) {
						j=degp->gene;
						(void)printf(" %d",j);
						degp1=degp->abs_list;
						while(degp1) {
							(void)printf("-%d",degp1->gene);
							degp1=degp1->abs_list;
						}
						degp=degp->next;
					}
					(void)fputc('\n',stdout);
				}
			}
		}
#endif
		/* Ok, here goes... */
		for(;;) {
			/* Pick first minimum degree node */
			for(i=0;i<no;i++) if(first[i]) break;
			if(i==no) break;
			j=first[i]->gene;
#ifdef DEBUG
			if(gflag[j]) {
				(void)fprintf(stderr,"Internal error when peeling gene %d from degree list %d\n",j,i);
				ABT_FUNC("Gene already peeled\n");
			}
#endif
			/* Get list of involved genes and to peel genes */
			n_peel=1;
			gp_inv[0]=j;
			gpen1[0]=gpen[j];
			gflag[j]|=4;
			degp=deglist[j].abs_list;
			while(degp) {
				k=degp->gene;
				gflag[k]|=2;
				gpen1[n_peel]=gpen[k];
				gp_inv[n_peel++]=k;
				degp=degp->abs_list;
			}
			n_inv=n_peel;
			p=adj_mat[j];
			while(p) {
				k=p->x;
				if(!(gflag[k]&2)) gp_inv[n_inv++]=k;
				p=p->next;
			}
			for(i=0;i<n_inv;i++) gflag[gp_inv[i]]|=8;
#ifdef TRACE
			if(CHK_TRACE(TRACE_LEVEL_2)) {
				(void)fputs("Peeling",stdout);
				for(i=0;i<n_peel;i++) {
					k=gp_inv[i];
					(void)fputc(i?',':' ',stdout);
					(void)printf("%d",k);
				}
				(void)fputs(" ->",stdout);
				if(n_peel==n_inv) (void)fputs(" *",stdout);
				else {
					for(;i<n_inv;i++) {
						k=gp_inv[i];
						(void)fputc(i>n_inv?',':' ',stdout);
						(void)printf("%d",k);
					}
				}
				fputs("  \n",stdout);
			}
#endif
			/* Get list of rfunctions */
			for(rfl=0,n_rf=i=0;i<n_peel;i++) {
				k=gp_inv[i];
				rfp1=rfuncs[k];
				while(rfp1) {
					rfp2=rfp1->next;
					rf1=rfp1->rf;
					if(!rf1->flag) {
						rf1->next=rfl;
						rfl=rf1;
						rf1->flag=1;
						n_rf++;
					}
#ifdef TRACE
					if(CHK_TRACE(TRACE_LEVEL_3)) (void)printf("[%s:%d] %s() - Adding %p to free_rfp_list\n",__FILE__,__LINE__,__func__,(void *)rfp1);
#endif
					rfp1->next=free_rfp_list;
					free_rfp_list=rfp1;
					rfp1=rfp2;
				}
				rfuncs[k]=0;
			}
			for(;i<n_inv;i++) {
				k=gp_inv[i];
				rfp3=rfuncs+k;
				rfp1=*rfp3;
				while(rfp1) {
					rf1=rfp1->rf;
					inv=rf1->inv;
					j1=rf1->n_inv;
					for(j=0;j<j1;j++) if(!(gflag[inv[j]]&8)) break;
					if(j==j1) {
						if(!rf1->flag) {
							rf1->next=rfl;
							rfl=rf1;
							rf1->flag=1;
							n_rf++;
						}
#ifdef TRACE
						if(CHK_TRACE(TRACE_LEVEL_3)) (void)printf("[%s:%d] %s() - Adding %p to free_rfp_list\n",__FILE__,__LINE__,__func__,(void *)rfp1);
#endif
						*rfp3=rfp1->next;
						rfp2=*rfp3;
						rfp1->next=free_rfp_list;
						free_rfp_list=rfp1;
					} else {
						rfp3=&rfp1->next;
						rfp2=*rfp3;
					}
					rfp1=rfp2;
				}
			}
			for(i=n_peel;i<n_inv;i++) gflag[gp_inv[i]]&=~8;
#ifdef TRACE
			if(CHK_TRACE(TRACE_LEVEL_2)) {
				rf1=rfl;
				if(rf1) {
					fputs("RF: ",stdout);
					while(rf1) {
						(void)fputc(' ',stdout);
						for(i=0;i<rf1->n_inv;i++) {
							(void)fputc(i?',':'(',stdout);
							(void)printf("%d",rf1->inv[i]);
						}
						(void)fputc(')',stdout);
						rf1=rf1->next;
					}
					(void)fputc('\n',stdout);
				}
			}
#endif
			/* Convert genes in rfunctions to references in gp_inv */
			for(i=0;i<n_inv;i++) gstate[gp_inv[i]]=i;
			rf1=rfl;
			while(rf1) {
				inv=rf1->inv;
				for(i=0;i<rf1->n_inv;i++) {
					j=inv[i];
					inv[i]=gstate[j];
				}
				rf1=rf1->next;
			}
			/* Setup output rfunction */
			n_out=n_inv-n_peel;
			if(n_out) {
				if(n_out>1) {
					if((rf=free_rf_list[n_out-2])) free_rf_list[n_out-2]=rf->next;
					else rf=get_new_rfunc(n_out,n_all);
					rf->next=0;
					rf->flag=0;
					xx=rf->p;
					nc=1;
					for(j=n_out,i=n_peel;i<n_inv;i++) {
						k1=gp_inv[i];
						nc*=n_all1[k1];
						rf->inv[--j]=k1;
						if((rfp1=free_rfp_list)) free_rfp_list=rfp1->next;
						else {
							if(!(rfp1=malloc(sizeof(struct gp_rfunc_ptr)))) ABT_FUNC(MMsg);
#ifdef TRACE
							if(CHK_TRACE(TRACE_LEVEL_3)) printf("Allocating rfp1 %p\n",rfp);
#endif
						}
						rfp1->next=rfuncs[k1];
						rfp1->rf=rf;
						rfuncs[k1]=rfp1;
						/* Fill in new adjacencies */
						for(i1=n_peel;i1<i;i1++) {
							k2=gp_inv[i1];
							p=add_node(k1,k2);
							if(!p->rf) (void)add_node(k2,k1);
						}
					}
				} else {
					nc=n_all1[gp_inv[n_peel]];
					xx=tmp_xx;
				}
				for(i=0;i<nc;i++) xx[i]=0.0;
			}
			/* Check we have enough space for rfunction indices */
			if(n_rf>rf_idx_n) {
				rf_idx_n=n_rf;
				if(!(fact[0]=realloc(fact[0],sizeof(int)*(2*rf_idx_n*ng)))) ABT_FUNC(MMsg);
				if(!(rf_ptr=realloc(rf_ptr,sizeof(void *)*rf_idx_n))) ABT_FUNC(MMsg);
				rfx[0]=fact[0]+rf_idx_n*ng;
				for(k1=1;k1<ng;k1++) {
					fact[k1]=fact[k1-1]+rf_idx_n;
					rfx[k1]=rfx[k1-1]+rf_idx_n;
				}
			}
			/* Perform peel operation */
			/* Initialize Gray code generator */
			for(i=0;i<n_inv;i++) {
				tt[i]=i;
				cc1[i]=n_all1[gp_inv[i]];
				cc[i]=cc1[i]-1;
				gstate[i]=0;
			}
			/* Initialize R-Function indices */
			k=0;
			rf1=rfl;
			for(k1=0;k1<n_inv;k1++) rfxn[k1]=0;
			while(rf1) {
				inv=rf1->inv;
				for(k2=1,k1=rf1->n_inv-1;k1>=0;k1--) {
					i2=inv[k1];
					fact[i2][k]=k2;
					k2*=cc1[i2];
					rfx[i2][rfxn[i2]++]=k;
				}
				rf_ptr[k++]=rf1->p;
				rf1=rf1->next;
			}
			/* Initialize output indices */
			for(k=1,i=n_peel;i<n_inv;i++) {
				ofact[i]=k;
				k*=cc1[i];
			}
			/* Initialize individual gene probs. */
			/* Loop through all combinations */
			i1=0;
			z1=0.0;
			if(n_out) {
				for(;;) {
					/* Get probs. from individual peeled gene functions and R-Functions */
					z=1.0;
					tp1=gpen1;
					for(i=0;i<n_peel;i++) z*=(**(tp1++));
					tp1=rf_ptr;
					for(i=0;i<n_rf;i++) z*=(**(tp1++));
					/* Store in output function */
					xx[i1]+=z;
					z1+=z;
					/* Get next combination */
					k=tt[0];
					tt[0]=0;
					if(!cc[k]) break;
					if(!(--cc[k])) {
						if(k<n_inv-1) {
							tt[k]=tt[k+1];
							tt[k+1]=k+1;
							cc[k]=cc1[k]-1;
						}
					}
					j=rfxn[k];
					inv=rfx[k];
					inv1=fact[k];
					gstate[k]++;
					if(gstate[k]<cc1[k]) {
						if(k<n_peel) gpen1[k]++;
						else i1+=ofact[k];
						for(k2=0;k2<j;k2++) {
							j1=*inv++;
							rf_ptr[j1]+=inv1[j1];
						}
					} else {
						gstate[k]=0;
						i=cc1[k]-1;
						if(k<n_peel) gpen1[k]-=i;
						else i1-=i*ofact[k];
						for(k2=0;k2<j;k2++) {
							j1=*inv++;
							rf_ptr[j1]-=inv1[j1]*i;
						}
					}
				}
				if(z1<=0.0) {
					clean_up(ng,rfl);
					return 0.0;
				}
				if(n_out==1) {
					k=gp_inv[n_inv-1];
					pd1=gpen[k];
					z1=0.0;
					for(i=0;i<nc;i++) z1+=(pd1[i]*=xx[i]);
					xx=pd1;
				}
				if(z1<RESCALE_LIMIT) {
					l+=log(z1);
					z1=1.0/z1;
					for(i=0;i<nc;i++) xx[i]*=z1;
				}
				for(i=n_peel;i<n_inv;i++) {
					k=gp_inv[i];
					rfp3=rfuncs+k;
					rfp1=*rfp3;
					while(rfp1) {
						rfp2=rfp1->next;
						if(rfp1->rf->flag) {
#ifdef TRACE
							if(CHK_TRACE(TRACE_LEVEL_3)) (void)printf("[%s:%d] %s() - Adding %p to free_rfp_list\n",__FILE__,__LINE__,__func__,(void *)rfp1);
#endif
							*rfp3=rfp1->next;
							rfp1->next=free_rfp_list;
							free_rfp_list=rfp1;
						} else rfp3=&rfp1->next;
						rfp1=rfp2;
					}
				}
			} else {
				for(;;) {
					/* Get probs. from individual peeled gene functions and R-Functions */
					z=1.0;
					tp1=gpen1;
					for(i=0;i<n_peel;i++) z*=(**(tp1++));
					tp1=rf_ptr;
					for(i=0;i<n_rf;i++) z*=(**(tp1++));
					z1+=z;
					/* Get next combination */
					k=tt[0];
					tt[0]=0;
					if(!cc[k]) break;
					if(!(--cc[k])) {
						if(k<n_inv-1) {
							tt[k]=tt[k+1];
							tt[k+1]=k+1;
							cc[k]=cc1[k]-1;
						}
					}
					k1=gstate[k]=(gstate[k]+1)%cc1[k];
					j=rfxn[k];
					inv=rfx[k];
					inv1=fact[k];
					if(k1) {
						gpen1[k]++;
						for(k2=0;k2<j;k2++) {
							j1=inv[k2];
							rf_ptr[j1]+=inv1[j1];
						}
					} else {
						i=cc1[k]-1;
						gpen1[k]-=i;
						for(k2=0;k2<j;k2++) {
							j1=inv[k2];
							rf_ptr[j1]-=inv1[j1]*i;
						}
					}
				}
				if(z1<=0.0) {
					clean_up(ng,rfl);
					return 0.0;
				}
				l+=log(z1);
			}
			/* Free used rfunctions */
			free_rfunc_list(rfl);
			/* Remove peeled nodes from degree list */
			for(i=0;i<n_peel;i++) {
				no--;
				k=gp_inv[i];
				if(!(gflag[k]&2)) {
#ifdef TRACE
					if(CHK_TRACE(TRACE_LEVEL_2)) (void)printf("Removing %d from degree list %d\n",k,deglist[k].deg);
#endif
					if(deglist[k].prev) deglist[k].prev->next=deglist[k].next;
					else {
						k1=deglist[k].deg;
						first[k1]=deglist[k].next;
					}
					if(deglist[k].next) deglist[k].next->prev=deglist[k].prev;
				}
				del_row1(k);
				/* Mark gene as peeled */
				gflag[k]|=2;
			}
			/* Recalculate degree for output nodes */
			for(i=n_peel;i<n_inv;i++) {
				k=gp_inv[i];
				if(!gflag[k]) {
					k1=calc_degree(k);
					if(k1>=0) {
						k2=0;
						for(j=1;j<=k1;j++) if(!gflag[gp_inv[j]]) k2++;
					} else k2=k1;
					k1=deglist[k].deg;
					if(k1!=k2) {
#ifdef TRACE
						if(CHK_TRACE(TRACE_LEVEL_2)) (void)printf("Moving %d from degree list %d to %d\n",k,k1,k2);
#endif
						if(deglist[k].prev) deglist[k].prev->next=deglist[k].next;
						else first[k1]=deglist[k].next;
						if(deglist[k].next) deglist[k].next->prev=deglist[k].prev;
						deglist[k].next=first[k2];
						deglist[k].prev=0;
						if(first[k2]) first[k2]->prev=deglist+k;
						first[k2]=deglist+k;
						deglist[k].deg=k2;
					}
				}
			}
		}
	}
/*	printf(" %d %d %d %d %d ",count[0],count[1],count[2],count[3],count[4]); */
/*	printf("l=%g\n",l);
	exit(0); */
	return l;
}
