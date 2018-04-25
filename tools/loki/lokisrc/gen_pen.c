/****************************************************************************
*                                                                          *
*     Loki - Programs for genetic analysis of complex traits using MCMC    *
*                                                                          *
*             Simon Heath - MSKCC - June 2001                              *
*                                                                          *
* gen_pen.c:                                                               *
*                                                                          *
* Calculate p(Y|S) for general penetrance function                         *
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
#include "seg_pen.h"
#include "gen_pen.h"
#include "handle_res.h"
#include "lk_malloc.h"

static struct gp_rfnode **rflist,*free_rfn_list;
static struct gp_rfunc_ptr **rfuncs,*free_rfp_list;
static struct gp_peel_op *free_ops_list[MAX_GP_RF_GENES];
static struct gp_rfunc *free_rf_list[MAX_GP_RF_GENES],**autozyg_list;
static struct deg_list *deglist,**first;
static int *gp_inv,*gp_inv1,*gp_flag,*pflag,*rf_idx,
**fact,*ofact,rf_idx_n,*rfxn,**rfx,
*cnv_gene,gp_xx_size,mem_list_size,mem_list_ptr;
static double *gp_tmp_p,*gp_z,*gp_xx;
static void **mem_list;

static void *realloc_and_remember(void *p,size_t size)
{
	int i;
	
	for(i=0;i<mem_list_ptr;i++) if(mem_list[i]==p) break;
	if(i<mem_list_ptr) {
		p=lk_realloc(p,size);
		mem_list[i]=p;
	} else p=0;
	return p;
}

static void *malloc_and_remember(size_t size)
{
	void *p;
	
	if(mem_list_ptr==mem_list_size) {
		if(mem_list_size) {
			mem_list_size*=1.5;
			mem_list=lk_realloc(mem_list,sizeof(void *)*mem_list_size);
		} else {
			mem_list_size=16;
			mem_list=lk_malloc(sizeof(void *)*mem_list_size);
		}
	}
	p=lk_malloc(size);
	mem_list[mem_list_ptr++]=p;
	return p;
}

static void free_rfunc(struct gp_rfunc *rf)
{
#ifdef TRACE_PEEL
	if(CHK_PEEL(TRACE_LEVEL_3)) (void)printf("In free_rfunc(%p)\n",(void *)rf);
#endif
	free(rf->inv);
	free(rf->p);
	free(rf); 
}

void free_gen_pen(void) 
{
	int i;
	struct gp_rfunc_ptr *rfp,*rfp1;
	struct gp_rfnode *rfn,*rfn1;
	struct gp_rfunc *rf,*rf1;
	struct gp_peel_op *ops,*ops1;
	
#ifdef TRACE_PEEL
	if(CHK_PEEL(TRACE_LEVEL_1)) (void)printf("In free_gen_pen()\n");
#endif
	rfp=free_rfp_list;
	while(rfp) {
		rfp1=rfp->next;
		free(rfp);
		rfp=rfp1;
	}
	free_rfp_list=0;
	rfn=free_rfn_list;
	while(rfn) {
		rfn1=rfn->next;
		free(rfn);
		rfn=rfn1;
	}
	free_rfn_list=0;
	if(mem_list) {
		for(i=0;i<mem_list_ptr;i++) if(mem_list[i]) free(mem_list[i]);
		free(mem_list);
		mem_list=0;
		mem_list_size=mem_list_ptr=0;
	}
	rf_idx_n=0;
	for(i=0;i<MAX_GP_RF_GENES;i++) {
		rf=free_rf_list[i];
		while(rf) {
			rf1=rf->next;
			free_rfunc(rf);
			rf=rf1;
		}
		free_rf_list[i]=0;
		ops=free_ops_list[i];
		while(ops) {
			ops1=ops->next;
			free(ops->inv);
			free(ops);
			ops=ops1;
		}
		free_ops_list[i]=0;
	}
}

void alloc_gen_pen(int ngenes)
{
	int i,n_all;
	
	prt_peel_trace(TRACE_LEVEL_1,"In alloc_gen_pen(%d)\n",ngenes);
	rflist=malloc_and_remember(sizeof(void *)*ngenes);
	autozyg_list=malloc_and_remember(sizeof(void *)*ngenes);
	rfuncs=malloc_and_remember(sizeof(void *)*ngenes);
	first=malloc_and_remember(sizeof(void *)*ngenes);
	fact=malloc_and_remember(sizeof(void *)*2*ngenes);
	rfx=fact+ngenes;
	deglist=malloc_and_remember(sizeof(struct deg_list)*ngenes);
	for(i=0;i<ngenes;i++) deglist[i].gene=i;
	gp_inv=malloc_and_remember(sizeof(int)*7*ngenes);
	gp_inv1=gp_inv+ngenes;
	gp_flag=gp_inv1+ngenes;
	pflag=gp_flag+ngenes;
	ofact=pflag+ngenes;
	rfxn=ofact+ngenes;
	cnv_gene=rfxn+ngenes;
	n_all=2;
	gp_tmp_p=malloc_and_remember(sizeof(double)*(ngenes+n_all*n_all));
	gp_z=gp_tmp_p+n_all*n_all;
	gp_xx_size=8;
	gp_xx=malloc_and_remember(sizeof(double)*gp_xx_size);
}

static struct gp_rfunc *get_new_rfunc(int n,int n_all)
{
	int nc;
	struct gp_rfunc *rf;
	
	rf=lk_malloc(sizeof(struct gp_rfunc));
	rf->n_inv=n;
	nc=(int)(.5+exp(log((double)n_all)*(double)n));
	rf->nc=nc;
	rf->p=lk_malloc(sizeof(double)*nc);
	rf->inv=lk_malloc(sizeof(int)*n);
	prt_peel_trace(TRACE_LEVEL_3,"In get_new_rfunc(%d,%d): allocating new rfunc %p\n",n,n_all,(void *)rf);
	return rf;
}

static struct gp_peel_op *get_new_peel_op(int n)
{
	struct gp_peel_op *ops;
	
	assert(n<=MAX_GP_RF_GENES && n>0);
	ops=free_ops_list[n-1];
	if(ops) {
		free_ops_list[n-1]=ops->next;
		prt_peel_trace(TRACE_LEVEL_3,"In get_new_peel_op(%d): getting peel_op %p from free list\n",n,(void *)ops);
	} else {
		ops=lk_malloc(sizeof(struct gp_peel_op));
		ops->n_inv=n;
		ops->inv=lk_malloc(sizeof(int)*n*2);
		ops->pflag=ops->inv+n;
		prt_peel_trace(TRACE_LEVEL_3,"In get_new_peel_op(%d): allocating new peel_op %p\n",n,(void *)ops);
	}
	ops->next=0;
	return ops;
}

static struct gp_rfnode *find_rfnode(int x,int y)
{
	struct gp_rfnode *p,**pp;
	
	pp=rflist+x;
	p=*pp;
	while(p) {
		if(p->y<=y) break;
		pp=&p->next;
		p=*pp;
	}
	if(!p || p->y!=y) {
		if((p=free_rfn_list)) free_rfn_list=p->next;
		else p=lk_malloc(sizeof(struct gp_rfnode));
		p->next=*pp;
		p->y=y;
		p->rf=0;
		*pp=p;
	}
	return p;
}

static void delete_rfnode(int x,int y)
{
	struct gp_rfnode *p,**pp;
	
	pp=rflist+x;
	p=*pp;
	while(p) {
		if(p->y<=y) break;
		pp=&p->next;
		p=*pp;
	}
	if(!p || p->y!=y) {
		(void)fprintf(stderr,"node %d,%d not found",x,y);
		ABT_FUNC("aborting\n");
	}
	*pp=p->next;
	p->next=free_rfn_list;
	free_rfn_list=p;
}

static int get_inv_list(int *inv,int i)
{
	int k=1;
	struct gp_rfnode *rfn;
	  
	inv[0]=i;
	rfn=rflist[i];
	while(rfn) {
		inv[k++]=rfn->y;
		rfn=rfn->next;
	}
	return k;
}

static void del_rfrow(int i)
{
	struct gp_rfnode *rfn;
	
	rfn=rflist[i];
	if(rfn) {
		while(rfn->next) {
			delete_rfnode(rfn->y,i);
			rfn=rfn->next;
		}
		delete_rfnode(rfn->y,i);
		rfn->next=free_rfn_list;
		free_rfn_list=rflist[i];
		rflist[i]=0;
#ifdef TRACE_PEEL
		if(CHK_PEEL(TRACE_LEVEL_3)) (void)printf("In del_rfrow(): Adding row %d to free_rfn_list\n",i);
#endif
	}
}

static int qfunc(const void *s1,const void *s2) 
{
	int i1,i2;
	
	i1=*(const int *)s1;
	i2=*(const int *)s2;
	if(i1<i2) return -1;
	if(i2<i1) return 1;
	return 0;
}

static int calc_degree(int i) 
{
	int i1,i2,j,j2,k=0,k1,k2;
	struct deg_list *p;
	struct gp_rfnode *rfn;
	
	if(rflist[i]) {
		k2=k=get_inv_list(gp_inv,i);
		gp_inv[k]=-1;
		for(j=1;j<k;j++) {
			i1=gp_inv[j];
			if(!gp_flag[i1]) {
				rfn=rflist[i1];
				if(rfn->next) {
					j2=1;
					while(rfn) {
						k1=rfn->y;
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
					else k1=(!rfn && j2==k)?1:0;
				} else {
					k1=k==2?1:0;
				}
				if(k1) {
#ifdef TRACE_PEEL
					if(CHK_PEEL(TRACE_LEVEL_2)) (void)printf("Absorbing gene %d into %d\n",i1,i);
#endif
					gp_flag[i1]|=1;
					p=deglist[i1].abs_list;
					if(p) {
						while(p->abs_list) p=p->abs_list;
						p->abs_list=deglist[i].abs_list;
					} else deglist[i1].abs_list=deglist[i].abs_list;
					deglist[i].abs_list=deglist+i1;
					if(deglist[i1].deg>=0) {
						i2=deglist[i1].deg;
#ifdef TRACE_PEEL
						if(CHK_PEEL(TRACE_LEVEL_2)) (void)printf("Removing %d from degree list %d (b)\n",i1,i2);
#endif
						if(deglist[i1].prev) deglist[i1].prev->next=deglist[i1].next;
						else first[i2]=deglist[i1].next;
						if(deglist[i1].next) deglist[i1].next->prev=deglist[i1].prev;
						deglist[i1].deg=-1;
					}
					rfn=rflist[i1];
					while(rfn) {
						k1=rfn->y;
						if(k1!=i && !gp_flag[k1]) {
							j2=deglist[k1].deg;
							if(j2>=0) {
#ifdef TRACE_PEEL
								if(CHK_PEEL(TRACE_LEVEL_2)) (void)printf("Moving %d from degree list %d to %d (b)\n",k1,j2,j2-1);
#endif
								if(!j2) ABT_FUNC("Internal error - degree shouldn't be zero\n");
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
	} else if(autozyg_list[i]) k=1;
	return k-1;
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
#ifdef TRACE_PEEL
			if(CHK_PEEL(TRACE_LEVEL_3)) (void)printf("[%s:%d] %s() - Adding %p to free_rfp_list\n",__FILE__,__LINE__,__func__,(void *)rfp1);
#endif
			rfp1->next=free_rfp_list;
			free_rfp_list=rfp1;
			rfp1=rfp2;
		}
		del_rfrow(i);
	}
	while(rfl) {
		rf=rfl->next;
		i=rfl->n_inv-1;
#ifdef TRACE_PEEL
		if(CHK_PEEL(TRACE_LEVEL_3)) (void)printf("[%s:%d] %s() - Adding %p to free_rf_list[%d]\n",__FILE__,__LINE__,__func__,(void *)rfl,i);
#endif
		rfl->next=free_rf_list[i];
		free_rf_list[i]=rfl;
		rfl=rf;
	}
}

static void free_rfunc_list(struct gp_rfunc *rfl)
{
	int k;
	struct gp_rfunc *rf;
	
	while(rfl) {
		rf=rfl->next;
		k=rfl->n_inv-1;
#ifdef TRACE_PEEL
		if(CHK_PEEL(TRACE_LEVEL_3)) (void)printf("[%s:%d] %s() - Adding %p to free_rf_list[%d]\n",__FILE__,__LINE__,__func__,(void *)rfl,k);
#endif
		rfl->next=free_rf_list[k];
		free_rf_list[k]=rfl;
		rfl=rf;
	}
}


double gen_pen(struct Locus *locus,int c,int *err,int flag,const struct loki *loki)
{ 
	int i,i1,i2,j,k,k1,k2,k3,a1,a2,n_all,n_peel,n_inv,n_out,*inv,*inv1,nc=0,*gt,**seg,ngenes1,si;
	int n_rf=0,*ff,*idx,out_idx,mtype,cpsize,ngenes,**genes,*grp=0,*pflag1,ids,idd,rec,nrec,n_genetic_groups;
	double *xx=0,like=0.0,z1,z,z3,**freq,*freq1=0,*eff,*count1,*rflp;
	struct gp_rfnode *rfnode,*rfnode1,**rfnodep;
	struct gp_rfunc *rf,*rf1,**rf2,*rfl;
	struct gp_rfunc_ptr *rfp,*rfp1,*rfp2,**rfp3;
	struct gp_peel_op *ops,*ops1,*opslist=0;
	struct deg_list *degp;
	struct Id_Record *id,*id_array;
	trait_pen_func *pen=0;
	
#ifdef TRACE_PEEL
	struct deg_list *degp1;
	if(CHK_PEEL(TRACE_LEVEL_1)) (void)printf("\nIn %s(%d,%d,%p,%d)\n",__func__,loc,c,(void *)err,flag);
#endif
	mtype=loki->models->models[0].var.type;
	if(!loki->data->n_nonid_records && !loki->models->use_student_t && (!(mtype&ST_CENSORED) || loki->models->censor_mode)) pen=s_penetrance1;
	else pen=s_penetrance;
	n_genetic_groups=loki->pedigree->n_genetic_groups;
	id_array=loki->pedigree->id_array;
	si=loki->params.si_mode;
	freq=locus->freq;
	gt=locus->gt;
	n_all=locus->n_alleles;
	if(n_all!=2) ABT_FUNC("Only for di-allelic QTLs\n");
	ngenes=loki->pedigree->comp_ngenes[c];
	genes=locus->genes;
	eff=locus->eff[0];
	*err=0;
	cpsize=loki->pedigree->comp_size[c];
	if(n_genetic_groups>1) grp=loki->peel->tl_group[c];
	else freq1=freq[0];
	if(loki->pedigree->singleton_flag[c]) {
#ifdef TRACE_PEEL
		if(CHK_PEEL(TRACE_LEVEL_1)) (void)fputs("\nHandling singletons\n",stdout);
#endif
		i=loki->pedigree->comp_start[c];
		if(n_genetic_groups==1) {
			xx=gp_xx+4;
			for(j=0;j<2;j++) {
				z1=freq1[j];
				for(k=j+1;k<2;k++) {
					z=z1*freq1[k];
					xx[j*2+k]=xx[k*2+j]=z;
				}
				xx[j*2+j]=z1*z1;
			}
			count1=loki->peel->seg_count[0];
			while(cpsize--) {
				if(id_array[i].res[0]) {
					for(j=0;j<4;j++) gp_xx[j]=xx[j];
					pen(gp_xx,i,locus,loki);
					for(z=0.0,j=0;j<4;j++) z+=gp_xx[j];
					like+=log(z);
					if(flag) {
						z1=safe_ranf()*z;
						z=0.0;
						for(j=0;j<4;j++) if(gp_xx[j]>0.0) {
							z+=gp_xx[j];
							if(z1<=z) break;
						}
							k1=1+(j&1);
						k2=1+(j>>1);
						k=k1>k2?k1*(k1-1)/2+k2:k2*(k2-1)/2+k1;
						assert(k>=0 && k<4);
						if(flag&4) {
							count1[k1-1]+=1.0;
							count1[k2-1]+=1.0;
						}
						nrec=id_array[i].n_rec;
						k1=(locus->flag&LOCUS_SAMPLED)?gt[i]:1;
						if(k1!=k) {
							z=(k1>1)?eff[k1-2]:0.0;
							if(k>1) z-=eff[k-2];
							for(rec=0;rec<nrec;rec++) {
								id_array[i].res[0][rec]+=z;
							}
						}
						gt[i]=k;
					}
				} else if(flag) {
					z1=safe_ranf();
					z=0.0;
					for(j=0;j<4;j++) {
						if(xx[j]>0.0) {
							z+=xx[j];
							if(z1<=z) break;
						}
					}
					k1=1+(j&1);
					k2=1+(j>>1);
					k=k1>k2?k1*(k1-1)/2+k2:k2*(k2-1)/2+k1;
					if(flag&4) {
						count1[k1-1]+=1.0;
						count1[k2-1]+=1.0;
					}
					assert(k>=0 && k<4);
					gt[i]=k;
				}
				i++;
			}
		} else {
			while(cpsize--) {
				if(id_array[i].res[0] || flag) {
					a1=genes[0][i]-1;
					k=grp[a1];
					freq1=freq[k];
					for(j=0;j<2;j++) {
						z1=freq1[j];
						for(k=j+1;k<2;k++) {
							z=z1*freq1[k];
							gp_xx[j*2+k]=gp_xx[k*2+j]=z;
						}
						gp_xx[j*2+j]=z1*z1;
					}
					if(id_array[i].res[0]) {
						pen(gp_xx,i,locus,loki);
						for(z=0.0,j=0;j<4;j++) z+=gp_xx[j];
						like+=log(z);
					} else z=1.0;
					if(flag) {
						z1=safe_ranf()*z;
						z=0.0;
						for(j=0;j<4;j++) {
							if(gp_xx[j]>0.0) {
								z+=gp_xx[j];
								if(z1<=z) break;
							}
						}
						k1=1+(j&1);
						k2=1+(j>>1);
						k=k1>k2?k1*(k1-1)/2+k2:k2*(k2-1)/2+k1;
						assert(k>=0 && k<4);
						if(id_array[i].res[0]) {
							nrec=id_array[i].n_rec;
							k1=(locus->flag&LOCUS_SAMPLED)?gt[i]:1;
							if(k1!=k) {
								z=(k1>1)?eff[k1-2]:0.0;
								if(k>1) z-=eff[k-2];
								for(rec=0;rec<nrec;rec++) {
									id_array[i].res[0][rec]+=z;
								}
							}
						}
						gt[i]=k;
						if(flag&4) {
							count1=loki->peel->seg_count[grp[a1]];
							count1[k1-1]+=1.0;
							count1[k2-1]+=1.0;
						}
					}
				}
				i++;				
			}
		}
#ifdef TRACE_PEEL
		if(CHK_PEEL(TRACE_LEVEL_1)) (void)printf("Returning like = %g, p = %g\n",like,exp(like));
#endif
		return like;
	}
	/* Build up initial (penetrance) R-Functions */
	for(i=0;i<ngenes;i++) {
		rflist[i]=0;
		rfuncs[i]=0;
		autozyg_list[i]=0;
	}
	i=loki->pedigree->comp_start[c];
	while(cpsize--) {
		if(id_array[i].res[0]) {
			a1=genes[0][i]-1;
			a2=genes[1][i]-1;
			assert(a1>=0 && a2>=0 && a1<ngenes && a2<ngenes);
#ifdef TRACE_PEEL
			if(CHK_PEEL(TRACE_LEVEL_2)) {
				(void)fputs("Ind ",stdout);
				print_orig_id(stdout,i+1);
				(void)printf(": genes %d %d\n",a1,a2);
			}
#endif
			if(a1==a2) {
				rf=autozyg_list[a1];
				if(!rf) {
					/* Allocate rfunction when a1==a2 */
					if((rf=free_rf_list[0])) free_rf_list[0]=rf->next;
					else rf=get_new_rfunc(1,2);
					rf->next=0;
					rf->flag=0;
					rf->inv[0]=a1;
					if((rfp=free_rfp_list)) free_rfp_list=rfp->next;
					else rfp=lk_malloc(sizeof(struct gp_rfunc_ptr));
					rfp->next=rfuncs[a1];
					rfp->rf=rf;
					rfuncs[a1]=rfp;
					xx=rf->p;
					for(j=0;j<2;j++) xx[j]=1.0;
					autozyg_list[a1]=rf;
				}
			} else {
				if(a1>a2) {
					k=a1;
					a1=a2;
					a2=k;
				}
				rfnode=find_rfnode(a1,a2);
				/* Is this a new combination? */
				if(!rfnode->rf) {
					/* Allocate rfunction */
					if((rf=free_rf_list[1])) free_rf_list[1]=rf->next;
					else rf=get_new_rfunc(2,2);
					rf->next=0;
					rf->flag=0;
					rf->inv[0]=a1;
					rf->inv[1]=a2;
					rfnode->rf=rf;
					if((rfp=free_rfp_list)) free_rfp_list=rfp->next;
					else rfp=lk_malloc(sizeof(struct gp_rfunc_ptr));
					rfp->next=rfuncs[a1];
					rfp->rf=rf;
					rfuncs[a1]=rfp;
					/* Add in reference to symmetrical node */
					if((rfp=free_rfp_list)) free_rfp_list=rfp->next;
					else rfp=lk_malloc(sizeof(struct gp_rfunc_ptr));
					rfp->next=rfuncs[a2];
					rfp->rf=rf;
					rfuncs[a2]=rfp;
					rfnodep=rflist+a2;
					if((rfnode=free_rfn_list)) free_rfn_list=rfnode->next;
					else rfnode=lk_malloc(sizeof(struct gp_rfnode));
					while((rfnode1=*rfnodep)) {
						if(rfnode1->y<a1) {
							rfnode->next=rfnode1;
							*rfnodep=rfnode;
							break;
						}
						rfnodep=&rfnode1->next;
					}
					if(!rfnode1) {
						*rfnodep=rfnode;
						rfnode->next=0;
					}
					rfnode->y=a1;
					rfnode->rf=rf;
					xx=rf->p;
					for(j=0;j<4;j++) xx[j]=1.0;
				} else rf=rfnode->rf;
			}
			xx=rf->p;
			/* Add in penetrance */
			if(a1!=a2) {
				pen(xx,i,locus,loki);
				z=xx[0]+xx[1]+xx[2]+xx[3];
				like+=log(z);
				z=1.0/z;
				for(j=0;j<4;j++) xx[j]*=z;
#ifdef TRACE_PEEL
				if(CHK_PEEL(TRACE_LEVEL_4)) {
					for(k1=j=0;j<2;j++) {
						for(k=0;k<2;k++) {
							(void)printf("%d %d : %g\n",j,k,xx[k1++]);
						}
					}
					printf("z=%g, like=%g\n",1.0/z,like);
				}
#endif
			} else {
				gp_tmp_p[0]=gp_tmp_p[3]=1.0;
				gp_tmp_p[1]=gp_tmp_p[2]=0.0;
				pen(gp_tmp_p,i,locus,loki);
				xx[0]=gp_tmp_p[0];
				xx[1]=gp_tmp_p[3];
				prt_peel_trace(TRACE_LEVEL_4,"%d : %g\n%d : %g\n",0,xx[0],1,xx[1]);
			}
		}
		i++;
	}
#ifdef TRACE_PEEL
	if(CHK_PEEL(TRACE_LEVEL_3)) {
		(void)fputs("Node adjacancies\n",stdout);
		printf("ngenes=%d\n",ngenes);
		for(i=0;i<ngenes;i++) {
			rfnode=rflist[i];
			if(rfnode) {
				(void)printf("%d:",i);
				while(rfnode) {
					(void)printf(" %d",rfnode->y);
					assert(rflist[rfnode->y] && rfnode->y<ngenes);
					rfnode=rfnode->next;
				}
				(void)fputc('\n',stdout);
			}
		}
	}
#endif
	/* Find peeling sequence using the minimum external degree approach */
	/* Find degree of each node and put into linked list */
	for(i=0;i<ngenes;i++) {
		first[i]=0;
		gp_flag[i]=0;
		deglist[i].abs_list=0;
		deglist[i].deg=-1;
	}
	for(i=0;i<ngenes;i++) {
		if(!gp_flag[i]) {
			k=calc_degree(i);
			if(k>=0) {
				j=0;
				for(k1=1;k1<=k;k1++) if(!gp_flag[gp_inv[k1]]) j++;
				deglist[i].deg=j;
				deglist[i].next=first[j];
				deglist[i].prev=0;
				if(first[j]) first[j]->prev=deglist+i;
				first[j]=deglist+i;
			}
		}
	}
#ifdef TRACE_PEEL
	if(CHK_PEEL(TRACE_LEVEL_3)) {
		(void)fputs("Degree lists\n",stdout);
		for(i=0;i<ngenes;i++) {
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
	/* Main peeling loop */
	ngenes1=ngenes;
	for(;;) {
		/* Pick a node of minimum degree */
		for(i=0;i<ngenes1;i++) if(first[i]) break;
		if(i==ngenes1) break;
		j=first[i]->gene;
		if(gp_flag[j]) {
			(void)fprintf(stderr,"Internal error when peeling gene %d from degree list %d\n",j,i);
			ABT_FUNC("Gene already peeled\n");
		}
		/* Get list of genes involved */
		n_inv=get_inv_list(gp_inv,j);
		/* Find out who we can peel */
		n_peel=1;
		gp_inv1[0]=j;
		gp_flag[j]|=4;
		degp=deglist[j].abs_list;
		while(degp) {
			j=degp->gene;
			gp_flag[j]|=2;
			gp_inv1[n_peel++]=j;
			degp=degp->abs_list;
		}
		k=n_peel;
		for(i=1;i<n_inv;i++) {
			j=gp_inv[i];
			if(!(gp_flag[j]&2)) {
				gp_inv1[k++]=j;
			}
		}
		if(n_inv>1) qsort(gp_inv1,n_inv,sizeof(int),qfunc);
		for(i=0;i<n_inv;i++) cnv_gene[gp_inv1[i]]=i;
#ifdef TRACE_PEEL
		if(CHK_PEEL(TRACE_LEVEL_2)) {
			(void)fputs("Peeling",stdout);
			for(i=j=0;i<n_inv;i++) {
				k=gp_inv1[i];
				if(gp_flag[k]&6) {
					(void)fputc(j++?',':' ',stdout);
					(void)printf("%d",k);
				}
			}
			(void)fputs(" ->",stdout);
			if(n_peel==n_inv) (void)fputs(" *",stdout);
			else {
				for(i=j=0;i<n_inv;i++) {
					k=gp_inv1[i];
					if(!(gp_flag[k]&6)) {
						(void)fputc(j++?',':' ',stdout);
						(void)printf("%d",k);
					}
				}
			}
			fputs("  ",stdout);
		}
#endif
		/* Get list of rfunctions for nodes to be peeled */
		for(rfl=0,n_rf=i=0;i<n_inv;i++) {
			k=gp_inv1[i];
			if(gp_flag[k]&6) {
				pflag[i]=1;
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
#ifdef TRACE_PEEL
					if(CHK_PEEL(TRACE_LEVEL_3)) (void)printf("\n[%s:%d] %s() - Adding %p to free_rfp_list\n",__FILE__,__LINE__,__func__,(void *)rfp1);
#endif
					rfp1->next=free_rfp_list;
					free_rfp_list=rfp1;
					rfp1=rfp2;
				}
				rfuncs[k]=0;
			} else {
				pflag[i]=0;
				rfp3=rfuncs+k;
				rfp1=*rfp3;
				while(rfp1) {
					rf1=rfp1->rf;
					if(!rf1->flag) {
						for(k1=0;k1<rf1->n_inv;k1++) {
							k3=rf1->inv[k1];
							for(k2=0;k2<n_inv;k2++) if(gp_inv1[k2]==k3) break;
							if(k2==n_inv) break;
						}
						if(k1==rf1->n_inv) {
							rf1->next=rfl;
							rfl=rf1;
							rf1->flag=1;
							n_rf++;
						}
					}
					if(rf1->flag) {
#ifdef TRACE_PEEL
						if(CHK_PEEL(TRACE_LEVEL_3)) (void)printf("\n[%s:%d] %s() - Adding %p to free_rfp_list\n",__FILE__,__LINE__,__func__,(void *)rfp1);
#endif
						*rfp3=rfp1->next;
						rfp1->next=free_rfp_list;
						free_rfp_list=rfp1;
					} else rfp3=&rfp1->next;
					rfp1=*rfp3;
				}
			}
		}
#ifdef TRACE_PEEL
		if(CHK_PEEL(TRACE_LEVEL_2)) {
			rf1=rfl;
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
#endif
		/* Convert genes in rfunctions to references in gp_inv1 */
		rf1=rfl;
		while(rf1) {
			for(i=0;i<rf1->n_inv;i++) {
				j=rf1->inv[i];
				rf1->inv[i]=cnv_gene[j];
			}
			rf1=rf1->next;
		}
		/* Setup output rfunction */
		n_out=n_inv-n_peel;
		if(n_out) {
			if((rf=free_rf_list[n_out-1])) free_rf_list[n_out-1]=rf->next;
			else rf=get_new_rfunc(n_out,2);
			rf->next=0;
			rf->flag=0;
			xx=rf->p;
			nc=rf->nc;
			for(i=0;i<nc;i++) xx[i]=0.0;
			for(j=0,i=0;i<n_inv;i++) {
				if(!pflag[i]) {
					k1=gp_inv1[i];
					rf->inv[j++]=k1;
					if((rfp1=free_rfp_list)) free_rfp_list=rfp1->next;
					else rfp1=lk_malloc(sizeof(struct gp_rfunc_ptr));
					rfp1->next=rfuncs[k1];
					rfp1->rf=rf;
					rfuncs[k1]=rfp1;
					/* Fill in new adjacencies */
					for(i1=0;i1<i;i1++) {
						if(!pflag[i1]) {
							k2=gp_inv1[i1];
							rfnode=find_rfnode(k1,k2);
							if(!rfnode->rf) (void)find_rfnode(k2,k1);
						}
					}
				}
			}
		}
		/* If sampling, store peel_op */
		if(flag) {
			/* Get new peel_op and add to opslist */
			ops=get_new_peel_op(n_inv);
			ops->next=opslist;
			opslist=ops;
			/* Fill in details */
			ops->n_peel=n_peel;
			ops->rfl=rfl;
			(void)memcpy(ops->inv,gp_inv1,n_inv*sizeof(int));
			(void)memcpy(ops->pflag,pflag,n_inv*sizeof(int));
		}
		/* Check we have enough space for rfunction indices */
		if(n_rf>rf_idx_n) {
			rf_idx_n=n_rf;
			if(rf_idx) rf_idx=realloc_and_remember(rf_idx,sizeof(int)*(rf_idx_n+2*rf_idx_n*ngenes));
			else rf_idx=malloc_and_remember(sizeof(int)*(rf_idx_n+2*rf_idx_n*ngenes));
			fact[0]=rf_idx+rf_idx_n;
			rfx[0]=fact[0]+rf_idx_n*ngenes;
			for(k1=1;k1<ngenes;k1++) {
				fact[k1]=fact[k1-1]+rf_idx_n;
				rfx[k1]=rfx[k1-1]+rf_idx_n;
			}
		}
		/* Skip operation for the moment if we will be sampling and this is a terminal op (n_out==0) */
		if(!flag || n_out) {
			/* Perform peelop */
			/* Initialize frequency probs. and allele indices */
			if(n_genetic_groups>1) {
				for(z=1.0,i=n_inv-1;i>=0;i--) {
					gp_inv[i]=0;
					gp_z[i]=z;
					if(pflag[i]) z*=freq[grp[gp_inv1[i]]][0]; 
				}
			} else {
				z1=freq1[0];
				for(z=1.0,i=n_inv-1;i>=0;i--) {
					gp_inv[i]=0;
					gp_z[i]=z;
					if(pflag[i]) z*=z1;
				}
			}
			z3=z;
			/* Initialize output index */
			for(k2=1,k1=0;k1<n_inv;k1++) {
				if(pflag[k1]) ofact[k1]=0;
				else {
					ofact[k1]=k2;
					k2<<=1;
				}
			}
			out_idx=0;
			z1=0.0;
			if(n_rf==1) {
				/* Initialize R-Function indices */
				k=0;
				inv=rfl->inv;
				for(k1=0;k1<n_inv;k1++) fact[k1][0]=0;
				for(k2=1,k1=0;k1<rfl->n_inv;k1++,k2<<=1) {
					i2=inv[k1];
					fact[i2][0]=k2;
				}
				k3=0;
				rflp=rfl->p;
				/* Loop through all combinations */
				do {
#ifdef TRACE_PEEL
					if(CHK_PEEL(TRACE_LEVEL_4)) for(j=0;j<n_inv;j++) (void)printf("%d ",gp_inv[j]);
#endif
					/* Insert frequencies for peeled genes and probabilities from the rfunction */
					/* Check for index out-of-bounds */
					assert(k3<rfl->nc);
					z=z3*rflp[k3];
					if(n_out) xx[out_idx]+=z;
					z1+=z;
#ifdef TRACE_PEEL
					if(CHK_PEEL(TRACE_LEVEL_4)) {
						(void)printf("[%d,%g]) : %g %g",k3,rfl->p[k3],z3,z);
					}
#endif
					/* Bump indices to get next combination */
					for(i=0;i<n_inv;i++) {
						gp_inv[i]^=1;
						if(gp_inv[i]) {
							k3+=fact[i][0];
							out_idx+=ofact[i];
							break;
						}
						k3-=fact[i][0];
						out_idx-=ofact[i];
					}
					/* Update frequency contribution from peeled genes */
					if(i<n_inv) {
						if(n_genetic_groups==1) {
							z3=gp_z[i];
							if(pflag[i]) z3*=freq1[gp_inv[i]];
							for(k=i-1;k>=0;k--) {
								gp_z[k]=z3;
								if(pflag[k]) {
									z3*=freq1[gp_inv[k]];
								}
							}
						} else {
							z3=gp_z[i];
							if(pflag[i]) {
								k=gp_inv[i];
								z3*=freq[grp[gp_inv1[i]]][k];
							}
							for(k=i-1;k>=0;k--) {
								gp_z[k]=z3;
								if(pflag[k]) {
									k1=gp_inv[k];
									z3*=freq[grp[gp_inv1[k]]][k1];
								}
							}
						}
					}
				} while(i<n_inv);
			} else {
				/* Initialize R-Function indices */
				k=0;
				rf1=rfl;
				for(k1=0;k1<n_inv;k1++) rfxn[k1]=0;
				while(rf1) {
					inv=rf1->inv;
					for(k2=1,k1=0;k1<rf1->n_inv;k1++,k2<<=1) {
						i2=inv[k1];
						fact[i2][k]=k2;
						rfx[i2][rfxn[i2]++]=k;
					}
					rf_idx[k++]=0;
					rf1=rf1->next;
				}
				/* Loop through all combinations */
				do {
#ifdef TRACE_PEEL
					if(CHK_PEEL(TRACE_LEVEL_4)) for(j=0;j<n_inv;j++) (void)printf("%d ",gp_inv[j]);
#endif
					/* Insert frequencies for peeled genes */
					z=z3;
					/* Insert probabilities from rfunctions */
					idx=rf_idx;
					rf1=rfl;
					while(rf1) {
						/* Check for out-of-bounds index */
						assert(*idx<rf1->nc);
						z*=rf1->p[*(idx++)];
						rf1=rf1->next;
					}
					if(n_out) xx[out_idx]+=z;
					z1+=z;
#ifdef TRACE_PEEL
					if(CHK_PEEL(TRACE_LEVEL_4)) {
						k1=0;
						rf1=rfl;
						while(rf1) {
							(void)fputc(k1?',':'(',stdout);
							(void)printf("[%d,%g]",rf_idx[k1],rf1->p[rf_idx[k1]]);
							k1++;
							rf1=rf1->next;
						}
						if(k1) (void)fputs(") ",stdout);
						(void)printf(": %g %g\n",z3,z);
					}
#endif
					/* Bump indices to get next combination */
					for(i=0;i<n_inv;i++) {
						gp_inv[i]^=1;
						ff=fact[i];
						k1=rfxn[i];
						inv=rfx[i];
						if(gp_inv[i]) {
							for(k=0;k<k1;k++) {
								k2=inv[k];
								rf_idx[k2]+=ff[k2];
							}
							out_idx+=ofact[i];
							break;
						}
						for(k=0;k<k1;k++) {
							k2=inv[k];
							rf_idx[k2]-=ff[k2];
						}
						out_idx-=ofact[i];
					}
					/* Update frequency contribution from peeled genes */
					if(i<n_inv) {
						if(n_genetic_groups==1) {
							z3=gp_z[i];
							if(pflag[i]) z3*=freq1[gp_inv[i]];
							for(k=i-1;k>=0;k--) {
								gp_z[k]=z3;
								if(pflag[k]) {
									z3*=freq1[gp_inv[k]];
								}
							}
						} else {
							z3=gp_z[i];
							if(pflag[i]) {
								k=gp_inv[i];
								z3*=freq[grp[gp_inv1[i]]][k];
							}
							for(k=i-1;k>=0;k--) {
								gp_z[k]=z3;
								if(pflag[k]) {
									k1=gp_inv[k];
									z3*=freq[grp[gp_inv1[k]]][k1];
								}
							}
						}
					}
				} while(i<n_inv);
			}
			if(z1<=0.0) {
				prt_peel_trace(TRACE_LEVEL_1,"Returning err=-1, p = 0.0 (zero prob.) (c)\n",stdout);
				clean_up(ngenes,rfl);
				*err=-1;
				return 0.0;
			}
			/* Rescale output rfunction */
			if(n_out) {
				z=1.0/z1;
				for(i=0;i<nc;i++) xx[i]*=z;
			}
			like+=log(z1);
			prt_peel_trace(TRACE_LEVEL_3,"Rescaling output function (z1 = %g)\n",z1);
		}
		/* Remove references to used rfunctions from *rfuncs */
		for(i=0;i<n_inv;i++) {
			if(!pflag[i]) {
				k=gp_inv1[i];
				rfp3=rfuncs+k;
				rfp1=*rfp3;
				while(rfp1) {
					rfp2=rfp1->next;
					if(rfp1->rf->flag) {
						prt_peel_trace(TRACE_LEVEL_3,"[%s:%d] %s() - Adding %p to free_rfp_list\n",__FILE__,__LINE__,__func__,(void *)rfp1);
						*rfp3=rfp1->next;
						rfp1->next=free_rfp_list;
						free_rfp_list=rfp1;
					} else rfp3=&rfp1->next;
					rfp1=rfp2;
				}
			}
		}
		/* Free used rfunctions (unless sampling) */
		if(!flag) free_rfunc_list(rfl);
		/* Remove peeled nodes from degree list */
		for(i=0;i<n_inv;i++) {
			if(pflag[i]) {
				ngenes1--;
				k=gp_inv1[i];
				if(!(gp_flag[k]&2)) {
					prt_peel_trace(TRACE_LEVEL_2,"Removing %d from degree list %d (a)\n",k,deglist[k].deg);
					if(deglist[k].prev) deglist[k].prev->next=deglist[k].next;
					else {
						k1=deglist[k].deg;
						first[k1]=deglist[k].next;
					}
					if(deglist[k].next) deglist[k].next->prev=deglist[k].prev;
				}
				del_rfrow(k);
				/* Mark gene as peeled */
				gp_flag[k]|=2;
			}
		}
		/* Recalculate degree for output nodes */
		for(i=0;i<n_inv;i++) {
			if(!pflag[i]) {
				k=gp_inv1[i];
				if(!gp_flag[k]) {
					k1=calc_degree(k);
					if(k1>=0) {
						k2=0;
						for(j=1;j<=k1;j++) if(!gp_flag[gp_inv[j]]) k2++;
					} else k2=k1;
					k1=deglist[k].deg;
					if(k1!=k2) {
						prt_peel_trace(TRACE_LEVEL_2,"Moving %d from degree list %d to %d (a)\n",k,k1,k2);
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
	if(flag && opslist) {
#ifdef TRACE_PEEL
		if(CHK_PEEL(TRACE_LEVEL_1)) (void)fputs("Starting sampling phase\n",stdout);
#endif
		/* Zero sampled alleles */
		(void)memset(gp_inv1,0,sizeof(int)*ngenes);
		ops=opslist;
		while(ops) {
			n_inv=ops->n_inv;
			n_peel=ops->n_peel;
			n_out=n_inv-n_peel;
			inv1=ops->inv;
			pflag1=ops->pflag;
			rfl=ops->rfl;
#ifdef TRACE_PEEL
			if(CHK_PEEL(TRACE_LEVEL_2)) {
				(void)fputs("Sampling",stdout);
				for(i=j=0;i<n_inv;i++) {
					k=inv1[i];
					if(pflag1[i]) {
						(void)fputc(j++?',':' ',stdout);
						(void)printf("%d",k);
					}
				}
				(void)fputs(" <-",stdout);
				if(n_peel==n_inv) (void)fputs(" *",stdout);
				else {
					for(i=j=0;i<n_inv;i++) {
						k=inv1[i];
						if(!(pflag1[i])) {
							(void)fputc(j++?',':' ',stdout);
							(void)printf("%d",k);
						}
					}
				}
				fputc('\n',stdout);
			}
#endif
			/* Set up temporary storage */
			nc=(int)(.5+exp(log((double)n_all)*(double)n_peel));
			if(nc>gp_xx_size) {
				if(gp_xx) {
					gp_xx=realloc_and_remember(gp_xx,sizeof(double)*nc);
				} else gp_xx=malloc_and_remember(sizeof(double)*nc);
				gp_xx_size=nc;
			}
			xx=gp_xx;
			for(i=0;i<nc;i++) xx[i]=0.0;
			/* Set allele indices and frequency probs. */
			if(n_genetic_groups>1) {
				for(z=1.0,i=n_inv-1;i>=0;i--) {
					gp_z[i]=z;
					if(pflag1[i]) {
						gp_inv[i]=0;
						z*=freq[grp[inv1[i]]][0];
					} else {
						k=gp_inv1[inv1[i]]-1;
						/* Check allele has been sampled */
						assert(k>=0 && k<n_all);
						gp_inv[i]=k;
					}
				}
			} else {
				z1=freq1[0];
				for(z=1.0,i=n_inv-1;i>=0;i--) {
					gp_z[i]=z;
					if(pflag1[i]) {
						gp_inv[i]=0;
						z*=z1;
					} else {
						k=gp_inv1[inv1[i]]-1;
						/* Check allele has been sampled */
						assert(k>=0 && k<n_all);
						gp_inv[i]=k;
					}
				}
			}
			/* Initialize output index */
			for(k2=1,k1=0;k1<n_inv;k1++) {
				if(!pflag1[k1]) ofact[k1]=0;
				else {
					ofact[k1]=k2;
					k2<<=1;
				}
			}
			out_idx=0;
			z1=0.0;
			z3=z;
			k=0;
			/* Set up R-Function indices */
			for(k1=0;k1<n_inv;k1++) rfxn[k1]=0;
			rf1=rfl;
			rf2=&rfl;
			while(rf1) {
				inv=rf1->inv;
				for(k2=1,k1=i=j=0;k1<rf1->n_inv;k1++,k2<<=1) {
					i2=inv[k1];
					if(pflag1[i2]) {
						j++;
						fact[i2][k]=k2;
						rfx[i2][rfxn[i2]++]=k;
					} else i+=k2*gp_inv[i2];
				}
				/* R-Functions only has already sampled genes, so remove from list */
				if(!j) {
					rf=rf1->next;
					k1=rf1->n_inv-1;
					prt_peel_trace(TRACE_LEVEL_3,"[%s:%d] %s() - Adding %p to free_rf_list[%d]\n",__FILE__,__LINE__,__func__,(void *)rf1,k);
					rf1->next=free_rf_list[k1];
					free_rf_list[k1]=rf1;
					*rf2=rf;
				} else {
					rf_idx[k++]=i;
					rf2=&rf1->next;
				}
				rf1=*rf2;
			}
			if(k==1) {
				k3=rf_idx[0];
				rflp=rfl->p;
				/* Loop through all combinations */	
				do {
					/* Insert frequencies for genes we are sampling*/
					/* Insert probabilities from rfunctions */
					/* Check for out-of-bounds index */
					assert(k3<rfl->nc);
					z=z3*rflp[k3];
					xx[out_idx]+=z;
					z1+=z;
					/* Bump indices to get next combination */
					for(i=0;i<n_inv;i++) {
						if(pflag1[i]) {
							gp_inv[i]^=1;
							if(gp_inv[i]) {
								k3+=fact[i][0];
								out_idx+=ofact[i];
								break;
							}
							k3-=fact[i][0];
							out_idx-=ofact[i];
						}
					}
					/* Update frequency contribution from the genes we are sampling */
					if(i<n_inv) {
						assert(pflag[i]);
						if(n_genetic_groups==1) {
							z3=gp_z[i];
							z3*=freq1[gp_inv[i]];
							for(k=i-1;k>=0;k--) {
								gp_z[k]=z3;
								if(pflag1[k]) {
									z3*=freq1[gp_inv[k]];
								}
							}
						} else {
							z3=gp_z[i];
							k=gp_inv[i];
							z3*=freq[grp[inv1[i]]][k];
							for(k=i-1;k>=0;k--) {
								gp_z[k]=z3;
								if(pflag1[k]) {
									k1=gp_inv[k];
									z3*=freq[grp[inv1[k]]][k1];
								}
							}
						}
					}
				} while(i<n_inv);
			} else {
				/* Loop through all combinations */
				do {
					/* Insert frequencies for genes we are sampling*/
					z=z3;
					/* Insert probabilities from rfunctions */
					idx=rf_idx;
					rf1=rfl;
					while(rf1) {
						/* Check for out-of-bounds index */
						assert(*idx<rf1->nc);
						z*=rf1->p[*(idx++)];
						rf1=rf1->next;
					}
					xx[out_idx]+=z;
					z1+=z;
					/* Bump indices to get next combination */
					for(i=0;i<n_inv;i++) {
						if(pflag1[i]) {
							gp_inv[i]^=1;
							ff=fact[i];
							k1=rfxn[i];
							inv=rfx[i];
							if(gp_inv[i]) {
								for(k=0;k<k1;k++) {
									k2=inv[k];
									rf_idx[k2]+=ff[k2];
								}
								out_idx+=ofact[i];
								break;
							}
							for(k=0;k<k1;k++) {
								k2=inv[k];
								rf_idx[k2]-=ff[k2];
							}
							out_idx-=ofact[i];
						}
					}
					/* Update frequency contribution from the genes we are sampling */
					if(i<n_inv) {
						assert(pflag1[i]);
						if(n_genetic_groups==1) {
							z3=gp_z[i];
							z3*=freq1[gp_inv[i]];
							for(k=i-1;k>=0;k--) {
								gp_z[k]=z3;
								if(pflag1[k]) {
									z3*=freq1[gp_inv[k]];
								}
							}
						} else {
							z3=gp_z[i];
							k=gp_inv[i];
							z3*=freq[grp[inv1[i]]][k];
							for(k=i-1;k>=0;k--) {
								gp_z[k]=z3;
								if(pflag1[k]) {
									k1=gp_inv[k];
									z3*=freq[grp[inv1[k]]][k1];
								}
							}
						}
					}
				} while(i<n_inv);
			}
			if(z1<=0.0) {
				/* If n_out > 0 then we're sampling so any problems should have been detected at the peeling stage */
				if(n_out) {
					printf("z1 = %g, n_rf=%d",z1,n_rf);
					ABT_FUNC("Internal error = zero prob. in sampling step\n");
				}
				/* If n_out==0 then this this technically still a peeling op, so a zero here might be valid (not a bug) */
#ifdef TRACE_PEEL
				if(CHK_PEEL(TRACE_LEVEL_1)) (void)fputs("Returning err=-1, p = 0.0 (zero prob.) (d)\n",stdout);
#endif
				clean_up(ngenes,rfl);
				*err=-1;
				return 0.0;
			}
			if(!n_out) like+=log(z1);
			/* Free used rfunctions */
			free_rfunc_list(rfl);
			/* Sample! */
			z=ranf()*z1;
			for(z3=0.0,i=0;i<nc;i++) {
				if(xx[i]>0.0) {
					z3+=xx[i];
					if(z<=z3) break;
				}
			}
			for(j=0;j<n_inv;j++) {
				if(pflag1[j]) {
					gp_inv1[inv1[j]]=(i&1)+1;
					i>>=1;
				}
			}
			ops1=ops->next;
			ops->next=free_ops_list[n_inv-1];
			free_ops_list[n_inv-1]=ops;
			ops=ops1;
		}
		/* Sample unobserved founder genes */
		for(i=0;i<ngenes;i++) {
			if(!gp_inv1[i]) {
				k=n_genetic_groups==1?0:grp[i];
				z=ranf();
				for(z1=0.0,j=0;j<2;j++) {
					z1+=freq[k][j];
					if(z<=z1) break;
				}
				gp_inv1[i]=j+1;
			}
		}
		/* Pass founder alleles down through pedigree */
		cpsize=loki->pedigree->comp_size[c];
		i=loki->pedigree->comp_start[c];
		seg=locus->seg;
		ff=loki->models->founder_flag;
		while(cpsize--) {
			id=id_array+i;
			a1=genes[X_MAT][i]-1;
			a2=genes[X_PAT][i]-1;
			k1=id->allele[X_MAT]=gp_inv1[a1];
			k2=id->allele[X_PAT]=gp_inv1[a2];
			if((flag&2) && !si && locus->flag&TL_LINKED) {
				/* Remove ambiguous segregations (cond. on ordered genotypes) */
				/* Don't bother if unlinked as it buys us nothing */
				if(ff[i]) seg[X_MAT][i]=seg[X_PAT][i]=-1;
				else {
					idd=id->dam;
					/* Check sampled allele in Dam */
					assert(k1==id_array[idd-1].allele[seg[X_MAT][i]]);
					if(id_array[idd-1].allele[X_MAT]==id_array[idd-1].allele[X_PAT]) seg[X_MAT][i]=-2;
					ids=id->sire;
					/* Check sampled allele in Sire */
					assert(k2==id_array[ids-1].allele[seg[X_PAT][i]]);
					if(id_array[ids-1].allele[X_MAT]==id_array[ids-1].allele[X_PAT]) seg[X_PAT][i]=-2;
				}
			}
			/* Correct residuals for new genotypes */
			if(loki->models->pruned_flag[i]) k=0;
			else {
				k=k1>k2?k1*(k1-1)/2+k2:k2*(k2-1)/2+k1;
				if(id->res[0]) {
					nrec=id->n_rec;
					k1=(locus->flag&LOCUS_SAMPLED)?gt[i]:1;
					if(k1!=k) {
						z=(k1>1)?eff[k1-2]:0.0;
						if(k>1) z-=eff[k-2];
						for(rec=0;rec<nrec;rec++) {
							id->res[0][rec]+=z;
						}
					}
				}
			}
			gt[i++]=k;
		}
		/* Get allele counts for frequency update */
		if(flag&4) {
			if(n_genetic_groups==1) {
				count1=loki->peel->seg_count[0];
				for(i=0;i<ngenes;i++) {
					k=gp_inv1[i]-1;
					count1[k]+=1.0;
				}
			} else {
				for(i=0;i<ngenes;i++) {
					count1=loki->peel->seg_count[grp[i]];
					k=gp_inv1[i]-1;
					count1[k]+=1.0;
				}
			}
		}
	}
	if(flag) locus->flag|=LOCUS_SAMPLED;
#ifdef TRACE_PEEL
	if(CHK_PEEL(TRACE_LEVEL_1)) (void)printf("Returning like = %g, p = %g\n",like,exp(like));
#endif
	return like;
}
