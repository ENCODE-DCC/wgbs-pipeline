/****************************************************************************
 *                                                                          *
 *     Loki - Programs for genetic analysis of complex traits using MCMC    *
 *                                                                          *
 *             Simon Heath - University of Washington                       *
 *                                                                          *
 *                        July 1997                                         *
 *                                                                          *
 * prep_do_peel_op.c:                                                       *
 *                                                                          *
 * Perform complex peeling operation                                        *
 *                                                                          *
 * Copyright (C) Simon C. Heath 1997, 2000, 2002                            *
 * This is free software.  You can distribute it and/or modify it           *
 * under the terms of the Modified BSD license, see the file COPYING        *
 *                                                                          *
 ****************************************************************************/

#include <config.h>
#include <stdlib.h>
#ifdef USE_DMALLOC
#include <dmalloc.h>
#endif
#include <math.h>
#include <stdio.h>

#include "utils.h"
#include "scan.h"
#include "prep_peel.h"

#define HASHTABLE_SIZE 2053  /* Should be prime */

static int hb_size=2048,n_bits,n_terms,hash_mode;
static struct bin_node *hashtable[HASHTABLE_SIZE];
static struct hash_block *first_hash_block=0,*hash_block;
static struct R_Func *rf_array;
int total_terms=0,total_comb=0;

static void get_gts(lk_ulong x,const int n,int *gt)
{
	int i;
	lk_ulong a;
	
	a=(1<<n_bits)-1;
	i=0;
	while(x)	{
		if(i>=n) ABT_FUNC("Internal error - invalid index value\n");
		gt[i++]=1+(int)(x&a);
		x>>=n_bits;
	}
	for(;i<n;i++) gt[i]=1;
}

static struct bin_node *get_new_element(lk_ulong idx)
{
	struct bin_node *element;
	lk_ulong *p;
	
	while(hash_block->ptr>=hash_block->size) {
		if(!hash_block->next) {
			if(!(hash_block->next=malloc(sizeof(struct hash_block)))) ABT_FUNC(MMsg);
			hash_block=hash_block->next;
			hash_block->next=0;
			if(!(hash_block->elements=malloc(sizeof(struct bin_node)*hb_size))) ABT_FUNC(MMsg);
			if(!(hash_block->idx=malloc(sizeof(lk_ulong)*hb_size))) ABT_FUNC(MMsg);
			hash_block->size=hb_size;
		} else hash_block=hash_block->next;
		hash_block->ptr=0;
	}
	element=hash_block->elements+hash_block->ptr;
	p=hash_block->idx+hash_block->ptr++;
	*p=idx;
	element->data=p;
	element->left=element->right=0;
	element->balance=0;
	n_terms++;
	return element;
}

static struct bin_node *insert_node(struct bin_node *node,lk_ulong idx,int *bal)
{
	int bb;
	lk_ulong idx1;
	
	idx1=*(lk_ulong *)node->data;
	if(idx!=idx1) {
		bb=node->balance;
		if(idx<idx1) {
			if(node->left) node->left=insert_node(node->left,idx,bal);
			else {
				node->left=get_new_element(idx);
				*bal=0;
			}
			if(!(*bal)) {
				switch(bb) {
				 case -1:
					node=rotate_left(node);
					*bal=1;
					break;
				 case 0:
					node->balance=-1;
					break;
				 case 1:
					node->balance=0;
					*bal=1;
				}
			}
		} else {
			if(node->right) node->right=insert_node(node->right,idx,bal);
			else {
				node->right=get_new_element(idx);
				*bal=0;
			}
			if(!(*bal)) {
				switch(bb) {
				 case -1:
					node->balance=0;
					*bal=1;
					break;
				 case 0:
					node->balance=1;
					break;
				 case 1:
					node=rotate_right(node);
					*bal=1;
				}
			}			
		}
	} else *bal=1;
	return node;
}

static lk_ulong get_index1(int n,int *gt)
{
	int i;
	lk_ulong x;
	
	x=0;
	for(i=n-1;i>=0;i--) {
		x<<=n_bits;
		x|=gt[i]-1;
	}
	return x;
}

static void get_nodes(struct bin_node *node,lk_ulong *tl,int *j)
{
	if(node->left)	get_nodes(node->left,tl,j);
	tl[(*j)++]=*(lk_ulong *)node->data;
	if(node->right) get_nodes(node->right,tl,j);
}

static void add_term(int n_out,int *gt_store)
{
	lk_ulong idx;
	int i,j;
	
	idx=get_index1(n_out,gt_store);
	if(hash_mode) {
		i=(int)idx;
		if(!hashtable[i])	hashtable[i]=get_new_element(idx);
	} else {
		i=(int)(idx%HASHTABLE_SIZE); /* Get hash index */
		if(hashtable[i]) hashtable[i]=insert_node(hashtable[i],idx,&j);
		else hashtable[i]=get_new_element(idx);
	}
}

/* Performs logical transmission check */
static int check_trans(int i,int par_flag,const int id,lk_ulong *req_set[])
{
	int par,al,al1,j;
	lk_ulong m,a,m1;
	
	if(par_flag==X_MAT) par=id_array[i-1].dam;
	else par=id_array[i-1].sire;
	j=ped_recode1[i-1]-1-id;
	al=id_array[i-1].allele[par_flag]-1;
	a=req_set[par_flag][j];
	m=1<<al;
	if(a&m) m=a;
	al1=id_array[par-1].allele[X_MAT]-1;
	m1=1<<al1;
	if(m1&m) return 0;
	al1=id_array[par-1].allele[X_PAT]-1;
	m1=1<<al1;
	if(m1&m) return 0;
	return 1;
}

void free_hash_blocks(void)
{
	struct hash_block *hb1;
	
	hash_block=first_hash_block;
	while(hash_block) {
		if(hash_block->elements) free(hash_block->elements);
		if(hash_block->idx) free(hash_block->idx);
		hb1=hash_block->next;
		free(hash_block);
		hash_block=hb1;
	}
	first_hash_block=0;
	total_terms=total_comb=0;
}

static int qs_func(const void *p1,const void *p2)
{
	int i1,i2,k1,k2;
	
	i1= rf_array[k1=*((int *)p1)].n_ind;
	i2= rf_array[k2=*((int *)p2)].n_ind;
	if(rf_array[k1].n_terms<rf_array[k2].n_terms) return 1;
	if(rf_array[k1].n_terms>rf_array[k2].n_terms) return -1;
	if(i1<i2) return 1;
	if(i1>i2) return -1; 
	return 0;
}

int do_peel_op(const struct Complex_Element *element,struct R_Func *r_func,const int n_all,const int id,lk_ulong **all_set,lk_ulong *req_set[])
{
	int i,j,k,k1,k2,k3,k4,n_out,n_peel,n_inv,*inv,n_rf,n_ind,n_other,ef,ef1,n_comb;
	int *gt_store,*gt_store1,*other_ptr,*other_list,*rf_ptr,*jnt_list,ht_size;
	double max_terms,fill;
	lk_ulong a,b,c,*tl;
	struct bin_node *elem;

	if(!first_hash_block) {
		if(!(first_hash_block=malloc(sizeof(struct hash_block)))) ABT_FUNC(MMsg);
		if(!(first_hash_block->elements=malloc(sizeof(struct bin_node)*hb_size))) ABT_FUNC(MMsg);
		if(!(first_hash_block->idx=malloc(sizeof(lk_ulong)*hb_size))) ABT_FUNC(MMsg);
		first_hash_block->next=0;
		first_hash_block->size=hb_size;
	}
	hash_block=first_hash_block;
	while(hash_block)	{
		hash_block->ptr=0;
		hash_block=hash_block->next;
	}
	hash_block=first_hash_block;
	n_terms=n_comb=0;
	n_bits=num_bits(n_all);
	n_inv=element->n_involved;
	n_peel=element->n_peel;
	n_out=n_inv-n_peel;
	inv=element->involved;
	if(n_bits*n_out>(int)LK_LONG_BIT) {
		(void)fprintf(stderr,"\nToo many individuals in output R-Function\nn_out = %d, n_all = %d, n_bits = %d, required size = %d, LONG_BIT = %d\n",n_out,n_all,n_bits,n_out*n_bits,(int)LK_LONG_BIT);
		ABT_FUNC(AbMsg);
	}
	max_terms=log((double)n_all)*(double)n_out;
	hash_mode=(log(2.0)*n_bits*n_out<log((double)HASHTABLE_SIZE));
	ht_size=hash_mode?1<<(n_bits*n_out):HASHTABLE_SIZE;
	for(i=0;i<ht_size;i++) hashtable[i]=0;
	n_rf=element->n_rfuncs;
	rf_array=r_func;
	if(n_rf>1) gnu_qsort(element->index,(size_t)n_rf,sizeof(int),qs_func);
#ifdef TRACE_PEEL
	if(CHK_PEEL(TRACE_LEVEL1)) (void)printf("In do_peel_op(), n_inv=%d, n_out=%d, n_rf=%d\n",n_inv,n_out,n_rf);
	if(CHK_PEEL(TRACE_LEVEL2)) {
		for(i=0;i<n_peel;i++) {
			if(i) fputc(',',stdout);
			print_orig_allele_id(stdout,inv[i]);
		}
		fputs("->",stdout);
		if(n_out) {
			for(;i<n_inv;i++) {
				if(i>n_peel) fputc(',',stdout);
				print_orig_allele_id(stdout,inv[i]);
			}
		} else fputc('*',stdout);
		fputc('\n',stdout);
	}
#endif
	if(!(gt_store=calloc((size_t)(5*n_inv+n_rf),sizeof(int)))) ABT_FUNC(MMsg);
	other_list=gt_store+n_inv;
	other_ptr=other_list+n_inv;
	jnt_list=other_ptr+n_inv;
	rf_ptr=jnt_list+n_inv;
	gt_store1=rf_ptr+n_rf;
	for(i=0;i<n_inv;i++) gt_store[i]=0;
	if(n_rf) {
		for(i=0;i<n_rf;i++) {
			j=element->index[i];
			n_ind=r_func[j].n_ind;
#ifdef TRACE_PEEL
			if(CHK_PEEL(TRACE_LEVEL2)) {
				for(k=0;k<n_ind;k++) {
					fputc(k?',':'(',stdout);
					k1=r_func[j].id_list[k];
					print_orig_allele_id(stdout,inv[k1]);
				}
				fputc(')',stdout);
			}
#endif
			a=c=0;
			b=(1<<n_bits)-1;
			for(k=n_ind-1;k>=0;k--)	{
				a<<=n_bits;
				k1=r_func[j].id_list[k];
				if(gt_store[k1]) a|=b;
				else {
					c|=1<<k1;
					gt_store[k1]=1;
				}
			}
			r_func[j].mask=a;
			r_func[j].mask1=c;
#ifdef TRACE_PEEL
			if(CHK_PEEL(TRACE_LEVEL2)) printf("[%d] ",r_func[j].n_terms);
#endif
		}	
#ifdef TRACE_PEEL
		if(CHK_PEEL(TRACE_LEVEL2)) fputc('\n',stdout);
#endif
	}
	for(n_other=i=0;i<n_inv;i++) if(!gt_store[i]) {
		if(element->flags[i]&(HAP_JNT|HAP_DAT))	{
			k= -inv[i];
			for(j=0;j<n_inv;j++) if(inv[j]==k) {
				jnt_list[n_other]=j;
				break;
			}
		}
		other_list[n_other++]=i;
	}
	ef1=0;
	for(i=0;i<n_inv;i++) gt_store[i]=0;
	i=0;
	for(;;) {
		for(;i<n_rf;i++) {
			j=element->index[i];
			ef=0;
			n_ind=r_func[j].n_ind;
			a=r_func[j].mask;
			if(a) {
				b=0;
				for(k=n_ind-1;k>=0;k--)	{
					b<<=n_bits;
					k1=r_func[j].id_list[k];
					if(gt_store[k1]) b|=gt_store[k1]-1;
				}
				b&=a;
				k1=r_func[j].n_terms;
				tl=r_func[j].index;
				for(k=rf_ptr[i];k<k1;k++)
				  if((tl[k]&a)==b) break;
				rf_ptr[i]=k;
			}
			if(rf_ptr[i]<r_func[j].n_terms) {
				get_gts(r_func[j].index[rf_ptr[i]],n_ind,gt_store1);
				for(k=0;k<n_ind;k++)	{
					k1=r_func[j].id_list[k];
					gt_store[k1]=gt_store1[k];
				}
			} else { /* Inconsistency found - try new combination */
				a=r_func[j].mask1;
				k=0;
				while(a)	{
					if(a&1) gt_store[k]=0;
					k++;
					a>>=1;
				}
				do	{
					if(!i) {
						ef1=1;
						break;
					}
					rf_ptr[i--]=0;
					rf_ptr[i]++;
					j=element->index[i];
					a=r_func[j].mask1;
					k=0;
					while(a)	{
						if(a&1) gt_store[k]=0;
						k++;
						a>>=1;
					}
				} while(rf_ptr[i]>=r_func[j].n_terms);
				i--;
			}
			if(ef1) break;
		}
		if(ef1) break;
		for(ef=i=0;i<=n_other;i++)	{
			if(i==n_other)	{
				if(!ef) {
					for(k=0;k<n_inv;k++)	{
						if(inv[k]>0) {
							id_array[inv[k]-1].allele[X_MAT]=gt_store[k];
							if(element->flags[k]&(HAP_JNT|HAP_DAT))	{
								j=ped_recode1[inv[k]-1]-1-id;
								for(k1=0;k1<n_inv;k1++) if(inv[k1]== -inv[k]) {
									if(!(all_set[gt_store[k]-1][j]&(1<<(gt_store[k1]-1)))) ef=1;
									break;
								}
							}
						} else id_array[-1-inv[k]].allele[X_PAT]=gt_store[k];
						if(ef) break;
					}
					/* Add Par. off transmission probs. */
					if(!ef) for(k=0;k<n_inv;k++) {
						if(element->flags[k]&HAD_P)
						  if((ef=check_trans(-inv[k],X_PAT,id,req_set))) break;
						if(element->flags[k]&HAD_M)
						  if((ef=check_trans(inv[k],X_MAT,id,req_set))) break;
						if(inv[k]>0 && element->flags[k]&(HAP_JNT|HAP_DAT))	{
							for(k1=0;k1<n_inv;k1++) if(inv[k1]== -inv[k]) {	
								j=ped_recode1[inv[k]-1]-1-id;
								if(!(all_set[gt_store[k]-1][j]&(1<<(gt_store[k1]-1)))) ef=1;
								break;
							}
							if(ef) break;
						}
					}
					if(!ef) {
						add_term(n_out,gt_store+n_peel);
						n_comb++;
						ef=1;
					}
				}
			}
			if(ef) {
				i--;
				while(i>=0)	{
					other_ptr[i]++;
					gt_store[other_list[i]]=0;
					if(other_ptr[i]<n_all) break;
					other_ptr[i--]=0;
				}
				i--;
				if(i<-1) break;
				ef=0;
				continue;
			}
			k1=other_list[i];
			if(gt_store[k1]) continue;
			j=ped_recode1[abs(inv[k1])-1]-1-id;
			if(element->flags[k1]&(HAP_JNT|HAP_DAT)) {
				k2=jnt_list[i];
				if(inv[k1]>0) {
					if(gt_store[k2]) {
						b=1<<(gt_store[k2]-1);
						for(k=other_ptr[i];k<n_all;k++) if(all_set[k][j]&b) {
							gt_store[k1]=k+1;
							break;
						}
						if((other_ptr[i]=k)==n_all) ef=1;
					} else {
						for(k=0;k<n_other;k++) if(other_list[k]==k2) break;
						if(k==n_other) ABT_FUNC("Internal error - no other?\n");
						if(i<k) {
							for(k3=other_ptr[i];k3<n_all;k3++) {
								a=all_set[k3][j];
								if(!a) continue;
								for(k4=other_ptr[k];k4<n_all;k4++) if(a&(1<<k4)) {
									gt_store[k1]=k3+1;
									gt_store[k2]=k4+1;
									break;
								}
								if(k4==n_all) other_ptr[k]=0;
								else {
									other_ptr[k]=k4;
									break;
								}
							}
							if((other_ptr[i]=k3)==n_all) ef=1;
						} else {
							for(k4=other_ptr[k];k4<n_all;k4++) {
								a=1<<k4;
								if(req_set[2][j]&a) {
									for(k3=other_ptr[i];k3<n_all;k3++) if(all_set[k3][j]&a) {
										gt_store[k1]=k3+1;
										gt_store[k2]=k4+1;
										break;
									}
									if(k3==n_all) other_ptr[i]=0;
									else {
										other_ptr[i]=k3;
										break;
									}
								} else other_ptr[i]=0;
							}
							if((other_ptr[k]=k4)==n_all) ef=1;
						}
					}
				} else {
					if(gt_store[k2]) {
						b=all_set[gt_store[k2]-1][j];
						for(k=other_ptr[i];k<n_all;k++) if(b&(1<<k)) {
							gt_store[k1]=k+1;
							break;
						}
						if((other_ptr[i]=k)==n_all) ef=1;
					} else {
						for(k=i+1;k<n_other;k++) if(other_list[k]==k2) break;
						if(k==n_other) ABT_FUNC("Internal error - why for this happen?\n");
						if(i<k) {
							for(k3=other_ptr[i];k3<n_all;k3++) {
								a=1<<k3;
								if(req_set[2][j]&a) {
									for(k4=other_ptr[k];k4<n_all;k4++) if(all_set[k4][j]&a) {
										gt_store[k1]=k3+1;
										gt_store[k2]=k4+1;
										break;
									}
									if(k4==n_all) other_ptr[k]=0;
									else {
										other_ptr[k]=k4;
										break;
									}
								} else other_ptr[k]=0;
							}
							if((other_ptr[i]=k3)==n_all) ef=1;
						} else {
							for(k4=other_ptr[k];k4<n_all;k4++) {
								a=all_set[k4][j];
								if(!a) continue;
								for(k3=other_ptr[i];k3<n_all;k3++) if(a&(1<<k3)) {
									gt_store[k1]=k3+1;
									gt_store[k2]=k4+1;
									break;
								}
								if(k3==n_all) other_ptr[i]=0;
								else {
									other_ptr[i]=k3;
									break;
								}
							}
							if((other_ptr[k]=k4)==n_all) ef=1;
						}
					}
				}
			} else {
				if(inv[k1]>0) {
					for(k=other_ptr[i];k<n_all;k++) if(all_set[k][j]) {
						gt_store[k1]=k+1;
						break;
					}
					if((other_ptr[i]=k)==n_all) ef=1;
				} else {
					a=req_set[2][j];
					for(k=other_ptr[i];k<n_all;k++) if(a&(1<<k))	{
						gt_store[k1]=k+1;
						break;
					}
					if((other_ptr[i]=k)==n_all) ef=1;
				}
			}
		}
		for(i=0;i<n_other;i++) gt_store[other_list[i]]=0;
		/* Cycle through R-Functions */
		if(n_rf)	{
			i=n_rf-1;
			do	{
				rf_ptr[i]++;
				j=element->index[i];
				a=r_func[j].mask1;
				if(a)	{
					k=0;
					while(a)	{
						if(a&1) gt_store[k]=0;
						k++;
						a>>=1;
					}
					if(rf_ptr[i]<r_func[j].n_terms) break;
				}
				rf_ptr[i--]=0;
			} while(i>=0);
			if(i>=0) ef=0;
		} else i=0;
		if(ef) break;
	}
	if(!n_terms) {
		free(gt_store);
		return 1;
	}
	i=element->out_index;
	if(i>=0) {
		r_func[i].n_terms=n_terms;
		for(j=0;j<n_out;j++) r_func[i].id_list[j]=inv[j+n_peel];
		if(!(tl=malloc(sizeof(lk_ulong)*n_terms))) ABT_FUNC(MMsg);
		r_func[i].index=tl;
		for(j=k=0;k<ht_size;k++) {
			elem=hashtable[k];
			if(elem)	{
				if(hash_mode) tl[j++]=(lk_ulong)k;
				else get_nodes(elem,tl,&j);
			}
		}
		fill=exp(log((double)n_terms)-max_terms);
	} else fill=1.0;
	total_terms+=n_terms;
	total_comb+=n_comb;
#ifdef TRACE_PEEL
	if(CHK_PEEL(TRACE_LEVEL1))
	  (void)printf("-> No. non-zero terms = %d (%g%% full), non-zero combs %d, %d, %d\n",n_terms,100.0*fill,n_comb,total_terms,total_comb);
#endif
	free(gt_store);
	return 0;
}
