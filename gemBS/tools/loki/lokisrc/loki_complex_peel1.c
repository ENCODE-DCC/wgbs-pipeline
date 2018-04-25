/****************************************************************************
 *                                                                          *
 *     Loki - Programs for genetic analysis of complex traits using MCMC    *
 *                                                                          *
 *             Simon Heath - University of Washington                       *
 *                         - Rockefeller University                         *
 *                                                                          *
 *                        August 1997                                       *
 *                                                                          *
 * loki_complex_peel.c:                                                     *
 *                                                                          *
 * Perform peeling calculations                                             *
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
#include <float.h>
#ifndef DBL_MAX
#define DBL_MAX MAXDOUBLE
#endif

#define IDX_PART_BIT 8
#define IDX_PART (1<<IDX_PART_BIT)
#define HASHTABLE_SIZE 2056

#include "ranlib.h"
#include "utils.h"
#include "loki.h"
#include "loki_peel.h"

static double **complex_freq,**complex_pen,*complex_out_p;
static int max_fnd,max_all,out_p_size,n_bits1;
static int n_terms,hash_mode,*complex_mem,max_pen,max_pen1;
static size_t hb_size=2048;
static struct bin_node *hashtable[HASHTABLE_SIZE];
static struct hash_block *first_hash_block,*hash_block;

/* Given index x (from R-Function), returns corresponding n allele types in gt */
static void get_gts(lk_ulong x,const int n,int *gt)
{
	int i=0;
	lk_ulong a;
	
	a=(1<<(n_bits1))-1;
	while(x)	{
		gt[i++]=1+(int)(x&a);
		x>>=n_bits1;
	}
	for(;i<n;i++) gt[i]=1;
}

/* The inverse of get_gts() - converts the n allele types in gt into an index */
lk_ulong get_index1(int n,int *gt,const int n_bits)
{
	int i;
	lk_ulong x;
	
	if(!n) return 0;
	x=gt[n-1]-1;
	for(i=n-2;i>=0;i--) {
		x<<=n_bits;
		x|=gt[i]-1;
	}
	return x;
}

/* Returns storage for storing a new non-zero element.  Elements are allocated
 * in blocks of size hb_size */
static struct bin_node *get_new_element(const lk_ulong idx,const double p)
{
	struct bin_node *element;
	struct hash_data *hd;
	
	while(hash_block->ptr>=hash_block->size) {
		if(!hash_block->next) {
			if(!(hash_block->next=malloc(sizeof(struct hash_block)))) ABT_FUNC(MMsg);
			hash_block=hash_block->next;
			hash_block->next=0;
			if(!(hash_block->elements=malloc(sizeof(struct bin_node)*hb_size))) ABT_FUNC(MMsg);
			if(!(hash_block->hd=malloc(sizeof(struct hash_data)*hb_size))) ABT_FUNC(MMsg);
			hash_block->size=hb_size;
		} else hash_block=hash_block->next;
		hash_block->ptr=0;
	}
	hd=hash_block->hd+hash_block->ptr;
	element=hash_block->elements+hash_block->ptr++;
	element->left=element->right=0;
	element->balance=0;
	hd->index=idx;
	hd->p=p;
	element->data=hd;
	n_terms++;
	return element;
}

struct peel_mem_block *get_new_memblock(size_t size,int flag)
{
	struct peel_mem_block *p;
	
	if(!(p=malloc(sizeof(struct peel_mem_block)))) ABT_FUNC(MMsg);
	p->index=0;
	p->val=0;
	p->size=size;
	p->ptr=0;
	p->next=0;
	if(flag==MRK_MBLOCK) if(!(p->index=malloc(sizeof(lk_ulong)*size))) ABT_FUNC(MMsg);
	if(!(p->val=malloc(sizeof(double)*size))) ABT_FUNC(MMsg);
	return p;
}

void get_rf_memory(struct R_Func *rf,size_t size,int flag,struct loki *loki)
{
	size_t i;
	struct peel_mem_block *p,*p1;
	
	p1=loki->peel->mem_block[flag];
	i=p1->size-p1->ptr;
	if(size>i) {
		p=p1;
		while(p->next)	{
			if(size<=p->next->size) break;
			p=p->next;
		}
		if(!(p->next)) {
			hb_size*=1.2;
			i=(size>hb_size)?size:hb_size;
			p=get_new_memblock(i,flag);
		} else {
			p1=p->next;
			p->next=p1->next;
			p=p1;
			p->ptr=0;
			p1=loki->peel->mem_block[flag];
		}
		p->next=p1->next;
		p1->next=p;
		p1=p;
	}
	i=p1->ptr;
	rf->index=p1->index+i;
	rf->p=p1->val+i;
	p1->ptr+=size;
	loki->peel->mem_block[flag]=p1;
}

/* Insert element with index idx into the binary tree hanging off *node */
static struct bin_node *insert_node(struct bin_node *node,const lk_ulong idx,const double p,int *bal)
{
	int bb;
	lk_ulong idx1;
	struct hash_data *hd;
	
	hd=node->data;
	idx1=hd->index;
	if(idx!=idx1) {
		bb=node->balance;
		if(idx<idx1) {
			if(node->left) node->left=insert_node(node->left,idx,p,bal);
			else {
				node->left=get_new_element(idx,p);
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
			if(node->right) node->right=insert_node(node->right,idx,p,bal);
			else {
				node->right=get_new_element(idx,p);
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
	} else {
		*bal=1;
		hd->p+=p;
	}
	return node;
}

/* Recursive routine for getting all elements from binary tree hanging off *node and storing
 * the values and indices from them in *tp and *tl */
static void get_nodes(struct bin_node *node,double *tp,lk_ulong *tl,int *j)
{
	struct hash_data *hd;
	
	if(node->left)	get_nodes(node->left,tp,tl,j);
	hd=node->data;
	tp[*j]=hd->p;
	tl[(*j)++]=hd->index;
	if(node->right) get_nodes(node->right,tp,tl,j);
}

/* Similar to get_nodes, but stops when the cumulative probability >= z.  Allows sampling
 * from distribution */
static int sample_comb(struct bin_node *node,const double z,double *p,lk_ulong *idx)
{
	struct hash_data *hd;
	
	if(node->left) if(sample_comb(node->left,z,p,idx)) return 1;
	hd=node->data;
	if(hd->p>0.0) {
		*p+=hd->p;
		if(*p>=z) {
			*idx=hd->index;
			return 1;
		}
	}
	if(node->right) if(sample_comb(node->right,z,p,idx)) return 1;
	return 0;
}

/* Clean up memory used */
static void free_hash_blocks(void)
{
	struct hash_block *hb1;
	
	hash_block=first_hash_block;
	while(hash_block) {
		if(hash_block->elements) free(hash_block->elements);
		if(hash_block->hd) free(hash_block->hd);
		hb1=hash_block->next;
		free(hash_block);
		hash_block=hb1;
	}
	first_hash_block=0;
}

void free_complex_mem(void)
{
	if(complex_mem) free(complex_mem);
	if(complex_out_p) free(complex_out_p);
	if(complex_freq) {
		if(complex_freq[0]) free(complex_freq[0]);
		free(complex_freq);
	}
	if(complex_pen) {
		if(complex_pen[0]) free(complex_pen[0]);
		free(complex_pen);
	}
	free_hash_blocks();
}

static void setup_complex_peel(const struct Complex_Element *element, const int sampling,int *temp[],int *n_other,int *n_jnt,int *n_trans,struct Id_Record *id_array)
{
	int i,j,k,k1,n_peel,n_inv,n_rf,*inv,*flags;
	static int complex_mem_size;
	
	n_inv=element->n_involved;
	n_peel=element->n_peel;
	n_rf=element->n_rfuncs;
	inv=element->involved;
	flags=element->flags;
	i=11*n_inv+n_rf;
	if(!complex_mem) {
		complex_mem_size=i;
		if(!(complex_mem=calloc((size_t)i,sizeof(int)))) ABT_FUNC(MMsg);
	} else {
		if(i>complex_mem_size) {
			complex_mem_size=i;
			if(!(complex_mem=realloc(complex_mem,i*sizeof(int)))) ABT_FUNC(MMsg);
		}
		(void)memset(complex_mem,0,i*sizeof(int));
	}
	temp[0]=complex_mem;
	for(i=1;i<12;i++) temp[i]=temp[i-1]+n_inv;
	/* 'other' alleles are alleles not in R-Functions and not already sampled
	 * (if on reverse sampling pass) */
	k1=sampling?n_peel:n_inv; /* When we are sampling, the output alleles have already been sampled */
	for(*n_other=i=0;i<n_inv;i++) {
		if(flags[i]&(HAP_JNT|HAP_DAT)) {
			k=-inv[i];
			for(j=0;j<n_inv;j++) if(inv[j]==k) {
				temp[3][i]=j;
				break;
			}
		}
		if(i<k1 && !(flags[i]&(IN_RF|FIXED_FLAG))) temp[1][(*n_other)++]=i;
	}
	*n_jnt=0;
	for(i=0;i<n_inv;i++) {
		if(flags[i]&(HAP_JNT|HAP_DAT)) {
			if(inv[i]>0) {
				k=-inv[i];
				for(j=0;j<n_inv;j++) if(inv[j]==k) {
					temp[4][(*n_jnt)++]=i;
					break;
				}
			}
		}
	}
	/* Make a list of alleles that we need transmission probs for */
	*n_trans=0;
	for(i=0;i<n_inv;i++) {
		if(flags[i]&HAD_P) {
			k=id_array[-inv[i]-1].sire;
			for(j=0;j<n_inv;j++) {
				if(inv[j]==k) temp[6+X_MAT][*n_trans]=j;
				else if(-inv[j]==k) temp[6+X_PAT][*n_trans]=j;
			}
			temp[8][(*n_trans)++]=i;
		}
		if(flags[i]&HAD_M) {
			k=id_array[inv[i]-1].dam;
			for(j=0;j<n_inv;j++) {
				if(inv[j]==k) temp[6+X_MAT][*n_trans]=j;
				else if(-inv[j]==k) temp[6+X_PAT][*n_trans]=j;
			}
			temp[8][(*n_trans)++]=i;
		}
	}
	/* Form list of alleles not in R-Functions */
	if(sampling) {
		for(i=0;i<n_peel;i++) temp[0][i]=0;
		for(;i<n_inv;i++) temp[0][i]=1;
	} else for(i=0;i<n_inv;i++) temp[0][i]=0;	
}

/* Perform general peeling operation defined in element (i.e., as opposed to a 'simple' 
 * nuclear family based peeling operation) */
double loki_complex_peelop(const struct Complex_Element *element,const int locus,const int s_flag,pen_func pen,const int n_all,struct R_Func *rf,double **freq,struct loki *loki)
{
	int i,j,k,k1,k2,k3,k4,n_out,n_peel,n_inv,*inv,n_rf,n_ind,n_other,ef,ef1,sampling=0,idx_shift=0;
	int *gt_store,*gt_store1,*other_ptr,*other_list,*rf_ptr,*jnt_list,ht_size,n_trans,*trans_idx,*flags;
	int *off_index[2],*fnd_list,n_fnd=0,all,*pen_list,n_pen,*temp_p[12],n_jnt,*jnt_idx,linktype,n_idx;
	double max_terms,*tp,prob=0.0,z,p1,Konst=0.0;
	lk_ulong a,b,m,*tl,**a_set;
	struct bin_node *elem;
	struct hash_data *hd;
	struct Id_Record *id_array;
	struct Marker *mark;
	
	n_bits1=num_bits(n_all);
	n_idx=1<<(n_bits1+n_bits1);
	mark=loki->markers->marker+locus;
 	a_set=mark->all_set;
	id_array=loki->pedigree->id_array;
	linktype=loki->markers->linkage[loki->markers->marker[locus].locus.link_group].type&LINK_TYPES;
	/* Get details about peeling operation */
	n_inv=element->n_involved; /* No. alleles involved in op */
	n_peel=element->n_peel; /* No. alleles to peel out (absorb) */
	inv=element->involved; /* List of involved alleles (alleles to peel come first) */
	flags=element->flags;
	n_out=element->n_out;
	n_rf=element->n_rfuncs; /* No. input R-Functions */
	/* Allocate first hash_block, if not already done so */
	if(!first_hash_block) {
		if(!(first_hash_block=malloc(sizeof(struct hash_block)))) ABT_FUNC(MMsg);
		if(!(first_hash_block->elements=malloc(sizeof(struct bin_node)*hb_size))) ABT_FUNC(MMsg);
		if(!(first_hash_block->hd=malloc(sizeof(struct hash_data)*hb_size))) ABT_FUNC(MMsg);
		first_hash_block->next=0;
		first_hash_block->size=hb_size;
	}
	/* Reset pointers to all hash_blocks */
	hash_block=first_hash_block;
	while(hash_block) {
		hash_block->ptr=0;
		hash_block=hash_block->next;
	}
	hash_block=first_hash_block;
	n_terms=0; /* No. non-zero terms in output R-Function */
	/* if s_flag is non-zero then we are doing a sampling run */
	/* if s_flag&OP_SAMPLING then we are on the reverse (sampling) phase */
	/* In any case, if !n_out and s_flag then we can sample */
	if(s_flag) {
		if(s_flag&OP_SAMPLING) sampling=1;
		else if(!n_out) sampling=1;
	}
	/* Check size of function we are assembling.  Check index for output function will
	 * fit into a lk_long.  Compute max_terms, the number of output terms if they are all non-zero.
	 * If sampling then we assemble function on the peeled alleles */
	if(sampling) {
		k=n_bits1*n_peel;
		max_terms=log((double)n_all)*(double)n_peel;
	} else {
		k=n_bits1*n_out;
		max_terms=log((double)n_all)*(double)n_out;
	}
 	if(k>(int)LK_LONG_BIT) {
		(void)fprintf(stderr,"\nToo many individuals in output R-Function for marker %s when %s\nn_peel = %d, n_all = %d, n_bits1 = %d, required size = %d, LONG_BIT = %d\n",loki->markers->marker[locus].name,sampling?"sampling":"peeling",n_peel,n_all,n_bits1,k,LK_LONG_BIT);
		ABT_FUNC(AbMsg);
	}
	hash_mode=(log(2.0)*k<log((double)HASHTABLE_SIZE));
	if(hash_mode) ht_size=1<<k;
	else {
		ht_size=IDX_PART;
		idx_shift=k-IDX_PART_BIT;
	}
	for(i=0;i<ht_size;i++) hashtable[i]=0;
	setup_complex_peel(element,sampling,temp_p,&n_other,&n_jnt,&n_trans,id_array);
	gt_store=temp_p[0];
	other_list=temp_p[1];
	other_ptr=temp_p[2];
	jnt_list=temp_p[3];
	jnt_idx=temp_p[4];
	gt_store1=temp_p[5];
	off_index[0]=temp_p[6];
	off_index[1]=temp_p[7];
	trans_idx=temp_p[8];
	pen_list=temp_p[9];
	fnd_list=temp_p[10];
	rf_ptr=temp_p[11];
	/* Compute masks (first run through only) - used for finding mutually 
	 * consistent terms from multiple input R-Functions */
 	if(n_rf) {
		for(i=0;i<n_rf;i++) {
			j=element->index[i];
			n_ind=rf[j].n_ind;
			a=0;
			k2=0;
			b=(1<<n_bits1)-1;
			for(k=n_ind-1;k>=0;k--) {
				a<<=n_bits1;
				k1=rf[j].id_list[k];
				if(gt_store[k1]) a|=b;
				else {
					k2|=1<<k1;
					gt_store[k1]=1;
				}
			}
			rf[j].mask[sampling]=a;
			rf[j].mask1[sampling]=k2;
	  	}
	}
	/* Make a list of alleles we need penetrances and/or founder probs for */
	n_pen=0;
	/* Make list of founder alleles being peeled */
	for(k=0;k<n_inv;k++) if(flags[k]&HAP_FND) n_fnd++;
	/* Get frequency info. for founder alleles (tricky because of recoding) */
	if(n_fnd) {
		if(!complex_freq) {
			if(!(complex_freq=malloc(sizeof(void *)*n_fnd))) ABT_FUNC(MMsg);
			max_all=n_all*n_fnd;
			if(!(complex_freq[0]=malloc(sizeof(double)*max_all))) ABT_FUNC(MMsg);
			max_fnd=n_fnd;
		} else {
			if(n_fnd>max_fnd) {
				max_fnd=n_fnd;
				tp=complex_freq[0];
				if(!(complex_freq=realloc(complex_freq,sizeof(void *)*n_fnd))) ABT_FUNC(MMsg);
				complex_freq[0]=tp;
			}
			if(n_fnd*n_all>max_all) {
				max_all=n_all*n_fnd;
				if(!(complex_freq[0]=realloc(complex_freq[0],sizeof(double)*max_all))) ABT_FUNC(MMsg);
			}
		}
		for(i=1;i<n_fnd;i++) complex_freq[i]=complex_freq[i-1]+n_all;
		n_fnd=0;
		for(k=0;k<n_inv;k++) if(flags[k]&HAP_FND) {
			j=abs(inv[k])-1;
			k3=id_array[j].group-1;
#ifdef DEBUG
			if(k3<0) ABT_FUNC("Internal error - bad group number\n");
#endif
			k2=n_all-1;
			a=mark->req_set[inv[k]<0?X_PAT:X_MAT][j];
			complex_freq[n_fnd][k2]=0.0;
			for(k1=0;k1<n_all;k1++) {
				if(a&(1<<k1)) complex_freq[n_fnd][k2]+=freq[k3][k1];
				else complex_freq[n_fnd][k1]=freq[k3][k1];
			}
			fnd_list[n_fnd++]=k;
		}
	}
	if(pen) {
		/* Find alleles which we need penetrances for */
		for(k=0;k<n_peel;k++) {
			k1=abs(inv[k])-1;
			if(linktype==LINK_AUTO || (linktype==LINK_X && id_array[k1].sex==2)) {
				for(k2=k+1;k2<n_inv;k2++) if(inv[k2]== -inv[k]) break;
				if(k2==n_inv) continue;
				if(id_array[k1].res[0]) pen_list[n_pen++]=k1;
			} else if(linktype==LINK_X && id_array[k1].res[0]) pen_list[n_pen++]=k1;
		}
		/* Pre-calculate penetrances */
		if(n_pen) {
			if(!complex_pen) {
				if(!(complex_pen=malloc(sizeof(void *)*n_pen))) ABT_FUNC(MMsg);
				max_pen1=n_pen*n_idx;
				if(!(complex_pen[0]=malloc(sizeof(double)*max_pen1))) ABT_FUNC(MMsg);
				max_pen=n_pen;
			} else {
				if(n_pen>max_pen)	{
					max_pen=n_pen;
					if(!(complex_pen=realloc(complex_pen,sizeof(void *)*n_pen))) ABT_FUNC(MMsg);
				}
				if(n_pen*n_idx>max_pen1) {
					max_pen1=n_pen*n_idx;
					if(!(complex_pen[0]=realloc(complex_pen[0],sizeof(double)*max_pen1))) ABT_FUNC(MMsg);
				}
			}
			for(k1=0;k1<n_pen*n_idx;k1++) complex_pen[0][k1]=1.0;
			for(k1=1;k1<n_pen;k1++) complex_pen[k1]=complex_pen[k1-1]+n_idx;
			for(k=0;k<n_pen;k++)	{
				pen(complex_pen[k],pen_list[k],&mark->locus,n_all,n_bits1,loki);
				z=0.0;
				for(k1=0;k1<n_idx;k1++) z+=complex_pen[k][k1];
				for(k1=0;k1<n_idx;k1++) complex_pen[k][k1]/=z;
				Konst+=log(z);
			}
		}
	}
	for(i=0;i<n_inv;i++) gt_store[i]=0;
	ef1=0;
	i=0;
	/* Main loop.
	 *
	 * (1) If sampling then fix allele types of already sampled alleles.
	 * (2) Find mutually conistent configurations from input R-Functions
	 * (3) Find consistent configurations for 'other' alleles
	 * (4) Compute transmission/founder/penetrance probs
	 * (5) If still possible then add to output function
	 * (6) Cycle until tired
	 */
	for(;;) {
		/* Fix allele types of sampled or fixed alleles */
		for(j=0;j<n_inv;j++) {
			if((sampling && j>=n_peel) || (flags[j]&FIXED_FLAG)) {
				k=inv[j];
				if(k>0) gt_store[j]=id_array[k-1].allele[X_MAT];
				else gt_store[j]=id_array[-k-1].allele[X_PAT];
			}
		}
		/* Cycle through R-Functions */
		for(;i<n_rf;i++) {
			j=element->index[i];
			ef=0;
			n_ind=rf[j].n_ind;
			/* Find term consistent with what's already been set up */
			a=rf[j].mask[sampling];
			if(a)	{
				b=0;
				for(k=n_ind-1;k>=0;k--)	{
					b<<=n_bits1;
					k1=rf[j].id_list[k];
					if(gt_store[k1]) b|=gt_store[k1]-1;
				}
				b&=a;
				k1=rf[j].n_terms;
				tl=rf[j].index;
				k=rf_ptr[i];
				for(;k<k1;k++) if((tl[k]&a)==b) break;
				rf_ptr[i]=k;
			}
			/* If consistent term found then get corresponding allele types */
			if(rf_ptr[i]<rf[j].n_terms) {
				get_gts(rf[j].index[rf_ptr[i]],n_ind,gt_store1);
				for(k=0;k<n_ind;k++)	{
					k1=rf[j].id_list[k];
					gt_store[k1]=gt_store1[k];
				}
			} else { /* Inconsistency found - try new combination */
				k2=rf[j].mask1[sampling];
				k=0;
				while(k2) {
					if(k2&1) gt_store[k]=0;
					k++;
					k2>>=1;
				}
				do {
					if(!i) {
						ef1=1;
						break;
					}
					rf_ptr[i--]=0;
					rf_ptr[i]++;
					j=element->index[i];
					k2=rf[j].mask1[sampling];
					k=0;
					while(k2) {
						if(k2&1) gt_store[k]=0;
						k++;
						k2>>=1;
					}
				} while(rf_ptr[i]>=rf[j].n_terms);
				i--;
			}
			if(ef1) break;
		}
		if(ef1) break;
		/* We have a consistent configuration from the R-Functions, now find 1 for the 'other' alleles */
		for(ef=i=0;i<=n_other;i++)	{
			if(i==n_other)	{
				if(!ef) {
					for(k=0;k<n_jnt;k++)	{
						k1=jnt_idx[k];
						k2=jnt_list[k1];
						j=abs(inv[k1])-1;
						if(!(a_set[j][gt_store[k1]-1]&(1<<(gt_store[k2]-1)))) {
							ef=1;
							break;
						}
					}
					/* Add Parent-off transmission probs. */
					p1=1.0;
					if(!ef) for(k=0;k<n_trans;k++) {
						k1=trans_idx[k];
						j=inv[k1];
						if(j<0) {
							j=-j-1;
							tp=id_array[j].tpp[X_PAT];
							m=1<<(gt_store[k1]-1);
							a=mark->req_set[X_PAT][j];
							if(a&m) m=a;
							z=0.0;
							all=gt_store[off_index[X_MAT][k]];
							if(m&(1<<(all-1))) z=tp[X_MAT];
							all=gt_store[off_index[X_PAT][k]];
							if(m&(1<<(all-1))) z+=tp[X_PAT];
							p1*=z;
							if(p1==0.0)	{
								ef=1;
								break;
							}
						} else {
							j--;
							tp=id_array[j].tpp[X_MAT];
							m=1<<(gt_store[k1]-1);
							a=mark->req_set[X_MAT][j];
							if(a&m) m=a;
							z=0.0;
							all=gt_store[off_index[X_MAT][k]];
							if(m&(1<<(all-1))) z=tp[X_MAT];
							all=gt_store[off_index[X_PAT][k]];
							if(m&(1<<(all-1))) z+=tp[X_PAT];
							p1*=z;
							if(p1==0.0)	{
								ef=1;
								break;
							}
						}
					}
#ifdef DEBUG
					if(p1<0.0) {
						fprintf(stderr,"Internal error - p1 = %g\n",p1);
						ABT_FUNC("Aborting\n");
					}
#endif
					/* Add penetrances and founder probs and store result */
					if(!ef) {
						for(k=0;k<n_rf;k++) p1*=rf[element->index[k]].p[rf_ptr[k]];
#ifdef DEBUG
						if(p1<0.0) {
							fprintf(stderr,"Internal error - p1 = %g\n",p1);
							ABT_FUNC("Aborting\n");
						}
#endif
						for(k1=0;k1<n_fnd;k1++) { /* Add founder frequencies */
							k=fnd_list[k1];
							p1*=complex_freq[k1][gt_store[k]-1];
						}
#ifdef DEBUG
						if(p1<0.0) {
							fprintf(stderr,"Internal error - p1 = %g\n",p1);
							ABT_FUNC("Aborting\n");
						}
#endif
						/* Add penetrances for those being peeled if required */
 						if(n_pen) {
							for(k=0;k<n_inv;k++)	{
								k1=inv[k];
								if(k1<0) id_array[-k1-1].allele[X_PAT]=gt_store[k];
								else id_array[k1-1].allele[X_MAT]=gt_store[k];
							}
							for(k=0;k<n_pen;k++)	{
								k1=pen_list[k];
								k2=((id_array[k1].allele[X_PAT]-1)<<n_bits1)|(id_array[k1].allele[X_MAT]-1);
								p1*=complex_pen[k][k2];
							}
						}
#ifdef DEBUG
						if(p1<0.0) {
							fprintf(stderr,"Internal error - p1 = %g\n",p1);
							ABT_FUNC("Aborting\n");
						}
#endif
						prob+=p1;
						if(sampling) a=get_index1(n_peel,gt_store,n_bits1);
						else a=get_index1(n_out,gt_store+n_inv-n_out,n_bits1);
						if(hash_mode) {
							k=(int)a;
							if(!hashtable[k]) hashtable[k]=get_new_element(a,p1);
							else {
								hd=hashtable[k]->data;
								hd->p+=p1;
							}
						} else {
							k=(a>>idx_shift);
							if(hashtable[k]) {
								hashtable[k]=insert_node(hashtable[k],a,p1,&k1);
							} else hashtable[k]=get_new_element(a,p1);
						}
						ef=1;
					}
				}
			}
			/* Find next combination */
			if(ef) {
				i--;
				while(i>=0) {
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
			j=abs(inv[k1])-1;
			if(flags[k1]&(HAP_JNT|HAP_DAT)) {
				k2=jnt_list[k1];
				if(inv[k1]>0) {
					if(gt_store[k2]) {
						b=1<<(gt_store[k2]-1);
						for(k=other_ptr[i];k<n_all;k++) if(a_set[j][k]&b) {
							gt_store[k1]=k+1;
							break;
						}
						if((other_ptr[i]=k)==n_all) ef=1;
					} else {
						for(k=0;k<n_other;k++) if(other_list[k]==k2) break;
#ifdef DEBUG
						if(k==n_other) ABT_FUNC("Internal error - no other?\n");
#endif
						if(i<k) {
							for(k3=other_ptr[i];k3<n_all;k3++) {
								a=a_set[j][k3];
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
								if(mark->temp[X_PAT][j]&a)	{
									for(k3=other_ptr[i];k3<n_all;k3++) if(a_set[j][k3]&a) {
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
						b=a_set[j][gt_store[k2]-1];
						for(k=other_ptr[i];k<n_all;k++) if(b&(1<<k))	{
							gt_store[k1]=k+1;
							break;
						}
						if((other_ptr[i]=k)==n_all) ef=1;
					} else {
						for(k=i+1;k<n_other;k++) if(other_list[k]==k2) break;
#ifdef DEBUG
						if(k==n_other) ABT_FUNC("Internal error - why for this happen?\n");
#endif
						if(i<k) {
							for(k3=other_ptr[i];k3<n_all;k3++) {
								a=1<<k3;
								if(mark->temp[X_PAT][j]&a)	{
									for(k4=other_ptr[k];k4<n_all;k4++) if(a_set[j][k4]&a)	{
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
								a=a_set[j][k4];
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
					for(k=other_ptr[i];k<n_all;k++) if(a_set[j][k]) {
						gt_store[k1]=k+1;
						break;
					}
					if((other_ptr[i]=k)==n_all) ef=1;
				} else {
					a=mark->temp[X_PAT][j];
					for(k=other_ptr[i];k<n_all;k++) if(a&(1<<k)) {
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
			do {
				rf_ptr[i]++;
				j=element->index[i];
				k2=rf[j].mask1[sampling];
				if(k2) {
					k=0;
					while(k2) {
						if(k2&1) gt_store[k]=0;
						k++;
						k2>>=1;
					}
					if(rf_ptr[i]<rf[j].n_terms) break;
				}
				rf_ptr[i--]=0;
			} while(i>=0);
			if(i>=0) ef=0;
		} else i=0;
		if(ef) break;
	}
	/* All valid combinations have been visited.  Have we found any ? */
	if(!n_terms) ABT_FUNC("Zero probability!\n");
#ifdef DEBUG
	if(prob<=0.0) {
		fprintf(stderr,"Error - probability = %g\n",prob);
		ABT_FUNC("Aborting\n");
	}
#endif
	/* If sampling then sample from output function */
	if(sampling) {
		if(n_peel) {
			do {
				z=safe_ranf()*prob;
				p1=0.0;
				for(j=k=0;k<ht_size;k++) {
					elem=hashtable[k];
					if(!elem) continue;
					if(hash_mode) {
						hd=elem->data;
						if(hd->p>0.0) {
							p1+=hd->p;
							if(z<=p1) {
								a=(lk_ulong)k;
								break;
							}
						}
					} else if(sample_comb(elem,z,&p1,&a)) break;
				}
#ifdef DEBUG
				if(k==ht_size) {
					ABT_FUNC("Internal error\n");
				}
#endif
			} while(k==ht_size);
			get_gts(a,n_peel,gt_store);
			for(i=0;i<n_peel;i++) {
				j=inv[i];
				if(j>0) {
					id_array[j-1].allele[X_MAT]=gt_store[i];
					id_array[j-1].flag|=SAMPLED_MAT;
				} else {
					id_array[-j-1].allele[X_PAT]=gt_store[i];
					id_array[-j-1].flag|=SAMPLED_PAT;
				}
			}
		} else for(i=0;i<n_inv;i++) {
			j=inv[i];
			k=abs(j)-1;
			id_array[k].flag|=(j<0?SAMPLED_PAT:SAMPLED_MAT);
		}
	} else { /* Otherwise normalize and store */
		i=element->out_index;
		if(i>=0)	{
			rf[i].n_terms=n_terms;
			get_rf_memory(rf+i,n_terms,MRK_MBLOCK,loki);
			tl=rf[i].index;
			tp=rf[i].p;
			if(hash_mode) {
				for(j=k=0;k<ht_size;k++) {
					elem=hashtable[k];
					if(!elem) continue;
					hd=elem->data;
					tp[j]=hd->p;
					tl[j++]=(lk_ulong)k;
				}
			} else {
				for(j=k=0;k<ht_size;k++) {
					elem=hashtable[k];
					if(elem) get_nodes(elem,tp,tl,&j);
				}
			}
			for(j=0;j<n_terms;j++) tp[j]/=prob;
		}
	}
	return log(prob)+Konst;
}

/* Similar to loki_complex_peelop, but for a trait locus */
double loki_trait_complex_peelop(const struct Complex_Element *element,const int locus,const int s_flag,struct R_Func *rf,trait_pen_func *trait_pen,double **freq,struct loki *loki)
{
	int i,j,k,k1,k2,k3,n_out,n_peel,n_inv,*inv,n_rf,n_ind,n_other,ef,ef1,sampling=0;
	int *gt_store,*other_ptr,*other_list,*rf_ptr,n_idx,n_all,*flags;
	int *off_index[2],*fnd_list,n_fnd=0,all,*pen_list,n_pen,*temp_p[12],n_trans,*trans_idx;
	double max_terms,*tp,prob=0.0,z,p1,Konst=0.0,*p_rf;
	struct Id_Record *id_array;
	struct Locus *loc;
	
	id_array=loki->pedigree->id_array;
	loc=&loki->models->tlocus[-1-locus];
	n_all=loc->n_alleles;
	n_bits1=num_bits(n_all);
	n_idx=n_all*n_all;
	flags=element->flags;
	n_terms=0; /* No. non-zero terms in output R-Function */
	/* Get details about peeling operation */
	n_inv=element->n_involved; /* No. alleles involved in op */
	n_peel=element->n_peel; /* No. allles to peel out (absorb) */
	n_out=n_inv-n_peel; /* No. to appear in output R-Function */
	n_rf=element->n_rfuncs; /* No. input R-Functions */
	inv=element->involved; /* List of involved alleles (alleles to peel come first) */
	/* if s_flag is non-zero then we are doing a sampling run */
	/* if s_flag&OP_SAMPLING then we are on the reverse (sampling) phase */
	/* In any case, if !n_out and s_flag then we can sample */
	if(s_flag) {
		if(s_flag&OP_SAMPLING) sampling=1;
		else if(!n_out) sampling=1;
	}
	/* In a sampling operation we sample the peeled alleles conditional on
	 * R-Functions and on the pivot (output) alleles which have all already been sampled */
	if(sampling) {
		for(i=0;i<n_out;i++)	{
			j=inv[i+n_peel];
			if(j<0) {
				if(!(id_array[-j-1].flag&SAMPLED_PAT)) break;
			} else if(!(id_array[j-1].flag&SAMPLED_MAT)) break;
		}	
#ifdef DEBUG
		if(i<n_out) ABT_FUNC("Internal error: Unsampled pivots\n");
#endif
	}
	if(sampling) max_terms=log((double)n_all)*(double)n_peel;
	else max_terms=log((double)n_all)*(double)n_out;
	i=(int)(.5+exp(max_terms));
	if(!complex_out_p) {
		out_p_size=i+n_rf;
		if(!(complex_out_p=malloc(sizeof(double)*out_p_size))) ABT_FUNC(MMsg);
	}
	if(i+n_rf>out_p_size) {
		out_p_size=i+n_rf;
		if(!(complex_out_p=realloc(complex_out_p,sizeof(double)*out_p_size))) ABT_FUNC(MMsg);
	}
	p_rf=complex_out_p+i;
	for(j=0;j<i;j++) complex_out_p[j]=0.0;
	setup_complex_peel(element,sampling,temp_p,&n_other,&k,&n_trans,id_array);
	gt_store=temp_p[0];
	other_list=temp_p[1];
	other_ptr=temp_p[2];
	off_index[0]=temp_p[6];
	off_index[1]=temp_p[7];
	trans_idx=temp_p[8];
	pen_list=temp_p[9];
	fnd_list=temp_p[10];
	rf_ptr=temp_p[11];
	/* Compute masks - used for finding mutually 
	 * consistent terms from multiple input R-Functions */
	if(n_rf)	{
		for(i=0;i<n_rf;i++) {
			j=element->index[i];
			n_ind=rf[j].n_ind;
			k2=0;
			k3=1;
			for(k=n_ind-1;k>=0;k--)	{
				k1=rf[j].id_list[k];
				if(!gt_store[k1])	{
					k2|=1<<k1;
					k3*=n_all;
					gt_store[k1]=1;
				}
			}
			rf[j].mask1[sampling]=k2;
			rf[j].n_terms=k3;
		}
	}
	/* Make a list of alleles we need penetrances and/or founder probs for */
	n_pen=0;
	for(k=0;k<n_inv;k++) if(flags[k]&HAP_FND) n_fnd++;
	if(n_fnd) {
		n_fnd=0;
		for(k=0;k<n_inv;k++) if(flags[k]&HAP_FND) {
			j=abs(inv[k])-1;
			fnd_list[n_fnd++]=k;
		}
	}
	for(k=0;k<n_peel;k++) {
		for(k1=k+1;k1<n_inv;k1++) if(inv[k1]== -inv[k]) break;
		if(k1==n_inv) continue;
		k1=abs(inv[k])-1;
		if(id_array[k1].res[0]) pen_list[n_pen++]=k1;
	}
	if(n_pen) {
		if(!complex_pen) {
			if(!(complex_pen=malloc(sizeof(void *)*n_pen))) ABT_FUNC(MMsg);
			max_pen1=n_pen*n_idx;
			if(!(complex_pen[0]=malloc(sizeof(double)*max_pen1))) ABT_FUNC(MMsg);
			max_pen=n_pen;
		} else {
			if(n_pen>max_pen)	{
				max_pen=n_pen;
				if(!(complex_pen=realloc(complex_pen,sizeof(void *)*n_pen))) ABT_FUNC(MMsg);
			}
			if(n_pen*n_idx>max_pen1) {
				max_pen1=n_pen*n_idx;
				if(!(complex_pen[0]=realloc(complex_pen[0],sizeof(double)*max_pen1))) ABT_FUNC(MMsg);
			}
		}
		for(k1=0;k1<n_pen*n_idx;k1++) complex_pen[0][k1]=1.0;
		for(k1=1;k1<n_pen;k1++) complex_pen[k1]=complex_pen[k1-1]+n_idx;
		for(k=0;k<n_pen;k++) {
			trait_pen(complex_pen[k],pen_list[k],loc,loki);
			z=0.0;
			for(k1=0;k1<n_idx;k1++) z+=complex_pen[k][k1];
			if(z<=0.0) {
				if(!(s_flag&1)) return -DBL_MAX;
				ABT_FUNC("Zero probability!\n");
			}
			for(k1=0;k1<n_idx;k1++) complex_pen[k][k1]/=z;
			Konst+=log(z);
		}
	}
	for(i=0;i<n_inv;i++) gt_store[i]=0;
	ef1=0;
	i=0;
	/* Main loop. */
	for(;;) {
		/* Fix allele types of sampled alleles */
		if(sampling) for(j=0;j<n_out;j++) {
			k=inv[j+n_peel];
			if(k>0) gt_store[j+n_peel]=id_array[k-1].allele[X_MAT];
			else gt_store[j+n_peel]=id_array[-k-1].allele[X_PAT];
		}
		/* Cycle through R-Functions */
		for(;i<n_rf;i++) {
			j=element->index[i];
			ef=0;
			n_ind=rf[j].n_ind;
			/* Find term consistent with what's already been set up */
			k3=rf_ptr[i];
			k2=0;
			for(k=n_ind-1;k>=0;k--)	{
				k1=rf[j].id_list[k];
				if(!gt_store[k1])	{
					gt_store[k1]=1+(k3%n_all);
					k3/=n_all;
				}
				k2=k2*n_all+gt_store[k1]-1;
			}
			p_rf[i]=rf[j].p[k2];
			if(p_rf[i]<=0.0) {/* Inconsistent - find the next 1 */
				for(;;) {
					k2=rf[j].mask1[sampling];
					k=0;
					while(k2) {
						if(k2&1) gt_store[k]=0;
						k++;
						k2>>=1;
					}
					if(++rf_ptr[i]<rf[j].n_terms) break;
					if(!i) {
						ef1=1;
						break;
					}
					rf_ptr[i--]=0;
					j=element->index[i];
				}
				i--;
			}
			if(ef1) break;
		}
		if(ef1) break;
		for(ef=i=0;i<=n_other;i++) {
			if(i==n_other) {
				if(!ef) {
					/* Add Parent-off transmission probs. */
					p1=1.0;
					if(!ef) for(k=0;k<n_trans;k++) {
						k1=trans_idx[k];
						j=inv[k1];
						if(j<0) {
							j=-j-1;
							tp=id_array[j].tpp[X_PAT];
							k3=gt_store[k1];
							z=0.0;
							all=gt_store[off_index[X_MAT][k]];
							if(k3==all) z=tp[X_MAT];
							all=gt_store[off_index[X_PAT][k]];
							if(k3==all) z+=tp[X_PAT];
							p1*=z;
							if(p1==0.0)	{
								ef=1;
								break;
							}
						} else {
							j--;
							tp=id_array[j].tpp[X_MAT];
							k3=gt_store[k1];
							z=0.0;
							all=gt_store[off_index[X_MAT][k]];
							if(k3==all) z=tp[X_MAT];
							all=gt_store[off_index[X_PAT][k]];
							if(k3==all) z+=tp[X_PAT];
							p1*=z;
							if(p1==0.0)	{
								ef=1;
								break;
							}
						}
					}
					/* Add penetrances and founder probs and store result */
					if(!ef) {
						for(k=0;k<n_rf;k++) p1*=p_rf[k];
						for(k1=0;k1<n_fnd;k1++) { /* Add founder frequencies */
							k=fnd_list[k1];
							j=abs(inv[k])-1;
							k2=id_array[j].group-1;
#ifdef DEBUG
							if(k2<0) ABT_FUNC("OOOK!\n");
#endif
							p1*=freq[k2][gt_store[k]-1];
						}
						/* Add penetrances for those being peeled if required */
 						if(n_pen) {
							for(k=0;k<n_inv;k++)	{
								k1=inv[k];
								if(k1<0) id_array[-k1-1].allele[X_PAT]=gt_store[k];
								else id_array[k1-1].allele[X_MAT]=gt_store[k];
							}
							for(k=0;k<n_pen;k++)	{
								k1=pen_list[k];
								k2=((id_array[k1].allele[X_PAT]-1)<<n_bits1)|(id_array[k1].allele[X_MAT]-1);
								p1*=complex_pen[k][k2];
							}
						}
						prob+=p1;
						k1=0;
						if(sampling) for(k=n_peel-1;k>=0;k--) k1=k1*n_all+gt_store[k]-1;
						else for(k=n_inv-1;k>=n_peel;k--) k1=k1*n_all+gt_store[k]-1;
						complex_out_p[k1]+=p1;
						k1=0;
						for(k=n_inv-1;k>=0;k--) k1=k1*n_all+gt_store[k]-1;
						ef=1;
					}
				}
			}
			/* Find next combination */
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
			} else gt_store[other_list[i]]=other_ptr[i]+1;
		}
		for(i=0;i<n_other;i++) gt_store[other_list[i]]=0;
		/* Cycle through R-Functions */
		if(n_rf)	{
			i=n_rf-1;
			do	{
				rf_ptr[i]++;
				j=element->index[i];
				k2=rf[j].mask1[sampling];
				if(k2) {
					k=0;
					while(k2) {
						if(k2&1) gt_store[k]=0;
						k++;
						k2>>=1;
					}
					if(rf_ptr[i]<rf[j].n_terms) break;
				}
				rf_ptr[i--]=0;
			} while(i>=0);
			if(i>=0) ef=0;
		} else i=0;
		if(ef) break;
	}
	/* All valid combinations have been visited.  Have we found any ? */
	if(prob<=0.0) {
		if(!(s_flag&1)) return -DBL_MAX;
		ABT_FUNC("Zero probability!\n");
	}
	n_terms=(int)(.5+exp(max_terms));
	/* If sampling then sample from output function */
	if(sampling) {
		do {
			z=ranf()*prob;
			p1=0.0;
			for(k=0;k<n_terms;k++) {
				if(complex_out_p[k]>0.0) {
					p1+=complex_out_p[k];
					if(z<=p1) break;
				}
			}
		} while(k==n_terms);
		for(k1=0;k1<n_peel;k1++) {
			gt_store[k1]=(k%n_all)+1;
			k/=n_all;
		}
		for(i=0;i<n_peel;i++) {
			j=inv[i];
			if(j>0) {
				id_array[j-1].allele[X_MAT]=gt_store[i];
				id_array[j-1].flag|=SAMPLED_MAT;
			} else {
				id_array[-j-1].allele[X_PAT]=gt_store[i];
				id_array[-j-1].flag|=SAMPLED_PAT;
			}
		}
	} else { /* Otherwise normalize and store */
		i=element->out_index;
		if(i>=0)	{
			get_rf_memory(rf+i,n_terms,TRT_MBLOCK,loki);
			for(j=0;j<n_terms;j++) rf[i].p[j]=complex_out_p[j]/prob;
		}
	}
	return log(prob)+Konst;
}

