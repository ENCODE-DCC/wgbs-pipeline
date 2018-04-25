/****************************************************************************
 *                                                                          *
 *     Loki - Programs for genetic analysis of complex traits using MCMC    *
 *                                                                          *
 *             Simon Heath - University of Washington                       *
 *                                                                          *
 *                       March 1997                                         *
 *                                                                          *
 * recode_factors.c:                                                        *
 *                                                                          *
 * Recodes factorial data.                                                  *
 *                                                                          *
 * Copyright (C) Simon C. Heath 1997, 2000, 2002                            *
 * This is free software.  You can distribute it and/or modify it           *
 * under the terms of the Modified BSD license, see the file COPYING        *
 *                                                                          *
 ****************************************************************************/

#include <config.h>
#include <stdlib.h>
#include <strings.h>
#ifdef USE_DMALLOC
#include <dmalloc.h>
#endif
#include <math.h>
#include <stdio.h>

#include "utils.h"
#include "scan.h"
#include "y.tab.h"

struct var_element **var_factors=0;
int n_factors;
static int factor,*tmp_rec;

static void check_for_factor(struct bin_node *node,int *i)
{
	int j;
	struct var_element *elem;
	struct scan_data *sd;
	
	sd=node->data;
	for(j=0;j<sd->n_elements;j++) {
		elem=sd->element+j;
		if((elem->type&ST_FACTOR) && !(elem->type&(ST_HAPLO|ST_MARKER|ST_TRAITLOCUS|ST_ID|ST_SIRE|ST_DAM|ST_FAMILY))) (*i)++;
	}
}

static void list_factors(struct bin_node *node,int *i)
{
	int j;
	struct var_element *elem;
	struct scan_data *sd;
	
	sd=node->data;
	for(j=0;j<sd->n_elements;j++) {
		elem=sd->element+j;
		if((elem->type&ST_FACTOR) && !(elem->type&(ST_HAPLO|ST_MARKER|ST_TRAITLOCUS|ST_ID|ST_SIRE|ST_DAM|ST_FAMILY)))
		  var_factors[(*i)++]=elem;
	}
}

static int tmp_cmp(const void *s1,const void *s2)
{
	int i1,i2;
	
	i1=tmp_rec[*(const int *)s1];
	i2=tmp_rec[*(const int *)s2];
	if(i1<i2) return 1;
	if(i2<i1) return -1;
	return 0;
}

static int cmp_factors(const void *s1,const void *s2)
{
	const struct label_data *n1,*n2;
	int i=0;
	
	n1=factor_recode[factor][*((const int *)s1)];
	n2=factor_recode[factor][*((const int *)s2)];
	
	switch(n1->type) {
	 case STRING:
		i=strcasecmp(n1->data.string,n2->data.string);
		break;
	 case INTEGER:
		i= n1->data.value - n2->data.value;
		break;
	 default:
		ABT_FUNC("Internal error - invalid type\n");
	}
	return i;
}

void recode_factors(int num_nodes,struct recode_table_tag *recode_table)
{
	int i,j,k,ncol,col,ind,*trans,*perm;
	struct InFile *infile;
	struct DataBlock *db;
	struct var_element *elem;
	
	n_factors=0;
	check_vars(root_var,&n_factors,check_for_factor);
	if(n_factors+n_markers) if(!(factor_recode=calloc((size_t)(n_factors+n_markers),sizeof(struct label_node **)))) ABT_FUNC(MMsg);
	if(!n_factors) return;
	(void)fputs("Recoding factors\n",stdout);
	if(!(var_factors=malloc(sizeof(void *)*n_factors))) ABT_FUNC(MMsg);
	n_factors=0;
	check_vars(root_var,&n_factors,list_factors);
	for(factor=0;factor<n_factors;factor++) {
		k=0;
		infile=Infiles;
		elem=var_factors[factor];
		while(infile) {
			if(infile->n_records) {
				ncol=infile->ncol;
				for(col=i=0;i<infile->nvar;i++) if(infile->element[i]) {
					if(elem==infile->element[i]) {
						db=infile->data;
						while(db) {
							for(j=0;j<db->record_ptr;j++) if(!check_missing(col,ncol,j,db)) {
								ind=db->records[j*ncol+col].node->index;
								if(!recode_table[ind].index) {
									recode_table[ind].index= ++k;
									recode_table[ind].node=db->records[j*ncol+col].node;
								}
								db->records[j*ncol+col].value=recode_table[ind].index;
							}
							db=db->next;
						}
						break;
					}
					col++;
				}
			}
			infile=infile->next;
		}
		elem->n_levels=k;
 		elem->index=factor;
		if(k)	{
			if(!(factor_recode[factor]=malloc(sizeof(struct label_node *)*k))) ABT_FUNC(MMsg);
			if(!(perm=malloc(sizeof(int)*k*2))) ABT_FUNC(MMsg);
			trans=perm+k;
			for(i=0;i<num_nodes;i++)	{
				j=recode_table[i].index;
				if(j) factor_recode[factor][j-1]=recode_table[i].node;
			}
			for(i=0;i<k;i++) perm[i]=i;
			gnu_qsort(perm,(size_t)k,sizeof(int),cmp_factors);
			for(i=0;i<k;i++) trans[perm[i]]=i;
			for(i=0;i<num_nodes;i++)	{
				j=recode_table[i].index;
				if(j)	{
					j=trans[j-1];
					factor_recode[factor][j]=recode_table[i].node;
					recode_table[i].index=0;
				}
			}
			infile=Infiles;
			while(infile) {
				if(infile->n_records) {
					ncol=infile->ncol;
					for(col=i=0;i<infile->nvar;i++) if(infile->element[i]) {
						if(elem==infile->element[i]) {
							db=infile->data;
							while(db) {
								for(j=0;j<db->record_ptr;j++) if(!check_missing(col,ncol,j,db)) {
									k=db->records[j*ncol+col].value;
									db->records[j*ncol+col].value=trans[k-1]+1;
								}
								db=db->next;
							}
							break;
						}
						col++;
					}
				}
				infile=infile->next;
			}
			free(perm);
		}
	}
}

void recode_marker_data(int num_nodes,struct recode_table_tag *recode_table)
{
	int i,j,k,hap[2],ncol,col,marker,nv,nvar,ind,id,haplo[2],haplo1[2],col_list[2],*tmp_perm,gtflag;
	struct InFile *infile;
	struct DataBlock *db;
	struct var_element *elem,*elem1,**elem_list;
	struct gt_data *gt;
	
	(void)fputs("Recoding marker data\n",stdout);
	for(marker=0;marker<n_markers;marker++) {
		k=0;
		infile=Infiles;
		elem=markers[marker].element;
		if(elem->type&ST_DATA) {
			elem1=elem;
			gtflag=1;
		} else {
			elem=markers[marker].hap_element[0];
			elem1=markers[marker].hap_element[1];
			gtflag=0;
		}
		if(!(elem || elem1)) continue;
		while(infile) {
			if(infile->n_records) {
				nvar=infile->nvar;
				elem_list=infile->element;
				for(nv=col=i=0;i<nvar;i++) if(elem_list[i])	{
					if(elem==elem_list[i]) {
						if(nv<2) {
							hap[nv]=0;
							col_list[nv++]=col;
						} else break;
					} else if(elem1==elem_list[i]) {
						if(nv<2) {
							hap[nv]=1;
							col_list[nv++]=col;
						} else break;
					}
					col++;
				}
				if(nv>2-gtflag) ABT_FUNC("Internal error: variable occurs too often in file list\n");
				if(nv) {
					ncol=infile->ncol;
					db=infile->data;
					while(db) {
						for(j=0;j<db->record_ptr;j++)	{
							id=(int)db->records[j*ncol+infile->id_col].value;
							if(!id) ABT_FUNC("Internal error - null id\n");
							if(id_array[id-1].haplo[0]) {
								for(i=0;i<2;i++) {
									haplo[i]=id_array[id-1].haplo[i][marker];
									id_array[id-1].haplo[i][marker]=0;
								}
							} else haplo[0]=haplo[1]=0;
							for(i=0;i<nv;i++) {
								col=col_list[i];
								if(!check_missing(col,ncol,j,db)) {
									if(gtflag) {
										gt=db->records[j*ncol+col].gt_data;
										if(gt->node1) {
											ind=gt->node1->index;
											if(!recode_table[ind].index) {
												recode_table[ind].index= ++k;
												recode_table[ind].node=gt->node1;
											}
											if(!id_array[id-1].haplo[0]) {
												if(!(id_array[id-1].haplo[0]=calloc(2*(size_t)n_markers,sizeof(int)))) ABT_FUNC(MMsg);
												id_array[id-1].haplo[1]=id_array[id-1].haplo[0]+n_markers;
											}
											id_array[id-1].haplo[0][marker]=recode_table[ind].index;
										}
										if(gt->node2) {
											ind=gt->node2->index;
											if(!recode_table[ind].index) {
												recode_table[ind].index= ++k;
												recode_table[ind].node=gt->node2;
											}
											if(!id_array[id-1].haplo[0]) {
												if(!(id_array[id-1].haplo[0]=calloc(2*(size_t)n_markers,sizeof(int)))) ABT_FUNC(MMsg);
												id_array[id-1].haplo[1]=id_array[id-1].haplo[0]+n_markers;
											}
											id_array[id-1].haplo[1][marker]=recode_table[ind].index;
										}
										free(gt);
										gt=db->records[j*ncol+col].gt_data=0;
									} else {
										ind=db->records[j*ncol+col].node->index;
										if(!recode_table[ind].index) {
											recode_table[ind].index= ++k;
											recode_table[ind].node=db->records[j*ncol+col].node;
										}
										db->records[j*ncol+col].node=0;
										if(!id_array[id-1].haplo[0]) {
											if(!(id_array[id-1].haplo[0]=calloc(2*(size_t)n_markers,sizeof(int)))) ABT_FUNC(MMsg);
											id_array[id-1].haplo[1]=id_array[id-1].haplo[0]+n_markers;
										}
										id_array[id-1].haplo[hap[i]][marker]=recode_table[ind].index;
									}
								}
							}
							if(id_array[id-1].haplo[0]) {
								for(i=0;i<2;i++) haplo1[i]=id_array[id-1].haplo[i][marker];
								for(i=0;i<2;i++) if(haplo[i] && haplo1[i] && haplo[i]!=haplo1[i]) break;
								if(i<2) for(i=0;i<2;i++) if(haplo[i] && haplo1[i^1] && haplo[i]!=haplo1[i^1]) break;
								if(i<2) {
									(void)fprintf(stderr,"Marker %s",markers[marker].var->name);
									if(markers[marker].index) (void)fprintf(stderr,"(%d)",markers[marker].index);
									(void)fputs(" - Error: individual ",stderr);
									print_orig_id(stderr,id,0);
									print_scan_err(" has inconsistent genotype information\n");
								} else {
									for(i=0;i<2;i++) id_array[id-1].haplo[i][marker]=haplo[i]?haplo[i]:haplo1[i];
								}
							}
						}
						db=db->next;
					}
				}
			}
			infile=infile->next;
		}
		markers[marker].element->index=marker;
		markers[marker].element->n_levels=k;
		if(k) {
			if(!(tmp_rec=malloc(sizeof(int)*k*2))) ABT_FUNC(MMsg);
			tmp_perm=tmp_rec+k;
			for(i=0;i<k;i++) {
				tmp_rec[i]=0;
				tmp_perm[i]=i;
			}
			for(id=0;id<ped_size;id++) {
				if(id_array[id].haplo[0]) {
					for(i=0;i<2;i++) {
						j=id_array[id].haplo[i][marker];
						if(j) tmp_rec[j-1]++;
					}
				}
			}
			gnu_qsort(tmp_perm,(size_t)k,sizeof(int),tmp_cmp);
			for(i=0;i<k;i++) {
				j=tmp_perm[i];
				tmp_rec[j]=i;
			}
			for(id=0;id<ped_size;id++) {
				if(id_array[id].haplo[0]) {
					for(i=0;i<2;i++) {
						j=id_array[id].haplo[i][marker];
						if(j) id_array[id].haplo[i][marker]=tmp_rec[j-1]+1;
					}
				}
			}
			if(!(factor_recode[n_factors+marker]=malloc(sizeof(struct label_node *)*k))) ABT_FUNC(MMsg);
			for(i=0;i<num_nodes;i++)	{
				j=recode_table[i].index;
				if(j)	{
					j=tmp_rec[j-1]+1;
					factor_recode[n_factors+marker][j-1]=recode_table[i].node;
					recode_table[i].index=0;
				}
			}
			free(tmp_rec);
		}
	}
}

