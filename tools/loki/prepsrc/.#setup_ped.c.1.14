/****************************************************************************
 *                                                                          *
 *     Loki - Programs for genetic analysis of complex traits using MCMC    *
 *                                                                          *
 *             Simon Heath - University of Washington                       *
 *                                                                          *
 *                       March 1997                                         *
 *                                                                          *
 * setup_ped.c:                                                             *
 *                                                                          *
 * Sets up pedigree structures and recodes ids                              *
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
#ifdef HAVE_LIMITS_H
#include <limits.h>
#endif
#include "utils.h"
#include "y.tab.h"
#include "scan.h"

int pruned_ped_size,n_comp,*comp_size,verbose=3,comp_sflag,n_orig_families;
int *rec_tab,*rec_tab1,n_genetic_groups=1;
struct label_data **group_recode;
struct Id_Record *id_array=0;

int print_orig_family(FILE *fptr,const int id,const int fg) 
{
	int fam,i=0;
	
	fam=id>0?id_array[id-1].fam_code:0;
	if(fg) {
		if(!fam || !family_recode[fam-1]) {
			(void)fputs("[*]:",fptr);
			i=4;
		} else if(family_recode[fam-1]->type==STRING) i=fprintf(fptr,"[%s]:",family_recode[fam-1]->data.string);
		else i=fprintf(fptr,"[%ld]:",family_recode[fam-1]->data.value);
	} else {
		if(!fam || !family_recode[fam-1]) {
			(void)fputs("*",fptr);
			i=1;
		} else if(family_recode[fam-1]->type==STRING) i=fprintf(fptr,"%s",family_recode[fam-1]->data.string);
		else i=fprintf(fptr,"%ld",family_recode[fam-1]->data.value);
	}
	return i;
}

void print_orig_id1(FILE *fptr,const int id,const int fg)
{
	if(!id) (void)fputc('*',fptr);
	else if(id<0) (void)fprintf(fptr,"BAD ID [%d]",id);
	else {
		if(!ped_recode[id-1]) (void)fputc('*',fptr);
		else if(ped_recode[id-1]->type==STRING) (void)fputs(ped_recode[id-1]->data.string,fptr);
		else (void)fprintf(fptr,"%ld",ped_recode[id-1]->data.value);
	}
	if(fg) (void)fputc(' ',fptr);
}

void print_orig_id(FILE *fptr,const int id,const int fg)
{
	if(family_id) print_orig_family(fptr,id,1);
	print_orig_id1(fptr,id,fg);
}

void print_orig_triple(FILE *fptr,const int id)
{
	if(family_id) {
		if(id) {
			print_orig_family(fptr,id,0);
			fputc(' ',fptr);
		} else (void)fputs("* ",fptr);
	}
	if(id) {
		print_orig_id1(fptr,id,1);
		print_orig_id1(fptr,id_array[id-1].sire,1);
		print_orig_id1(fptr,id_array[id-1].dam,1);
	} else (void)fputs("* * * ",fptr);
}

static void do_recode(int i,int *id)
{
	int id1;
	
	if(ped_recode1[i-1]<0) {
		(void)fprintf(stderr,"Pedigree Error - Circular pedigree found.  Ids in list follow:\n");
		print_orig_triple(stderr,i);
		(void)fputc('\n',stderr);
		*id= -1;
		return;
	}     
	ped_recode1[i-1]= -1;
	id1=id_array[i-1].sire;
	if(id1 && !ped_recode1[id1-1]) do_recode(id1,id);
	if(*id>0) {
		id1=id_array[i-1].dam;
		if(id1 && !ped_recode1[id1-1]) do_recode(id1,id);
	}
	if(*id<0) {
		print_orig_triple(stderr,i);
		(void)fputc('\n',stderr);
		return;
	}
	ped_recode1[i-1]= (*id)++;
	return;
}

static void infect(int i,int *group,int *equiv,int *ng)
{
	int j,ids,idd,i1,i2;

	ids=id_array[i].sire;
	idd=id_array[i].dam;
	if(group[i]<0) {
		(void)fprintf(stderr,"Pedigree Error - Circular pedigree found.  Ids in cycle follow:\n");
		print_orig_triple(stderr,i+1);
		(void)fputc('\n',stderr);
		*ng= -1;
		return;
	}     
	group[i]=-1;
	if(ids && group[ids-1]<1) infect(ids-1,group,equiv,ng);
	if(*ng>=0 &&idd && group[idd-1]<1) infect(idd-1,group,equiv,ng);
	if(*ng>=0) {
		if(ids && group[ids-1]>0) {
			if(idd && group[idd-1]>0)  {
				i1=equiv[group[ids-1]-1];
				i2=equiv[group[idd-1]-1];
				if(i1!=i2) for(j=0;j<(*ng);j++) if(equiv[j]==i2) equiv[j]=i1;
			}
			group[i]=group[ids-1];
		} else if(idd && group[idd-1]>0) group[i]=group[idd-1];
		else {
			equiv[*ng]= *ng;
			group[i]= ++(*ng);
		}
	} else {
		print_orig_triple(stderr,i+1);
		(void)fputc('\n',stderr);
		return;
	}
}

static int cmp_comp_size(const void *s1,const void *s2)
{
	int sz1,sz2;
	
	sz1=comp_size[*((const int *)s1)];
	sz2=comp_size[*((const int *)s2)];
	if(sz1<sz2) return 1;
	if(sz1>sz2) return -1;
	return 0;
}

void count_components(char *LogFile)
{
	int i,j,k,*group,*equiv,*code,*perm,*nbits,ids,idd;
	FILE *flog;
	
	n_comp=k=0;
	if(!(group=calloc((size_t)(3*ped_size),sizeof(int)))) ABT_FUNC(MMsg);
	equiv=group+ped_size;
	perm=equiv+ped_size;
	n_comp=0;
	for(j=i=0;i<ped_size;i++) {
		if(ped_recode1[i] && !group[i]) infect(i,group,equiv,&j);
		if(j<0) ABT_FUNC(AbMsg);
	}
	if(!j) {
		free(group);
		return;
	}
	if(!(code=calloc((size_t)j,sizeof(int)))) {
		(void)fprintf(stderr,"*** i = %d, j = %d ***\n",i,j);
		ABT_FUNC(MMsg);
	}
	for(i=0;i<j;i++) if(!code[equiv[i]]) code[equiv[i]]= ++n_comp;
	for(i=0;i<ped_size;i++) if(ped_recode1[i]) id_array[i].component=code[equiv[group[i]-1]];
	if(!(comp_size=calloc((size_t)n_comp,sizeof(int)))) ABT_FUNC(MMsg);
	for(i=0;i<ped_size;i++)	{
		j=ped_recode1[i];
		if(j)	{
			comp_size[id_array[i].component-1]++;
			equiv[j-1]=i;
		}
	}
	if(n_comp>1) {
		for(i=0;i<n_comp;i++) {
			perm[i]=i;
			code[i]=comp_size[i];
		}
		gnu_qsort(perm,(size_t)n_comp,sizeof(int),cmp_comp_size);
		for(i=0;i<n_comp;i++) comp_size[i]=code[perm[i]];
		for(i=n_comp-1;i>=0;i--) if(comp_size[i]>1) break;
		if(i<n_comp-2) {
			comp_sflag=1;
			k=i+2;
			comp_size[k-1]=1+n_comp-k;
		}
		for(i=0;i<n_comp;i++) code[perm[i]]=i;
		for(i=0;i<ped_size;i++) if(ped_recode1[i]) {
			j=code[id_array[i].component-1];
			if(comp_sflag && comp_size[j]==1) id_array[i].component=k;
			else id_array[i].component=j+1;
		}
		if(comp_sflag) n_comp=k;
		group[0]=1;
		for(i=1;i<n_comp;i++) group[i]=group[i-1]+comp_size[i-1];
		for(j=0;j<pruned_ped_size;j++) {
			i=equiv[j];
			ped_recode1[i]=group[id_array[i].component-1]++;
		}
		for(i=0;i<ped_size;i++) {
			j=ped_recode1[i];
			if(j) equiv[j-1]=i;
		}
	}
	nbits=group;
	for(i=0;i<n_comp;i++) nbits[i]=0;
	for(j=0;j<pruned_ped_size;j++) {
		i=equiv[j];
		ids=id_array[i].sire;
		if(ids && !ped_recode1[ids-1]) ids=0;
		idd=id_array[i].dam;
		if(idd && !ped_recode1[idd-1]) idd=0;
		k=id_array[i].component-1;
		if(ids) nbits[k]+=2;
		else nbits[k]--;
	}
	if(verbose) {
		if(comp_sflag) {
			j=comp_size[n_comp-1];
			k=n_comp+comp_size[n_comp-1]-1;
		} else {
			j=0;
			k=n_comp;
		}
		if(verbose&1) (void)printf("No. components = %d\n",k);
		if((verbose&2)) {
			if(LogFile && (LogFile=add_file_dir(LogFile))) {
				flog=fopen(LogFile,"a");
				free(LogFile);
			} else flog=0;
			if(flog)	{
				(void)fputs("\n****************** Counting components ******************\n\n",flog);
				(void)fprintf(flog,"     No. components = %d\n",k);
				if(j) (void)fprintf(flog,"     No. components with single individuals (lumped) = %d\n",j);
				(void)fputc('\n',flog);
				for(i=0;i<n_comp-comp_sflag;i++) {
					fprintf(flog,"     Component %d, size = %d, nbits = %d\n",i+1,comp_size[i],nbits[i]<0?0:nbits[i]);
				}
				if(comp_sflag) fprintf(flog,"     Component %d, size = %d, singletons\n",n_comp,comp_size[n_comp-1]);
				fputc('\n',flog);
				(void)fclose(flog);
			}
		}
	}
	free(code);
	free(group);
}

void setup_pedigree(int num_nodes,struct recode_table_tag *recode_table,char *LogFile)
{
	int i,j,k,k1,k2,ncol,col,id_col,realcol,id,fam=-1,n_dummy,fam_col,s;
	struct InFile *infile;
	struct DataBlock *db;
	struct var_element *elem;
	FILE *flog;
	
	if(verbose&1) (void)printf("Initializing and recoding pedigree\n");
	infile=Infiles;
	k=k1=0;
	while(infile) {
		if(infile->n_records) {
			db=infile->data;
			ncol=infile->ncol;
			while(db) {
				for(j=0;j<db->record_ptr;j++)	{
					for(col=realcol=0;realcol<infile->nvar;realcol++) {
						elem=infile->element[realcol];
						if(elem)	{
							if(elem->type&(ST_FAMILY|ST_ID|ST_SIRE|ST_DAM)) if(!check_missing(col,ncol,j,db))	{
								i=db->records[j*ncol+col].node->index;
								if(!recode_table[i].index) {
									if(elem->type&ST_FAMILY) recode_table[i].index= --k1;
									else recode_table[i].index= ++k;
									recode_table[i].node=db->records[j*ncol+col].node;
								}
								if((elem->type&ST_FAMILY) || !family_id) db->records[j*ncol+col].value=recode_table[i].index;
							}
							col++;
						}
					}
				}
				db=db->next;
			}
		}
		infile=infile->next;
	}
	if(!k) ABT_FUNC("Internal error: No. id levels = 0\n");
	/* Count pedigree size and recode again if ids are repeated within families */
	if(family_id) {
		if(!k1) ABT_FUNC("Internal error: No. family levels = 0\n");
		k1=-k1;
		n_orig_families=k1;
		if(!(family_recode=malloc(sizeof(struct label_node *)*k1))) ABT_FUNC(MMsg);
		if(!(rec_tab1=calloc((size_t)(k1*k),sizeof(int)))) ABT_FUNC(MMsg);
		k=0;
		infile=Infiles;
		while(infile) {
			if(infile->n_records) {
				db=infile->data;
				ncol=infile->ncol;
				fam_col=infile->family_col;
				while(db) {
					for(j=0;j<db->record_ptr;j++)	{
						fam=-(int)db->records[j*ncol+fam_col].value-1;
						if(fam<0 || fam>=k1) ABT_FUNC("Internal error - bad family code\n");
						for(col=realcol=0;realcol<infile->nvar;realcol++) {
							elem=infile->element[realcol];
							if(elem)	{
								if(elem->type&(ST_ID|ST_SIRE|ST_DAM)) if(!check_missing(col,ncol,j,db))	{
									i=db->records[j*ncol+col].node->index;
									id=recode_table[i].index;
									s=(id-1)*k1+fam;
									if(!rec_tab1[s]) rec_tab1[s]=++k;
									db->records[j*ncol+col].value=rec_tab1[s];
								}
								col++;
							}
						}
					}
					db=db->next;
				}
			}
			infile=infile->next;
		}
	} 
	ped_size=k;
	if(!(ped_recode=malloc(sizeof(struct label_node *)*ped_size))) ABT_FUNC(MMsg);
	if(!(rec_tab=malloc(sizeof(int *)*num_nodes))) ABT_FUNC(MMsg);
	for(i=0;i<num_nodes;i++) {
		j=recode_table[i].index;
		if(j>0) {
			if(family_id) {
				s=(j-1)*k1;
				for(k2=0;k2<k1;k2++,s++) if(rec_tab1[s]) {
					ped_recode[rec_tab1[s]-1]=recode_table[i].node;
				}
			} else ped_recode[j-1]=recode_table[i].node;
		} else if(j<0) family_recode[-j-1]=recode_table[i].node;
		rec_tab[i]=j;
	}
	if(!(id_array=calloc((size_t)ped_size,sizeof(struct Id_Record)))) ABT_FUNC(MMsg);
	if(!(id_array[0].mkflag=malloc(sizeof(int)*(n_markers+1)*ped_size))) ABT_FUNC(MMsg);
	for(i=1;i<ped_size;i++) id_array[i].mkflag=id_array[i-1].mkflag+n_markers+1;
	/* Put sire/dam information into id_array */
	infile=Infiles;
	while(infile) {
		if(infile->n_records) {
			/* Find (one of) the column(s) the id is in, and check if file contains other pedigree data */
			for(j=i=0;i<infile->nvar;i++)	{
				elem=infile->element[i];
				if(elem && (elem->type&(ST_SIRE|ST_DAM))) {
					j=1;
					break;
				}
			}
			if(j) {/* Any sire/dam columns in data? */
				db=infile->data;
				ncol=infile->ncol;
				id_col=infile->id_col;
				fam_col=infile->family_col;
				while(db) {
					for(j=0;j<db->record_ptr;j++)	{
						id=(int)db->records[j*ncol+id_col].value;
						if(!id) ABT_FUNC("Internal error - ID is listed as missing\n");
						if(family_id) {
							fam=-(int)db->records[j*ncol+fam_col].value;
							id_array[id-1].fam_code=fam;
						}
						for(col=realcol=0;realcol<infile->nvar;realcol++) {
							elem=infile->element[realcol];
							if(elem)	{
								if(col!=id_col && col!=fam_col && (elem->type&(ST_ID|ST_SIRE|ST_DAM|ST_FAMILY))) {
									i=(int)db->records[j*ncol+col].value;
									if(i) switch(elem->type&(ST_ID|ST_SIRE|ST_DAM|ST_FAMILY))	{
									 case ST_FAMILY:
										if(fam!=-i) {
											print_orig_id(stderr,id,0);
											print_scan_err(" changes family code within 1 record\n");
										}
										break;
									 case ST_SIRE:
										if(id_array[id-1].sire) {
											if(id_array[id-1].sire!=i)	{
												print_orig_id(stderr,id,0);
												print_scan_err(" listed with multiple fathers\n");
											}
										} else id_array[id-1].sire=i;
										break;
									 case ST_DAM:
										if(id_array[id-1].dam) {
											if(id_array[id-1].dam!=i) {
												print_orig_id(stderr,id,0);
												print_scan_err(" listed with multiple mothers\n");
											}
										} else id_array[id-1].dam=i;
										break;
									 case ST_ID:
										if(i!=id) {
											print_orig_id(stderr,id,0);
											print_scan_err(" changes id within 1 record\n");
										}
									 default:
										ABT_FUNC("Internal error - weird variable type\n");
									}
								}
								col++;
							}
						}
					}
					db=db->next;
				}
			}
		}
		infile=infile->next;
	}
	for(i=0;i<ped_size;i++) {
		if(id_array[i].sire && !id_array[i].dam) id_array[i].dam= ++k;
		else if(!id_array[i].sire && id_array[i].dam) id_array[i].sire= ++k;
	}
	n_dummy=k-ped_size;
	if(n_dummy) {
		if(!(ped_recode=realloc(ped_recode,sizeof(struct label_node *)*k))) ABT_FUNC(MMsg);
		if(!(id_array=realloc(id_array,sizeof(struct Id_Record)*k))) ABT_FUNC(MMsg);
		for(;ped_size<k;ped_size++) {
			ped_recode[ped_size]=0;
			id_array[ped_size].sire=id_array[ped_size].dam=0;
			id_array[ped_size].haplo[0]=id_array[ped_size].haplo[1]=0;
			id_array[ped_size].flag=0;
			id_array[ped_size].sex=0;
		}
		if(!(id_array[0].mkflag=realloc(id_array[0].mkflag,sizeof(int)*(n_markers+1)*k))) ABT_FUNC(MMsg);
		for(i=1;i<ped_size;i++) id_array[i].mkflag=id_array[i-1].mkflag+n_markers+1;
	}
	/* Now recode ids (again!) so that parents are coded before children */
	/* This is only necessary for calculating the inverse NRM matrix, but */
     /* it can be useful at other places as well */
	if(!(ped_recode1=calloc((size_t)ped_size,sizeof(int)))) ABT_FUNC(MMsg);
	for(i=id=1;i<=ped_size;i++) if(!ped_recode1[i-1]) do_recode(i,&id);
	/* Check a few things with the pedigree */
	for(id=0;id<ped_size;id++) {
		i=id_array[id].sire;
		j=id_array[id].dam;
		if(i==(id+1) || j==(id+1)) {
			if(scan_error_n<max_scan_errors)	{
				(void)fputs("Individual ",stderr);
				print_orig_id(stderr,id+1,0);
				(void)fputs(" selfs!\n",stderr);
			}
			scan_error_n++;
		}
		if(i==j && i) {
			(void)fputs("Individual ",stderr);
			print_orig_id(stderr,id+1,0);
			(void)fputs(" has identical father and mother\n",stderr);
		}
	}
	if(verbose&2) {
		if(LogFile && (LogFile=add_file_dir(LogFile))) {
			flog=fopen(LogFile,"a");
			free(LogFile);
		} else flog=0;
		if(flog) {
			(void)fprintf(flog,"\n******************* Recoding pedigree *******************\n\n");
			(void)fprintf(flog,"     No. individuals found = %d\n",ped_size-n_dummy);
			(void)fprintf(flog,"     No. 'dummy' parents created = %d\n",n_dummy);
			(void)fprintf(flog,"     Pedigree size = %d\n",ped_size);
			if(flog!=stdout) (void)fclose(flog);
		}
	}
}

void handle_groups(char *LogFile)
{
	int i,j,k,k1,sire,dam,er=0,*ct,*ct1;
	FILE *flog;
	
 	for(i=0;i<ped_size;i++) id_array[i].group=0;
	if(!group_elem) return;
	if(verbose&1) (void)printf("Handling genetic groups\n");
	for(k=0;k<n_id_records;k++) if(id_elements[k]==group_elem) break;
	if(k==n_id_records) ABT_FUNC("Internal error - can't find group records\n");
	if(!group_elem->n_levels) ABT_FUNC("No data for genetic group found\n");
	if(!(ct=malloc(sizeof(int)*2*group_elem->n_levels))) ABT_FUNC(MMsg);
	ct1=ct+group_elem->n_levels;
	for(i=0;i<group_elem->n_levels;i++) ct[i]=0;
 	for(i=0;i<ped_size;i++)	{
		sire=id_array[i].sire;
		dam=id_array[i].dam;
		if(sire || dam) {
			if(!(sire && dam)) ABT_FUNC("Internal error - one parent unknown\n");
			if(!(ped_recode[sire-1] && ped_recode[dam-1])) {
				(void)fputs("Individual ",stderr);
				print_orig_id(stderr,i+1,0);
				(void)fputs(" has only one parent - illegal when using genetic groups\n",stderr);
				er=1;
			}
		}
	}
	if(er) ABT_FUNC("Errors - aborting\n");
	/* Flag which individuals have trait and marker data */
 	for(i=0;i<ped_size;i++)	{
		if(id_array[i].data && id_array[i].data[k].flag)
		  id_array[i].group=(int)id_array[i].data[k].data.value;
		sire=id_array[i].sire;
		dam=id_array[i].dam;
		if(sire || dam) id_array[i].group=0;
		else {
			if(!id_array[i].group) {
				if(ped_recode[i]) {
					(void)fputs("Founder individual ",stderr);
					print_orig_id(stderr,i+1,0);
					(void)fprintf(stderr," has no genetic group specified\n");
					er=1;
				} else ABT_FUNC("Internal error - dummy founder with no genetic information\n");
			} else ct[id_array[i].group-1]++;
		}
	}
	if(er) ABT_FUNC("Errors - aborting\n");
	for(j=0;j<n_factors;j++) if(group_elem==var_factors[j]) break;
	if(!(group_recode=malloc(sizeof(struct label_node *)*group_elem->n_levels))) ABT_FUNC(MMsg);
	for(k1=k=i=0;i<group_elem->n_levels;i++) if(ct[i]) {
		group_recode[k]=factor_recode[j][i];	
		ct[k++]=ct[i];
		ct1[i]=k;
		k1+=ct[i];
	}
	n_genetic_groups=k;
	for(i=0;i<ped_size;i++) if(id_array[i].group) id_array[i].group=ct1[id_array[i].group-1];
	if(verbose&1) (void)printf("No. genetic groups = %d\n",k);
	if((verbose&2) && LogFile) flog=fopen(LogFile,"a");
	if((verbose&2) && LogFile && (LogFile=add_file_dir(LogFile))) {
		flog=fopen(LogFile,"a");
		free(LogFile);
	} else flog=0;
	if(flog)	{
		(void)fprintf(flog,"\n********************* Genetic groups ********************\n\n");
		(void)fputs("         No. genetics groups       No. Records\n",flog);
		(void)fputs("         -------------------       -----------\n",flog);
		(void)fprintf(flog,"             %5d                  %5d\n\n",k,k1);
		if(k>1) for(i=0;i<k;i++) {
			if(factor_recode[j][i]->type==STRING) {
				k1=36-fprintf(flog,"                 %s",factor_recode[j][i]->data.string);
				while((k1--)>0) (void)fputc(' ',flog);
			} else (void)fprintf(flog,"                 %-12ld",factor_recode[j][i]->data.value);
			(void)fprintf(flog,"%5d\n",ct[i]);
		}
		(void)fclose(flog);
	}
	free(ct);
}

void prune_pedigree(char *LogFile)
{
	int i,i1,j,k,tj,*n_kids=0,sire,dam,*perm;
	FILE *flog;
	
	if((verbose&1) && syst_var[PRUNE_OPTION]) (void)printf("Pruning pedigree\n");
	/* Flag which individuals have trait and marker data */
 	for(i=0;i<ped_size;i++)	{
		id_array[i].flag=0;
		if(id_array[i].data) for(k=0;k<n_id_records;k++) {
			if((id_elements[k]->type&ST_TRAIT) && id_array[i].data[k].flag) {
				id_array[i].flag|=1;
				break;
			}
		}
		if(!id_array[i].flag && id_array[i].data1) for(j=0;j<id_array[i].nrec;j++) {
			if(id_array[i].data1[j]) for(k=0;k<n_nonid_records;k++) {
				if((nonid_elements[k]->type&ST_TRAIT) && id_array[i].data1[j][k].flag) {
					id_array[i].flag|=1;
					break;
				}
			}
		}
		if(id_array[i].group) id_array[i].flag|=1;
		if(id_array[i].haplo[0]) {
			for(k=j=0;j<n_markers;j++) {
				if(id_array[i].haplo[0][j] && id_array[i].haplo[1][j]) {
					id_array[i].flag|=2;
					k=1;
				} else id_array[i].haplo[0][j]=id_array[i].haplo[1][j]=0;
			}
			if(!k) {
				free(id_array[i].haplo[0]);
				id_array[i].haplo[0]=id_array[i].haplo[1]=0;
			}
		}
	}
	if(!syst_var[PRUNE_OPTION]) {
		pruned_ped_size=ped_size;
		return;
	}
	if(!(n_kids=calloc((size_t)ped_size*2,sizeof(int)))) ABT_FUNC(MMsg);
	perm=n_kids+ped_size;
	for(i=0;i<ped_size;i++) perm[ped_recode1[i]-1]=i;
	tj=0;
	do {
		for(i1=0;i1<ped_size;i1++) {
			i=perm[ped_size-i1-1];
			if(!(id_array[i].flag&8)) {
				if(id_array[i].flag&19) {
					sire=id_array[i].sire;
					dam=id_array[i].dam;
					if(sire)	{
						n_kids[sire-1]++;
						id_array[sire-1].flag|=16;
					}
					if(dam) {
						n_kids[dam-1]++;
						id_array[dam-1].flag|=16;
					}
				}
			}
		}
		for(i=j=0;i<ped_size;i++) {
			if(!(id_array[i].flag&8)) {
				if(!(id_array[i].flag&3)) {
					if(!n_kids[i])	{
						id_array[i].flag|=8;
						j++;
					} else if(n_kids[i]==1 && !id_array[i].sire) {
						id_array[i].flag|=8;
						j++;
					}
				}
			}
		}
		if(j)	{
			for(i=0;i<ped_size;i++)	{
				if(!(id_array[i].flag&8)) {
					sire=id_array[i].sire;
					dam=id_array[i].dam;
					if(sire && dam) {
						if(id_array[sire-1].flag&8) {
							if(id_array[dam-1].flag&8) id_array[i].sire=id_array[i].dam=0;
							else {
								id_array[sire-1].flag&=3;
								j--;
							}
						} else if(id_array[dam-1].flag&8) {
							id_array[dam-1].flag&=3;
							j--;
						}
				   } 
				}
				n_kids[i]=0;
				id_array[i].flag&=15;
			}
		}
		tj+=j;
	}
	while(j);
	if(tj) {
		if(verbose&1) (void)printf("No. pruned = %d\n",tj);
		if((verbose&2) && LogFile && (LogFile=add_file_dir(LogFile))) {
			flog=fopen(LogFile,"a");
			free(LogFile);
		} else flog=0;
		if(flog) {
			(void)fprintf(flog,"\n******************** Pruning pedigree *******************\n\n");
			(void)fprintf(flog,"     No. individuals pruned = %d, IDs follow:\n",tj);
		}
		for(i1=0,k=i=0;i<ped_size;i++) {
			j=perm[i];
			if(id_array[j].flag&8) {
				ped_recode1[j]=0;
				if(ped_recode[j] && flog)	{
					if(!((i1++)%10)) (void)fputs("\n     ",flog);
					else (void)fputc(',',flog);
					print_orig_id(flog,j+1,0);
				}
				k++;
			} else ped_recode1[j]-=k;
		}
		pruned_ped_size=ped_size-k;
		if(flog)	{
			(void)fputc('\n',flog);
			(void)fprintf(flog,"\n     New pedigree size = %d\n",pruned_ped_size);
		}
		if(verbose&1) (void)printf("Pruned pedigree size = %d\n",pruned_ped_size);
		if(flog && flog!=stdout) (void)fclose(flog);
	} else pruned_ped_size=ped_size;
	free(n_kids);
}

void check_sex(void)
{
	int i,j,k,k1,k2,ids,idd,er=0;
	const struct label_data *node;
	struct sex_def *sd;
	struct var_element *sex_elem;
	struct express **sex_exp;
	
	for(i=0;i<ped_size;i++)	id_array[i].sex=0;
	if(sex_def) {
		sd=sex_def;
		if(verbose&1) (void)printf("Checking sex\n");
		while(sd) {
			sex_elem=sd->sex_elem;
			sex_exp=sd->sex_exp;
			sd=sd->next;
			for(k=0;k<n_id_records;k++) if(id_elements[k]==sex_elem) break;
			if(k==n_id_records) ABT_FUNC("Couldn't find sex record\n");
			j=sex_elem->index;
			if(sex_elem->type&ST_INTTYPE)	{
				if(sex_exp[0]->type!=ST_INTEGER) ABT_FUNC("Type mismatch with sex command\n");
			} else if(sex_exp[0]->type!=ST_STRING) ABT_FUNC("Type mismatch with sex command\n");
			for(i=0;i<ped_size;i++)	{
				if(id_array[i].data && id_array[i].data[k].flag) {
					k1=(int)id_array[i].data[k].data.value;
					node=factor_recode[j][k1-1];
					if(node->type==INTEGER) {
						for(k2=0;k2<2;k2++) if(node->data.value==sex_exp[k2]->arg.value) break;
						if(k2==2 && id_array[i].sex<=0) id_array[i].sex=-1;
						else if(id_array[i].sex>0 && id_array[i].sex!=k2+5) id_array[i].sex=-2;
						else id_array[i].sex=k2+5;
					} else if(node->type==STRING) {
						for(k2=0;k2<2;k2++) if(!mystrcmp(node->data.string,sex_exp[k2]->arg.string)) break;
						if(k2==2 && id_array[i].sex<=0) id_array[i].sex=-1;
						else if(id_array[i].sex>0 && id_array[i].sex!=k2+5) id_array[i].sex=-2;
						else id_array[i].sex=k2+5;
					}
				}
			}
		}
	}
	for(i=0;i<ped_size;i++)	if(id_array[i].sex<0) {
		print_orig_id(stderr,i+1,0);
		if(id_array[i].sex==-1) 
		  (void)fputs(" has invalid sex information\n",stderr);
		else
		  (void)fputs(" has conflicting sex information\n",stderr);
		er=1;
	}
	if(er) ABT_FUNC("Errors - aborting\n");
 	for(i=0;i<ped_size;i++) {
		ids=id_array[i].sire;
		idd=id_array[i].dam;
		if(ids && id_array[ids-1].sex>=0) {
			if(id_array[ids-1].sex && (id_array[ids-1].sex&3)!=1) {
				print_orig_id(stderr,ids,0);
				if(id_array[ids-1].sex&4) (void)fprintf(stderr," declared as female yet used as a father\n");
				else (void)fprintf(stderr," used as both father and mother\n");
				er=1;
				id_array[ids-1].sex= -1;
			} else id_array[ids-1].sex=1;
		}
		if(idd && id_array[idd-1].sex>=0) {
			if(id_array[idd-1].sex && (id_array[idd-1].sex&3)!=2)	{
				print_orig_id(stderr,idd,0);
				if(id_array[idd-1].sex&4) (void)fprintf(stderr," declared as male yet used as a mother\n");
				else (void)fprintf(stderr," used as both father and mother\n");
				er=1;
				id_array[idd-1].sex= -2;
			} else id_array[idd-1].sex=2;
		}
	}
	if(er) ABT_FUNC("Pedigree errors - aborting\n");
	for(i=0;i<ped_size;i++) id_array[i].sex&=3;
}

void check_sex2(void)
{
	int i,sex_reqd=0,er=0;
	struct Link *linkp;
	
	linkp=links;
	while(linkp && !sex_reqd) {
		if(linkp->type!=LINK_AUTO) sex_reqd=1;
		linkp=linkp->next;
	}
	if(sex_reqd) for(i=0;i<ped_size;i++) if(ped_recode1[i] && !id_array[i].sex) {
		print_orig_id(stderr,i+1,0);
		fputs(" has no sex information (required for sex linked loci)\n",stderr);
		er=1;
	}
	if(er) ABT_FUNC("Pedigree errors - aborting\n");
}

static void inbr_rep(FILE *fptr,int i,int ids,int idd,char *inbr[])
{
	static int fg;
	
	if(!fg) {
		fg=1;
		(void)fputs("\n**************** Close inbreeding report *****************\n\n",fptr);
	}
	(void)fputs("     ",fptr);
	if(family_id) print_orig_family(fptr,ids,0);
	fputc(' ',fptr);
	print_orig_id1(fptr,ids,1);
	print_orig_id1(fptr,idd,0);
	(void)fprintf(fptr," - %s marriage\n",inbr[i]);
}

static int check_siblings(int ids,int idd)
{
	int j=0;

	if(ids&&idd) {
		if(id_array[ids-1].sire && id_array[ids-1].sire==id_array[idd-1].sire) j++;
		if(id_array[ids-1].dam && id_array[ids-1].dam==id_array[idd-1].dam) j++;
	}
	return j;
}

static int check_ancestor(int i,int j,int *k,FILE *fptr)
{
	if(i && j) {
		if(id_array[j-1].sire==i || id_array[j-1].dam==i) {
			if(fptr)	{
				print_orig_id(fptr,i,0);
				(void)fputc(',',fptr);
				print_orig_id(fptr,j,0);
			}
			return 1;
		}
		(*k)++;
		if(check_ancestor(i,id_array[j-1].sire,k,fptr))	{
			if(fptr)	{
				(void)fputc(',',fptr);
				print_orig_id(fptr,j,0);
			}
			return 1;
		}
		if(check_ancestor(i,id_array[j-1].dam,k,fptr)) {
			if(fptr)	{
				(void)fputc(',',fptr);
				print_orig_id(fptr,j,0);
			}	
			return 1;
		}
		(*k)--;
	}
	return 0;
}

/* Produce report on various types of close inbreeding which could indicate pedigree errors */
void check_inbreeding(char *LogFile)
{
	int i,j,ids,idd,cc[7],gp1,gp2,gp3,gp4,k;
	static char *inbr[]={"full sib","half sib","offspring parent",
		"grandchild grandparent","niece/nephew aunt/uncle","first cousin","descendent ancestor"};
	FILE *flog;

	for(i=0;i<7;i++) cc[i]=0;
	if(LogFile && (LogFile=add_file_dir(LogFile))) {
		flog=fopen(LogFile,"a");
		free(LogFile);
	} else flog=0;
	for(i=0;i<n_families;i++) {
		ids=family[i].sire;
		idd=family[i].dam;
		if(ids && idd)	{
			j=check_siblings(ids,idd);
			if(j) {
				cc[2-j]++;
				if(flog) inbr_rep(flog,2-j,ids,idd,inbr);
			} else {
				k=1;
				if(ped_recode1[ids-1]<ped_recode1[idd-1]) j=check_ancestor(ids,idd,&k,0);
				else j=check_ancestor(idd,ids,&k,0);
				if(j) {
					if(k>2) k=5;
					cc[1+k]++;
					if(flog)	{
						inbr_rep(flog,1+k,ids,idd,inbr);
						if(k==5)	{
							(void)fputs("       -> ",flog);
							if(ped_recode1[ids-1]<ped_recode1[idd-1]) (void)check_ancestor(ids,idd,&k,flog);
							else (void)check_ancestor(idd,ids,&k,flog);
							(void)fputc('\n',flog);
						}
					}
				}
			}
			if(!j) {
				gp1=id_array[ids-1].sire;
				gp2=id_array[ids-1].dam;
				gp3=id_array[idd-1].sire;
				gp4=id_array[idd-1].dam;
				j=check_siblings(gp1,gp3);
				if(j!=2) j=check_siblings(gp1,gp4);
				if(j!=2) j=check_siblings(gp2,gp3);
				if(j!=2) j=check_siblings(gp2,gp4);
				if(j==2)	{
					cc[5]++;
					if(flog) inbr_rep(flog,5,ids,idd,inbr);
				} else {
					j=check_siblings(gp1,idd);
					if(j!=2) j=check_siblings(gp2,idd);
					if(j!=2) j=check_siblings(gp3,ids);
					if(j!=2) j=check_siblings(gp3,ids);
					if(j==2)	{
						cc[4]++;
						if(flog) inbr_rep(flog,4,ids,idd,inbr);
					}
				}
			}
		}
	}
	if(flog) (void)fclose(flog);
	for(i=0;i<7;i++) if(cc[i]) break;
	if(i<7) (void)printf("Close inbreeding found - more details in logfile\n");
	for(i=0;i<7;i++) if(cc[i]) (void)printf("%d %s marriage%c\n",cc[i],inbr[i],cc[i]==1?' ':'s');
}

void count_loops(char *LogFile)
{
	int i,j,k1,k2,k3,k4,pivot,kid,fam,ids,idd,*famlist,nf,*compflag;
	FILE *flog;
	
	if(LogFile && (LogFile=add_file_dir(LogFile))) {
		flog=fopen(LogFile,"a");
		free(LogFile);
	} else flog=0;
	if(flog) (void)fputs("\n********************* Counting loops *********************\n",flog);
	for(i=0;i<ped_size;i++) id_array[i].order=0;
	if(!(compflag=malloc(sizeof(int)*(n_comp+n_families)))) ABT_FUNC(MMsg);
	famlist=compflag+n_comp;
	for(i=0;i<n_comp;i++) compflag[i]=0;
	for(fam=0;fam<n_families;fam++) {
		famlist[fam]=fam;
		ids=family[fam].sire;
		idd=family[fam].dam;
		if(ids) {
			id_array[ids-1].order++;
			id_array[idd-1].order++;
			for(j=0;j<family[fam].nkids;j++) {
				kid=family[fam].kids[j];
				id_array[kid].order++;
			}
		} else {
			kid=family[fam].kids[0];
			id_array[kid].order++;
		}
	}
	nf=n_families;
	k4=0;
	while(nf) {
		for(i=0;i<nf;i++) {
			fam=famlist[i];
			ids=family[fam].sire;
			idd=family[fam].dam;
			pivot=0;
			if(ids) {
				if(id_array[ids-1].order>1) pivot=ids;
				if(id_array[idd-1].order>1) {
					if(pivot) pivot= -1;
					else pivot=idd;
				}
			}
			if(pivot<0) continue;
			for(j=0;j<family[fam].nkids;j++) {
				kid=family[fam].kids[j];
				if(id_array[kid].order>1) {
					if(pivot) break;
					pivot=kid+1;
				}
			}
			if(j<family[fam].nkids) continue;
			if(pivot>0) id_array[pivot-1].order--;
			famlist[i]=famlist[--nf];
			i= -1;
		}
		if(nf) {
			for(k3=INT_MAX,k2=i=0;i<nf;i++) {
				fam=famlist[i];
				ids=family[fam].sire;
				idd=family[fam].dam;
				k1=0;
				if(ids) {
					if(id_array[ids-1].order>1) k1++;
					if(id_array[idd-1].order>1) k1++;
				}
				for(j=0;j<family[fam].nkids;j++) {
					kid=family[fam].kids[j];
					if(id_array[kid].order>1) k1++;
				}
				if(k1<k3) {
					k3=k1;
					k2=i;
				}
			}
			fam=famlist[k2];
			ids=family[fam].sire;
			idd=family[fam].dam;
			if(ids) {
				if(id_array[ids-1].order>1) {
					k4+=id_array[ids-1].order-1;
					compflag[id_array[ids-1].component-1]=ids;
					id_array[ids-1].order=1;
				}
				if(id_array[idd-1].order>1) {
					k4+=id_array[idd-1].order-1;
					compflag[id_array[idd-1].component-1]=idd;
					id_array[idd-1].order=1;
				}
			}
			for(j=0;j<family[fam].nkids;j++) {
				kid=family[fam].kids[j];
				if(id_array[kid].order>1) {
					k4+=id_array[kid].order-1;
					id_array[kid].order=1;
					compflag[id_array[kid].component-1]=kid+1;
				}
			}
			k4--;
			famlist[k2]=famlist[--nf];
		}
	}
	if(k4) (void)printf("No. loops = %d\n",k4);
	else (void)fputs("Pedigree has no loops\n",stdout);
	if(flog) {
		if(k4) {
			for(i=j=0;i<n_comp;i++) if(compflag[i]) j++;
			k1=fprintf(flog,"\n     No. loops = %d.  Loops found in componen%s:",k4,j>1?"ts":"t");
			for(i=0;i<n_comp;i++) if(compflag[i]) {
				if(!k1 || family_id) k1=fputs("\n    ",flog);
				k1+=fprintf(flog," %d",i+1);
				if(family_id) {
					fputs(" (",flog);
					print_orig_family(flog,compflag[i],0);
					fputc(')',flog);
				}
				if(k1>=60) k1=0;
			}
			if(k1) (void)fputc('\n',flog);
		} else (void)fputs("\n     Pedigree has no loops\n",flog);
		(void)fclose(flog);
	}
	free(compflag);
}

static void do_yprint(int i,const int *list,const int ny)
{
	int j,k,ch,locus;
	
	print_orig_id(stdout,i+1,1);
	if(id_array[i].haplo[0]) {
		for(j=0;j<ny;j++)	{
			if(j) (void)fputc('_',stdout);
			locus=list[j];
			ch=id_array[i].haplo[0][j];
			if(ch) {
				if(factor_recode[n_factors+locus][ch-1]->type==STRING)
				  (void)fputs(factor_recode[n_factors+locus][ch-1]->data.string,stdout);
				else (void)fprintf(stdout,"%ld",factor_recode[n_factors+locus][ch-1]->data.value);
			} else (void)fputc('*',stdout);
		}
	}
	(void)fputc('\n',stdout);
	for(j=0;j<id_array[i].nkids;j++) {
		k=id_array[i].kids[j];
		if(id_array[k].sex==1) do_yprint(k,list,ny);
	}
}

void check_ymark(void)
{
	int i,j,locus,ny=0,*list;
	struct Link *pl;
	
	for(locus=0;locus<n_markers;locus++) {
		pl=links;
		i=markers[locus].link-1;
		while(i && pl) {
			pl=pl->next;
			i--;
		}
		if(!pl) ABT_FUNC("Invalid linkage group\n");
		if(pl->type==LINK_Y) ny++;
	}
	if(ny) (void)printf("No. y markers = %d\n",ny);
	else return;
	if(!(list=malloc(sizeof(int)*ny))) ABT_FUNC(MMsg);
	ny=0;
	free(list);
	for(locus=0;locus<n_markers;locus++) {
		pl=links;
		i=markers[locus].link-1;
		while(i && pl)	{
			pl=pl->next;
			i--;
		}
		if(pl->type==LINK_Y) list[ny++]=locus;
	}
	for(i=0;i<ped_size;i++) {
		j=ped_recode1[i];
		if(!j || (id_array[i].sex!=1 || id_array[i].sire)) continue;
		(void)fputs("\n***********\n",stdout);
		do_yprint(i,list,ny);
	}
}

