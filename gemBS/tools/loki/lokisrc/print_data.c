/****************************************************************************
 *                                                                          *
 *     Loki - Programs for genetic analysis of complex traits using MCMC    *
 *                                                                          *
 *                      Simon Heath - MSKCC                                 *
 *                                                                          *
 *                          August 2000                                     *
 *                                                                          *
 * print_data.c:                                                            *
 *                                                                          *
 * Auxillary routines for writing out coded data                            *
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
#include <stdio.h>

#include "utils.h"
#include "loki.h"
#include "loki_peel.h"
#include "loki_utils.h"
#include "version.h"
#include "print_data.h"

void Print_Data(const char *filename,struct loki *loki)
{
	int i,k,k1,k2,rec,nrec,type;
	struct id_data *data;
	struct Variable *var;
	double x;
	FILE *fptr;
	char *s;
	struct Id_Record *id_array;
	struct Model *mod;
	
	if(!(fptr=fopen(filename,"w"))) {
		(void)fprintf(stderr,"Couldn't open output file '%s': ",filename);
		perror(0);
		return;
	}
	id_array=loki->pedigree->id_array;
	mod=loki->models->models;
	for(i=0;i<loki->pedigree->ped_size;i++) {
		nrec=id_array[i].n_rec;
		if(!nrec) nrec=1;
		for(rec=0;rec<nrec;rec++) {
			(void)fprintf(fptr,"%d ",id_array[i].comp+1);
			print_orig_triple(fptr,i+1);
			(void)fputc(' ',fptr);
			if(id_array[i].sex)	(void)fprintf(fptr,"%d ",id_array[i].sex);
			else (void)fputs("* ",fptr);
			for(k=0;k<mod->n_terms;k++) {
				type=mod->term[k].vars[0].type;
				if(type&(ST_TRAITLOCUS|ST_SEX)) continue;
				k1=mod->term[k].vars[0].var_index;
				if(type&ST_MARKER) {
					if(loki->markers->marker[k1].haplo[i]) (void)fprintf(fptr,"%d %d ",loki->markers->marker[k1].haplo[i]&65535,loki->markers->marker[k1].haplo[i]>>16);
					else (void)fputs("* * ",fptr);
				} else {
					data=0;
					if(type&ST_CONSTANT) {
						if(id_array[i].data) data=id_array[i].data+k1;
					} else if(id_array[i].data1) data=id_array[i].data1[rec]+k1;
					if(!data || !data->flag) (void)fputs("* ",fptr);
					else {
						if(type&ST_FACTOR) {
							k2=(int)data->data.value;
							(void)fprintf(fptr,"%d ",k2);
						} else {
							if(data->flag&ST_INTTYPE) x=(double)data->data.value;
							else x=data->data.rvalue;
							(void)fprintf(fptr,"%g ",x);
						}
					}
				}
			}
			if(id_array[i].res[0]) {
				type=mod->var.type;
				k1=mod->var.var_index;
				if(type&ST_CONSTANT) data=id_array[i].data+k1;
				else data=id_array[i].data1[rec]+k1;
				if(type&ST_FACTOR) {
					k2=(int)data->data.value;
					(void)fprintf(fptr,"%d\n",k2);
				} else {
					if(data->flag&ST_INTTYPE) x=(double)data->data.value;
					else x=data->data.rvalue;
					(void)fprintf(fptr,"%g\n",x);
				}
			} else (void)fputs("*\n",fptr);
		}
	}
	(void)fclose(fptr);
	i=strlen(filename)+5;
	if(!(s=malloc((size_t)i))) ABT_FUNC(MMsg);
	(void)snprintf(s,i,"%s_idx",filename);
	if(!(fptr=fopen(s,"w"))) {
		(void)fprintf(stderr,"Couldn't open output index file '%s': ",s);
		perror(0);
	} else {
		(void)(void)fputs("******************* Phenotype Data **********************\n\n",fptr);
		(void)fprintf(fptr,"     %s: %s",LOKI_NAME,ctime(&loki->sys.lktime.start_time));
		(void)fputs("\n\n     1  component\n     2  id\n     3  father\n     4  mother\n     5  sex\n",fptr);
		i=5;
		for(k=0;k<mod->n_terms;k++) {
			type=mod->term[k].vars[0].type;
			k1=mod->term[k].vars[0].var_index;
			if(type&(ST_TRAITLOCUS|ST_SEX)) continue;
			if(type&ST_MARKER) {
				(void)fprintf(fptr,"   %3d  ",++i);
				print_marker_name(fptr,k);
			} else {
				if(type&ST_CONSTANT) var=loki->data->id_variable+k1;
				else var=loki->data->nonid_variable+k1;
				(void)fprintf(fptr,"   %3d  %s",++i,var->name);
				if(var->index) (void)fprintf(fptr,"(%d)",var->index);
			}
			(void)fputc('\n',fptr);
		}
		type=mod->var.type;
		k1=mod->var.var_index;
		if(type&ST_CONSTANT) var=loki->data->id_variable+k1;
		else var=loki->data->nonid_variable+k1;
		(void)fprintf(fptr,"   %3d  %s",++i,var->name);
		if(var->index) (void)fprintf(fptr,"(%d)",var->index);
		(void)fputc('\n',fptr);
		(void)fclose(fptr);
	}
	free(s);
}

void Print_Genotypes(const struct output_gen *og,struct loki *loki)
{
	int i,j,k,k1,k2,l;
	FILE *fptr;
	char *s;
	struct Locus **locilist;
	struct Link *linkage;
	struct Marker *marker;
	
	if(!loki->markers->n_markers) return;
	if(!(fptr=fopen(og->file,"w"))) {
		(void)fprintf(stderr,"Couldn't open output file '%s': ",og->file);
		perror(0);
		return;
	}
	if(!(locilist=malloc(sizeof(void *)*loki->markers->n_markers))) ABT_FUNC(MMsg);
	marker=loki->markers->marker;
	linkage=loki->markers->linkage;
	for(j=0;j<loki->pedigree->ped_size;j++) {
		if(og->link_group) i=k=og->link_group-1;
		else {
			i=0;
			k=loki->markers->n_links-1;
		}
		for(;i<=k;i++) if(linkage[i].n_markers) {
			for(k1=0;k1<linkage[i].n_markers;k1++) {
				k2=linkage[i].mk_index[k1];
				if(marker[k2].haplo[j]) break;
			}
			if(k1<linkage[i].n_markers) break;
		}
		if(i>k) continue;
		if(loki->pedigree->family_id) {
			print_orig_family(fptr,j+1,0);
			(void)fputc(' ',fptr);
		}
		print_orig_id1(fptr,j+1);
		(void)fputc(' ',fptr);
		if(og->link_group) i=k=og->link_group-1;
		else {
			i=0;
			k=loki->markers->n_links-1;
		}
		for(;i<=k;i++) if(linkage[i].n_markers) {
			for(k1=0;k1<linkage[i].n_markers;k1++) locilist[k1]=&marker[linkage[i].mk_index[k1]].locus;
			gnu_qsort(locilist,(size_t)linkage[i].n_markers,(size_t)sizeof(void *),cmp_loci);
			for(k1=0;k1<linkage[i].n_markers;k1++) {
				k2=locilist[k1]->index;
				(void)fprintf(fptr,"%d %d ",marker[k2].haplo[j]>>16,marker[k2].haplo[j]&65535);
			}
		}
		(void)fputc('\n',fptr);
	}
	(void)fclose(fptr);
	i=strlen(og->file)+5;
	if(!(s=malloc((size_t)i))) ABT_FUNC(MMsg);
	(void)snprintf(s,i,"%s_idx",og->file);
	if(!(fptr=fopen(s,"w"))) {
		(void)fprintf(stderr,"Couldn't open output index file '%s': ",s);
		perror(0);
	} else {
		(void)(void)fputs("******************** Genotype Data **********************\n\n",fptr);
		(void)fprintf(fptr,"     %s: %s",LOKI_NAME,ctime(&loki->sys.lktime.start_time));
		if(og->link_group) (void)fprintf(fptr,"\n     Linkage group '%s'\n",linkage[og->link_group-1].name);
		(void)fputs("\n     1      id\n",fptr);
		if(og->link_group) i=k=og->link_group-1;
		else {
			i=0;
			k=loki->markers->n_links-1;
		}
		j=2;
		for(;i<=k;i++) if(linkage[i].n_markers) {
			for(k1=0;k1<linkage[i].n_markers;k1++) locilist[k1]=&marker[linkage[i].mk_index[k1]].locus;
			gnu_qsort(locilist,(size_t)linkage[i].n_markers,(size_t)sizeof(void *),cmp_loci);
			for(k1=0;k1<linkage[i].n_markers;k1++) {
				k2=locilist[k1]->index;
				(void)fprintf(fptr,"   %3d-%-3d  %s\n",j,j+1,marker[k2].name);
				j+=2;
			}
		}
		(void)fputs("\n   Allele frequencies\n\n",fptr);
		if(og->link_group) i=k=og->link_group-1;
		else i=0;
		for(;i<=k;i++) if(linkage[i].n_markers) {
			for(k1=0;k1<linkage[i].n_markers;k1++) locilist[k1]=&marker[linkage[i].mk_index[k1]].locus;
			gnu_qsort(locilist,(size_t)linkage[i].n_markers,(size_t)sizeof(void *),cmp_loci);
			for(k1=0;k1<linkage[i].n_markers;k1++) {
				(void)fputs("   ",fptr);
				k2=locilist[k1]->index;
				for(l=0;l<marker[k2].locus.n_alleles;l++)	{
					if(marker[k2].freq_set[0][l]) (void)fprintf(fptr,"%g ",marker[k2].locus.freq[0][l]);
					else (void)fputs("* ",fptr);
				}
				(void)fputc('\n',fptr);
			}
		}
		(void)fclose(fptr);
	}
	free(s);
	free(locilist);
}

