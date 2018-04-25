/****************************************************************************
 *                                                                          *
 *     Loki - Programs for genetic analysis of complex traits using MCMC    *
 *                                                                          *
 *             Simon Heath - University of Washington                       *
 *                                                                          *
 *                      June 1997                                           *
 *                                                                          *
 * write_report.c:                                                          *
 *                                                                          *
 * Write report to logfile                                                  *
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
#include <float.h>
#ifndef DBL_MAX
#define DBL_MAX MAXDOUBLE
#endif

#include "utils.h"
#include "scan.h"
#include "y.tab.h"

void WriteReport(char *LogFile)
{
	int i,j,k,k1,k2,k3,i1,rec,fg,*tmp,nlink,cens=0,comp,temp_size;
	FILE *flog=0;
	struct model *model;
	struct model_list *mlist,**model_list;
	struct Link *linkp,**llist;
	struct var_element *elem;
	struct scan_data *sd;
	double mean,min,max,y,het,hom;
	
	if(LogFile && (LogFile=add_file_dir(LogFile))) flog=fopen(LogFile,"a");
	if(LogFile) free(LogFile);
	if(!flog) return;
	(void)fputs("\n*********************** Final report ********************\n\n",flog);
	if(Models) (void)fprintf(flog,"     Model%s\n\n       ",Models->next?"s:":":");
	model=Models;
	k2=0;
	while(model) {
		if(model->trait->vtype&ST_ARRAY) (void)fprintf(flog,"%s(%d) = ",model->trait->name,model->index);
		else (void)fprintf(flog,"%s = ",model->trait->name);
		j=0;
		mlist=model->model_list;
		while(mlist) {
			j++;
			mlist=mlist->next;
		}
		if(!(model_list=malloc(sizeof(void *)*j))) ABT_FUNC(MMsg);
		k=j;
		mlist=model->model_list;
		while(mlist) {
			model_list[--j]=mlist;
			mlist=mlist->next;
		}
		for(j=0;j<k;j++) {
			mlist=model_list[j];
			for(i=0;i<mlist->nvar;i++) {
				elem=mlist->element[i];
				sd=elem->arg.var->data;
				if(sd->vtype&ST_ARRAY) {
					if(elem->type&ST_MARKER) k1=markers[elem->index].index;
					else k1=elem->oindex;
					(void)fprintf(flog,"%s(%d)",sd->name,k1);
				} else (void)fprintf(flog,"%s",sd->name);
				if(elem->type&(ST_RANDOM|ST_ID|ST_DAM|ST_SIRE)) {
					(void)fputc('\'',flog);
					k2=1;
				}
				if(i<mlist->nvar-1) (void)fputs("*",flog);
			}
			if(j<k-1) (void)fputs(" + ",flog);
		}
		free(model_list);
		model=model->next;
		if(model) (void)fputs("\n       ",flog);
	}
	if(k2) (void)fputs("\n\n       ( \' indicates a random effect )",flog);
	fg=0;
	for(i=0;i<n_id_records;i++) if(id_elements[i]->type&(ST_MODEL|ST_TRAIT)) {
		elem=id_elements[i];
		if(elem->type&ST_CENSORED) cens=1;
		if(elem->type&ST_FACTOR) continue;
		if(!fg) {
			(void)fputs("\n\n       Name           Mean        Min      Max     No. Records\n",flog);
			(void)fputs("       ----           ----        ---      ---     -----------\n",flog);
			fg=1;
		}
		mean=0.0;
		max= -DBL_MAX;
		min=DBL_MAX;
		for(k=j=0;j<ped_size;j++) if(ped_recode1[j] && id_array[j].data) {
			if(id_array[j].data[i].flag&1) {
				k++;
				if(elem->type&ST_INTTYPE) y=(double)id_array[j].data[i].data.value;
				else y=id_array[j].data[i].data.rvalue;
				mean+=y;
				if(y>max) max=y;
				if(y<min) min=y;
			}
		}
		sd=elem->arg.var->data;
		if(sd->vtype&ST_ARRAY) k1=20-fprintf(flog,"       %s(%d)",sd->name,elem->oindex);
		else k1=20-fprintf(flog,"       %s",sd->name);
		while((k1--)>0) (void)fputc(' ',flog);
		if(k) {
			mean/=(double)k;
 			if(elem->type&ST_INTTYPE) (void)fprintf(flog,"%8g %8d %8d       %6d\n",mean,(int)min,(int)max,k);
			else (void)fprintf(flog,"%8g %8g %8g       %6d\n",mean,min,max,k);
		} else (void)fprintf(flog,"%8g %8s %8s       %6d\n",mean," "," ",k);
	}
	for(i=0;i<n_nonid_records;i++) if(nonid_elements[i]->type&(ST_MODEL|ST_TRAIT)) {
		elem=nonid_elements[i];
		if(elem->type&ST_CENSORED) cens=1;
		if(elem->type&ST_FACTOR) continue;
		if(!fg) {
			(void)fputs("\n\n       Name           Mean        Min      Max     No. Records\n",flog);
			(void)fputs("       ----           ----        ---      ---     -----------\n",flog);
			fg=1;
		}
		mean=0.0;
		max= -DBL_MAX;
		min=DBL_MAX;
		for(k=j=0;j<ped_size;j++) if(ped_recode1[j] && id_array[j].data1)	{
			for(k1=0;k1<id_array[j].nrec;k1++) if(ped_recode1[j] && id_array[j].data1[k1]) {
				if(id_array[j].data1[k1][i].flag&1)	{
					k++;
					if(elem->type&ST_INTTYPE) y=(double)id_array[j].data1[k1][i].data.value;
					else y=id_array[j].data1[k1][i].data.rvalue;
					mean+=y;
					if(y>max) max=y;
					if(y<min) min=y;
				}
			}
		}
		sd=elem->arg.var->data;
		if(sd->vtype&ST_ARRAY) k1=20-fprintf(flog,"       %s(%d)",sd->name,elem->oindex);
		else k1=20-fprintf(flog,"       %s",sd->name);
		while((k1--)>0) (void)fputc(' ',flog);
		if(k) {
			mean/=(double)k;
			if(elem->type&ST_INTTYPE) (void)fprintf(flog,"%8g %8d %8d       %6d\n",mean,(int)min,(int)max,k);
			else (void)fprintf(flog,"%8g %8g %8g       %6d\n",mean,min,max,k);
		} else (void)fprintf(flog,"%8g %8s %8s       %6d\n",mean," "," ",k);
	}
	fg=0;
	tmp=0;
	temp_size=0;
	for(i=0;i<n_id_records;i++) if(id_elements[i]->type&(ST_MODEL|ST_TRAIT)) {
		elem=id_elements[i];
		if(!(elem->type&ST_FACTOR)) continue;
		if(!fg) {
			(void)fputs("\n\n       Name        No. levels      No. Records\n",flog);
			(void)fputs("       ----        ----------      -----------\n",flog);
			fg=1;
		}
		if(elem->n_levels>temp_size) {
			if(tmp) {
				if(!(tmp=realloc(tmp,sizeof(int)*elem->n_levels))) ABT_FUNC(MMsg);
			} else if(!(tmp=malloc(sizeof(int)*elem->n_levels))) ABT_FUNC(MMsg);
			temp_size=elem->n_levels;
		}
		for(j=0;j<elem->n_levels;j++) tmp[j]=0;
		for(k=j=0;j<ped_size;j++) if(ped_recode1[j] && id_array[j].data) {
			if(id_array[j].data[i].flag&1) {
				k++;
				tmp[id_array[j].data[i].data.value-1]++;
			}
		}
		for(k2=j=0;j<elem->n_levels;j++) if(tmp[j]) k2++;
		sd=elem->arg.var->data;
		if(sd->vtype&ST_ARRAY) k1=20-fprintf(flog,"       %s(%d)",sd->name,elem->oindex);
		else k1=20-fprintf(flog,"       %s",sd->name);
		while((k1--)>0) (void)fputc(' ',flog);
		(void)fprintf(flog,"%5d             %5d\n\n",k2,k);
		for(k=0;k<n_factors;k++) if(elem==var_factors[k]) break;
		for(j=0;j<elem->n_levels;j++) {
			if(tmp[j]) {
				if(factor_recode[k][j]->type==STRING) {
					k1=38-fprintf(flog,"                        %s",factor_recode[k][j]->data.string);
					while((k1--)>0) (void)fputc(' ',flog);
				} else (void)fprintf(flog,"                        %-14ld",factor_recode[k][j]->data.value);
				(void)fprintf(flog,"%5d\n",tmp[j]);
			}
		}
		if(k2) (void)fputc('\n',flog);
		if(k2<j) {
			for(k1=j=0;j<elem->n_levels;j++) {
				if(tmp[j]) {
					factor_recode[k][k1]=factor_recode[k][j];
					tmp[j]=k1;
					k1++;
				}
			}
			for(k=j=0;j<ped_size;j++) if(ped_recode1[j] && id_array[j].data) {
				if(id_array[j].data[i].flag&1) {
					k1=id_array[j].data[i].data.value-1;
					id_array[j].data[i].data.value=tmp[k1]+1;
				}
			}
			elem->n_levels=k2;
		}
	}
	for(i=0;i<n_nonid_records;i++) if(nonid_elements[i]->type&(ST_MODEL|ST_TRAIT)) {
		elem=nonid_elements[i];
		if(!(elem->type&ST_FACTOR)) continue;
		if(!fg) {
			(void)fputs("\n\n       Name        No. levels      No. Records\n",flog);
			(void)fputs("       ----        ----------      -----------\n",flog);
			fg=1;
		}
		if(elem->n_levels>temp_size) {
			if(tmp) {
				if(!(tmp=realloc(tmp,sizeof(int)*elem->n_levels))) ABT_FUNC(MMsg);
			} else if(!(tmp=malloc(sizeof(int)*elem->n_levels))) ABT_FUNC(MMsg);
			temp_size=elem->n_levels;
		}
		for(j=0;j<elem->n_levels;j++) tmp[j]=0;
		for(k=j=0;j<ped_size;j++) if(ped_recode1[j] && id_array[j].data1) {
			for(k1=0;k1<id_array[j].nrec;k1++) if(id_array[j].data1[k1]) {
				if(id_array[j].data1[k1][i].flag&1)	{
					k++;
					tmp[id_array[j].data1[k1][i].data.value-1]++;
				}
			}
		}
		for(k2=j=0;j<elem->n_levels;j++) if(tmp[j]) k2++;
		sd=elem->arg.var->data;
		if(sd->vtype&ST_ARRAY) k1=20-fprintf(flog,"       %s(%d)",sd->name,elem->oindex);
		else k1=20-fprintf(flog,"       %s",sd->name);
		while((k1--)>0) (void)fputc(' ',flog);
		(void)fprintf(flog,"%5d             %5d\n\n",k2,k);
		for(k=0;k<n_factors;k++) if(elem==var_factors[k]) break;
		for(j=0;j<elem->n_levels;j++) {
			if(factor_recode[k][j]->type==STRING) {
				k1=38-fprintf(flog,"                        %s",factor_recode[k][j]->data.string);
				while((k1--)>0) (void)fputc(' ',flog);
			} else (void)fprintf(flog,"                        %-14ld",factor_recode[k][j]->data.value);
			(void)fprintf(flog,"%5d\n",tmp[j]);
		}
		if(k2) (void)fputc('\n',flog);
		if(k2<j) {
			for(k1=j=0;j<elem->n_levels;j++) {
				if(tmp[j]) {
					factor_recode[k][k1]=factor_recode[k][j];
					tmp[j]=k1;
					k1++;
				}
			}
			for(k=j=0;j<ped_size;j++) if(ped_recode1[j] && id_array[j].data1) {
				for(k1=0;k1<id_array[j].nrec;k1++) if(id_array[j].data1[k1]) {
					if(id_array[j].data1[k1][i].flag&1) {
						k3=id_array[j].data1[k1][i].data.value-1;
						id_array[j].data1[k1][i].data.value=tmp[k3]+1;
					}
				}
			}
			elem->n_levels=k2;
		}
	}
	if(tmp) free(tmp);
	if(cens) {
		(void)fputs("\n\n     Censored traits:\n\n       Name           Records   Cens.   Non Cens.\n",flog);
		(void)fputs("       ----           -------   ----    ---------\n",flog);
		for(i=0;i<n_id_records;i++) if(id_elements[i]->type&ST_CENSORED) {
			elem=id_elements[i];
			if(!elem->type&(ST_MODEL|ST_TRAIT)) continue;
			for(k=k1=k2=j=0;j<ped_size;j++) if(ped_recode1[j] && id_array[j].data) {
				if(id_array[j].data[i].flag&1) k++;
				if(id_array[j].data[i].flag&2) k1++;
				if(id_array[j].data[i].flag&4) k2++;
			}
			sd=elem->arg.var->data;
			if(sd->vtype&ST_ARRAY) i1=20-fprintf(flog,"       %s(%d)",sd->name,elem->oindex);
			else i1=20-fprintf(flog,"       %s",sd->name);
			while((i1--)>0) (void)fputc(' ',flog);
			(void)fprintf(flog,"  %5d    %5d    %5d\n",k,k1,k2);
		}
		for(i=0;i<n_nonid_records;i++) if(nonid_elements[i]->type&ST_CENSORED) {
			elem=nonid_elements[i];
			if(!elem->type&(ST_MODEL|ST_TRAIT)) continue;
			for(k=k1=k2=j=0;j<ped_size;j++) if(ped_recode1[j] && id_array[j].data1) {
				for(rec=0;rec<id_array[j].nrec;rec++) if(id_array[j].data1[rec]) {
					if(id_array[j].data1[rec][i].flag&1) k++;
					if(id_array[j].data1[rec][i].flag&2) k1++;
					if(id_array[j].data1[rec][i].flag&4) k2++;
				}
			}
			sd=elem->arg.var->data;
			if(sd->vtype&ST_ARRAY) i1=20-fprintf(flog,"       %s(%d)",sd->name,elem->oindex);
			else i1=20-fprintf(flog,"       %s",sd->name);
			while((i1--)>0) (void)fputc(' ',flog);
			(void)fprintf(flog,"  %5d    %5d    %5d\n",k,k1,k2);
		}
	}
	if(n_markers) {
		j=0;
		linkp=links;
		while(linkp) {
			j++;
			linkp=linkp->next;
		}
		if(!(llist=malloc(sizeof(void *)*j))) ABT_FUNC(MMsg);
		nlink=j;
		linkp=links;
		while(linkp) {
			llist[--j]=linkp;
			linkp=linkp->next;
		}
		if(!(tmp=calloc((size_t)n_markers,sizeof(int)))) ABT_FUNC(MMsg);
		for(i=0;i<ped_size;i++) if(ped_recode1[i] && id_array[i].haplo[0]) {
			for(j=0;j<n_markers;j++)
			  if(id_array[i].haplo[0][j] || id_array[i].haplo[1][j]) tmp[j]++;
		}
		for(k1=0;k1<nlink;k1++)	{
			linkp=llist[k1];
			(void)fprintf(flog,"\n\n     Linkage Group - %s:\n\n       Marker      No. Alleles     No. Records      Het.\n",linkp->name?linkp->name:"__Not_Named__");
			(void)fprintf(flog,"       ------      -----------     -----------      ----\n");
			for(i=0;i<linkp->n_loci;i++) {
				j=linkp->element[i]->index;
				if(markers[j].var->vtype&ST_ARRAY) k=20-fprintf(flog,"       %s(%d)",markers[j].var->name,markers[j].index);
				else k=20-fprintf(flog,"       %s",markers[j].var->name);
				while((k--)>0) (void)fputc(' ',flog);
				(void)fprintf(flog,"%5d             %5d",linkp->element[i]->n_levels,tmp[j]);
				hom=het=0.0;
				for(k2=0;k2<ped_size;k2++) if(id_array[k2].haplo[0]) {
					if(id_array[k2].haplo[0][j] && id_array[k2].haplo[1][j]) {
						if(id_array[k2].haplo[0][j]==id_array[k2].haplo[1][j]) hom++;
						else het++;
					}
				}
				if(hom+het) het=het/(hom+het);
				(void)fprintf(flog,"       %7.5f\n",het);
			}
		}
		free(tmp);	
		for(k2=j=0;j<n_markers;j++) {
			k1=markers[j].element->n_levels;
			if(k1>k2) k2=k1;
		}
		if(k2) {
			if(!(tmp=malloc(k2*sizeof(int)))) ABT_FUNC(MMsg);
			for(j=0;j<n_markers;j++) if(markers[j].element->n_levels) {
				if(markers[j].var->vtype&ST_ARRAY) (void)fprintf(flog,"\n\n     Marker - %s(%d):\n\n       Allele       No. Records\n",markers[j].var->name,markers[j].index);
				else (void)fprintf(flog,"\n\n     Marker - %s:\n\n       Allele       No. Records\n",markers[j].var->name);
				(void)fputs("       ------       -----------\n",flog);
				k1=markers[j].element->n_levels;
				for(i=0;i<k1;i++) tmp[i]=0;
				i1=0;
				for(i=0;i<ped_size;i++) if(ped_recode1[i] && id_array[i].haplo[0]) {
					comp=id_array[i].component-1;
					for(k3=0;k3<2;k3++) {
						k2=id_array[i].haplo[k3][j];
						if(k2>i1) i1=k2;
						if(k2) {
							k2=markers[j].allele_trans[comp][k2-1];
							if(k2<0 || k2>=k1) ABT_FUNC("Internal error - illegal recoded allele value\n");
							tmp[k2]++;
						}
					}
				}
				for(i=0;i<k1;i++) {
					if(factor_recode[n_factors+j][i]->type==STRING)	{
						k=23-fprintf(flog,"        %s",factor_recode[n_factors+j][i]->data.string);
						while((k--)>0) (void)fputc(' ',flog);
					} else (void)fprintf(flog,"        %-15ld",factor_recode[n_factors+j][i]->data.value);
					(void)fprintf(flog,"%5d\n",tmp[i]);
				}
			}
			free(tmp);
		}
		free(llist);
	}
	if(flog!=stdout) (void)fclose(flog);
}
