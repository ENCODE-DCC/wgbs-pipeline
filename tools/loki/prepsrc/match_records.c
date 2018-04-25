/****************************************************************************
 *                                                                          *
 *     Loki - Programs for genetic analysis of complex traits using MCMC    *
 *                                                                          *
 *             Simon Heath - University of Washington                       *
 *                                                                          *
 *                       March 1997                                         *
 *                                                                          *
 * match_records.c:                                                         *
 *                                                                          *
 * Matches records to id.                                                   *
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
#include "y.tab.h"

int n_id_records,n_nonid_records;
struct var_element **id_elements=0,**nonid_elements=0;

void match_records(void)
{
	int i,i1,j,k,k1,k2,k3,id,ncol;
	struct InFile *infile,*infile1;
	struct var_element *elem;
	struct id_data *data,**dataptr=0;
	struct DataBlock *db;
	struct scan_data *sd;
	
	n_id_records=n_nonid_records=0;
	/* Count constant and inconstant records */
	infile=Infiles;
	while(infile) {
		for(i=0;i<infile->nvar;i++) {
			elem=infile->element[i];
			if(elem && !(elem->type&(ST_ID|ST_FAMILY|ST_SIRE|ST_DAM|ST_HAPLO|ST_FLAG|ST_MARKER)))	{
				if(elem->type&ST_CONSTANT) n_id_records++;
				else n_nonid_records++;
				elem->type|=ST_FLAG;
			}
		}
		infile=infile->next;
	}
	/* Allocate space for list, and put element pointers for
	 * constant and inconstant records in list */
	if(!(n_id_records+n_nonid_records)) return;
	if(!(id_elements=malloc((n_id_records+n_nonid_records)*sizeof(void *)))) ABT_FUNC(MMsg);
	nonid_elements=id_elements+n_id_records;
	infile=Infiles;
	k=j=0;
	while(infile) {
		for(i=0;i<infile->nvar;i++) {
			elem=infile->element[i];
			if(elem && !(elem->type&(ST_ID|ST_FAMILY|ST_SIRE|ST_DAM|ST_HAPLO))) {
				if(elem->type&ST_FLAG) {
					if(elem->type&ST_CONSTANT) id_elements[j++]=elem;
					else nonid_elements[k++]=elem;
					elem->type&=~ST_FLAG;
				}
			}
		}
		infile=infile->next;
	}
	/* Find out which individuals have data */
	for(i=0;i<ped_size;i++) id_array[i].flag=id_array[i].nrec=0;
	infile=Infiles;
	while(infile) {
		db=infile->data;
		ncol=infile->ncol;
		while(db) {
			for(j=0;j<db->record_ptr;j++)	{
				id=(int)db->records[j*ncol+infile->id_col].value;
				k1=n_nonid_records?0:1;
				i=-1;
				for(i1=0;i1<infile->nvar;i1++) {
					elem=infile->element[i1];
					if(!elem) continue;
					i++;
					if(check_missing(i,ncol,j,db)) continue;
					if(!id_array[id-1].flag) for(k=0;k<n_id_records;k++) {
						if(elem==id_elements[k]) {
							id_array[id-1].flag=1;
							break;
						}
					}
					if(!k1) for(k=0;k<n_nonid_records;k++) {
						if(elem==nonid_elements[k]) {
							id_array[id-1].nrec++;
							k1=1;
							break;
						}
					}
					if(k1 && id_array[id-1].flag) break;
				}
			}
			db=db->next;
		}
		infile=infile->next;
	}
	/* Count how many data records we have */
	for(i=j=k=0;i<ped_size;i++) {
		if(id_array[i].flag) j++;
		k+=id_array[i].nrec;
	}
	/* Allocate space for data records */
	if(!(j+k)) return;
	if(!(data=calloc((size_t)(k*n_nonid_records+j*n_id_records),sizeof(struct id_data)))) ABT_FUNC(MMsg);
	RemBlock=AddRemem(data,RemBlock);
	if(k) {
		if(!(dataptr=malloc(k*sizeof(void *)))) ABT_FUNC(MMsg);
		RemBlock=AddRemem(dataptr,RemBlock);
	}
	for(i=0;i<ped_size;i++)	{
		if(id_array[i].flag)	{
			id_array[i].data=data;
			data+=n_id_records;
		} else id_array[i].data=0;
		if(id_array[i].nrec)	{
			id_array[i].data1=dataptr;
			for(k=0;k<id_array[i].nrec;k++) {
				dataptr[k]=data;
				data+=n_nonid_records;
			}
			dataptr+=id_array[i].nrec;
			id_array[i].nrec=0;
		} else id_array[i].data1=0;
	}
	/* Go through data again and assign to individuals */
	infile=Infiles;
	while(infile) {
		db=infile->data;
		ncol=infile->ncol;
		while(db) {
			for(j=0;j<db->record_ptr;j++)	{
				id=(int)db->records[j*ncol+infile->id_col].value;
				i= -1;
				for(k1=i1=0;i1<infile->nvar;i1++) {
					elem=infile->element[i1];
					if(!elem) continue;
					i++;
					if(check_missing(i,ncol,j,db)) continue;
					for(k=0;k<n_id_records;k++) if(elem==id_elements[k]) {
						if(id_array[id-1].data[k].flag) {
							k2=0;
							if(elem->type&(ST_INTTYPE|ST_FACTOR)) {
								if(id_array[id-1].data[k].data.value!=db->records[j*ncol+i].value) k2=1;
							} else if(id_array[id-1].data[k].data.rvalue!=db->records[j*ncol+i].rvalue) k2=1;
							if(k2) {
								sd=elem->arg.var->data;
								if(sd->vtype&ST_ARRAY) (void)fprintf(stderr,"Error: File %s (col %d), variable '%s(%d)', id ",infile->name,i1+1,sd->name,elem->oindex);
								else (void)fprintf(stderr,"Error: File %s, variable '%s', id ",infile->name,sd->name);
								print_orig_id(stderr,id,0);
								(void)fputs(" - constant type isn't\n",stderr);
								if(elem->type&ST_INTTYPE) (void)fprintf(stderr,"(old: %ld, new: %ld)\n",id_array[id-1].data[k].data.value,db->records[j*ncol+i].value);
								else if(elem->type&ST_FACTOR) {
									for(k2=0;k2<n_factors;k2++) if(elem==var_factors[k2]) break;
									if(k2==n_factors) ABT_FUNC("Internal error - factor not found\n");
									k3=id_array[id-1].data[k].data.value-1;
									if(factor_recode[k2][k3]->type==STRING) (void)fprintf(stderr,"(old: '%s', ",factor_recode[k2][k3]->data.string);
									else (void)fprintf(stderr,"(old: %ld, ",factor_recode[k2][k3]->data.value);
									k3=db->records[j*ncol+i].value-1;
									if(factor_recode[k2][k3]->type==STRING) (void)fprintf(stderr,"new: '%s')\n",factor_recode[k2][k3]->data.string);
									else (void)fprintf(stderr,"new: %ld)\n",factor_recode[k2][k3]->data.value);
								} else (void)fprintf(stderr,"(old: %g, new: %g)\n",id_array[id-1].data[k].data.rvalue,db->records[j*ncol+i].rvalue);
								if((++scan_error_n)>=max_scan_errors) abt(__FILE__,__LINE__,"Too many errors - aborting\n");
							}
						}
						id_array[id-1].data[k].flag=1|(elem->type&(ST_INTTYPE|ST_REALTYPE));
						if(elem->type&ST_INTTYPE) id_array[id-1].data[k].data.value=db->records[j*ncol+i].value;
						else id_array[id-1].data[k].data.rvalue=db->records[j*ncol+i].rvalue;
						break;
					}
					for(k=0;k<n_nonid_records;k++) if(elem==nonid_elements[k]) {
						k2=id_array[id-1].nrec;
						if(id_array[id-1].data1[k2][k].flag) ABT_FUNC("Internal error - record already used\n");
						id_array[id-1].data1[k2][k].flag=1|(elem->type&(ST_INTTYPE|ST_REALTYPE));
						if(elem->type&ST_INTTYPE) id_array[id-1].data1[k2][k].data.value=db->records[j*ncol+i].value;
						else id_array[id-1].data1[k2][k].data.rvalue=db->records[j*ncol+i].rvalue;
						k1=1;
						break;
					}
				}
				if(k1) id_array[id-1].nrec++;
			}
			db=db->next;
		}
		infile1=infile->next;
		free_infile(infile);
		infile=infile1;
	}
	Infiles=0;
	for(i=0;i<ped_size;i++) id_array[i].flag=0;
}
