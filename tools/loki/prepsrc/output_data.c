/****************************************************************************
 *                                                                          *
 *     Loki - Programs for genetic analysis of complex traits using MCMC    *
 *                                                                          *
 *                     Simon Heath - CNG, Paris                             *
 *                                                                          *
 *                           August 2002                                    *
 *                                                                          *
 * output_data.c:                                                           *
 *                                                                          *
 * Output the data in text format (not used by Loki)                        *
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
#include <stdio.h>

#include "utils.h"
#include "scan.h"
#include "y.tab.h"

void Output_Data(void)
{
	int i,j,ids,idd,*perm;
	FILE *fptr;
	
	if(!pruned_ped_size) return;
	if((fptr=fopen(OutputFile,"w"))) {
		if(!(perm=malloc(pruned_ped_size*sizeof(int)))) ABT_FUNC(MMsg);
		for(i=0;i<ped_size;i++)	{
			j=ped_recode1[i];
			if(j) perm[j-1]=i;
		}
		for(j=0;j<pruned_ped_size;j++) {
			i=perm[j];
			(void)fprintf(fptr,"%d ",id_array[i].component);
			print_orig_triple(fptr,i+1);
			ids=id_array[i].sire;
			idd=id_array[i].dam;
			if(ids) ids=ped_recode1[ids-1];
			if(idd) idd=ped_recode1[idd-1];
			(void)fprintf(fptr,"%d %d %d %d",id_array[i].sex,j+1,ids,idd);
			if(Affected) (void)fprintf(fptr," %d",id_array[i].affected);
			(void)fputc('\n',fptr);
		}
		(void)fclose(fptr);
		free(perm);
	}
}

void Output_Raw_Data(void)
{
	int i,j,k,x,ch,*perm;
	FILE *fptr;
	
	if(!pruned_ped_size) return;
	if((fptr=fopen(OutputRawFile,"w"))) {
		if(!(perm=malloc(pruned_ped_size*sizeof(int)))) ABT_FUNC(MMsg);
		for(i=0;i<ped_size;i++)	{
			j=ped_recode1[i];
			if(j) perm[j-1]=i;
		}
		for(j=0;j<pruned_ped_size;j++) {
			i=perm[j];
			print_orig_triple(fptr,i+1);
			if(id_array[i].haplo[0]) {
				for(x=0;x<n_markers;x++) {
					for(k=0;k<2;k++) {
						ch=id_array[i].haplo[k][x];
						if(ch--) {
							if(factor_recode[n_factors+x][ch]->type==STRING) {
								(void)fputc(' ',fptr);
								(void)fputs(factor_recode[n_factors+x][ch]->data.string,fptr);
							} else (void)fprintf(fptr," %ld",factor_recode[n_factors+x][ch]->data.value);
						} else (void)fputs(" *",fptr);
					}
				}
			} else for(x=0;x<n_markers;x++) fputs(" * *",fptr);
			(void)fputc('\n',fptr);
		}
		(void)fclose(fptr);
		free(perm);
	}
}
