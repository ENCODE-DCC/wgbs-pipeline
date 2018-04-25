/****************************************************************************
 *                                                                          *
 *     Loki - Programs for genetic analysis of complex traits using MCMC    *
 *                                                                          *
 *             Simon Heath - University of Washington                       *
 *                                                                          *
 *                      March/April 1997                                    *
 *                                                                          *
 * check_inbr.c:                                                            *
 *                                                                          *
 * Calculate inbreeding coefficients (info only)                            *
 *                                                                          *
 * Copyright (C) Simon C. Heath 1997, 2000, 2002, 2003                      *
 *                                                                          *
 * Gutted version of calc_nrm.c, which has been moved to loki               *
 * June 2003 (SCH)                                                          *
 *                                                                          *
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
#include "prep_utils.h"
#include "libhdr.h"
#include "scan.h"

static double *u;

static int cmp_inbr(const void *s1,const void *s2)
{
	double x1,x2;
	int i;
	
	i=*((const int *)s1);
	x1=u[i];
	i=*((const int *)s2);
	x2=u[i];
	if(x1<x2) return 1;
	if(x1>x2) return -1;
	return 0;
}

void Check_Inbr(char *LogFile)
{
	int id,i,j,k,k1,idd,ids,idd1,ids1,comp,max_comp=0,n_inbr;
	int *sire_list,*dam_list,*temp,*perm;
	FILE *flog;
	double d,*v,xx,Detl,TDetl,avg_inbr;

	if(!pruned_ped_size || !LogFile || !(LogFile=add_file_dir(LogFile)))	return;
	/* Find maximum component size */
	for(i=0;i<n_comp;i++) if(comp_size[i]>max_comp) max_comp=comp_size[i];
	if(!(sire_list=calloc((size_t)(pruned_ped_size+1+max_comp*3),sizeof(int)))) ABT_FUNC(MMsg);
	dam_list=sire_list+max_comp;
	perm=dam_list+max_comp;
	temp=perm+pruned_ped_size;
	if(!(u=malloc(max_comp*2*sizeof(double)))) ABT_FUNC(MMsg);
	v=u+max_comp;
	for(i=0;i<ped_size;i++)	{
		j=ped_recode1[i];
		if(j) perm[j-1]=i;
	}
	TDetl=0.0;
	flog=fopen(LogFile,"a");
	free(LogFile);
	(void)fputs("\n********* Calculation of inbreeding coefficients ********\n\n",flog);
	(void)fprintf(flog,"     Pedigree size = %d\n\n",pruned_ped_size);
	/* Do the calculations one component at a time */
	for(id=comp=0;comp<n_comp;comp++) {
		(void)fflush(stdout);
		for(j=0;j<comp_size[comp];j++)	{
			i=perm[j+id];
			ids=id_array[i].sire;
			idd=id_array[i].dam;
			if(ids) ids=ped_recode1[ids-1]-id;
			if(idd) idd=ped_recode1[idd-1]-id;
			sire_list[j]=ids;
			dam_list[j]=idd;
		}
		/* Zero inbreeding count */
		avg_inbr=0.0;
		n_inbr=0;
		for(j=0;j<comp_size[comp];j++) u[j]=0.0;
		Detl=0.0;
		for(i=0;i<comp_size[comp];i++) {
			xx=0.0;
			ids=sire_list[i]-1;
			idd=dam_list[i]-1;
			if(ids>=0) xx=u[ids];
			if(idd>=0) xx+=u[idd];
			xx=1.0-.25*xx;
			d=1.0/xx;
			u[i]+=xx;
			xx=sqrt(xx);
			v[i]=xx;
			Detl+=log(xx);
			for(k=i+1;k<comp_size[comp];k++) {
				xx=0.0;
				ids1=sire_list[k]-1;
				idd1=dam_list[k]-1;
				if(ids1>=i) xx=v[ids1];
				if(idd1>=i) xx+=v[idd1];
				if(xx>0.0) {
					xx=.5*xx;
					u[k]+=xx*xx;
				}
				v[k]=xx;
			}
		}
		for(j=0;j<comp_size[comp];j++) {
			if(u[j]>1.0) {
				xx=u[j]-1.0;
				avg_inbr+=xx;
				n_inbr++;
			}
		}
		Detl+=Detl;
		TDetl+=Detl;
		(void)fprintf(flog,"     Component %d\n",comp+1);
		(void)fprintf(flog,"       Component size = %d\n",comp_size[comp]);
		(void)fprintf(flog,"       Log Determinant of NRM matrix = %f\n",Detl);
		(void)fprintf(flog,"       No. inbred individuals = %d\n",n_inbr);
		if(n_inbr) {
			(void)fprintf(flog,"       Average inbreeding coefficient (overall) = %g\n",avg_inbr/(double)comp_size[comp]);
			(void)fprintf(flog,"       Average inbreeding coefficient (inbred individuals only) = %g\n",avg_inbr/(double)n_inbr);
			(void)fputs("\n       Ids and inbreeding coefficients of inbred individuals:\n\n",flog);
			for(k=j=0;j<comp_size[comp];j++) if(u[j]>1.0) temp[k++]=j;
			gnu_qsort(temp,(size_t)k,sizeof(int),cmp_inbr); 
			for(j=0;j<n_inbr;j++) {
				k1=temp[j];
				xx=u[k1]-1.0;
				k=perm[id+k1];
				(void)fputs("          ",flog);
				print_orig_id(flog,k+1,0);
				(void)fprintf(flog," %g\n",xx);
			}
		}
		id+=comp_size[comp];
	}
	free(sire_list);
	free(u);
	if(n_comp>1) (void)fprintf(flog,"\n     Log Determinant of combined NRM matrix = %f\n",TDetl);
	(void)fclose(flog);
}
