/****************************************************************************
 *                                                                          *
 *     Loki - Programs for genetic analysis of complex traits using MCMC    *
 *                                                                          *
 *             Simon Heath - University of Washington                       *
 *                                                                          *
 *                        July 1997                                         *
 *                                                                          *
 * init_fam.c:                                                              *
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

int n_families=0;
struct Family *family=0;

void InitFamilies(char *LogFile)
{
	int i,j,k,kid,ids,idd,*tp,nk,nf;
	FILE *flog;
	
	for(i=0;i<ped_size;i++) if(ped_recode1[i]) {
		id_array[i].flag=id_array[i].nfam=id_array[i].nkids=id_array[i].family=0;
		id_array[i].kids=0;
	}
	/* Count kids */
	for(i=j=0;i<ped_size;i++) if(ped_recode1[i]) {
		ids=id_array[i].sire;
		idd=id_array[i].dam;
		if(ids) {
			id_array[ids-1].nkids++;
			j++;
		}
		if(idd) {
			id_array[idd-1].nkids++;
			j++;
		}
	}
	if(j)	{
		if(!(tp=malloc(sizeof(int)*j))) ABT_FUNC(MMsg);
		RemBlock=AddRemem(tp,RemBlock);
		for(i=0;i<ped_size;i++) if(ped_recode1[i]) {
			j=id_array[i].nkids;
			if(j)	{
				id_array[i].kids=tp;
				tp+=j;
				id_array[i].nkids=0;
			}
		}
		for(i=0;i<ped_size;i++) if(ped_recode1[i]) {
			ids=id_array[i].sire;
			idd=id_array[i].dam;
			if(ids) {
				j=id_array[ids-1].nkids++;
				id_array[ids-1].kids[j]=i;
			}
			if(idd) {
				j=id_array[idd-1].nkids++;
				id_array[idd-1].kids[j]=i;
			}
		}
	}
	for(nf=i=0;i<ped_size;i++) if(ped_recode1[i] && !id_array[i].flag) {
		ids=id_array[i].sire;
		idd=id_array[i].dam;
		if(ids||idd) {
			j=ids?ids:idd;
			for(k=0;k<id_array[j-1].nkids;k++) {
				kid=id_array[j-1].kids[k];
				if(id_array[kid].sire==ids && id_array[kid].dam==idd) id_array[kid].flag=1;
			}
			n_families++;
			nf++;
			if(ids) id_array[ids-1].nfam++;
			if(idd) id_array[idd-1].nfam++;
		} else if(!id_array[i].nkids) n_families++;
	}
	if(LogFile && (LogFile=add_file_dir(LogFile))) {
		if((flog=fopen(LogFile,"a"))) {
			(void)fputs("\n******************* Family Initialization ***************\n\n",flog);
			(void)fprintf(flog,"     No. Full Sib Families = %d\n",nf);
			if(n_families-nf) (void)fprintf(flog,"     No. Single individuals = %d\n",n_families-nf);
			(void)fputc('\n',flog);
			(void)fclose(flog);
		}
		free(LogFile);
	}
	(void)printf("No. Full Sib Families = %d\n",nf);
	if(n_families-nf) (void)printf("No. Single individuals = %d\n",n_families-nf);
	if(nf) {
		if(!(family=calloc((size_t)nf,sizeof(struct Family)))) ABT_FUNC(MMsg);
		if(!(tp=malloc(sizeof(int)*2*nf))) ABT_FUNC(MMsg);
		RemBlock=AddRemem(tp,RemBlock);
	} else tp=0;
	for(i=0;i<ped_size;i++) if(ped_recode1[i]) {
		if(id_array[i].nfam) {
			id_array[i].famlist=tp;
			tp+=id_array[i].nfam;
			id_array[i].nfam=0;
		} else id_array[i].famlist=0;
	}
	n_families=0;
	for(nk=i=0;i<ped_size;i++) if(ped_recode1[i] && !id_array[i].family) {
		ids=id_array[i].sire;
		idd=id_array[i].dam;
		if(ids||idd) {
			j=ids?ids:idd;
			for(k=0;k<id_array[j-1].nkids;k++) {
				kid=id_array[j-1].kids[k];
				if(id_array[kid].sire==ids && id_array[kid].dam==idd) {
					id_array[kid].family=n_families+1;
					family[n_families].nkids++;
				}
			}
			nk+=family[n_families].nkids;
			if(ids) id_array[ids-1].famlist[id_array[ids-1].nfam++]=n_families;
			if(idd) id_array[idd-1].famlist[id_array[idd-1].nfam++]=n_families;
			family[n_families].sire=ids;
			family[n_families++].dam=idd;
		} else if(!id_array[i].nkids)	{
/*			nk++;
			family[n_families].nkids=1;
			family[n_families].sire=0;
			family[n_families++].dam=0;
			id_array[i].family=n_families; */
		}
	}
	if(nk) {
		if(!(tp=malloc(sizeof(int)*nk))) ABT_FUNC(MMsg);
		RemBlock=AddRemem(tp,RemBlock);
	}
	for(i=0;i<n_families;i++) {
		family[i].kids=tp;
		tp+=family[i].nkids;
		family[i].nkids=0;
	}
	for(i=0;i<ped_size;i++) if(ped_recode1[i] && (j=id_array[i].family))
	  family[j-1].kids[family[j-1].nkids++]=i;
}

