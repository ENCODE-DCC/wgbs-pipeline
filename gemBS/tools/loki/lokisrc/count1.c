/****************************************************************************
 *                                                                          *
 *     Loki - Programs for genetic analysis of complex traits using MCMC    *
 *                                                                          *
 *                     Simon Heath - CNG                                    *
 *                                                                          *
 *                       April 2002                                         *
 *                                                                          *
 * count_relationships.c:                                                   *
 *                                                                          *
 * Counts different types of relationships                                  *
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
#include "control_parse.h"
#include "scan.h"

static struct relate_off **relate;

static double kin(int a,int b) 
{
	int i,j,ids,ids1;
	
	if(!(a&&b)) return 0.0;
	i=ped_recode1[a-1];
	j=ped_recode1[b-1];
	if(!(i&&j)) return 0.0;
	if(i==j) return .5*(1.0+kin(id_array[a-1].sire,id_array[a-1].dam));
	ids=id_array[a-1].sire;
	if(ids) ids=ped_recode1[ids-1];
	ids1=id_array[b-1].sire;
	if(ids1) ids1=ped_recode1[ids1-1];
	if(!ids) {
		if(j<i || !ids1) return 0.0;
		return .5*(kin(a,id_array[b-1].sire)+kin(a,id_array[b-1].dam));
	}
	if(!ids1) {
		if(i<j) return 0.0;
		return .5*(kin(b,id_array[a-1].sire)+kin(b,id_array[a-1].dam));
	}
	if(i<j) return .5*(kin(a,id_array[b-1].sire)+kin(a,id_array[b-1].dam));
	return .5*(kin(b,id_array[a-1].sire)+kin(b,id_array[a-1].dam));
}

void count_relationships(void)
{
	int i,j,j1,k,nr,i1,comp,cs,rel;
	int *counts,*counts1,*perm;
	int ids,ids1,ids2,ids3,ids4,ids5,idd,idd1,idd2,idd3,idd4,idd5;
	double z;
	char *relate[]={"full-sibs","half-sibs","parent:offspring","grandparent:grandchild","great-grandparent:great-grandchild",
		  "avuncular","half-avuncular",
		  "half first-cousins","first-cousins","first-cousins","double first-cousins",0};
	double coeff[]={.5,.25,.5,.25,.125,
		  .25,.125,
		  .0625,.125,.125,.25};
	char *msg,*unknown="unknown";
	
	nr=0;
	while(relate[nr]) nr++;
	if(!nr) return;
	if(!(counts=calloc((size_t)(2*nr+1),sizeof(int)))) ABT_FUNC(MMsg);
	counts1=counts+nr;
	if(!(perm=malloc(sizeof(int)*pruned_ped_size))) ABT_FUNC(MMsg);
	for(i=0;i<ped_size;i++)	{
		j=ped_recode1[i];
		if(j) perm[j-1]=i;
	}
	if(!(relate=malloc(sizeof(void *)*pruned_ped_size))) ABT_FUNC(MMsg);
	comp=id_array[perm[0]].component;
	cs=0;
	if(!comp_sflag || comp<n_comp) for(i1=1;i1<pruned_ped_size;i1++) {
		relate[i1]=0;
		i=perm[i1];
		if(id_array[i].component!=comp) {
			cs=i1;
			comp=id_array[i].component;
			if(comp_sflag && comp==n_comp) break;
			continue;
		}
		if(!id_array[i].sire || !ped_recode1[id_array[i].sire-1]) continue;
		for(j1=cs;j1<i1;j1++) {
			j=perm[j1];
			z=2.0*kin(i+1,j+1);
			if(z>0.0) {
				rel=-1;
				ids=id_array[i].sire;
				ids1=id_array[j].sire;
				idd=id_array[i].dam;
				idd1=id_array[j].dam;
				if(ids==ids1) {
					if(idd==idd1) rel=0; /* Full sibs */
					else rel=1; /* Paternal half sibs */
				} else if(idd==idd1) rel=1; /* Maternal half sibs */
				else if(ids==j+1 || idd==j+1) rel=2; /* Parent-offspring */
				else {
					if(ids>j+1) {
						ids2=id_array[ids-1].sire;
						idd2=id_array[ids-1].dam;
						if(ids2==j+1 || idd2==j+1) rel=3; /* grandparents */
						else if(ids2>j+1 && (id_array[ids2-1].sire==j+1 || id_array[ids2-1].dam==j+1)) rel=4; /* great-grandparents */
						else if(idd2>j+1 && (id_array[idd2-1].sire==j+1 || id_array[idd2-1].dam==j+1)) rel=4; /* great-grandparents */
					}
					if(rel<0 && idd>j+1) {
						ids2=id_array[idd-1].sire;
						idd2=id_array[idd-1].dam;
						if(ids2==j+1 || idd2==j+1) rel=3; /* grandparents */
						else if(ids2>j+1 && (id_array[ids2-1].sire==j+1 || id_array[ids2-1].dam==j+1)) rel=4; /* great-grandparents */
						else if(idd2>j+1 && (id_array[idd2-1].sire==j+1 || id_array[idd2-1].dam==j+1)) rel=4; /* great-grandparents */
					}
				}
				if(rel<0) { /* Paternal avuncular */
					ids2=id_array[ids-1].sire;
					idd2=id_array[ids-1].dam;
					if(ids2==ids1) {
						if(idd2==idd1) rel=5; /* avuncular */
						else rel=6; /* half-avuncular */
					} else if(idd2==idd1) rel=6;
				}
				if(rel<0) { /* Maternal avuncular */
					ids2=id_array[idd-1].sire;
					idd2=id_array[idd-1].dam;
					if(ids2==ids1) {
						if(idd2==idd1) rel=5; /* avuncular */
						else rel=6; /* half-avuncular */
					} else if(idd2==idd1) rel=6;
				}
				if(rel<0) { /* First cousins */
					ids2=id_array[ids-1].sire;
					idd2=id_array[ids-1].dam;
					ids3=id_array[idd-1].sire;
					idd3=id_array[idd-1].dam;
					ids4=id_array[ids1-1].sire;
					idd4=id_array[ids1-1].dam;
					ids5=id_array[idd1-1].sire;
					idd5=id_array[idd1-1].dam;
					k=0;
					if(ids2 && (ids2==ids4 || ids2==ids5)) k++; /* Cousins */
					if(ids3 && (ids3==ids4 || ids3==ids5)) k++;
					if(idd2 && (idd2==idd4 || idd2==idd5)) k++;
					if(idd3 && (idd3==idd4 || idd3==idd5)) k++;
					if(k) rel=6+k;
				}
				print_orig_id(stdout,j+1,1);
				print_orig_id(stdout,i+1,1);
				if(rel<0) msg=unknown;
				else msg=relate[rel];
				printf("%g %s",z,msg);
				if(rel>=0 && z>coeff[rel]) fputc('+',stdout);
				fputc('\n',stdout);
			}
		}
	}
	free(perm);
	free(counts);
}
