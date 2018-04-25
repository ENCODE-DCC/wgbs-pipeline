/****************************************************************************
 *                                                                          *
 *     Loki - Programs for genetic analysis of complex traits using MCMC    *
 *                                                                          *
 *             Simon Heath - University of Washington                       *
 *                                                                          *
 *                      March/April 1997                                    *
 *                                                                          *
 * calc_nrm.c:                                                              *
 *                                                                          *
 * Calculate Inverse of NRM matrix using algorithm of Quaas (1976)          *
 *                                                                          *
 * June 2003 (SCH) - Moved from prep to loki                                * 
 *                                                                          *
 * Copyright (C) Simon C. Heath 1997, 2000, 2002, 2003                      *
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
#include "libhdr.h"
#include "loki.h"
#include "loki_utils.h"
#include "sparse.h"

static struct loki *lk;

static void free_nrm(void)
{
	int i;
	
	if(lk && lk->models->AIMatrix) {
		for(i=0;i<lk->pedigree->n_comp;i++) {
			if(lk->models->AIMatrix[i]) free(lk->models->AIMatrix[i]);
		}
		free(lk->models->AIMatrix);
	}
}

/* Calculate inverse of NRM (G) matrix using Quaas/Henderson method 
 * See Numerical Recipes in C 2nd edition pps. 78-79 for a description 
 * of the sparse matrix storage. Note the matrix is half stored.
 * 
 * The strategy is to build up the matrix directly, which means knowing
 * what non-zero diagonal elements will be present.  This is relatively
 * simple for the inverse NRM matrix as it has a very simple structure */
void Calculate_NRM(struct loki *loki)
{
	int id,i,j,k,k1,idd,ids,idd1,ids1,cs,nz=0,comp,max_comp=0;
	int *sire_list,*dam_list,*temp,*spouses=0,n_comp;
	double d,*v,*u,xx,d2,d4,Detl,TDetl;
	struct SparseMatRec *AIM=0;
	struct Id_Record *id_array;
	
	message(INFO_MSG,"Calculating NRM matrix for polygenic effect\n");
	lk=loki;
	/* Find maximum component size */
	n_comp=loki->pedigree->n_comp;
	for(i=0;i<n_comp;i++) if(loki->pedigree->comp_size[i]>max_comp) max_comp=loki->pedigree->comp_size[i];
	if(!(sire_list=calloc((size_t)(loki->pedigree->ped_size+max_comp),sizeof(int)))) ABT_FUNC(MMsg);
	dam_list=sire_list+max_comp;
	if(!(u=malloc(max_comp*2*sizeof(double)))) ABT_FUNC(MMsg);
	v=u+max_comp;
	TDetl=0.0;
	if(!(loki->models->AIMatrix=malloc(sizeof(void *)*n_comp))) ABT_FUNC(MMsg);
	id_array=loki->pedigree->id_array;
	/* Do this one component at a time */
	for(id=comp=0;comp<n_comp;comp++) {
		cs=loki->pedigree->comp_size[comp];
		/* Calculate storage requirements */
		if(!(temp=calloc((size_t)cs+1,sizeof(int)))) ABT_FUNC(MMsg);
		/* Count non-zero (off diagonal) elements in matrix */
		for(j=k=0;j<cs;j++) {
			ids=id_array[id+j].sire;
			idd=id_array[id+j].dam;
			if(ids) ids-=id;
			if(idd) idd-=id;
			sire_list[j]=ids;
			dam_list[j]=idd;
			if(ids!=idd) {
				k++;
				if(ids>idd) temp[ids-1]++;
				else temp[idd-1]++;
			}
		}
		nz=0;
		spouses=0;
		if(k) {
			if(!(spouses=calloc((size_t)k,sizeof(int)))) ABT_FUNC(MMsg);
			for(j=k=0;j<cs;j++)	{
				k1=k;
				k+=temp[j];
				temp[j]=k1;
				ids=sire_list[j];
				idd=dam_list[j];
				if(ids && ids!=(j+1)) nz++;
				if(idd && idd!=ids && idd!=(j+1)) nz++;
				if(ids!=idd) {
					if(ids>idd)	{
						for(k1=temp[ids-1];k1<temp[ids];k1++) {
							if(spouses[k1]==idd) break;
							if(!spouses[k1]) {
								spouses[k1]=idd;
								nz++;
								break;
							}
						}
					} else {
						for(k1=temp[idd-1];k1<temp[idd];k1++) {
							if(spouses[k1]==ids) break;
							if(!spouses[k1]) {
								spouses[k1]=ids;
								nz++;
								break;
							}
						}
					}
				}
			}
		}
		if(!(AIM=malloc((nz+1+cs)*sizeof(struct SparseMatRec)))) ABT_FUNC(MMsg);
		/* Set up coordinates for non-zero elements */
		nz=cs+1;
		for(j=0;j<cs;j++) {
			u[j]=0.0;
			AIM[j].val=0.0;
			AIM[j].x=nz;
			ids=sire_list[j];
			idd=dam_list[j];
			if(ids && ids!=(j+1)) {
				AIM[nz].x=ids-1;
				AIM[nz++].val=0.0;
			}
			if(idd && idd!=ids && idd!=(j+1)) {
				AIM[nz].x=idd-1;
				AIM[nz++].val=0.0;
			}
			if(!spouses) continue;
			for(k=temp[j];k<temp[j+1];k++) {
				if(!spouses[k]) break;
				AIM[nz].x=spouses[k]-1;
				AIM[nz++].val=0.0;
			}
		}
		AIM[cs].x=nz;
		if(spouses) free(spouses);
		Detl=0.0;
		for(i=0;i<cs;i++) {
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
			for(k=i+1;k<cs;k++) {
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
			AIM[i].val+=d;
			d2= -.5*d;
			d4=.25*d;
			j=AIM[i].x;
			if(ids>=0) {
				AIM[j++].val+=d2;
				AIM[ids].val+=d4;
			}
			if(idd>=0) {
				AIM[j++].val+=d2;
				AIM[idd].val+=d4;
				if(ids>=0) {
					if(ids==idd) AIM[ids].val+=d4;
					else if(ids>idd) {
						for(j=AIM[ids].x;j<AIM[ids+1].x;j++) if(AIM[j].x==idd) {
							AIM[j].val+=d4;
							break;
						}
					} else {
						for(j=AIM[idd].x;j<AIM[idd+1].x;j++) if(AIM[j].x==ids) {
							AIM[j].val+=d4;
							break;
						}
					}
				}
			}
		}
		Detl+=Detl;
		TDetl+=Detl;
		AIM[cs].val=Detl;
		loki->models->AIMatrix[comp]=AIM;
		message(DEBUG_MSG," Component %d, non-zero off-diagonals = %d, L(Det) NRM Matrix = %g\n",comp+1,nz-1-cs,Detl);
		free(temp);
		id+=cs;
	}
	free(sire_list);
	free(u);
	if(atexit(free_nrm)) message(WARN_MSG,"Unable to register exit function free_nrm()\n");
}
