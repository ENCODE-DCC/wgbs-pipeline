/****************************************************************************
 *                                                                          *
 *     Loki - Programs for genetic analysis of complex traits using MCMC    *
 *                                                                          *
 *                    Simon Heath - MSKCC                                   *
 *                                                                          *
 *                       October 2001                                       *
 *                                                                          *
 * count_dbr.c:                                                             *
 *                                                                          *
 * Routines for couting and displaying double recombinant probabilities     *
 *                                                                          *
 * Copyright (C) Simon C. Heath 1997, 2000, 2002                            *
 * This is free software.  You can distribute it and/or modify it           *
 * under the terms of the Modified BSD license, see the file COPYING        *
 *                                                                          *
 ****************************************************************************/

#include <config.h>
#include <stdlib.h>
#include <string.h>
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif
#ifdef USE_DMALLOC
#include <dmalloc.h>
#endif
#include <stdio.h>
#include <math.h>
#include <sys/types.h>
#include <sys/ipc.h>
#include <sys/shm.h>

#include "utils.h"
#include "loki.h"
#include "loki_peel.h"
#include "loki_monitor.h"
#include "count_dbr.h"

static unsigned int **dbr_count[2],**dbr_countn[2],*dbr_mem;
static int *perm;

int init_dbr_shm(void)
{
	int i,j,k;
	
	int er=-1;
	
	if(n_links && n_markers) {
		for(k=i=0;i<n_links;i++) {
			j=linkage[i].n_markers;
			if(j>2) k+=j-2;
		}
		if(k) {
			lpar->dbr_mem_size=k*ped_size*4*sizeof(int);
			lpar->dbr_shm_id=shmget(IPC_PRIVATE,lpar->dbr_mem_size,IPC_CREAT|0600);
			if(lpar->dbr_shm_id<0) {
				perror("Failed to get shared memory segment");
				ABT_FUNC(AbMsg);
			}
		} 
		er=0;
	}
	return er;
}

int init_dbr_count(void)
{
	int i,j,k,k1,k2,er=-1;
	unsigned int *tmp;
	
	if(lpar->dbr_shm_id>=0) {
		dbr_mem=shmat(lpar->dbr_shm_id,0,0);
		if(dbr_mem==(void *)-1) {
			perror("Failed to attach shared memory segment");
			ABT_FUNC(AbMsg);
		}
		if(!(perm=malloc(sizeof(int)*n_markers))) ABT_FUNC(MMsg);
		if(!(dbr_count[0]=malloc(sizeof(void *)*4*n_markers))) ABT_FUNC(MMsg);
		dbr_count[1]=dbr_count[0]+n_markers;
		dbr_countn[0]=dbr_count[1]+n_markers;
		dbr_countn[1]=dbr_countn[0]+n_markers;
		for(i=0;i<n_markers*4;i++) dbr_count[0][i]=0;
		tmp=dbr_mem;
		for(i=0;i<n_links;i++) {
			get_locuslist(perm,i,&j,1);
			if(j>2) {
				gnu_qsort(perm,(size_t)j,sizeof(int),cmp_loci_pos);
				for(k=1;k<j-1;k++) {
					k1=perm[k];
					for(k2=0;k2<2;k2++) {
						dbr_count[k2][k1]=tmp;
						tmp+=ped_size;
						dbr_countn[k2][k1]=tmp;
						tmp+=ped_size;
					}
				}
			}
		}
		er=0;
	}
	return er;
}

void free_dbr_count(void)
{
	if(lpar->dbr_shm_id>=0) {
		(void)shmctl(lpar->dbr_shm_id,IPC_RMID,0);
		lpar->dbr_shm_id=-1;
		lpar->dbr_flag=0;
	}
	if(perm) free(perm);
	if(dbr_count[0]) free(dbr_count[0]);
	dbr_count[0]=0;
	dbr_mem=0;
	perm=0;
}

void zero_dbr_count(void) 
{
	if(dbr_mem) memset(dbr_mem,0,lpar->dbr_mem_size);
}

void count_dbr(void) 
{
	int i,j,k,k1,k2,*kk,id,seg[2][3];
	
	if(dbr_mem) {
		for(i=0;i<n_links;i++) {
			get_locuslist(perm,i,&j,1);
			if(j>2) {
				gnu_qsort(perm,(size_t)j,sizeof(int),cmp_loci_pos);
				for(id=0;id<ped_size;id++) if(id_array[id].sire) {
					for(k=0;k<2;k++) {
						seg[k][0]=marker[perm[0]].locus.seg[k][id];
						seg[k][1]=marker[perm[1]].locus.seg[k][id];
					}
					for(k=1;k<j-1;k++) {
						for(k1=0;k1<2;k1++) seg[k1][2]=marker[perm[k+1]].locus.seg[k1][id];
						for(k2=0;k2<2;k2++) {
							kk=seg[k2];
							for(k1=0;k1<3;k1++) if(kk[k1]<0) break;
							if(k1==3) {
								k1=perm[k];
								if(kk[0]!=kk[1] && kk[1]!=kk[2]) dbr_count[k2][k1][id]++;
								dbr_countn[k2][k1][id]++;
							}
						}
						for(k1=0;k1<2;k1++) {
							seg[k1][0]=seg[k1][1];
							seg[k1][1]=seg[k1][2];
						}
					}
				}
			}
		}
	}
}

void write_dbr_count(int fd)
{
	int i,j,k;
	FILE *fptr;
	
	if(dbr_mem) {
		if(!(fptr=fdopen(fd,"w"))) {
			perror("write_dbr_count(): can't fdopen()");
			return;
		}
		for(i=0;i<n_markers;i++) {
			if(dbr_count[0][i]) {
				print_marker_name(fptr,i);
				(void)fputc('\n',fptr);
				for(j=0;j<ped_size;j++) if(id_array[j].sire) {
					for(k=0;k<2;k++) if(dbr_count[k][i][j]) {
						print_orig_id(fptr,j+1);
						fprintf(fptr," : %c %d %d\n",k==X_MAT?'m':'p',dbr_count[k][i][j],dbr_countn[k][i][j]);
					}
				}
			}
		}
		fflush(fptr);
	}
}
