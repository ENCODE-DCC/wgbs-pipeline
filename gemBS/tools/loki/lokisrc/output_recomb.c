/****************************************************************************
 *                                                                          *
 *     Loki - Programs for genetic analysis of complex traits using MCMC    *
 *                                                                          *
 *                        Simon Heath - MSKCC                               *
 *                                                                          *
 *                           January 2002                                   *
 *                                                                          *
 * output_recomb.c:                                                         *
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
#include "loki.h"
#include "loki_peel.h"
#include "output_recomb.h"

int output_recomb;
char *Recombfile;

/* Print recombinations forced by the data */
void print_recomb(void)
{
	int id,i,j,link,k,k1,s;
	int *perm,sg[2],flag[2];
	size_t max,sz;
	FILE *fptr;
	
	if(!n_markers) return;
	if(!Recombfile) ABT_FUNC("Internal error - null file name\n");
	if(!(fptr=fopen(Recombfile,"w"))) {
		fputs("Couldn't open file ",stderr);
		perror(Recombfile);
		return;
	}
	if(!(perm=malloc(sizeof(int)*n_markers))) ABT_FUNC(MMsg);
	max=get_max_idlen();
	for(link=0;link<n_links;link++) {
		get_locuslist(perm,link,&k,0);
		if(k) {
			for(i=0;i<k;i++) {
				j=perm[i];
				if(i) fputc(',',fptr);
				fprintf(fptr,"%s",marker[j].name);
				if(marker[j].index) fprintf(fptr,"(%d)",marker[j].index);
			}
			fputc('\n',fptr);
			for(id=0;id<ped_size;id++) if(id_array[id].sire) {
				sg[0]=sg[1]=-1;
				flag[0]=flag[1]=0;
				for(i=0;i<k;i++) {
					j=perm[i];
					for(k1=0;k1<2;k1++) {
						s=marker[j].locus.seg[k1][id];
						if(s==0 || s==1) {
							if(sg[k1]>=0) {
								if(s!=sg[k1]) {
									flag[k1]=1;
								}
							} else sg[k1]=s;
						}
					}
				}
				for(k1=0;k1<2;k1++) if(flag[k1]) {
					sz=print_orig_id(fptr,id+1);
					for(;sz<max;sz++) fputc(' ',fptr);
					fputs(k1?" m: ":" p: ",fptr);
					for(i=0;i<k;i++) {
						j=perm[i];
						s=marker[j].locus.seg[k1][id];
						if(s==0 || s==1) fputc(s?'1':'0',fptr);
						else fputc('.',fptr);
					}
					fputc('\n',fptr);
				}
			}
		}
	}
	fclose(fptr);
	free(perm);
}
