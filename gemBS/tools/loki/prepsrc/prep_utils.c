/****************************************************************************
 *                                                                          *
 *     Loki - Programs for genetic analysis of complex traits using MCMC    *
 *                                                                          *
 *                   Simon Heath - CNG, Paris                               *
 *                                                                          *
 *                       August 2002                                        *
 *                                                                          *
 * prep_utils.c:                                                            *
 *                                                                          *
 * Copyright (C) Simon C. Heath 1997, 2000, 2002                            *
 * This is free software.  You can distribute it and/or modify it           *
 * under the terms of the Modified BSD license, see the file COPYING        *
 *                                                                          *
 ****************************************************************************/

#include <config.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>

#include "utils.h"
#include "libhdr.h"
#include "scan.h"
#include "y.tab.h"
#include "min_deg.h"
#include "prep_utils.h"

void print_orig_allele_id(FILE *fptr,const int i)
{
	if(i>0) {
		print_orig_id(fptr,i,0);
		(void)fputc('m',fptr);
	} else {
		print_orig_id(fptr,-i,0);
		(void)fputc('p',fptr);
	}
}

void New_DFE(const char *sfile,const int line,const char *file)
{
	(void)fprintf(stderr,"[%s:%d] Error writing to file '%s'\n",sfile,line,file);
	if(errno) perror("loki");
	exit(EXIT_FAILURE);
}

void add_to_list(const int i,int *j,int *list)
{
	int k;
	
	for(k=0;k<(*j);k++) if(list[k]==i) break;
	if(k==(*j)) list[(*j)++]=i;
}

int find_id_code(char *buf,int type,int fam)
{
	struct label_data *node=0;
	int i,flag;
	char *p;
	
	flag=fam<0?1:0;
	if(type==INTEGER) {
		i=(int)strtol(buf,&p,10);
		if(*p) (void)fprintf(stderr,"Garbage after id code '%s'\n",buf);
		else node=find_node(&i,type,flag);
	} else node=find_node(buf,type,flag);
	if(node) i=rec_tab[node->index];
	else {
		if(strcmp("*",buf)) {
			(void)fprintf(stderr,"Id code '%s' not found\n",buf);
			i=-1;
		} else i=0;
	}
	if(fam && i>=0) i=rec_tab1[(i-1)*n_orig_families+fam-1];
	else if(flag) i=-i;
	return i;
}

void cat_file(FILE *in,FILE *out,char *fname)
{
	char buf[1024];
	size_t l;
	
	if(fseek(in,0,SEEK_SET)<0) ABT_FUNC("Couldn't seek in input file\n");
	do {
		l=fread(buf,1,1024,in);
		if(fwrite(buf,1,l,out)!=l) DataFileError(fname);
	} while(l);
}

static void blank_ind(const int i,int *blank,const int locus,const int fam)
{
	int j,k;
	
	j=ped_recode1[i-1];
	if(j && blank[j-1]>=0) {
		blank[j-1]=fam+1;
		for(k=0;k<2;k++) id_array[i-1].haplo[k][locus]=0;
	}
}

void blank_fam(const int fam,int *blank,const int locus)
{
	int i,k;
	
	i=family[fam-1].sire;
	if(i) blank_ind(i,blank,locus,fam);
	i=family[fam-1].dam;
	if(i) blank_ind(i,blank,locus,fam);
	for(k=0;k<family[fam-1].nkids;k++) {
		i=family[fam-1].kids[k]+1;
		blank_ind(i,blank,locus,fam);
	}
}

