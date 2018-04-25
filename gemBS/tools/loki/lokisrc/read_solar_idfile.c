/****************************************************************************
 *                                                                          *
 *     Loki - Programs for genetic analysis of complex traits using MCMC    *
 *                                                                          *
 *                   Simon Heath - CNG, Evry                                *
 *                                                                          *
 *                        October 2002                                      *
 *                                                                          *
 * read_solar_idfile.c:                                                     *
 *                                                                          *
 * Read in pedindex file from SOLAR with translations between external IDs  *
 * and SOLAR IBDIDs                                                         *
 *                                                                          *
 * Copyright (C) Simon C. Heath 2002                                        *
 * This is free software.  You can distribute it and/or modify it           *
 * under the terms of the Modified BSD license, see the file COPYING        *
 *                                                                          *
 ****************************************************************************/

#include <config.h>
#include <stdlib.h>
#include <string.h>
#ifdef USE_DMALLOC
#include <dmalloc.h>
#endif
#include <stdio.h>
#include <ctype.h>

#include "utils.h"
#include "loki.h"
#include "read_solar_idfile.h"

#define BUFSIZE 256

static struct bin_node **id_root,*fam_root;

static struct bin_node *new_node(union arg_type *p,int i)
{
	struct bin_node *p1;
	struct rs_data *rs;
	
	if(!(p1=malloc(sizeof(struct bin_node)))) ABT_FUNC(MMsg);
	p1->left=p1->right=0;
	p1->balance=0;
	if(!(rs=malloc(sizeof(struct rs_data)))) ABT_FUNC(MMsg);
	rs->data=p;
	rs->id=i;
	p1->data=rs;
	return p1;
}

static struct bin_node *add_node(struct bin_node *node,int flag,union arg_type *data,int id,int *bal)
{
	int k;
	struct rs_data *rs;
	
	rs=node->data;
	if(flag==ST_STRING) k=strcmp(rs->data->string,data->string);
	else k=rs->data->value-data->value;
	if(k<0) {
		if(node->left) node->left=add_node(node->left,flag,data,id,bal);
		else {
			node->left=new_node(data,id);
			*bal=0;
		}
		if(!(*bal)) {
			switch(node->balance) {
			 case -1:
				node=rotate_left(node);
				*bal=1;
				break;
			 case 0:
				node->balance=-1;
				break;
			 case 1:
				node->balance=0;
				*bal=1;
			}
		}
	} else if(k>0) {
		if(node->right) node->right=add_node(node->right,flag,data,id,bal);
		else {
			node->right=new_node(data,id);
			*bal=0;
		}
		if(!(*bal)) {
			switch(node->balance) {
			 case -1:
				node->balance=0;
				*bal=1;
				break;
			 case 0:
				node->balance=1;
				break;
			 case 1:
				node=rotate_right(node);
				*bal=1;
			}
		}
	} else *bal=1;
	return node;
}

static int find_node(struct bin_node *node,int flag,union arg_type *data)
{
	int i=-1,k;
	struct rs_data *rs;
	
	rs=node->data;
	if(flag==ST_STRING) k=strcmp(rs->data->string,data->string);
	else k=rs->data->value-data->value;
	if(k<0) {
		if(node->left) i=find_node(node->left,flag,data);
	} else if(k>0) {
		if(node->right) i=find_node(node->right,flag,data);
	} else i=rs->id;
	return i;
}

void read_solar_idfile(int *trans,struct loki *loki)
{
	FILE *fptr;
	int i,j,k,ibdid_start,ibdid_width,famid_start,famid_width,id_start,id_width,line;
	char buf[BUFSIZE],buf1[BUFSIZE],*fname,*p,*p1;
	union arg_type data;
	size_t l;
	struct Id_Record *id_array;
	
	fputs("Making translation table from Loki IDs to Solar IDs...",stdout);
	fflush(stdout);
	id_array=loki->pedigree->id_array;
	for(i=0;i<loki->pedigree->ped_size;i++) trans[i]=0;
	if(loki->pedigree->family_id) {
		fam_root=new_node(loki->pedigree->fam_recode.recode,0);
		for(i=1;i<loki->pedigree->n_orig_families;i++) {
			fam_root=add_node(fam_root,loki->pedigree->fam_recode.flag,loki->pedigree->fam_recode.recode+i,i,&j);
		}
		if(!(id_root=malloc(sizeof(void *)*loki->pedigree->n_orig_families))) ABT_FUNC(MMsg);
		for(i=0;i<loki->pedigree->n_orig_families;i++) id_root[i]=0;
		for(i=0;i<loki->pedigree->ped_size;i++) {
			j=id_array[i].fam_code-1;
			if(!id_root[j]) id_root[j]=new_node(loki->pedigree->id_recode.recode+i,i);
			else id_root[j]=add_node(id_root[j],loki->pedigree->id_recode.flag,loki->pedigree->id_recode.recode+i,i,&k);
		}
	} else {
		if(!(id_root=malloc(sizeof(void *)))) ABT_FUNC(MMsg);
		*id_root=new_node(loki->pedigree->id_recode.recode,0);
		for(i=1;i<loki->pedigree->ped_size;i++) {
			*id_root=add_node(*id_root,loki->pedigree->id_recode.flag,loki->pedigree->id_recode.recode+i,i,&k);
		} 
	}
	if(!(fptr=fopen(SOLAR_PEDINDEX_FILE,"r"))) abt(__FILE__,__LINE__,"%s(): Couldn't open file '%s' for input\n",__func__,SOLAR_PEDINDEX_FILE);
	if(!(fgets(buf,BUFSIZE,fptr))) ABT_FUNC("Error reading from file\n");
	qstrip(buf);
	l=strlen(buf);
	if(!l) ABT_FUNC("Null filename read in\n");
	if(!(fname=malloc(l+1))) ABT_FUNC(MMsg);
	memcpy(fname,buf,l+1);
	i=0;
	ibdid_start=famid_start=id_start=-1;
	ibdid_width=famid_width=id_width=0;
	while(fgets(buf,BUFSIZE,fptr)) {
		l=strlen(buf);
		if(l>=54) {
			j=atoi(buf);
			p1=p=buf+25;
			while(*p1 && !isspace((int)*p1)) p1++;
			*p1=0;
			if(!strcmp(p,"IBDID")) {
				ibdid_start=i;
				ibdid_width=j;
				if(buf[53]!='I') ABT_FUNC("IBDID field not integer\n");
				if(j>=BUFSIZE) ABT_FUNC("IBDID field width too large\n");
			} else if(!strcmp(p,"FAMID")) {
				famid_start=i;
				famid_width=j;
				if(j>=BUFSIZE) ABT_FUNC("FAMID field width too large\n");
			} else if(!strcmp(p,"ID")) {
				id_start=i;
				id_width=j;
				if(j>=BUFSIZE) ABT_FUNC("ID field width too large\n");
			}
			i+=j;
		}
	}
	if(ibdid_start<0 || id_start<0) ABT_FUNC("Error reading pedindex file\n");
	if(loki->pedigree->family_id) {
		if(famid_start<0) ABT_FUNC("No FAMID in pedindex file\n");
	} else if(famid_start>=0) ABT_FUNC("Unexpected FAMID found in pedindex file\n");
	fclose(fptr);
	if(ibdid_width+famid_width+id_width+1>=BUFSIZE) ABT_FUNC("BUFSIZE too small\n");
	if(!(fptr=fopen(fname,"r"))) abt(__FILE__,__LINE__,"%s(): Couldn't open file '%s' for input\n",__func__,fname);
	line=0;
	while(fgets(buf,BUFSIZE,fptr)) {
		line++;
		l=strlen(buf);
		if((l>=(size_t)(id_start+id_width))) {
			memcpy(buf1,buf+ibdid_start,ibdid_width);
			buf1[ibdid_width]=0;
			i=atoi(buf1);
			if(i<1 || i>loki->pedigree->ped_size) {
				fprintf(stderr,"Invalid IBDID read in at line %d of file %s\n",line,fname);
				ABT_FUNC("Aborting\n");
			}
			if(famid_start>=0) {
				memcpy(buf1,buf+famid_start,famid_width);
				buf1[famid_width]=0;
				qstrip(buf1);
				l=strlen(buf1);
				if(!l) {
					fprintf(stderr,"Null FAMID read in at line %d of file %s\n",line,fname);
					ABT_FUNC("Aborting\n");
				}
			}
			p=buf1;
			p1=buf1+famid_width+1;
			memcpy(p1,buf+id_start,id_width);
			p1[id_width]=0;
			qstrip(p1);
			l=strlen(buf1);
			if(!l) {
				fprintf(stderr,"Null ID read in at line %d of file %s\n",line,fname);
				ABT_FUNC("Aborting\n");
			}
			if(loki->pedigree->family_id) {
				if(loki->pedigree->fam_recode.flag==ST_STRING) data.string=p;
				else data.value=atoi(p);
				j=find_node(fam_root,loki->pedigree->fam_recode.flag,&data);
				if(j<0) {
					fprintf(stderr,"Family ID '%s' not found at line %d of file %s\n",p,line,fname);
					ABT_FUNC("Aborting\n");
				}
			} else j=0;
			if(loki->pedigree->id_recode.flag==ST_STRING) data.string=p1;
			else data.value=atoi(p1);
			k=find_node(id_root[j],loki->pedigree->id_recode.flag,&data);
			if(k<0) {
				fprintf(stderr,"ID '%s' not found at line %d of file %s\n",p1,line,fname);
				ABT_FUNC("Aborting\n");
			}
			trans[k]=i;
		}
	}
	fclose(fptr);
	free(fname);
	if(loki->pedigree->family_id) {
		free_bin_tree(fam_root,free);
		for(i=0;i<loki->pedigree->n_orig_families;i++) free_bin_tree(id_root[i],free);
	} else free_bin_tree(*id_root,free);
	free(id_root);
	fputs("Done\n",stdout);
}
