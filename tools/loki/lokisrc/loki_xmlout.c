/****************************************************************************
 *                                                                          *
 *     Loki - Programs for genetic analysis of complex traits using MCMC    *
 *                                                                          *
 *                        Simon Heath - CNG                                 *
 *                                                                          *
 *                           June 2004                                      *
 *                                                                          *
 * write_xml_dump:                                                          *
 *                                                                          *
 * Routines for write XML format dump file                                  *
 *                                                                          *
 * Copyright (C) Simon C. Heath 2004                                        *
 * This is free software.  You can distribute it and/or modify it           *
 * under the terms of the Modified BSD license, see the file COPYING        *
 *                                                                          *
 ****************************************************************************/

#include <config.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <sys/times.h>
#include <errno.h>
#include <sys/wait.h>
#include <signal.h>
#include <assert.h>

#include "bzlib.h"
#include "version.h"
#include "utils.h"
#include "string_utils.h"
#include "lk_malloc.h"
#include "loki.h"
#include "libhdr.h"
#include "loki_peel.h"
#include "loki_output.h"
#include "loki_utils.h"
#include "sample_rand.h"
#include "snprintf.h"
#include "loki_ibd.h"
#include "loki_tlmoves.h"
#include "loki_xmlout.h"

static bz_stream bz_str;
static char *out_buf=0,*in_buf=0;
static size_t buf_size=16384;
  
FILE *open_xmloutfile(struct loki *loki)
{
	int i=0,j;
	char *name;
	FILE *fout=0;
	
	name=make_file_name("_out.xml");
	j=loki->sys.syst_var[SYST_BACKUPS].flag?loki->sys.syst_var[SYST_BACKUPS].data.value:1;
	if(j) i=mkbackup(name,j);
	if(!i) {
		fout=fopen(name,"w");
		if(!fout || i) (void)fprintf(stderr,"[%s:%d] %s(): File Error.  Couldn't open '%s' for output\n",__FILE__,__LINE__,__func__,loki->names[LK_DUMPFILE]);
	}
	bz_str.bzalloc=0;
	bz_str.bzfree=0;
	bz_str.opaque=0;
	if(!i && BZ2_bzCompressInit(&bz_str,9,0,0)!=BZ_OK) i=1;
	if(!i) {
		out_buf=lk_malloc(buf_size);
		in_buf=lk_malloc(buf_size);
	}
	bz_str.next_out=out_buf;
	bz_str.next_in=in_buf;
	bz_str.avail_out=(unsigned int)buf_size;
	bz_str.avail_in=0;
	if(i && fout) {
		fclose(fout);
		fout=0;
	}
	return fout;
}

static int write_xmlout_buffer(FILE *fptr)
{
	int i=0;
	size_t s;
	
	s=buf_size-bz_str.avail_out;
	if(s && fwrite(out_buf,(size_t)1,s,fptr)!=s) i=1;
	bz_str.next_out=out_buf;
	bz_str.avail_out=(unsigned int)buf_size;
	return i;
}

void close_xmloutfile(FILE *fptr)
{
	int i;
	char *p;
	
	p=bz_str.next_in+bz_str.avail_in;
	assert(buf_size-bz_str.avail_in>=8);
	memcpy(p,"</loki>\n",(size_t)8);
	bz_str.avail_in+=8;
	do {
		i=BZ2_bzCompress(&bz_str,BZ_FINISH);
		if(buf_size>bz_str.avail_out) write_xmlout_buffer(fptr);
	} while(i!=BZ_STREAM_END);
	BZ2_bzCompressEnd(&bz_str);
	fclose(fptr);
}

int output_xmlheader(FILE *fptr,struct loki *loki)
{
	int i=0,j,s,n_comp,fam=0,pedsize,comp,cs,ng_flag,lnk;
	struct Id_Record *id_array;
	string *ss=0;
	char *p;
	
	p=bz_str.next_in+bz_str.avail_in;
	s=snprintf(p,buf_size-bz_str.avail_in,"<?xml version='1.0' encoding='UTF-8'?>\n<loki program='%s'>\n",LOKI_NAME);
	bz_str.avail_in+=s;
	p+=s;
	id_array=loki->pedigree->id_array;
	n_comp=loki->pedigree->n_comp;
	ng_flag=loki->pedigree->n_genetic_groups>1?1:0;
	if(loki->pedigree->family_id) {
		fam=id_array[0].fam_code;
		ss=str_print_orig_family(ss,1);
		s=snprintf(p,buf_size-bz_str.avail_in," <family orig='%s'>\n",get_cstring(ss));
	} else s=snprintf(p,buf_size-bz_str.avail_in," <family>\n");
	bz_str.avail_in+=s;
	for(i=comp=0;comp<n_comp;comp++) {
		p=bz_str.next_in+bz_str.avail_in;
		s=0;
		if(loki->pedigree->family_id) {
			j=id_array[i].fam_code;
			if(j!=fam) {
				fam=j;
				ss=str_print_orig_family(ss,i+1);
				s=snprintf(p,buf_size-bz_str.avail_in," </family>\n<family orig='%s'>\n",get_cstring(ss));
			}
		}
		s+=snprintf(p+s,buf_size-bz_str.avail_in-s,"  <comp>\n");
		bz_str.avail_in+=s;
		cs=loki->pedigree->comp_size[comp];
		for(j=0;j<cs;j++,i++) {
			p=bz_str.next_in+bz_str.avail_in;
			s=snprintf(p,buf_size-bz_str.avail_in,"   <ind id='%d'",i+1);
			bz_str.avail_in+=s;
			p+=s;
			if(id_array[i].sire) {
				s=snprintf(p,buf_size-bz_str.avail_in," father='%d'",id_array[i].sire);
				bz_str.avail_in+=s;
				p+=s;
			}
			if(id_array[i].dam) {
				s=snprintf(p,buf_size-bz_str.avail_in," mother='%d'",id_array[i].dam);
				bz_str.avail_in+=s;
				p+=s;
			}
			if(id_array[i].sex) {
				s=snprintf(p,buf_size-bz_str.avail_in," sex='%d'",id_array[i].sex);
				bz_str.avail_in+=s;
				p+=s;
			}
			if(id_array[i].affected) {
				s=snprintf(p,buf_size-bz_str.avail_in," aff='%d'",id_array[i].affected);
				bz_str.avail_in+=s;
				p+=s;
			}
			if(id_array[i].proband) {
				s=snprintf(p,buf_size-bz_str.avail_in," proband='%d'",id_array[i].proband);
				bz_str.avail_in+=s;
				p+=s;
			}
			if(ng_flag && id_array[i].group) {
				s=snprintf(p,buf_size-bz_str.avail_in," group='%d'",id_array[i].group);
				bz_str.avail_in+=s;
				p+=s;
			}
			if(loki->pedigree->id_recode.recode[i].string) {
				ss=str_print_orig_id(ss,i+1);
				s=snprintf(p,buf_size-bz_str.avail_in," orig='%s'",get_cstring(ss));
				bz_str.avail_in+=s;
				p+=s;
			}
			assert(buf_size-bz_str.avail_in>=2);
			*(p++)='>';
			*(p++)='\n';
			bz_str.avail_in+=2;
			BZ2_bzCompress(&bz_str,BZ_RUN);
			if(!bz_str.avail_in) bz_str.next_in=in_buf;
			if(buf_size>bz_str.avail_out) write_xmlout_buffer(fptr);
		}
		p=bz_str.next_in+bz_str.avail_in;
		s=snprintf(p,buf_size-bz_str.avail_in,"  </comp>\n");
		bz_str.avail_in+=s;
		p+=s;
	}
	s=snprintf(p,buf_size-bz_str.avail_in," </family>\n <seg_ids>\n");
	bz_str.avail_in+=s;
	pedsize=loki->pedigree->ped_size;
	for(i=j=0;i<pedsize;i++) if(id_array[i].sire) {
		p=bz_str.next_in+bz_str.avail_in;
		s=0;
		if(!j++) {
			assert(buf_size-bz_str.avail_in>=1);
			p[s++]=' ';
		}
		s+=snprintf(p+s,buf_size-bz_str.avail_in-s," %d",i+1);
		if(!(j%8)) {
			j=0;
			assert(buf_size-bz_str.avail_in>=1);
			p[s++]='\n';
		}
		bz_str.avail_in+=s;
		BZ2_bzCompress(&bz_str,BZ_RUN);
		if(!bz_str.avail_in) bz_str.next_in=in_buf;
		if(buf_size>bz_str.avail_out) write_xmlout_buffer(fptr);
	}
	s=0;
	p=bz_str.next_in+bz_str.avail_in;
	if(j) {
		assert(buf_size-bz_str.avail_in>=1);
		p[s++]='\n';
	}
	s+=snprintf(p+s,buf_size-bz_str.avail_in-s," </seg_ids>\n");
	bz_str.avail_in+=s;
	BZ2_bzCompress(&bz_str,BZ_RUN);
	if(!bz_str.avail_in) bz_str.next_in=in_buf;
	if(buf_size>bz_str.avail_out) write_xmlout_buffer(fptr);
	for(lnk=0;lnk<loki->markers->n_links;lnk++) {
		p=bz_str.next_in+bz_str.avail_in;
		s=snprintf(p,buf_size-bz_str.avail_in," <link");
		if(loki->markers->linkage[lnk].name) {
			s+=snprintf(p+s,buf_size-bz_str.avail_in-s," name='%s'",loki->markers->linkage[lnk].name);
		}
		switch(loki->markers->linkage[lnk].type) {
		 case LINK_X:
			s+=snprintf(p+s,buf_size-bz_str.avail_in-s," type='x'");
			break;
		 case LINK_Y:
			s+=snprintf(p+s,buf_size-bz_str.avail_in-s," type='y'");
			break;
		}
		s+=snprintf(p+s,buf_size-bz_str.avail_in-s,">\n");
		bz_str.avail_in+=s;
		BZ2_bzCompress(&bz_str,BZ_RUN);
		if(!bz_str.avail_in) bz_str.next_in=in_buf;
		if(buf_size>bz_str.avail_out) write_xmlout_buffer(fptr);
		for(i=0;i<loki->markers->linkage[lnk].n_markers;i++) {
			j=loki->markers->linkage[lnk].mk_index[i];
			p=bz_str.next_in+bz_str.avail_in;
			s=snprintf(p,buf_size-bz_str.avail_in,"  <locus name='%s'",loki->markers->marker[j].name);
			if(loki->markers->sex_map) s+=snprintf(p+s,buf_size-bz_str.avail_in-s," male='%g' female='%g'/>\n",loki->markers->marker[j].locus.pos[X_PAT],loki->markers->marker[j].locus.pos[X_MAT]);
			else s+=snprintf(p+s,buf_size-bz_str.avail_in-s," avg='%g'/>\n",loki->markers->marker[j].locus.pos[0]);
			bz_str.avail_in+=s;
			BZ2_bzCompress(&bz_str,BZ_RUN);
			if(!bz_str.avail_in) bz_str.next_in=in_buf;
			if(buf_size>bz_str.avail_out) write_xmlout_buffer(fptr);
		}
		p=bz_str.next_in+bz_str.avail_in;
		s=snprintf(p,buf_size-bz_str.avail_in-s," </link>\n");
		bz_str.avail_in+=s;
		BZ2_bzCompress(&bz_str,BZ_RUN);
		if(!bz_str.avail_in) bz_str.next_in=in_buf;
		if(buf_size>bz_str.avail_out) write_xmlout_buffer(fptr);
	}
	if(ss) free_string(ss);
	return 0;
}

int output_xmlout(FILE *fptr,int lp,int flag,struct loki *loki)
{
	int i=0,j,k,k1,k2,k3,n,lnk,s[2],ss,pedsize,skip;
	struct Locus **list;
	struct Id_Record *id_array;
	char *p;
	
	id_array=loki->pedigree->id_array;
	pedsize=loki->pedigree->ped_size;
	p=bz_str.next_in+bz_str.avail_in;
	if(flag) {
		k2=snprintf(p,buf_size-bz_str.avail_in," <iter val='%d' diff='1'>\n",lp);
	} else k2=snprintf(p,buf_size-bz_str.avail_in," <iter val='%d'>\n",lp);
	bz_str.avail_in+=k2;
	BZ2_bzCompress(&bz_str,BZ_RUN);
	if(!bz_str.avail_in) bz_str.next_in=in_buf;
	if(buf_size>bz_str.avail_out) write_xmlout_buffer(fptr);
	for(lnk=0;!i && lnk<loki->markers->n_links;lnk++) {
		list=get_sorted_locuslist(lnk,&n,1);
		if(!n) continue;
		p=bz_str.next_in+bz_str.avail_in;
		k2=snprintf(p,buf_size-bz_str.avail_in,"  <seg_rec");
		if(loki->markers->n_links>1 && loki->markers->linkage[lnk].name) {
			k2+=snprintf(p+k2,buf_size-bz_str.avail_in-k2," link='%s'",loki->markers->linkage[lnk].name);
		}
		k2+=snprintf(p+k2,buf_size-bz_str.avail_in-k2,">\n");
		bz_str.avail_in+=k2;
		BZ2_bzCompress(&bz_str,BZ_RUN);
		if(!bz_str.avail_in) bz_str.next_in=in_buf;
		skip=0;
		for(j=0;!i && j<pedsize;j++) if(id_array[j].sire) {
			if(flag) {
				for(k1=0;k1<2;k1++) {
					for(k=0;k<n;k++) if(list[k]->seg[k1][j]!=list[k]->seg_bk[k1][j]) break;
					if(k<n) break;
				}
				if(k1==2) {
					skip++;
					continue; /* Not changed from last time - don't print */
				}
			}
			if(skip) {
				p=bz_str.next_in+bz_str.avail_in;
				if(skip>15) k1=snprintf(p,buf_size-bz_str.avail_in,"   </skip n='%d'>\n",skip);
				else for(k1=0;k1<skip;k1++) p[k1]='\n';
				bz_str.avail_in+=k1;
				BZ2_bzCompress(&bz_str,BZ_RUN);
				if(!bz_str.avail_in) bz_str.next_in=in_buf;
				if(buf_size>bz_str.avail_out) write_xmlout_buffer(fptr);
				skip=0;
			}
			for(k1=0;k1<2;k1++) s[k1]=list[0]->seg[k1][j];
			p=bz_str.next_in+bz_str.avail_in;
			k2=snprintf(p,buf_size-bz_str.avail_in,"   %d",(s[0]<<1)|s[1]);
			k3=0;
			for(k=1;k<n;k++) {
				ss=list[k]->seg[0][j];
				if(ss!=s[0]) {
					p[k2++]=' ';
					k2+=snprintf(p+k2,buf_size-bz_str.avail_in-k2,"%d",k);
					s[0]=ss;
					k3=1;
				}
			}
			for(k=1;k<n;k++) {
				ss=list[k]->seg[1][j];
				if(ss!=s[1]) {
					if(k3!=1) p[k2++]=' ';
					if(k3!=2) {
						p[k2++]=',';
						k3=2;
					}
					k2+=snprintf(p+k2,buf_size-bz_str.avail_in-k2,"%d",k);
					s[1]=ss;
				}
			}
			p[k2++]='\n';
			bz_str.avail_in+=k2;
			BZ2_bzCompress(&bz_str,BZ_RUN);
			if(!bz_str.avail_in) bz_str.next_in=in_buf;
			if(buf_size>bz_str.avail_out) write_xmlout_buffer(fptr);
		}
		p=bz_str.next_in+bz_str.avail_in;
		assert(buf_size-bz_str.avail_in>=13);
		memcpy(p,"  </seg_rec>\n",(size_t)13);
		bz_str.avail_in+=13;
		BZ2_bzCompress(&bz_str,BZ_RUN);
		if(!bz_str.avail_in) bz_str.next_in=in_buf;
		if(buf_size>bz_str.avail_out) write_xmlout_buffer(fptr);
		for(k=0;k<n;k++) memcpy(list[k]->seg_bk[0],list[k]->seg[0],sizeof(int)*2*pedsize);
	}
	p=bz_str.next_in+bz_str.avail_in;
	assert(buf_size-bz_str.avail_in>=9);
	memcpy(p," </iter>\n",(size_t)9);
	bz_str.avail_in+=9;
	BZ2_bzCompress(&bz_str,BZ_RUN);
	if(!bz_str.avail_in) bz_str.next_in=in_buf;
	if(buf_size>bz_str.avail_out) write_xmlout_buffer(fptr);
	return i;
}
