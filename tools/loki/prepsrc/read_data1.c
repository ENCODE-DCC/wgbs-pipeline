/****************************************************************************
 *                                                                          *
 *     Loki - Programs for genetic analysis of complex traits using MCMC    *
 *                                                                          *
 *             Simon Heath - University of Washington                       *
 *                                                                          *
 *                       March 1997                                         *
 *                                                                          *
 * read_data.c:                                                             *
 *                                                                          *
 * Reads in data from all input files and recodes factorial data.           *
 * Also calls setup_pedigree() <setup_ped.c> to initialize pedigree         *
 * data structures.                                                         *
 *                                                                          *
 * Copyright (C) Simon C. Heath 1997, 2000, 2002                            *
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
#include <math.h>
#include <stdio.h>
#include <ctype.h>
#include <errno.h>
#if HAVE_REGCOMP
#include <sys/types.h>
#include <regex.h>
#endif
#include <sys/wait.h>
#include <assert.h>

#include "config.h"
#include "utils.h"
#include "parser.h"
#include "parse.tab.h"
#include "string_utils.h"
#include "prep.h"
#include "prep_input.h"
#include "lk_malloc.h"
#include "libhdr.h"

#define INIT_BLOCK_SIZE 128 /* Start allocating memory in blocks of INIT_BLOCK_SIZE records, doubling */
#define MAX_BLOCK_SIZE 512  /* the size if more space required up to MAX_BLOCK_SIZE */

#define S_BLOCK_SIZE 512 /* For string allocation: should be at least as big as BUFFER_SIZE */
#define LINE_COUNT 5000 /* How often to print 'At line' */

static char *StringData=0;
static size_t StringPos;

static int get_scope(char *p)
{
	int i,j,scope;
	
	scope=j=0;
	while(*p) {
		i=toupper((int)*p);
		switch(i) {
		 case '!':
			j=1;
			break;
		 case 'P':
			scope|=j?~ST_PED:ST_PED;
			j=0;
			break;
		 case 'F':
			scope|=j?~ST_FACTOR:ST_FACTOR;
			j=0;
			break;
		 case 'G':
			scope|=j?~(ST_MARKER|ST_HAPLO):(ST_MARKER|ST_HAPLO);
			j=0;
			break;
		 case 'I':
			scope|=j?~ST_INTTYPE:ST_INTTYPE;
			break;
		 case 'C':
		 case 'R':
			scope|=j?~ST_REALTYPE:ST_REALTYPE;
			j=0;
			break;
		 default:
			ABT_FUNC("Illegal missing scope\n");
		}
		p++;
	}
	return scope;
}

struct miss_var_tag
{
	struct express **Missing;
	int nmiss;
};

static struct miss_var_tag *process_missing(struct prep_file *file,struct loki *loki)
{
	int i,k,scope=0,ncol,col,total_miss=0,match,ae_flag;
	struct miss_var_tag *miss_var;
	struct Missing *ms;
	struct parse_var_list *vl,*vl1;
	struct parse_term *v,*v1;
	struct parse_var *vv,*vv1;
	
	ncol=abs(file->ncol);
	assert(ncol);
	miss_var=lk_malloc(sizeof(struct miss_var_tag)*ncol);
	for(i=0;i<ncol;i++) miss_var->nmiss=0;
	ms=loki->prep->data->missing;
	while(ms) {
		if(ms->scope) scope=get_scope(ms->scope);
		col=0;
		vl=file->varlist;
		while(vl) {
			v=&vl->term;
			if((v->type&VT_TYPES)==VARIABLE) {
				k=1;
				if(v->type&VT_ARRAY_ELEMENT) {
					ae_flag=1;
					vv=get_array_var(v->elem.vv1->head,&v->elem.vv1->ex);
				} else {
					ae_flag=0;
					vv=v->elem.vv;
					if((vv->type&(VT_ARRAY|VT_SIZE_SET))==(VT_ARRAY|VT_SIZE_SET)) k=vv->size;
				}
				if((vl1=ms->vl)) {
					while(vl1) {
						v1=&vl1->term;
						if((v1->type&VT_TYPES)==VARIABLE) {
							match=0;
							if(ae_flag) {
								if(v1->type&VT_ARRAY_ELEMENT) {
									vv1=get_array_var(v1->elem.vv1->head,&v1->elem.vv1->ex);
									if(vv1==vv) match=1;
								} else {
									vv1=v1->elem.vv;
									if((vv1->type&VT_ARRAY) && (v1->elem.vv1->head==vv1)) match=1;
								}
							} else if(!(v1->type&VT_ARRAY_ELEMENT)) {
								if(v1->elem.vv==vv) match=1;
							}
							if(match) {
								for(i=0;i<k;i++) miss_var[col+i].nmiss++;
								total_miss+=k;
							}
						}
						vl1=vl1->next;
					}
				} else if(ms->scope) {
					if(vv->type&scope) {
 						for(i=0;i<k;i++) miss_var[col+i].nmiss++;
						total_miss+=k;
					}
				} else {
					for(i=0;i<k;i++) miss_var[col+i].nmiss++;
					total_miss+=k;
				}
				col+=k;
			}
			vl=vl->next;
		}
		ms=ms->next;
	}
	return miss_var;
}

int ReadData1(char *lfile,struct loki *loki)
{
	int i,j,skip,fs_reg=0,error_n=0,popen_flag,rec;
	struct prep_file *infile;
	char *rs,*fs,*gs,*p;
	FILE *fptr,*flog;
	string *buf=0;
	tokens *tok=0;
	void *tbuf=0;
#ifdef HAVE_REGCOMP	
	regex_t preg;
#endif

	infile=loki->prep->data->infiles;
	if(!(StringData=malloc(S_BLOCK_SIZE))) ABT_FUNC(MMsg);
	loki->sys.RemBlock=AddRemem(StringData,loki->sys.RemBlock);
	StringPos=0;
	if(lfile && (p=add_file_dir(lfile))) {
		flog=fopen(p,"a");
		free(p);
	} else flog=0;
	if(flog) i=fputs("\n******************** Reading in data ********************\n\n",flog);
	while(infile) {
		/* Handle field/record separators */
		fs_reg=0;
		fs=infile->fs?get_cstring(infile->fs):0;
		rs=infile->rs?get_cstring(infile->rs):0;
		gs=infile->gs?get_cstring(infile->gs):0;
		skip=infile->skip;
		if(fs) {
			if(strlen(fs)>1) {
#if HAVE_REGCOMP
				if((i=regcomp(&preg,fs,REG_EXTENDED))) fs=0;
				else fs_reg=1;
#endif				  
			} else if(fs[0]==' ') fs=0;
		}
		if(rs && *rs==0) rs=0;
		/* Set up missing data handlers */
		if(loki->prep->data->missing) process_missing(infile,loki); 
		if(!error_n) {
			/* Check for '|' as last (non-space) character in file name */
			p=infile->name;
			i=(int)strlen(p);
			while(i && isspace((int)p[i--]));
			assert(i); /* Filename just whitespace? */
			if(p[i]=='|') {
				/* Use popen on this */
				p[i]=0;
				if(!(fptr=popen(p,"r"))) {
					(void)fprintf(stderr,"read_data(): Can't execute '%s': %s\n",infile->name,strerror(errno));
					error_n++;
					break;
				}
				message(INFO_MSG,"Reading output from shell command '%s'\n",p);
			} else {
				/* Open file (checking for compressed files) */
				fptr=open_readfile_and_check(p,&popen_flag,loki);
				message(INFO_MSG,"Reading data from file '%s'\n",p);
			}
			rec=0;
			do {
				if(!rs) buf=fget_string(fptr,buf,&tbuf);
				else buf=fget_string_gen(fptr,buf,&tbuf,rs[0]);
				rec++;
				if(infile->fixed_start) {
					/* Fixed format read */
					p=get_cstring(buf);
					for(i=0;i<infile->ncol;i++) {
						fputc('{',stdout);
						for(j=infile->fixed_start[i];j<infile->fixed_end[i];j++) {
							if(j>=buf->len) break;
							fputc(p[j],stdout);
						}
						fputs("} ",stdout);
					}
					fputc('\n',stdout);
				} else {
					/* Free format read */
					/* Tokenize line */
					if(fs_reg) tok=reg_tokenize(get_cstring(buf),&preg,tok);
					else tok=tokenize(get_cstring(buf),fs?fs[0]:0,tok);
					if(tok) {
						for(i=0;i<tok->n_tok;i++) {
							printf("<%s> ",tok->toks[i]);
						}
						fputc('\n',stdout);
					}
				}
			} while(buf->len);
			fclose(fptr);
		}
		infile=infile->next;
	}
	if(tbuf) free_fget_buffer(&tbuf);
	if(buf) free_string(buf);
	if(tok) free_tokens(tok);
	return 0;
}
