/****************************************************************************
 *                                                                          *
 *     Loki - Programs for genetic analysis of complex traits using MCMC    *
 *                                                                          *
 *                     Simon Heath - CNG, Evry                              *
 *                                                                          *
 *                         February 2003                                    *
 *                                                                          *
 * loki_compress.c:                                                         *
 *                                                                          *
 * Copyright (C) Simon C. Heath 2003                                        *
 * This is free software.  You can distribute it and/or modify it           *
 * under the terms of the Modified BSD license, see the file COPYING        *
 *                                                                          *
 ****************************************************************************/

#include <config.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#ifdef USE_DMALLOC
#include <dmalloc.h>
#endif
#include <stdio.h>
#include <sys/param.h>
#include <sys/stat.h>
#include <locale.h>
#include <assert.h>

#include "utils.h"
#include "libhdr.h"
#include "lk_malloc.h"
#include "loki_compress.h"
#include "snprintf.h"
#include "string_utils.h"

#define DEFAULT_PATH "/bin:/usr/bin:/usr/local/bin:/opt/local/bin";

static struct lk_compress *lkc;

char *find_prog(const char *prog,const char *path)
{
	char *p,*p1,*path1,*prog1,name[MAXPATHLEN];
	int sz,sz1,found,i;
	struct stat buf;
	tokens *tok;
	
	prog1=strdup(prog);
	found=0;
	tok=tokenize(prog1,':',0);
	for(i=0;!found && i<tok->n_tok;i++) {
		sz1=(int)strlen(tok->toks[i]);
		if(!(p1=path1=strdup(path))) ABT_FUNC(MMsg);
		while((p=my_strsep(&path1,":"))) {
			if(!*p) {
				p=".";
				sz=1;
			} else {
				sz=(int)strlen(p);
				while(p[sz-1]=='/') p[--sz]=0;
			}
			assert(sz+sz1+1<MAXPATHLEN);
			(void)snprintf(name,MAXPATHLEN,"%s/%s",p,tok->toks[i]);
			if(!stat(name,&buf) && S_ISREG(buf.st_mode) && !access(name,X_OK)) {
				found=1;
				break;
			}
		}
		(void)free(p1);
	}
	free(prog1);
	if(tok) free_tokens(tok);
	if(found) {
		if(!(p=strdup(name))) ABT_FUNC(MMsg);
		return p;
	}
	return 0;
}

static void free_compress(void)
{
  int i,j;
	
  if(lkc) {
    for(i=0;i<COMPRESS_NONE;i++) {
      free(lkc->compress_suffix[i]);
      for(j=0;j<2;j++)
	if(lkc->comp_path[i][j]) free(lkc->comp_path[i][j]);
    }
	  free(lkc);
	  lkc=0;
  }
}

struct lk_compress *init_compress(void)
{
  int i,j;
  char *pnames[COMPRESS_NONE][2]={{"gzip",0},{"bzip2",0},{"zip","funzip"},{"compress",0}};
  char *suff[]={"gz","bz2","zip","Z"};
  char *path;

	if(!lkc) {
		lkc=lk_malloc(sizeof(struct lk_compress));
		(void)setlocale(LC_ALL,"");
		if(!(path=getenv("PATH"))) path=DEFAULT_PATH;
		for(i=0;i<COMPRESS_NONE;i++) {
			lkc->compress_suffix[i]=strdup(suff[i]);
			for(j=0;j<2;j++) 
			  lkc->comp_path[i][j]=pnames[i][j]?find_prog(pnames[i][j],path):0;
		}
		if(!lkc->comp_path[COMPRESS_ZIP][1]) {
			if(lkc->comp_path[COMPRESS_ZIP][0]) {
				free(lkc->comp_path[COMPRESS_ZIP][0]);
				lkc->comp_path[COMPRESS_ZIP][0]=0;
			}
		}
		for(i=0;i<COMPRESS_NONE;i++) if(lkc->comp_path[i][0]) break;
		lkc->default_compress=i;
	
		if(atexit(free_compress)) ABT_FUNC("Unable to register exit function free_compress()\n");
	}
  return lkc;
}

