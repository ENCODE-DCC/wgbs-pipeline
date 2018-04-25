#include <config.h>

#include <stdlib.h>
#include <string.h>
#ifdef USE_DMALLOC
#include <dmalloc.h>
#endif
#include <stdio.h>
#include <math.h>
#include <glob.h>
#include <sys/stat.h>

#include "config.h"
#include "lk_malloc.h"
#include "libhdr.h"
#include "snprintf.h"

struct bkfile {
	char *name;
	int num;
};

static int cmp_func(const void *s1,const void *s2) 
{
	int i,j,k;
	
	i=((struct bkfile *)s1)->num;
	j=((struct bkfile *)s2)->num;
	if(i<j) k=1;
	else if(i>j) k=-1;
	else k=0;
	return k;
}

int mkbackup(const char *fname,int nback) 
{
	int i=0,j,k,k1,k2;
	char *buf,*p,*p1;
	struct stat sbuf;
	struct bkfile *bkfiles;
	glob_t globbuf;
	
	if(!fname) {
		(void)fprintf(stderr,"mkbackup(): NULL filename\n");
		return -1;
	}
	if(nback<1) {
		(void)fprintf(stderr,"mkbackup(): invalid backup number\n");
		return -1;
	}
	if(nback==1) {
		if(!stat(fname,&sbuf)) {
			i=1;
			j=(int)strlen(fname);
			buf=(char *)lk_malloc((size_t)j+2);
			if(buf) {
				(void)strncpy(buf,fname,j);
				buf[j]='~';
				buf[j+1]='\0';
				i=rename(fname,buf);
				free(buf);
			}
			if(i) (void)fprintf(stderr,"mkbackup(): Couldn't rename old seedfile\n");
		}
	} else {
		j=(int)strlen(fname);
		buf=lk_malloc((size_t)j+2);
		i=1;
		if(buf) {
			(void)strncpy(buf,fname,j);
			buf[j]='*';
			buf[j+1]='\0';
			globbuf.gl_offs=2;
			(void)glob(buf,0,NULL,&globbuf);
			free(buf);
			if(globbuf.gl_pathc) {
				bkfiles=lk_malloc(sizeof(struct bkfile)*globbuf.gl_pathc);
				if(bkfiles) {
					for(k=k1=0;k<(int)globbuf.gl_pathc;k++) {
						p=globbuf.gl_pathv[k];
						if(p[j]=='~') {
							k2=(int)strtol(p+j+1,&p1,10);
							if(k2>0 && k2<nback && *p1++=='~' && *p1=='\0') {
								bkfiles[k1].name=p;
								bkfiles[k1++].num=k2;
							}
						} else if(p[j]=='\0') {
							bkfiles[k1].name=p;
							bkfiles[k1++].num=0;
						}
					}
					if(k1) {
						qsort(bkfiles,(size_t)k1,sizeof(struct bkfile),cmp_func);
						k2=1+(int)(log((double)nback)/log(10.0));
						buf=lk_malloc((size_t)(k2+3+j));
						if(buf) {
							for(k=0;k<k1;k++) {
								(void)snprintf(buf,k2+3+j,"%s~%d~",fname,bkfiles[k].num+1);
								i=rename(bkfiles[k].name,buf);
								if(i) break;
							}
							free(buf);
						}
					} else i=0;
					free(bkfiles);
				}
			} else i=0;
 			globfree(&globbuf);
		}
	}
	return i;
}

