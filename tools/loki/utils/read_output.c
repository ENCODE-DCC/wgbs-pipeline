#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <float.h>

#include "utils.h"
#include "string_utils.h"
#include "lk_malloc.h"
#include "read_output.h"

#ifndef DBL_MAX
#define DBL_MAX MAXDOUBLE
#endif

int read_header(FILE *fptr,struct ro_data *dat,string **err_st)
{
	int i,j,k,err=0,lfg=0,fg,line=0,nlink=0,ncov=0;
	void *tbuf=0;
	string *s=0;
	char *p,*p1;
	char *match[]={"-----","Model: ","Linkage groups:","Total Map Length: ","Output format: ",
		  "No. genetic groups: ","Output columns:","QTL data blocks:","Number of QTL: ","No. QTL: ",0};
	struct ro_link *link=0,*link1;
	struct ro_mark *mark=0,*mark1;
	struct ro_cov *cov=0,*cov1;
	double v[2];
	
	dat->link=0;
	dat->model=0;
	dat->cov=0;
	dat->format=dat->sex_map=dat->n_cov_col=0;
	dat->n_genetic_groups=1;
	dat->total_map[0]=dat->total_map[1]=-DBL_MAX;
	dat->n_qtl[0]=dat->n_qtl[1]=-1;
	*err_st=0;
	do {
		line++;
		s=fget_string(fptr,s,&tbuf);
		if(!s->len) break;
		p=get_cstring(s);
		i=0;
		while(match[i] && strncmp(p,match[i],strlen(match[i]))) i++;
		if(match[i]) {
			if(!i) break;
			switch(i) {
			 case 1:
				p+=7;
				j=(int)strlen(p);
				dat->model=p1=lk_malloc(j);
				while(--j) *(p1++)=*(p++);
				*p1=0;
				break;
			 case 2:
				lfg=1;
				break;
			 case 3:
				p+=18;
				while(*p && isspace((int)*p)) p++;
				if(!*p) err=1;
				else {
					dat->total_map[0]=strtod(p,&p1);
					if(p==p1) err=1;
				}
				if(!err) {
					if(strncmp(p1,"cM",2)) err=1;
					else p1+=2;
				}
				if(!err) {
					while(*p1 && isspace((int)*p1)) p1++;
					if(*p1) {
						dat->total_map[1]=strtod(p1,&p);
						if(p==p1) err=1;
						else dat->sex_map=1;
						if(!err && strncmp(p,"cM",2)) err=1;
					}
				}
				lfg=0;
				break;
			 case 4:
				dat->format=atoi(p+15);
				break;
			 case 5:
				dat->n_genetic_groups=atoi(p+20);
				break;
			 case 6:
				lfg=2;
				break;
			 case 7:
				lfg=0;
				break;
			 case 8:
				dat->n_qtl[0]=strtol(p+15,&p1,10);
				if(p1==p+15 || strncmp(p1," to ",4)) err=1;
				else {
					dat->n_qtl[1]=strtol(p1+4,&p,10);
					if(p==p1+4) err=1;
				}
				break;
			 case 9:
				dat->n_qtl[0]=dat->n_qtl[1]=atoi(p+9);
				break;
			}
		} else if(lfg==1) {
			while(*p && isspace((int)*p)) p++;
			fg=0;
			if(isdigit((int)*p)) {
				j=strtol(p,&p1,10);
				if(*(p1++)==':') {
					fg=1;
					if(j!=++nlink) err=2;
					else {
						link1=lk_malloc(sizeof(struct ro_link));
						link1->next=0;
						link1->name=0;
						link1->markers=mark=0;
						for(j=0;j<2;j++) link1->map_r1[j]=link1->map_r2[j]=-DBL_MAX;
						if(link) link->next=link1;
						else dat->link=link1;
						link=link1;
						while(*p1 && isspace((int)*p1)) p1++;
						j=(int)strlen(p1)-1;
						while(j>0 && isspace((int)p1[j])) j--;
						if(j>0 && p1[j]==')') {
							if(p1[--j]!='M' || p1[--j]!='c') err=1;
							if(!err) {
								j--;
								while(j>0 && !isspace((int)p1[j])) j--;
								v[1]=strtod(p1+j,&p);
								if(p==p1 || *p!='c') err=1;
							}
							if(!err) {
								if(j>4) {
									j-=5;
									if(strncmp(p1+j,"cM to ",6)) err=1;
								}
							}
							if(!err) {
								while(j>0 && !isspace((int)p1[j])) j--;
								if(!j || p1[j+1]!='(') err=1;
								else {
									v[0]=strtod(p1+j+2,&p);
									if(p==p1 || *p!='c') err=1;
								}
							}
						} else err=1;
						if(!err) {
							while(j>0 && isspace((int)p1[j])) j--;
							if(j>8 && p1[j]==')') {
								link->map_r1[1]=v[0];
								link->map_r2[1]=v[1];
								dat->sex_map=1;
								if(p1[--j]!='M' || p1[--j]!='c') err=1;
								if(!err) {
									j--;
									while(j>0 && !isspace((int)p1[j])) j--;
									v[1]=strtod(p1+j,&p);
									if(p==p1 || *p!='c') err=1;
								}
								if(!err) {
									if(j>4) {
										j-=5;
										if(strncmp(p1+j,"cM to ",6)) err=1;
									}
								}
								if(!err) {
									while(j>0 && !isspace((int)p1[j])) j--;
									if(!j || p1[j+1]!='(') err=1;
									else {
										v[0]=strtod(p1+j+2,&p);
										if(p==p1 || *p!='c') err=1;
									}
								}
								while(j>0 && isspace((int)p1[j])) j--;
							}
							link->map_r1[0]=v[0];
							link->map_r2[0]=v[1];
							if(j>8 && p1[j]==':') {
								j-=9;
								if(strncmp(p1+j,"Map range:",10)) err=1;
								else {
									j--;
									while(j>0 && isspace((int)p1[j])) j--;
									link->name=lk_malloc(j+2);
									memcpy(link->name,p1,j+1);
									link->name[j+1]=0;
								}
							}
						}
					}
				}
			}
			if(!err && !fg && link) {
				mark1=lk_malloc(sizeof(struct ro_mark));
				mark1->next=0;
				mark1->name=0;
				mark1->pos[0]=mark1->pos[1]=-DBL_MAX;
				if(mark) mark->next=mark1;
				else link->markers=mark1;
				mark=mark1;
				j=(int)strlen(p)-1;
				for(k=dat->sex_map;!err && k>=0;k--) {
					while(j && isspace((int)p[j])) j--;
					while(j && !isspace((int)p[j])) j--;
					mark->pos[k]=strtod(p+j+1,&p1);
					if(p1==p+j+1 || !isspace((int)*p1)) err=1;
				}
				if(!err) {
					if(j>2 && p[j-1]=='-' && p[j-2]==' ') {
						j-=1;
						mark->name=lk_malloc((size_t)j);
						memcpy(mark->name,p,j-1);
						mark->name[j-1]=0;
					} else err=1;
				}
			}
		} else if(lfg==2) {
			while(*p && isspace((int)*p)) p++;
			if(isdigit((int)*p)) {
				j=strtol(p,&p1,10);
				if(*(p1++)==':') {
					if(j!=++ncov) err=3;
					else {
						cov1=lk_malloc(sizeof(struct ro_cov));
						cov1->next=0;
						cov1->name=0;
						if(cov) cov->next=cov1;
						else dat->cov=cov1;
						cov=cov1;
						while(*p1 && isspace((int)*p1)) p1++;
						j=strlen(p1);
						cov->name=lk_malloc((size_t)j);
						memcpy(cov->name,p1,j-1);
						cov->name[j-1]=0;
					}
				}
			}
		}
	} while(!err);
	if(!err) {
		if(dat->sex_map) {
			link=dat->link;
			while(link) {
				if(link->map_r1[1]==-DBL_MAX) err=4;
				else {
					mark=link->markers;
					while(mark) {
						if(mark->pos[1]==-DBL_MAX) err=4;
						mark=mark->next;
					}
				}
				link=link->next;
			}
			if(dat->total_map[1]==-DBL_MAX) err=4;
		}
	}
	switch(err) {
	 case 1:
		*err_st=strprintf(0,"Line %d: format error reading line:\n>> %s",line,get_cstring(s));
		break;
	 case 2:
		*err_st=strprintf(0,"Line %d: missing linkage group - expecting %d when reading line:\n>> %s",line,nlink,get_cstring(s));
		break;
	 case 3:
		*err_st=strprintf(0,"Line %d: missing effect - expecting %d when reading line:\n>> %s",line,ncov,get_cstring(s));
		break;
	 case 4:
		*err_st=strprintf(0,"Error - not all map positions have two values\n");
		break;
	}
	if(s) free_string(s);
	if(tbuf) free_fget_buffer(&tbuf);
	return err;
}
