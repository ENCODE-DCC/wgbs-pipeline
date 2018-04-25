/****************************************************************************
 *                                                                          *
 *     Loki - Programs for genetic analysis of complex traits using MCMC    *
 *                                                                          *
 *             Simon Heath - University of Washington                       *
 *                                                                          *
 *                       March 1997                                         *
 *                                                                          *
 * utils.c:                                                                 *
 *                                                                          *
 * Small utility routines used by both prep and loki                        *
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
#include <math.h>
#ifdef HAVE_IEEEFP_H
#include <ieeefp.h>
#endif
#include <stdarg.h>
#include <stdio.h>
#include <ctype.h>
#include <signal.h>
#include <time.h>
#include <sys/times.h>
#include <sys/stat.h>
#include <errno.h>
#include <float.h>
#include <limits.h>
#ifdef HAVE_SYS_SYSTEMINFO_H
#include <sys/systeminfo.h>
#endif
#include <assert.h>

#include "ranlib.h"
#include "lk_malloc.h"
#include "utils.h"

#ifndef MAXHOSTNAMELEN
#define MAXHOSTNAMELEN 64
#endif

static char b64digits[64]="ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
static unsigned char b64conv[256];

const char *FMsg="File Error - aborting\n";
const char *IntErr="Fatal Internal Error - aborting\n";
const char *MMsg="Out of Memory Error - aborting\n";
const char *AbMsg="Aborting\n";
int from_abt;
static char *file_prefix,*file_dir;

static char *lfile;
static loki_time *ltime;
static int set_exit;

#if !HAVE_ISINF
#if HAVE_FPCLASS
int isinf(double x)
{
  fpclass_t c;
	
  c=fpclass(x);
  return (c==FP_NINF || c==FP_PINF)?1:0;
}
#elif HAVE_FINITE
int isinf (double x)
{
  return !finite(x) && x==x; 
}
#endif
#endif

void init_lib_utils(void)
{
  int i;
	
  set_strcmpfunc(strcasecmp);
  memset(b64conv,0xff,(size_t)256);
  for(i=0;i<64;i++) b64conv[(int)b64digits[i]]=i;
  return;
}

#if 0
#define SEG_LOG2 0.69314718055994530942

double addlog(double x1,double x2)
{
  double yy, zz, r;

  yy=.5*(x1+x2);
  zz=cosh(x1-yy);
  if(isinf(zz)) {
    r=x1>x2?x1:x2;
  } else r=(log(zz)+yy+SEG_LOG2);
  return r;
}
#endif

/* Returns log(exp(x1)+exp(x2)) */
double addlog(const double x1,const double x2)
{
  if(x1>x2) {
    register const double diff=x2-x1;
    return (diff<-745 ? x1 : x1+log1p(exp(diff)));
  } else {
    register const double diff=x1-x2;
    return (diff<-745 ? x2 : x2+log1p(exp(diff)));
  }
}

/* Returns log(exp(x1)-exp(x2)); x1>x2 */
double sublog(const double x1,const double x2)
{
  if(x1>x2) {
    register const double diff=x2-x1;
    return (diff<-745 ? x1 : x1+log1p(-exp(diff)));
  } else {
    return NAN;
  }
}

int proc_verbose_level(const int level,const int op)
{
	static int verbose_level;
	
	if(op) verbose_level=level;
	return verbose_level;
}

int check_output(const int level) 
{
  int i;
  static int mask[4]={11,15,10,8};
	
  i=get_verbose_level();
  return (mask[i]&level);
}

void message(const int level, const char *fmt, ...)
{
  va_list args;
	
  if(check_output(level)) {
    va_start(args,fmt);
	 if(level==ERROR_MSG) (void)vfprintf(stderr,fmt,args);
    else (void)vfprintf(stdout,fmt,args);
    va_end(args);
  }
}

int count_bits(unsigned long x)
{
  int n=0;
	
  if(x) do n++; while(x&=(x-1));
  return n;
}

static void free_utils(void) 
{
  if(file_prefix) {
    free(file_prefix);
    file_prefix=0;
  }
  if(file_dir) {
    free(file_dir);
    file_dir=0;
  }
}

char *add_file_dir(const char *p)
{
  char *p1=0,*s2;
  const char *s1;
  size_t s;
	
  if(p) {
    if(file_dir) {
      s=strlen(p)+2+strlen(file_dir);
      p1=lk_malloc(s);
      if(p1) {
	s1=file_dir;
	s2=p1-1;
	while((*++s2=*s1++));
	if(s2[-1]!='/') *s2='/';
	else s2--;
	s1=p;
	while((*++s2=*s1++));
      }
    } else p1=strdup(p);
  }
  return p1;
}

int set_file_prefix(const char *p)
{
  int err=0;
  char *p1;
	
  if(*p) {
    p1=strdup(p);
    if(p1) {
      if(file_prefix) free(file_prefix);
      file_prefix=p1;
      if(!set_exit) {
	if(!atexit(free_utils)) set_exit=1;
      }
    } else err=2;
  } else err=1;
  return err;
}
			
int set_file_dir(const char *p) 
{
  int err=0,s;
  char *p1;
  struct stat sbuf;
	
  if(*p) {
    p1=strdup(p);
    if(p1) {
      s=stat(p,&sbuf);
      if(!s) {
	if(sbuf.st_mode&S_IFDIR) {
	  if((sbuf.st_mode&(S_IXUSR|S_IWUSR))==(S_IXUSR|S_IWUSR)) {
	    if(file_dir) free(file_dir);
	    file_dir=p1;
	    if(!set_exit) {
	      if(!atexit(free_utils)) set_exit=1;
	    }
	  } else err=5;
	} else err=4;
      } else err=3;
    } else err=2;
  } else err=1;
  return err;
}

char *utl_error(int i)
{
  static char *errs[]={
    "No Error",
    "Out of memory",
    "Null Pointer",
    "Couldn't stat() directory",
    "Not a directory",
    "Insufficient permissions",
    "Bad error number"
  };
	
  if(i>UTL_MAX_ERR) i=UTL_MAX_ERR+1;
  return errs[i];
}

void abt(const char *file, const int line, const char *fmt, ...)
{
  va_list args;

  if(stdout) (void)fflush(stdout);
  (void)fprintf(stderr,"[%s:%d] ",file,line);
  va_start(args,fmt);
  (void)vfprintf(stderr,fmt,args);
  va_end(args);
  from_abt=1; /* Avoid certain cleanup routines if aborting */
  exit(EXIT_FAILURE); 	
}

struct _lk_list {
  struct _lk_list *next;
};

void *reverse_list(void *s) 
{
  struct _lk_list *p,*p1,*p2;

  p=s;
  p2=0;
  while(p) {
    p1=p->next;
    p->next=p2;
    p2=p;
    p=p1;
  }
  return p2;
}

void free_list(void *s,void (*f)(void *))
{
  struct _lk_list *p,*p1;
	
  p=s;
  while(p) {
    p1=p->next;
    if(f) f(p);
    free(p);
    p=p1;
  }
}

int list_length(void *s)
{
  int c=0;
  struct _lk_list *p;
	
  p=s;
  while(p) {
    c++;
    p=p->next;
  }
  return c;
}

static const struct {
	char *errmsg;
	int err;
} errs[]= {
	{NULL,0},
	{"Invalid limits",ERANGE},
	{"Out of range",ERANGE},
	{"Invalid characters",EINVAL},
	{"Index out of range",ERANGE},
	{"Blank field",EINVAL},
};

double strtodnum(const char *s,double low,double high,char **errstr)
{
	int er=0,old_errno;
	double x=0.0;
	char *s1;
	
	assert(s);	
	old_errno=errno;
	errno=0;
	if(low>high) er=1;
	else {
		x=strtod(s,&s1);
		if(errno==EINVAL) er=3;
		else if(errno==ERANGE) er=2;
		else if(*s1 || s==s1) er=3;
		else if(x<low || x>high) er=2;
	}
	if(errstr) *errstr=errs[er].errmsg;
	if(er) {
		errno=errs[er].err;
		x=0.0;
	} else errno=old_errno;
	return x;
}

long strtolnum(const char *s,long low,long high,char **errstr)
{
	int er=0,old_errno;
	long x=0;
	char *s1;
	
	assert(s);	
	old_errno=errno;
	errno=0;
	if(low>high) er=1;
	else {
		x=strtol(s,&s1,10);
		if(errno==EINVAL) er=3;
		else if(errno==ERANGE) er=2;
		else if(*s1 || s==s1) er=3;
		else if(x<low || x>high) er=2;
	}
	if(errstr) *errstr=errs[er].errmsg;
	if(er) {
		errno=errs[er].err;
		x=0;
	} else errno=old_errno;
	return x;
}

char *getcolumn(const tokens *tok,const int i,const char *missing,char **errstr)
{
	char *s=0;
	int er=0,old_errno;
	
	assert(tok);
	old_errno=errno;
	errno=0;
	if(i>=0 && tok->n_tok>i) {
		s=tok->toks[i];
		if(!*s) s=0;
		else if(missing && !strcmp(missing,s)) s=0;
	} else er=4;
	if(errstr) *errstr=errs[er].errmsg;
	if(er) errno=errs[er].err;
	else errno=old_errno;
	return s; 
}

char *getstrcolumn(const tokens *tok,const int i,const char *missing,char **errstr)
{
	char *p;
	
	p=getcolumn(tok,i,missing,errstr);
	if(p) p=strdup(p);
	return p;
}

double getdnumcolumn(const tokens *tok,const int i,const char *missing,char **errstr)
{
	double x=0.0;
	char *p,*err;
	
	p=getcolumn(tok,i,missing,&err);
	if(!err) {
		if(p) x=strtodnum(p,DBL_MIN,DBL_MAX,&err);
		else {
			err=errs[5].errmsg;
			errno=errs[5].err;
		}
	}
	if(errstr) *errstr=err;
	return x;
}

long getlnumcolumn(const tokens *tok,const int i,const char *missing,char **errstr)
{
	long x=0;
	char *p,*err;
	
	p=getcolumn(tok,i,missing,&err);
	if(!err) {
		if(p) x=strtolnum(p,LONG_MIN,LONG_MAX,&err);
		else {
			err=errs[5].errmsg;
			errno=errs[5].err;
		}
	}
	if(errstr) *errstr=err;
	return x;
}

int mystrcmp(const char *p1, const char *p2)
{
  if(!p1) {
    if(!p2) return 0;
    return 1;
  }
  if(!p2) return 1;
  return strcmp(p1,p2);
}

void qstrip(char *s1)
{
  char *p,*p1;
	
  p=s1;
  p1=s1-1;
  while(*s1) {
    if(!isspace((int)*s1)) break;
    s1++;
  }
  while(*s1) {
    if(!isspace((int)*s1)) p1=p;
    *(p++)= *(s1++);
  }
  *(++p1)='\0';
}

tokens *copy_ntokens(tokens *tok,int n)
{
  int i;
  tokens *tok1=0;
  char *p,*p1;
  size_t sz;
	
  if(tok && tok->n_tok) {
    tok1=malloc(sizeof(tokens));
    if(tok1) {
      if(!n || n>tok->n_tok) n=tok->n_tok;
      tok1->size=tok1->n_tok=n;
      if(!(tok1->toks=malloc(sizeof(void *)*tok1->size))) {
	free(tok1);
	tok1=0;
      } else {
	p=tok->toks[n-1];
	sz=(p-tok->toks[0])+strlen(p)+1;
	if(!(p1=malloc(sz))) {
	  free_tokens(tok1);
	  tok1=0;
	} else {
	  for(i=0;i<n;i++) {
	    p=tok->toks[i];
	    tok1->toks[i]=p1;
	    while(*p) *p1++=*p++;
	    *p1++=*p;
	  }
	}
      }
    }
  }
  return tok1;
}

tokens *tokenize(char *s,const int ch,tokens *tok)
{
  int n_toks=0;
  char **p=0,*p1;
	
  if(!tok) {
    tok=lk_malloc(sizeof(tokens));
    tok->size=16;
    if(tok) {
      if(!(tok->toks=lk_malloc(sizeof(void *)*tok->size))) {
	free(tok);
	tok=0;
      }
    }
  }
  if(tok) {
    p=tok->toks;
    if((p1=s)) {
      if(!ch) { /* Split on white space */
	for(;;) {
	  while(*s && isspace((int)*s)) s++;
	  if(!*s) break;
	  if(n_toks==tok->size) {
	    tok->size<<=1;
	    if(!(p=lk_realloc(p,sizeof(void *)*tok->size))) {
	      free_tokens(tok);
	      tok=0;
	      break;
	    }
	    tok->toks=p;
	  }
	  p[n_toks++]=p1;
	  while(*s && !isspace((int)*s)) {
	    *p1++=*s++;
	  }
	  if(*s) s++;
	  *p1++=0;
	}
      } else { /* Split on token */
	for(;;) {
	  if(!*s) break;
	  if(n_toks==tok->size) {
	    tok->size<<=1;
	    if(!(p=lk_realloc(p,sizeof(void *)*tok->size))) {
	      free_tokens(tok);
	      tok=0;
	      break;
	    }
	    tok->toks=p;
	  }
	  p[n_toks++]=p1;
	  while(*s && *s!=ch) {
 	    *p1++=*s++;
	  }
	  if(*s) s++;
	  *p1++=0;
	  qstrip(p[n_toks-1]);
	}
      }
    }
  }
  if(tok) {
    if(n_toks==1 && !*p[0]) n_toks--;
    tok->n_tok=n_toks;
  }
  return tok;
}

#ifdef HAVE_REGCOMP	
tokens *reg_tokenize(char *s,regex_t *reg,tokens *tok)
{
  int n_toks=0,res;
  char **p=0;
  regmatch_t match;
	
  if(!tok) {
    tok=lk_malloc(sizeof(tokens));
    tok->size=16;
    if(tok) {
      if(!(tok->toks=lk_malloc(sizeof(void *)*tok->size))) {
	free(tok);
	tok=0;
      }
    }
  }
  if(tok) {
    p=tok->toks;
    if(s) {
      do {
	if(n_toks==tok->size) {
	  tok->size<<=1;
	  if(!(p=lk_realloc(p,sizeof(void *)*tok->size))) {
	    free_tokens(tok);
	    tok=0;
	    break;
	  }
	  tok->toks=p;
	}
	p[n_toks++]=s;
	if(!(res=regexec(reg,s,1,&match,0))) {
	  s[match.rm_so]=0;
	  s+=match.rm_eo;
	}
	qstrip(p[n_toks-1]);
      } while(!res);
    }
  }
  if(tok) {
    if(n_toks==1 && !*p[0]) n_toks--;
    tok->n_tok=n_toks;
  }
  return tok;
}
#endif

char *make_file_name(const char *s)
{
  size_t l;
  char *s1,*s2,*s3;
	
  assert(s);
  if(!file_prefix) {
    if(set_file_prefix(DEFAULT_FILE_PREFIX)) return 0;
  }
  l=strlen(s)+strlen(file_prefix)+2;
  if(file_dir) l+=strlen(file_dir)+1;
  assert(l);
  if(!(s1=lk_malloc(l))) return 0;
  s2=s1-1;
  s3=file_dir;
  if(s3) {
    while((*++s2 = *s3++));
    if(s2[-1]!='/') *s2='/';
    else s2--;
  }
  s3=file_prefix;
  while((*++s2 = *s3++));
  while((*s2++ = *s++));
  return s1;
}

void print_start_time(const char *progname,const char *mode,char *logfile,loki_time *lt)
{
  FILE *flog=0;
  char *buf;
  struct stat sbuf;
  int i,j;
	
  ltime=lt;
  if(logfile) logfile=add_file_dir(logfile);
  if(logfile && mode[0]=='w') {
    if(!stat(logfile,&sbuf)) {
      i=1;
      j=(int)strlen(logfile);
      buf=lk_malloc((size_t)j+2);
      if(buf) {
	(void)strncpy(buf,logfile,(size_t)j);
	buf[j]='~';
	buf[j+1]='\0';
	i=rename(logfile,buf);
	free(buf);
      }
      if(i) (void)fprintf(stderr,"print_start_time(): Couldn't rename old logfile\n");
    }
  }
  if(logfile) flog=fopen(logfile,mode);
  else flog=stdout;
  if(flog) {
    (void)fprintf(flog,"\n********************** Starting *************************\n\n");
    (void)fprintf(flog,"     %s: %s",progname,ctime(&lt->start_time));
    if(logfile) {
      (void)fclose(flog);
      lfile=logfile;
    }
  } else {
    (void)fprintf(stderr,"Couldn't write to log file %s\n",logfile?logfile:"<NULL>");
    if(logfile) free(logfile);
  }
}

void print_end_time(void)
{
  FILE *flog;
  time_t end_time;
  double l;
  struct tms tbuf;
  long tps;
  int t=0,flag=0;
  char hname[MAXHOSTNAMELEN+1];
	
  if(from_abt || !ltime) return;
  if(lfile) flog=fopen(lfile,"a");
  else flog=stdout;
  if(flog)	{
#ifdef HAVE_SYS_SYSTEMINFO_H
    if(sysinfo(SI_HOSTNAME,hname,MAXHOSTNAMELEN)<0)
#else
      if(gethostname(hname,MAXHOSTNAMELEN)<0)
#endif
	(void)strncpy(hname,"UNKNOWN",MAXHOSTNAMELEN+1);
		
    (void)fprintf(flog,"\n*********************** Exiting *************************\n\n");
    (void)fprintf(flog,"     Hostname:     %s\n", hname);
    tps=sysconf (_SC_CLK_TCK);
    errno=0;
    (void)times(&tbuf);
    if(errno) perror("print_end_time():");
    else {
      (void)fprintf (flog,"     System time:  %.4f seconds\n",ltime->extra_stime+(double)tbuf.tms_stime/(double)tps);
      (void)fprintf (flog,"     User time:    %.4f seconds\n",ltime->extra_utime+(double)tbuf.tms_utime/(double)tps);
    }
    end_time=time(0);
    l=ltime->extra_time+difftime(end_time,ltime->start_time);    
    (void)fputs("     Elapsed time: ",flog);
    if(l>86400.0) {
      t=(int)(l/86400.0);
      l-=(double)t*86400.0;
      (void)fprintf(flog,"%d day%s",t,t!=1?"s, ":", ");
      flag=1;
    }
    if(l>3600.0) {
      t=(int)(l/3600.0);
      l-=(double)t*3600.0;
      flag=1;
    }
    if(flag) (void)fprintf(flog,"%d hour%s",t,t!=1?"s, ":", ");
    if(l>60.0) {
      t=(int)(l/60.0);
      l-=(double)t*60.0;
      flag=1;
    }
    if(flag) (void)fprintf(flog,"%d minute%s",t,t!=1?"s, ":", ");
    (void)fprintf(flog,"%.4f seconds\n",l);
    if(lfile) (void)fclose(flog);
  }
  if(lfile) {
    free(lfile);
    lfile=0;
  }
}

void gen_perm(int *x,int n)
{
  int i,j;
	
  while(n) {
    j=(int)(safe_ranf()*(double)(n--));
    i=x[j];
    x[j]=x[n];
    x[n]=i;
  }
}

int txt_print_double(double x,FILE *fptr)
{
  int i,er=0;
  double y,z;

  y=frexp(x,&i);
  if(fprintf(fptr,"%d:",i)<0) er=1;
  if(!er && y<0.0) {
    y=-y;
    if(fputc('-',fptr)==EOF) er=1;
  }
  if(!er) {
    if(y) {
      while(y>0.0 && !er) {
	y=frexp(y,&i);
	y=ldexp(y,i+6);
	y=modf(y,&z);
	if(fputc(b64digits[(int)z],fptr)==EOF) er=1;
      }
    } else if(fputc('0',fptr)==EOF) er=1;
  }
  return er;
}

int txt_get_double(char *p,char **p1,double *x)
{
  int e,j,mf=0;
  char *p2;
  double y;
	
  e=(int)strtol(p,p1,10);
  p2=*p1;
  if(*p2++!=':') return 1;
  if(*p2=='-') {
    mf=1;
    p2++;
  }
  p=p2;
  while(b64conv[(int)*p2++]<64);
  y=0.0;
  *p1=--p2;
  while(p2-->p) {
    y+=(double)b64conv[(int)*p2];
    y=frexp(y,&j);
    y=ldexp(y,j-6);
  }
  if(mf) y=-y;
  *x=ldexp(y,e);
  return 0;
}
