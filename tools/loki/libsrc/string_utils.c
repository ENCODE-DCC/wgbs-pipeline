/****************************************************************************
 *                                                                          *
 *     Loki - Programs for genetic analysis of complex traits using MCMC    *
 *                                                                          *
 *                     Simon Heath - CNG, Evry                              *
 *                                                                          *
 *                         February 2003                                    *
 *                                                                          *
 * string_utils.c:                                                          *
 *                                                                          *
 * Copyright (C) Simon C. Heath 2003                                        *
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
#include <stdarg.h>
#include <ctype.h>
#include <math.h>
#include <assert.h>
#include <pthread.h>
#include "snprintf.h"

#define string_srep srep
#define _fget_bptr struct f_buffer **

#define __FG_BUF_SIZE 4096
#define __FG_STR_SIZE 256

struct f_buffer {
  size_t ptr,ptr1,size;
  char buf[__FG_BUF_SIZE];
};


#include "config.h"
#include "lk_malloc.h"
#include "string_utils.h"

#define SBUF(s) (s)->srep->buf

static int (*strcmpfunc)(const char *,const char *);

/* Max length of a printed real with %g format */
#define real_length(x) 13

struct srep {
  struct srep *next;
  size_t slen;
  size_t rfc;
  char *buf;
};

static pthread_mutex_t free_string_mutex = PTHREAD_MUTEX_INITIALIZER;
static pthread_mutex_t free_srep_mutex = PTHREAD_MUTEX_INITIALIZER;

static string *free_string_list;
static struct srep *free_srep_list;

static struct srep *new_srep(void)
{
  struct srep *sr;
	
  pthread_mutex_lock(&free_srep_mutex);
  if((sr=free_srep_list)) {
    free_srep_list=sr->next;
    pthread_mutex_unlock(&free_srep_mutex);
  } else {
    pthread_mutex_unlock(&free_srep_mutex);
    sr=lk_malloc(sizeof(struct srep));
  }
  if(sr) {
    sr->slen=0;
    sr->buf=0;
    sr->next=0;
    sr->rfc=0;
  }
  return sr;
}

static string *new_string(void)
{
  string *s;
	
  pthread_mutex_lock(&free_string_mutex);
  if((s=free_string_list)) {
    free_string_list=s->next;
    pthread_mutex_unlock(&free_string_mutex);
  } else {
    pthread_mutex_unlock(&free_string_mutex);
    s=lk_malloc(sizeof(string));
  }
  s->next=0;
  s->len=0;
  return s;
}

static char *alloc_sbuf(struct srep *sr,size_t sz)
{
  size_t ln;
	
  ln=sz+1;
  sr->slen=ln;
  sr->rfc=1;
  sr->buf=lk_malloc(ln);
  return sr->buf;
}

static char *realloc_sbuf(struct srep *sr,size_t sz)
{
  size_t ln;
	
  ln=sz+1;
  sr->buf=lk_realloc(sr->buf,ln);
  sr->slen=ln;
  return sr->buf;
}

static void free_sbuf(struct srep *sr)
{
  if(!(--(sr->rfc))) {
    if(sr->buf) {
      free(sr->buf);
      sr->buf=0;
    }
    sr->slen=0;
    pthread_mutex_lock(&free_srep_mutex);
    sr->next=free_srep_list;
    free_srep_list=sr;
    pthread_mutex_unlock(&free_srep_mutex);
  }
}

void free_string(string *s)
{
  if(s) {
    free_sbuf(s->srep);
    s->srep=0;
    pthread_mutex_lock(&free_string_mutex);
    s->next=free_string_list;
    free_string_list=s;
    pthread_mutex_unlock(&free_string_mutex);
  }
}

void free_string_lib(void)
{
  pthread_mutex_lock(&free_string_mutex);
  free_list(free_string_list,0);
  free_string_list=0;
  pthread_mutex_unlock(&free_string_mutex);
  pthread_mutex_lock(&free_srep_mutex);
  free_list(free_srep_list,0);
  free_srep_list=0;
  pthread_mutex_lock(&free_srep_mutex);
}

static struct srep *clean_copy(struct srep *sr)
{
  struct srep *sr1=0;
	
  sr1=new_srep();
  if(sr1) {
    alloc_sbuf(sr1,sr->slen);
    memcpy(sr1->buf,sr->buf,sr->slen);
    sr->rfc--;
  }
  return sr1;
}

void trim_comment(string *s,int c)
{
  size_t i;
  char *bf;
  struct srep *sr;
	
  if(s) {
    sr=s->srep;
    if(sr->rfc>1) s->srep=sr=clean_copy(sr);
    bf=sr->buf;
    for(i=0;i<s->len;i++) if(bf[i]==c) {
	bf[i]=0;
	s->len=i-1;
	break;
      }
  }
}

void *get_fget_buffer(struct f_buffer **buffer,size_t *sz)
{
  void *p;
  struct f_buffer *buf;
	
  buf=*buffer;
  if(buf->ptr>buf->ptr1) {
    p=buf->buf+buf->ptr1;
    *sz=buf->ptr-buf->ptr1;
  } else {
    p=0;
    *sz=0;
  }
  return p;    			 
}

string *fget_string(FILE *fptr,string *s,struct f_buffer **buffer)
{
  size_t i;
  int fg;
  struct srep *sr;
  struct f_buffer *buf;
  char *p;
	
  buf=*buffer;
  if(!buf) {
    buf=lk_malloc(sizeof(struct f_buffer));
    buf->ptr=buf->ptr1=0;
    buf->size=__FG_BUF_SIZE;
    *buffer=buf;
  }
  if(!s) {
    s=new_string();
    s->srep=sr=new_srep();
    sr->buf=lk_malloc((size_t)__FG_STR_SIZE);
    sr->slen=__FG_STR_SIZE;
  } else {
    sr=s->srep;
    if(sr->rfc>1) {
      sr->rfc--;
      s->srep=sr=new_srep();
      sr->slen=__FG_STR_SIZE;
      sr->buf=lk_malloc((size_t)__FG_STR_SIZE);
    }
    s->len=0;
  }
  sr->buf[0]=0;
  assert(buf->ptr>=buf->ptr1 && buf->ptr<=buf->size);
  p=buf->buf;
  if(buf->ptr==buf->ptr1) {
    buf->ptr1=0;
    /*    n=fread(p,1,buf->size,fptr);
    buf->ptr=n;
    if(!n) {
      s->len=0;
      return s;
      } */
    char *p1=fgets(p,buf->size,fptr);
    if(!p1) {
      s->len=0;
      return s;
    }
    buf->ptr=strlen(p1);
  }
  do {
    for(i=buf->ptr1;i<buf->ptr;i++) if(p[i]=='\n' || p[i]=='\r') break;
    if(i<buf->ptr) {
      fg=0;
      if(p[i]=='\r') {
	if(i+1==buf->ptr) fg=-1;
	else {
	  if(p[i+1]=='\n') fg=1;
	  p[i]='\n';
	}
      }
      if(fg>=0) {
	addl_cstring_to_string(s,p+buf->ptr1,i-buf->ptr1+1);
	buf->ptr1=i+1+fg;
	break;
      } else {
	if(i-buf->ptr1) addl_cstring_to_string(s,p+buf->ptr1,i-buf->ptr1);
	p[0]='\r';
	buf->ptr1=0;
	/*	n=fread(p+1,1,buf->size-1,fptr);
		buf->ptr=n+1; */
	char *p1=fgets(p+1,buf->size-1,fptr);
	
	buf->ptr=1+(p1?strlen(p1):0);
      }
    } else {
      addl_cstring_to_string(s,p+buf->ptr1,i-buf->ptr1);
      buf->ptr1=0;
      /*      n=fread(p,1,buf->size,fptr); */
      char *p1=fgets(p,buf->size,fptr);
      buf->ptr=p1?strlen(p1):0;
    }
  } while(buf->ptr>buf->ptr1);
  return s;
}

string *fget_string_gen(FILE *fptr,string *s,struct f_buffer **buffer,int c)
{
  size_t i,n;
  struct srep *sr;
  struct f_buffer *buf;
  char *p;
	
  buf=*buffer;
  if(!buf) {
    buf=lk_malloc(sizeof(struct f_buffer));
    buf->ptr=buf->ptr1=0;
    buf->size=__FG_BUF_SIZE;
    *buffer=buf;
  }
  if(!s) {
    s=new_string();
    s->srep=sr=new_srep();
    sr->buf=lk_malloc((size_t)__FG_STR_SIZE);
    sr->slen=__FG_STR_SIZE;
  } else {
    sr=s->srep;
    if(sr->rfc>1) {
      sr->rfc--;
      s->srep=sr=new_srep();
      sr->slen=__FG_STR_SIZE;
      sr->buf=lk_malloc((size_t)__FG_STR_SIZE);
    }
    s->len=0;
  }
  sr->buf[0]=0;
  assert(buf->ptr>=buf->ptr1 && buf->ptr<=buf->size);
  p=buf->buf;
  if(buf->ptr==buf->ptr1) {
    buf->ptr1=0;
    n=fread(p,1,buf->size,fptr);
    buf->ptr=n;
    if(!n) {
      s->len=0;
      return s;
    }
  }
  do {
    for(i=buf->ptr1;i<buf->ptr;i++) if(p[i]==c) break;
    if(i<buf->ptr) {
      addl_cstring_to_string(s,p+buf->ptr1,i-buf->ptr1+1);
      buf->ptr1=i+1;
      break;
    }
    addl_cstring_to_string(s,p+buf->ptr1,i-buf->ptr1);
    buf->ptr1=0;
    n=fread(p,1,buf->size,fptr);
    buf->ptr=n;
  } while(buf->ptr>buf->ptr1);
  return s;
}

void set_strcmpfunc(int (*f)(const char *,const char *))
{
  strcmpfunc=f;
}

int string_cmp(string *s1,string *s2)
{
  return strcmpfunc(SBUF(s1),SBUF(s2));
}

string *add_to_string(string *s,char c)
{
  char *bf;
  struct srep *sr;
	
  if(s) {
    sr=s->srep;
    if(sr->rfc>1) s->srep=sr=clean_copy(sr);
    bf=sr->buf;
    if(++s->len==sr->slen) {
      bf=realloc_sbuf(sr,sr->slen<<1);
    }
    bf[s->len-1]=c;
    bf[s->len]=0;
  } else {
    s=make_string(7,0);
    s->len=1;
    bf=SBUF(s);
    bf[0]=c;
    bf[1]=0;
  }
  return s;
}

string *addl_cstring_to_string(string *s,char *p,size_t sz)
{
  struct srep *sr;
  size_t i;
	
  if(p) {
    if(s) {
      sr=s->srep;
      if(sr->rfc>1) s->srep=sr=clean_copy(sr);
      if(sr->slen<s->len+sz+1) realloc_sbuf(sr,s->len+sz);
      memcpy(sr->buf+s->len,p,sz);
      s->len+=sz;
      sr->buf[s->len]=0;
    } else {
      i=(sz<7?7:sz);
      s=make_string(i,0);
      s->len=sz;
      memcpy(SBUF(s),p,sz);
      SBUF(s)[sz]=0;
    }
  } 
  return s;
}

string *add_cstring_to_string(string *s,char *p)
{
  return p?addl_cstring_to_string(s,p,strlen(p)):s;
}

string *init_string(string *s,char *p)
{
  struct srep *sr;
  size_t sz;
	
  sz=p?strlen(p):0;
  if(s) {
    sr=s->srep;
    if(sr->rfc>1) s->srep=sr=clean_copy(sr);
    if(sr->slen<sz+1) realloc_sbuf(sr,sz);
  } else s=make_string(sz<7?7:sz+1,0);
  s->len=sz;
  memcpy(SBUF(s),p,sz);
  SBUF(s)[sz]=0;
  return s;
}

char *copy_cstring(string *s)
{
  char *p=0;
	
  if(s) {
    p=s->srep->buf;
    s->srep->rfc++;
  }
  return p;
}

char *extract_cstring(string *s)
{
  char *p=0;
  struct srep *sr;
	
  if(s) {
    sr=s->srep;
    if(sr->rfc>1) {
      p=lk_malloc((size_t)s->len+1);
      strncpy(p,sr->buf,s->len+1);
    } else {
      p=sr->buf;
      sr->buf=0;
    }
    free_string(s);
  }
  return p;
}

string *copy_string1(string *s1,string *s)
{
  if(s1) free_string(s1);
  if(s) return copy_string(s);
  else return make_string(1,0);
}

string *strip_quotes(string *s)
{
  size_t l;
  char *bf;
  struct srep *sr;
	
  if(s) {
    l=s->len;
    if(l>=2) {
      sr=s->srep;
      if(sr->rfc>1) s->srep=sr=clean_copy(sr);
      bf=sr->buf;
      if(bf[0]=='\"' && bf[l-1]=='\"') {
	if(l>2) memmove(bf,bf+1,l-2);
	bf[l-2]=0;
	s->len-=2;
      }
    }
  }
  return s;
}

string *copy_string(string *s)
{
  string *s1=0;
	
  if(s) {
    s1=new_string();
    s1->srep=s->srep;
    s->srep->rfc++;
    s1->len=s->len;
  }
  return s1;
}

static size_t int_length(int i)
{
  size_t len=0;
	
  if(!i) len=1;
  else {
    if(i<0) {
      i=-i;
      len++; /* add space for unary - */
    }	
    len+=1+(int)(log((double)i)/LOG_TEN+1.0e-16);
  }	
  return len;
}

static string *lk_addn_to_string(string *s,char *p,size_t n,int fg)
{
  char *bf;
  struct srep *sr;
	
  if(s) {
    if(n) {
      sr=s->srep;
      if(sr->rfc>1) s->srep=sr=clean_copy(sr);
      if(s->len+n>=sr->slen) {
	realloc_sbuf(sr,s->len+n);
      }
      bf=sr->buf;
      if(fg) {
	if(s->len) memmove(bf+n,bf,s->len);
	memmove(bf,p,n);
      } else memmove(bf+s->len,p,n);
      s->len+=n;
      bf[s->len]=0;
    }
  } else {
    s=new_string();
    s->srep=new_srep();
    bf=alloc_sbuf(s->srep,n<7?7:n);
    s->len=n;
    if(n)	memmove(bf,p,n);
    bf[n]=0;
  }
  return s;
}

string *addn_to_string(string *s,char *p,size_t n)
{
  return lk_addn_to_string(s,p,n,0);
}

string *paddn_to_string(string *s,char *p,size_t n)
{
  return lk_addn_to_string(s,p,n,1);
}

string *int_to_string(int i)
{
  string *s;
  size_t j;
	
  j=int_length(i);
  s=make_string(j,1);
  snprintf(SBUF(s),j,"%d",i);
  return s;
}

string *real_to_string(double x)
{
  string *s;
  int j;
	
  j=real_length(x);
  s=make_string(j,1);
  s->len=snprintf(SBUF(s),j,"%g",x);
  return s;
}

string *add_strings(string *s,string *s1)
{
  struct srep *sr;
	
  if(!s) {
    return s1;
  }		
  if(!s1) {
    return s;
  }
  sr=s->srep;
  if(sr->rfc>1) s->srep=sr=clean_copy(sr);
  if(s->len+s1->len>=sr->slen) {
    realloc_sbuf(sr,s->len+s1->len);
  }
  strncpy(sr->buf+s->len,SBUF(s1),s1->len);
  s->len+=s1->len;
  free_string(s1);
  return s;
}

string *make_string(size_t i,int fg)
{
  string *s;
  char *p;
	
  if(!i) return 0;
  s=new_string();
  s->srep=new_srep();
  p=alloc_sbuf(s->srep,i);
  if(fg) {
    s->len=i;
    p[i--]=0;
    while(i) p[i--]=' ';
  } else {
    s->len=0;
    *p=0;
  }
  return s;
}

static string *check_str(string *s,size_t sz)
{
  struct srep *sr;
	
  if(s) {
    sr=s->srep;
    if(sr->rfc>1) {
      sr->rfc--;
      sr=new_srep();
      (void)alloc_sbuf(sr,sz);
      s->srep=sr;
    } else if(sz>=sr->slen) (void)realloc_sbuf(sr,sz);
  } else s=make_string(sz,0);
  return s;
}

string *strputc(const char c,string *s)
{
  char *p;
	
  s=check_str(s,1);
  p=SBUF(s);
  p[0]=c;
  p[1]=0;
  s->len=1;
  return s;
}

string *strputs(char *s1,string *s)
{
  size_t sz;
	
  assert(s1);
  sz=strlen(s1);
  s=check_str(s,sz);
  s->len=sz;
  strncpy(SBUF(s),s1,sz+1);
  return s;
}

string *strprintf(string *s,const char *fmt, ...)
{
  int sz;
	
  va_list args;
  /* get size */
  va_start(args,fmt);
  sz=vsnprintf(0,(size_t)0,fmt,args);
  va_end(args);
  assert(sz>=0);
  s=check_str(s,sz);
  va_start(args,fmt);
  (void)vsnprintf(get_cstring(s),(size_t)sz+1,fmt,args);
  s->len=sz;
  va_end(args);
  return s;
}

string *straprintf(string *s,const char *fmt, ...)
{
  size_t sz,off=0;
	
  va_list args;
  /* get size */
  va_start(args,fmt);
  sz=vsnprintf(0,(size_t)0,fmt,args);
  va_end(args);
  if(s) off=s->len;
  s=check_str(s,sz+off);
  va_start(args,fmt);
  (void)vsnprintf(get_cstring(s)+off,(size_t)sz+1,fmt,args);
  s->len=sz+off;
  va_end(args);
  return s;
}
