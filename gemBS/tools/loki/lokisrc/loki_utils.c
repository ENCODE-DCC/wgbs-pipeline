/****************************************************************************
 *                                                                          *
 *     Loki - Programs for genetic analysis of complex traits using MCMC    *
 *                                                                          *
 *             Simon Heath - University of Washington                       *
 *                                                                          *
 *                      August 1997                                         *
 *                                                                          *
 * loki_utils.c:                                                            *
 *                                                                          *
 * Utility routines for loki                                                *
 *                                                                          *
 * Copyright (C) Simon C. Heath 1997, 2000, 2002                            *
 * This is free software.  You can distribute it and/or modify it           *
 * under the terms of the Modified BSD license, see the file COPYING        *
 *                                                                          *
 ****************************************************************************/

#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <string.h>
#include <math.h>
#include <stdarg.h>
#include <assert.h>

#include "utils.h"
#include "loki.h"
#include "loki_peel.h"
#include "string_utils.h"
#include "loki_utils.h"

static int sort_sex;
static const struct loki *loki;

void init_utils(const struct loki *lk)
{
  loki=lk;
  init_lib_utils();
  return;
}

#ifdef TRACE_PEEL
void prt_peel_trace(const int level, const char *fmt, ...)
{
  va_list args;
	
  if(CHK_PEEL(level)) {
    va_start(args,fmt);
    (void)vfprintf(stdout,fmt,args);
    va_end(args);
  }
}
#else 
void prt_peel_trace(const int level _U_, const char *fmt _U_, ...)
{
}
#endif
size_t print_orig_family(FILE *fptr,const int id,const int fg)
{
  int fam;
  size_t sz;
	
  fam=id>0?loki->pedigree->id_array[id-1].fam_code:0;
  if(!fam) sz=1;
  else if(loki->pedigree->fam_recode.flag==ST_STRING) sz=strlen(loki->pedigree->fam_recode.recode[fam-1].string);
  else sz=(size_t)(1.000001+log((double)loki->pedigree->fam_recode.recode[fam-1].value)/log(10.0));
  if(fg) sz+=3;
  if(fptr) {
    if(fg) {
      if(!fam) (void)fputs("[*]:",fptr);
      else if(loki->pedigree->fam_recode.flag==ST_STRING) (void)fprintf(fptr,"[%s]:",loki->pedigree->fam_recode.recode[fam-1].string);
      else (void)fprintf(fptr,"[%d]:",loki->pedigree->fam_recode.recode[fam-1].value);
    } else {
      if(!fam) (void)fputs("*",fptr);
      else if(loki->pedigree->fam_recode.flag==ST_STRING) (void)fprintf(fptr,"%s",loki->pedigree->fam_recode.recode[fam-1].string);
      else (void)fprintf(fptr,"%d",loki->pedigree->fam_recode.recode[fam-1].value);
    }
  }
  return sz;
}

void print_marker_name(FILE *fptr,const int i)
{
  fputs(loki->markers->marker[i].name,fptr);
}

size_t print_orig_id1(FILE *fptr,int i)
{
  size_t sz;
  int j;
  double z;
	
  assert(i>=0 && i<=loki->pedigree->ped_size);
  if(!i) sz=1;
  else if(!loki->pedigree->id_recode.recode[i-1].string) {
    i=0;
    sz=1;
  } else if(loki->pedigree->id_recode.flag==ST_STRING) sz=strlen(loki->pedigree->id_recode.recode[i-1].string);
  else {
    j=loki->pedigree->id_recode.recode[i-1].value;
    z=(double)abs(j);
    sz=(size_t)(1.000001+log(z)/log(10.0));
    if(j<0) sz++;
  }
  if(fptr) {
    if(!i) (void)fputc('*',fptr);
    else {
      if(loki->pedigree->id_recode.flag==ST_STRING) (void)fprintf(fptr,"%s",loki->pedigree->id_recode.recode[i-1].string);
      else (void)fprintf(fptr,"%d",loki->pedigree->id_recode.recode[i-1].value);
    }
  }
  return sz;
}

size_t print_orig_id(FILE *fptr,int i)
{
  size_t sz=0;
	
  if(loki->pedigree->family_id) sz=print_orig_family(fptr,i,1);
  sz+=print_orig_id1(fptr,i);
  return sz;
}

string *str_print_orig_family(string *s,const int i)
{
  int fam;

  assert(i>=0 && i<loki->pedigree->ped_size && loki->pedigree->family_id);
  fam=i?loki->pedigree->id_array[i-1].fam_code:0;
  if(!fam) s=strputc('*',s);
  else {
    if(loki->pedigree->fam_recode.flag==ST_STRING) s=strputs(loki->pedigree->fam_recode.recode[fam-1].string,s);
    else s=strprintf(s,"%d",loki->pedigree->fam_recode.recode[fam-1].value);
  }
  return s;
}

string *str_print_orig_id(string *s,int i)
{
  assert(i>=0 && i<=loki->pedigree->ped_size);
  if(i && !loki->pedigree->id_recode.recode[i-1].string) i=0;
  if(!i) s=strputc('*',s);
  else {
    if(loki->pedigree->id_recode.flag==ST_STRING) s=strputs(loki->pedigree->id_recode.recode[i-1].string,s);
    else s=strprintf(s,"%d",loki->pedigree->id_recode.recode[i-1].value);
  }
  return s;
}

size_t get_max_idlen(void)
{
  int i;
  size_t sz;
  static size_t max;
	
  if(!max) {
    for(i=0;i<loki->pedigree->ped_size;i++) {
      sz=(int)print_orig_id(0,i+1);
      if(sz>max) max=sz;
    }
  }
  return max;
}

void print_orig_triple(FILE *fptr,const int i)
{
  if(loki->pedigree->family_id) {
    (void)print_orig_family(fptr,i,0);
    (void)fputc(' ',fptr);
  }
  if(i)	{
    (void)print_orig_id1(fptr,i);
    (void)fputc(' ',fptr);
    (void)print_orig_id1(fptr,loki->pedigree->id_array[i-1].sire);
    (void)fputc(' ',fptr);
    (void)print_orig_id1(fptr,loki->pedigree->id_array[i-1].dam);
  }
  else (void)fputs("* * *",fptr);
}

void print_orig_allele_id(FILE *fptr,const int i)
{
  if(i>0) {
    (void)print_orig_id(fptr,i);
    (void)fputc('m',fptr);
  } else {
    (void)print_orig_id(fptr,-i);
    (void)fputc('p',fptr);
  }
}

void print_allele_type1(FILE *fptr,const int locus,const int j)
{
  struct Marker *mk;
	
  mk=loki->markers->marker+locus;
  if(mk->recode[j]) (void)fputs(mk->recode[j],fptr);
  else (void)fputs(LUMPED_ALLELE,fptr);
}

/* Print the original allele code for (recoded) allele i */
void print_allele_type(FILE *fptr,const int locus,const int comp,const int i)
{
  print_allele_type1(fptr,locus,loki->markers->marker[locus].allele_trans[comp][i]);
}

/* Print original code for maternal or paternal allele (depending on flag) for
 * individual id.  Takes account of set recoding */
void print_allele_name(FILE *fptr,const int id,const int locus,const int flag)
{
  lk_ulong c;
  int i,j,allele,comp;
	
  comp=loki->pedigree->id_array[id].comp;
  allele=loki->pedigree->id_array[id].allele[flag]-1;
  c=1<<allele;
  if(c&loki->markers->marker[locus].req_set[flag][id]) {
    i=j=0;
    c=loki->markers->marker[locus].req_set[flag][id];
    while(c)	{
      if(c&1) {
	(void)fputc(j++?',':'[',fptr);
	print_allele_type(fptr,locus,comp,i);
      }
      i++;
      c>>=1;
    }
    (void)fputc(']',fptr);
  } else print_allele_type(fptr,locus,comp,allele);
}

void set_sort_sex(const int s)
{
  sort_sex=s;
}

/* Comparison function for qsort(), used to sort loci into position
 * order within a linkage group */
int cmp_loci(const void *s1,const void *s2)
{
  const struct Locus *l1,*l2;
  double x1,x2;
	
  l1=*((const struct Locus **)s1);
  l2=*((const struct Locus **)s2);
  x1=l1->pos[sort_sex];
  x2=l2->pos[sort_sex];
  if(x1<x2) return -1;
  if(x1>x2) return 1;
  return 0;
}

/* Get list of loci (marker + trait) in linkage group.  Returns number found
 * in count */
void get_locuslist(struct Locus **list,const int link,int *count,int flag)
{
  int i,j,k;
	
  for(i=k=0;i<loki->markers->linkage[link].n_markers;i++) {
    j=loki->markers->linkage[link].mk_index[i];
    list[i]=&loki->markers->marker[j].locus;
  }
  if(!flag) for(j=0;j<loki->params.n_tloci;j++) if((loki->models->tlocus[j].flag&TL_LINKED) && loki->models->tlocus[j].link_group==link)
    list[i++]=&loki->models->tlocus[j];
  *count=i;
}

void print_allele_type2(FILE *fptr,struct Marker *mark,int comp,int i)
{
  int j;
  
  if(i) {
    if(i<=0 || i>mark->locus.n_alleles) {
      fprintf(fptr,"<OOPS %d>",i);
    } else {
		 assert(i>0 && i<=mark->locus.n_alleles); 
		 j=mark->allele_trans[comp][i-1];
		 if(j<0 || j>=mark->locus.n_alleles) {
			 fprintf(fptr,"<ACK! %d %d>",i,j);
		 } else {
			 assert(j>=0 && j<mark->locus.n_alleles);
			 if(mark->recode[j]) (void)fputs(mark->recode[j],fptr);
			 else (void)fputs(LUMPED_ALLELE,fptr);
		 }
	 }
  } else (void)fputc('*',fptr);
}

static int sxdef[][2]={{0,0}, /* LINK_AUTO */
		       {1,2}, /* LINK_X */
		       {2,1}, /* LINK_Y */
		       {2,2}, /* LINK_MIT */
		       {1,1}, /* LINK_Z */
		       {2,2}  /* LINK_W */
};

void print_derived_gtypes(FILE *fptr,struct Id_Record *id,struct Marker *mark,int *ngens,gtype **gtypes)
{
  int i,c,g,k,idx,ng,linktype,sx_ix,sx;

  linktype=loki->markers->linkage[mark->locus.link_group].type&LINK_TYPES_MASK;
  sx_ix=sxdef[linktype][0];
  sx=sxdef[linktype][1];
  if(sx_ix==2 && id->sex!=sx) return;
  idx=id->idx;
  ng=ngens[idx];
  if(!ng) {
    if(!sx_ix || (sx_ix==1 && id->sex==sx)) fputs(" [*,*]",fptr);
    else fputs(" [*]",fptr);
  } else {
    c=id->comp;
    for(i=0;i<ng;i++) {
      fputs(" [",fptr);
      k=0;
      if(sx!=1 || id->sex==1) {
	g=gtypes[idx][i].mat;
	if(g) print_allele_type2(fptr,mark,c,g);
	else fputc('*',fptr);
	k=1;
      }
      if(sx!=2 || id->sex==2) {
	if(k) fputc(',',fptr);
	g=gtypes[idx][i].pat;
	if(g) print_allele_type2(fptr,mark,c,g);
	else fputc('*',fptr);
      }
      fputc(']',fptr);
    }
  }
}

void print_orig_gtypes(FILE *fptr,struct Id_Record *id,struct Marker *mark)
{
  int k,c,g1,g2,linktype,sx,sx_ix;
	
  linktype=loki->markers->linkage[mark->locus.link_group].type&LINK_TYPES_MASK;
  sx_ix=sxdef[linktype][0];
  sx=sxdef[linktype][1];
  if(sx_ix==2 && id->sex!=sx) return;
  k=mark->haplo[id->idx];
  if(!k) {
    if(!sx_ix || (sx_ix==1 && id->sex==sx)) fputs("*,*",fptr);
    else fputs("*",fptr);
  } else {
    g1=(int)(sqrt(.25+2.0*(double)k)-.49999);
    g2=k-(g1*(g1+1)/2);
    c=id->comp;
    print_allele_type2(fptr,mark,c,g1);
    if(!sx_ix || (sx_ix==1 && id->sex==sx)) {
      fputc(',',fptr);
      print_allele_type2(fptr,mark,c,g2);
    }
  }
}

void print_family_gtypes(FILE *fptr,nuc_fam *fam,struct Marker *mark,int *ngens,gtype **gtypes)
{
  int i=0;
  struct Id_Record *par,*kid;
	
  par=fam->father;
  if(!par) ABT_FUNC("Internal error - no father\n");
  fputs("Father: ",fptr);
  print_orig_id(fptr,par->idx+1);
  fputs(" (",fptr);
  print_orig_gtypes(fptr,par,mark);
  fputc(')',fptr);
  if(ngens && gtypes) print_derived_gtypes(fptr,par,mark,ngens,gtypes);
  fputc('\n',fptr);
  par=fam->mother;
  if(!par) ABT_FUNC("Internal error - no father\n");
  fputs("Mother: ",fptr);
  print_orig_id(fptr,par->idx+1);
  fputs(" (",fptr);
  print_orig_gtypes(fptr,par,mark);
  fputc(')',fptr);
  if(ngens && gtypes) print_derived_gtypes(fptr,par,mark,ngens,gtypes);
  fputc('\n',fptr);
  while((kid=fam->kids[i++])) {
    fputs("  ",fptr);
    print_orig_id(fptr,kid->idx+1);
    fputs(" (",fptr);
    print_orig_gtypes(fptr,kid,mark);
    fputc(')',fptr);
    if(ngens && gtypes) print_derived_gtypes(fptr,kid,mark,ngens,gtypes);
    fputc('\n',fptr);
  }
}

static void print_flags(FILE *fptr,int flags)
{
  int i,k=0,tst[]={IN_RF,HAP_JNT,HAD_M,HAD_P,FIXED_FLAG,HAP_FND};
  int c[]={'R','J','M','P','F','A'};
	
  for(i=0;i<6;i++) if(flags&tst[i]) {
    if(!k++) fputc('[',fptr);
    fputc(c[i],fptr);
  }
  if(k) fputc(']',fptr);
}

void print_peel_elem(FILE *fptr,void *elem,int type,int sample)
{ 
  int j,k,k1;
  struct Simple_Element *s_elem;
  struct Complex_Element *c_elem;
	
  if(type==PEEL_SIMPLE) {
    s_elem=elem;
    if(s_elem->sire) {
      if(sample) {
	(void)fputs("--> Sampling family: ",fptr);
	if(loki->pedigree->family_id) print_orig_family(fptr,abs(s_elem->sire),1);
	if(s_elem->pivot==-2) {
	  fputs("[Parents]",fptr);
	} else if(s_elem->pivot==-3) {
	  fputs("[Kids]",fptr);
	} else if(!s_elem->pivot) (void)fputc('*',fptr);
	else print_orig_id1(fptr,s_elem->pivot);
	(void)fputs(" ->",fptr);
	print_orig_id1(fptr,abs(s_elem->sire));
	(void)fputc(',',fptr);
	print_orig_id1(fptr,abs(s_elem->dam));
	(void)fputc(' ',fptr);
	for(j=0;j<s_elem->n_off;j++) {
	  (void)fputc(j?',':'(',fptr);
	  print_orig_id1(fptr,abs(s_elem->off[j]));
	}
	fputc(')',fptr);
      } else {
	(void)fputs("--> Peeling family: ",fptr);
	if(loki->pedigree->family_id) print_orig_family(fptr,abs(s_elem->sire),1);
	print_orig_id1(fptr,abs(s_elem->sire));
	(void)fputc(',',fptr);
	print_orig_id1(fptr,abs(s_elem->dam));
	(void)fputc(' ',fptr);
	for(j=0;j<s_elem->n_off;j++) {
	  (void)fputc(j?',':'(',fptr);
	  print_orig_id1(fptr,abs(s_elem->off[j]));
	}
	(void)fputs(") -> ",fptr);
	if(s_elem->pivot==-2) {
	  fputs("[Parents]",fptr);
	} else if(s_elem->pivot==-3) {
	  fputs("[Kids]",fptr);
	} else if(!s_elem->pivot) (void)fputc('*',fptr);
	else print_orig_id1(fptr,s_elem->pivot);
      }
    } else {
      if(sample) (void)fputs("--> Sampling singletons:",fptr);
      else (void)fputs("--> Peeling singletons:",fptr);
      for(j=0;j<s_elem->n_off;j++) {
	fputs("\n  ",fptr);
	print_orig_id1(fptr,s_elem->off[j]);
      }
    }
    (void)fputc('\n',fptr);
  } else if(type==PEEL_COMPLEX) {
    c_elem=elem;
    if(sample) (void)fputs("--> Complex sample:",fptr);
    else (void)fputs("--> Complex peel:",fptr);
    if(loki->pedigree->family_id) {
      fputc(' ',fptr);
      print_orig_family(fptr,c_elem->involved[0],1);
    }
    for(j=0;j<c_elem->n_peel;j++) {
      fputc(j?',':' ',fptr);
      k=c_elem->involved[j];
      print_orig_id1(fptr,abs(k));
      fputc(k<0?'p':'m',fptr);
      print_flags(fptr,c_elem->flags[j]);
    }
    for(k1=0;j<c_elem->n_involved-c_elem->n_out;j++) {
      if(!k1++) fputs(" {",fptr);
      else fputc(' ',fptr);
      k=c_elem->involved[j];
      print_orig_id1(fptr,abs(k));
      fputc(k<0?'p':'m',fptr);
      print_flags(fptr,c_elem->flags[j]);
    }
    if(k1) fputc('}',fptr);
    fputs(" ->",fptr);
    for(k1=0;j<c_elem->n_involved;j++) {
      fputc(k1++?',':' ',fptr);
      k=c_elem->involved[j];
      print_orig_id1(fptr,abs(k));
      fputc(k<0?'p':'m',fptr);
      print_flags(fptr,c_elem->flags[j]);
    }
    if(!k1) fputs(" *",fptr);
    else fprintf(fptr," [%d]",c_elem->out_index);
    for(j=0;j<c_elem->n_rfuncs;j++) {
      k=c_elem->index[j];
      if(!j) fputs(" (",fptr);
      else fputc(',',fptr);
      fprintf(fptr,"%d",k);
    }
    if(j) fputc(')',fptr);
    fputc('\n',fptr);
  } else {
    fprintf(fptr,"print_peel_elem() not implemented for this type\n");
  }
}

void print_peelseq(FILE *fptr,struct Peelseq_Head *peel)
{
  struct Simple_Element *simple_em;
  struct Complex_Element *complex_em;

  int i;
  struct Peelseq_Head *pp;
  for(i=0;i<loki->pedigree->n_comp;i++) {
    pp=peel+i;
    while(pp->type) {
      if(pp->type==PEEL_SIMPLE) {
	simple_em=pp->ptr.simple;
	print_peel_elem(fptr,simple_em,PEEL_SIMPLE,0);
	pp=&simple_em->next;
      } else {
	complex_em=pp->ptr.complex;
	print_peel_elem(fptr,complex_em,PEEL_COMPLEX,0);
	pp=&complex_em->next;
      }
    }
  }
}

