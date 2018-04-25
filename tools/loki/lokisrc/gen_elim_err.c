/****************************************************************************
 *                                                                          *
 *     Loki - Programs for genetic analysis of complex traits using MCMC    *
 *                                                                          *
 *                     Simon Heath - CNG, Evry                              *
 *                                                                          *
 *                         February 2003                                    *
 *                                                                          *
 * gen_elim_err.c:                                                          *
 *                                                                          *
 * Error reporting for genotype elimination code                            *
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
#include <math.h>
#include <stdio.h>

#include "utils.h"
#include "libhdr.h"
#include "version.h"
#include "loki.h"
#include "loki_peel.h"
#include "loki_utils.h"
#include "lk_malloc.h"
#include "gen_elim.h"

static struct loki *lk;

struct mark_err {
  struct Marker *mark;
  int n_err;
};

gen_elim_err *add_gen_elim_err(gen_elim_err *err,nuc_fam *fam,struct Id_Record *id,int type)
{
  gen_elim_err *new_err;
	
  new_err=lk_malloc(sizeof(gen_elim_err));
  new_err->fam=fam;
  new_err->id=id;
  new_err->err_type=type;
  new_err->next=err;
  return new_err;
}

gen_elim_err *invert_err_list(gen_elim_err *err)
{
  gen_elim_err *err1=0,*err2;
	
  while(err) {
    err2=err->next;
    err->next=err1;
    err1=err;
    err=err2;
  }
  return err1;
}

static int cmp_mark_err(const void *s1,const void *s2)
{
  int i;
  struct Marker *m1,*m2;

  i=((struct mark_err *)s2)->n_err-((struct mark_err *)s1)->n_err;
  if(!i) {
    m1=((struct mark_err *)s1)->mark;
    m2=((struct mark_err *)s2)->mark;
    i=(lk->markers->linkage[m1->locus.link_group].type&LINK_TYPES_MASK)-
      (lk->markers->linkage[m2->locus.link_group].type&LINK_TYPES_MASK);
  }
  return i;
}

void handle_errors(struct loki *loki,gen_elim_err **errors)
{
  static FILE *fptr;
  int i,i1,j,k,k1,kd,sz,nf,*perm,n_mark,linktype,nmk,fam_flag=0,*ngens,*sex_errs;
  nuc_fam *fam,*f1;
  gen_elim_err *err,*err1;
  struct Marker *mark;
  struct mark_err *errs;
  struct Id_Record *kid;
  unsigned long mask;
  char *ltype[]={"AUTO","X","Y","MT","Z","W"};
  gtype **gtypes;
	
  n_mark=loki->markers->n_markers;
  fam=loki->family->families;
  nf=loki->family->cp_start[loki->pedigree->n_comp];
  sex_errs=lk_calloc((size_t)loki->pedigree->ped_size,sizeof(int));
  errs=lk_malloc(sizeof(struct mark_err)*n_mark);
  for(i=0;i<nf;i++) {
    fam[i].n_err=0;
    for(j=0;j<N_LINK_TYPES;j++) fam[i].n_err1[j]=0;
  }
  if(!loki->names[LK_GENERRFILE]) loki->names[LK_GENERRFILE]=make_file_name(".err");
  j=loki->sys.syst_var[SYST_BACKUPS].flag?loki->sys.syst_var[SYST_BACKUPS].data.value:1;
  if(j) (void)mkbackup(loki->names[LK_GENERRFILE],j);
  if(!(fptr=fopen(loki->names[LK_GENERRFILE],"w"))) abt(__FILE__,__LINE__,"%s(): File Error.  Couldn't open '%s' for output\n",__func__,loki->names[LK_GENERRFILE]);
  (void)fprintf(fptr,"Created by %s: %s",LOKI_NAME,ctime(&loki->sys.lktime.start_time));
  fputs("\nDetected Mendelian errors\n",fptr);
  for(sz=nmk=i=0;i<n_mark;i++) if((err=errors[i])) {
    mark=loki->markers->marker+i;
    ngens=mark->ngens;
    gtypes=mark->gtypes;
    errs[nmk].mark=mark;
    errs[nmk].n_err=0;
    linktype=loki->markers->linkage[mark->locus.link_group].type&LINK_TYPES_MASK;
    if(linktype==LINK_AUTO || linktype==LINK_X) fam_flag=1;
    k=(int)strlen(mark->name)+1;
    if(k>sz) sz=k;
    j=k1=0;
    while(err) {
      errs[nmk].n_err++;
      if(err->err_type<GEN_ELIM_PASS2_KID) k=1;
      else k=2;
      if(j<k) {
	fputs("\n************************************************************\n",fptr);
	switch(linktype) {
	case LINK_AUTO:
	  fprintf(fptr,"* Marker %s, Pass %d\n",mark->name,k);
	  break;
	case LINK_X:
	  fprintf(fptr,"* Marker %s (X-linked), Pass %d\n",mark->name,k);
	  break;
	case LINK_Y:
	  fprintf(fptr,"* Marker %s (Y-linked), Pass %d\n",mark->name,k);
	  break;
	case LINK_MIT:
	  fprintf(fptr,"* Marker %s (Mitochondrial), Pass %d\n",mark->name,k);
	  break;
	case LINK_Z:
	  fprintf(fptr,"* Marker %s (Z-linked), Pass %d\n",mark->name,k);
	  break;
	case LINK_W:
	  fprintf(fptr,"* Marker %s (Z-linked), Pass %d\n",mark->name,k);
	  break;
	}
	fputs("************************************************************\n",fptr);
	j=k;
	k1=0;
      }
      if(!k1) fputc('\n',fptr);
      switch(err->err_type) {
      case GEN_ELIM_X_HET_MALE:
	(void)fputs("Individual ",fptr);
	print_orig_id(fptr,err->id->idx+1);
	sex_errs[err->id->idx]++;
	if(linktype==LINK_X) (void)fprintf(fptr," is male and heterozygous (X-linked marker)\n");
	else if(linktype==LINK_Z) (void)fprintf(fptr," is female and heterozygous (Z-linked marker)\n");
	else (void)fputs(" <BAD LINKTYPE ERROR>\n",fptr);
	k1=1;
	break;
      case GEN_ELIM_Y_HET_MALE:
	sex_errs[err->id->idx]++;
	(void)fputs("Individual ",fptr);
	print_orig_id(fptr,err->id->idx+1);
	(void)fprintf(fptr," is heterozygous (Y-linked marker)\n");
	k1=1;
	break;
      case GEN_ELIM_MIT_HET_FEMALE:
	sex_errs[err->id->idx]++;
	(void)fputs("Individual ",fptr);
	print_orig_id(fptr,err->id->idx+1);
	if(linktype==LINK_MIT) (void)fprintf(fptr," is heterozygous (Mitochondrial marker)\n");
	else if(linktype==LINK_W) (void)fprintf(fptr," is heterozygous (W-linked marker)\n");
	else (void)fputs(" <BAD LINKTYPE ERROR>\n",fptr);
	k1=1;
	break;
      case GEN_ELIM_Y_OBS_FEMALE:
	sex_errs[err->id->idx]++;
	(void)fputs("Individual ",fptr);
	print_orig_id(fptr,err->id->idx+1);
	(void)fprintf(fptr," is female and observed (Y-linked marker)\n");
	k1=1;
	break;
      case GEN_ELIM_MIT_OBS_MALE:
	sex_errs[err->id->idx]++;
	(void)fputs("Individual ",fptr);
	print_orig_id(fptr,err->id->idx+1);
	if(linktype==LINK_MIT) (void)fprintf(fptr," is male and observed (Mitochondrial marker)\n");
	else if(linktype==LINK_W) (void)fprintf(fptr," is male and observed (W-linked marker)\n");
	k1=1;
	break;
      case GEN_ELIM_HALF_OBS:
	(void)fputs("Individual ",fptr);
	print_orig_id(fptr,err->id->idx+1);
	(void)fprintf(fptr," is half observed\n");
	k1=1;
	break;
      case GEN_ELIM_PASS1:
	err->fam->n_err++;
	err->fam->n_err1[linktype]++;
	fputs("Error when adding child ",fptr);
	print_orig_id(fptr,err->id->idx+1);
	fputc('\n',fptr);
	print_family_gtypes(fptr,err->fam,mark,0,0);
	break;
      case GEN_ELIM_PASS2_KID:
	err->fam->n_err++;
	err->fam->n_err1[linktype]++;
	fputs("Error when processing child ",fptr);
	print_orig_id(fptr,err->id->idx+1);
	fputc('\n',fptr);
	print_family_gtypes(fptr,err->fam,mark,ngens,gtypes);
	break;
      case GEN_ELIM_PASS2_PAR:
	err->fam->n_err++;
	err->fam->n_err1[linktype]++;
	fprintf(fptr,"Error when processing %s\n",err->id->sex==2?"mother":"father");
	print_family_gtypes(fptr,err->fam,mark,ngens,gtypes);
	break;
      case GEN_ELIM_PASS2_YM:
	fprintf(fptr,"Mismatch in %s line descended from ",linktype==LINK_Y?"male":"female");
	print_orig_id(fptr,err->id->idx+1);
	fputs("\n\n",fptr);
	kd=-1;
	while((kid=err->fam->kids[++kd])) {
	  if(ngens[kid->idx]) {
	    print_orig_id(fptr,kid->idx+1);
	    fputs(" (",fptr);
	    print_orig_gtypes(fptr,kid,mark);
	    fputs(")\n",fptr);
	  }
	}
	kid=err->fam->kids[0];
			if(kid->family) {
				kid->family->n_err++;
				kid->family->n_err1[linktype]++;
			}
	free(err->fam->kids);
	free(err->fam);
	break;
      default:
	ABT_FUNC("Internal error - bad error code\n");
      }
      err1=err->next;
      free(err);
      err=err1;
    }
    nmk++;
  }
  fputs("\n***************************************************\nBreak down of errors by marker\n\n",fptr);
  lk=loki;
  gnu_qsort(errs,nmk,sizeof(struct mark_err),cmp_mark_err);
  fputs("Marker",fptr);
  if(sz<6) sz=6;
  for(i=0;i<sz+2-6;i++) fputc(' ',fptr);
  fputs("Linkage  No. errors\n\n",fptr); 
  for(i=0;i<nmk;i++) {
    mark=errs[i].mark;
    linktype=loki->markers->linkage[mark->locus.link_group].type&LINK_TYPES_MASK;
    fprintf(fptr,"%-*s  %-4s  %d\n",sz,mark->name,ltype[linktype],errs[i].n_err);
  }
  if(fam_flag) {
    fputs("\n***************************************************\nBreak down of errors by nuclear family\n\n",fptr);
    sz=get_max_idlen()+1;
    if(sz<8) sz=8;
    fputs("Father",fptr);
    for(i=0;i<sz-6;i++) fputc(' ',fptr);
    fputs("Mother",fptr);
    for(i=0;i<sz-6;i++) fputc(' ',fptr);
    perm=lk_malloc(sizeof(int)*nf);
    mask=0L;
    for(i=j=0;j<nf;j++) if(fam[j].n_err) {
      perm[i++]=j;
      for(k=0;k<N_LINK_TYPES;k++) if(fam[j].n_err1[k]) mask|=(1L<<k);
    }
    if(mask==1L) {
      fputs("No. errors\n\n",fptr);
      mask=0L;
    } else {
      for(k=0;k<N_LINK_TYPES;k++) if(mask&(1L<<k)) {
	fprintf(fptr,"%-8s",ltype[k]);
      }
      fputs("Total errors\n\n",fptr);
    }
    sort_fam_error(perm,i,loki);
    for(j=0;j<i;j++) {
      f1=fam+perm[j];
      k=print_orig_id(fptr,f1->father->idx+1);
      for(i1=0;i1<sz-k;i1++) fputc(' ',fptr);
      k=print_orig_id(fptr,f1->mother->idx+1);
      for(i1=0;i1<sz-k;i1++) fputc(' ',fptr);
      if(mask) {
	for(k=0;k<N_LINK_TYPES;k++) if(mask&(1L<<k)) {
	  fprintf(fptr,"%-8d",f1->n_err1[k]);
	}
      }
      fprintf(fptr,"%d\n",f1->n_err);
    }
    free(perm);
  }
  for(i=j=0;i<loki->pedigree->ped_size;i++) if(sex_errs[i]) j++;
  if(j) {
    perm=lk_malloc(sizeof(int)*j);
    for(i=j=0;i<loki->pedigree->ped_size;i++) if(sex_errs[i]) perm[j++]=i;
    sort_ind_error(perm,j,sex_errs);
    fputs("\n***************************************************\nErrors in sex-limited markers by individual\n\n",fptr);
    sz=get_max_idlen()+2;
    if(sz<12) sz=12;
    fprintf(fptr,"%-*sSex  No. Errors\n\n",sz,"Individual");
    for(k=0;k<j;k++) {
      i=perm[k];
      k1=print_orig_id(fptr,i+1);
      for(i1=0;i1<sz-k1;i1++) fputc(' ',fptr);
      fprintf(fptr,"%c    %d\n",loki->pedigree->id_array[i].sex==2?'F':'M',sex_errs[i]);
    }
    free(perm);
  }
  free(errs);
  fclose(fptr);
}

