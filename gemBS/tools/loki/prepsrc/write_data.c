/***************************************************************************
 *                                                                          *
 *     Loki - Programs for genetic analysis of complex traits using MCMC    *
 *                                                                          *
 *             Simon Heath - University of Washington                       *
 *                                                                          *
 *                      June 1997                                           *
 *                                                                          *
 * write_data.c:                                                            *
 *                                                                          *
 * Write out data files for use by Loki                                     *
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
#include <stdio.h>
#include <math.h>
#include <errno.h>
#include <sys/wait.h>
#ifdef HAVE_LIMITS_H
#include <limits.h>
#endif
#include <ctype.h>
#include <assert.h>

#include "version.h"
#include "prep_utils.h"
#include "utils.h"
#include "libhdr.h"
#include "lk_malloc.h"
#include "scan.h"
#include "y.tab.h"
#include "write_data.h"

static void escape_space(const char *s,FILE *fp)
{
  while(*s) {
    if(isspace((int)*s)) fputc('\\',fp);
    fputc(*s++,fp);
  }
}

static void write_phen_var(int count,struct var_element **list,char *label,FILE *fptr)
{
  int i,j,k,k1;
  struct scan_data *sd;
	
  (void)fprintf(fptr," <group name='%s'>\n",label);
  for(i=0;i<count;i++) {
    sd=list[i]->arg.var->data;
    fprintf(fptr,"  <var name='%s",sd->name);
    if(sd->vtype&ST_ARRAY) fprintf(fptr,"(%d)",list[i]->oindex);
    j=list[i]->type;
    if(j&ST_CENSORED) fputs("' flags='censored",fptr);
    fputs("' type=",fptr);
    if(j&ST_FACTOR) {
      fputs("'factor'",fptr);
      k=list[i]->n_levels;
      if(k) {
	fputs(">\n   <levels>",fptr);
	for(j=0;j<n_factors;j++) if(var_factors[j]==list[i]) break;
	assert(j<n_factors);
	if(factor_recode[j][0]->type==STRING) {
	  for(k1=0;k1<k;k1++) {
	    if(k1) fputc(' ',fptr);
	    escape_space(factor_recode[j][k1]->data.string,fptr);
	  }
	} else {
	  for(k1=0;k1<k;k1++) {
	    if(k1) fputc(' ',fptr);
	    fprintf(fptr,"%ld",factor_recode[j][k1]->data.value);
	  }
	}
	fputs("</levels>\n  </var>\n",fptr);
      } else fputs("/>\n",fptr);
    } else {
      if(j&ST_INTEGER) fputs("'integer'/>\n",fptr);
      else fputs("'real'/>\n",fptr);
    }
  }
  fputs(" </group>\n",fptr);
}

static void write_phen_data(int count,int flag,struct var_element **list,int *list1,char *label,FILE *fptr)
{
  int i,j,k,fg,rec,nrec,type;
  struct id_data *s;
	
  (void)fprintf(fptr," <data group='%s'>\n",label);
  for(i=0;i<ped_size;i++) if(ped_recode1[i]) {
    nrec=flag?id_array[i].nrec:1;
    for(rec=0;rec<nrec;rec++) {
      fg=0;
      for(j=0;j<count;j++) {
	k=list1[j];
	type=list[j]->type;
	if(k<n_id_records) {
	  if(id_array[i].data) {
	    if(!fg++) fprintf(fptr,"  <rec id='%d'>",ped_recode1[i]);
	    else fputc(',',fptr);
	    s=id_array[i].data+k;
	    if(s->flag&2)fputs("<cens>",fptr);
	    if((type&ST_FACTOR)|(s->flag&ST_INTTYPE)) fprintf(fptr,"%ld",s->data.value);
	    else fprintf(fptr,"%g",s->data.rvalue);
	    if(s->flag&2)fputs("</cens>",fptr);
	  } 
	} else {
	  k-=n_id_records;
	  if(id_array[i].data1 && id_array[i].data1[rec]) {
	    if(!fg++) fprintf(fptr,"  <rec id='%d'>",ped_recode1[i]);
	    else fputc(',',fptr);
	    s=id_array[i].data1[rec]+k;
	    if(s->flag&2)fputs("<cens>",fptr);
	    if((type&ST_FACTOR)|(s->flag&ST_INTTYPE)) fprintf(fptr,"%ld",s->data.value);
	    else fprintf(fptr,"%g",s->data.rvalue);
	    if(s->flag&2)fputs("</cens>",fptr);
	  }
	}
      }
      if(fg) fputs("</rec>\n",fptr);
    }
  }
  (void)fputs(" </data>\n",fptr);
}

void WriteXMLData(char *LogFile,struct lk_compress *lkc)
{
  int i,j,k,k1,k2,n_orig_fam1,id,id1,ids,idd,link,m_id;
  int id_count=0,nonid_count=0,*trans;
  int *family_recode1=0,*perm,comp,*comp_to_family=0,*elem_list1=0,*comp_start;
  char *fname=0,*map_str[]={"gen_pos_female","gen_pos_male","gen_pos"};
  struct Link *linkp;
  struct var_element **elem_list=0,*elem;
  struct scan_data *sd;
  struct model *model;
  struct model_list *mlist;
  static char *labels[]={"markers","single_rec","multiple_rec","system"};
  char *filter;
  unsigned long tm;
  FILE *fptr;
	
  errno=0;
  filter=Filter;
  if(!filter && (lkc->default_compress<COMPRESS_NONE)) {
    j=lkc->default_compress;
    filter=lkc->comp_path[j][j==COMPRESS_ZIP?1:0];
  }
  fname=make_file_name(".xml");
  if(filter) {
    i=child_open(WRITE,fname,filter);
    if(!(fptr=fdopen(i,"w"))) DataFileError(fname);
    if(errno && errno!=ESPIPE) DataFileError(fname);
    errno=0;
  } else if(!(fptr=fopen(fname,"w"))) abt(__FILE__,__LINE__,"%s(): File Error.  Couldn't open '%s' for writing\n",__func__,fname);
  (void)fputs("<?xml version='1.0' encoding='UTF-8' standalone='no'?>\n",fptr);
  (void)fputs("<!DOCTYPE loki SYSTEM 'http://loki.homeunix.net/loki.dtd'>\n",fptr);
  (void)fprintf(fptr,"<loki program='%s'>\n",PREP_NAME); /* Start tag */
  if(LogFile && *LogFile) {
    (void)fprintf(fptr," <log file='%s'/>\n",LogFile);
  }
  /* Pedigree structure */
  perm=lk_malloc(sizeof(int)*pruned_ped_size);
  comp_to_family=lk_calloc((size_t)n_comp*2,sizeof(int));
  comp_start=comp_to_family+n_comp;
  for(id=comp=0;comp<n_comp;comp++) {
    comp_start[comp]=id;
    id+=comp_size[comp];
  }
  for(i=0;i<ped_size;i++) {
    j=ped_recode1[i];
    if(!j) continue;
    perm[j-1]=i;
  }
  if(family_id) {
    family_recode1=lk_calloc((size_t)n_orig_families,sizeof(int));
    for(i=0;i<ped_size;i++) if(ped_recode1[i]) {
      k=id_array[i].fam_code;
      if(k) family_recode1[k-1]=1;
    }
    n_orig_fam1=0;
    for(i=0;i<n_orig_families;i++) if(family_recode1[i]) family_recode1[i]=++n_orig_fam1;
    id=0;
    for(comp=0;comp<n_comp;comp++) {
      j=comp_size[comp];
      i=perm[id+j-1];
      k=id_array[i].fam_code;
      comp_to_family[comp]=k-1;
      id+=j;
    }
  } 
  j=family_id?n_orig_families:1;
  id=0;
  for(i=0;i<j;i++) {
    if(family_id) {
      if(!family_recode1[i]) continue;
      (void)fputs(" <family orig='",fptr);
      if(family_recode[i]->type==STRING) (void)fprintf(fptr,"%s'>\n",family_recode[i]->data.string);
      else (void)fprintf(fptr,"%d'>\n",(int)family_recode[i]->data.value);
    } else (void)fputs(" <family>\n",fptr);
    for(k1=comp=0;comp<n_comp;comp++) if(comp_to_family[comp]==i) k1++;
    k1++; /* Temporary until Loki can cope with missing comp */
    for(comp=0;comp<n_comp;comp++) if(comp_to_family[comp]==i) {
      if(k1>1) (void)fputs("  <comp>\n",fptr);
      id=comp_start[comp];
      for(k=0;k<comp_size[comp];k++,id++) {
	(void)fprintf(fptr,"%s<ind id='%d'",k1>1?"   ":"  ",id+1);
	id1=perm[id];
	ids=id_array[id1].sire;
	idd=id_array[id1].dam;
	if(ids) (void)fprintf(fptr," father='%d'",ped_recode1[ids-1]);
	if(idd) (void)fprintf(fptr," mother='%d'",ped_recode1[idd-1]);
	if(id_array[id1].sex) (void)fprintf(fptr," sex='%d'",id_array[id1].sex);
	if(Affected && id_array[id1].affected) (void)fprintf(fptr," aff='%d'",id_array[id1].affected);
	if(Proband && id_array[id1].proband) (void)fprintf(fptr," proband='%d'",id_array[id1].proband);
	if(n_genetic_groups>1 && id_array[id1].group) (void)fprintf(fptr," group='%d'",id_array[id1].group);
	if(ped_recode[id1]) {
	  (void)fputs(" orig='",fptr);
	  print_orig_id1(fptr,id1+1,0);
	  (void)fputc('\'',fptr);
	}
	(void)fputs("/>\n",fptr);
      }
      if(k1>1) (void)fputs("  </comp>\n",fptr);
    }
    (void)fputs(" </family>\n",fptr);
  }
  if(family_recode1) free(family_recode1);
  free(perm);
  free(comp_to_family);
  /* System var information */
  j=0;
  if(traitlocus) {
    (void)fprintf(fptr," <group name='%s'>\n  <var name='%s",labels[3],traitlocus->var->name);
    if(traitlocus->var->vtype&ST_ARRAY) (void)fprintf(fptr,"(%d)",traitlocus->index);
    (void)fputs("' type='loki_qtl'/>\n",fptr);
    j=1;
  }
  if(n_genetic_groups>1) {
    if(!j++) (void)fprintf(fptr," <group name='%s'>\n",labels[3]);
    (void)fputs("  <var name='__GROUP__' type='group'>\n   <levels>",fptr);
    if(group_recode[0]->type==STRING) {
      escape_space(group_recode[0]->data.string,fptr);
      for(k=1;k<n_genetic_groups;k++) {
	fputc(' ',fptr);
	escape_space(group_recode[k]->data.string,fptr);
      }
    } else {
      fprintf(fptr,"%ld",group_recode[0]->data.value);
      for(k=1;k<n_genetic_groups;k++) fprintf(fptr," %ld",group_recode[k]->data.value);
    }
    fputs("</levels>\n  </var>\n",fptr);
  }
  model=Models;
  while(model) {
    mlist=model->model_list;
    while(mlist) {
      for(i=0;i<mlist->nvar;i++) if(mlist->element[i]->type&ST_ID) {
	elem=mlist->element[i];
	sd=elem->arg.var->data;
	if(!j++) (void)fprintf(fptr," <group name='%s'>\n",labels[3]);
	(void)fprintf(fptr,"  <var name='%s",sd->name);
	if(sd->vtype&ST_ARRAY) {
	  if(elem->type&ST_MARKER) k=markers[elem->index].index;
	  else k=elem->oindex;
	  (void)fprintf(fptr,"(%d)",k);
	}
	(void)fputs("' type='ind'/>\n",fptr);
	break;
      }
      if(i<mlist->nvar) break;
      mlist=mlist->next;
    }
    if(mlist) break;
    model=model->next;
  }
  for(i=0;i<n_markers;i++) if(markers[i].n_sub_elements) {
    if(!j++) (void)fprintf(fptr," <group name='%s'>\n",labels[3]);
    (void)fputs("  <var name=",fptr);
    if(markers[i].var->vtype&ST_ARRAY) fprintf(fptr,"'%s(%d)'",markers[i].var->name,markers[i].index);
    else fprintf(fptr,"'%s'",markers[i].var->name);
    (void)fputs(" type='super_locus'/>\n",fptr);
  }
  if(j) fputs(" </group>\n",fptr);
  /* Marker information */
  (void)fprintf(fptr," <group name='%s' type='microsat'>\n",labels[0]);
  for(i=0;i<n_markers;i++) if(markers[i].allele_trans) {
    (void)fputs("  <var name=",fptr);
    if(markers[i].var->vtype&ST_ARRAY) fprintf(fptr,"'%s(%d)'",markers[i].var->name,markers[i].index);
    else fprintf(fptr,"'%s'",markers[i].var->name);
    if(markers[i].element->type&ST_LUMPED) {
      for(j=0;j<n_markers;j++)  if(markers[i].element->arg.element==markers[j].element) {
	if(markers[j].var->vtype&ST_ARRAY) fprintf(fptr," parent='%s(%d)'",markers[j].var->name,markers[j].index);
	else fprintf(fptr," parent='%s'",markers[j].var->name);
	break;
      }
      assert(j<n_markers);
    }
    (void)fputs(">\n",fptr);
    k=markers[i].element->n_levels;
    if(k) {
      fputs("   <levels>",fptr);
      if(factor_recode[n_factors+i][0]->type==STRING) {
	escape_space(factor_recode[n_factors+i][0]->data.string,fptr);
	for(j=1;j<k;j++) {
	  fputc(' ',fptr);
	  escape_space(factor_recode[n_factors+i][j]->data.string,fptr);
	}
      } else {
	fprintf(fptr,"%ld",factor_recode[n_factors+i][0]->data.value);
	for(j=1;j<k;j++) fprintf(fptr," %ld",factor_recode[n_factors+i][j]->data.value);
      }
      fputs("</levels>\n",fptr);
    }
    (void)fputs("  </var>\n",fptr);
  }
  (void)fputs(" </group>\n",fptr);
  /* Phenotype information */
  if(n_id_records+n_nonid_records) {
    elem_list=lk_malloc(sizeof(void *)*(n_id_records+n_nonid_records));
    elem_list1=lk_malloc(sizeof(int)*(n_id_records+n_nonid_records));
    id_count=nonid_count=0;
    tm=ST_MODEL|ST_TRAIT;
    for(i=0;i<n_id_records;i++) if(id_elements[i]->type&tm) {
      elem_list1[id_count]=i;
      elem_list[id_count++]=id_elements[i];
    }
    j=id_count;
    for(i=0;i<n_nonid_records;i++) if(nonid_elements[i]->type&tm) {
      elem_list1[j]=i+n_id_records;
      elem_list[j++]=nonid_elements[i];
    }
    for(id=0;id<ped_size;id++) if(ped_recode1[id]) {
      if(id_array[id].nrec>1 && id_array[id].data1) break;
    }
    if(id<ped_size) nonid_count=j-id_count;
    else id_count=j;
    if(id_count) write_phen_var(id_count,elem_list,labels[1],fptr);
    if(nonid_count) write_phen_var(nonid_count,elem_list+id_count,labels[2],fptr);
  }
  /* Linkage information */
  /* Map information */
  linkp=links;
  link=0;
  while(linkp) {
    link++;
    (void)fputs(" <link",fptr);
    if(linkp->name) (void)fprintf(fptr," name='%s'",linkp->name);
    switch(linkp->type) {
    case LINK_AUTO:
      break;
    case LINK_X:
      fputs(" type='x'",fptr);
      break;
    case LINK_Y:
      fputs(" type='y'",fptr);
      break;
    default:
      ABT_FUNC("Unknown linkage type\n");
    }
    fputs(">\n",fptr);
    for(i=0;i<n_markers;i++) if(markers[i].link==link) {
      if(markers[i].element->type&ST_LUMPED) continue;
      if(markers[i].allele_trans || markers[i].n_sub_elements) {
	fputs("  <locus name='",fptr);
	if(markers[i].var->vtype&ST_ARRAY) fprintf(fptr,"%s(%d)",markers[i].var->name,markers[i].index);
	else fputs(markers[i].var->name,fptr);
	fputc('\'',fptr);
	for(j=0;j<3;j++) if(markers[i].pos_set[j]) fprintf(fptr," %s='%g'",map_str[j],markers[i].pos[j]);
	fputs("/>\n",fptr);
      }
    }
    (void)fputs(" </link>\n",fptr);
    linkp=linkp->next;
  }
  /* Genotype data */
  (void)fprintf(fptr," <data group='%s'>\n",labels[0]);
  for(i=0;i<ped_size;i++) if(ped_recode1[i] && id_array[i].haplo[0]) {
    for(j=0;j<n_markers;j++) {
      k=id_array[i].haplo[0][j];
      k1=id_array[i].haplo[1][j];
      if(k&&k1) break;
    }
    if(j==n_markers) continue;
    comp=id_array[i].component-1;
    fprintf(fptr,"   <rec id='%d'>",ped_recode1[i]);
    for(j=0;j<n_markers;j++) if(markers[j].allele_trans) {
      trans=markers[j].allele_trans[comp];
      k=id_array[i].haplo[0][j];
      k1=id_array[i].haplo[1][j];
      if(k || k1) {
	if(k) k=trans[k-1]+1;
	if(k1) k1=trans[k1-1]+1;
	if(k<k1) {
	  k2=k;
	  k=k1;
	  k1=k2;
	}
	k=k*(k+1)/2+k1;
      }
      if(j) fputc(',',fptr);
      if(k) fprintf(fptr,"%d",k);
    }
    fputs("</rec>\n",fptr);
  }
  (void)fputs(" </data>\n",fptr);
  /* Phenotype data */
  if(id_count) write_phen_data(id_count,0,elem_list,elem_list1,labels[1],fptr);
  if(nonid_count) write_phen_data(nonid_count,1,elem_list+id_count,elem_list1+id_count,labels[2],fptr);
  /* Model information */
  model=Models;
  if(model) {
    m_id=model->next?1:0;
    while(model) {
      if(m_id) fprintf(fptr," <model id='%d'>\n",m_id++);
      else fputs(" <model>\n",fptr);
      fprintf(fptr,"  <trait>%s",model->trait->name);
      if(model->trait->vtype&ST_ARRAY) fprintf(fptr,"(%d)",model->index);
      fputs("</trait>\n",fptr);
      mlist=model->model_list;
      while(mlist) {
	k=0;
	for(j=0;j<mlist->nvar;j++) {
	  elem=mlist->element[j];
	  sd=elem->arg.var->data;
	  if(elem->type&(ST_ID|ST_SIRE|ST_DAM|ST_RANDOM)) k=1;
	}
	fputs("  <term",fptr);
	if(k) fputs(" random='yes'",fptr);
	fputc('>',fptr);
	for(j=0;j<mlist->nvar;j++) {
	  if(j) fputc('*',fptr);
	  elem=mlist->element[j];
	  sd=elem->arg.var->data;
	  (void)fputs(sd->name,fptr);
	  if(sd->vtype&ST_ARRAY) {
	    if(elem->type&ST_MARKER) k=markers[elem->index].index;
	    else k=elem->oindex;
	    (void)fprintf(fptr,"(%d)",k);
	  }
	}
	mlist=mlist->next;
	fputs("</term>\n",fptr);
      }
      fputs(" </model>\n",fptr);
      model=model->next;
    }
  }
  (void)fputs("</loki>\n",fptr); /* Closing tag */
  if(fclose(fptr)) DataFileError(fname);
  if(elem_list) free(elem_list);
  if(elem_list1) free(elem_list1);
  free(fname);
  if(Filter) do i=wait(&j); while(i>0);
}
