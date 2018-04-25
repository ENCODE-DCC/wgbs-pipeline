/****************************************************************************
 *                                                                          *
 *     Loki - Programs for genetic analysis of complex traits using MCMC    *
 *                                                                          *
 *                   Simon Heath - CNG, Evry                                *
 *                                                                          *
 *                        October 2005                                      *
 *                                                                          *
 * loki_ibs_check.c:                                                        *
 *                                                                          *
 * Routines for checking pedigree IBS relationships                         *
 *                                                                          *
 * Copyright (C) Simon C. Heath 2005                                        *
 * This is free software.  You can distribute it and/or modify it           *
 * under the terms of the Modified BSD license, see the file COPYING        *
 *                                                                          *
 ****************************************************************************/

#include <config.h>
#include <stdlib.h>
#ifdef USE_DMALLOC
#include <dmalloc.h>
#endif
#include <math.h>
#include <stdio.h>

#include "lk_malloc.h"
#include "utils.h"
#include "loki.h"
#include "loki_peel.h"
#include "loki_utils.h"
#include "seg_pen.h"
#include "loki_ibd.h"
#include "ibs_check.h"
#include "kinship.h"

static char *rtd[]={"PO","HS","FS","1CS","DCS","2CS","--","AV","UN"};
static int pflag;

static void print_id(FILE *fptr,int i,struct loki *loki)
{
  if(loki->pedigree->family_id) (void)print_orig_family(fptr,i,0);
  fputc('\t',fptr);
  (void)print_orig_id1(fptr,i);
}

static void print_res(FILE *fptr,struct result *res,struct loki *loki,int self_flag)
{
  double z,zz,sd;

  print_id(fptr,res->id+1,loki);
  z=res->ibs_obs;
  zz=res->ibs_exp;
  sd=res->ibs_sd;
  if(!self_flag) {
    fputc('\t',fptr);
    print_id(fptr,res->id1+1,loki);
    fprintf(fptr,"\t%.2f\t%g\t%.2f\t%.2f\t%d\t%d\t%d\t%d\t%s\n",zz,z,z-zz,(z-zz)/sd,res->n,res->n_disc,res->n1,res->n2,rtd[res->rel]);
  } else fprintf(fptr,"\t%.2f\t%g\t%.2f\t%.2f\t%d\n",zz,z,z-zz,(z-zz)/sd,res->n);
}

static void print_res_head(FILE *fptr,int self_flag)  
{
  if(self_flag) {
    if(pflag) fputs("Ped\t",fptr);
    fputs("Id",fptr);
		
  } else {
    if(pflag) fputs("Ped 1\tId 1\tPed 2\tId 2",fptr);
    else fputs("Id 1\tId 2",fptr);
  }
  fputs("\tExp\tObs\tExp-Obs\tZ\tn",fptr);
  if(!self_flag) fputs("\tn_disc\tn1\tn2\tRel",fptr);
  fputc('\n',fptr);
}

static int cmp_res(const void *s1,const void *s2)
{
  const struct result *res1,*res2;
  double x1,x2;
	
  res1=s1;
  res2=s2;
  x1=(res1->ibs_obs-res1->ibs_exp)/res1->ibs_sd;
  x2=(res2->ibs_obs-res2->ibs_exp)/res2->ibs_sd;
  if(x1<x2) return -1;
  if(x1>x2) return 1;
  return 0;
}

struct sample {
  int id;
  int rec;
  int n_longs;
  unsigned long *fnd;
};

void ibs_check(int no_unrel,struct loki *loki)
{
  int *n_longs,*inbr,n_data,ped_size,nl,n_self,n_rel,n_diff_fam,n_genetic_groups;
  unsigned long *founders,*tp,*tp1;
  int xx[4],xx1[4],np,n_pairs,rel_idx=0,max_all,n1,n2,n_comp,*comp_size,rec,rec1;
  int **ibs_tab,kk,old_id,n_pairs_fam,mod_val,n_disc;
  int i,i1,i2,j,k,k1,k2,k3,id,id1,comp,cs,g0,g1,g2,n,nlong,rel_flag=0,n_all;
  double zz,z,sd,*ff,rz[7],**lcf1,**lcf2,**lcf3,*cf1,*cf2,*cf3,rzz,rzz1,zk1,zk2,w[9],ww[5];
  double tz1=0,tz2=0;
  double zz_un,zz_po,zz_fs,zz_hs,zz_cs,*inb,zk1a,zk2a;
  int states[][4]={{0,0,0,0},{0,0,1,1},{0,0,0,1},{0,0,1,2},
		   {0,1,0,0},{0,1,2,2},{0,1,0,1},{0,1,0,2},{0,1,2,3}};
  int coeff[]={1,1,2,1,2,1,2,4,1};
  double rt[][3]={{0,1.0,0},{0,.5,.5},{.25,.5,.25},{0,.25,.75},
		  {.0625,.3125,.5625},{0,.0625,.9375}};
  struct result *results=0,*results_self=0,*results_same_fam=0,*results_diff_fam=0,*res;
  struct sample *sample_list;
  FILE *fptr;
  struct Marker **mklist,*mark;
  struct Id_Record *id_array;

  get_founder_params(&founders,&n_longs,0,&inbr,loki);
  calc_kin(inbr,loki);
  pflag=loki->pedigree->family_id?1:0;
  ped_size=loki->pedigree->ped_size;
  id_array=loki->pedigree->id_array;
  n_genetic_groups=loki->pedigree->n_genetic_groups;
  for(n_data=i=0;i<ped_size;i++) n_data+=id_array[i].n_gt_sets;
  if(!n_data || !loki->markers->n_markers) {
    fputs("ibs_check(): No data!\n",stderr);
  }
  sample_list=lk_malloc(sizeof(struct sample)*n_data);
  n_comp=loki->pedigree->n_comp;
  comp_size=loki->pedigree->comp_size;
  tp=founders;
  for(n_data=i=comp=0;comp<n_comp;comp++) {
    nlong=n_longs[comp];
    cs=comp_size[comp];
    for(j=0;j<cs;j++,tp+=nlong,i++) {
      for(k=0;k<id_array[i].n_gt_sets;k++) {
	sample_list[n_data].id=i;
	sample_list[n_data].n_longs=n_longs[comp];
	sample_list[n_data].fnd=tp;
	sample_list[n_data++].rec=id_array[i].gt_idx[k];
      }
    }
  }
  mklist=lk_malloc(sizeof(void *)*loki->markers->n_markers);
  for(nl=i=0;i<loki->markers->n_markers;i++) {
    j=loki->markers->marker[i].locus.link_group;
    if(loki->markers->linkage[j].type==LINK_AUTO) {
      mark=loki->markers->marker+i;
      /* Check for unobserved or monomorphic markers */
      k1=0;
      for(k=0;k<n_data;k++) {
	if((k2=mark->orig_gt[sample_list[k].rec])) {
	  if(!k1) k1=k2;
	  else if(k1!=k2) break;
	}
      }
      if(k==n_data && k1) { /* Observed, but only 1 genotype */
	/* Check if this is a homozygous or heterozygous genotype */
	g1=(int)(sqrt(.25+2.0*(double)k1)-.49999);
	g2=k1-(g1*(g1+1)/2);
	if(g1!=g2) k=0;
      }
      if(k<n_data) mklist[nl++]=mark;
    }
  }
  for(i=k=0;i<n_data;i++) {
    id=sample_list[i].rec;
    for(j=0;j<nl;j++) if(mklist[j]->orig_gt[id]) break;
    if(j<nl) sample_list[k++]=sample_list[i];
  }
  n_data=k;
  if(!n_data) return;
  n_pairs=n_data*(n_data-1)/2;
  lcf1=lk_malloc(sizeof(void *)*nl*3);
  lcf1[0]=lk_malloc(sizeof(double)*nl*15);
  lcf2=lcf1+nl;
  lcf3=lcf2+nl;
  for(i=1;i<nl;i++) lcf1[i]=lcf1[i-1]+8;
  lcf2[0]=lcf1[nl-1]+8;
  lcf3[0]=lcf2[0]+4;
  max_all=mklist[0]->locus.n_alleles;
  for(i=1;i<nl;i++) {
    if(mklist[i]->locus.n_alleles>max_all) max_all=mklist[i]->locus.n_alleles;
    lcf2[i]=lcf2[i-1]+7;
    lcf3[i]=lcf2[i]+4;
  }
  k=(max_all+1)*(max_all+2)/2;
  if(!(ibs_tab=malloc(sizeof(void *)*k))) exit(EXIT_FAILURE);
  if(!(ibs_tab[0]=malloc(sizeof(int)*k*k))) exit(EXIT_FAILURE);
  for(i=1;i<k;i++) ibs_tab[i]=ibs_tab[i-1]+k;
  for(i=k=0;k<=max_all;k++) {
    for(k1=0;k1<=k;k1++,i++) {
      for(j=k2=0;k2<=max_all;k2++) {
	for(k3=0;k3<=k2;k3++,j++) {
	  kk=0;
	  if(k && k==k2) {
	    if(k1 && k1==k3) kk=2;
	    else kk=1;
	  } else if((k && k==k3) || (k1 && (k1==k2 || k1==k3))) kk=1;
	  ibs_tab[i][j]=kk;
	}
      }
    }
  }
  inb=lk_malloc(sizeof(double)*ped_size);
  fptr=fopen("ibs.inbr","w");
  old_id=n_pairs_fam=k=0;
  for(i=0;i<n_data;i++) {
    id=sample_list[i].id;
    if(id_array[id].fam_code!=id_array[old_id].fam_code) {
      n_pairs_fam+=k*(k-1)/2;
      k=0;
      old_id=id;
    }
    k++;
    if(inbr[id]) {
      xx[0]=xx[1]=id;
      xx1[0]=xx1[1]=0;
      z=calc_gen(2,xx,xx1,0);
      inb[id]=2.0*(z-.5);
    } else inb[id]=0.0;
    if(fptr) {
      print_id(fptr,id+1,loki);
      fprintf(fptr," %g",inb[id]);
      fputc('\n',fptr);
    }
  }
  n_pairs_fam+=k*(k-1)/2;
  if(fptr) fclose(fptr);
  zz_un=zz_po=zz_fs=zz_hs=zz_cs=0.0;
  ww[0]=ww[1]=ww[2]=ww[3]=ww[4]=0.0;
  n=0;
  for(k=0;k<nl;k++) {
    cf1=lcf1[k];
    cf2=lcf2[k];
    cf3=lcf3[k];
    if(n_genetic_groups==1) ff=mklist[k]->locus.freq[0];
    else {
      ff=0;
      ABT_FUNC("Not working with multiple genetic groups (yet)\n");
    }
    n_all=mklist[k]->locus.n_alleles;
    for(k1=0;k1<7;k1++) rz[k1]=0.0;
    for(k1=0;k1<n_all;k1++) {
      zk1=ff[k1];
      rz[0]+=(rzz=zk1*zk1);
      rz[1]+=(rzz1=rzz*zk1);
      rz[2]+=rzz1*zk1;
      for(k2=0;k2<n_all;k2++) if(k1!=k2) {
	zk2=ff[k2];
	rz[3]+=rzz*zk2;
	rz[4]+=rzz1*zk2;
	rz[5]+=rzz*zk2*zk2;
	for(k3=0;k3<n_all;k3++) if(k3!=k1 && k3!=k2) rz[6]+=rzz*zk2*ff[k3];
      }
    }
    cf1[0]=1.0-rz[0];
    cf1[1]=2.0*rz[3];
    cf1[2]=4.0*(rz[4]+rz[6]);
    cf2[0]=rz[0];
    cf2[1]=rz[1];
    cf2[2]=rz[2]+2.0*rz[5];
    cf1[3]=cf1[2]+2.0*cf2[2]; /* Unrelated  - mean */
    cf2[3]=cf1[2]+4.0*cf2[2]; /* Unrelated - variance */
    cf1[4]=cf1[0]+2.0*cf2[0]; /* Parent offspring mean */
    cf1[5]=.5*cf1[0]+.25*cf1[2]+.5+cf2[0]+.5*cf2[2]; /* Full sibs mean */
    cf1[6]=.5*(cf1[0]+cf1[2])+cf2[0]+cf2[2]; /* Half sibs mean */
    cf1[7]=.25*cf1[0]+.75*cf1[2]+.5*cf2[0]+1.5*cf2[2]; /* 1st cousins mean */
    cf3[0]=cf1[0]+2.0*cf2[0];
    cf3[1]=cf1[1]+2.0*cf2[1];
    cf3[2]=.5*cf1[2]+cf2[2];
    if(n_all>1) {
      zz_un+=cf1[3];
      zz_po+=cf1[4];
      zz_fs+=cf1[5];
      zz_hs+=cf1[6];
      zz_cs+=cf1[7];
      n++;
    }
  }
  if(n) {
    zz_un/=(double)n;
    zz_po/=(double)n;
    zz_fs/=(double)n;
    zz_hs/=(double)n;
    zz_cs/=(double)n;
  }
  printf("%d %g %g %g %g %g\n",n,zz_un,zz_po,zz_fs,zz_hs,zz_cs);
  k=no_unrel?n_pairs_fam:n_pairs;
  mod_val=(int)(1.0+(double)k/100.0);
  if(mod_val<1000) mod_val=1000;
  k1=k+n_data;
  if(k1) {
    results=lk_malloc(sizeof(struct result)*k1);
    results_self=results+k;
    results_diff_fam=results+n_pairs_fam;
    results_same_fam=results_diff_fam; 
  }
  n_rel=n_self=n_diff_fam=0;
  old_id=i2=0;
  for(np=i=0;i<n_data;i++) {
    id=sample_list[i].id;
    tp=sample_list[i].fnd;
    nlong=sample_list[i].n_longs;
    rec=sample_list[i].rec;
    zk1=inb[id];
    zk1a=1.0-zk1;
    if(no_unrel && id_array[id].fam_code!=id_array[old_id].fam_code) {
      old_id=id;
      i2=i;
    }
    comp=id_array[id].comp;
    for(i1=i2;i1<=i;i1++) {
      id1=sample_list[i1].id;
      if(id!=id1) {
	np++;
	if(!(np%mod_val)) {
	  fputs("At: ",stdout);
	  k=no_unrel?n_pairs_fam:n_pairs;
	  printf(" (%d/%d) - %2.0f%%       \r",np,k,100.0*(double)np/(double)k);
	  fflush(stdout);
	}
	if(id_array[id1].comp!=comp) rel_flag=0;
	else {
	  tp1=sample_list[i1].fnd;
	  for(k=0;k<nlong;k++) if(tp[k]&tp1[k]) break;
	  rel_flag=k<nlong?1:0;
	}
	if(rel_flag) {
	  rel_idx=6;
	  xx[0]=xx[1]=id;
	  xx[2]=xx[3]=id1;
	  if(inbr[id]||inbr[id1]) {
	    rzz=0.0;
	    for(k1=0;k1<9;k1++) {
	      for(k2=0;k2<4;k2++) xx1[k2]=states[k1][k2];
	      w[k1]=(double)coeff[k1]*calc_gen(4,xx,xx1,0);
	      rzz+=w[k1];
	    }
	    w[0]=w[0]-.5*(w[2]+w[4]-w[6])+.25*w[7];
	    w[1]=w[1]-.5*(w[2]+w[4]-w[6])-w[3]-w[5]+.75*w[7]+w[8];
	    w[2]=2.0*(w[2]-w[6])-w[7];
	    w[3]=2.0*(w[3]-w[8])-w[7];
	    w[4]=2.0*(w[4]-w[6])-w[7];
	    w[5]=2.0*(w[5]-w[8])-w[7];
	  } else {
	    for(k1=0;k1<6;k1++) w[k1]=0.0;
	    for(;k1<9;k1++) {
	      for(k2=0;k2<4;k2++) xx1[k2]=states[k1][k2];
	      w[k1]=(double)coeff[k1]*calc_gen(4,xx,xx1,0);
	    }
	  }
	  w[6]=4.0*w[6];
	  w[7]=4.0*w[7];
	  w[8]=4.0*w[8];
	  if(fptr) {
	    print_id(fptr,id+1,loki);
	    (void)fputc(' ',fptr);
	    print_id(fptr,id1+1,loki);
	  }
	  for(k1=0;k1<9;k1++) fprintf(fptr," %g",w[k1]);
	  (void)fputc('\n',fptr);
	  ww[0]=w[0]+w[6];
	  ww[1]=w[2]+w[4]+w[7];
	  ww[2]=w[1]+ww[1];
	  ww[3]=w[3]+w[5];
	  ww[4]=w[8];
	  for(k1=0;k1<6;k1++) {
	    for(k2=0;k2<3;k2++) if(fabs(w[6+k2]-rt[k1][k2])>1.0e-12) break;
	    if(k2==3) break;
	  }
	  if(k1==1) {
	    k1=7;
	    if(id_array[id].sire && id_array[id].sire==id_array[id1].sire) k1=1;
	    else if(id_array[id].dam && id_array[id].dam==id_array[id1].dam) k1=1;
	  }
	  rel_idx=k1;
	} else {
	  rel_idx=8;
	  zk2=inb[id1];
	  zk2a=1.0-zk2;
	  tz1=zk1a*zk2a;
	  tz2=zk2a*zk1+zk2*zk1a;
	  w[0]=zk1*zk2;
	}
      }
      z=zz=sd=0.0;
      n=n1=n2=n_disc=0;
      if(i==i1) {
	for(k=0;k<nl;k++) {
	  k1=mklist[k]->orig_gt[rec];
	  if(k1) {
	    g0=(int)(sqrt(.25+2.0*(double)k1)-.49999);
	    g1=k1-(g0*(g0+1)/2);
	    if(g0) {
	      if(g0==g1) z++;
	      cf2=lcf2[k];
	      zz+=zk1+(1.0-zk1)*cf2[0];
	      n++;
	    }
	  }
	}
      } else {
	rec1=sample_list[i1].rec;
	if(rel_flag) {
	  for(k=0;k<nl;k++) {
	    g0=mklist[k]->orig_gt[rec];
	    g2=mklist[k]->orig_gt[rec1];
	    if(g0) n1++;
	    if(g2) n2++;
	    if(g0&&g2) { /* Both typed for this locus */
	      z+=(double)ibs_tab[g0][g2];
	      cf1=lcf1[k];
	      cf2=lcf2[k];
	      rzz=ww[1]*cf1[0]+ww[3]*cf1[1]+ww[4]*cf1[2];
	      rzz1=ww[0]+ww[2]*cf2[0]+ww[3]*cf2[1]+ww[4]*cf2[2];
	      zz+=rzz+2.0*rzz1;
	      sd+=rzz+4.0*rzz1;
	      n++;
	      if(g0!=g2) n_disc++;
	    }
	  }
	} else {
	  for(k=0;k<nl;k++) {
	    g0=mklist[k]->orig_gt[rec];
	    g2=mklist[k]->orig_gt[rec1];
	    if(g0) n1++;
	    if(g2) n2++;
	    if(g0&&g2) { /* Both typed for this locus */
	      z+=(double)ibs_tab[g0][g2];
	      cf1=lcf1[k];
	      cf2=lcf2[k];
	      rzz=tz2*cf1[1]+tz1*cf1[2];
	      rzz1=w[0]*cf2[0]+tz2*cf2[1]+tz1*cf2[2];
	      zz+=rzz+2.0*rzz1;
	      sd+=rzz+4.0*rzz1;
	      n++;
	      if(g0!=g2) n_disc++;
	    }
	  }
	}
      }
      if(n>1) {
	rzz=(double)n;
	z/=rzz;
	zz/=rzz;
	if(i==i1) {
	  res=results_self+(n_self++);
	  if(zz==1.0) {
	    sd=1.0;
	  } else sd=sqrt(zz*(1.0-zz)/rzz);
	  rel_idx=-1;
	} else {
	  if(id==id1) {
	    sd=1.0;
	    zz=2.0;
	  } else sd=sqrt((sd-zz*zz*rzz)/(rzz*(rzz-1.0)));
	  if(!no_unrel && id_array[id].fam_code!=id_array[id1].fam_code) {
	    res=results_diff_fam+(n_diff_fam++);
	  } else if(!rel_flag) {
	    res=--results_same_fam;
	  } else {
	    res=results+(n_rel++);
	  }
	}
	res->id=id;
	res->id1=id1;
	res->ibs_obs=z;
	res->ibs_exp=zz;
	res->ibs_sd=sd;
	res->n=n;
	res->n1=n1;
	res->n2=n2;
	res->n_disc=n_disc;
	res->rel=rel_idx;
      }
    }
  }
  fclose(fptr);
  fputs("At: ",stdout);
  k=no_unrel?n_pairs_fam:n_pairs;
  printf(" (%d/%d) - %2.0f%%\n",np,k,100.0*(double)np/(double)k);
  k=results_diff_fam-results_same_fam;
  if(n_self) {
    fputs("Sorting results: self\r",stdout);
    fflush(stdout);
    qsort(results_self,n_self,sizeof(struct result),cmp_res);
  }
  if(n_rel) {
    fputs("Sorting results: related pairs\r",stdout);
    fflush(stdout);
    qsort(results,n_rel,sizeof(struct result),cmp_res);
  }
  if(k) {
    fputs("Sorting results: unrelated pairs (within families)\r",stdout);
    fflush(stdout);
    qsort(results_same_fam,k,sizeof(struct result),cmp_res); 
  }
  if(n_diff_fam) {
    fputs("Sorting results: unrelated pairs (between families)\r",stdout);
    fflush(stdout);
    qsort(results_diff_fam,n_diff_fam,sizeof(struct result),cmp_res);
  }
  fputc('\n',stdout);
  if(n_self) {
    fputs("Saving results: self\r",stdout);
    fflush(stdout);
    if(!(fptr=fopen("self.dat","w"))) exit(EXIT_FAILURE);
    print_res_head(fptr,1);
    for(i=0;i<n_self;i++) print_res(fptr,results_self+i,loki,1);
    fclose(fptr);
  }
  if(n_rel) {
    fputs("Saving results: related pairs\r",stdout);
    fflush(stdout);
    if(!(fptr=fopen("related.dat","w"))) exit(EXIT_FAILURE);
    print_res_head(fptr,0);
    for(i=0;i<n_rel;i++) print_res(fptr,results+i,loki,0);
    fclose(fptr);
  }
  if(k) {
    fputs("Saving results: unrelated pairs (within families)\r",stdout);
    fflush(stdout);
    if(!(fptr=fopen("within_fam.dat","w"))) exit(EXIT_FAILURE);
    print_res_head(fptr,0);
    for(i=0;i<k;i++) print_res(fptr,results_same_fam+i,loki,0);
    fclose(fptr);
  }
  if(n_diff_fam) {
    fputs("Saving results: unrelated pairs (between families)\r",stdout);
    fflush(stdout);
    qsort(results_diff_fam,n_diff_fam,sizeof(struct result),cmp_res);
    if(!(fptr=fopen("between_fam.dat","w"))) exit(EXIT_FAILURE);
    print_res_head(fptr,0);
    for(i=0;i<n_diff_fam;i++) print_res(fptr,results_diff_fam+i,loki,0);
    fclose(fptr);
  }
  fputc('\n',stdout);
  if(results) free(results);
  free(lcf1[0]);
  free(lcf1);
  free(ibs_tab[0]);
  free(ibs_tab);
  free(inb);
  free(mklist);
  free(sample_list);
}
