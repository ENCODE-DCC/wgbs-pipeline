/****************************************************************************
 *                                                                          *
 *     Loki - Programs for genetic analysis of complex traits using MCMC    *
 *                                                                          *
 *                   Simon Heath - MSKCC                                    *
 *                                                                          *
 *                       August 2000                                        *
 *                                                                          *
 * kinship.c:                                                               *
 *                                                                          *
 * Calculate kinship coefficients                                           *
 *                                                                          *
 * Copyright (C) Simon C. Heath 1997, 2000, 2002                            *
 * This is free software.  You can distribute it and/or modify it           *
 * under the terms of the Modified BSD license, see the file COPYING        *
 *                                                                          *
 ****************************************************************************/

#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "utils.h"
#include "lk_malloc.h"
#include "loki.h"
#include "loki_utils.h"
#include "kinship.h"

/* Simple kinship coefficient */
double kinship(const int a, const int b,const struct Id_Record *id_array)
{
  if(!(a&&b)) return 0.0;
  if (a==b) return 0.5+0.5*kinship(id_array[a-1].sire,id_array[a-1].dam,id_array);
  if(!id_array[a-1].sire) {
    if(b<a) return 0.0;
    if(!id_array[b-1].sire) return 0.0;
    return (kinship(a,id_array[b-1].sire,id_array)+kinship(a,id_array[b-1].dam,id_array))*0.5;
  }
  if(!id_array[b-1].sire) {
    if(a<b) return 0.0;
    return (kinship(b,id_array[a-1].sire,id_array)+kinship(b,id_array[a-1].dam,id_array))*0.5;
  }
  if(a<b) return (kinship(a,id_array[b-1].sire,id_array)+kinship(a,id_array[b-1].dam,id_array))*0.5;
  return (kinship(b,id_array[a-1].sire,id_array)+kinship(b,id_array[a-1].dam,id_array))*0.5;
}

static struct kin **kin;
static struct Id_Record *id_array;
static int ped_size;

double calc_gen(int n,int *g,int *b,int xfl)
{
  int i,j,k,k1,bk,ids,idd,t[2],*b1,*g1,n1,sex;
  double p,z,z1;
  static int *tmp[2];
	
  /* Clean up if requested */
  if(!n) {
    if(tmp[0]) {
      free(tmp[0]);
      tmp[0]=0;
    }
    return 0.0;
  }
  if(!tmp[0]) {
    tmp[0]=lk_malloc(sizeof(int)*2*ped_size);
    tmp[1]=tmp[0]+ped_size;
  }
  /* Boundary condition 1 */
  /* Check for an individual being present in >2 blocks (or >1 if male and x-linked) */
  for(i=0;i<n;i++) {
    k=g[i];
    tmp[0][k]=tmp[1][k]=-1;
  }
  for(i=0;i<n;i++) {
    k=g[i];
    for(j=0;j<2;j++) {
      if(tmp[j][k]==b[i]) break;
      if(tmp[j][k]<0) {
	tmp[j][k]=b[i];
	break;
      }
    }
    if(xfl) {
      sex=id_array[k].sex;
      k1=sex;
    } else k1=2;
    if(j==k1) return 0.0;
  }
  /* Boundary condition 2 */
  /* Check for >1 founder in 1 block */
  j=bk=0;
  while(j<n) {
    for(k=-1,i=0;i<n;i++) if(b[i]==bk) {
      if(!id_array[g[i]].sire) {
	if(k<0) k=g[i];
	else if(k!=g[i]) break;
      }
      j++;
    }
    if(i<n) return 0.0;
    bk++;
  }
  /* If sex-linked, remove multiple copies of male genes */
  if(xfl) {
    for(i=0;i<n;i++)  {
      k=g[i];
      if(id_array[k].sex==1) {
	for(j=i+1;j<n;j++) if(g[j]==k) break;
	if(j<n) break;
      }
    }
    if(i<n) {
      g1=lk_malloc(sizeof(int)*n*2);
      b1=g1+n;
      for(i=0;i<n;i++) b1[i]=0;
      for(i=j=0;i<n;i++) if(!b1[i]) {
	k=g[i];
	g1[j]=k;
	b1[j++]=b[i];
	if(id_array[k].sex==1) {
	  for(k1=i+1;k1<n;k1++) if(g[k1]==k) b1[k1]=-1;
	}
      }
      p=calc_gen(j,g1,b1,xfl);
      free(g1);
      return p;
    }
  }
  /* Boundary condition 3 */
  /* Check if only founder genes are involved */
  for(i=0;i<n;i++) {
    k=g[i];
    if(id_array[k].sire) break;
  }
  if(i==n) {
    for(i=0;i<n;i++) {
      k=g[i];
      tmp[0][k]=0;
    }
    for(j=i=0;i<n;i++) {
      k=g[i];
      if(!tmp[0][k]) {
	j++;
	tmp[0][k]=1;
      }
    }
    k=1<<(n-j);
    p=1.0/(double)k;
    return p;
  }
  /* Not a boundary condition - use recurrence rules */
  /* Pick highest numbered non-founder gene */
  for(j=i=0;i<n;i++) if(id_array[g[i]].sire && g[i]>j) j=g[i];
  /* Find out how many times it is present */
  t[0]=t[1]=-1;
  for(k=i=0;i<n;i++) if(g[i]==j) {
    k++;
    if(t[0]<0) t[0]=b[i];
    else if(t[0]!=b[i]) t[1]=b[i];
  }
  if(xfl && id_array[j].sex==1 && k>1) {
    fputs("Internal error! Multiple male genes at x-linked locus\n",stderr);
    exit(EXIT_FAILURE);
  }
  ids=id_array[j].sire-1;
  idd=id_array[j].dam-1;
  if(k==1) {
    g1=lk_malloc(sizeof(int)*n);
    if(xfl && id_array[j].sex==1) {
      for(i=0;i<n;i++) {
	if(g[i]==j) g1[i]=idd;
	else g1[i]=g[i];
      }
      p=calc_gen(n,g1,b,xfl);
    } else {
      for(i=0;i<n;i++) {
	if(g[i]==j) g1[i]=ids;
	else g1[i]=g[i];
      }
      z1=calc_gen(n,g1,b,xfl);
      p=.5*z1;
      for(i=0;i<n;i++) if(g[i]==j) {
	g1[i]=idd;
	break;
      }
      z1=calc_gen(n,g1,b,xfl);
      p+=.5*z1;
    }
    free(g1);
  } else {
    n1=n-k+2;
    g1=lk_malloc(sizeof(int)*n1*2);
    b1=g1+n1;
    if(t[1]>=0) {
      g1[0]=ids;
      g1[1]=idd;
      b1[0]=t[0];
      b1[1]=t[1];
      for(k1=2,i=0;i<n;i++) if(g[i]!=j) {
	g1[k1]=g[i];
	b1[k1++]=b[i];
      }
      z=1.0/(double)(1<<k);
      p=calc_gen(n1,g1,b1,xfl);
      b1[0]=t[1];
      b1[1]=t[0];
      z1=calc_gen(n1,g1,b1,xfl);
      p=z*(p+z1);
    } else {
      g1[0]=ids;
      g1[n1-1]=idd;
      b1[0]=b1[n1-1]=t[0];
      z=1.0/(double)(1<<k);
      for(k1=1,i=0;i<n;i++) if(g[i]!=j) {
	g1[k1]=g[i];
	b1[k1++]=b[i];
      }
      if(id_array[ids].sire || id_array[idd].sire) {
	z1=calc_gen(n1,g1,b1,xfl);
	p=(1.0-2.0*z)*z1;
      } else p=0.0;
      z1=calc_gen(n1-1,g1,b1,xfl);
      p+=z*z1;
      g1[0]=idd;
      z1=calc_gen(n1-1,g1,b1,xfl);
      p+=z*z1;
    }
    free(g1);
  }
  return p;
}
				
void calc_kin(int *inbr,struct loki *loki)
{
  int i,id,id1,j,k,cs,comp,xx[4],xx1[4],k1,k2,n_comp,*comp_size;
  struct kin *kk1;
  double w[9],z;
  int states[][4]={{0,0,0,0},{0,0,1,1},{0,0,0,1},{0,0,1,2},
		   {0,1,0,0},{0,1,2,2},{0,1,0,1},{0,1,0,2},{0,1,2,3}};
  int coeff[]={1,1,2,1,2,1,2,4,1};
  FILE *fptr;
	
  message(INFO_MSG,"Calculating detailed kinship coefficients\n");
  ped_size=loki->pedigree->ped_size;
  id_array=loki->pedigree->id_array;
  n_comp=loki->pedigree->n_comp;
  comp_size=loki->pedigree->comp_size;
  fptr=fopen("loki_ibs.kin","w");
  kin=lk_malloc(sizeof(void *)*n_comp);
  for(j=i=0;i<n_comp;i++) {
    cs=comp_size[i];
    j+=cs*(cs+1)/2;
  }
  kin[0]=lk_malloc(sizeof(struct kin)*j);
  for(i=1;i<n_comp;i++) {
    cs=comp_size[i-1];
    j=cs*(cs+1)/2;
    kin[i]=kin[i-1]+j;
  }
  for(i=comp=0;comp<n_comp;comp++) {
    cs=comp_size[comp];
    for(j=1;j<cs;j++) {
      id=j+i;
      for(k=0;k<j;k++) {
	id1=k+i;
	xx[0]=xx[1]=id;
	xx[2]=xx[3]=id1;
	if(inbr[id]>0.0 || inbr[id1]>0.0) {
	  z=0.0;
	  for(k1=0;k1<9;k1++) {
	    for(k2=0;k2<4;k2++) xx1[k2]=states[k1][k2];
	    w[k1]=(double)coeff[k1]*calc_gen(4,xx,xx1,0);
	    z+=w[k1];
	  }
	  if(fabs(z-1.0)>1.0e-16) printf("Warning: z=%g\n",z);
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
	kk1=kin[comp]+j*(j-1)/2+k;
	for(k1=0;k1<9;k1++) {
	  kk1->p[k1]=w[k1];
	}
	if(fptr) {
	  print_orig_id(fptr,id+1);
	  (void)fputc(' ',fptr);
	  print_orig_id(fptr,id1+1);
	  for(k1=0;k1<9;k1++) fprintf(fptr," %-6g",w[k1]);
	  (void)fputc('\n',fptr);
	}
      }
    }
    i+=cs;
  }
  free(kin[0]);
  free(kin);
  if(fptr) fclose(fptr);
}
