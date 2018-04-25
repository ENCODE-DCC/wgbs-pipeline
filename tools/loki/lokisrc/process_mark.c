/****************************************************************************
*                                                                          *
*     Loki - Programs for genetic analysis of complex traits using MCMC    *
*                                                                          *
*                     Simon Heath - CNG, Evry                              *
*                                                                          *
*                            May 2003                                      *
*                                                                          *
* process_mark.c:                                                          *
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
#include <assert.h>

#include "utils.h"
#include "loki.h"
#include "loki_peel.h"
#include "loki_utils.h"
#include "lk_malloc.h"
#include "gen_elim.h"

/* Add individuals genotype to list, reallocating array if required. */
int add_ind_gtype(gtype **gt,int k,int *j,const int *g)
{
	int i,g1[2],k1,k2,sx,del=0;
	int map[]={0,1,0,1,1,2,0,2,0};
	gtype *gg;
	
	gg=*gt;
	for(i=0;i<k;i++) {
		g1[X_MAT]=gg[i].mat;
		g1[X_PAT]=gg[i].pat;
		k2=0;
		for(sx=0;sx<2;sx++) {
			k1=0;
			if(g[sx]) {
				if(g1[sx]) {
					if(g[sx]==g1[sx]) k1=1;
				} else k1=3;
			} else k1=g1[sx]?2:1;
			if(!k1) break;
			k2=k2*3+k1-1;
		}
		if(sx==2) {
			k2=map[k2];
			if(k2==1) {
				if(!del) {
					gg[i].mat=g[0];
					gg[i].pat=g[1];
				} else gg[i].mat=-1;
				del=1;
			} else if(!k2) break;
		}
	}
	if(i==k) {
		if(del) {
			for(i=k1=k2=0;i<k;i++) if(gg[i].mat>=0) {
				gg[k1].mat=gg[i].mat;
				gg[k1++].pat=gg[i].pat;
			}
				k=k1;
		} else {
			if(k==*j) {
				(*j)<<=1;
				*gt=lk_realloc(*gt,sizeof(gtype)*(*j));
			}
			gg=*gt+k;
			gg->mat=g[0];
			gg->pat=g[1];
			k++;
		}
	}
	return k;
}

/* Assign initial genotype lists based on microsat marker data */
static int assign_gtypes(int ix,int *hap,gtype **gts,int *ngens)
{
	int k,a1,a2,fix=0;
	gtype *gt;
	
	if((k=hap[ix])) {
		a1=(int)(sqrt(.25+2.0*(double)k)-.49999);
		a2=k-(a1*(a1+1)/2);
		if(a1!=a2) {
			gt=lk_malloc(sizeof(gtype)*2);
			gt[0].mat=gt[1].pat=a1;
			gt[0].pat=gt[1].mat=a2;
			ngens[ix]=2;
		} else {
			gt=lk_malloc(sizeof(gtype));
			gt[0].mat=gt[0].pat=a1;
			ngens[ix]=1;
			fix=1;
		}
		gts[ix]=gt;
	} else ngens[ix]=0;
	return fix;
}

/* Assign initial genotype lists based on individual microsat/snp marker data */
static void get_family_gtypes(nuc_fam *fam,int *hap,gtype **gts,int *ngens)
{
	int i,sx;
	struct Id_Record *par[2],**kids,*kid;
	
	par[X_PAT]=fam->father;
	par[X_MAT]=fam->mother;
	for(sx=0;sx<2;sx++) if(par[sx] && !par[sx]->flag) {
		par[sx]->flag=assign_gtypes(par[sx]->idx,hap,gts,ngens)?-1:1;
	}
		i=-1;
	kids=fam->kids;
	while((kid=kids[++i])) {
		if(!kid->flag) {
			kid->flag=assign_gtypes(kid->idx,hap,gts,ngens)?-1:1;
		}
	}
}

/* Add parental genotypes to family set.
* Check for duplicates, and deal correctly with wildcards */
static int add_gt_famlist(fam_gtype **fgt,int k,int *sz,int *g)
{
	int i,j,k1,k2,sx,ng[2],match=0,g1[2],found=0;
	int op[]={1,1,2,1,1,0,2,0,2},deleted=0;
	fam_gtype *fg;
	
	fg=*fgt;
	for(j=sx=0;sx<2;sx++,j+=2) {
		for(i=0;i<2;i++) if(!g[j+i]) break;
		ng[sx]=i;
		if(i==2 && g[j]<g[j+1]) {
			k1=g[j];
			g[j]=g[j+1];
			g[j+1]=k1;
		}
	}
	for(i=0;i<k;i++) {
		match=0;
		for(j=sx=0;sx<2;sx++,j+=2) {
			match*=3;
			k1=ng[sx]*3;
			if((g1[0]=fg[i].par[sx].mat)) k1++;
			if((g1[1]=fg[i].par[sx].pat)) k1++;
			k2=0;
			switch(k1) {
				case 0: /* (0,0) (0,0) */
					k2=1;
					break;
				case 3: /* (A,0) (0,0) */
				case 6: /* (A,B) (0,0) */
					k2=2;
					break;
				case 1: /* (0,0) (X,0) */
				case 2: /* (0,0) (X,Y) */
					k2=3;
					break;
				case 4: /* (A,0) (X,0) */
					if(g1[0]==g[j]) k2=1;
					break;
				case 5: /* (A,0) (X,Y) */
					if(g1[0]==g[j] || g1[1]==g[j]) k2=3;
					break;
				case 7: /* (A,B) (X,0) */
					if(g1[0]==g[j] || g1[0]==g[j+1]) k2=2;
					break;
				case 8: /* (A,B) (X,Y) */
					if(g1[0]==g[j] && g1[1]==g[j+1]) k2=1;
					break;
			}
			if(!k2) break;
			match+=(k2-1);
		}
		match=op[match];
		/* Match found */
		if(sx==2 && match) {
			found++;
			if(match==2) {
				if(found>1) {
					fg[i].par[0].mat=-1;
					deleted++;
				} else {
					for(j=sx=0;sx<2;sx++) {
						fg[i].par[sx].mat=g[j++];
						fg[i].par[sx].pat=g[j++];
					}
				}
			} else break;
		}
	}
	if(!found) {
		if(k==*sz) {
			(*sz)<<=1;
			*fgt=lk_realloc(*fgt,sizeof(fam_gtype)*(*sz));
			fg=*fgt;
		}
		for(sx=i=0;sx<2;sx++) {
			fg[k].par[sx].mat=g[i++];
			fg[k].par[sx].pat=g[i++];
		}
		k++;
	} else if(deleted) {
		for(i=j=0;i<k;i++) {
			if(fg[i].par[0].mat>=0) {
				if(j<i) memcpy(fg+j,fg+i,sizeof(fam_gtype));
				j++;
			}
		}
		k=j;
	}
	return k;
}

/* Add parental genotypes to family set.
* Check for duplicates, but don't do anything
* special with wildcards */
static int simple_add_gt_famlist(fam_gtype **fgt,int n,int *size,int *gg)
{
	int i,t,sx;
	fam_gtype *ff;
	
	ff=*fgt;
	/* We're not bothered with order for parents */
	for(i=0;i<=2;i+=2) {
		if(gg[i]<gg[i+1]) {
			t=gg[i];
			gg[i]=gg[i+1];
			gg[i+1]=t;
		}
	}
	for(i=0;i<n;i++) {
		for(sx=t=0;sx<2;sx++,t+=2) {
			if(!(gg[t]==ff[i].par[sx].mat && gg[t+1]==ff[i].par[sx].pat)) break;
		}
		if(sx==2) break;
	}
	if(i==n) {
		if(n==*size) {
			(*size)<<=1;
			*fgt=lk_realloc(*fgt,sizeof(fam_gtype)*(*size));
			ff=*fgt;
		}
		for(sx=t=0;sx<2;sx++) {
			ff[n].par[sx].mat=gg[t++];
			ff[n].par[sx].pat=gg[t++];
		}
		n++;
	}
	return n;
}

/* Add kid's genotypes to family set. */
static int add_kid_gt_plist(gtype *gt,int n,fam_gtype *fgt,int n1,fam_gtype **fgt1,int *size)
{
	int i,sx,g[2],gg[4],a,k,k1,k2;
	
	g[0]=gt->mat;
	g[1]=gt->pat;
	if(n) {
		for(i=0;i<n;i++) {
			gg[0]=fgt[i].par[0].mat;
			gg[1]=fgt[i].par[0].pat;
			gg[2]=fgt[i].par[1].mat;
			gg[3]=fgt[i].par[1].pat;
			for(a=1,k=k1=k2=0;k1<2;k1++,k2+=2,a<<=1) {
				if(g[k1]) {
					for(sx=0;sx<2;sx++) if(g[k1]==gg[k2+sx] || !gg[k2+sx]) {
						if(!gg[k2+sx]) gg[k2+sx]=g[k1];
						k|=a;
						break;
					}
				} else k|=a;
			}
				if(k==3) n1=add_gt_famlist(fgt1,n1,size,gg);
		}
	} else {
		gg[0]=g[0];
		gg[2]=g[1];
		gg[1]=gg[3]=0;
		n1=add_gt_famlist(fgt1,n1,size,gg);
	}
		return n1;
}

void del_gtypes(struct Id_Record *id,gtype **gtypes,int *ngens)
{
	int idx;
	
	idx=id->idx;
	if(gtypes[idx]) {
		free(gtypes[idx]);
		gtypes[idx]=0;
	}
	ngens[idx]=0;
	id->flag++;
}

static void update_parents(struct Id_Record *par[],gtype **gtypes,int *ngens,int n_famgt,fam_gtype *fam_gt)
{
	int j,k,k1,k2,k4,k5,idx,sx;
	int gg[4],g[2],g1[2],sz,sz1=0;
	gtype *gt,*gt1;
	
	gt1=0;
	for(sx=0;sx<2;sx++) {
		gt=0;
		idx=par[sx]->idx;
		if(par[sx]->flag<1) continue; /* Parent is fixed.  Nothing more to add */
		k=ngens[idx];
		if(k==2) {
			if(gtypes[idx]->mat && gtypes[idx]->pat) {
				if(gtypes[idx][0].mat && gtypes[idx][0].mat==gtypes[idx][1].pat &&
					gtypes[idx][0].pat && gtypes[idx][0].pat==gtypes[idx][1].mat)
					continue; /* Parent is known apart from phase, which can't be resolved using info. from kids */
			}
		}
		k=0;
		if((sz=ngens[idx])) {
			gt=lk_malloc(sizeof(gtype)*sz);
			for(j=0;j<ngens[idx];j++) {
				k2=0;
				if((g[0]=gtypes[idx][j].mat)) k2++;
				if((g[1]=gtypes[idx][j].pat)) k2++;
				k2*=3;
				for(k1=0;k1<n_famgt;k1++) {
					gg[0]=fam_gt[k1].par[sx].mat;
					gg[1]=fam_gt[k1].par[sx].pat;
					k4=k2;
					if(gg[0]) k4++;
					if(gg[1]) k4++;
					g1[0]=g[0];
					g1[1]=g[1];
					switch(k4) {
						case 3: /* (A,0) (0,0) */
						case 6: /* (A,B) (0,0) */
							k4=1;
							break;
						case 4: /* (A,0) (X,0) */
							k4=1;
							if(g[0]) {
								if(g[0]!=gg[0]) g1[1]=gg[0];
							} else {
								if(g[1]!=gg[0]) g1[0]=gg[0];
							}
								break;
						case 5: /* (A,0) (X,Y) */
							for(k4=0;k4<2;k4++) {
								if(g[k4]) {
									for(k5=0;k5<2;k5++) if(g[k4]==gg[k5]) {
										g1[k4^1]=gg[k5^1];
										break;
									}
										if(k5<2) break;
								}
							}
							k4=(k4<2)?1:0;
							break;
						case 7: /* (A,B) (X,0) */
							if(g[0]==gg[0] || g[1]==gg[0]) k4=1;
							else k4=0;
							break;
						case 8: /* (A,B) (X,Y) */
							if((g[0]==gg[0] && g[1]==gg[1]) || (g[0]==gg[1] && g[1]==gg[0])) k4=1;
							else k4=0;
							break;
						default:
							ABT_FUNC("Internal error\n");
							break;
					}
					if(k4) k=add_ind_gtype(&gt,k,&sz,g1);
				}
			}
			if(k==1 && !(gt->mat || gt->pat)) {
				k=0;
				free(gt);
				gt=0;
			}
		} else if(n_famgt) {
			if(gt1) {
				k=sz1;
				if(sz1) {
					gt=lk_malloc(sizeof(gtype)*k);
					memcpy(gt,gt1,sizeof(gtype)*k);
				} else gt=0;
			} else {
				sz=n_famgt*2;
				gt=lk_malloc(sizeof(gtype)*sz);
				for(k1=0;k1<n_famgt;k1++) {
					gg[0]=fam_gt[k1].par[sx].mat;
					gg[1]=fam_gt[k1].par[sx].pat;
					k=add_ind_gtype(&gt,k,&sz,gg);
					if(gg[0]!=gg[1]) {
						k2=gg[0];
						gg[0]=gg[1];
						gg[1]=k2;
						k=add_ind_gtype(&gt,k,&sz,gg);
					}
				}
				if(k==1 && !(gt->mat || gt->pat)) {
					k=0;
					free(gt);
					gt=0;
				}
				gt1=gt;
				sz1=k;
			}
		}
#ifdef DEBUG
		if(k==1 && !(gt->mat || gt->pat)) ABT_FUNC("Internal error\n");
#endif		
		j=0;
		if(k==ngens[idx]) {
			for(k1=0;k1<k;k1++) {
				g[0]=gtypes[idx][k1].mat;
				g[1]=gtypes[idx][k1].pat;
				for(k2=0;k2<k;k2++) if(gt[k2].mat==g[0] && gt[k2].pat==g[1]) break;
				if(k2==k) break;
			}
			if(k1==k) j=1;
		}
		if(!j) {
			par[sx]->flag++;
			if(gtypes[idx]) free(gtypes[idx]);
			gtypes[idx]=gt;
			ngens[idx]=k;
			if(par[sx]->family) {
				if(!par[sx]->family->flag) par[sx]->family->flag=1;
			}
		} else if(gt) free(gt);
	}
}

static void update_kids(struct Id_Record **kids,gtype **gtypes,int *ngens,int n_famgt,fam_gtype *fam_gt)
{
	int j,k,k1,k2,k3,k4,k5,idx;
	int i1,i2,i3,i4;
	int gg[4],g[2],g1[2],g2[2],sz,sz1=0;
	gtype *gt,*gt1;
	struct Id_Record *kid;
	
	j=-1;
	gt1=0;
	while((kid=kids[++j])) {
		idx=kid->idx;
		if(kid->flag<0) continue; /* Kid is fixed */
		k=0;
		if((sz=ngens[idx])) {
			gt=lk_malloc(sizeof(gtype)*sz);
			for(k5=0;k5<ngens[idx];k5++) {
				k2=0;
				if((g[0]=gtypes[idx][k5].mat)) k2++;
				if((g[1]=gtypes[idx][k5].pat)) k2+=2;
				for(k1=0;k1<n_famgt;k1++) {
					gg[0]=fam_gt[k1].par[0].mat;
					gg[1]=fam_gt[k1].par[0].pat;
					gg[2]=fam_gt[k1].par[1].mat;
					gg[3]=fam_gt[k1].par[1].pat;
					i3=(gg[0]==gg[1])?1:2;
					i4=(gg[2]==gg[3])?3:4;
					for(i1=0;i1<i3;i1++) {
						g2[0]=gg[i1]; 
						for(i2=2;i2<i4;i2++) {
							g2[1]=gg[i2];
							k3=0;
							if(g2[0]) k3++;
							if(g2[1]) k3+=2;
							k3=k3*4+k2;
							g1[0]=g[0];
							g1[1]=g[1];
							k4=0;
							if(k3>3)	switch(k3) {
								case 4: /* (0,0) (A,0) */
									k4=1;
									g1[0]=g2[0];
									break;
								case 5: /* (A,0) (A,0) */
								case 7: /* (A,B) (A,0) */
									if(g[0]==g2[0]) k4=1;
									break;
								case 6: /* (0,A) (A,0) */
									g1[0]=g2[0];
									k4=1;
									break;
								case 8: /* (0,0) (0,A) */
									k4=1;
									g1[1]=g2[1];
									break;
								case 9: /* (A,0) (0,A) */
									g1[1]=g2[1];
									k4=1;
									break;
								case 10: /* (0,A) (0,A) */
								case 11: /* (A,B) (0,A) */
									if(g[1]==g2[1]) k4=1;
									break;
								case 12: /* (0,0) (A,B) */
									k4=1;
									g1[0]=g2[0];
									g1[1]=g2[1];
									break;
								case 13: /* (A,0) (A,B) */
									if(g[0]==g2[0]) {
										k4=1;
										g1[1]=g2[1];
									}
									break;
								case 14: /* (0,A) (A,B) */
									if(g[1]==g2[1]) {
										k4=1;
										g1[0]=g2[0];
									}
									break;
								case 15: /* (A,B) (A,B) */
									if(g[0]==g2[0] && g[1]==g2[1]) k4=1;
									break;
								default:
									ABT_FUNC("Internal error\n");
									break;
							} else k4=1;
							if(k4) k=add_ind_gtype(&gt,k,&sz,g1);
						}
					}
				}
			}
			if(k==1 && !(gt->mat || gt->pat)) {
				k=0;
				free(gt);
				gt=0;
			}
		} else {
			if(gt1) {
				k=sz1;
				if(sz1) {
					gt=lk_malloc(sizeof(gtype)*k);
					memcpy(gt,gt1,sizeof(gtype)*k);
				} else gt=0;
			} else {
				sz=8;
				gt=lk_malloc(sizeof(gtype)*sz);
				for(k1=0;k1<n_famgt;k1++) {
					gg[0]=fam_gt[k1].par[0].mat;
					gg[1]=fam_gt[k1].par[0].pat;
					gg[2]=fam_gt[k1].par[1].mat;
					gg[3]=fam_gt[k1].par[1].pat;
					k4=gg[0]==gg[1]?1:2;
					k5=gg[2]==gg[3]?3:4;
					for(k2=0;k2<k4;k2++) {
						g1[0]=gg[k2];
						for(k3=2;k3<k5;k3++) {
							g1[1]=gg[k3];
							k=add_ind_gtype(&gt,k,&sz,g1);
						}
					}
				}
				gt1=gt;
				if(k==1 && !(gt->mat || gt->pat)) {
					k=0;
					free(gt);
					gt=0;
				}
				sz1=k;
			}
		}
		k1=0;
		if(k<=ngens[idx]) {
			for(k1=0;k1<ngens[idx];k1++) {
				g[0]=gtypes[idx][k1].mat;
				g[1]=gtypes[idx][k1].pat;
				for(k2=0;k2<k;k2++) if((gt[k2].mat==g[0] && gt[k2].pat==g[1]) ||
											  (gt[k2].mat==g[1] && gt[k2].pat==g[0])) break;
				if(k2==k) break;
			}
			if(k1==ngens[idx]) k1=1;
			if(k<ngens[idx]) k1=2;
		}
		if(k1!=1) {
			if(gtypes[idx]) free(gtypes[idx]);
			gtypes[idx]=gt;
			ngens[idx]=k;
			kid->flag++;
			if(!k1) {
				k2=kid->nkids;
				for(k3=0;k3<k2;k3++) {
					if(!kid->kids[k3]->family->flag) kid->kids[k3]->family->flag=1;
				}
			}
		} else if(gt) free(gt);
	}
}

static gen_elim_err *check_family_pass1(nuc_fam *fam,fam_gtype **fam_gt,int *fam_gt_size,fam_gtype **fam_gt1,int *fam_gt1_size,int *n_famgt,int *ngens,gtype **gtypes)
{
	int j,k,k1,k2,k3,k4,gg[4],idx,midx,sidx;
	struct Id_Record **kids,*kid;
	gtype *gt,*gt1;
	gen_elim_err *gerr=0;
	
	/* Make initial family set */
	*n_famgt=0;
	kids=fam->kids;
	midx=fam->mother->idx;
	sidx=fam->father->idx;
	k1=ngens[midx];
	k2=ngens[sidx];
	if(!k1) {
		if(k2) {
			gt=gtypes[sidx];
			gg[0]=gg[1]=0;
			for(k=0;k<k2;k++) {
				gg[2]=gt[k].mat;
				gg[3]=gt[k].pat;
				*n_famgt=simple_add_gt_famlist(fam_gt,*n_famgt,fam_gt_size,gg);
			}
		}
	} else if(!k2) {		
		gt=gtypes[midx];
		gg[2]=gg[3]=0;
		for(k=0;k<k1;k++) {
			gg[0]=gt[k].mat;
			gg[1]=gt[k].pat;
			*n_famgt=simple_add_gt_famlist(fam_gt,*n_famgt,fam_gt_size,gg);
		}
	} else {
		gt=gtypes[midx];
		gt1=gtypes[sidx];
		for(k3=0;k3<k1;k3++) {
			gg[0]=gt[k3].mat;
			gg[1]=gt[k3].pat;
			for(k4=0;k4<k2;k4++) {
				gg[2]=gt1[k4].mat;
				gg[3]=gt1[k4].pat;
				*n_famgt=simple_add_gt_famlist(fam_gt,*n_famgt,fam_gt_size,gg);
			}
		}
	}
	/* Add in kids' information to family set; check for consistency */
	j=-1;
	while((kid=kids[++j])) {
		idx=kid->idx;
		if(!ngens[idx]) continue;
		k=*n_famgt;
		if(k>*fam_gt1_size) {
			free(*fam_gt1);
			*fam_gt1=lk_malloc(sizeof(fam_gtype)*k);
			*fam_gt1_size=k;
		}
		if(k) memcpy(*fam_gt1,*fam_gt,sizeof(fam_gtype)*k);
		*n_famgt=0;
		for(k1=0;k1<ngens[idx];k1++) *n_famgt=add_kid_gt_plist(gtypes[idx]+k1,k,*fam_gt1,*n_famgt,fam_gt,fam_gt_size);
		if(!*n_famgt) gerr=add_gen_elim_err(0,fam,kid,GEN_ELIM_PASS1);
	}
	return gerr;
}

gen_elim_err *process_mark(struct Marker *mark,struct loki *loki,int comp)
{
	int i,j,k,*prune,fam1,fam2,*haplo,*ngens,n_famgt,f1,nc;
	int nf,*perm,k1,k2,k3,k4,idx,sx,*id_status,*id_stat1,err_check;
	int gg[4],g[2],fam_gt_size=6,fam_gt1_size=6,trial,prev_ct,*err_list=0;
	int linktype;
	struct Id_Record *par[2],**kids,*kid,*id;
	gtype **gtypes;
	nuc_fam *fam;
	fam_gtype *fam_gt,*fam_gt1;
	gen_elim_err *gerr=0,*gerr1;
	
	linktype=loki->markers->linkage[mark->locus.link_group].type&LINK_TYPES_MASK;
	k=0;
	switch(linktype) {
		case LINK_AUTO:
			k=1;
			break;
		case LINK_X:
			gerr=process_mark_x(mark,loki,LINK_X,-1);
			break;
			/* Note - no point doing set recoding on Y or mitochondrial markers */
		case LINK_Y:
			gerr=process_mark_y(mark,loki,LINK_Y,-1);
			break;
		case LINK_MIT:
			gerr=process_mark_y(mark,loki,LINK_MIT,-1);
			break;
		case LINK_W:
			gerr=process_mark_y(mark,loki,LINK_W,-1);
			break;
		case LINK_Z:
			gerr=process_mark_x(mark,loki,LINK_Z,-1);
			break;
		default:
			ABT_FUNC("Unknown link type\n");
	}
	if(!k) return gerr;
	prune=mark->locus.pruned_flag;
	haplo=mark->haplo;
	gtypes=mark->gtypes;
	ngens=mark->ngens;
	id=loki->pedigree->id_array;
	fam=loki->family->families;
	err_check=loki->params.error_correct_type;
	/* Future expansion */
	if(err_check && err_check!=1) ABT_FUNC("Unknown error correction method\n");
	nc=loki->pedigree->n_comp;
	if(err_check) {
		err_list=lk_calloc((size_t)nc,sizeof(int));
	}
	/* Do all components */
	if(comp<0) {
		i=0;
		j=loki->pedigree->ped_size;
		k=j+2*loki->family->cp_start[nc];
		fam1=0;
		fam2=loki->family->cp_start[nc];
	} else {
		/* Do a particular component */
		i=loki->pedigree->comp_start[comp];
		j=i+loki->pedigree->comp_size[comp];
		fam1=loki->family->cp_start[comp];
		fam2=loki->family->cp_start[comp+1];
		k=j-i+2*(fam2-fam1);
	}
	for(;i<j;i++,id++) {
		ngens[i]=0;
		gtypes[i]=0;
		id->flag=0;
	}
	id_status=lk_malloc(sizeof(int)*k);
	id_stat1=id_status;
	fam_gt=lk_malloc(sizeof(fam_gtype)*fam_gt_size);
	fam_gt1=lk_malloc(sizeof(fam_gtype)*fam_gt1_size);
	perm=lk_malloc(sizeof(int)*(fam2-fam1));
	/* Check for singletons and pruned families */
	for(i=fam1;i<fam2;i++) {
		fam[i].flag=0;
		kids=fam[i].kids;
		/* Singleton component ? */
		if(!fam[i].father) {
			fam[i].status=0;
			get_family_gtypes(fam+i,haplo,gtypes,ngens);
		} else {
			/* Pruned ? */
			par[X_PAT]=fam[i].father;
			par[X_MAT]=fam[i].mother;
			k=0;
			if(!(prune[par[X_MAT]->idx] || prune[par[X_PAT]->idx])) {
				k1=0;
				while((kid=kids[k1++])) {
					if(!prune[kid->idx]) {
						k++;
						break;
					}
				}
			}
			if(k) {
				fam[i].flag=1;
			} else fam[i].flag=-1;
		}
	}
	message(DEBUG_MSG,"Pass 1\n");
	for(nf=0,i=0;i<fam2;i++) if(fam[i].flag>0) perm[nf++]=i;
	for(i=0;i<nf;i++) {
		f1=perm[i];
		if(err_list && err_list[fam[f1].comp]) continue;
		par[X_PAT]=fam[f1].father;
		par[X_MAT]=fam[f1].mother;
		/* Make a note of who has already been processed */
		prev_ct=0;
		for(sx=0;sx<2;sx++) if((id_stat1[sx]=ngens[par[sx]->idx])) prev_ct++;
		j=-1;
		kids=fam[f1].kids;
		while((kid=kids[++j])) if((id_stat1[2+j]=ngens[kid->idx])) prev_ct++;
		get_family_gtypes(fam+f1,haplo,gtypes,ngens);
		for(gerr1=0,trial=0;trial<2;trial++) {
			/* Assign initial genotype sets */
			gerr1=check_family_pass1(fam+f1,&fam_gt,&fam_gt_size,&fam_gt1,&fam_gt1_size,&n_famgt,ngens,gtypes);
			if(gerr1) {
				/* We hit a problem.  Try deleting passed on information to see if this resolves anything */
				if(err_check) break;
				k=0;
				for(sx=0;sx<2;sx++) if(id_stat1[sx]) {
					del_gtypes(par[sx],gtypes,ngens);
					(void)assign_gtypes(par[sx]->idx,haplo,gtypes,ngens);	
				}
					j=-1;
				kids=fam[f1].kids;
				while((kid=kids[++j])) if(id_stat1[2+j]) {
					del_gtypes(kid,gtypes,ngens);
					(void)assign_gtypes(kid->idx,haplo,gtypes,ngens);
				}
					/* Been here before...? */
					if(!(trial || !prev_ct)) {
						/*	  fam[f1].flag=-1; */
						free(gerr1);
						gerr1=0;
						continue;
					} else break;
			}
			/* Parents and kids are consistent, so now recreate sets for non-fixed individuals */
			/* Parents first */
			update_parents(par,gtypes,ngens,n_famgt,fam_gt);
			/* Now kids */
			update_kids(kids,gtypes,ngens,n_famgt,fam_gt);
			break; 
		}
		if(gerr1) {
			if(err_check) {
				free(gerr1);
				gerr1=0;
				message(DEBUG_MSG,"Error found in component %d - skipping\n",fam[f1].comp+1);
				err_list[fam[f1].comp]=1;
			} else {
				gerr1->next=gerr;
				gerr=gerr1;
				get_family_gtypes(fam+f1,haplo,gtypes,ngens);
			}
		}
		if(n_famgt) {
			fam[f1].gtypes=lk_malloc(sizeof(fam_gtype)*n_famgt);
			memcpy(fam[f1].gtypes,fam_gt,sizeof(fam_gtype)*n_famgt);
		} else fam[f1].gtypes=0;
		fam[f1].n_gt=n_famgt;
		fam[f1].status=id_stat1;
		for(sx=0;sx<2;sx++) {
			idx=par[sx]->idx;
			if(ngens[idx]==1 && gtypes[idx]->mat && gtypes[idx]->pat) par[sx]->flag=-1;
			*id_stat1++=par[sx]->flag;
		}
		j=-1;
		while((kid=kids[++j])) {
			idx=kid->idx;
			if(ngens[idx]==1 && gtypes[idx]->mat && gtypes[idx]->pat) kid->flag=-1;
			*id_stat1++=kid->flag;
		}
	}
	message(DEBUG_MSG,"Pass 2\n");
	do {
		if(err_list) {
			for(nf=0,i=0;i<fam2;i++) if(fam[i].flag>0 && !err_list[fam[i].comp]) perm[nf++]=i;
		} else {
			for(nf=0,i=0;i<fam2;i++) if(fam[i].flag>0) perm[nf++]=i;
		}
		for(i=0;i<nf;i++) {
			f1=perm[i];
			if(err_list && err_list[fam[f1].comp]) continue;
			par[X_PAT]=fam[f1].father;
			par[X_MAT]=fam[f1].mother;
			n_famgt=fam[f1].n_gt;
			if(n_famgt) memcpy(fam_gt,fam[f1].gtypes,sizeof(fam_gtype)*n_famgt);
			kids=fam[f1].kids;
			/* Add in kids' information (if changed) to family set; check for consistency */
			j=-1;
			while((kid=kids[++j])) {
				if(fam[f1].status[2+j]!=kid->flag) {
					idx=kid->idx;
					if(!ngens[idx]) continue;
					k=n_famgt;
					if(k>fam_gt1_size) {
						free(fam_gt1);
						fam_gt1=lk_malloc(sizeof(fam_gtype)*k);
						fam_gt1_size=k;
					}
					if(k) memcpy(fam_gt1,fam_gt,sizeof(fam_gtype)*k);
					n_famgt=0;
					for(k1=0;k1<ngens[idx];k1++) n_famgt=add_kid_gt_plist(gtypes[idx]+k1,k,fam_gt1,n_famgt,&fam_gt,&fam_gt_size);
					/* Inconsistency */
					if(!n_famgt) {
						if(err_check) break;
						gerr=add_gen_elim_err(gerr,fam+f1,kid,GEN_ELIM_PASS2_KID);
						/* Go back to previous state */
						n_famgt=k;
						if(k) memcpy(fam_gt,fam_gt1,sizeof(fam_gtype)*k);
						fam[f1].flag=-1;
					}
				}
			}
			if(kids[j]) {
				assert(err_check);
				message(DEBUG_MSG,"Error found in component %d - skipping\n",fam[f1].comp+1);
				err_list[fam[f1].comp]=1;
				continue;
			}
			/* Add in parents' information (if changed) to family set; check for consistency */
			for(k3=sx=0;sx<2;sx++,k3+=2) {
				if(fam[f1].status[sx]!=par[sx]->flag) {
					idx=par[sx]->idx;
					k=n_famgt;
					if(!k) k=1;
					if(k>fam_gt1_size) {
						free(fam_gt1);
						fam_gt1=lk_malloc(sizeof(fam_gtype)*k);
						fam_gt1_size=k;
					}
					if(n_famgt) memcpy(fam_gt1,fam_gt,sizeof(fam_gtype)*n_famgt);
					else {
						fam_gt1[0].par[X_MAT].mat=fam_gt1[0].par[X_MAT].pat=0;
						fam_gt1[0].par[X_PAT].mat=fam_gt1[0].par[X_PAT].pat=0;
					}
					n_famgt=0;
					for(j=0;j<ngens[idx];j++) {
						k2=0;
						if((g[0]=gtypes[idx][j].mat)) k2++;
						if((g[1]=gtypes[idx][j].pat)) k2++;
						if(g[0]!=g[1]) {
							for(k1=0;k1<j;k1++) {
								if(g[0]==gtypes[idx][k1].pat && g[1]==gtypes[idx][k1].mat) break;
							}
							if(k1<j) continue; /* Already done mirror pattern */
						}
						if(g[0]<g[1]) {
							k4=g[0];
							g[0]=g[1];
							g[1]=k4;
						}
						k2*=3;
						for(k1=0;k1<k;k1++) {
							gg[0]=fam_gt1[k1].par[X_MAT].mat;
							gg[1]=fam_gt1[k1].par[X_MAT].pat;
							gg[2]=fam_gt1[k1].par[X_PAT].mat;
							gg[3]=fam_gt1[k1].par[X_PAT].pat;
							k4=k2;
							if(gg[k3]) k4++;
							if(gg[k3+1]) k4++;
							switch(k4) {
								case 3: /* (A,0) (0,0) */
								case 6: /* (A,B) (0,0) */
									gg[k3]=g[0];
									gg[k3+1]=g[1];
									k4=2;
									break;
								case 4: /* (A,0) (X,0) */
									if(g[0]==gg[k3]) k4=1;
									else {
										k4=2;
										gg[k3+1]=g[0];
									}
										break;
								case 5: /* (A,0) (X,Y) */
									if(g[0]==gg[k3] || g[0]==gg[k3+1]) k4=1;
									else k4=0;
									break;
								case 7: /* (A,B) (X,0) */
									if(g[0]==gg[k3] || g[1]==gg[k3]) {
										k4=2;
										gg[k3]=g[0];
										gg[k3+1]=g[1];
									} else k4=0;
									break;
								case 8: /* (A,B) (X,Y) */
									if(g[0]==gg[k3] && g[1]==gg[k3+1]) k4=1;
									else k4=0;
									break;
								default:
									ABT_FUNC("Internal error\n");
									break;
							}
							if(k4) {
								n_famgt=add_gt_famlist(&fam_gt,n_famgt,&fam_gt_size,gg);
							}
						}
					}
					if(!n_famgt) {
						if(err_check) break;
						gerr=add_gen_elim_err(gerr,fam+f1,par[sx],GEN_ELIM_PASS2_PAR);
						/* Go back to previous state */
						n_famgt=k;
						if(k) memcpy(fam_gt,fam_gt1,sizeof(fam_gtype)*k);
						fam[f1].flag=-1;
					}
				}
			}
			if(sx<2) {
				assert(err_check);
				message(DEBUG_MSG,"Error found in component %s - skipping\n",fam[f1].comp+1);
				err_list[fam[f1].comp]=1;
				continue;
			}
			if(fam[f1].flag>=0) {
				update_parents(par,gtypes,ngens,n_famgt,fam_gt);
				update_kids(kids,gtypes,ngens,n_famgt,fam_gt);
			}
			if(n_famgt) {
				if(fam[f1].gtypes) free(fam[f1].gtypes);
				fam[f1].gtypes=lk_malloc(sizeof(fam_gtype)*n_famgt);
				memcpy(fam[f1].gtypes,fam_gt,sizeof(fam_gtype)*n_famgt);
			} else fam[f1].gtypes=0;
			fam[f1].n_gt=n_famgt;
			id_stat1=fam[f1].status;
			for(sx=0;sx<2;sx++) {
				idx=par[sx]->idx;
				if(ngens[idx]==1 && gtypes[idx]->mat && gtypes[idx]->pat) par[sx]->flag=-1;
				*id_stat1++=par[sx]->flag;
			}
			j=-1;
			while((kid=kids[++j])) {
				idx=kid->idx;
				if(ngens[idx]==1 && gtypes[idx]->mat && gtypes[idx]->pat) kid->flag=-1;
				*id_stat1++=kid->flag;
			}
			if(fam[f1].flag>=0) fam[f1].flag=0;
		}
	} while(nf);
	if(err_check) {
		/*    for(i=0;i<nc;i++) if(err_list[i]) process_mark_errors(mark,loki,i); */
		free(err_list);
	}
	free(perm);
	free(id_status);
	free(fam_gt);
	free(fam_gt1);
	return invert_err_list(gerr);
}

