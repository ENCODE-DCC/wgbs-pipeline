/****************************************************************************
*                                                                          *
*     Loki - Programs for genetic analysis of complex traits using MCMC    *
*                                                                          *
*                     Simon Heath - CNG, Evry                              *
*                                                                          *
*                            May 2003                                      *
*                                                                          *
* process_mark_x.c:                                                        *
*                                                                          *
* Check inheritance of x or z linked markers                               *
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
#include "lk_malloc.h"
#include "loki_utils.h"
#include "gen_elim.h"

static int het_sex; /* heterogametic sex (male in mammals, female in birds) */

/* Add individual male's genotype (female if Z linked) to list, reallocating array if required. */
static int add_ind_gtype_xm(gtype **gt,int k,int *j,const int g)
{
	int i,g1;
	gtype *gg;
	
	gg=*gt;
	if(het_sex==1) {
		if(!g) {
			k=1;
			gg->mat=0;
		} else if(!k || gg->mat) {
			for(i=0;i<k;i++) {
				g1=gg[i].mat;
				assert(g1);
				if(g==g1) break;
			}
			if(i==k) {
				if(k==*j) {
					(*j)<<=1;
					*gt=lk_realloc(*gt,sizeof(gtype)*(*j));
				}
				gg=*gt+k;
				gg->mat=g;
				k++;
			}
		}
	} else {
		if(!g) {
			k=1;
			gg->pat=0;
		} else if(!k || gg->pat) {
			for(i=0;i<k;i++) {
				g1=gg[i].pat;
				assert(g1);
				if(g==g1) break;
			}
			if(i==k) {
				if(k==*j) {
					(*j)<<=1;
					*gt=lk_realloc(*gt,sizeof(gtype)*(*j));
				}
				gg=*gt+k;
				gg->pat=g;
				k++;
			}
		}
	}
	return k;
}

/* Add parental genotypes to family set.
* Check for duplicates, but don't do anything
* special with wildcards */
static int simple_add_gt_famlist_x(fam_gtype **fgt,int n,int *size,int *gg)
{
	int i,t;
	fam_gtype *ff;
	
	ff=*fgt;
	if(het_sex==1) {
		/* We're not bothered with order for mother */
		if(gg[0]<gg[1]) {
			t=gg[0];
			gg[0]=gg[1];
			gg[1]=t;
		}
		for(i=0;i<n;i++) {
			if(gg[0]==ff[i].par[0].mat && gg[1]==ff[i].par[0].pat && gg[2]==ff[i].par[1].mat) break;
		}
	} else {
		/* or father for Z linkage */
		if(gg[2]<gg[3]) {
			t=gg[2];
			gg[2]=gg[3];
			gg[3]=t;
		}
		for(i=0;i<n;i++) {
			if(gg[1]==ff[i].par[0].pat && gg[2]==ff[i].par[1].mat && gg[3]==ff[i].par[1].pat) break;
		}
	}
	if(i==n) {
		if(n==*size) {
			(*size)<<=1;
			*fgt=lk_realloc(*fgt,sizeof(fam_gtype)*(*size));
		}
	}
	if(het_sex==1) {
		ff[n].par[0].mat=gg[0];
		ff[n].par[0].pat=gg[1];
		ff[n].par[1].mat=gg[2];
	} else {
		ff[n].par[0].pat=gg[1];
		ff[n].par[1].mat=gg[2];
		ff[n].par[1].pat=gg[3];
	}
	n++;
	return n;
}

/* Add parental genotypes to family set.
* Check for duplicates, and deal correctly with wildcards */
static int add_gt_famlist_x(fam_gtype **fgt,int k,int *sz,int *g)
{
	int par,par1,i,j,k1,k2,s,ng[2],match=0,g1[2],found=0;
	int op[]={1,1,2,1,1,0,2,0,2},deleted=0;
	fam_gtype *fg;
	
	fg=*fgt;
	s=het_sex-1; /* Homogametic parent */
	par=s*2; /* Genotype index for homogametic parent */
	par1=3-het_sex; /* Genotype index for heterogametic parent */
	for(i=0;i<2;i++) if(!g[par+i]) break;
	ng[s]=i;
	if(i==2 && g[par]<g[par+1]) {
		k1=g[par];
		g[par]=g[par+1];
		g[par+1]=k1;
	}
	ng[s^1]=g[par1]?1:0;
	for(i=0;i<k;i++) {
		/* First do homogametic parent */
		match=0;
		k1=ng[s]*3;
		if((g1[0]=fg[i].par[s].mat)) k1++;
		if((g1[1]=fg[i].par[s].pat)) k1++;
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
				if(g1[0]==g[par]) k2=1;
				break;
			case 5: /* (A,0) (X,Y) */
				if(g1[0]==g[par] || g1[1]==g[par]) k2=3;
				break;
			case 7: /* (A,B) (X,0) */
				if(g1[0]==g[par] || g1[0]==g[par+1]) k2=2;
				break;
			case 8: /* (A,B) (X,Y) */
				if(g1[0]==g[par] && g1[1]==g[par+1]) k2=1;
				break;
		}
		/* Then heterogametic parent */
		if(k2) {
			match=3*(k2-1);
			k1=ng[s^1]*2;
			g1[0]=(het_sex==1)?fg[i].par[1].mat:fg[i].par[0].pat;
			if(g1[0]) k1++;
			k2=0;
			switch(k1) {
				case 0: /* (0) (0) */
					k2=1;
					break;
				case 1: /* (0) (X) */
					k2=3;
					break;
				case 2: /* (A) (0) */
					k2=2;
					break;
				case 3: /* (A) (X) */
					if(g1[0]==g[par1]) k2=1;
					break;
			}
			if(k2) match+=k2-1;
		}
		if(k2) {
			match=op[match];
			if(match) /* Match found */ {
				found++;
				if(match==2) {
					if(found>1) {
						fg[i].par[s].mat=-1;
						deleted++;
					} else {
						fg[i].par[s].mat=g[par];
						fg[i].par[s].pat=g[par+1];
						if(het_sex==1) fg[i].par[1].mat=g[2];
						else fg[i].par[0].pat=g[1];
					}
				} else break;
			}
		}
	}
	if(!found) {
		if(k==*sz) {
			(*sz)<<=1;
			*fgt=lk_realloc(*fgt,sizeof(fam_gtype)*(*sz));
			fg=*fgt;
		}
		fg[k].par[s].mat=g[par];
		fg[k].par[s].pat=g[par+1];
		if(het_sex==1) {
			fg[k].par[1].mat=g[2];
			fg[k].par[1].pat=0;
		} else {
			fg[k].par[0].pat=g[2];
			fg[k].par[0].mat=0;
		}
		k++;
	} else if(deleted) {
		for(i=j=0;i<k;i++) {
			if(fg[i].par[s].mat>=0) {
				if(j<i) memcpy(fg+j,fg+i,sizeof(fam_gtype));
				j++;
			}
		}
		k=j;
	}
	return k;
}

/* Add kid's genotypes to family set. */
static int add_kid_gt_plist_x(gtype *gt,int n,fam_gtype *fgt,int n1,fam_gtype **fgt1,int *size,int sex)
{
	int i,j,s,par,sx,g[2],gg[4],k;
	
	s=het_sex-1; /* Homogametic parent */
	par=s*2; /* Genotype index for homogametic parent */
	if(sex==het_sex) {
		if(het_sex==1) {
			g[0]=gt->mat;
			g[1]=0;
		} else {
			g[0]=0;
			g[1]=gt->pat;
		}
	} else {
		g[0]=gt->mat;
		g[1]=gt->pat;
	}
	if(n) {
		for(i=0;i<n;i++) {
			gg[0]=fgt[i].par[0].mat;
			gg[1]=fgt[i].par[0].pat;
			gg[2]=fgt[i].par[1].mat;
			gg[3]=fgt[i].par[1].pat;
			k=0;
			if(g[s]) {
				for(sx=0;sx<2;sx++) {
					if(g[s]==gg[par+sx] || !gg[par+sx]) {
						if(!gg[par+sx]) gg[par+sx]=g[s];
						k|=1;
						break;
					}
				}
			} else k|=1;
			if(g[s^1]) {
				j=3-het_sex;
				if(g[s^1]==gg[j] || !gg[j]) {
					if(!gg[j]) gg[j]=g[s^1];
					k|=2;
				}
			} else k|=2;
			if(k==3) n1=add_gt_famlist_x(fgt1,n1,size,gg);
		}
	} else {
		if(het_sex==1) {
			gg[0]=g[0];
			gg[2]=g[1];
			gg[1]=gg[3]=0;
		} else {
			gg[1]=g[0];
			gg[2]=g[1];
			gg[0]=gg[3]=0;
		}
		n1=add_gt_famlist_x(fgt1,n1,size,gg);
	}
	return n1;
}

/* Assign initial genotype lists based on X-linked microsat marker data */
static int assign_gtypes_x(int ix,int *hap,gtype **gts,int *ngens,int sex,gen_elim_err **gerr,int flag,struct loki *loki)
{
	int k,a1,a2,fix=0,z=0;
	gtype *gt=0;
	
	if((k=hap[ix])) {
		a1=(int)(sqrt(.25+2.0*(double)k)-.49999);
		a2=k-(a1*(a1+1)/2);
		if(sex==het_sex) {
			if(a2) {
				if(a1!=a2) {
					if(!flag) {
						*gerr=add_gen_elim_err(*gerr,0,loki->pedigree->id_array+ix,GEN_ELIM_X_HET_MALE);
					}
					hap[ix]=0;
					z=1;
				}
			}
			if(!z) {
				gt=lk_malloc(sizeof(gtype));
				if(het_sex==1) {
					gt[0].mat=a1;
					gt[0].pat=0;
				} else {
					gt[0].pat=a1;
					gt[0].mat=0;
				}
				ngens[ix]=1;
				fix=1;
			}
		} else {
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
		}
		gts[ix]=gt;
	} else ngens[ix]=0;
	return fix;
}

/* Assign initial genotype lists based on microsat X-linked marker data */
static gen_elim_err *get_family_gtypes_x(nuc_fam *fam,int *hap,gtype **gts,int *ngens,int flag,struct loki *loki,gen_elim_err *gerr)
{
	int i,sx;
	struct Id_Record *par[2],**kids,*kid;
	
	par[X_PAT]=fam->father;
	par[X_MAT]=fam->mother;
	for(sx=0;sx<2;sx++) if(par[sx] && !par[sx]->flag) {
		par[sx]->flag=assign_gtypes_x(par[sx]->idx,hap,gts,ngens,2-sx,&gerr,flag,loki)?-1:1;
	}
		i=-1;
	kids=fam->kids;
	while((kid=kids[++i])) {
		if(!kid->flag) {
			kid->flag=assign_gtypes_x(kid->idx,hap,gts,ngens,kid->sex,&gerr,flag,loki)?-1:1;
		}
	}
	return gerr;
}

static void update_parents_x(struct Id_Record *par[],gtype **gtypes,int *ngens,int n_famgt,fam_gtype *fam_gt)
{
	int j,k,k1,k2,k4,k5,idx,sx,s;
	int gg[3],g[2],g1[2],sz;
	gtype *gt;
	
	s=het_sex-1;
	for(sx=0;sx<2;sx++) {
		gt=0;
		idx=par[sx]->idx;
		if(par[sx]->flag<1) continue; /* Parent is fixed.  Nothing more to add */
		k=ngens[idx];
		if(k==2 && sx==s) {
			if(gtypes[idx]->mat && gtypes[idx]->pat) {
				if(gtypes[idx][0].mat && gtypes[idx][0].mat==gtypes[idx][1].pat &&
					gtypes[idx][0].pat && gtypes[idx][0].pat==gtypes[idx][1].mat)
					continue; /* Parent is known apart from phase, which can't be resolved using info. from kids */
			}
		}
		k=0;
		if((sz=ngens[idx])) {
			gt=lk_malloc(sizeof(gtype)*sz);
			if(sx==s) {
				for(j=0;j<ngens[idx];j++) {
					k2=0;
					if((g[0]=gtypes[idx][j].mat)) k2++;
					if((g[1]=gtypes[idx][j].pat)) k2++;
					k2*=3;
					for(k1=0;k1<n_famgt;k1++) {
						gg[0]=fam_gt[k1].par[0].mat;
						gg[1]=fam_gt[k1].par[0].pat;
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
			} else {
				for(j=0;j<ngens[idx];j++) {
					g[0]=het_sex==1?gtypes[idx][j].mat:gtypes[idx][j].pat;
					assert(g[0]);
					for(k1=0;k1<n_famgt;k1++) {
						k4=0;
						gg[0]=(het_sex==1)?fam_gt[k1].par[1].mat:fam_gt[k1].par[0].pat;
						if(gg[0]) {
							if(g[0]==gg[0]) k4=1;
						} else k4=1;
						if(k4) k=add_ind_gtype_xm(&gt,k,&sz,g[0]);
					}
				}
				j=het_sex==1?gt->mat:gt->pat;
				if(k==1 && !j ) {
					k=0;
					free(gt);
					gt=0;
				}
			}
		} else if(n_famgt) {
			sz=n_famgt*2;
			gt=lk_malloc(sizeof(gtype)*sz);
			if(sx==s) {
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
			} else {
				for(k1=0;k1<n_famgt;k1++) {
					gg[0]=het_sex==1?fam_gt[k1].par[1].mat:fam_gt[k1].par[0].pat;
					k=add_ind_gtype_xm(&gt,k,&sz,gg[0]);
				}
				k1=het_sex==1?gt->mat:gt->pat;
				if(k==1 && !k1) {
					k=0;
					free(gt);
					gt=0;
				}
			}
		}
		j=0;
		if(k==ngens[idx]) {
			if(sx==s) {
				for(k1=0;k1<k;k1++) {
					g[0]=gtypes[idx][k1].mat;
					g[1]=gtypes[idx][k1].pat;
					for(k2=0;k2<k;k2++) if(gt[k2].mat==g[0] && gt[k2].pat==g[1]) break;
					if(k2==k) break;
				}
			} else {
				if(het_sex==1) {
					for(k1=0;k1<k;k1++) {
						g[0]=gtypes[idx][k1].mat;
						for(k2=0;k2<k;k2++) if(gt[k2].mat==g[0]) break;
						if(k2==k) break;
					}
				} else {
					for(k1=0;k1<k;k1++) {
						g[0]=gtypes[idx][k1].pat;
						for(k2=0;k2<k;k2++) if(gt[k2].pat==g[0]) break;
						if(k2==k) break;
					}
				}
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

static void update_kids_x(struct Id_Record **kids,gtype **gtypes,int *ngens,int n_famgt,fam_gtype *fam_gt)
{
	int j,k,k1,k2,k3,k4,k5,idx,sex,i1,i2,i3,i4;
	int gg[4],g[2],g1[2],g2[2],sz,sz1[2]={0,0};
	gtype *gt,*gt1[2];
	struct Id_Record *kid;
	
	j=-1;
	while((kid=kids[++j])) {
		idx=kid->idx;
		if(kid->flag<0) continue; /* Kid is fixed */
		k=0;
		sex=kid->sex;
		if((sz=ngens[idx])) {
			gt=lk_malloc(sizeof(gtype)*sz);
			if(sex!=het_sex) {
				for(k5=0;k5<ngens[idx];k5++) {
					k2=0;
					if((g[0]=gtypes[idx][k5].mat)) k2++;
					if((g[1]=gtypes[idx][k5].pat)) k2+=2;
					for(k1=0;k1<n_famgt;k1++) {
						gg[0]=fam_gt[k1].par[0].mat;
						gg[1]=fam_gt[k1].par[0].pat;
						gg[2]=fam_gt[k1].par[1].mat;
						gg[3]=fam_gt[k1].par[1].pat;
						i3=(het_sex==2 || gg[0]==gg[1])?1:2;
						i4=(het_sex==1 || gg[2]==gg[3])?3:4;
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
				for(k5=0;k5<ngens[idx];k5++) {
					g[0]=het_sex==1?gtypes[idx][k5].mat:gtypes[idx][k5].pat;
					assert(g[0]);
					for(k1=0;k1<n_famgt;k1++) {
						gg[0]=fam_gt[k1].par[het_sex-1].mat;
						gg[1]=fam_gt[k1].par[het_sex-1].pat;
						i3=(gg[0]==gg[1])?1:2;
						for(i1=0;i1<i3;i1++) {
							k4=0;
							if(gg[i1]) {
								if(gg[i1]==g[0]) k4=1;
							} else k4=1;
							if(k4) k=add_ind_gtype_xm(&gt,k,&sz,g[0]);
						}
					}
				}
				k1=het_sex==1?gt->mat:gt->pat;
				if(k==1 && !k1) {
					k=0;
					free(gt);
					gt=0;
				}
			}
		} else {
			k1=2-sex;
			k=sz1[k1];
			if(k) {
				gt=lk_malloc(sizeof(gtype)*k);
				memcpy(gt,gt1[k1],sizeof(gtype)*k);
			} else {
				sz=8;
				gt=lk_malloc(sizeof(gtype)*sz);
				if(sex!=het_sex) {
					for(k1=0;k1<n_famgt;k1++) {
						gg[0]=fam_gt[k1].par[0].mat;
						gg[1]=fam_gt[k1].par[0].pat;
						gg[2]=fam_gt[k1].par[1].mat;
						gg[3]=fam_gt[k1].par[1].pat;
						k4=(het_sex==2 || gg[0]==gg[1])?1:2;
						k5=(het_sex==1 || gg[2]==gg[3])?3:4;
						for(k2=0;k2<k4;k2++) {
							g1[0]=gg[k2];
							for(k3=2;k3<k5;k3++) {
								g1[1]=gg[k3];
								k=add_ind_gtype(&gt,k,&sz,g1);
							}
						}
					}
					if(k==1 && !(gt->mat || gt->pat)) {
						k=0;
						free(gt);
						gt=0;
					}
				} else {
					for(k1=0;k1<n_famgt;k1++) {
						gg[0]=fam_gt[k1].par[het_sex-1].mat;
						gg[1]=fam_gt[k1].par[het_sex-1].pat;
						k4=gg[0]==gg[1]?1:2;
						for(k2=0;k2<k4;k2++) {
							k=add_ind_gtype_xm(&gt,k,&sz,gg[k2]);
						}
					}
					k1=het_sex==1?gt->mat:gt->pat;
					if(k==1 && !k1) {
						k=0;
						free(gt);
						gt=0;
					}
				}
				gt1[2-sex]=gt;
				sz1[2-sex]=k;
			}
		}
		k1=0;
		if(sex!=het_sex) {
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
		} else {
			if(k==ngens[idx]) {
				if(het_sex==1) {
					for(k1=0;k1<k;k1++) {
						g[0]=gtypes[idx][k1].mat;
						for(k2=0;k2<k;k2++) if(gt[k2].mat==g[0]) break;
						if(k2==k) break;
					}
				} else {
					for(k1=0;k1<k;k1++) {
						g[0]=gtypes[idx][k1].pat;
						for(k2=0;k2<k;k2++) if(gt[k2].pat==g[0]) break;
						if(k2==k) break;
					}
				}
				if(k1==ngens[idx]) k1=1;
			}
		} 
		if(k1!=1) {
			if(gtypes[idx]) free(gtypes[idx]);
			gtypes[idx]=gt;
			ngens[idx]=k;
			kid->flag++;
			if(!k1) {
				k2=kid->nkids;
				for(k3=0;k3<k2;k3++) {
					if(sex==het_sex && kid->kids[k3]->sex==het_sex) continue; /* Only flag female kids of males */
					if(!kid->kids[k3]->family->flag) kid->kids[k3]->family->flag=1;
				}
			}
		} else if(gt) free(gt);
	}
}

gen_elim_err *process_mark_x(struct Marker *mark,struct loki *loki,int ltype,int comp)
{
	int i,j,k,*prune,fam1,fam2,*haplo,*ngens,n_famgt,f1;
	int nf,*perm,k1,k2,k3,k4,idx,sx,*id_status,*id_stat1;
	int gg[3],g[2],fam_gt_size=6,fam_gt1_size=6,trial,prev_ct;
	struct Id_Record *par[2],**kids,*kid,*id;
	gtype **gtypes,*gt,*gt1;
	nuc_fam *fam;
	fam_gtype *fam_gt,*fam_gt1;
	gen_elim_err *gerr=0;
	
	prune=mark->locus.pruned_flag;
	haplo=mark->haplo;
	gtypes=mark->gtypes;
	ngens=mark->ngens;
	id=loki->pedigree->id_array;
	fam=loki->family->families;
	het_sex=ltype==LINK_X?1:2;
	if(comp<0) {
		/* Do all components */
		i=0;
		j=loki->pedigree->ped_size;
		k=j+2*loki->family->cp_start[loki->pedigree->n_comp];
		fam1=0;
		fam2=loki->family->cp_start[loki->pedigree->n_comp];
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
	if(!(id_status=malloc(sizeof(int)*k))) ABT_FUNC(MMsg);
	id_stat1=id_status;
	if(!(fam_gt=malloc(sizeof(fam_gtype)*fam_gt_size))) ABT_FUNC(MMsg);
	if(!(fam_gt1=malloc(sizeof(fam_gtype)*fam_gt1_size))) ABT_FUNC(MMsg);
	if(!(perm=malloc(sizeof(int)*(fam2-fam1)))) ABT_FUNC(MMsg);
	for(i=fam1;i<fam2;i++) {
		/* Check for singletons and pruned families */
		fam[i].flag=0;
		kids=fam[i].kids;
		if(!fam[i].father) {
			/* Singleton component ? */
			fam[i].status=0;
			gerr=get_family_gtypes_x(fam+i,haplo,gtypes,ngens,0,loki,gerr);
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
			/* Family is pruned for this locus */
		}
	}
	for(nf=0,i=0;i<fam2;i++) if(fam[i].flag>0) perm[nf++]=i;
	message(DEBUG_MSG,"Pass 1\n");
	for(i=0;i<nf;i++) {
		f1=perm[i];
		par[X_PAT]=fam[f1].father;
		par[X_MAT]=fam[f1].mother;
		/* Make a note of who has already been processed */
		prev_ct=0;
		for(sx=0;sx<2;sx++) if((id_stat1[sx]=ngens[par[sx]->idx])) prev_ct++;
		j=-1;
		kids=fam[f1].kids;
		while((kid=kids[++j])) if((id_stat1[2+j]=ngens[kid->idx])) prev_ct++;
		for(trial=0;trial<2;trial++) { 
			/* We go around the second time if we hit a problem! */
			/* Assign initial genotype sets */
			gerr=get_family_gtypes_x(fam+f1,haplo,gtypes,ngens,trial,loki,gerr);
			/* Make initial family set */
			n_famgt=0;
			k1=ngens[par[0]->idx];
			k2=ngens[par[1]->idx];
			if(!k1) {
				/* Mother not typed */
				gt=gtypes[par[1]->idx];
				gg[0]=gg[1]=0;
				for(k=0;k<k2;k++) {
					gg[2]=gt[k].mat;
					n_famgt=simple_add_gt_famlist_x(&fam_gt,n_famgt,&fam_gt_size,gg);
				}
			} else if(!k2) {
				/* Father not typed */
				gt=gtypes[par[0]->idx];
				gg[2]=0;
				for(k=0;k<k1;k++) {
					gg[0]=gt[k].mat;
					gg[1]=gt[k].pat;
					n_famgt=simple_add_gt_famlist_x(&fam_gt,n_famgt,&fam_gt_size,gg);
				}
			} else {
				/* Both parents typed */
				gt=gtypes[par[0]->idx];
				gt1=gtypes[par[1]->idx];
				for(k3=0;k3<k1;k3++) {
					gg[0]=gt[k3].mat;
					gg[1]=gt[k3].pat;
					for(k4=0;k4<k2;k4++) {
						gg[2]=gt1[k4].mat;
						n_famgt=simple_add_gt_famlist_x(&fam_gt,n_famgt,&fam_gt_size,gg);
					}
				}
			}
			/* Add in kids' information to family set; check for consistency */
			j=-1;
			while((kid=kids[++j])) {
				idx=kid->idx;
				if(!ngens[idx]) continue;
				k=n_famgt;
				if(k>fam_gt1_size) {
					free(fam_gt1);
					if(!(fam_gt1=malloc(sizeof(fam_gtype)*k))) ABT_FUNC(MMsg);
					fam_gt1_size=k;
				}
				if(k) memcpy(fam_gt1,fam_gt,sizeof(fam_gtype)*k);
				n_famgt=0;
				for(k1=0;k1<ngens[idx];k1++) n_famgt=add_kid_gt_plist_x(gtypes[idx]+k1,k,fam_gt1,n_famgt,&fam_gt,&fam_gt_size,kid->sex);
				if(!n_famgt) {
					/* This error is due to this family alone */
					if(trial || !prev_ct) 
						gerr=add_gen_elim_err(gerr,fam+f1,kid,GEN_ELIM_PASS1);
					break;
				}
			}
			if(kid) { 
				/* We hit a problem.  Try deleting passed on information to see if this resolves anything */
				k=0;
				for(sx=0;sx<2;sx++) if(id_stat1[sx]) del_gtypes(par[sx],gtypes,ngens);
				j=-1;
				kids=fam[f1].kids;
				while((kid=kids[++j])) if(id_stat1[2+j]) del_gtypes(kid,gtypes,ngens);
				fam[f1].flag=-1;
				if(!(trial || !prev_ct)) { 
					/* Been here before... */
					continue;
				} else {
					(void)get_family_gtypes_x(fam+f1,haplo,gtypes,ngens,1,loki,0);
					break;
				}
			}
			/* Parents and kids are consistent, so now recreate sets for non-fixed individuals */
			/* Parents first */
			update_parents_x(par,gtypes,ngens,n_famgt,fam_gt);
			/* Now kids */
			update_kids_x(kids,gtypes,ngens,n_famgt,fam_gt);
			break; 
		}
		if(n_famgt) {
			if(!(fam[f1].gtypes=malloc(sizeof(fam_gtype)*n_famgt))) ABT_FUNC(MMsg);
			memcpy(fam[f1].gtypes,fam_gt,sizeof(fam_gtype)*n_famgt);
		} else fam[f1].gtypes=0;
		fam[f1].n_gt=n_famgt;
		fam[f1].status=id_stat1;
		for(sx=0;sx<2;sx++) {
			idx=par[sx]->idx;
			if(ngens[idx]==1 && gtypes[idx]->mat && (sx==1 || gtypes[idx]->pat)) par[sx]->flag=-1;
			*id_stat1++=par[sx]->flag;
		}
		j=-1;
		while((kid=kids[++j])) {
			idx=kid->idx;
			if(ngens[idx]==1 && gtypes[idx]->mat && (kid->sex==het_sex || gtypes[idx]->pat)) kid->flag=-1;
			*id_stat1++=kid->flag;
		}
	}
	message(DEBUG_MSG,"Pass 2\n");
	do {
		for(nf=i=0;i<fam2;i++) if(fam[i].flag>0) perm[nf++]=i;
		for(i=0;i<nf;i++) {
			f1=perm[i];
			par[X_PAT]=fam[f1].father;
			par[X_MAT]=fam[f1].mother;
			n_famgt=fam[f1].n_gt;
			/*      printf("A: \n");
			print_family_gtypes(stdout,fam+f1,mark,ngens,gtypes); */
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
						if(!(fam_gt1=malloc(sizeof(fam_gtype)*k))) ABT_FUNC(MMsg);
						fam_gt1_size=k;
					}
					if(k) memcpy(fam_gt1,fam_gt,sizeof(fam_gtype)*k);
					n_famgt=0;
					for(k1=0;k1<ngens[idx];k1++) n_famgt=add_kid_gt_plist_x(gtypes[idx]+k1,k,fam_gt1,n_famgt,&fam_gt,&fam_gt_size,kid->sex);
					if(!n_famgt) {
						/* Inconsistency */
						gerr=add_gen_elim_err(gerr,fam+f1,kid,GEN_ELIM_PASS2_KID);
						/* Go back to previous state */
						n_famgt=k;
						if(k) memcpy(fam_gt,fam_gt1,sizeof(fam_gtype)*k);
						fam[f1].flag=-1;
					}
				}
			}
			/* Add in parents' information (if changed) to family set; check for consistency */
			for(sx=0;sx<2;sx++) {
				if(fam[f1].status[sx]!=par[sx]->flag) {
					idx=par[sx]->idx;
					k=n_famgt;
					if(!k) k=1;
					if(k>fam_gt1_size) {
						free(fam_gt1);
						if(!(fam_gt1=malloc(sizeof(fam_gtype)*k))) ABT_FUNC(MMsg);
						fam_gt1_size=k;
					}
					if(n_famgt) memcpy(fam_gt1,fam_gt,sizeof(fam_gtype)*n_famgt);
					else {
						fam_gt1[0].par[X_MAT].mat=fam_gt1[0].par[X_MAT].pat=0;
						fam_gt1[0].par[X_PAT].mat=fam_gt1[0].par[X_PAT].pat=0;
					}
					n_famgt=0;
					if(!sx) {
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
								k4=k2;
								if(gg[0]) k4++;
								if(gg[1]) k4++;
								switch(k4) {
									case 3: 
										/* (A,0) (0,0) */
									case 6: 
										/* (A,B) (0,0) */
										gg[0]=g[0];
										gg[1]=g[1];
										k4=2;
										break;
									case 4: 
										/* (A,0) (X,0) */
										if(g[0]==gg[0]) k4=1;
										else {
											k4=2;
											gg[1]=g[0];
										}
										break;
									case 5: 
										/* (A,0) (X,Y) */
										if(g[0]==gg[0] || g[0]==gg[1]) k4=1;
										else k4=0;
										break;
									case 7: 
										/* (A,B) (X,0) */
										if(g[0]==gg[0] || g[1]==gg[0]) {
											k4=2;
											gg[0]=g[0];
											gg[1]=g[1];
										} else k4=0;
										break;
								case 8: 
									/* (A,B) (X,Y) */
									if(g[0]==gg[0] && g[1]==gg[1]) k4=1;
									else k4=0;
									break;
								default:
									ABT_FUNC("Internal error\n");
									break;
								}
								if(k4) n_famgt=add_gt_famlist_x(&fam_gt,n_famgt,&fam_gt_size,gg);
							}
						}
					} else {
						for(j=0;j<ngens[idx];j++) {
							g[0]=gtypes[idx][j].mat;
							for(k1=0;k1<k;k1++) {
								k4=0;
								gg[0]=fam_gt1[k1].par[X_MAT].mat;
								gg[1]=fam_gt1[k1].par[X_MAT].pat;
								if((gg[2]=fam_gt1[k1].par[X_PAT].mat)) {
									if(g[0]==gg[2]) k4=1;
								} else {
									gg[2]=g[0];
									k4=2;
								}
								if(k4) n_famgt=add_gt_famlist_x(&fam_gt,n_famgt,&fam_gt_size,gg);
							}
						}
					}
					if(!n_famgt) {
						gerr=add_gen_elim_err(gerr,fam+f1,par[sx],GEN_ELIM_PASS2_PAR);
						/* Go back to previous state */
						n_famgt=k;
						if(k) memcpy(fam_gt,fam_gt1,sizeof(fam_gtype)*k);
						fam[f1].flag=-1;
					}
				}
			}
			if(fam[f1].flag>=0) {
				update_parents_x(par,gtypes,ngens,n_famgt,fam_gt);
				update_kids_x(kids,gtypes,ngens,n_famgt,fam_gt);
			}
			if(n_famgt) {
				if(fam[f1].gtypes) free(fam[f1].gtypes);
				if(!(fam[f1].gtypes=malloc(sizeof(fam_gtype)*n_famgt))) ABT_FUNC(MMsg);
				memcpy(fam[f1].gtypes,fam_gt,sizeof(fam_gtype)*n_famgt);
			} else fam[f1].gtypes=0;
			fam[f1].n_gt=n_famgt;
			id_stat1=fam[f1].status;
			for(sx=0;sx<2;sx++) {
				idx=par[sx]->idx;
				if(ngens[idx]==1 && gtypes[idx]->mat && (sx==1 || gtypes[idx]->pat)) par[sx]->flag=-1;
				*id_stat1++=par[sx]->flag;
			}
			j=-1;
			while((kid=kids[++j])) {
				idx=kid->idx;
				if(ngens[idx]==1 && gtypes[idx]->mat && (kid->sex==het_sex || gtypes[idx]->pat)) kid->flag=-1;
				*id_stat1++=kid->flag;
			}
			if(fam[f1].flag>=0) fam[f1].flag=0;
			/*      printf("B: \n");
			print_family_gtypes(stdout,fam+f1,mark,ngens,gtypes); */
		}
	} while(nf);
	free(perm);
	free(id_status);
	free(fam_gt);
	free(fam_gt1);
	return invert_err_list(gerr);
}
