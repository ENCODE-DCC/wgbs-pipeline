/****************************************************************************
 *                                                                          *
 *     Loki - Programs for genetic analysis of complex traits using MCMC    *
 *                                                                          *
 *                     Simon Heath - CNG, Evry                              *
 *                                                                          *
 *                         February 2003                                    *
 *                                                                          *
 * gen_elim.c:                                                              *
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
#include "loki.h"
#include "loki_peel.h"
#include "loki_utils.h"
#include "ped_utils.h"
#include "gen_elim.h"

static gt_node *free_gt_node_list;

void check_fixed(struct Id_Record *id,gtype **gtypes,int linktype)
{
	int j;
	
	j=id->idx;
	switch(linktype) {
	 case LINK_AUTO:
		if(gtypes[j]->mat && gtypes[j]->pat) id->flag=-1;
		break;
	 case LINK_Y:
		if(gtypes[j]->pat) id->flag=-1;
		break;
	 case LINK_X:
		if(id->sex==1) {
			if(gtypes[j]->mat) id->flag=-1;
		} else if(gtypes[j]->mat && gtypes[j]->pat) id->flag=-1;
		break;
	 default:
		ABT_FUNC("Unknown linktype\n");
	}
}

/* Assign initial genotype lists based on marker data */
static int assign_gtypes(int ix,int *hap,gtype **gts,int *ngens,int linktype,int sex,char *name)
{
	int k,a1,a2,err=0;
	gtype *gt;
	
	if((k=hap[ix])) {
		a1=(int)(.50001+sqrt(.25+2.0*(double)(k-1)));
		a2=k-(a1*(a1-1)/2);
		if(linktype==LINK_X) {
			if(sex==1 && a1!=a2) {
				(void)fprintf(stderr,"Error - X-linked marker %s: Individual ",name);
				print_orig_id(stderr,ix+1);
				(void)fputs(" is male and heterozygous",stderr);
				hap[ix]=0;
				err++;
			}
		} else if(linktype==LINK_Y) {
			if(sex==2 || a1!=a2) {
				(void)fprintf(stderr,"Error - Y-linked marker %s: Individual ",name);
				print_orig_id(stderr,ix+1);
				if(sex==2) (void)fputs(" is female and has genetic data\n",stderr);
				else (void)fputs( " is heterozygous\n",stderr);
				a1=a2=0;
				hap[ix]=0;
				err++;
			}
		}
		if(a1) { /* If typed, set gtype list */
			if(a1!=a2) {
				if(!(gt=malloc(sizeof(gtype)*2))) ABT_FUNC(MMsg);
				gt[0].mat=gt[1].pat=a1;
				gt[0].pat=gt[1].mat=a2;
				ngens[ix]=2;
			} else {
				if(!(gt=malloc(sizeof(gtype)))) ABT_FUNC(MMsg);
				if(linktype==LINK_Y) {
					gt[0].pat=a1;
					gt[0].mat=0;
				} else if (linktype==LINK_X && sex==1) {
					gt[0].mat=a1;
					gt[0].pat=0;
				} else gt[0].mat=gt[0].pat=a1;
				ngens[ix]=1;
			}
			gts[ix]=gt;
		} else ngens[ix]=0;
	} else ngens[ix]=0;
	return err;
}

static void add_gt_x(gtype *gt,fam_gtype *fgt,int *ct,int sex)
{
	int i,j,k,k1,ct1,sx;
	int ff[3][3];
	unsigned int g[2],m[2],p;
	
	g[0]=gt->mat;
	g[1]=gt->pat;
	if(!*ct) {
		if(sex==1) {
			fgt[0].par[X_MAT].mat=g[0];
			fgt[0].par[X_MAT].pat=fgt[0].par[X_PAT].mat=fgt[0].par[X_PAT].pat=0;
			*ct=1;
		} else {
			for(sx=0;sx<2;sx++) {
				fgt[0].par[sx].mat=g[sx];
				fgt[1].par[sx].mat=g[sx^1];
				fgt[0].par[sx].pat=fgt[1].par[sx].pat=0;
			}
			*ct=2;
		}
	} else {
		ct1=0;
		for(i=0;i<*ct;i++) {
			m[0]=fgt[i].par[X_MAT].mat;
			m[1]=fgt[i].par[X_MAT].pat;
			p=fgt[i].par[X_PAT].mat;
			for(j=0;j<2;j++) if(!m[j]) break;
			if(sex==2) {
				if(j&&p) {
					for(k1=0;k1<j;k1++) if((g[0]==m[k1] && g[1]==p)||(g[1]==m[k1] && g[0]==p)) break;
				} else k1=j;
			} else {
				if(j) {
					for(k1=0;k1<j;k1++) if(g[0]==m[k1]) break;
				} else k1=j;
			}
			if(k1<j) { /* Match to previous set */
				for(k=0;k<2;k++) ff[ct1][k]=m[k];
				ff[ct1][2]=p;
				ct1++;
			} else { /* No match, check if unused slots */
				if(j<2) {
					if(sex==1) {
						if(g[0]!=m[0]) {
							m[1]=g[0];
							for(k=0;k<2;k++) ff[ct1][k]=m[k];
							ff[ct1++][2]=p;
						}
					} else {
						if(!p) {
							if(g[0]!=m[0]) {
								m[1]=g[0];
								for(k=0;k<2;k++) ff[ct1][k]=m[k];
								ff[ct1++][2]=g[1];
							}
							if(g[0]!=g[1] && m[0]!=g[1]) {
								m[1]=g[1];
								for(k=0;k<2;k++) ff[ct1][k]=m[k];
								ff[ct1++][2]=g[0];
							}
						} else {
							if(p==g[1]) {
								m[1]=g[0];
								for(k=0;k<2;k++) ff[ct1][k]=m[k];
								ff[ct1++][2]=p;
							}
							if(g[0]!=g[1] && p==g[0]) {
								m[1]=g[1];
								for(k=0;k<2;k++) ff[ct1][k]=m[k];
								ff[ct1++][2]=p;
							}
						}
					}
				} else if(!p) {
#ifdef DEBUG
					if(sex==1) ABT_FUNC("Internal error\n");
#endif
					if(g[0]==m[0] || g[0]==m[1]) {
						for(k=0;k<2;k++) ff[ct1][k]=m[k];
						ff[ct1++][2]=g[1];
					}
					if(g[0]!=g[1] && (g[1]==m[0] || g[1]==m[1])) {
						for(k=0;k<2;k++) ff[ct1][k]=m[k];
						ff[ct1++][2]=g[0];
					}
				}
			}
		}
		for(i=0;i<ct1;i++) {
			if(ff[i][0]>ff[i][1]) {
				fgt[i].par[X_MAT].mat=ff[i][0];
				fgt[i].par[X_MAT].pat=ff[i][1];
			} else {
				fgt[i].par[X_MAT].mat=ff[i][1];
				fgt[i].par[X_MAT].pat=ff[i][0];
			}
			fgt[i].par[X_PAT].mat=ff[i][2];
			fgt[i].par[X_PAT].pat=0;
		}
		*ct=ct1;
	}
#ifdef DEBUG
	if(*ct>3) ABT_FUNC("Internal error - too many combinations\n");
#endif
}

static void add_gt(gtype *gt,fam_gtype *fgt,int *ct)
{
	int i,j,k,k1,k2,ct1,sx;
	int ff[3][4];
	unsigned int g[2],m[2],p[2];
	
	g[0]=gt->mat;
	g[1]=gt->pat;
	if(!*ct) {
		for(sx=0;sx<2;sx++) {
			fgt[0].par[sx].mat=g[sx];
			fgt[0].par[sx].pat=0;
		}
		(*ct)++;
	} else {
		ct1=0;
		for(i=0;i<*ct;i++) {
			m[0]=fgt[i].par[X_MAT].mat;
			m[1]=fgt[i].par[X_MAT].pat;
			p[0]=fgt[i].par[X_PAT].mat;
			p[1]=fgt[i].par[X_PAT].pat;
			for(j=0;j<2;j++) if(!m[j]) break;
			for(k=0;k<2;k++) if(!p[k]) break;
			if(j&&k) {
				for(k1=0;k1<j;k1++) {
					for(k2=0;k2<k;k2++) if((g[0]==m[k1] && g[1]==p[k2]) || (g[1]==m[k1] && g[0]==p[k2])) break;
					if(k2<k) break;
				}
			} else k1=j;
			if(k1<j) { /* Match to previous set */
				for(k=0;k<2;k++) {
					ff[ct1][k]=m[k];
					ff[ct1][k+2]=p[k];
				}
				ct1++;
			} else { /* No Match, check if unused slots */
				if(j<2) {
					if(k<2) {
						if(g[0]!=m[0]) m[1]=g[0];
						if(g[1]!=p[0]) p[1]=g[1];
						for(k=0;k<2;k++) {
							ff[ct1][k]=m[k];
							ff[ct1][k+2]=p[k];
						}
						ct1++;
						if(g[0]!=g[1] && m[0]!=p[0]) {
							m[1]=(g[1]!=m[0])?g[1]:0;
							p[1]=(g[0]!=p[0])?g[0]:0;
							for(k=0;k<2;k++) {
								ff[ct1][k]=m[k];
								ff[ct1][k+2]=p[k];
							}
							ct1++;
						}
					} else {
						if(g[1]==p[0] || g[1]==p[1]) {
							m[j]=g[0];
							for(k1=0;k1<2;k1++) {
								ff[ct1][k1]=m[k1];
								ff[ct1][k1+2]=p[k1];
							}
							ct1++;
						}
						if(g[0]!=g[1] && (g[0]==p[0] || g[0]==p[1])) {
							m[j]=g[1];
							for(k1=0;k1<2;k1++) {
								ff[ct1][k1]=m[k1];
								ff[ct1][k1+2]=p[k1];
							}
							ct1++;
						}
					}
				} else if(k<2) {
					if(g[0]==m[0] || g[0]==m[1]) {
						p[k]=g[1];
						for(k1=0;k1<2;k1++) {
							ff[ct1][k1]=m[k1];
							ff[ct1][k1+2]=p[k1];
						}
						ct1++;
					}
					if(g[0]!=g[1] && (g[1]==m[0] || g[1]==m[1])) {
						p[k]=g[0];
						for(k1=0;k1<2;k1++) {
								ff[ct1][k1]=m[k1];
							ff[ct1][k1+2]=p[k1];
						}
						ct1++;
					}
				}
			}
		}
		for(i=0;i<ct1;i++) {
			for(k=sx=0;sx<2;sx++,k+=2) {
				if(ff[i][k]>ff[i][k+1]) {
					fgt[i].par[sx].mat=ff[i][k];
					fgt[i].par[sx].pat=ff[i][k+1];
				} else {
					fgt[i].par[sx].mat=ff[i][k+1];
					fgt[i].par[sx].pat=ff[i][k];
				}
			}
		}
		*ct=i;
	}
#ifdef DEBUG
	if(*ct>3) ABT_FUNC("Internal error - too many combinations\n");
#endif
}

static void add_gt_y(gtype *gt,fam_gtype *fgt,int *ct)
{
	int g;
	
	g=gt->pat;
	if(!*ct) {
		fgt[0].par[X_PAT].pat=g;
		fgt[0].par[X_MAT].mat=fgt[0].par[X_MAT].pat=fgt[0].par[X_PAT].pat=0;
		(*ct)++;
	} else if(fgt[0].par[X_PAT].pat!=g) *ct=0;
}

gt_node *new_gt_node(void)
{
	gt_node *p;
	
	if(free_gt_node_list) {
		p=free_gt_node_list;
		free_gt_node_list=p->next;
	} else if(!(p=malloc(sizeof(gt_node)))) ABT_FUNC(MMsg);
	return p;
}

gt_node *find_gt_node(gt_node *p,int i)
{
	while(p) {
		if(p->i>=i) break;
		p=p->next;
	}
	if(p && p->i==i) return p;
	return 0;
}

void add_gt_node(gt_node **p,int i)
{
	gt_node *p1;
	
	while((p1=*p)) {
		if(p1->i>=i) break;
		p=&p1->next;
	}
	if(!p1 || p1->i!=i) {
		*p=new_gt_node();
		(*p)->i=i;
		(*p)->next=p1;
	}
}

/* Find intersect of 2 genotype lists in *gt1 *gt2.  Returns length of output list, and
 * puts output list in *gt_out.  Needs to account for wildcard genotypes (0) */
static int get_intersect(gtype *gt1,gtype *gt_out,int n1,gtype *gt2,int n2,int n_all,gt_node **work,int *iwork)
{
	int i,j,k,g1,g2,*iwork1;
	gt_node **work1,*p1,*p;

	work1=work+n_all;
	iwork1=iwork+n_all;
	for(i=0;i<n_all*2;i++) work[i]=0;
	memset(iwork,0,sizeof(int)*2*n_all);
 	for(i=0;i<n1;i++) {
		g1=gt1[i].mat;
		g2=gt1[i].pat;
		if(!g1) for(j=0;j<n_all;j++) add_gt_node(work+j,g2);
		else if(!g2) for(j=0;j<n_all;j++) add_gt_node(work+g1-1,j+1);
		else add_gt_node(work+g1-1,g2);
	}
 	for(i=0;i<n2;i++) {
		g1=gt2[i].mat;
		g2=gt2[i].pat;
		if(!g1) for(j=0;j<n_all;j++) add_gt_node(work1+j,g2);
		else if(!g2) for(j=0;j<n_all;j++) add_gt_node(work1+g1-1,j+1);
		else add_gt_node(work1+g1-1,g2);
	}
	/* Get intersect of the two sets */
	for(i=0;i<n_all;i++) {
		p=work[i];
		p1=work1[i];
		k=0;
		while(p&&p1) {
			if(p->i==p1->i) {
				k++;
				iwork1[p->i-1]++;
				p=p->next;
				p1=p1->next;
			} else {
				if(p->i>p1->i) p1=p1->next;
				else p=p->next;
			}
		}
		iwork[i]=k;
	}
	/* Store */
	for(i=j=0;i<n_all;i++) {
		if(iwork1[i]==n_all) {
			gt_out[j].mat=0;
			gt_out[j++].pat=i+1;
		}
		if(iwork[i]==n_all) {
			gt_out[j].mat=i+1;
			gt_out[j++].pat=0;
		} else {
			p=work[i];
			p1=work1[i];
			while(p&&p1) {
				if(p->i==p1->i) {
					if(iwork1[p->i-1]!=n_all) {
						gt_out[j].mat=i+1;
						gt_out[j++].pat=p->i;
					}
					p=p->next;
					p1=p1->next;
				} else {
					if(p->i>p1->i) p1=p1->next;
					else p=p->next;
				}
			}
		}
		if((p=work[i])) {
			while(p->next) p=p->next;
			p->next=free_gt_node_list;
			free_gt_node_list=work[i];
		}
		if((p=work1[i])) {
			while(p->next) p=p->next;
			p->next=free_gt_node_list;
			free_gt_node_list=work1[i];
		}
	}
	return j;
}

int gen_elim_marker(int ix,struct loki *loki) 
{
	int i,j,k,comp,*prune,err=0,fam1,fam2,linktype,*haplo,*ngens,ngt,n_famgt;
	int nf,*perm,ordered,k1,k2,k3,k4,n_all,*iwork,p_idx[2],sx,nkid[2],*id_status,*id_stat1;
	int gg[4];
	struct Marker *mark;
	struct Locus *loc;
	char *name;
	struct Id_Record *par[2],**kids,*kid,*id;
	gtype **gtypes,gtlist[12],gtlist1[12],*gt;
	gt_node **work;
	nuc_fam *fam;
	fam_gtype fam_gt[6];
	
	mark=loki->markers.marker+ix;
	name=mark->name;
	loc=&mark->locus;
	prune=loc->pruned_flag;
	haplo=mark->haplo;
	gtypes=mark->gtypes;
	ngens=mark->ngens;
	message(DEBUG_MSG,"Genotype elimination for microsat locus %s\n",name);
	id=loki->pedigree.id_array;
	for(i=0;i<loki->pedigree.ped_size;i++,id++) {
		ngens[i]=0;
		gtypes[i]=0;
		id->flag=0;
	}
	i+=2*loki->family->cp_start[loki->pedigree.n_comp];
	if(!(id_status=malloc(sizeof(int)*i))) ABT_FUNC(MMsg);
	id_stat1=id_status;
	n_all=mark->locus.n_alleles;
	i=n_all*(n_all+1)/2;
	if(!(iwork=malloc(sizeof(int)*(n_all+i)*2))) ABT_FUNC(MMsg);
	if(!(work=malloc(sizeof(void *)*(i*2+n_all)))) ABT_FUNC(MMsg);
	/* First pass - check for errors within each nuclear family */
	/* We will *not* pass information from one family to another at
	 * this stage to make error reporting easier */
	message(DEBUG_MSG,"Pass 1\n",name);
	linktype=loki->markers.linkage[loc->link_group].type&LINK_TYPES;
	fam1=0;
	fam=loki->family->families;
	for(comp=0;comp<loki->pedigree.n_comp;comp++) {
		n_all=mark->n_all1[comp];
		fam2=loki->family->cp_start[comp+1];
		for(;fam1<fam2;fam1++) {
			fam[fam1].flag=0;
			kids=fam[fam1].kids;
			if(!fam[fam1].father) {
				/* This component is just unrelated individuals */
				fam[fam1].status=0;
				i=ngt=0;
				while((kid=kids[i++])) {
					j=kid->idx;
					if(!prune[j]) err+=assign_gtypes(j,haplo,gtypes,ngens,linktype,kid->sex,name);
				}
			} else { /* Normal component */
				par[X_PAT]=fam[fam1].father;
				par[X_MAT]=fam[fam1].mother;
				/* Check if family is pruned for this locus */
				k=0;
				if(!prune[par[X_MAT]->idx] || !prune[par[X_PAT]->idx]) k=1;
				else {
					i=0;
					while((kid=kids[i++])) {
						if(!prune[kid->idx]) {
							k++;
							break;
						}
					}
				}
				if(k) {
					n_famgt=0;
					/* Get genotypes of parents.  Start construction of family genotype list */
					for(sx=0;sx<2;sx++) {
						p_idx[sx]=par[sx]->idx;
						if(!par[sx]->flag) {
							k=assign_gtypes(p_idx[sx],haplo,gtypes,ngens,linktype,1,name);
							par[sx]->flag=1;
							err+=k;
						} else k=0;
						if(!k && haplo[p_idx[sx]]) {
							gt=gtypes[p_idx[sx]];
							fam_gt[0].par[sx].mat=gt->mat;
							fam_gt[0].par[sx].pat=gt->pat;
							n_famgt=1;
						} else fam_gt[0].par[sx].mat=fam_gt[0].par[sx].pat=0;
					}
					/* If at least one parent is typed then the parents are not symmetric */
					ordered=n_famgt;
					ngt=0;
					i=-1;
					/* Get genotypes of kids; construct family genotype list, merging information
					 from parents */
					while((kid=kids[++i])) {
						j=kid->idx;
						if(!prune[j]) {
							if(!kid->flag) {
								k=assign_gtypes(j,haplo,gtypes,ngens,linktype,kid->sex,name);
								kid->flag=1;
								err+=k;
							} else k=0;
							if(!k && haplo[j]) {
								gt=gtypes[j];
								for(k=0;k<ngt;k++) if(gtlist[k].mat==gt->mat && gtlist[k].pat==gt->pat) break;
								if(k==ngt) {
									if(ngt==4) n_famgt=0;
									else {
										gtlist[ngt].mat=gt->mat;
										gtlist[ngt].pat=gt->pat;
										switch(linktype) {
										 case LINK_X:
											add_gt_x(gtlist+ngt,fam_gt,&n_famgt,kid->sex);
											break;
										 case LINK_Y:
											add_gt_y(gtlist+ngt,fam_gt,&n_famgt);
											break;
										 case LINK_AUTO:
											add_gt(gtlist+ngt,fam_gt,&n_famgt);
											break;
										 default:
											ABT_FUNC("Bad link type\n");
										}
										ngt++;
									}
									if(!n_famgt) {
										fprintf(stderr,"\nMendelian error in nuclear family for marker %s when adding child ",mark->name);
										print_orig_id(stderr,j+1);
										fputc('\n',stderr);
#ifdef DEBUG
										print_family_gtypes(stderr,fam+fam1,mark,ngens,gtypes);
#else
										print_family_gtypes(stderr,fam+fam1,mark,0,0);
#endif
										fputc('\n',stderr);
										err++;
										exit(0);
										break;
									}
								}
							}
						}
					}
					/* If parents not typed, create mirrors of family gtypes */
					if(!ordered && linktype==LINK_AUTO) {
						j=n_famgt;
						for(i=0;i<j;i++) {
							if(fam_gt[i].par[X_PAT].mat!=fam_gt[i].par[X_MAT].mat || fam_gt[i].par[X_PAT].pat !=fam_gt[i].par[X_MAT].pat) {
								fam_gt[n_famgt].par[X_PAT].mat=fam_gt[i].par[X_MAT].mat;
								fam_gt[n_famgt].par[X_PAT].pat=fam_gt[i].par[X_MAT].pat;
								fam_gt[n_famgt].par[X_MAT].mat=fam_gt[i].par[X_PAT].mat;
								fam_gt[n_famgt++].par[X_MAT].pat=fam_gt[i].par[X_PAT].pat;
							}
						}
					}
					fam[fam1].n_gt=n_famgt;
					if(n_famgt) {
						/* Create genotype lists for untyped parents */
						for(sx=0;sx<2;sx++) if(!haplo[p_idx[sx]]) {
							if(linktype==LINK_Y) {
								if(sx==X_MAT) continue;
								if(n_famgt>1) ABT_FUNC("This shouldn't happen\n");
								gtlist[0].mat=0;
								gtlist[0].pat=fam_gt[0].par[X_PAT].pat;
								j=1;
							} else if(linktype==LINK_X && sx==X_PAT) {
								for(i=j=0;i<n_famgt;i++) {
									gt=&fam_gt[i].par[sx];
									if(!gt->mat) break;
									for(k=0;k<j;k++) if(gt->mat==gtlist[k].mat) break;
									if(k==j) {
										gtlist[j].mat=gt->mat;
										gtlist[j++].pat=0;
									}
								}
								if(i<n_famgt) j=0;
							} else {
								for(j=i=0;i<n_famgt;i++) { /* First pick up combinations with wild cards */
									gt=&fam_gt[i].par[sx];
									if(!gt->pat) for(k=0;k<j;k++) if(gtlist[k].mat==gt->mat) break;
									if(k==j) {
										gtlist[j].mat=gt->mat;
										gtlist[j++].pat=0;
									}
								}
								/* Now pick up remaining combinations */
								for(i=0;i<n_famgt;i++) {
									gt=&fam_gt[i].par[sx];
									if(gt->pat) for(k=0;k<j;k++) {
										if(!gtlist[k].pat) {
											if(gtlist[k].mat==gt->mat || gtlist[k].mat==gt->pat) break;
										} else if(gtlist[k].mat==gt->mat && gtlist[k].pat==gt->pat) break;
									}
									if(k==j) {
#ifdef DEBUG
										if(j==6) ABT_FUNC("Internal error - too many genotypes\n");
#endif
										gtlist[j].mat=gt->mat;
										gtlist[j++].pat=gt->pat;
									}
								}
								k=j;
								for(i=0;i<k;i++) {
									gt=gtlist+i;
									if(gt->mat!=gt->pat) {
										gtlist[j].mat=gt->pat;
										gtlist[j++].pat=gt->mat;
									}
								}
							}
							/* Remove special case with (0,0) genotype */
							if(j==1 && !(gtlist[0].mat||gtlist[0].pat)) j=0;
							if(j) {
								gt=0;
								if(ngens[p_idx[sx]]) {
									k=get_intersect(gtlist,gtlist1,j,gtypes[p_idx[sx]],ngens[p_idx[sx]],n_all,work,iwork);
									if(k) {
										if(!(gt=malloc(sizeof(gtype)*k))) ABT_FUNC(MMsg);
										memcpy(gt,gtlist1,sizeof(gtype)*k);
									} else {
										/* Ignore inconsistencies for now - they will be picked up in stage 2, and
										 will be easier to report there */
										message(DEBUG_MSG,"Inconsistency noted.  Will be reported later\n");
									}
								}
								if(!gt) {
									k=j;
									if(!(gt=malloc(sizeof(gtype)*k))) ABT_FUNC(MMsg);
									memcpy(gt,gtlist,sizeof(gtype)*k);
								}
								par[sx]->flag++;
								ngens[p_idx[sx]]=k;
								gtypes[p_idx[sx]]=gt;
							}
						}
						/* Create genotype lists for unfixed kids */
						/* First check if there are any... */
						i=0;
						nkid[0]=nkid[1]=0;
						while((kid=kids[++i])) {
							if(kid->flag>=0) nkid[2-kid->sex]++;
						}
						if(linktype==LINK_Y) {
							if(nkid[1] && fam_gt[0].par[X_PAT].pat) {
								i=0;
								gg[0]=fam_gt[0].par[X_PAT].pat;
								while((kid=kids[++i])) {
									if(kid->flag>=0 && kid->sex==1) {
										if(!(gt=malloc(sizeof(gtype)))) ABT_FUNC(MMsg);
										gt->pat=gg[0];
										gt->mat=0;
										ngens[kid->idx]=1;
										gtypes[kid->idx]=gt; 	
									}
								}
							}
						} else if(linktype==LINK_X && nkid[1]) {
							/* First we look for genotype possibilities with wildcards */
							for(i=0;i<n_famgt;i++) {
								if(!fam_gt[i].par[X_MAT].pat) break;
							}
							if(i==n_famgt) for(i=j=0;i<n_famgt;i++) {
								gg[0]=fam_gt[i].par[X_MAT].mat;
								gg[1]=fam_gt[i].par[X_MAT].pat;
								for(k3=0;k3<2;k3++) {
									for(k=0;k<j;k++) if(gtlist[k].mat==gg[k3]) break;
									if(k==j) {
										gtlist[j].mat=gg[k3];
										gtlist[j++].pat=0;
									}
								}
								i=0;
								while((kid=kids[++i])) {
									if(kid->flag>=0 && kid->sex==1) {
										if(ngens[kid->idx]) {
											k1=0;
											gt=gtypes[kid->idx];
											for(k=0;k<ngens[kid->idx];k++) {
												gg[0]=gt[k].mat;
												for(k2=0;k2<j;k2++) if(gg[0]==gtlist[k2].mat) break;
												if(k2<j) gt[k1++].mat=gg[0];
											}
											ngens[kid->idx]=k1;
											if(!k1) {
												free(gtypes[kid->idx]);
												/* Ignore inconsistencies for now - they will be picked up in stage 2, and
												 will be easier to report there */
												message(DEBUG_MSG,"Inconsistency noted.  Will be reported later\n");
											}
										}
										if(!ngens[kid->idx]) {
											if(!(gt=malloc(sizeof(gtype)*j))) ABT_FUNC(MMsg);
											ngens[kid->idx]=j;
											gtypes[kid->idx]=gt; 
										}
									}
								}
							}
							
						} else if((fam_gt[0].par[X_MAT].pat || fam_gt[0].par[X_PAT].pat) && (nkid[0] || (nkid[1] && linktype==LINK_AUTO))) {
							/* First we look for genotype possibilities with wildcards */
							for(j=i=0;i<n_famgt;i++) {
								if(!fam_gt[i].par[X_MAT].pat || !fam_gt[i].par[X_PAT].pat) {
									gg[0]=fam_gt[i].par[X_MAT].mat;
									gg[1]=fam_gt[i].par[X_MAT].pat;
									gg[2]=fam_gt[i].par[X_PAT].mat;
									gg[3]=fam_gt[i].par[X_PAT].pat;
									k1=(gg[0]==gg[1])?1:2;
									if(linktype==LINK_AUTO) k2=(gg[2]==gg[3])?3:4;
									else k2=3;
									for(k3=0;k3<k1;k3++) for(k4=2;k4<k2;k4++) {
										if(gg[k3] && gg[k4]) continue;
										for(k=0;k<j;k++) {
											if(gtlist[k].mat==gg[k3]) {
												if(gtlist[k].pat==gg[k4] || !gtlist[k].pat) break;
											} else if(!gtlist[k].mat && gtlist[k].pat==gg[k4]) break;
										}
										if(k==j) {
#ifdef DEBUG
											if(j==4) {
												ABT_FUNC("Internal error - too many genotypes\n");
											}
#endif
											gtlist[j].mat=gg[k3];
											gtlist[j++].pat=gg[k4];
										}
									}
								}
							}
							for(i=0;i<n_famgt;i++) {
								gg[0]=fam_gt[i].par[X_MAT].mat;
								gg[1]=fam_gt[i].par[X_MAT].pat;
								gg[2]=fam_gt[i].par[X_PAT].mat;
								gg[3]=fam_gt[i].par[X_PAT].pat;
								k1=(gg[0]==gg[1])?1:2;
								if(linktype==LINK_AUTO) k2=(gg[2]==gg[3])?3:4;
								else k2=3;
								for(k3=0;k3<k1;k3++) for(k4=2;k4<k2;k4++) {
									if(!gg[k3] || !gg[k4]) continue;
									for(k=0;k<j;k++) {
										if(gtlist[k].mat==gg[k3]) {
											if(gtlist[k].pat==gg[k4] || !gtlist[k].pat) break;
										} else if(!gtlist[k].mat && gtlist[k].pat==gg[k4]) break;
									}
									if(k==j) {
#ifdef DEBUG
										if(j==12) {
											ABT_FUNC("Internal error - too many genotypes\n");
										}
#endif
										gtlist[j].mat=gg[k3];
										gtlist[j++].pat=gg[k4];
									}
								}
							}
							i=-1;
							while((kid=kids[++i])) {
								if(kid->flag>=0 && (kid->sex==2 || linktype==LINK_AUTO)) {
									gt=0;
									if(ngens[kid->idx]) {
										k=get_intersect(gtlist,gtlist1,j,gtypes[kid->idx],ngens[kid->idx],n_all,work,iwork);
										if(k) {
											if(!(gt=malloc(sizeof(gtype)*k))) ABT_FUNC(MMsg);
											memcpy(gt,gtlist1,sizeof(gtype)*k);
										} else {
											/* Ignore inconsistencies for now - they will be picked up in stage 2, and
											 will be easier to report there */
											message(DEBUG_MSG,"Inconsistency noted.  Will be reported later\n");
										}
									}
									if(!gt) {
										k=j;
										if(!(gt=malloc(sizeof(gtype)*k))) ABT_FUNC(MMsg);
										memcpy(gt,gtlist,sizeof(gtype)*k);
									}
									kid->flag++;
									ngens[kid->idx]=k;
									gtypes[kid->idx]=gt; 
								}
							}
						}
						if(!(fam[fam1].gtypes=malloc(sizeof(fam_gtype)*n_famgt))) ABT_FUNC(MMsg);
						memcpy(fam[fam1].gtypes,fam_gt,sizeof(fam_gtype)*n_famgt);
					} else fam[fam1].gtypes=0;
					fam[fam1].status=id_stat1;
					for(sx=0;sx<2;sx++) {
						if(ngens[p_idx[sx]]==1) check_fixed(par[sx],gtypes,linktype);
						*id_stat1++=par[sx]->flag;
					}
					i=-1;
					while((kid=kids[++i])) {
						if(ngens[kid->idx]==1) check_fixed(kid,gtypes,linktype);
						*id_stat1++=kid->flag;
					}
					fam[fam1].flag=1;
				} else {
					fam[fam1].status=0;
					fam[fam1].flag=-1;
					printf("Pruned: ");
					print_orig_id(stdout,fam[fam1].father->idx+1);
					fputc(' ',stdout);
					print_orig_id(stdout,fam[fam1].mother->idx+1);
					fputc('\n',stdout);
				}
			}
		}
	}
	if(!err) {
		/* Second pass - iterative genotype elimination */
		message(DEBUG_MSG,"Pass 2\n",name);
		fam1=0;
		fam=loki->family->families;
		for(comp=0;comp<loki->pedigree.n_comp;comp++) {
			fam2=loki->family->cp_start[comp+1];
			if(!(perm=malloc(sizeof(int)*(fam2-fam1)))) ABT_FUNC(MMsg);
			do {
				for(nf=0,i=fam1;i<fam2;i++) if(fam[i].flag>0) perm[nf++]=i;
				printf("comp %d, nf=%d\n",comp,nf);
				switch(linktype) {
				 case LINK_AUTO:
					for(i=0;i<nf && !err;i++) {
						j=perm[i];
						err+=recheck_fam(fam+j,gtypes,ngens,haplo);
					}
					break;
				 case LINK_X:
					for(i=0;i<nf && !err;i++) {
						j=perm[i];
						err+=recheck_fam_x(fam+j,gtypes,ngens,haplo);
					}
					break;
				 case LINK_Y:
				 default:
					ABT_FUNC("Unknown link type\n");
				}
			} while(nf && !err);
			free(perm);
			fam1=fam2;
		}
	}
	free(work);
	free(iwork);
	free(id_status);	
	return err;
}

int Gen_Elim(struct loki *loki) 
{
	int i,err=0;
	gt_node *p;
	
	message(INFO_MSG,"Genotype elimination\n");
	for(i=0;i<loki->markers.n_markers;i++) {
		recode_alleles(i,loki);
		err+=gen_elim_marker(i,loki);
	}
	while(free_gt_node_list) {
		p=free_gt_node_list->next;
		free(free_gt_node_list);
		free_gt_node_list=p;
	}
	return err;
}
