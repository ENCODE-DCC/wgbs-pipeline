/****************************************************************************
 *                                                                          *
 *     Loki - Programs for genetic analysis of complex traits using MCMC    *
 *                                                                          *
 *             Simon Heath - University of Washington                       *
 *                                                                          *
 *                        October 1997                                      *
 *                                                                          *
 * meiosis_scan.c:                                                          *
 *                                                                          *
 * Perform update of segregation indicators using meiosis sampler           *
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
#include <math.h>
#include <stdio.h>
#include <assert.h>

#include "ranlib.h"
#include "utils.h"
#include "loki.h"
#include "lk_malloc.h"
#include "loki_peel.h"
#include "loki_utils.h"
#include "seg_pen.h"
#include "gen_pen.h"
#include "meiosis_scan.h"

static int *nnfd_st,*par_st,*gpfam_st,*fam_st,mem_list_size,mem_list_ptr;
static int *nnfd_list,*par_list,*fam_list,*gpfam_list,*temp_list;
static struct nuc_family *families;
static void **mem_list;

static double mscan_prob[]={MSCAN_INDIVIDUAL,
	  MSCAN_HS_FAMILY,
	  MSCAN_GP_FAMILY,
	  MSCAN_FS_FAMILY};

int set_mscan_probs(double *p)
{
	int i,err=0;
	double z=0.0;
	char *emsg[]={"Negative probabilies","Probabilities all zero"};
	
	for(i=0;i<4;i++) {
		if(p[i]<0.0) {
			err=1;
			break;
		}
		z+=p[i];
		
	}
	if(!err) {
		if(z==0.0) err=2;
		else {
			for(i=0;i<4;i++) mscan_prob[i]=p[i]/z;
			if(check_output(INFO_MSG)) {
				fputs("Setting mscan probabilities:\n",stdout);
				printf("  p(Individual updates) = %g\n",mscan_prob[0]);
				printf("  p(Halfsib family updates) = %g\n",mscan_prob[1]);
				printf("  p(Fullsib family updates) = %g\n",mscan_prob[2]);
				printf("  p(Grand parental family updates) = %g\n",mscan_prob[3]);
			}
		}
	}
	if(err) message(WARN_MSG,"set_mscan_probs(): Illegal probabilities specified (%s)",emsg[err-1]);
	return err;
}

static void *realloc_and_remember(void *p,size_t size)
{
	int i;
	
	for(i=0;i<mem_list_ptr;i++) if(mem_list[i]==p) break;
	if(i<mem_list_ptr) {
		p=realloc(p,size);
		mem_list[i]=p;
	} else p=0;
	return p;
}

static void *malloc_and_remember(size_t size)
{
	void *p;
	
	if(mem_list_ptr==mem_list_size) {
		if(mem_list_size) {
			mem_list_size*=1.5;
			mem_list=lk_realloc(mem_list,sizeof(void *)*mem_list_size);
		} else {
			mem_list_size=16;
			mem_list=lk_malloc(sizeof(void *)*mem_list_size);
		}
	}
	p=malloc(size);
	if(p) mem_list[mem_list_ptr++]=p;
	return p;
}

static void init_families(const struct loki *loki)
{
	int j,i1,cs,comp,nnfd,npar,nfam,ngpfam;
	int k,k1,k2,nk,ids,idd,*kid_list,n_comp,ped_size;
	struct Id_Record *id_array,*kid;
	
	n_comp=loki->pedigree->n_comp;
	ped_size=loki->pedigree->ped_size;
	if(!(nnfd_st=malloc_and_remember(sizeof(int)*(n_comp+1)*4))) ABT_FUNC(MMsg);
	par_st=nnfd_st+n_comp+1;
	fam_st=par_st+n_comp+1;
	gpfam_st=fam_st+n_comp+1;
	if(!(nnfd_list=malloc_and_remember(sizeof(int)*3*ped_size))) ABT_FUNC(MMsg);
	par_list=nnfd_list+ped_size;
	temp_list=par_list+ped_size;
	id_array=loki->pedigree->id_array;
	for(j=0;j<ped_size;j++) id_array[j].flag=0;
	for(nk=nfam=ngpfam=nnfd=npar=j=comp=0;comp<n_comp;comp++) {
		nnfd_st[comp]=nnfd;
		par_st[comp]=npar;
		cs=loki->pedigree->comp_size[comp];
		for(i1=0;i1<cs;i1++,j++) {
			if(id_array[j].sire) {
				nnfd_list[nnfd++]=j;
				if(!id_array[j].flag) {
					ids=id_array[j].sire;
					idd=id_array[j].dam;
					for(k1=k=0;k<id_array[ids-1].nkids;k++) {
						kid=id_array[ids-1].kids[k];
						if(kid->dam==idd) {
							kid->flag=1;
							nk++;
							if(1 ||kid->nkids) k1=1;
						}
					}
					nfam++;
					ngpfam+=k1;
				}
			}
			if(id_array[j].nkids) par_list[npar++]=j;
		}
	}
	nnfd_st[comp]=nnfd;
	par_st[comp]=npar;
	if(!nfam) {
		for(comp=0;comp<=n_comp;comp++) fam_st[comp]=gpfam_st[comp]=0;
		return;
	}
	if(!(families=malloc_and_remember(sizeof(struct nuc_family)*nfam))) ABT_FUNC(MMsg);
	if(!(fam_list=malloc_and_remember(sizeof(int)*(nfam+ngpfam+nk)))) ABT_FUNC(MMsg);
	gpfam_list=fam_list+nfam;
	kid_list=gpfam_list+ngpfam;
	for(nk=nfam=ngpfam=j=comp=0;comp<n_comp;comp++) {
		fam_st[comp]=nfam;
		gpfam_st[comp]=ngpfam;
		cs=loki->pedigree->comp_size[comp];
		for(i1=0;i1<cs;i1++,j++) {
			if(id_array[j].sire) {
				if(id_array[j].flag<2) {
					ids=id_array[j].sire;
					idd=id_array[j].dam;
					families[nfam].kids=kid_list+nk;
					for(k2=k1=k=0;k<id_array[ids-1].nkids;k++) {
						kid=id_array[ids-1].kids[k];
						if(kid->dam==idd) {
							kid->flag=2;
							kid_list[nk++]=kid->idx;
							k2++;
							if(1 || kid->nkids) k1=1;
						}
					}
					fam_list[nfam]=nfam;
 					if(k1) gpfam_list[ngpfam++]=nfam;
					families[nfam++].nkids=k2;
				}
			}
		}
	}
	fam_st[comp]=nfam;
	gpfam_st[comp]=ngpfam;
}

static void free_meiosis_scan(void)
{
	meiosis_scan(-1,0);
}

void meiosis_scan(int link,const struct loki *loki)
{
	int i,i1,j,k,k1,k2,k3,kk,n_loci,s,s1,s2,s3,s1a=0,s2a=0,ss,locus,comp;
	int idd,ids,mf,ffg[4],scl,update_type,cs,nkids,state;
	int sex,**seg,**seg1,kid1,kid,nkids1,*kids,lffg,symflag,n_qtl;
	double x,x1,xx,xx1,xx2,rr[2],z,z1,z2,zz,zz1[4];
	static double *recom[2],*pp[4],*pen[4],*lks[4],**lk_store;
	static int *loc_fg,***seg_list,ctr,*tmp_arr=0,tmp_arr_size=0;
	static struct Locus **loci1,*qtl;
	struct Id_Record *id_array,**kids1;
	struct Marker *marker;

	if(link<0) { /* Free space and prevent reuse */
		if(mem_list) {
			for(i=0;i<mem_list_ptr;i++) if(mem_list[i]) free(mem_list[i]);
			free(mem_list);
			mem_list=0;
			mem_list_size=mem_list_ptr=0;
		}
		return;
	}
	n_loci=loki->markers->linkage[link].n_markers+loki->params.max_tloci;
	if(!n_loci) return;
	id_array=loki->pedigree->id_array;
	marker=loki->markers->marker;
	if(!pp[0]) { /* Allocate space on first time through */
		for(j=i=0;i<loki->markers->n_links;i++) {
			k=loki->markers->linkage[i].n_markers;
			if(k>j) j=k;
		}
		j+=loki->params.max_tloci;
		if(!(pp[0]=malloc_and_remember(sizeof(double)*(14*j-1)))) ABT_FUNC(MMsg);
		pen[0]=pp[0]+4*j;
		lks[0]=pen[0]+4*j;
		recom[0]=lks[0]+4*j;
		recom[1]=recom[0]+j;
		for(i=1;i<4;i++) {
			pp[i]=pp[i-1]+j;
			pen[i]=pen[i-1]+j;
			lks[i]=lks[i-1]+j;
		}
		if(!(seg_list=malloc_and_remember(sizeof(void *)*j*2))) ABT_FUNC(MMsg);
		lk_store=(void *)(seg_list+j);
		if(!(loc_fg=malloc_and_remember(sizeof(int)*j))) ABT_FUNC(MMsg);
		if(!(loci1=malloc_and_remember(sizeof(void *)*j))) ABT_FUNC(MMsg);
		init_families(loki);
		/* Normalize mscan type probabilities */
		for(z1=0.0,update_type=0;update_type<4;update_type++) z1+=mscan_prob[update_type];
		if(z1<=0.0) ABT_FUNC("Invalid probabilities set for mscan_prob[]\n");
		z1=1.0/z1;
		for(update_type=0;update_type<4;update_type++) mscan_prob[update_type]*=z1;
		tmp_arr_size=32;
		if(!(tmp_arr=malloc_and_remember(sizeof(int)*tmp_arr_size))) ABT_FUNC(MMsg);
		if(atexit(free_meiosis_scan)) message(WARN_MSG,"Unable to register exit function free_meiosis_scan()\n");
	}
	z=safe_ranf();
	for(z1=0.0,update_type=0;update_type<4;update_type++) {
		z1+=mscan_prob[update_type];
		if(z<=z1) break;
	}
	switch(update_type) {
	 case 0:
		memcpy(temp_list,nnfd_list,sizeof(int)*nnfd_st[loki->pedigree->n_comp]);
		break;
	 case 1:
		memcpy(temp_list,par_list,sizeof(int)*par_st[loki->pedigree->n_comp]);
		break;
	 case 2:
		memcpy(temp_list,gpfam_list,sizeof(int)*gpfam_st[loki->pedigree->n_comp]);
		break;
	 case 3:
		memcpy(temp_list,fam_list,sizeof(int)*fam_st[loki->pedigree->n_comp]);
		break;
	}
	for(comp=0;comp<loki->pedigree->n_comp;comp++) if(!(loki->pedigree->singleton_flag[comp])) {
		qtl=0;
		get_locuslist(loci1,link,&k,0);
		/* Pick loci in this component where no. alleles >=2 */
		n_qtl=k1=0;
		for(k2=0;k2<k;k2++) {
			if(loci1[k2]->type&ST_MARKER) {
				if(loki->markers->marker[loci1[k2]->index].n_all1[comp]>=2) loci1[k1++]=loci1[k2];
			} else {
				n_qtl++;
				qtl=loci1[k2];
			}
		}
		if(!k1) continue;
		n_loci=k1;
		/* If we have >1 linked QTL, then pick one to update conditional on trait data alone,
		 * the others will be updated conditional on sampled genotypes for observed individuals */
		if(n_qtl>1) {
			k1=(int)(safe_ranf()*(double)n_qtl);
			for(i1=k2=0;k2<k;k2++) if(loci1[k2]->type&ST_TRAITLOCUS) {
				if(i1++==k1) {
					qtl=loci1[k2];
					break;
				}
			}
			for(k2=0;k2<k;k2++) if((loci1[k2]->type&ST_TRAITLOCUS) && qtl!=loci1[k2]) {
				z=seg_pen(loci1[k2],comp,&i1,0,loki);
				if(i1) {
					(void)fprintf(stderr,"%d %d %d %d\n",k1,comp,i,(int)z);
					ABT_FUNC("Illegal segregation pattern\n");
				}
				loci1[k2]->lk_store[comp]=z;
			}
		} 
		gnu_qsort(loci1,(size_t)n_loci,sizeof(void *),cmp_loci);
		for(s=0;s<2;s++) {
			x=loci1[0]->pos[s];
			for(k=1;k<n_loci;k++) {
				x1=loci1[k]->pos[s];
				recom[s][k-1]=.5*(1.0-exp(0.02*(x-x1)));
				x=x1;
			}
		}
		for(k=0;k<n_loci;k++) {
			seg_list[k]=loci1[k]->seg;
			lk_store[k]=loci1[k]->lk_store;
		}
		if(!update_type) {
			/* Update both segs of each non-founder individual */
			cs=nnfd_st[comp+1];
			i1=nnfd_st[comp];
			zz=(double)(cs-i1);
			for(;i1<cs;i1++) {
				k=i1+(int)(safe_ranf()*(zz--));
				i=temp_list[k];
				temp_list[k]=temp_list[i1];
				temp_list[i1]=i;
				ids=id_array[i].sire;
				idd=id_array[i].dam;
				for(k=0;k<n_loci;k++) { /* Calculate penetrances */
					locus=(loci1[k]->type&ST_MARKER)?loci1[k]->index:-1-loci1[k]->index;
					if(loci1[k]->pruned_flag[i]) {
						/* Individual pruned for this locus - 
						 * seg has no effect on penetrance */
						for(k1=0;k1<4;k1++) pen[k1][k]=0.25;
						continue;
					}
					seg=seg_list[k];
					if(locus>=0) {
						mf=marker[locus].m_flag[i]&6;
						/* marker[locus].m_flag[i] is a flag where the following bits mean:
						 *  0  i is homozygous at locus
						 *  1  maternal seg is forced by the data
						 *  2  paternal seg is forced by the data
						 * 
						 * Note: this is all based on the original marker data */
						if(marker[locus].m_flag[idd-1]&1) mf|=8;
						if(marker[locus].m_flag[ids-1]&1) mf|=16;
					} else mf=0;
					if(!id_array[idd-1].sire && id_array[idd-1].nkids==1) mf|=8; 
					if(!id_array[ids-1].sire && id_array[ids-1].nkids==1) mf|=16;
					s1=s1a=seg[X_MAT][i];
					s2=s2a=seg[X_PAT][i];
					assert(s1>=0 && s2>=0);
					s=(s1<<1)|s2;
					scl=0;
					xx2=xx1=lk_store[k][comp];
					if(!mf) { /* Default case - both segs must be checked */
						/* Current state */
						lks[s][k]=xx1;
						ffg[s]=0;
						ss=s;
						for(k2=k1=0;k1<3;k1++) {
							/* Flip maternal bit once */
							if(k1&1) {
								seg[X_MAT][i]=(s1a^=1);
								ss^=2;
 								if(pass_founder_genes1(locus,i,X_MAT,loki)) {
									xx=xx1=seg_pen1(loci1[k],qtl,comp,ffg+ss,0,loki);
								} else {
									xx=xx1;
									ffg[ss]=ffg[ss^2];
								}
							} else {
								/* Flip paternal bit twice */
								seg[X_PAT][i]=(s2a^=1);
								ss^=1;
 								if(pass_founder_genes1(locus,i,X_PAT,loki)) {
									xx=xx1=seg_pen1(loci1[k],qtl,comp,ffg+ss,0,loki);
								} else {
									xx=xx1;
									ffg[ss]=ffg[ss^1];
								}
							}
							lks[ss][k]=xx;
							/* Get maximum value(for scaling purposes) */
							if(!ffg[ss] && xx>xx2) xx2=xx;
						}
						scl=1;
					} else if(mf==6) { /* Both segs fixed */
						for(k1=0;k1<4;k1++) pen[k1][k]=0.0;
						pen[s][k]=1.0;
						lks[s][k]=xx1;
					} else if(mf==24) { /* Both segs indetermined by data 
												* (because parents are homozygous) */
						for(k1=0;k1<4;k1++) pen[k1][k]=0.25;
						for(k1=0;k1<4;k1++) lks[k1][k]=xx1;
					} else if(mf==18) { /* Mat seg fixed, pat seg indetermined */
						s=s1<<1;
						for(k1=0;k1<2;k1++) {
							pen[s|k1][k]=0.5;
							lks[s|k1][k]=xx1;
						}
						s=(s1^1)<<1;
						for(k1=0;k1<2;k1++) pen[s|k1][k]=0.0;
					} else if(mf==12) { /* Pat seg fixed, mat seg indetermined */
						for(k1=0;k1<4;k1+=2) {
							pen[s2|k1][k]=0.5;
							lks[s2|k1][k]=xx1;
							pen[(s2^1)|k1][k]=0.0;
						}
					} else if(mf&10) { /* Don't need to vary mat seg */
						/* Current state */
						xx2=xx1=lk_store[k][comp];
						lks[s][k]=xx1;
						ffg[s]=0;
						seg[X_PAT][i]=(s2a^=1);
						ss=(s1<<1)|s2a;
						if(pass_founder_genes1(locus,i,X_PAT,loki)) xx=xx1=seg_pen1(loci1[k],qtl,comp,ffg+ss,0,loki);
						else xx=xx1;
						lks[ss][k]=xx;
						if(!ffg[ss] && xx>xx2) xx2=xx;
						k1=(s1^1)<<1;
						if(mf&2) {
							for(k2=0;k2<2;k2++) ffg[k1|k2]=1.0;
						} else {
							for(k2=0;k2<2;k2++) {
								lks[k1|k2][k]=lks[(s1<<1)|k2][k];
								ffg[k1|k2]=ffg[(s1<<1)|k2];
							}
						}
						scl=1;
					} else if(mf&20) { /* Don't need to vary pat seg */
						/* Current state */
						xx2=xx1=lk_store[k][comp];
						lks[s][k]=xx1;
						ffg[s]=0;
						seg[X_MAT][i]=(s1a^=1);
						ss=(s1a<<1)|s2;
						if(pass_founder_genes1(locus,i,X_MAT,loki)) xx=xx1=seg_pen1(loci1[k],qtl,comp,ffg+ss,0,loki);
						else xx=xx1;
						lks[ss][k]=xx;
						if(!ffg[ss] && xx>xx2) xx2=xx;
						k2=s2^1;
						if(mf&4) {
							for(k1=0;k1<4;k1+=2) ffg[k1|k2]=-1;
						} else {
							for(k1=0;k1<4;k1+=2) {
								lks[k1|k2][k]=lks[k1|s2a][k];
								ffg[k1|k2]=ffg[k1|s2a];
							}
						}
						scl=1;
					} else ABT_FUNC("Internal error - illegal flag\n");
					if(scl) {
						/* Convert from log to normal scale */
						for(z=0.0,k1=0;k1<4;k1++) {
							if(!ffg[k1]) pen[k1][k]=exp(lks[k1][k]-xx2);
							else pen[k1][k]=0.0;
							z+=pen[k1][k];
						}
						assert(z>0.0);
						z=1.0/z;
						for(k1=0;k1<4;k1++) pen[k1][k]*=z;
					}
				}
				/* Forward phase */
				for(k1=0;k1<4;k1++) pp[k1][0]=pen[k1][0];
				for(k=1;k<n_loci;k++) {
					rr[0]=recom[0][k-1];
					rr[1]=recom[1][k-1];
					xx=0.0;
					for(k1=0;k1<4;k1++) {
						if((z2=pen[k1][k])>0.0) {
							z1=0.0;
							for(k2=0;k2<4;k2++) {
								if((z=pp[k2][k-1])>0.0) {
									z*=((k1&1)!=(k2&1))?rr[1]:1.0-rr[1];
									z*=((k1&2)!=(k2&2))?rr[0]:1.0-rr[0];
									z1+=z;
								}
							}
							z2*=z1;
						}
						xx+=(pp[k1][k]=z2);
					}
					xx=1.0/xx;
					for(k1=0;k1<4;k1++) pp[k1][k]*=xx;
				}
				/* Backwards phase */
				k=n_loci-1;
				locus=(loci1[k]->type&ST_MARKER)?loci1[k]->index:-1-loci1[k]->index;
				seg=seg_list[k];
				if(!(loci1[k]->pruned_flag[i])) {
					s1=seg[X_MAT][i];
					s2=seg[X_PAT][i];
					s=(s1<<1)|s2;
					z=safe_ranf();
					for(z1=0.0,k1=0;k1<4;k1++) if(pp[k1][k]>0.0) {
						z1+=pp[k1][k];
						if(z<=z1) break;
					}
					if(k1!=s) {
						s1a=seg[X_MAT][i]=(k1&2)>>1;
						s2a=seg[X_PAT][i]=k1&1;
						if(s1a!=s1) (void)pass_founder_genes1(locus,i,X_MAT,loki);
						if(s2a!=s2) (void)pass_founder_genes1(locus,i,X_PAT,loki);
					}
					lk_store[k][comp]=lks[k1][k];
				} else {
					s1a=seg[X_MAT][i]=ranf()<.5?0:1;
					s2a=seg[X_PAT][i]=ranf()<.5?0:1;
				}
				for(k=n_loci-2;k>=0;k--) {
					rr[0]=recom[0][k];
					rr[1]=recom[1][k];
					for(xx=0.0,k1=0;k1<4;k1++) {
						if((z2=pp[k1][k])) {
							z2*=((k1&1)!=s2a)?rr[1]:1.0-rr[1];
							z2*=(((k1&2)>>1)!=s1a)?rr[0]:1.0-rr[0];
						}
						xx+=(pp[k1][k]=z2);
					}
					z=safe_ranf()*xx;
					for(z1=0.0,k1=0;k1<4;k1++) if(pp[k1][k]>0.0) {
						z1+=pp[k1][k];
						if(z<=z1) break;
					}
					s1a=(k1&2)>>1;
					s2a=k1&1;	
					locus=(loci1[k]->type&ST_MARKER)?loci1[k]->index:-1-loci1[k]->index;
					seg=seg_list[k];
					if(!(loci1[k]->pruned_flag[i])) {
						s1=seg[X_MAT][i];
						s2=seg[X_PAT][i];
						s=(s1<<1)|s2;
						if(k1!=s) {
							seg[X_MAT][i]=s1a;
							seg[X_PAT][i]=s2a;
							if(s1a!=s1) (void)pass_founder_genes1(locus,i,X_MAT,loki);
							if(s2a!=s2) (void)pass_founder_genes1(locus,i,X_PAT,loki);
						}
						lk_store[k][comp]=lks[k1][k];
					} else {
						seg[X_MAT][i]=s1a;
						seg[X_PAT][i]=s2a;
					}
				}
			}
		} else if(update_type==1) {
			ctr++;
			/* Try to flip all maternal or all paternal segs in a half sib family */
			i1=par_st[comp];
			cs=par_st[comp+1];
			zz=(double)(cs-i1);
			for(;i1<cs;i1++) {
				k=i1+(int)(safe_ranf()*(zz--)); /* Pick random parent */
				i=temp_list[k];
				temp_list[k]=temp_list[i1];
				temp_list[i1]=i;
				nkids=id_array[i].nkids;
				sex=id_array[i].sex;
				if(nkids>tmp_arr_size) {
					tmp_arr_size=nkids;
					if(!(tmp_arr=realloc_and_remember(tmp_arr,sizeof(int)*nkids))) ABT_FUNC(MMsg);
				}
				kids=tmp_arr;
				for(k=0;k<nkids;k++) kids[k]=id_array[i].kids[k]->idx;
				for(k=0;k<n_loci;k++) { /* Calculate penetrances */
					locus=(loci1[k]->type&ST_MARKER)?loci1[k]->index:-1-loci1[k]->index;
					lks[0][k]=lk_store[k][comp];
					seg=seg_list[k];
					loc_fg[k]=0;
					if(!id_array[i].sire || (locus>=0 && marker[locus].m_flag[i]&1)) {
						pen[0][k]=pen[1][k]=0.5;
						lks[1][k]=lks[0][k];
					} else {
						if(locus>=0) {
							k2=sex==1?4:2;
							for(k1=0;k1<nkids;k1++) {
								kid=kids[k1];
								if(marker[locus].m_flag[kid]&k2) break;
							}
						} else k1=nkids;
						if(k1<nkids) { /* If one kid is fixed, can't change */
							pen[0][k]=1.0;
							pen[1][k]=0.0;
						} else {
							k2=2-sex;
							for(k1=0;k1<nkids;k1++) {
								kid=kids[k1];
								seg[k2][kid]^=1;
							}
							loc_fg[k]=1; /* Flag that this locus has been flipped */
							j=pass_founder_genes1a(locus,kids,nkids,k2,loki);
							if(j) {
								lks[1][k]=seg_pen1(loci1[k],qtl,comp,&k2,0,loki);
								if(!k2) {
									z=lks[0][k];
									z1=lks[1][k];
									if(z>z1) {
										z1=exp(z1-z);
										z=1.0;
									} else {
										z=exp(z-z1);
										z1=1.0;
									}
									z2=z+z1;
									pen[0][k]=z/z2;
									pen[1][k]=z1/z2;
								} else {
									pen[0][k]=1.0;
									pen[1][k]=0.0;
								}
							} else {
								lks[1][k]=lks[0][k];
								pen[0][k]=pen[1][k]=0.5;
							}
						}
					}
				}
				/* Forward phase */
				k2=2-sex;
				for(k1=0;k1<2;k1++) pp[k1][0]=pen[k1][0];
				seg=seg_list[0];
				for(k=1;k<n_loci;k++) {
					seg1=seg;
					seg=seg_list[k];
					for(k1=0;k1<2;k1++) pp[k1][k]=pen[k1][k];
					if(pen[1][k]>0.0) {
						rr[1]=recom[k2][k-1];
						rr[0]=1.0-rr[1];
						z1=z2=1.0;
						for(k1=0;k1<nkids;k1++) {
							kid=kids[k1];
							s1=seg1[k2][kid]^loc_fg[k-1];
							s2=seg[k2][kid]^loc_fg[k];
							ss=s1^s2;
							z1*=rr[ss];
							z2*=rr[ss^1];
						}
						pp[0][k]*=(pp[0][k-1]*z1+pp[1][k-1]*z2);
						pp[1][k]*=(pp[0][k-1]*z2+pp[1][k-1]*z1);
						z=pp[0][k]+pp[1][k];
						pp[0][k]/=z;
						pp[1][k]/=z;
					}
				}
				/* Backwards phase */
				k=n_loci-1;
				seg=seg_list[k];
				k1=(safe_ranf()<=pp[0][k])?0:1;
				if(k1!=loc_fg[k]) {
					locus=(loci1[k]->type&ST_MARKER)?loci1[k]->index:-1-loci1[k]->index;
					for(j=0;j<nkids;j++) {
						kid=kids[j];
						seg[k2][kid]^=1;
					}
					pass_founder_genes1a(locus,kids,nkids,k2,loki);
				}
				lk_store[k][comp]=lks[k1][k];
				for(k=n_loci-2;k>=0;k--) {
					seg1=seg;
					seg=seg_list[k];
					if(pp[1][k]>0.0) {
						rr[1]=recom[k2][k];
						rr[0]=1.0-rr[1];
						z1=z2=1.0;
						for(k1=0;k1<nkids;k1++) {
							kid=kids[k1];
							s1=seg1[k2][kid];
							s2=seg[k2][kid]^loc_fg[k];
							ss=s1^s2;
							z1*=rr[ss];
							z2*=rr[ss^1];
						}
						pp[0][k]*=z1;
						pp[1][k]*=z2;
						z=pp[0][k]+pp[1][k];
						k1=(ranf()*z<=pp[0][k])?0:1;
					} else k1=0;
					if(k1!=loc_fg[k]) {
						locus=(loci1[k]->type&ST_MARKER)?loci1[k]->index:-1-loci1[k]->index;
						for(j=0;j<nkids;j++) {
							kid=kids[j];
							seg[k2][kid]^=1;
						}
						(void)pass_founder_genes1a(locus,kids,nkids,k2,loki);
					}
					lk_store[k][comp]=lks[k1][k];
				}
			}
		} else if(update_type==2) {
			/* Try to flip segs in a 3 generation family */
			i1=gpfam_st[comp];
			cs=gpfam_st[comp+1];
			zz=(double)(cs-i1);
			for(;i1<cs;i1++) {
				k=i1+(int)(safe_ranf()*(zz--));
				i=temp_list[k];
				temp_list[k]=temp_list[i1];				
				temp_list[i1]=i;
				nkids=families[i].nkids;
				kids=families[i].kids;
				ids=id_array[kids[0]].sire-1;
				idd=id_array[kids[0]].dam-1;
				symflag=0;
				if(!(id_array[ids].sire || id_array[idd].sire) &&
					(id_array[ids].nkids==nkids && id_array[idd].nkids==nkids)) symflag=1;
				for(k=0;k<n_loci;k++) { /* Calculate penetrances */
					locus=(loci1[k]->type&ST_MARKER)?loci1[k]->index:-1-loci1[k]->index;
					lks[0][k]=lk_store[k][comp];
					seg=seg_list[k];
					loc_fg[k]=0;
					k1=0;
					if(symflag) {
						if(locus>=0) {
							if(!marker[locus].haplo[ids] && !marker[locus].haplo[idd]) {
								if(marker[locus].mterm && marker[locus].mterm[0]) {
									if(!(id_array[ids].res[0] || id_array[idd].res[0])) k1=1;
								} else k1=1;
							}
						} else {
							if(!(id_array[ids].res[0] || id_array[idd].res[0])) k1=1;
						}
					}
					if(k1) {
						pen[0][k]=pen[1][k]=0.5;
						lks[1][k]=lks[0][k];
						continue;
					}
					if(locus>=0) {
						/* Check if we can switch at this locus */
						for(k1=0;k1<nkids;k1++) {
							/* Check kids */
							kid=kids[k1];
							/* If either of the kid's seg are fixed, we can only
							 * perform the move if both segs are the same */
							if((marker[locus].m_flag[kid]&6)&&(seg[0][kid]!=seg[1][kid])) break;
							/* Check grandkids */
							nkids1=id_array[kid].nkids;
							kids1=id_array[kid].kids;
							k2=id_array[kid].sex==1?4:2;
							for(k3=0;k3<nkids1;k3++) {
								kid1=kids1[k3]->idx;
								if(marker[locus].m_flag[kid1]&k2) break;
							}
							if(k3<nkids1) break;
						} 
					} else k1=nkids;
					if(k1<nkids) { /* Can't switch at this locus */
						pen[0][k]=1.0;
						pen[1][k]=0.0;
					} else {
						for(k1=0;k1<nkids;k1++) {
							kid=kids[k1];
							/* Swap segs of kids */
							k3=seg[0][kid];
							seg[0][kid]=seg[1][kid];
							seg[1][kid]=k3;
							/* Flip segs of grandkids */
							nkids1=id_array[kid].nkids;
							kids1=id_array[kid].kids;
							k2=2-id_array[kid].sex;
							for(k3=0;k3<nkids1;k3++) {
								kid1=kids1[k3]->idx;
								seg[k2][kid1]^=1;
							}
						}
						loc_fg[k]=1; /* Flag that this locus has been flipped */
						j=pass_founder_genes1b(locus,kids,nkids,loki);
						if(j) {
							lks[1][k]=seg_pen1(loci1[k],qtl,comp,&k2,0,loki);
							if(!k2) {
								z=lks[0][k];
								z1=lks[1][k];
								if(z>z1) {
									z1=exp(z1-z);
									z=1.0;
								} else {
									z=exp(z-z1);
									z1=1.0;
								}
								z2=z+z1;
								pen[0][k]=z/z2;
								pen[1][k]=z1/z2;
							} else {
								pen[0][k]=1.0;
								pen[1][k]=0.0;
							}
						} else {
							lks[1][k]=lks[0][k];
							pen[0][k]=pen[1][k]=0.5;
						}
					}
				}
				/* Forward phase */
				for(k1=0;k1<2;k1++) pp[k1][0]=pen[k1][0];
				seg=seg_list[0];
				for(k=1;k<n_loci;k++) {
					seg1=seg;
					seg=seg_list[k];
					for(k1=0;k1<2;k1++) pp[k1][k]=pen[k1][k];
					if(pen[1][k]>0.0) {
						rr[0]=recom[0][k-1];
						rr[1]=recom[1][k-1];
						z1=z2=1.0;
						for(k1=0;k1<nkids;k1++) {
							kid=kids[k1];
							/* Handle swaps in kids' segs */
							s2=seg[loc_fg[k]][kid]; /* Maternal seg of current locus if no swap occurs */
							s3=seg[loc_fg[k]^1][kid]; /* Maternal seg of current locus if swap occurs */
							if(s2!=s3) {
								s1=seg1[loc_fg[k-1]][kid]; /* Maternal seg of previous locus */
								z1*=(s1==s2)?1.0-rr[0]:rr[0]; /* Prob. of swap */
								z2*=(s1==s3)?1.0-rr[0]:rr[0]; /* Prob. of no swap */
								s1=seg1[loc_fg[k-1]^1][kid]; /* Paternal seg of previous locus */
								z1*=(s1==s3)?1.0-rr[1]:rr[1]; /* Prob. of no swap */
								z2*=(s1==s2)?1.0-rr[1]:rr[1]; /* Prob. of swap */
							}
							/* Handle flips in grandkids */
							nkids1=id_array[kid].nkids;
							kids1=id_array[kid].kids;
							k2=2-id_array[kid].sex;
							for(k3=0;k3<nkids1;k3++) {
								kid1=kids1[k3]->idx;
								s1=seg1[k2][kid1]^loc_fg[k-1];
								s2=seg[k2][kid1]^loc_fg[k];
								if(s1==s2) {
									z1*=1.0-rr[k2];
									z2*=rr[k2];
								} else {
									z2*=1.0-rr[k2];
									z1*=rr[k2];
								}
							}
						}
						pp[0][k]*=(pp[0][k-1]*z1+pp[1][k-1]*z2);
						pp[1][k]*=(pp[0][k-1]*z2+pp[1][k-1]*z1);
						z=pp[0][k]+pp[1][k];
						pp[0][k]/=z;
						pp[1][k]/=z;
					}
				}
				/* Backwards phase */
				k=n_loci-1;
				locus=(loci1[k]->type&ST_MARKER)?loci1[k]->index:-1-loci1[k]->index;
				seg=seg_list[k];
				k1=(ranf()<=pp[0][k])?0:1;
				if(k1!=loc_fg[k]) {
					for(j=0;j<nkids;j++) {
						kid=kids[j];
						k2=seg[0][kid];
						seg[0][kid]=seg[1][kid];
						seg[1][kid]=k2;
						nkids1=id_array[kid].nkids;
						kids1=id_array[kid].kids;
						k2=2-id_array[kid].sex;
						for(k3=0;k3<nkids1;k3++) {
							kid1=kids1[k3]->idx;
							seg[k2][kid1]^=1;
						}
					}
					(void)pass_founder_genes1b(locus,kids,nkids,loki);
				}
				lk_store[k][comp]=lks[k1][k];
				for(k=n_loci-2;k>=0;k--) {
					locus=(loci1[k]->type&ST_MARKER)?loci1[k]->index:-1-loci1[k]->index;
					seg1=seg;
					seg=seg_list[k];
					if(pp[1][k]>0.0) {
						rr[0]=recom[0][k];
						rr[1]=recom[1][k];
						z1=z2=1.0;
						for(k1=0;k1<nkids;k1++) {
							kid=kids[k1];
							s2=seg[loc_fg[k]][kid];
							s3=seg[loc_fg[k]^1][kid];
							if(s2!=s3) {
								s1=seg1[0][kid];
								z1*=(s1==s2)?1.0-rr[0]:rr[0];
								z2*=(s1==s3)?1.0-rr[0]:rr[0];
								s1=seg1[1][kid];
								z1*=(s1==s3)?1.0-rr[1]:rr[1];
								z2*=(s1==s2)?1.0-rr[1]:rr[1];
							}
							nkids1=id_array[kid].nkids;
							kids1=id_array[kid].kids;
							k2=2-id_array[kid].sex;
							for(k3=0;k3<nkids1;k3++) {
								kid1=kids1[k3]->idx;
								s1=seg1[k2][kid1];
								s2=seg[k2][kid1]^loc_fg[k];
								if(s1==s2) {
									z1*=1.0-rr[k2];
									z2*=rr[k2];
								} else {
									z2*=1.0-rr[k2];
									z1*=rr[k2];
								}
							}
						}
						pp[0][k]*=z1;
						pp[1][k]*=z2;
						z=pp[0][k]+pp[1][k];
						k1=(ranf()*z<=pp[0][k])?0:1;
					} else k1=0;
					if(k1!=loc_fg[k]) {
						for(j=0;j<nkids;j++) {
							kid=kids[j];
							k2=seg[0][kid];
							seg[0][kid]=seg[1][kid];
							seg[1][kid]=k2;
							nkids1=id_array[kid].nkids;
							kids1=id_array[kid].kids;
							k2=2-id_array[kid].sex;
							for(k3=0;k3<nkids1;k3++) {
								kid1=kids1[k3]->idx;
								seg[k2][kid1]^=1;
							}
						}
						(void)pass_founder_genes1b(locus,kids,nkids,loki);
					}
					lk_store[k][comp]=lks[k1][k];
				}
			} 
		} else if(update_type==3) {
			/* Try to flip maternal and paternal segs in a nuclear family */
			i1=fam_st[comp];
			cs=fam_st[comp+1];
			zz=(double)(cs-i1);
			for(;i1<cs;i1++) {
				k=i1+(int)(safe_ranf()*(zz--));
				i=temp_list[k];
				temp_list[k]=temp_list[i1];
				temp_list[i1]=i;
				nkids=families[i].nkids;
				kids=families[i].kids;
				kid=kids[0];
				ids=id_array[kid].sire;
				idd=id_array[kid].dam;
				for(k=0;k<n_loci;k++) { /* Calculate penetrances */
					locus=(loci1[k]->type&ST_MARKER)?loci1[k]->index:-1-loci1[k]->index;
					seg=seg_list[k];
					loc_fg[k]=0;
					if(locus>=0) {
						/* Check if we can switch at this locus */
						/* If any kid has a seg fixed, then we can't switch that seg */
						for(mf=0,k1=0;k1<nkids;k1++) {
							/* Check kids */
							kid=kids[k1];
							mf|=(marker[locus].m_flag[kid]&6);
						}
						if(marker[locus].m_flag[idd-1]&1) mf|=8;
						if(marker[locus].m_flag[ids-1]&1) mf|=16;
					} else mf=0;
					if(!id_array[idd-1].sire && id_array[idd-1].nkids==nkids) mf|=8;
					if(!id_array[ids-1].sire && id_array[ids-1].nkids==nkids) mf|=16;
					lffg=ffg[0]=0;
					xx=xx2=lks[0][k]=lk_store[k][comp];
					loc_fg[k]=state=0;
					for(kk=3,k1=0;k1<kk;k1++) {
						if(k1&1) { /* Flip maternal segs */
							state^=2;
							if(!(mf&10)) { 
								for(k3=0;k3<nkids;k3++) { /* Flip all maternal segs in the kids */
									kid=kids[k3];
									seg[X_MAT][kid]^=1;
								}
								k2=pass_founder_genes1a(locus,kids,nkids,X_MAT,loki);
								loc_fg[k]^=2; /* Update state for this locus */
								if(k2) {
									xx=lks[state][k]=seg_pen1(loci1[k],qtl,comp,ffg+state,0,loki);
									lffg=ffg[state];
									if(!ffg[state] && xx>xx2) xx2=xx;
								} else { /* Nothing changed */
									lks[state][k]=xx;
									ffg[state]=lffg;
								}
							} else if(mf&2) { /* Maternal seg can not be switched */
								lks[state][k]=lks[state^1][k]=16;
								ffg[state]=ffg[state^1]=-1;
								kk=2;
							} else if(mf&8) { /* Maternal seg has equal likelihood if switched */
								lks[state][k]=lks[state^2][k];
								ffg[state]=ffg[state^2];
							}
						} else { /* Flip paternal segs */
							state^=1;
							if(!(mf&20)) {
								for(k3=0;k3<nkids;k3++) { /* Flip all paternal segs in the kids */
									kid=kids[k3];
									seg[X_PAT][kid]^=1;
								}
								k2=pass_founder_genes1a(locus,kids,nkids,X_PAT,loki);
								loc_fg[k]^=1; /* Update state for this locus */
								if(k2) {
									xx=lks[state][k]=seg_pen1(loci1[k],qtl,comp,ffg+state,0,loki);
									lffg=ffg[state];
									if(!ffg[state] && xx>xx2) xx2=xx;
								} else {
									lks[state][k]=xx;
									ffg[state]=lffg;
								}
							} else if(mf&4) { /* Paternal seg can not be switched */
								lks[state][k]=lks[state^2][k]=17;
								ffg[state]=ffg[state^2]=-1;
								state^=1;
								kk=2;
							} else if(mf&16) { /* Paternal seg has equal likelihood if switched */
								lks[state][k]=lks[state^1][k];
								ffg[state]=ffg[state^1];
							}
						}
					}
					/* Convert from log to normal scale */
					for(z=0.0,k1=0;k1<4;k1++) {
						if(!ffg[k1]) pen[k1][k]=exp(lks[k1][k]-xx2);
						else pen[k1][k]=0.0;
						z+=pen[k1][k];
					}
					if(!(z>0.0)) {
						printf("ZERO: %g\n",z);
						for(k1=0;k1<4;k1++) {
							printf("%d ",ffg[k1]);
							if(!ffg[k1]) {
								printf("%g %g\n",lks[k1][k],pen[k1][k]);
							}
							printf("\n");
						}
						assert(z>0.0);
					}
					z=1.0/z;
					for(k1=0;k1<4;k1++) pen[k1][k]*=z;
				}
				/* Forward phase */
				for(k1=0;k1<4;k1++) pp[k1][0]=pen[k1][0];
				seg=seg_list[0];
				for(k=1;k<n_loci;k++) {
					seg1=seg;
					seg=seg_list[k];
					rr[0]=recom[0][k-1]; /* Male recom  fraction */
					rr[1]=recom[1][k-1]; /* Female recom fraction */
					xx=0.0;
					for(k1=0;k1<4;k1++) {
						if((z2=pen[k1][k])>0.0) {
							for(k2=0;k2<4;k2++) {
								zz1[k2]=pp[k2][k-1];
								if(zz1[k2]>0.0) {
									for(z=1.0,k3=0;k3<nkids;k3++) {
										kid=kids[k3];
										s1=seg[X_MAT][kid];
										s2=seg[X_PAT][kid];
										mf=loc_fg[k]^k1;
										if(mf&1) s2^=1;
										if(mf&2) s1^=1;
										s1a=seg1[X_MAT][kid];
										s2a=seg1[X_PAT][kid];
										mf=loc_fg[k-1]^k2;
										if(mf&1) s2a^=1;
										if(mf&2) s1a^=1;
										z*=(s1a==s1)?1.0-rr[X_MAT]:rr[X_MAT];
										z*=(s2a==s2)?1.0-rr[X_PAT]:rr[X_PAT];
									}
									zz1[k2]*=z;
								}
							}
							for(z1=0.0,k2=0;k2<4;k2++) z1+=zz1[k2];
							z2*=z1;
						}
						xx+=(pp[k1][k]=z2);
					}
					xx=1.0/xx;
					for(k1=0;k1<4;k1++) pp[k1][k]*=xx;
				}
				/* Backwards phase */	
				k=n_loci-1;
				locus=(loci1[k]->type&ST_MARKER)?loci1[k]->index:-1-loci1[k]->index;
				seg=seg_list[k];
				z=safe_ranf();
				for(z1=0.0,k1=0;k1<4;k1++) if(pp[k1][k]>0.0) {
					z1+=pp[k1][k];
					if(z<=z1) break;
				}
				if(k1!=loc_fg[k]) {
					mf=loc_fg[k]^k1;
					if(mf) {
						for(k3=0;k3<nkids;k3++) {
							kid=kids[k3];
							if(mf&1) seg[X_PAT][kid]^=1;
							if(mf&2) seg[X_MAT][kid]^=1;
						}
						if(mf&1) (void)pass_founder_genes1a(locus,kids,nkids,X_PAT,loki);
						if(mf&2) (void)pass_founder_genes1a(locus,kids,nkids,X_MAT,loki);
					}
				}
				lk_store[k][comp]=lks[k1][k];
				for(k=n_loci-2;k>=0;k--) {
					seg1=seg;
					locus=(loci1[k]->type&ST_MARKER)?loci1[k]->index:-1-loci1[k]->index;
					seg=seg_list[k];
					rr[0]=recom[0][k];
					rr[1]=recom[1][k];
					for(xx=0.0,k1=0;k1<4;k1++) {
						if((z2=pp[k1][k])) {
							mf=loc_fg[k]^k1;
							for(k3=0;k3<nkids;k3++) {
								kid=kids[k3];
								s1a=seg1[X_MAT][kid];
								s2a=seg1[X_PAT][kid];
								s1=seg[X_MAT][kid];
								s2=seg[X_PAT][kid];
								if(mf&1) s2^=1;
								if(mf&2) s1^=1;
								z2*=(s1==s1a)?1.0-rr[X_MAT]:rr[X_MAT];
								z2*=(s2==s2a)?1.0-rr[X_PAT]:rr[X_PAT];
							}
						}
						xx+=(pp[k1][k]=z2);
					}
					z=safe_ranf()*xx;
					for(z1=0.0,k1=0;k1<4;k1++) if(pp[k1][k]>0.0) {
						z1+=pp[k1][k];
						if(z<=z1) break;
					}
					if(k1!=loc_fg[k]) {
						mf=loc_fg[k]^k1;
						if(mf) {
							for(k3=0;k3<nkids;k3++) {
								kid=kids[k3];
								if(mf&1) seg[X_PAT][kid]^=1;
								if(mf&2) seg[X_MAT][kid]^=1;
							}
							if(mf&1) (void)pass_founder_genes1a(locus,kids,nkids,X_PAT,loki);
							if(mf&2) (void)pass_founder_genes1a(locus,kids,nkids,X_MAT,loki);
						}
					}
					lk_store[k][comp]=lks[k1][k];
				}
			}
		}
		if(qtl) {
			z=gen_pen(qtl,comp,&i,1,loki);
#ifndef NDEBUG
			if(i) {
				(void)fprintf(stderr,"%d %d %d %d\n",k,comp,i,(int)z);
				ABT_FUNC("Illegal segregation pattern\n");
			}
			if(fabs(z-qtl->lk_store[comp])>1.0e-8) {
				(void)fprintf(stderr,"Seg. prob mismatch %g %g (comp=%d, locus=%d, type=%d, n_qtl=%d)\n",z,loki->models->tlocus[k].lk_store[comp],comp,k,update_type,n_qtl);
				ABT_FUNC("Aborting\n");
			}
#endif
			for(k2=0;k2<loki->params.n_tloci;k2++) if(loki->models->tlocus[k2].flag && k2!=k) {
				z=gen_pen(loki->models->tlocus+k2,comp,&i,0,loki);
#ifndef NDEBUG
				if(i) {
					(void)fprintf(stderr,"%d %d %d %d\n",k2,comp,i,(int)z);
					ABT_FUNC("Illegal segregation pattern\n");
				}
#endif
				loki->models->tlocus[k2].lk_store[comp]=z;
			}
		}
	}
}
