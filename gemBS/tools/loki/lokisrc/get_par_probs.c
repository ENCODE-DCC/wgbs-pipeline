/****************************************************************************
 *                                                                          *
 *     Loki - Programs for genetic analysis of complex traits using MCMC    *
 *                                                                          *
 *             Simon Heath - CNG                                            *
 *                                                                          *
 *                        April 2002                                        *
 *                                                                          *
 * get_par_probs.c:                                                         *
 *                                                                          *
 * Get parental probabilities (for peeling)                                 *
 *                                                                          *
 * Copyright (C) Simon C. Heath 1997, 2000, 2002                            *
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
#include <assert.h>

#include "utils.h"
#include "loki.h"
#include "loki_peel.h"
#include "get_par_probs.h"

/* Get parental probability dist., put into *val */
double get_par_probs(double *val,const int i,struct Marker *mark,pen_func pen,lk_ulong **a_set,double **freq,
							struct R_Func *rf,struct loki *loki)
{
	int j,k,k1,fflag=0,nb1,n_idx,comp,n_all,n_bits;
	double f,f1,lf[2],*fq=0,*tmp,p=0.0,z;
	lk_ulong mask,cm[2],l,l1,l2,a;
	struct Id_Record *id_array;

#ifdef TRACE_PEEL
		if(CHK_PEEL(TRACE_LEVEL_2)) (void)printf("In %s(%p,%d,%s,%p)\n",__func__,(void *)val,i,mark->name,(void *)pen);
#endif
	id_array=loki->pedigree->id_array;
	comp=id_array[i].comp;
	n_all=mark->n_all1[comp];
	n_bits=num_bits(n_all);
 	nb1=1<<n_bits;
	mask=nb1-1;
	n_idx=1<<(n_bits+n_bits);
 	if(!mark->locus.founder_flag[i]) fflag=1;
	else {
		k=id_array[i].group-1;
		assert(k>=0);
		fq=freq[k];
		/* If found add founder probs */
#ifdef TRACE_PEEL
		if(CHK_PEEL(TRACE_LEVEL_2)) {
			(void)fputs("get_par_probs(): inserting founder probs. for ",stdout);
			print_orig_id(stdout,i+1);
			(void)fputc('\n',stdout);
			if(CHK_PEEL(TRACE_LEVEL_3)) {
				for(j=0;j<n_all;j++) (void)printf("%g ",fq[j]);
				(void)fputc('\n',stdout);
			}
		}
#endif
		for(j=0;j<2;j++) {
			lf[j]=0.0;
			cm[j]=mark->req_set[j][i];
			if(cm[j]) {
				a=cm[j];
				k=0;
				while(a) {
					if(a&1) {
						lf[j]+=fq[k];
					}
					k++;
					a>>=1;
				}
			}
		}
	}
	if(pen) {
		tmp=val;
		for(j=0;j<n_idx;j++) *(tmp++)=0.0;
	}
	if((k=id_array[i].rfp)>=0) { /* Insert Previously computed R_Func */
		if(!pen) for(j=0;j<n_all;j++) {
			a=a_set[i][j];
			tmp=val+j;
			while(a) {
				if(a&1) *tmp=0.0;
				tmp+=nb1;
				a>>=1;
			}
		}
#ifdef TRACE_PEEL
		if(CHK_PEEL(TRACE_LEVEL_2)) {
			(void)fputs("get_par_probs(): inserting previously computed R-Function for ",stdout);
			print_orig_id(stdout,i+1);
			(void)fputc('\n',stdout);
		}
#endif
		k1=rf[k].n_terms;
		if(fflag) {
			for(j=0;j<k1;j++) {
				z=rf[k].p[j];
				val[rf[k].index[j]]=z;
				p+=z;
			}
		} else {
			for(j=0;j<k1;j++) {
				l=rf[k].index[j];
				l1=l&mask;
				l2=(l>>n_bits)&mask;
				if((cm[X_MAT])&(1<<l1)) f=lf[X_MAT];
				else f=fq[l1];
				if((cm[X_PAT])&(1<<l2)) f1=lf[X_PAT];
				else f1=fq[l2];
				val[l]=rf[k].p[j]*f*f1;
				p+=val[l];
			}
		}
	} else { /* No R-Function - remove illegal genotypes */
		if(fflag) {
			for(j=0;j<n_all;j++) {
				a=a_set[i][j];
				tmp=val+j;
				while(a) {
					if(a&1) {
						*tmp=1.0;
						p+=1.0;
					}
					tmp+=nb1;
					a>>=1;
				}
			}
		} else {
			for(j=0;j<n_all;j++)	{
				a=a_set[i][j];
				if(a) {
					l=j;
					if((cm[X_MAT])&(1<<j)) f=lf[X_MAT];
					else f=fq[j];
					k=0;
					while(a) {
						if(a&1) {
							if((cm[X_PAT])&(1<<k)) f1=lf[X_PAT];
							else f1=fq[k];
/*							printf("l=%d, f=%g, f1=%g (%d,%d) (%g,%g)\n",(int)l,f,f1,j,k,fq[j],fq[k]); */
							val[l]=f*f1;
							p+=f*f1;
						}
						l+=nb1;
						k++;
						a>>=1;
					}
				}
			}
		}
	}
	/* Add penetrances */
	if(pen) {
		pen(val,i,&mark->locus,n_all,n_bits,loki);
		p=0.0;
		for(j=0;j<n_all;j++) {
			a=a_set[i][j];
			l=j;
			while(a) {
				if(a&1) p+=val[l];
				a>>=1;
				l+=nb1;
			}
		}
	}
	assert(p>0.0);
	z=1.0/p;
	for(j=0;j<n_all;j++) {
		a=a_set[i][j];
		l=j;
		while(a) {
			if(a&1) val[l]*=z;
			a>>=1;
			l+=nb1;
		}
	}
#ifdef TRACE_PEEL
	if(CHK_PEEL(TRACE_LEVEL_2)) {
		(void)printf("get_par_probs(): Returning with %g for ",p);
		print_orig_id(stdout,i+1);
		(void)fputc('\n',stdout);
	}
#endif
	return p;
}

/* Get parental probability dist. for x linked markers, put into *val */
double get_par_probs_x(double *val,const int i,struct Marker *mark,pen_func pen,lk_ulong **a_set,double **freq,
							struct R_Func *rf,struct loki *loki)
{
	int j,k,k1,fflag=0,nb1,sex,n_idx,comp,n_bits,n_all;
	double f,f1,lf[2],*fq=0,*tmp,p=0.0,z;
	lk_ulong mask,cm[2],l,l1,l2,a;
	struct Id_Record *id_array;
	
#ifdef TRACE_PEEL
		if(CHK_PEEL(TRACE_LEVEL_2)) (void)printf("In %s(%p,%d,%s,%p)\n",__func__,(void *)val,i,mark->name,(void *)pen);
#endif
	id_array=loki->pedigree->id_array;
	comp=id_array[i].comp;
	n_all=mark->n_all1[comp];
	n_bits=num_bits(n_all);
 	nb1=1<<n_bits;
	mask=nb1-1;
	n_idx=1<<(n_bits+n_bits);
	sex=id_array[i].sex;
 	if(!mark->locus.founder_flag[i]) fflag=1;
	else {
		k=id_array[i].group-1;
		assert(k>=0);
		fq=freq[k];
		/* If found add founder probs */
#ifdef TRACE_PEEL
		if(CHK_PEEL(TRACE_LEVEL_2)) {
			(void)fputs("get_par_probs(): inserting founder probs. for ",stdout);
			print_orig_id(stdout,i+1);
			(void)fputc('\n',stdout);
			if(CHK_PEEL(TRACE_LEVEL_3)) {
				for(j=0;j<n_all;j++) (void)printf("%g ",fq[j]);
				(void)fputc('\n',stdout);
			}
		}
#endif
		for(j=0;j<sex;j++) {
			lf[j]=0.0;
			cm[j]=mark->req_set[j][i];
			if(cm[j]) {
				a=cm[j];
				k=0;
				while(a) {
					if(a&1) lf[j]+=fq[k];
					k++;
					a>>=1;
				}
			}
		}
	}
	if(pen) {
		tmp=val;
		for(j=0;j<n_idx;j++) *(tmp++)=0.0;
	}
	if((k=id_array[i].rfp)>=0) { /* Insert Previously computed R_Func */
		if(!pen) {
			if(sex==1) {
				for(j=0;j<n_all;j++) {
					a=a_set[i][j];
					if(a) val[j]=0.0;
				}
			} else {
				for(j=0;j<n_all;j++) {
					a=a_set[i][j];
					tmp=val+j;
					while(a) {
						if(a&1) *tmp=0.0;
						tmp+=nb1;
						a>>=1;
					}
				}
			}
		}
#ifdef TRACE_PEEL
		if(CHK_PEEL(TRACE_LEVEL_2)) {
			(void)fputs("get_par_probs(): inserting previously computed R-Function for ",stdout);
			print_orig_id(stdout,i+1);
			(void)fputc('\n',stdout);
		}
#endif
		k1=rf[k].n_terms;
		if(fflag) {
			for(j=0;j<k1;j++) {
				z=rf[k].p[j];
				val[rf[k].index[j]]=z;
				p+=z;
			}
		} else {
			for(j=0;j<k1;j++) {
				l=rf[k].index[j];
				l1=l&mask;
				l2=(l>>n_bits)&mask;
				if((cm[X_MAT])&(1<<l1)) f=lf[X_MAT];
				else f=fq[l1];
				if(sex==2) {
					if((cm[X_PAT])&(1<<l2)) f1=lf[X_PAT];
					else f1=fq[l2];
				} else f1=1.0;
				val[l]=rf[k].p[j]*f*f1;
				p+=val[l];
			}
		}
	} else { /* No R-Function - remove illegal genotypes */
		if(fflag) {
			if(sex==1) {
				for(j=0;j<n_all;j++) {
					if(a_set[i][j]) {
						val[j]=1.0;
						p+=1.0;
					}
				}
			} else {
				for(j=0;j<n_all;j++) {
					a=a_set[i][j];
					tmp=val+j;
					while(a) {
						if(a&1) {
							*tmp=1.0;
							p+=1.0;
						}
						tmp+=nb1;
						a>>=1;
					}
				}
			}
		} else {
			if(sex==2) {
				for(j=0;j<n_all;j++)	{
					a=a_set[i][j];
					if(a) {
						l=j;
						if((cm[X_MAT])&(1<<j)) f=lf[X_MAT];
						else f=fq[j];
						k=0;
						while(a) {
							if(a&1) {
								if((cm[X_PAT])&(1<<k)) f1=lf[X_PAT];
								else f1=fq[k];
								val[l]=f*f1;
								p+=f*f1;
							}
							l+=nb1;
							k++;
							a>>=1;
						}
					}
				}
			} else {
				for(j=0;j<n_all;j++)	{
					if(a_set[i][j]) {
						if((cm[X_MAT])&(1<<j)) f=lf[X_MAT];
						else f=fq[j];
						val[j]=f;
						p+=f;
					}
				}
			}
		}
	}
	/* Add penetrances */
	if(pen) {
		pen(val,i,&mark->locus,n_all,n_bits,loki);
		p=0.0;
		if(sex==1) {
			for(j=0;j<n_all;j++) {
				if(a_set[i][j]) p+=val[j];
			}
			z=1.0/p;
			for(j=0;j<n_all;j++) {
				if(a_set[i][j]) val[j]*=z;
			}
		} else {
			for(j=0;j<n_all;j++) {
				a=a_set[i][j];
				l=j;
				while(a) {
					if(a&1) p+=val[l];
					a>>=1;
					l+=nb1;
				}
			}
			z=1.0/p;
			for(j=0;j<n_all;j++) {
				a=a_set[i][j];
				l=j;
				while(a) {
					if(a&1) val[l]*=z;
					a>>=1;
					l+=nb1;
				}
			}
		}
	}
#ifdef TRACE_PEEL
	if(CHK_PEEL(TRACE_LEVEL_2)) {
		(void)printf("get_par_probs(): Returning with %g for ",p);
		print_orig_id(stdout,i+1);
		(void)fputc('\n',stdout);
	}
#endif
	return p;
}

/* Similar to get_par_probs, but for a trait locus */
double get_trait_par_probs(double *val,const int i,const int locus,trait_pen_func *pen,double **freq,struct R_Func *rf,struct loki *loki)
{
	int j,k,fflag=0,l1,l2,n_idx,n_all,n_markers;
	double f,*fq=0,z,z1,*tmp,*tmp1,*tmp2;
	struct Id_Record *id_array;
	struct Locus *loc;
	
#ifdef TRACE_PEEL
	if(CHK_PEEL(TRACE_LEVEL_2)) (void)printf("In %s(%p,%d,%d)\n",__func__,(void *)val,i,locus);
#endif
	id_array=loki->pedigree->id_array;
	n_markers=loki->markers->n_markers;
	loc=&loki->models->tlocus[-1-locus];
	n_all=loc->n_alleles;
	n_idx=n_all*n_all;
	j=id_array[i].sire;
	if(j && !(loki->models->pruned_flag[j-1])) fflag=1;
	else {
		k=id_array[i].group-1;
		assert(k>=0);
		fq=freq[k];
	}
	tmp=val;
	if((k=id_array[i].rfp)>=0) {/* Insert Previously computed R_Func */
		tmp1=rf[k].p;
		if(fflag) for(j=0;j<n_idx;j++) *(tmp++)=*(tmp1++);
		else {
			for(j=l1=0;l1<n_all;l1++) {
				f=fq[l1];
				tmp2=fq;
				for(l2=0;l2<n_all;l2++) *(tmp++)=*(tmp1++)*(*(tmp2++))*f;
			}
		}
	} else {
		if(fflag) {
			for(j=0;j<n_idx;j++) *(tmp++)=1.0;
		} else {
			for(l1=0;l1<n_all;l1++) {
				f=fq[l1];
				tmp1=fq;
				for(l2=0;l2<n_all;l2++)	*(tmp++)=*(tmp1++)*f;
			}
		}
	}
	if(id_array[i].res[0]) pen(val,i,loc,loki);
	tmp=val;
	z=0.0;
	for(j=0;j<n_idx;j++) z+=*(tmp++);
	if(z>0.0)	{
		z1=1.0/z;
		tmp=val;
		for(j=0;j<n_idx;j++) *(tmp++)*=z1;
	}
	return z;
}

