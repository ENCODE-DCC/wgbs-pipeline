/****************************************************************************
*                                                                          *
*     Loki - Programs for genetic analysis of complex traits using MCMC    *
*                                                                          *
*             Simon Heath - University of Washington                       *
*                         - Rockefeller University                         *
*                                                                          *
*                        August 1997                                       *
*                                                                          *
* loki_peel.c:                                                             *
*                                                                          *
* Perform peeling calculations                                             *
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
#include <float.h>
#include <assert.h>

#ifndef DBL_MAX
#define DBL_MAX MAXDOUBLE
#endif

#include "ranlib.h"
#include "utils.h"
#include "loki.h"
#include "loki_peel.h"
#include "loki_utils.h"
#include "seg_pen.h"
#include "get_par_probs.h"
#include "lk_malloc.h"
#include "loki_simple_peel.h"
#include "loki_trait_simple_peel.h"

static double **freq;
static struct Peelseq_Head **peel_list;
static int *peel_alls;
static struct loki *lk;

static void peel_dealloc(void)
{
	int k;
	struct peel_mem *work;
	struct peel_mem_block *p,*p1;
	
#ifdef TRACE_PEEL
	if(CHK_PEEL(TRACE_LEVEL_1)) (void)printf("In %s()\n",__func__);
#endif
	work=&lk->peel->workspace;
	if(freq) {
		if(freq[0]) free(freq[0]);
		free(freq);
	}
	if(work) {
		if(work->r_funcs) free(work->r_funcs);
		if(work->s0) free(work->s0);
		if(work->s1) free(work->s1);
		if(work->s2) free(work->s2);
		if(work->s3) free(work->s3);
		if(work->s4) free(work->s4);
		if(work->s5) free(work->s5);
		if(work->s6) free(work->s6);
		if(work->s7) free(work->s7);
		if(work->s8) free(work->s8);
	}
	for(k=0;k<2;k++) {
		p1=lk->peel->first_mem_block[k];
		while(p1) {
			p=p1->next;
			if(p1->index) free(p1->index);
			if(p1->val) free(p1->val);
			free(p1);
			p1=p;
		}
	}
	seg_dealloc(lk);
	free_complex_mem();
}

void peel_alloc(struct loki *loki)
{
	int i,j,i1,k,k1,k2;
	struct peel_mem *work;
	
#ifdef TRACE_PEEL
	if(CHK_PEEL(TRACE_LEVEL_1)) (void)printf("In %s()\n",__func__);
#endif
	lk=loki;
	work=&loki->peel->workspace;
	k=loki->models->tlocus?2:0;
	for(i=0;i<loki->markers->n_markers;i++) {
		j=loki->markers->marker[i].locus.n_alleles;
		if(j>k) k=j;
	}
	if(!k) return;
	freq=lk_malloc(sizeof(void *)*2*loki->pedigree->n_genetic_groups);
	freq[0]=lk_malloc(sizeof(double)*k*2*loki->pedigree->n_genetic_groups);
	for(i=1;i<loki->pedigree->n_genetic_groups*2;i++) freq[i]=freq[i-1]+k;
	k1=k*k;
	work->s0=lk_malloc(sizeof(struct fset)*k1*k1);
	work->s1=lk_malloc(sizeof(int)*(k+k1*2));
	peel_alls=work->s1+2*k1;
	j=num_bits(k);
	j=1<<(j+j);
	/* Note the 16 below should change if more than diallelic trait loci are fitted */
	if(loki->peel->max_rfuncs) {
		work->r_funcs=lk_malloc(sizeof(struct R_Func)*loki->peel->max_rfuncs);
		work->s8=lk_malloc(sizeof(int)*loki->peel->max_out);
	} else {
		work->r_funcs=0;
		work->s8=0;
	}
	work->s2=lk_malloc(sizeof(double)*(16+k1+4*j));
	if(loki->peel->max_peel_off) {
		work->s3=lk_malloc(sizeof(lk_ulong)*loki->peel->max_peel_off*(k+2));
		work->s4=lk_malloc(sizeof(int)*loki->peel->max_peel_off);
	}
	k2=loki->peel->max_peel_ops;
	work->s5=lk_malloc(sizeof(void *)*k2);
	i1=loki->markers->n_markers+loki->params.max_tloci;
	work->s6=lk_malloc(sizeof(void *)*i1);
	work->s7=lk_malloc(sizeof(double)*2*i1);
	peel_list=(struct Peelseq_Head **)work->s5;
	loki->peel->first_mem_block[MRK_MBLOCK]=get_new_memblock(MB_SIZE,MRK_MBLOCK);
	loki->peel->first_mem_block[TRT_MBLOCK]=get_new_memblock(MB_SIZE,TRT_MBLOCK);
	seg_alloc(loki);
	if(atexit(peel_dealloc)) message(WARN_MSG,"Unable to register exit function peel_dealloc()\n");
}

static void rearrange_rfuncs(struct Complex_Element *em,struct R_Func *rf)
{
	int i,j,k,*inv,*ids,n_inv,n_ids,n_rf,*idx;
	struct R_Func *rf1;
	
	n_rf=em->n_rfuncs;
	if(n_rf) {
		idx=em->index;
		inv=em->involved;
		n_inv=em->n_involved;
		for(i=0;i<n_rf;i++) {
			rf1=rf+idx[i];
			ids=rf1->id_list;
			n_ids=rf1->n_ind;
			for(j=0;j<n_ids;j++) {
				for(k=0;k<n_inv;k++) {
					if(inv[k]==ids[j]) {
						ids[j]=k;
						break;
					}
				}
				assert(k<n_inv);
			}
		}
	}
}

/* Peels locus perm[idx] where the n_loci entries in perm are all in 1 linkage group.
* If sample_flag then the locus is to be sampled */
double peel_locus(struct Locus **llist,int idx,int n_loci,int sample_flag,struct loki *loki)
{
	struct Peelseq_Head *pp,*pp_head;
	struct Simple_Element *simple_em;
	struct Complex_Element *complex_em; 
	struct Id_Record *id,*id_array;
	struct Marker *mark=0;
	struct peel_mem *work;
	trait_pen_func *trait_pen=0;
	int comp,i,i1,j,k,k1,k2,k3,nn_all,n_peel_ops,ids,idd,allele,grp,cs,*a_trans;
	int rec,nrec,mtype,locus,sample_freq=0,trn[2],locus1,**seg,unlinked=0,*f_flag;
	int ***seglist,n_all,linktype,n_markers,si,*ind_list;
	double *recom1,*recom2;
	double like=0.0,like1,z,z1,theta,Mtp[2],Ptp[2],*pos,*eff,**count,*tpp,*freq1,*count1;
	signed char *freq_set;
	lk_ulong a,b,c,lump,**a_set=0;
	pen_func *pen=0;
	struct Locus *loc,*loc1;
	struct R_Func *rf;
	
#ifdef TRACE_PEEL
	if(CHK_PEEL(TRACE_LEVEL_1)) (void)printf("In %s(%p,%d,%d,%d)\n",__func__,(void *)llist,idx,n_loci,sample_flag);
#endif
	id_array=loki->pedigree->id_array;
	n_markers=loki->markers->n_markers;
	work=&loki->peel->workspace;
	si=loki->params.si_mode;
	seglist=(int ***)work->s6;
	recom1=work->s7;
	recom2=recom1+n_markers+loki->params.n_tloci;
	loc=llist[idx];
	sample_flag&=~OP_SAMPLING;
	if(loc->type&ST_TRAITLOCUS) { 
		/* Trait locus */
		locus1=n_markers;
		locus=-1-loc->index;
		nn_all=loc->n_alleles;
		if(loc->flag&TL_LINKED) {
			j=loc->link_group;
			linktype=loki->markers->linkage[j].type&LINK_TYPES_MASK;
		} else {
			linktype=LINK_AUTO;
			unlinked=1;
		}
		mtype=loki->models->models[0].var.type;
		if(!loki->data->n_nonid_records && !loki->models->use_student_t && (loc->flag&LOCUS_SAMPLED) && (!(mtype&ST_CENSORED) || loki->models->censor_mode)) trait_pen=s_penetrance1;
		else trait_pen=s_penetrance;
		sample_freq=sample_flag==1?1:0;
		pp_head=(linktype==LINK_AUTO?loki->peel->peelseq_head_gen:loki->peel->peelseq_head_gen_x);
	} else { 
		/* Marker locus */
		locus1=locus=loc->index;
		mark=loki->markers->marker+locus;
		nn_all=loc->n_alleles;
		a_set=mark->all_set;
		sample_freq=0;
		if(sample_flag==1) {
			for(k=0;k<loki->pedigree->n_genetic_groups;k++) {
				for(i=j=0;i<nn_all;i++) if(mark->freq_set[k][i]!=1) j++;
				sample_freq|=(j>1)?1:0;
			}
		}
		j=loc->link_group;
		linktype=loki->markers->linkage[j].type&LINK_TYPES_MASK;
		pp_head=loki->peel->peelseq_head[locus];
	}
	pos=loc->pos;
	seg=loc->seg;
#ifdef TRACE_PEEL
	if(CHK_PEEL(TRACE_LEVEL_2)) {
		(void)printf("locus %s, sample_flag=%d, nn_all=%d, linktype=%d\n",locus<0?"QTL":mark->name,sample_flag,nn_all,linktype);
	}
#endif
	if(nn_all<2) return 0.0;
	if(n_loci>1) {
		k2=0;
		for(i=0;i<idx;i++) {
			loc1=llist[i];
			seglist[k2]=loc1->seg;
			recom1[k2]=.5*(1.0-exp(-0.02*(pos[X_MAT]-loc1->pos[X_MAT])));
			recom2[k2++]=.5*(1.0-exp(-0.02*(pos[X_PAT]-loc1->pos[X_PAT])));
			if(recom1[k2-1]<=1.0e-8) recom1[k2-1]=1.0e-8;
			if(recom2[k2-1]<=1.0e-8) recom2[k2-1]=1.0e-8;
		}
		for(i=idx+1;i<n_loci;i++) {
			loc1=llist[i];
			seglist[k2]=loc1->seg;
			recom1[k2]=.5*(1.0-exp(0.02*(pos[X_MAT]-loc1->pos[X_MAT])));
			recom2[k2++]=.5*(1.0-exp(0.02*(pos[X_PAT]-loc1->pos[X_PAT])));
			if(recom1[k2-1]<=1.0e-8) recom1[k2-1]=1.0e-8;
			if(recom2[k2-1]<=1.0e-8) recom2[k2-1]=1.0e-8;
		}
	}
	count=freq+loki->pedigree->n_genetic_groups;
	if(sample_freq) {
		for(j=0;j<loki->pedigree->n_genetic_groups;j++) {
			z=1.0/(double)nn_all;
			if(locus>=0 && mark->count_flag[j]) for(i=0;i<nn_all;i++) count[j][i]=mark->counts[j][i]+z;
			else for(i=0;i<nn_all;i++) count[j][i]=z;
		}
	}
	if(locus<0)	{
		for(i=0;i<nn_all;i++) {
			peel_alls[i]=i;
			for(j=0;j<loki->pedigree->n_genetic_groups;j++) freq[j][i]=loc->freq[j][i];
		}
	}
	eff=0;
	mtype=0;
	if(locus<0 || (mark->mterm && mark->mterm[0])) {
		mtype=loki->models->models[0].var.type;
		eff=loc->eff[0];
	}
	for(j=comp=0;comp<loki->pedigree->n_comp;comp++) {
		for(k=0;k<2;k++) {
			loki->peel->mem_block[k]=loki->peel->first_mem_block[k];
			loki->peel->mem_block[k]->ptr=0;
		}
		rf=loki->peel->workspace.r_funcs;
		ind_list=loki->peel->workspace.s8;
		if(locus<0)	{
			n_all=nn_all;
			lump=0;
		} else {
			a_trans=mark->allele_trans[comp];
			n_all=mark->n_all1[comp];
			if(n_all<2)	{
				j+=loki->pedigree->comp_size[comp];
				continue;
			}
			for(i=0;i<nn_all;i++) {
				for(grp=0;grp<loki->pedigree->n_genetic_groups;grp++) freq[grp][i]=0.0;
				peel_alls[i]=n_all-1;
			}
			for(i=0;i<nn_all;i++) if((k=a_trans[i])>=0) peel_alls[k]=i;
			for(grp=0;grp<loki->pedigree->n_genetic_groups;grp++)
				for(i=0;i<nn_all;i++) freq[grp][peel_alls[i]]+=loc->freq[grp][i];
			lump=0;
			for(i=0;i<nn_all;i++) if(peel_alls[i]==n_all-1) lump|=(1<<i);
		}
		n_peel_ops=0;
		cs=loki->pedigree->comp_size[comp];
		id=id_array+j;
		i=cs;
		i1=j;
		f_flag=loc->founder_flag+j;
		while(i--) {
			id->flag=0;
			id->allele[X_MAT]=id->allele[X_PAT]= -1;
			k3=*(f_flag++);
			if(k3==2) {
				id++;
				i1++;
				continue;
			}
			/* Set up transmission probs. */
			if(!k3 && n_loci>1) {
				for(k=idx-1;k>=0;k--) {
					k2=seglist[k][X_MAT][i1];
					if(k2<0) continue;
					theta=recom1[k];
					if(theta<1.0e-16) {
						fprintf(stderr,"Warning: [1] theta very low (%g) for locus %d\n",theta,locus);
					}
					Mtp[k2]=1.0-theta;
					Mtp[1-k2]=theta;
					break;
				}
				if(k<0) Mtp[0]=Mtp[1]=1.0;
				for(k=idx-1;k>=0;k--) {
					k2=seglist[k][X_PAT][i1];
					if(k2<0) continue;
					theta=recom2[k];
					if(theta<1.0e-16) {
						fprintf(stderr,"Warning: [2] theta very low (%g) for locus %d\n",theta,locus);
					}
					Ptp[k2]=1.0-theta;
					Ptp[1-k2]=theta;
					break;
				}
				if(k<0) Ptp[0]=Ptp[1]=1.0;
				for(k=idx;k<n_loci-1;k++) {
					k2=seglist[k][X_MAT][i1];
					if(k2<0) continue;
					theta=recom1[k];
					if(theta<1.0e-16) theta=1.0e-16;
					Mtp[k2]*=1.0-theta;
					Mtp[1-k2]*=theta;
					break;
				}
				for(k=idx;k<n_loci-1;k++) {
					k2=seglist[k][X_PAT][i1];
					if(k2<0) continue;
					theta=recom2[k];
					if(theta<1.0e-16) theta=1.0e-16;
					Ptp[k2]*=1.0-theta;
					Ptp[1-k2]*=theta;
					break;
				}
				z=1.0/(Mtp[0]+Mtp[1]);
				z1=1.0/(Ptp[0]+Ptp[1]);
				for(k=0;k<2;k++) {
					Mtp[k]*=z;
					Ptp[k]*=z1;
					id->tpp[X_MAT][k]=Mtp[k];
					id->tpp[X_PAT][k]=Ptp[k];
				}
				id->tp[X_MM_PM]=Mtp[X_MAT]*Ptp[X_MAT];
				id->tp[X_MM_PP]=Mtp[X_MAT]*Ptp[X_PAT];
				id->tp[X_MP_PM]=Mtp[X_PAT]*Ptp[X_MAT];
				id->tp[X_MP_PP]=Mtp[X_PAT]*Ptp[X_PAT];
			} else {	
				for(k=0;k<2;k++) id->tpp[X_MAT][k]=id->tpp[X_PAT][k]=0.5;
				for(k=0;k<4;k++) id->tp[k]=0.25;
			}
			id->rfp= -1;
			if(locus>=0) {
				for(k=0;k<2;k++) {
					if(mark->nhaps[k][i1]==1) {
						a=mark->temp[k][i1];
						k1=1;
						while(!(a&1)) {
							a>>=1;
							k1++;
						}
						id->allele[k]=k1;
					}
				}
				if(mark->ngens[i1]==1) id->flag|=(SAMPLED_MAT|SAMPLED_PAT);
			}
			id++;
			i1++;
		}
		if(locus>=0 && mark->mterm && mark->mterm[0]) pen=&penetrance;
		sample_flag&=~OP_SAMPLING;
#ifdef TRACE_PEEL
		if(CHK_PEEL(TRACE_LEVEL_2)) (void)printf("Peeling component %d\n",comp+1);
#endif
		/* Go through all operations in peeling sequence */
		pp=pp_head+comp;
		like1=0.0;
		while(pp->type) {
			if(pp->type==PEEL_SIMPLE) {
				simple_em=pp->ptr.simple;
				if(sample_flag && simple_em->pivot) peel_list[n_peel_ops++]=pp;
				if(simple_em->pivot>0) {
					k=simple_em->out_index;
					k1=simple_em->pivot;
					rf[k].id_list=ind_list;
					rf[k].n_ind=2;
					*ind_list++=k1;
					*ind_list++=-k1;
				} else if(simple_em->pivot<0) {
					k=simple_em->out_index;
					rf[k].id_list=ind_list;
					rf[k].n_ind=2;
					k1=simple_em->sire;
					*ind_list++=k1;
					*ind_list++=-k1;
					rf[k+1].id_list=ind_list;
					rf[k+1].n_ind=2;
					k1=simple_em->dam;
					*ind_list++=k1;
					*ind_list++=-k1;
				}
				if(locus<0)	{
					z=loki_trait_simple_peelop(simple_em,locus,sample_flag,freq,rf,trait_pen,loki);
					if(z== -DBL_MAX) {
						if(sample_flag&1) ABT_FUNC("Internal error - shouldn't be here\n");
						/* Error return on likelihood (not sampling) run.  Clean up and return error */
						return z;
					}
					like1+=z;
				} else {
					switch(linktype) {
						case LINK_AUTO:
							z=loki_simple_peelop(simple_em,locus,sample_flag,*pen,a_set,freq,rf,loki);
							break;
						case LINK_X:
							z=loki_simple_peelop_x(simple_em,locus,sample_flag,*pen,a_set,freq,rf,loki);
							break;
						default:
							ABT_FUNC("Link type not implemented\n");
					}
					if(z== -DBL_MAX) {
						if(sample_flag&1) {
							ABT_FUNC("Internal error - shouldn't be here\n");
						}
						/* Error return on likelihood (not sampling) run.  Clean up and return error */
						return z;
					}
					like1+=z;
				}
				pp= &simple_em->next;
			} else {
				complex_em=pp->ptr.complex;
				k=complex_em->out_index;
				if(k>=0) {
					rf[k].id_list=ind_list;
					k2=complex_em->n_involved-complex_em->n_out;
					for(k1=0;k2<complex_em->n_involved;k2++) {
						*ind_list++=complex_em->involved[k2];
						k1++;
					}
					rf[k].n_ind=k1;
				}
				if(sample_flag && k>=0) peel_list[n_peel_ops++]=pp;
				rearrange_rfuncs(complex_em,rf);
				if(locus<0)	{
					z=loki_trait_complex_peelop(complex_em,locus,sample_flag,rf,trait_pen,freq,loki);
					if(z== -DBL_MAX) {
						if(sample_flag&1) ABT_FUNC("Internal error - shouldn't be here\n");
						/* Error return on likelihood (not sampling) run.  Clean up and return error */
						return z;
					}
					like1+=z;
				} else {	
					z=loki_complex_peelop(complex_em,locus,sample_flag,*pen,n_all,rf,freq,loki);
					if(z== -DBL_MAX) {
						if(sample_flag&1) ABT_FUNC("Internal error - shouldn't be here\n");
						/* Error return on likelihood (not sampling) run.  Clean up and return error */
						return z;
					}
					like1+=z;
				}
				pp= &complex_em->next;
			}
		}
		/* If we're sampling, do them again in reverse */
		if(sample_flag) {
			sample_flag|=OP_SAMPLING;
			for(i=n_peel_ops-1;i>=0;i--) {
				pp=peel_list[i];
				if(pp->type==PEEL_SIMPLE) {
					simple_em=pp->ptr.simple;
					if(simple_em->pivot== -2) {
						k1=simple_em->out_index;
					} else {
						if(locus<0) (void)loki_trait_simple_sample(simple_em,locus,sample_flag,freq,rf,trait_pen,loki);
						else {
							if(linktype==LINK_AUTO) (void)loki_simple_sample(simple_em,locus,*pen,a_set,freq,rf,loki);
							else (void)loki_simple_sample(simple_em,locus,*pen,a_set,freq,rf,loki);
						}
					}
				} else {
					complex_em=pp->ptr.complex;	
					if(locus<0) (void)loki_trait_complex_peelop(complex_em,locus,sample_flag,rf,trait_pen,freq,loki);
					else (void)loki_complex_peelop(complex_em,locus,sample_flag,*pen,n_all,rf,freq,loki);
				}
			}
			/* Get allele counts for frequency update */
			if(sample_freq) {
				f_flag=loc->founder_flag+j;
				if(locus>=0) {
					a_trans=mark->allele_trans[comp];
					for(id=id_array+j,i=0;i<cs;i++,id++) {
						if(*(f_flag++)!=1) continue;
						i1=i+j;
						grp=id->group-1;
						freq1=loc->freq[grp];
						count1=count[grp];
						for(k1=0;k1<2;k1++) {
							allele=id->allele[k1]-1;
							assert(allele>=0 && allele<n_all); /* Check if sampled allele is in range */
							a=mark->req_set[k1][i1];
							if(a&(1<<allele)) {
								z=0.0;
								k=0;
								while(a) {
									if(a&1) {
										k2=a_trans[k];
										if(k2<0) {
											b=lump;
											k3=0;
											while(b) {
												if(b&1) z+=freq1[k3];
												b>>=1;
												k3++;
											}
										} else z+=freq1[k2];
									}
									a>>=1;
									k++;
								}
								z=1.0/z;
								k=0;
								a=mark->req_set[k1][i1];
								while(a) {
									if(a&1) {
										k2=a_trans[k];
										if(k2<0) {
											b=lump;
											k3=0;
											while(b) {
												if(b&1) count1[k3]+=z*freq1[k3];
												b>>=1;
												k3++;
											}
										} else count1[k2]+=z*freq1[k2];
									}
									a>>=1;
									k++;
								}
							} else {
								k2=a_trans[allele];
								if(k2<0) {
									z=0.0;
									k3=0;
									b=lump;
									while(b) {
										if(b&1) z+=freq1[k3];
										b>>=1;
										k3++;
									}
									z=1.0/z;
									k3=0;
									b=lump;
									while(b) {
										if(b&1) count1[k3]+=z*freq1[k3];
										b>>=1;
										k3++;
									}
								} else count1[k2]+=1.0;
							}
						}
					}
				} else {
					for(id=id_array+j,i=0;i<cs;i++,id++) {
						if(*(f_flag++)!=1) continue;
						grp=id->group-1;
						for(k1=0;k1<2;k1++) {
							allele=id->allele[k1]-1;
							count[grp][allele]+=1.0;
						}
					}
				}
			}
			/* Correct residuals, get segregation pattern */
			if(locus<0) { 
				/* For trait loci */
				id=id_array+j;
				for(i=0;i<cs;i++,id++) {
					i1=i+j;
					if(loc->pruned_flag[i1]) continue;
					k=id->allele[X_MAT];
					k1=id->allele[X_PAT];
					if(k>k1) k=k*(k-1)/2+k1;
					else k=k1*(k1-1)/2+k;
					/* Correct residuals for new genotypes */
					if(id->res[0]) {
						nrec=id->n_rec;
						k1=(loc->flag&LOCUS_SAMPLED)?loc->gt[i1]:1;
						if(k1!=k) {
							z=(k1>1)?eff[k1-2]:0.0;
							if(k>1) z-=eff[k-2];
							for(rec=0;rec<nrec;rec++) {
								id->res[0][rec]+=z;
							}
						}
					} 
					loc->gt[i1]=k;
					idd=id->dam;
					if(idd && !loc->pruned_flag[idd-1])	{
						k=id->allele[X_MAT];
						for(k1=0;k1<2;k1++)
							trn[k1]=(k==id_array[idd-1].allele[k1])?1:0;
						assert(trn[0] || trn[1]);
						if(trn[0] && trn[1])	{
							if(si || unlinked)	{
								tpp=id->tpp[X_MAT];
								z=ranf()*(tpp[0]+tpp[1]);
								seg[X_MAT][i1]=(z<=tpp[0])?0:1;
							} else seg[X_MAT][i1]= -2;
						} else if(trn[0]) seg[X_MAT][i1]=0;
						else seg[X_MAT][i1]=1;
					} else seg[X_MAT][i1]= -1;
					ids=id->sire;
					if(ids && !loc->pruned_flag[ids-1])	{
						k=id->allele[X_PAT];
						for(k1=0;k1<2;k1++)
							trn[k1]=(k==id_array[ids-1].allele[k1])?1:0;
						assert(trn[0] || trn[1]);
						if(trn[0] && trn[1])	{
							if(si || unlinked)	{
								tpp=id->tpp[X_PAT];
								z=ranf()*(tpp[0]+tpp[1]);
								seg[X_PAT][i1]=(z<=tpp[0])?0:1;
							} else seg[X_PAT][i1]= -2;
						} else if(trn[0]) seg[X_PAT][i1]=0;
						else seg[X_PAT][i1]=1;
					} else seg[X_PAT][i1]= -1;
				}
			} else { 
				/* For marker loci */
				id=id_array+j;
				for(i=0;i<cs;i++,id++) {
					i1=i+j;
					if(loc->pruned_flag[i1]) continue;
					k=id->allele[X_MAT];
					k1=id->allele[X_PAT];
					if(k>k1) k=k*(k-1)/2+k1;
					else k=k1*(k1-1)/2+k;
					if(eff) {
						if(id->res[0]) {
							nrec=id->n_rec;
							k1=(loc->flag&LOCUS_SAMPLED)?loc->gt[i1]:1;
							if(k1!=k) {
								z=(k1>1)?eff[k1-2]:0.0;
								if(k>1) z-=eff[k-2];
								for(rec=0;rec<nrec;rec++) {
									id->res[0][rec]+=z;
								}
							}
						}
					}
					loc->gt[i1]=k;
					idd=id->dam;
					if(idd && !loc->pruned_flag[idd-1]) {
						c=mark->req_set[X_MAT][i1];
						k=id->allele[X_MAT]-1;
						assert(k>=0 && k<n_all); /* Check if sampled allele is in range */
						a=1<<k;
						if(a&c) a=c;
						for(k1=0;k1<2;k1++) {
							k=id_array[idd-1].allele[k1]-1;
							b=1<<k;
							trn[k1]=(a&b)?1:0;
						}
						assert(trn[0] || trn[1]);
						if(trn[0] && trn[1])	{
							if(si)	{
								tpp=id->tpp[X_MAT];
								z=ranf()*(tpp[0]+tpp[1]);
								seg[X_MAT][i1]=(z<=tpp[0])?0:1;
							} else seg[X_MAT][i1]= -2;
						} else {
							seg[X_MAT][i1]=trn[0]?0:1;
						}
					} else seg[X_MAT][i1]= -1;
					ids=id->sire;
					if(ids && !loc->pruned_flag[ids-1]) {
						c=mark->req_set[X_PAT][i1];
						k=id->allele[X_PAT]-1;
						assert(k>=0 && k<n_all); /* Check if sampled allele is in range */
						a=1<<k;
						if(a&c) a=c;
						for(k1=0;k1<2;k1++) {
							k=id_array[ids-1].allele[k1]-1;
							b=1<<k;
							trn[k1]=(a&b)?1:0;
						}
						assert(trn[0] || trn[1]);
						if(trn[0] && trn[1]) {
							if(si)	{
								tpp=id->tpp[X_PAT];
								z=ranf()*(tpp[0]+tpp[1]);
								seg[X_PAT][i1]=(z<=tpp[0])?0:1;
							} else seg[X_PAT][i1]= -2;
						} else {
							seg[X_PAT][i1]=trn[0]?0:1;
						}
					} else seg[X_PAT][i1]= -1;
				}
			}
		}
		like+=like1;
		j+=cs;
	}
	if(sample_freq) {
		for(grp=0;grp<loki->pedigree->n_genetic_groups;grp++) {
			z1=0.0;
			z=1.0;
			count1=count[grp];
			freq1=loc->freq[grp];	
			if(locus<0) {
				for(i=0;i<nn_all;i++) z1+=count1[i];
				for(i=0;i<nn_all-1;i++) {
					z1-=count1[i];
					freq1[i]=z*genbet(count1[i],z1);
					z-=freq1[i];
				}
				freq1[i]=z;
			} else {
				freq_set=mark->freq_set[grp];
				for(i=i1=0;i<nn_all;i++)	{
					if(freq_set[i]!=1) {
						z1+=count1[i];
						i1=i;
					} else z-=freq1[i];
				}
				for(i=0;i<i1;i++) if(freq_set[i]!=1) {
					z1-=count1[i];
					freq1[i]=z*genbet(count1[i],z1);
					z-=freq1[i];
				}
					freq1[i]=z;
			}
		}
	}
	if(sample_flag) loc->flag|=LOCUS_SAMPLED|RFMASK_OK;
	return like;
}
