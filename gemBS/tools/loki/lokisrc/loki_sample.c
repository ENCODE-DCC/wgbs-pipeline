/****************************************************************************
*                                                                          *
*     Loki - Programs for genetic analysis of complex traits using MCMC    *
*                                                                          *
*             Simon Heath - Rockefeller University                         *
*                                                                          *
*                       October 1997                                       *
*                                                                          *
* loki_sample.c:                                                           *
*                                                                          *
* Main sampling loop                                                       *
*                                                                          *
* Copyright (C) Simon C. Heath 1997, 2000, 2002                            *
* This is free software.  You can distribute it and/or modify it           *
* under the terms of the Modified BSD license, see the file COPYING        *
*                                                                          *
****************************************************************************/

#include <config.h>
#include <stdlib.h>
#include <string.h>
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif
#ifdef USE_DMALLOC
#include <dmalloc.h>
#endif
#include <math.h>
#include <stdio.h>
#include <float.h>
#include <errno.h>
#include <sys/time.h>
#include <sys/times.h>
#include <sys/types.h>
#include <sys/wait.h>

#include "ranlib.h"
#include "utils.h"
#include "libhdr.h"
#include "loki.h"
#include "loki_peel.h"
#include "loki_utils.h"
#include "seg_pen.h"
#include "gen_pen.h"
#include "loki_ibd.h"
#include "loki_dump.h"
#include "loki_tlmoves.h"
#include "sample_cens.h"
#include "handle_res.h"
#include "sample_effects.h"
#include "sample_rand.h"
#include "print_data.h"
#include "calc_var_locus.h"
#include "loki_output.h"
#include "loki_output_stat5.h"
#include "kinship.h"
#include "genedrop.h"
#include "sample_nu.h"
#include "meiosis_scan.h"
#include "update_segs.h"
#include "loki_npl.h"
#include "read_solar_idfile.h"
#include "pseudo_chrom.h"
#include "loki_haplo.h"
#include "loki_xmlout.h"

#define SEGS_COMPLETE 1
#define MK_GENES_OK 2
#define QTL_GENES_OK 4
#define GENES_OK (MK_GENES_OK|QTL_GENES_OK)

void SampleLoop(int read_dump_flag,int append_output_flag,struct loki *loki)
{
	int lp,i,j,i1,j2,k,k1,k2,flag=0,comp,*naffect=0,**affs=0,dumped=0,ibdflag,*ntl_linked=0,*ntl_linked1=0,lp1;
	int **haplo_store[2]={0,0},analysis,polygenic_flag,num_iter;
	FILE *fptr=0,*ffreq=0,*fpos=0,*fmpos=0,*hapfile=0,*qptr=0;
	double z,*ss=0,*ss2=0,**pairs=0,*trpos=0,*trpos1=0;
	struct output_gen *og,*og1;
	struct Locus **locilist=0,**mpos_perm=0,*loc;
	double addlog(double x1,double x2);
	char *mposname=0,*qname=0;
	char *gtt[]={"0 0","1 1","1 2","2 2"};
	int *solar_trans,n_ibd=0,n_markers,ps_sampled=0;
	long old_pos=-1;
	struct Marker *marker;
	FILE *xmloutfile=0; 
	int xmloutflag=1;
	int blocking=1;
	
	message(INFO_MSG,"Setting up sampler\n");
	/* Output copy of phenotype data if required */
	if(loki->names[LK_PHENFILE]) {
		Print_Data(loki->names[LK_PHENFILE],loki);
		free(loki->names[LK_PHENFILE]);
	}
	/* Output copy of genotype data if required */
	og=loki->params.Output_Gen;
	while(og) {
		Print_Genotypes(og,loki);
		og1=og->next;
		free(og->file);
		free(og);
		og=og1;
	}
	analysis=loki->params.analysis;
	n_markers=loki->markers->n_markers;
	marker=loki->markers->marker;
	/* Check for IBD analysis, and set up memory structures if required */
	ibdflag=SetupIBD(loki);
	if(ibdflag) analysis=(IBD_ANALYSIS|ESTIMATE_IBD);
	else if(analysis&IBD_ANALYSIS) SetupNPL(loki);
	message(INFO_MSG,"Analysis = %d\n",analysis,ibdflag);
	set_sort_sex(0);
	lp=lp1=0;
	/* If standard (quantative) analysis, define normal constants */
	if(!analysis) {
		z=RES_PRIOR_V0*.5;
		loki->sys.res_prior_konst=log(z)*z-lgamma(z)+log(RES_PRIOR_S0)*z;
	}
	if(analysis&NULL_ANALYSIS) loki->params.est_aff_freq=0;
	/* Restarting ? */
	if(read_dump_flag) {
		/* Yes - read in dump file */
		if(check_output(INFO_MSG)) {
			(void)fputs("Retrieving program state: ",stdout);
			(void)fflush(stdout);
		}
		j=read_dump(&lp,&lp1,&n_ibd,&old_pos,&flag,analysis,loki);
		if(j<0) {
			message(WARN_MSG,"read_dump() FAILED\n",stdout);
			ABT_FUNC(AbMsg);
		} else message(INFO_MSG,"OK\n");
		flag&=~GENES_OK;
	} else {
		/* No - Get initial genotype samples for all markers */
		if(!(analysis&NULL_ANALYSIS)) {
			for(k=0;k<n_markers;k++) {
				loc=&loki->markers->marker[k].locus;
				(void)peel_locus(&loc,0,1,1,loki);
				/*				printf("%d %g\n",k,z); */
			}
		}
		if(!analysis) {
			/* Quantitative analysis */
			if(!loki->models->tlocus) loki->params.start_tloci=0;
			/* If we have a starting number of QTL, create them now */
			else {
				for(i=0;i<loki->params.start_tloci;i++)	{
					k=get_new_traitlocus(2,loki);
					loc=loki->models->tlocus+k;
					loc->model_flag=1;
					if(k<0) ABT_FUNC("Internal error - Couldn't get trait locus\n");
					for(j=0;j<loki->pedigree->n_genetic_groups;j++) loc->freq[j][0]=loc->freq[j][1]=0.5;
					for(k1=0;k1<loki->models->n_models;k1++) {
						loc->eff[k1][0]=loc->eff[k1][1]=0.0;
					}
					loc->link_group=get_tl_position(loc->pos,loki);
					loc->flag=loc->link_group<0?TL_UNLINKED:TL_LINKED;
					loc->eff[0][1]=1.0;
					loki->models->grand_mean[0]=1.0;
					loki->models->residual_var[0]=0.0025;
					loc->freq[0][0]=0.9999;
					loc->freq[0][1]=0.0001;
					loc->flag=TL_UNLINKED;
					z=peel_locus(&loc,0,1,1,loki);
					printf("YOYOYO! z=%g\n",z);
					for(i=0;i<loki->pedigree->ped_size;i++) {
					  print_orig_triple(stdout,i+1);
					  printf("\t%s\n",gtt[loc->gt[i]]);
					}
					exit(0);
					/* Handle weights for t model */
					if(loki->models->use_student_t) {
						init_sample_nu(loki);
						sample_nu(loki);
					}
				}
				flag=0;
			}
		} 
		/* Set up M-Sampler stuff */
		if((loki->params.lm_ratio>0.0 || loki->params.est_aff_freq || (analysis&IBD_ANALYSIS)) && !(analysis&NULL_ANALYSIS)) {
			sample_segs(loki);
			for(k=0;k<n_markers;k++) {
				pass_founder_genes(&marker[k].locus,loki);
				for(comp=0;comp<loki->pedigree->n_comp;comp++) {
					marker[k].locus.lk_store[comp]=seg_pen(&marker[k].locus,comp,&i,0,loki);
					if(i) {
						(void)fprintf(stderr,"seg_pen returned error code %d for marker %s",(int)marker[k].locus.lk_store[comp],marker[k].name);
						(void)fprintf(stderr," comp %d\n",comp+1);
						ABT_FUNC("Illegal segregation pattern\n");
					}
				}
			}
			loc=loki->models->tlocus;
			for(k=0;k<loki->params.n_tloci;k++) {
				if(loc[k].flag) {
					pass_founder_genes(loc+k,loki);
					for(comp=0;comp<loki->pedigree->n_comp;comp++) {
						loc[k].lk_store[comp]=gen_pen(loc+k,comp,&i,0,loki);
					}
				}
			}
			flag=(loki->params.si_mode)?1:0;
			if(!(analysis&NULL_ANALYSIS)) {
				if(analysis&IBD_ANALYSIS) {
					if(!flag) {
						sample_segs(loki);
						flag=1;
					}
					for(k=0;k<n_markers;k++) pass_founder_genes(&marker[k].locus,loki);
				}
			}
		}
	}
	polygenic_flag=0;
	for(i=0;i<loki->models->n_models;i++) polygenic_flag|=loki->models->models[i].polygenic_flag;
	/* Setup outputfiles */
	if(ibdflag) {
		if(!loki->names[LK_IBDFILE]) loki->names[LK_IBDFILE]=make_file_name(".ibd");
	} else {
		if(!loki->names[LK_OUTPUTFILE]) loki->names[LK_OUTPUTFILE]=make_file_name(".out");
		if(!(analysis&(IBD_ANALYSIS))) {
			if(!loki->names[LK_POSFILE]) loki->names[LK_POSFILE]=make_file_name(".pos");
		}
	}
	if(loki->params.genv_out) qname=make_file_name(".qtl");
	for(i=0;i<n_markers;i++) if(marker[i].pos_set==2) break;
	if(i<n_markers) {
		mposname=make_file_name(".mpos");
		if(!(locilist=malloc(sizeof(void *)*n_markers))) ABT_FUNC(MMsg);
		if(!(mpos_perm=malloc(sizeof(void *)*n_markers))) ABT_FUNC(MMsg);
		for(k2=k1=0;k1<loki->markers->n_links;k1++) {
			get_locuslist(locilist,k1,&k,1);
			gnu_qsort(locilist,(size_t)k,sizeof(void *),cmp_loci);
			for(j=0;j<k;j++) mpos_perm[k2++]=locilist[j];
		}
		free(locilist);
		locilist=0;
	}
	if(loki->params.output_haplo && !loki->names[LK_HAPLOFILE]) loki->names[LK_HAPLOFILE]=make_file_name(".hap");
	if(append_output_flag) {
		if(qname) if(!(qptr=fopen(qname,"a"))) abt(__FILE__,__LINE__,"%s(): File Error.  Couldn't open '%s' for output\n",__func__,qname);
		if(mposname) if(!(fmpos=fopen(mposname,"a"))) abt(__FILE__,__LINE__,"%s(): File Error.  Couldn't open '%s' for output\n",__func__,mposname);
		if(loki->names[LK_FREQFILE]) if(!(ffreq=fopen(loki->names[LK_FREQFILE],"a"))) abt(__FILE__,__LINE__,"%s(): File Error.  Couldn't open '%s' for output\n",__func__,loki->names[LK_FREQFILE]);
		if(!ibdflag) {
			if(!(fptr=fopen(loki->names[LK_OUTPUTFILE],"a"))) abt(__FILE__,__LINE__,"%s(): File Error.  Couldn't open '%s' for output\n",__func__,loki->names[LK_OUTPUTFILE]);
			if(loki->names[LK_POSFILE]) if(!(fpos=fopen(loki->names[LK_POSFILE],"r+"))) abt(__FILE__,__LINE__,"%s(): File Error.  Couldn't open '%s' for output\n",__func__,loki->names[LK_POSFILE]);
		}
		if(loki->params.output_haplo) if(!(hapfile=fopen(loki->names[LK_HAPLOFILE],"a"))) abt(__FILE__,__LINE__,"%s(): File Error.  Couldn't open '%s' for output\n",__func__,loki->names[LK_HAPLOFILE]);
	} else {
		j=loki->sys.syst_var[SYST_BACKUPS].flag?loki->sys.syst_var[SYST_BACKUPS].data.value:1;
		if(j) {
			if(loki->names[LK_OUTPUTFILE]) i=mkbackup(loki->names[LK_OUTPUTFILE],j);
		}
		if(qname) {
			if(!(qptr=fopen(qname,"w"))) abt(__FILE__,__LINE__,"%s(): File Error.  Couldn't open '%s' for output\n",__func__,qname);
			OutputQTLHead(qptr,loki);
		}
		if(mposname) if(!(fmpos=fopen(mposname,"w"))) abt(__FILE__,__LINE__,"%s(): File Error.  Couldn't open '%s' for output\n",__func__,mposname);
		if(loki->names[LK_FREQFILE]) if(!(ffreq=fopen(loki->names[LK_FREQFILE],"w"))) abt(__FILE__,__LINE__,"%s(): File Error.  Couldn't open '%s' for output\n",__func__,loki->names[LK_FREQFILE]);
		if(!ibdflag) {
			if(!(fptr=fopen(loki->names[LK_OUTPUTFILE],"w"))) abt(__FILE__,__LINE__,"%s(): File Error.  Couldn't open '%s' for output\n",__func__,loki->names[LK_OUTPUTFILE]);
			if(loki->names[LK_POSFILE]) if(!(fpos=fopen(loki->names[LK_POSFILE],"w"))) abt(__FILE__,__LINE__,"%s(): File Error.  Couldn't open '%s' for output\n",__func__,loki->names[LK_POSFILE]);
			OutputHeader(fptr,loki);
		}
		if(loki->params.output_haplo) if(!(hapfile=fopen(loki->names[LK_HAPLOFILE],"w"))) abt(__FILE__,__LINE__,"%s(): File Error.  Couldn't open '%s' for output\n",__func__,loki->names[LK_HAPLOFILE]);
		if(!(xmloutfile=open_xmloutfile(loki))) exit(EXIT_FAILURE);
		output_xmlheader(xmloutfile,loki);
		if(ffreq) OutputFreqHeader(ffreq,loki);
	}
	if(qname) free(qname);
	if(mposname) free(mposname);
	OutputHeader(stdout,loki);
	if((analysis&IBD_ANALYSIS)||loki->params.output_haplo) {
		if(!(naffect=malloc(sizeof(int)*loki->pedigree->n_comp))) ABT_FUNC(MMsg);
		if(!(affs=malloc(sizeof(void *)*loki->pedigree->n_comp))) ABT_FUNC(MMsg);
		for(i=k1=comp=0;comp<loki->pedigree->n_comp;comp++) {
			naffect[comp]=0;
			for(j=0;j<loki->pedigree->comp_size[comp];j++) if(loki->pedigree->id_array[i+j].affected==2) naffect[comp]++;
			if(naffect[comp]) {
				if(!(affs[comp]=malloc(sizeof(int)*naffect[comp]))) ABT_FUNC(MMsg);
				for(k=j=0;j<loki->pedigree->comp_size[comp];j++) if(loki->pedigree->id_array[i+j].affected==2) affs[comp][k++]=i+j;
				k1+=k;
			} else affs[comp]=0;
			i+=loki->pedigree->comp_size[comp];
		}
		if(loki->params.output_haplo && k1) {
			if(!(haplo_store[0]=malloc(sizeof(void *)*2*k1))) ABT_FUNC(MMsg);
			haplo_store[1]=haplo_store[0]+k1;
			if(!(haplo_store[0][0]=malloc(sizeof(int)*2*k1*n_markers))) ABT_FUNC(MMsg);
			for(i=1;i<2*k1;i++) haplo_store[0][i]=haplo_store[0][i-1]+n_markers;
		}
	}
	if(analysis&IBD_ANALYSIS) {
		/* Allocate memory for affected only analyses */
		if(n_markers) {
			k=2;
			if(!(ss=malloc(sizeof(double)*n_markers*k))) ABT_FUNC(MMsg);
			if(!(locilist=malloc(sizeof(void *)*n_markers))) ABT_FUNC(MMsg);
			ss2=ss+n_markers;
		}
		for(i=0;i<n_markers;i++) ss[i]=ss2[i]=0.0;
		if(!(pairs=malloc(sizeof(void *)*loki->pedigree->n_comp))) ABT_FUNC(MMsg);
		for(i=comp=0;comp<loki->pedigree->n_comp;comp++) {
			if(naffect[comp]) {
				j=naffect[comp]*(naffect[comp]+1)/2;
				if(!(pairs[comp]=malloc(sizeof(double)*j))) ABT_FUNC(MMsg);
				for(k1=j=0;j<naffect[comp];j++) {
					for(k=0;k<=j;k++) {
						pairs[comp][k1++]=kinship(affs[comp][j]+1,affs[comp][k]+1,loki->pedigree->id_array);
					}
				}
			} else pairs[comp]=0;
			i+=loki->pedigree->comp_size[comp];
		}
	} else if(loki->models->tlocus && loki->params.max_tloci) {
		/* Allocate memory for quantitative analyses */
		if(!(ntl_linked=malloc(sizeof(int)*(loki->markers->n_links+1)*2))) ABT_FUNC(MMsg);
		ntl_linked1=ntl_linked+loki->markers->n_links+1;
		for(i=0;i<=loki->markers->n_links;i++) ntl_linked1[i]=-1;
		if(!(trpos=malloc(sizeof(double)*(1+loki->markers->sex_map)*2*loki->params.max_tloci))) ABT_FUNC(MMsg);
		trpos1=trpos+(1+loki->markers->sex_map)*loki->params.max_tloci;
		loc=loki->models->tlocus;
		if(read_dump_flag) {
			for(k=0;k<=loki->markers->n_links;k++) ntl_linked1[k]=0;
			for(k1=0;k1<loki->params.n_tloci;k1++) if(loc[k1].flag && !(loc[k1].flag&TL_LINKED)) ntl_linked1[0]++;
			for(k2=k=0;k<loki->markers->n_links;k++) {
				for(k1=0;k1<loki->params.n_tloci;k1++) {
					if((loc[k1].flag&TL_LINKED) && loc[k1].link_group==k) {
						ntl_linked1[k+1]++;
						trpos1[k2++]=loc[k1].pos[0];
						if(loki->markers->sex_map) trpos1[k2++]=loc[k1].pos[1];
					}
				}
			}
		}
	}
	loki->sys.catch_sigs=1;
	num_iter=loki->params.num_iter;
	for(++lp;!loki->sys.sig_caught && (!num_iter || lp<=num_iter);lp++) {
#ifdef USE_DMALLOC
		if(dmalloc_verify(0)==DMALLOC_ERROR) {
			(void)fprintf(stderr,"[%s:%d] %s(): Error returned from dmalloc_verify().\nAttempting to abort nicely.\n",__FILE__,__LINE__,__func__);
			break;
		}
#endif
		if(!(analysis&NULL_ANALYSIS)) {
			/* Sample genetic portion of model */
			if(loki->params.pseudo_flag) {
				if(lp>=loki->params.pseudo_start) {
					k=0;
					if(!ps_sampled) {
						sample_pseudo_chromosomes(loki);
						k=1;
					} else if(loki->params.pseudo_freq) {
						if(!(lp%loki->params.pseudo_freq)) {
							sample_pseudo_chromosomes(loki);
							k=1;
						}
					}
					if(k) {
						sample_segs(loki);
						update_seg_probs(0,3,loki);
						flag|=GENES_OK;
					}
					ps_sampled=1;
				}
			}
			if(ranf()<loki->params.lm_ratio) {
				/* M sample */
				if(!(flag&SEGS_COMPLETE)) {
					sample_segs(loki);
					flag|=SEGS_COMPLETE;
				}
				if(!(flag&MK_GENES_OK)) {
					update_seg_probs(0,(flag&QTL_GENES_OK)?1:3,loki);
					flag|=GENES_OK;
				}
				for(k=0;k<loki->markers->n_markers;k++) {
					for(comp=0;comp<loki->pedigree->n_comp;comp++) {
						z=loki->markers->marker[k].locus.lk_store[comp];
#ifndef NDEBUG
						if(isnan(z) || isinf(z)) {
							ABT_FUNC("Floating point errors\n");
						}
#endif
					}
				}
				for(k=0;k<loki->markers->n_links;k++) meiosis_scan(k,loki);
				update_seg_probs(5,0,loki);
				flag&=~MK_GENES_OK;
			} else { /* L sample */
				if(blocking) {
					if(!(flag&SEGS_COMPLETE)) {
						sample_segs(loki);
						flag|=SEGS_COMPLETE;
					}
					k=((flag&MK_GENES_OK)?0:1)|((flag&QTL_GENES_OK)?0:2);
					update_seg_probs(7,k,loki);
					
				} else {
					if((!loki->params.si_mode) && (flag&SEGS_COMPLETE)) {
						k=((flag&MK_GENES_OK)?0:1)|((flag&QTL_GENES_OK)?0:2);
						update_seg_probs(7,k,loki);
						reprune_segs(loki);
						flag&=~SEGS_COMPLETE;
					}
					for(k=0;k<loki->markers->n_links;k++) Sample_LinkageGroup(k,loki);
					k=GENES_OK;
					if(!loki->params.si_mode) k|=SEGS_COMPLETE;
					flag&=~k;
				}
			}
			if(loki->params.output_haplo) {
				if(loki->params.sample_freq[0] && lp>=loki->params.sample_from[0] && !(lp%loki->params.sample_freq[0])) {
					if(!(flag&SEGS_COMPLETE)) {
						sample_segs(loki);
						flag|=SEGS_COMPLETE;
					}
					/*					sample_haplo(hapfile,naffect,affs,haplo_store,loki); */
				}
			}
		}
		/* Estimate allele frequencies in affecteds only */
		if(loki->params.est_aff_freq) {
			/* Get full segregation pattern and founder gene info */
			if(!(flag&SEGS_COMPLETE)) {
				sample_segs(loki);
				flag|=SEGS_COMPLETE;
			}
			update_seg_probs(8,(flag&QTL_GENES_OK)?1:3,loki);
			flag|=GENES_OK;
		}
		if(analysis&ESTIMATE_IBD) {
			if(loki->params.sample_freq[1] && lp>=loki->params.sample_from[1] && !(lp%loki->params.sample_freq[1]))
				(void)printf("At: %d\n",lp);
			if(loki->params.sample_freq[0] && lp>=loki->params.sample_from[0] && !(lp%loki->params.sample_freq[0])) {
				n_ibd++;
				Handle_IBD(loki);
			}
		} else if(analysis&IBD_ANALYSIS) {
			if(analysis&NULL_ANALYSIS) genedrop(0,loki);
			for(k=0;k<n_markers;k++) {
				z=0.0;
				k2=0;
				for(comp=0;comp<loki->pedigree->n_comp;comp++) {
					if(naffect[comp]>1) {
						z+=score_ibd(k,&k1,comp,naffect[comp],affs[comp],pairs[comp],loki);
						k2+=k1;
					}
				}
				if(k2) z=sqrt(z/(double)k2);
				ss[k]=z;
				if(analysis&NULL_ANALYSIS) break;
			}
			if(loki->params.sample_freq[0] && lp>=loki->params.sample_from[0] && !(lp%loki->params.sample_freq[0])) {
				(void)fprintf(fptr,"%d",lp);
				for(k1=0;k1<loki->markers->n_links;k1++) {
					get_locuslist(locilist,k1,&k2,1);
					gnu_qsort(locilist,(size_t)k2,sizeof(void *),cmp_loci);
					for(k=0;k<k2;k++) (void)fprintf(fptr," %g",ss[locilist[k]->index]);
				}
				(void)fputc('\n',fptr);
			}
			if(loki->params.sample_freq[1] && lp>=loki->params.sample_from[1] && !(lp%loki->params.sample_freq[1])) {
				(void)fprintf(stdout,"%d",lp);
				for(k1=0;k1<loki->markers->n_links;k1++) {
					get_locuslist(locilist,k1,&k2,1);
					gnu_qsort(locilist,(size_t)k2,sizeof(int),cmp_loci);
					for(k=0;k<k2;k++) (void)fprintf(stdout," %g",ss[locilist[k]->index]);
				}
				(void)fputc('\n',stdout);
			}
		} else {
			if(!(loki->params.si_mode)) {
				reprune_segs(loki);
				flag&=~(GENES_OK|SEGS_COMPLETE);
			} else flag&=~GENES_OK;
			/* Quantitative analysis */
			if(loki->models->tlocus) {
				TL_Birth_Death(loki);
				loc=loki->models->tlocus;
				for(z=0.0,k=0;k<loki->params.n_tloci;k++) if(loc[k].flag) z+=1.0;
				for(k=0;k<loki->params.n_tloci;k++) {
					if(loc[k].flag) {
						if(ranf()<.5) Flip_TL_Alleles(k,loki);
						if(ranf()<.25/z) Flip_TL_Mode(k,loki);
						Sample_TL_Position(k,loki);
					}
				}
			}
			if(loki->models->use_student_t) sample_nu(loki);
			if(loki->models->censored_flag) Sample_Censored(loki);
			if(loki->models->models) {
				sample_effects(loki);
				if(loki->models->res_var_set[0]!=1) (void)Sample_ResVar(loki);
				if(loki->models->n_random) sample_rand(loki);
				if(polygenic_flag) sample_additive_var(loki);
				switch(loki->models->tau_mode) {
					case 0:
					case 1:
						loki->models->tau[0]=loki->models->tau_beta[0];
						break;
					case 2:
						loki->models->tau[0]=loki->models->residual_var[0]*loki->models->tau_beta[0];
						break;
				}
			}
#ifndef NDEBUG
			if((loki->params.debug_level)&4) (void)fputc('\n',stdout);
			z=Recalc_Res(0,loki);
			if(z>1.0e-8) printf("Warning: err=%g\n",z);
#endif
			k=0;
			loc=loki->models->tlocus;
			if(loki->params.sample_freq[0] && lp>=loki->params.sample_from[0] && !(lp%loki->params.sample_freq[0])) {
				for(i=0;i<loki->params.n_tloci;i++) if(loc[i].flag) calc_var_locus(loc+i,loki);
				OutputSample(fptr,lp,loki);
				if(qptr) OutputQTLvect(qptr,lp,loki);
				if(ffreq) OutputFreq(ffreq,lp,loki);
				k=1;
			}
			if(loki->params.sample_freq[1] && lp>=loki->params.sample_from[1] && !(lp%loki->params.sample_freq[1])) {
				if(!k) {
					for(i=0;i<loki->params.n_tloci;i++) if(loc[i].flag) calc_var_locus(loc+i,loki);
				}
				OutputSample(stdout,lp,loki);
			}
			if(ntl_linked && fpos && lp>=loki->params.sample_from[0]) {
				for(k=0;k<=loki->markers->n_links;k++) ntl_linked[k]=0;
				for(k1=0;k1<loki->params.n_tloci;k1++) if(loc[k1].flag && !(loc[k1].flag&TL_LINKED)) ntl_linked[0]++;
				for(k2=k=0;k<loki->markers->n_links;k++) {
					for(k1=0;k1<loki->params.n_tloci;k1++) {
						if((loc[k1].flag&TL_LINKED) && loc[k1].link_group==k) {
							ntl_linked[k+1]++;
							trpos[k2++]=loc[k1].pos[0];
							if(loki->markers->sex_map) trpos[k2++]=loc[k1].pos[1];
						}
					}
				}
				for(j=k=0;k<=loki->markers->n_links;k++) if(ntl_linked[k]!=ntl_linked1[k]) break;
				if(k==loki->markers->n_links+1) {
					for(k=0;k<k2;k++) if(trpos[k]!=trpos1[k]) break;
					if(k==k2) j=1;
				}
				if(old_pos>=0) (void)fseek(fpos,old_pos,SEEK_SET);
				else {
					(void)fseek(fpos,0,SEEK_END);
					old_pos=ftell(fpos);
				}
				if(lp1) (void)fprintf(fpos,":%d\n",lp-lp1);
				if(!j) {
					for(i1=k2=k=0;k<=loki->markers->n_links;k++) {
						ntl_linked1[k]=ntl_linked[k];
						if(ntl_linked[k]) {
							if(i1++) (void)fputc(' ',fpos);
							(void)fprintf(fpos,"%d %d",k,ntl_linked[k]);
							if(k) {
								for(k1=0;k1<ntl_linked[k];k1++) {
									for(j2=0;j2<1+loki->markers->sex_map;j2++) {
										trpos1[k2]=trpos[k2];
										(void)fprintf(fpos," %g",trpos[k2++]);
									}
								}
							}
						}
					}
					lp1=lp;
					old_pos=ftell(fpos);
					(void)fputs(":1\n",fpos);
				}
				(void)fflush(fpos);
			}
		}
		if(loki->params.sample_freq[0] && lp>=loki->params.sample_from[0] && !(lp%loki->params.sample_freq[0])) {
			if(!(flag&SEGS_COMPLETE)) {
				sample_segs(loki);
				flag|=SEGS_COMPLETE;
			}
			output_xmlout(xmloutfile,lp,xmloutflag++,loki);
		} 
		if(loki->params.dump_freq && lp && !(lp%loki->params.dump_freq)) {
			write_dump(lp,lp1,n_ibd,old_pos,flag,analysis,loki);
			dumped=1;
		} else dumped=0;
		if(fmpos && lp>=loki->params.sample_from[0]) {
			if(lp==1) {
				for(k=k1=0;k1<n_markers;k1++) {
					i=mpos_perm[k1]->index;
					if(marker[i].pos_set==2) {
						if(k++) (void)fputc(' ',fmpos);
						(void)fprintf(fmpos,"%s",marker[i].name);
					}
				}
				(void)fputc('\n',fmpos);
			}
			(void)fprintf(fmpos,"%d ",lp);
			for(k=k1=0;k1<n_markers;k1++) {
				i=mpos_perm[k1]->index;
				if(marker[i].pos_set==2) {
					(void)fprintf(fmpos,"%g ",mpos_perm[k1]->pos[0]);
					if(loki->markers->sex_map) (void)fprintf(fmpos,"%g ",mpos_perm[k1]->pos[1]);
				}
			}
			(void)fputc('\n',fmpos);
			(void)fflush(fmpos);
		}
	}
	if(analysis&ESTIMATE_IBD) {
		switch((loki->params.ibd_mode)&IBD_MODE_MASK) {
			case DEFAULT_IBD_MODE:
				Output_Sample_IBD(lp-1,n_ibd,loki);
				break;
			case MERLIN_IBD_MODE:
				Output_Merlin_IBD(n_ibd,loki);
				break;
			case SOLAR_IBD_MODE:
				if(!(solar_trans=malloc(sizeof(int)*loki->pedigree->ped_size))) ABT_FUNC(MMsg);
				read_solar_idfile(solar_trans,loki);
				Output_Solar_IBD(n_ibd,solar_trans,loki);
				free(solar_trans);
				break;
		}
	}
	if(!dumped) {
		write_xml_dump(lp-1,lp1,n_ibd,old_pos,flag,analysis,loki);
		/*		write_dump(lp-1,lp1,n_ibd,old_pos,flag,analysis,loki);*/
	}
	if(fptr) (void)fclose(fptr);
	if(fpos) (void)fclose(fpos);
	if(ffreq) (void)fclose(ffreq);
	if(hapfile) (void)fclose(hapfile);
	if(qptr) (void)fclose(qptr);
	if(xmloutfile) (void)close_xmloutfile(xmloutfile);
	if(polygenic_flag && loki->names[LK_POLYFILE]) {
		if(!(fptr=fopen(loki->names[LK_POLYFILE],"w"))) abt(__FILE__,__LINE__,"%s(): File Error.  Couldn't open '%s' for output\n",__func__,loki->names[LK_POLYFILE]);
		Output_BV(fptr,loki);
		(void)fclose(fptr);
	}
	if(mpos_perm) free(mpos_perm);
	if(ss) free(ss);
	if(locilist) free(locilist);
	if(naffect) free(naffect);
	if(affs) free(affs);
	if(pairs) free(pairs);
	if(ntl_linked) free(ntl_linked);
	if(trpos) free(trpos);
	if(haplo_store[0]) {
		if(haplo_store[0][0]) free(haplo_store[0][0]);
		free(haplo_store[0]);
	}
}
