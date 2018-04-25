/****************************************************************************
 *                                                                          *
 *     Loki - Programs for genetic analysis of complex traits using MCMC    *
 *                                                                          *
 *                      Simon Heath - MSKCC                                 *
 *                                                                          *
 *                          August 2000                                     *
 *                                                                          *
 * loki_dump.c:                                                             *
 *                                                                          *
 * Routines for reading and writing dump files                              *
 *                                                                          *
 * Copyright (C) Simon C. Heath 1997, 2000, 2002                            *
 * This is free software.  You can distribute it and/or modify it           *
 * under the terms of the Modified BSD license, see the file COPYING        *
 *                                                                          *
 ****************************************************************************/

#include <config.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/times.h>
#include <errno.h>
#include <sys/wait.h>
#include <signal.h>
#include <assert.h>

#include "utils.h"
#include "string_utils.h"
#include "loki.h"
#include "libhdr.h"
#include "loki_peel.h"
#include "loki_output.h"
#include "sample_rand.h"
#include "loki_ibd.h"
#include "loki_dump.h"

int read_dump(int *lp,int *lp1,int *n_ibd,long *old_pos,int *flag,int analysis,struct loki *loki)
{
	int i=0,j,k,k1,k2,nl=0,type,nrec,v2[4],pflag=0,polygenic_flag=0;
	char *tmp,*tmp1;
	double v1[4],z;
	struct Model *mod;
	struct Id_Record *id_array;
	struct Locus *loc=0;
	string *s=0;
	void *tbuf=0;
	FILE *fdump;

	signal(SIGCHLD,SIG_IGN);
	id_array=loki->pedigree->id_array;
	k=k1=k2=0;
	if(loki->names[LK_FILTER]) {
		errno=0;
		i=child_open(READ,loki->names[LK_DUMPFILE],loki->names[LK_FILTER]);
		fdump=fdopen(i,"r");
		if(errno && errno!=ESPIPE) i=__LINE__;
		else i=0;
		errno=0;
	} else fdump=fopen(loki->names[LK_DUMPFILE],"r");
	if(!fdump || i) (void)fprintf(stderr,"[%s:%d] %s(): File Error.  Couldn't open '%s' for input\n",__FILE__,__LINE__,__func__,loki->names[LK_DUMPFILE]);
	else {
		s=fget_string(fdump,s,&tbuf);
		if(!s->len) i=__LINE__;
		tmp=get_cstring(s);
		if(!i && strncmp(tmp,"Loki.dump:",10)) i=__LINE__;
		if(!i) {
			tmp1=tmp+10;
			for(j=0;!i && j<2;j++) {
				if(j && (*tmp1++!=',')) i=__LINE__;
				else v2[j]=strtol(tmp1,&tmp,16);
				tmp1=tmp;
			}
		}
		if(!i) {
			*lp=(int)v2[0];
			*lp1=(int)v2[1];
		}
		if(!i && (analysis&IBD_ANALYSIS)) {
			if(analysis&ESTIMATE_IBD) i=read_ibd_dump(fdump,n_ibd,tmp,loki);
			else if(*tmp!='\n') i=__LINE__;
		}
		(void)fputs("[Parameters] ",stdout);
		(void)fflush(stdout);
		if(!(analysis&IBD_ANALYSIS)) {
			for(polygenic_flag=k=0;k<loki->models->n_models;k++) if(loki->models->models[k].polygenic_flag) polygenic_flag++;
			for(j=0;!i && j<3;j++) {
				if(*tmp1++!=',') i=__LINE__;
				else v2[j]=strtol(tmp1,&tmp,16);
				tmp1=tmp;
			}
			if(!i) {
				nl=(int)v2[0];
				loki->sys.n_cov_columns=(int)v2[2];
				k=3;
				if(polygenic_flag) k++;
				for(j=0;!i && j<k;j++) {
					if(*tmp1++!=',') i=__LINE__;
					if(txt_get_double(tmp1,&tmp,v1+j)) i=__LINE__;
					tmp1=tmp;
				}
				for(j=0;!i && j<loki->models->n_random;j++) {
					if(*tmp1++!=',') i=__LINE__;
					if(txt_get_double(tmp1,&tmp,loki->models->c_var[j])) i=__LINE__;
					tmp1=tmp;
				}
			}
			if(!i) {
				loki->models->residual_var[0]=v1[0];
				loki->models->tau[0]=v1[1];
				loki->models->grand_mean[0]=v1[2];
				if(polygenic_flag) loki->models->additive_var[0]=v1[3];
				if(*tmp!='\n') i=__LINE__;
				else {
					s=fget_string(fdump,s,&tbuf);
					if(!s->len) i=__LINE__;
					tmp=get_cstring(s);
				}
			}
			mod=loki->models->models;
			for(j=0;!i && j<mod->n_terms;j++) {
				type=mod->term[j].vars[0].type;
				if(type&(ST_TRAITLOCUS|ST_ID)) continue;
				for(k=0;!i && k<mod->term[j].df;k++) {
					if(txt_get_double(tmp,&tmp1,mod->term[j].eff+k)) i=__LINE__;
					if(*tmp1=='\n') {
						s=fget_string(fdump,s,&tbuf);
						if(!s->len) i=__LINE__;
						tmp=get_cstring(s);
					} else {
						if(*tmp1++!=',') i=__LINE__;
						tmp=tmp1;
					}
				}
			}
			*old_pos=-1;
			*flag=0;
			if(!i) {
				if(strncmp(tmp,"LKST:",5)) {
					*old_pos=strtol(++tmp1,&tmp,10);
					if(*tmp1==',') {
						*flag=(int)strtol(++tmp1,&tmp,16);
						tmp1=tmp;
					}
					if(*tmp1!='\n') i=__LINE__;
					if(!i) {
						s=fget_string(fdump,s,&tbuf);
						if(!s->len) i=__LINE__;
						tmp=get_cstring(s);
					}
				}
			}
			if(!i && strncmp(tmp,"LKST:",5)) i=__LINE__;
			else if(!i) {
				j=(int)strtol(tmp+5,&tmp1,16);
				if(j!=N_MOVE_STATS) i=__LINE__;
			}
			for(j=0;!i && j<N_MOVE_STATS;j++) {
				if(*tmp1++!=',') i=__LINE__;
				else loki->sys.move_stats[j].success=(int)strtol(tmp1,&tmp,16);
				if(!i && *tmp++!=',') i=__LINE__;
				else loki->sys.move_stats[j].n=(int)strtol(tmp,&tmp1,16);
			}
			if(!i) {
				if(*tmp1!='\n') i=__LINE__;
				else {
					s=fget_string(fdump,s,&tbuf);
					if(!s->len) i=__LINE__;
					tmp=get_cstring(s);
				}
			}
			if(!i) {
				if(strncmp(tmp,"LKQT:",5)) i=__LINE__;
				else {
					j=(int)strtol(tmp+5,&tmp1,16);
					if(*tmp1!='\n' || j!=nl) i=__LINE__;
				}
			}
			for(j=0;!i && j<nl;j++) {
				s=fget_string(fdump,s,&tbuf);
				if(!s->len) i=__LINE__;
				tmp=get_cstring(s);
				if(!i) {
					k1=(int)strtol(tmp,&tmp1,16);
					assert(k1==2);
					k=get_new_traitlocus(k1,loki);
					loc=loki->models->tlocus+k;
					loc->model_flag=1;
					if(*tmp1++!=',') i=__LINE__;
					else {
						loc->link_group=(int)strtol(tmp1,&tmp,10);
						if(*tmp++!=',') i=__LINE__;
						else loc->flag=(int)strtol(tmp,&tmp1,16);
					}
					assert(k>=0);
				}
				for(k2=0;!i && k2<2;k2++) {
					if(*tmp1++!=',') i=__LINE__;
					else {
						if(txt_get_double(tmp1,&tmp,loc->pos+k2)) i=__LINE__;
						else tmp1=tmp;
					}
				}
				if(!i) {
					k1=k1*(k1+1)/2-1;
					for(k2=0;!i && k2<k1;k2++) {
						if(*tmp1++!=',') i=__LINE__;
						else {
							if(txt_get_double(tmp1,&tmp,loc->eff[0]+k2)) i=__LINE__;
							else tmp1=tmp;
						}
					}
				}
				for(k1=0;!i && k1<loki->pedigree->n_genetic_groups;k1++) {
					z=1.0;
					for(k2=0;!i && k2<loc->n_alleles-1;k2++) {
						if(*tmp1++!=',') i=__LINE__;
						else {
							if(txt_get_double(tmp1,&tmp,loc->freq[k1]+k2)) i=__LINE__;
							else {
								z-=loc->freq[k1][k2];
								tmp1=tmp;
							}
						}
					}
					if(!i) loc->freq[k1][k2]=z;
				}
				if(!i && *tmp1!='\n') i=__LINE__;
			}
		}
		if(!i) {
			s=fget_string(fdump,s,&tbuf);
			if(!s->len) i=__LINE__;
			tmp=get_cstring(s);
			if(!i) {
				if(strncmp(tmp,"LKMK:",5)) i=__LINE__;
				else {
					j=(int)strtol(tmp+5,&tmp1,16);
					if(*tmp1!='\n' || j!=loki->markers->n_markers) i=__LINE__;
				}
			}
		}
		for(j=0;!i && j<loki->markers->n_markers;j++) {
			loc=&loki->markers->marker[j].locus;
			s=fget_string(fdump,s,&tbuf);
			if(!s->len) i=__LINE__;
			tmp=get_cstring(s);
			k1=(int)strtol(tmp,&tmp1,16);
			if(*tmp1++!=',' || k1!=loc->n_alleles) i=__LINE__;
			else {
				loc->link_group=(int)strtol(tmp1,&tmp,10);
				if(*tmp++!=',') i=__LINE__;
				else loc->flag=((int)strtol(tmp,&tmp1,16))&~RFMASK_OK;
			}
			for(k=0;!i && k<2;k++) {
				if(*tmp1++!=',') i=__LINE__;
				else {
					if(txt_get_double(tmp1,&tmp,loc->pos+k)) i=__LINE__;
					else tmp1=tmp;
				}
			}
			for(k=0;!i && k<loki->pedigree->n_genetic_groups;k++) {
				z=1.0;
				for(k2=0;!i && k2<k1-1;k2++) {
					if(*tmp1++!=',') i=__LINE__;
					else {
						if(txt_get_double(tmp1,&tmp,loc->freq[k]+k2)) i=__LINE__;
						else {
							z-=loc->freq[k][k2];
							tmp1=tmp;
						}
					}
				}
				if(!i) loc->freq[k][k2]=z;
			}
			if(!i && *tmp1!='\n') i=__LINE__;
		}
		if(!i && loki->markers->n_markers+nl) {
			(void)fputs("[Genotypes] ",stdout);
			(void)fflush(stdout);
			if(nl) {
				s=fget_string(fdump,s,&tbuf);
				tmp=get_cstring(s);
				if(!s->len) i=__LINE__;
				else if(strncmp(tmp,"LKQS:",5)) i=__LINE__;
				else k2=(int)strtol(tmp+5,&tmp1,16);
				if(!i && *tmp1++!=',') i=__LINE__;
				else {
					k=(int)strtol(tmp1,&tmp,16);
					if(*tmp!='\n' || k!=loki->pedigree->ped_size || k2>=k) i=__LINE__;
				}
				if(k2) for(j=0;!i && j<loki->params.n_tloci;j++) if(loki->models->tlocus[j].flag) {
 					s=fget_string(fdump,s,&tbuf);
					tmp=get_cstring(s);
					if(!s->len) i=__LINE__;
					loc=loki->models->tlocus+j;
					for(k=0;!i && k<loki->pedigree->ped_size;k++) if(id_array[k].sire || id_array[k].dam) {
						k1=(int)*tmp++;
						if(k1>='0'&&k1<='9') k1-='0';
						else if(k1>='a'&&k1<='f') k1-='a'-10;
						else i=__LINE__;
						if(!i) {
							loc->seg[1][k]=(k1&3)-2;
							loc->seg[0][k]=(k1>>2)-2;
						}
					}
					if(!i && *tmp!='\n') i=__LINE__;
				}
				s=fget_string(fdump,s,&tbuf);
				tmp=get_cstring(s);
				if(!s->len) i=__LINE__;
				else if(strncmp(tmp,"LKQG\n",5)) i=__LINE__;
				for(k=0;!i && k<loki->pedigree->ped_size;k++) {
					s=fget_string(fdump,s,&tbuf);
					tmp=get_cstring(s);
					if(!s->len) i=__LINE__;
					for(k1=j=0;!i && j<loki->params.n_tloci;j++) if(loki->models->tlocus[j].flag) {
						if(k1++ && *tmp++!=',') i=__LINE__;
						else {
							loki->models->tlocus[j].gt[k]=(int)strtol(tmp,&tmp1,16);
							tmp=tmp1;
						}
					}
					if(!i && *tmp1!='\n') i=__LINE__;
				}
			}
			if(!i && loki->markers->n_markers) {
				s=fget_string(fdump,s,&tbuf);
				tmp=get_cstring(s);
				if(!s->len) i=__LINE__;
				else if(strncmp(tmp,"LKMS:",5)) i=__LINE__;
				else k2=(int)strtol(tmp+5,&tmp1,16);
				if(!i && *tmp1++!=',') i=__LINE__;
				else {
					k=(int)strtol(tmp1,&tmp,16);
					if(*tmp!='\n' || k!=loki->markers->n_markers) i=__LINE__;
				}
				if(k2) for(j=0;!i && j<loki->markers->n_markers;j++) {
					loc=&loki->markers->marker[j].locus;
					s=fget_string(fdump,s,&tbuf);
					tmp=get_cstring(s);
					if(!s->len) i=__LINE__;
					for(k=0;!i && k<loki->pedigree->ped_size;k++) if(id_array[k].sire || id_array[k].dam) {
						k1=(int)*tmp++;
						if(k1>='0'&&k1<='9') k1-='0';
						else if(k1>='a'&&k1<='f') k1-='a'-10;
						else i=__LINE__;
						if(!i) {
							loc->seg[1][k]=(k1&3)-2;
							loc->seg[0][k]=(k1>>2)-2;
						}
					}
					if(!i && *tmp!='\n') i=__LINE__;
				}
				for(k=j=0;!i && j<loki->markers->n_markers;j++) if(loki->markers->marker[j].mterm[0]) k++;
				if(k) {
					s=fget_string(fdump,s,&tbuf);
					tmp=get_cstring(s);
					if(!s->len) i=__LINE__;
					else if(strncmp(tmp,"LKMG:",5)) i=__LINE__;
					else k2=(int)strtol(tmp+5,&tmp1,16);
					if(*tmp1++!=',' || k2!=k) i=__LINE__;
					else {
						k2=(int)strtol(tmp1,&tmp,16);	
						if(k2!=loki->pedigree->ped_size) i=__LINE__;
					}
					for(k1=0;!i && k1<k;k1++) {
						if(*tmp++!=',') i=__LINE__;
						else k2=(int)strtol(tmp,&tmp1,16);
						if(k2<0 || k2>=loki->markers->n_markers || !loki->markers->marker[k2].mterm[0]) i=__LINE__;
						else tmp=tmp1;
					}
					if(!i && *tmp++!='\n') i=__LINE__;
					for(k=0;!i && k<loki->pedigree->ped_size;k++) {
						s=fget_string(fdump,s,&tbuf);
						tmp=get_cstring(s);
						if(!s->len) i=__LINE__;
						for(k1=j=0;!i && j<loki->markers->n_markers;j++) if(loki->markers->marker[j].mterm[0]) {
							if(k1++ && *tmp++!=',') i=__LINE__;
							else {
								loki->markers->marker[j].locus.gt[k]=(int)strtol(tmp,&tmp1,16);
								tmp=tmp1;
							}
						}
						if(!i && *tmp!='\n') i=__LINE__;
					}
				}
			}
		}
		if(!i && !(analysis&IBD_ANALYSIS)) {
			(void)fputs("[Residuals] ",stdout);
			(void)fflush(stdout);
			s=fget_string(fdump,s,&tbuf);
			tmp=get_cstring(s);
			if(!s->len) i=__LINE__;
			else if(strncmp(tmp,"LKRS:",5)) i=__LINE__;
			else k2=(int)strtol(tmp+5,&tmp1,16);
			if(*tmp1++!='\n' || k2!=loki->models->censored_flag) i=__LINE__;
			for(j=0;!i && j<loki->pedigree->ped_size;j++) {
				if(!i && polygenic_flag) {
					s=fget_string(fdump,s,&tbuf);
					tmp=get_cstring(s);
					if(!s->len) i=__LINE__;
					if(txt_get_double(tmp,&tmp1,&id_array[j].bv[0])) i=__LINE__;
					if(!i) {
						if(*tmp1=='\n') {
							if(pflag==1) i=__LINE__;
							else pflag=2;
						} else {
							if(!i && *tmp1++!=',') i=__LINE__;
							if(pflag==2) i=__LINE__;
							pflag=1;
							if(txt_get_double(tmp1,&tmp,&id_array[j].bvsum[0])) i=__LINE__;
							if(!i && *tmp++!=',') i=__LINE__;
							if(txt_get_double(tmp,&tmp1,&id_array[j].bvsum2[0])) i=__LINE__;
							if(!i && *tmp1!='\n') i=__LINE__;
						}
					}
				}
				if(!id_array[j].res[0]) continue;
				s=fget_string(fdump,s,&tbuf);
				tmp=get_cstring(s);
				if(!s->len) i=__LINE__;
				nrec=id_array[j].n_rec;
				if(!nrec) nrec=1;
				for(k=0;!i && k<nrec;k++) {
					if(k && *tmp++!=',') i=__LINE__;
					else if(txt_get_double(tmp,&tmp1,id_array[j].res[0]+k)) i=__LINE__;
					if(!i && loki->models->use_student_t) {
						if(*tmp1++!=',') i=__LINE__;
						else if(txt_get_double(tmp1,&tmp,id_array[j].vv[0]+k)) i=__LINE__;
						tmp1=tmp;
					}
					if(!i && loki->models->censored_flag && id_array[j].pseudo_qt[0]) {
						if(*tmp1++!=',') i=__LINE__;
						else if(txt_get_double(tmp1,&tmp,id_array[j].pseudo_qt[0]+k)) i=__LINE__;
					} else tmp=tmp1;
				}
				if(!i && *tmp!='\n') i=__LINE__;
			}
			if(!i && pflag==1) {
				s=fget_string(fdump,s,&tbuf);
				tmp=get_cstring(s);
				if(!s->len) i=__LINE__;
				else loki->params.bv_iter=(int)strtol(tmp,&tmp1,16);
				if(!i && *tmp1!='\n') i=__LINE__;
			}
		}
		if(!i) {
			(void)fputs("[Seed] ",stdout);
			(void)fflush(stdout);
			s=fget_string(fdump,s,&tbuf);
			tmp=get_cstring(s);
			if(!s->len) i=__LINE__;
			else if(strncmp(tmp,"LKSD:",5)) i=__LINE__;
			if(!i && binreadseed(fdump,tmp+5)) i=__LINE__;
		}
		if(!i) {
			(void)fputs("[Time",stdout);
			(void)fflush(stdout);
			s=fget_string(fdump,s,&tbuf);
			tmp=get_cstring(s);
			if(!s->len) i=__LINE__;
			else {
				if(strncmp(tmp,"LKTM:",5)) {
					(void)fputs(" !Not found!] ",stdout);
					(void)fflush(stdout);
				} else {
					(void)fputs("] ",stdout);
					(void)fflush(stdout);
					if(txt_get_double(tmp+5,&tmp1,&loki->sys.lktime.extra_time)) i=__LINE__;
					if(!i && *tmp1++!=',') i=__LINE__;
					if(!i && txt_get_double(tmp1,&tmp,&loki->sys.lktime.extra_stime)) i=__LINE__;
					if(!i && *tmp++!=',') i=__LINE__;
					if(!i && txt_get_double(tmp,&tmp1,&loki->sys.lktime.extra_utime)) i=__LINE__;
					if(!i && *tmp1!='\n') i=__LINE__;
					if(!i) {
						s=fget_string(fdump,s,&tbuf);
						tmp=get_cstring(s);
						if(!s->len) i=__LINE__;
					}
				}
				if(!i && strncmp(tmp,"Ldmp.end\n",9)) i=__LINE__;
			}
		}
	}
	if(s) free_string(s);
	if(tbuf) free_fget_buffer(&tbuf);
	if(fdump) (void)fclose(fdump);
	if(i) {
		if(i>0) (void)fprintf(stderr,"[%s:%d] Error in %s()\n",__FILE__,i,__func__);
		else {
			(void)fprintf(stderr,"[%s:%d] Error in %s()\n",__FILE__,-i,"read_ibd_dump");
			i=-i;
		}
	}
	if(!i && pflag==2) (void)fputs("Polygenic summary information not in file\n",stdout);
	signal(SIGCHLD,SIG_DFL);
	while(waitpid(-1,&i,WNOHANG)>0);
	return -i;
}

void write_dump(int lp,int lp1,int n_ibd,long old_pos,int flag,int analysis,struct loki *loki)
{
	int i=0,j,k,k1,k2,nx,nl=0,type,nrec,polygenic_flag=0;
	double ex1,ex2;
	struct tms tbuf;
	long tps;
	FILE *fdump;
	struct Model *mod=0;
	struct Id_Record *id_array;
	struct Locus *loc;
	
	signal(SIGCHLD,SIG_IGN);
	id_array=loki->pedigree->id_array;
	(void)fputs("Dumping program state: ",stdout);
	(void)fflush(stdout);
	if(!loki->names[LK_DUMPFILE]) loki->names[LK_DUMPFILE]=make_file_name(".dump");
	if(loki->names[LK_DUMPFILE]) {
		j=loki->sys.syst_var[SYST_BACKUPS].flag?loki->sys.syst_var[SYST_BACKUPS].data.value:1;
		if(j) i=mkbackup(loki->names[LK_DUMPFILE],j);
		if(!i) {
			if(loki->names[LK_FILTER]) {
				errno=0;
				i=child_open(WRITE,loki->names[LK_DUMPFILE],loki->names[LK_FILTER]);
				fdump=fdopen(i,"w");
				if(errno && errno!=ESPIPE) i=1;
				else i=0;
				errno=0;
			} else fdump=fopen(loki->names[LK_DUMPFILE],"w");
			if(!fdump || i) (void)fprintf(stderr,"[%s:%d] %s(): File Error.  Couldn't open '%s' for output\n",__FILE__,__LINE__,__func__,loki->names[LK_DUMPFILE]);
			else {
				if(!i && fprintf(fdump,"Loki.dump:%x,%x",lp,lp1)<0) i=1;
				if(!i && (analysis&IBD_ANALYSIS)) {
					if(analysis&ESTIMATE_IBD) i=write_ibd_dump(fdump,n_ibd,loki);
					else if(fputc('\n',fdump)<0) i=1;
				}
				(void)fputs("[Parameters] ",stdout);
				(void)fflush(stdout);
				if(!(analysis&IBD_ANALYSIS)) {
					for(polygenic_flag=k=0;k<loki->models->n_models;k++) if(loki->models->models[k].polygenic_flag) polygenic_flag++;
					nx=0;
					if(loki->models->models) {
						mod=loki->models->models;
						for(j=0;j<mod->n_terms;j++) {
							type=mod->term[j].vars[0].type;
							if(type&(ST_TRAITLOCUS|ST_ID)) continue;
							nx+=mod->term[j].df;
						}
					}
					for(nl=j=0;j<loki->params.n_tloci;j++) if(loki->models->tlocus[j].flag) nl++;
					if(!i && fprintf(fdump,",%x,%x,%x,",nl,nx,loki->sys.n_cov_columns)<0) i=1;
					if(loki->models->models) {
						if(!i) i=txt_print_double(loki->models->residual_var[0],fdump);
						if(!i && fputc(',',fdump)<0) i=1;
						if(!i) i=txt_print_double(loki->models->tau[0],fdump);
						if(!i && fputc(',',fdump)<0) i=1;
						if(!i) i=txt_print_double(loki->models->grand_mean[0],fdump);
						if(!i && polygenic_flag) {
							if(!i && fputc(',',fdump)<0) i=1;
							i=txt_print_double(loki->models->additive_var[0],fdump);
						}
						for(j=0;!i && j<loki->models->n_random;j++) {
							if(!i && fputc(',',fdump)<0) i=1;
							i=txt_print_double(loki->models->c_var[j][0],fdump);
						}
						if(!i && fputc('\n',fdump)<0) i=1;
						for(j=0;!i && j<mod->n_terms;j++) {
							type=mod->term[j].vars[0].type;
							if(type&(ST_TRAITLOCUS|ST_ID)) continue;
							for(k=0;!i && k<mod->term[j].df;k++) {
								i=txt_print_double(mod->term[j].eff[k],fdump);	
								if(!i && fputc('\n',fdump)<0) i=1;
							}
						}
					}
					if(!i && fprintf(fdump,"%ld,%x",old_pos,flag)<0) i=1;
					if(!i && fprintf(fdump,"\nLKST:%x",N_MOVE_STATS)<0) i=1;
					for(j=0;!i && j<N_MOVE_STATS;j++) if(fprintf(fdump,",%x,%x",loki->sys.move_stats[j].success,loki->sys.move_stats[j].n)<0) i=1;
					if(!i && fprintf(fdump,"\nLKQT:%x\n",nl)<0) i=1;
					for(j=0;!i && j<loki->params.n_tloci;j++) {
						loc=loki->models->tlocus+j;
						if(loc->flag) {
							if(fprintf(fdump,"%x,%d,%x",loc->n_alleles,loc->link_group,loc->flag)<0) i=1;
							for(k=0;!i && k<2;k++) {
								if(fputc(',',fdump)<0) i=1;
								else i=txt_print_double(loc->pos[k],fdump);
							}
							k=loc->n_alleles;
							k1=k*(k+1)/2-1;
							for(k=0;!i && k<k1;k++) {
								if(fputc(',',fdump)<0) i=1;
								else i=txt_print_double(loki->models->tlocus[j].eff[0][k],fdump);
							}
							for(k=0;!i && k<loki->pedigree->n_genetic_groups;k++) {
								for(k1=0;!i && k1<loc->n_alleles-1;k1++) {
									if(fputc(',',fdump)<0) i=1;
									else i=txt_print_double(loc->freq[k][k1],fdump);
								}
							}
							if(!i && fputc('\n',fdump)<0) i=1;
						}
					}
				}
				if(!i && fprintf(fdump,"LKMK:%x\n",loki->markers->n_markers)<0) i=1;
				for(j=0;!i && j<loki->markers->n_markers;j++) {
					loc=&loki->markers->marker[j].locus;
					if(fprintf(fdump,"%x,%d,%x",loc->n_alleles,loc->link_group,loc->flag)<0) i=1;
					for(k=0;!i && k<2;k++) {
						if(fputc(',',fdump)<0) i=1;
						else i=txt_print_double(loc->pos[k],fdump);	
					}
					for(k=0;!i && k<loki->pedigree->n_genetic_groups;k++) {
						for(k1=0;!i && k1<loc->n_alleles-1;k1++) {
							if(fputc(',',fdump)<0) i=1;
							else i=txt_print_double(loc->freq[k][k1],fdump);
						}
					}
					if(!i && fputc('\n',fdump)<0) i=1;
				}
				if(!i && loki->markers->n_markers+nl) {
					(void)fputs("[Genotypes] ",stdout);
					(void)fflush(stdout);
					for(k2=k=0;k<loki->pedigree->ped_size;k++) if(id_array[k].sire || id_array[k].dam) k2++;
					for(k1=j=0;j<loki->params.n_tloci;j++) {
						loc=loki->models->tlocus+j;
						if(loc->flag) k1++;
					}
					if(nl && k1) {
						if(k2) {
							if(fprintf(fdump,"LKQS:%x,%x\n",k2,loki->pedigree->ped_size)<0) i=1;
							for(j=0;!i && j<loki->params.n_tloci;j++) {
								loc=loki->models->tlocus+j;
								if(loc->flag) {
									for(k=0;!i && k<loki->pedigree->ped_size;k++) if(id_array[k].sire || id_array[k].dam) {
										if(fprintf(fdump,"%x",((loc->seg[0][k]+2)<<2)|(loc->seg[1][k]+2))<0) i=1;
									}
									if(!i && fputc('\n',fdump)<0) i=1;
								}
							}
						}
						if(!i && fputs("LKQG\n",fdump)<0) i=1;
						for(k=0;!i && k<loki->pedigree->ped_size;k++) {
							for(k1=j=0;!i && j<loki->params.n_tloci;j++) {
								loc=loki->models->tlocus+j;
								if(loc->flag) {
									if(k1++ && fputc(',',fdump)<0) i=1;
									if(!i && fprintf(fdump,"%x",loc->gt[k])<0) i=1;
								}
							}
							if(!i && fputc('\n',fdump)<0) i=1;
						}
					}
					if(loki->markers->n_markers) {
						if(k2) {
							if(fprintf(fdump,"LKMS:%x,%x\n",k2,loki->markers->n_markers)<0) i=1;
							for(j=0;!i && j<loki->markers->n_markers;j++) {
								loc=&loki->markers->marker[j].locus;
								for(k=0;!i && k<loki->pedigree->ped_size;k++) if(id_array[k].sire || id_array[k].dam) {
									if(fprintf(fdump,"%x",((loc->seg[0][k]+2)<<2)|(loc->seg[1][k]+2))<0) i=1;
								}
								if(!i && fputc('\n',fdump)<0) i=1;
							}
						}
						if(loki->models->models) {
							for(k=j=0;!i && j<loki->markers->n_markers;j++) if(loki->markers->marker[j].mterm[0]) k++;
							if(!i && k) {
								if(fprintf(fdump,"LKMG:%x,%x",k,loki->pedigree->ped_size)<0) i=1;
								for(j=0;!i && j<loki->markers->n_markers;j++) if(loki->markers->marker[j].mterm[0] && fprintf(fdump,",%x",j)<0) i=1;
								if(!i && fputc('\n',fdump)<0) i=1;
								for(k=0;!i && k<loki->pedigree->ped_size;k++) {
									for(k1=j=0;!i && j<loki->markers->n_markers;j++) if(loki->markers->marker[j].mterm[0]) {
										if(k1++ && fputc(',',fdump)<0) i=1;
										if(!i && fprintf(fdump,"%x",loki->markers->marker[j].locus.gt[k])<0) i=1;
									}
									if(!i && fputc('\n',fdump)<0) i=1;
								}
							}
						}
					}
				}
				if(!i && !(analysis&IBD_ANALYSIS)) {
					(void)fputs("[Residuals] ",stdout);
					(void)fflush(stdout);
					if(fprintf(fdump,"LKRS:%x\n",loki->models->censored_flag)<0) i=1;
					if(loki->models->models) for(j=0;!i && j<loki->pedigree->ped_size;j++) {
						if(polygenic_flag) {
							i=txt_print_double(id_array[j].bv[0],fdump);
							if(!i && fputc(',',fdump)<0) i=1;
							if(!i) i=txt_print_double(id_array[j].bvsum[0],fdump);
							if(!i && fputc(',',fdump)<0) i=1;
							if(!i) i=txt_print_double(id_array[j].bvsum2[0],fdump);
							if(!i && fputc('\n',fdump)<0) i=1;
						}
						if(!id_array[j].res[0]) continue;
						nrec=id_array[j].n_rec;
						if(!nrec) nrec=1;
						for(k=0;!i && k<nrec;k++) {
							if(k && fputc(',',fdump)<0) i=1;
							else i=txt_print_double(id_array[j].res[0][k],fdump);
							if(!i && loki->models->use_student_t) {
								if(fputc(',',fdump)<0) i=1;
								else i=txt_print_double(id_array[j].vv[0][k],fdump);
							}
							if(!i && loki->models->censored_flag && id_array[j].pseudo_qt[0]) {
								if(fputc(',',fdump)<0) i=1;
								else i=txt_print_double(id_array[j].pseudo_qt[0][k],fdump);
							}
						}
						if(!i && fputc('\n',fdump)<0) i=1;
					}
					if(!i && polygenic_flag) if(fprintf(fdump,"%x\n",loki->params.bv_iter)<0) i=1;
				}
				if(!i) {
					(void)fputs("[Seed] ",stdout);
					(void)fflush(stdout);
					if(fputs("LKSD:",fdump)<0) i=1;
					if(!i) i=bindumpseed(fdump);
				}
				if(!i) {
					(void)fputs("[Time] ",stdout);
					(void)fflush(stdout);
					if(fputs("LKTM:",fdump)<0) i=1;
					if(!i) i=txt_print_double(loki->sys.lktime.extra_time+difftime(time(0),loki->sys.lktime.start_time),fdump);
					if(!i) {
						ex1=loki->sys.lktime.extra_stime;
						ex2=loki->sys.lktime.extra_utime;
						tps=sysconf (_SC_CLK_TCK);
						errno=0;
						(void)times(&tbuf);
						if(errno) perror("print_end_time():");
						else {
							ex1+=(double)tbuf.tms_stime/(double)tps;
							ex2+=(double)tbuf.tms_utime/(double)tps;
						}
						if(!i && fputc(',',fdump)<0) i=1;
						if(!i) i=txt_print_double(ex1,fdump);
						if(!i && fputc(',',fdump)<0) i=1;
						if(!i) i=txt_print_double(ex2,fdump);
						if(!i && fputc('\n',fdump)<0) i=1;
					}
				}
				if(!i && fputs("Ldmp.end\n",fdump)<0) i=1;
			}
			if(fdump) {
				(void)fclose(fdump);
				if(i) (void)unlink(loki->names[LK_DUMPFILE]);
			}
		}
	} else i=1;
	if(i) (void)fputs("FAILED\n",stdout);
	else (void)fputs("OK\n",stdout);
	signal(SIGCHLD,SIG_DFL);
	while(waitpid(-1,&i,WNOHANG)>0);
}

