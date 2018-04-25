/****************************************************************************
 *                                                                          *
 *     Loki - Programs for genetic analysis of complex traits using MCMC    *
 *                                                                          *
 *                        Simon Heath - CNG                                 *
 *                                                                          *
 *                           June 2004                                      *
 *                                                                          *
 * write_xml_dump:                                                          *
 *                                                                          *
 * Routines for write XML format dump file                                  *
 *                                                                          *
 * Copyright (C) Simon C. Heath 2004                                        *
 * This is free software.  You can distribute it and/or modify it           *
 * under the terms of the Modified BSD license, see the file COPYING        *
 *                                                                          *
 ****************************************************************************/

#include <config.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <sys/times.h>
#include <errno.h>
#include <sys/wait.h>
#include <signal.h>
#include <assert.h>

#include "version.h"
#include "utils.h"
#include "string_utils.h"
#include "loki.h"
#include "libhdr.h"
#include "loki_peel.h"
#include "loki_output.h"
#include "sample_rand.h"
#include "loki_ibd.h"
#include "loki_compress.h"
#include "loki_dump.h"

void write_xml_dump(int lp,int lp1,int n_ibd _U_,long old_pos _U_,int flag _U_,int analysis,struct loki *loki)
{
  int i=0,j;
  struct lk_compress *lkc;
  struct lk_param *lkp;
  char *filter;
  FILE *fdump;
	
  lkc=loki->compress;
  filter=loki->names[LK_FILTER];
  if(!filter && (lkc->default_compress<COMPRESS_NONE)) {
    j=lkc->default_compress;
    filter=lkc->comp_path[j][j==COMPRESS_ZIP?1:0];
  }
  if(!loki->names[LK_DUMPFILE]) loki->names[LK_DUMPFILE]=make_file_name(".dump");
  if(loki->names[LK_DUMPFILE]) {
    j=loki->sys.syst_var[SYST_BACKUPS].flag?loki->sys.syst_var[SYST_BACKUPS].data.value:1;
    if(j) i=mkbackup(loki->names[LK_DUMPFILE],j);
    if(!i) {
      if(filter) {
	errno=0;
	i=child_open(WRITE,loki->names[LK_DUMPFILE],filter);
	fdump=fdopen(i,"w");
	if(errno && errno!=ESPIPE) i=1;
	else i=0;
	errno=0;
      } else fdump=fopen(loki->names[LK_DUMPFILE],"w");
      if(!fdump || i) (void)fprintf(stderr,"[%s:%d] %s(): File Error.  Couldn't open '%s' for output\n",__FILE__,__LINE__,__func__,loki->names[LK_DUMPFILE]);
      else {
	lkp=&loki->params;
	if(fputs("<?xml version='1.0' encoding='UTF-8'?>\n",fdump)==EOF) i=1;
	if(!i && fprintf(fdump,"<loki program='%s'>\n",LOKI_NAME)<0) i=1; /* Start tag */
	/* Run Parameters */
	if(!i && fprintf(fdump," <group type='run_parameters'>\n")<0) i=1;
	if(!(analysis&IBD_ANALYSIS)) {
	  if(lp) {
	    if(!i && fprintf(fdump,"  <var id='iter'>%d",lp)<0) i=1;
	    if(!i && lp1 && fprintf(fdump,",%d",lp1)<0) i=1;
	    if(!i && lkp->bv_iter && fprintf(fdump,",%d",lkp->bv_iter)<0) i=1;
	    if(!i && fputs("</var>>\n",fdump)==EOF) i=1;
	  }
	  if(lkp->num_iter && !i && fprintf(fdump,"  <var id='num_iter'>%d</var>\n",lkp->num_iter)<0) i=1;
	  if(!i && fprintf(fdump,"  <var id='sample_from'>%d,%d</var>\n",lkp->sample_from[0],lkp->sample_from[1])<0) i=1;
	  if(!i && fprintf(fdump,"  <var id='sample_freq'>%d,%d</var>\n",lkp->sample_freq[0],lkp->sample_freq[1])<0) i=1;
	  if(analysis && !i && fprintf(fdump,"  <var id='analysis'>%d</var>\n",analysis)<0) i=1;
	  if(lkp->dump_freq && !i && fprintf(fdump,"  <var id='dump_freq'>%d</var>\n",lkp->dump_freq)<0) i=1;
	  if(lkp->output_type!=DEFAULT_OUTPUT_TYPE && !i && fprintf(fdump,"  <var id='output_type'>%d</var>\n",lkp->output_type)<0) i=1;
	  if(lkp->output_haplo && !i && fprintf(fdump,"  <var id='output_haplo'>%d</var>\n",lkp->output_haplo)<0) i=1;
	  if(lkp->si_mode && !i && fprintf(fdump,"  <var id='si_mode'>%d</var>\n",lkp->si_mode)<0) i=1;
	  if(lkp->debug_level && !i && fprintf(fdump,"  <var id='debug_limit'>%d</var>\n",lkp->debug_level)<0) i=1;
	  if(lkp->prune_option!=DEFAULT_PRUNE_OPTION && !i && fprintf(fdump,"  <var id='prune_option'>%d</var>\n",lkp->prune_option)<0) i=1;
	  if(lkp->peel_trace && !i && fprintf(fdump,"  <var id='peel_trace'>%d</var>\n",lkp->peel_trace)<0) i=1;
	  if(lkp->verbose_level!=OUTPUT_NORMAL && !i && fprintf(fdump,"  <var id='verbose_level'>%d</var>\n",lkp->verbose_level)<0) i=1;
	  if(lkp->ibd_mode && !i && fprintf(fdump,"  <var id='ibd_mode'>%d</var>\n",lkp->ibd_mode)<0) i=1;
	  if(lkp->compress_ibd && !i && fprintf(fdump,"  <var id='compress_ibd'>%d</var>\n",lkp->compress_ibd)<0) i=1;
	  if(lkp->est_aff_freq && !i && fprintf(fdump,"  <var id='est_aff_freq'>%d</var>\n",lkp->est_aff_freq)<0) i=1;
	  if(lkp->ranseed_set && !i && fprintf(fdump,"  <var id='ranseed_set'>%d</var>\n",lkp->ranseed_set)<0) i=1;
	  if(lkp->map_function && !i && fprintf(fdump,"  <var id='ranseed_set'>%d</var>\n",lkp->map_function)<0) i=1;
	  if(lkp->pseudo_flag) {
	    if(!i && fprintf(fdump,"  <var id='pseudo_flag'>%d</var>\n",lkp->pseudo_flag)<0) i=1;
	    if(lkp->pseudo_freq && !i && fprintf(fdump,"  <var id='pseudo_freq'>%d</var>\n",lkp->pseudo_freq)<0) i=1;
	    if(lkp->pseudo_start!=1 && !i && fprintf(fdump,"  <var id='pseudo_start'>%d</var>\n",lkp->pseudo_start)<0) i=1;
	  }
	  if(lkp->error_correct_type && !i && fprintf(fdump,"  <var id='error_correct_type'>%d</var>\n",lkp->error_correct_type)<0) i=1;
	  if(lkp->genv_out && !i && fprintf(fdump,"  <var id='stat5_out'>%d</var>\n",lkp->genv_out)<0) i=1;
	  if(lkp->n_tloci) {
	    if(lkp->max_tloci!=DEFAULT_MAX_TLOCI && !i && fprintf(fdump,"  <var id='max_tloci'>%d</var>\n",lkp->max_tloci)<0) i=1;
	    if(lkp->min_tloci && !i && fprintf(fdump,"  <var id='max_tloci'>%d</var>\n",lkp->min_tloci)<0) i=1;
	    if(lkp->start_tloci!=lkp->min_tloci && !i && fprintf(fdump,"  <var id='start_tloci'>%d</var>\n",lkp->start_tloci)<0) i=1;
	  }
	  if(!i && lkp->limit_time>0.0) {
	    if(fprintf(fdump,"  <var id='limit_time'>%g</var>\n",lkp->limit_time)<0) i=1;
	    if(!i && fprintf(fdump,"  <var id='limit_time_type'>%d</var>\n",lkp->limit_timer_type)<0) i=1;
	  }
	}
	if(!i && fputs(" </group>\n",fdump)==EOF) i=1;
	if(!i && fputs("</loki>\n",fdump)==EOF) i=1; /* Closing tag */
	fputc('\n',fdump);
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
