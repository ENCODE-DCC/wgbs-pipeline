/****************************************************************************
 *                                                                          *
 *     Loki - Programs for genetic analysis of complex traits using MCMC    *
 *                                                                          *
 *                     Simon Heath - CNG, France                            *
 *                                                                          *
 *                           August 2002                                    *
 *                                                                          *
 * fenris.c:                                                                *
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
#include <sys/types.h>
#ifdef USE_MONITOR
#include <sys/ipc.h>
#include <sys/shm.h>
#endif
#include <sys/time.h>
#include <signal.h>

#include "version.h"
#include "lkgetop.h"
#include "ranlib.h"
#include "utils.h"
#include "libhdr.h"
#include "loki.h"
#include "loki_peel.h"
#include "fenris_peel.h"
#include "fenris.h"

unsigned int RunID;
int use_student_t,analysis,tau_mode,no_overdominant;
double *tau_beta,*tau,lm_ratio;
int catch_sigs,sig_caught,multiple_rec;
int map_function,singleton_flag,est_aff_freq,*ibd_mode;
char *Filter,*Seedfile;
struct remember *RemBlock;
static struct remember *FirstRemBlock;
loki_time lt;
static int ranseed_set,*n_bits,*fixed;
static char *LogFile;

int trace_level=0;

/* Called when program receives signal */
static void int_handler(int i)
{
	char *sigs[]={"SIGTERM","SIGQUIT","SIGHUP","SIGINT","SIGALRM","SIGVTALRM","SIGPROF","UNKNOWN"};
	int j,signals[]={SIGTERM,SIGQUIT,SIGHUP,SIGINT,SIGALRM,SIGVTALRM,SIGPROF};
	
	for(j=0;j<7;j++) if(i==signals[j]) break;
	if(catch_sigs && !sig_caught) {
		sig_caught=i;
		if(j>=4 && j<7) (void)fprintf(stderr,"Received signal %s: %s closing down.\n",sigs[j],FENRIS_NAME);
		else (void)fprintf(stderr,"Received signal %s: please wait while %s closes down.\nRepeat signal to abort immediately.\n",sigs[j],FENRIS_NAME);
		return;
	}
	(void)fprintf(stderr,"Exiting on signal %s (%d)\n",sigs[j],i);
	exit(i);
}

void FreeStuff(void)
{
	int i;

	if(Seedfile) free(Seedfile);
	if(LogFile) free(LogFile);
	if(Filter) free(Filter);
	if(Freqfile) free(Freqfile);
	if(Haplofile) free(Haplofile);
	if(Outputfile) free(Outputfile);
	if(OutputIBDfile) free(OutputIBDfile);
	if(OutputPosfile) free(OutputPosfile);
	if(id_array) free(id_array);
	if(founder_flag) {
 		if(founder_flag[0]) free(founder_flag[0]);
		free(founder_flag);
	}
	if(marker) {
		if(marker[0].mterm) free(marker[0].mterm);
		for(i=0;i<n_markers;i++) {
			if(marker[i].locus.variance) free(marker[i].locus.variance);
		}
		free(marker);
	}
	if(allele_trans) {
		if(allele_trans[0]) {
			if(allele_trans[0][0]) free(allele_trans[0][0]);
			free(allele_trans[0]);
		}
		free(allele_trans);
	}
	if(linkage) {
		for(i=0;i<n_links;i++) {
			if(linkage[i].ibd_list) {
				free(linkage[i].ibd_list->pos);
				free(linkage[i].ibd_list);
			}
		}
		free(linkage);
	}
	if(comp_size) free(comp_size);
	if(id_variable) free(id_variable);
	else if(nonid_variable) free(nonid_variable);
	if(n_bits) free(n_bits);
	if(fixed) free(fixed);
	FreeRemem(FirstRemBlock);
	return;
}

void ignore_handler(int i)
{
	i=i;
	/* Do nothing */
}

int main(int argc,char *argv[])
{
	FILE *fptr;
	int i=0,j,max_bits=22,c,errcode=0,fx;
	struct sigaction s_action;
	struct itimerval itimerval;
	double sec,frac,error_probs[5],*p;

	lt.start_time=time(0);
	while((c=getopt(argc,argv,"b:t:v"))!=-1) switch(c) {
	 case 'b':
		j=atoi(optarg);
		if(j<4) abt(__FILE__,__LINE__,"Not going to do much with %d bi%s...\n",j,j==1?"t":"ts");
		max_bits=j;
		break;
	 case 't':
		j=atoi(optarg);
		trace_level=j;
		break;
	 case 'v':
		(void)printf("%s\n",FENRIS_NAME);
		return 0;
	}
	if(optind>=argc) abt(__FILE__,__LINE__,"No parameter file specified\n");
	/* If program terminates normally or via a signal I want it to go through
	 * FreeStuff() so that Memory can be cleared and I can check for unfreed blocks
	 */
	if(!(FirstRemBlock=malloc(sizeof(struct remember)))) ABT_FUNC(MMsg);
	RemBlock=FirstRemBlock;
	RemBlock->pos=0;
	RemBlock->next=0;
	s_action.sa_handler=ignore_handler;
	s_action.sa_flags=0;
	(void)sigaction(SIGPIPE,&s_action,0L);
	ReadBinFiles(&LogFile,1);
	if(models) {
		for(i=0;i<n_models;i++) if(models[i].term) free(models[i].term);
		free(models);
		models=0;
	}
	AllocFenrisStruct();
	print_start_time(FENRIS_NAME,"a",LogFile,&lt);
 	if((fptr=fopen(argv[optind],"r"))) i=ReadParam(fptr,argv[optind],&ranseed_set);
	else {
		(void)printf("Couldn't open '%s' for input as parameter file\nAborting...\n",argv[optind]);
		exit(EXIT_FAILURE);
	}
	(void)fclose(fptr);
	if(Dumpfile) free(Dumpfile);
	if(Polyfile) free(Polyfile);
	if(i) exit(EXIT_FAILURE);
	if(atexit(FreeStuff)) ABT_FUNC("Unable to register exit function FreeStuff()\n");
	if(atexit(print_end_time)) ABT_FUNC("Unable to register exit function print_end_time()\n");
	s_action.sa_handler=int_handler;
	s_action.sa_flags=0;
	(void)sigaction(SIGINT,&s_action,0L);
	(void)sigaction(SIGHUP,&s_action,0L);
	(void)sigaction(SIGQUIT,&s_action,0L);
	(void)sigaction(SIGTERM,&s_action,0L);
	(void)sigaction(SIGALRM,&s_action,0L);
	(void)sigaction(SIGVTALRM,&s_action,0L);
	(void)sigaction(SIGPROF,&s_action,0L);
	if(!(ranseed_set&2))	{
		if(!Seedfile) {
			if(getseed("seedfile",0)) init_ranf(135421);
		} else if(getseed(Seedfile,0)) init_ranf(135421);
		ranseed_set|=2;
	}
	if(LogFile)	{
		fptr=fopen(LogFile,"a");
		if(fptr) {
			(void)fputs("\n     Starting state of RNG (seedfile):\n\n",fptr);
			(void)dumpseed(fptr,1);
			(void)fputc('\n',fptr);
			(void)fclose(fptr);
		}
	}
	FenrisSetup();
	if(n_comp) {
		if(limit_time>0.0) {
			frac=modf(limit_time,&sec);
			itimerval.it_value.tv_sec=(time_t)sec;
			itimerval.it_value.tv_usec=(time_t)(frac*1000000.0);
			/* Set it to trigger another timer event after 5 minutes in case
			 * we get stuck in the peel operation */
			itimerval.it_interval.tv_sec=300;
			itimerval.it_interval.tv_usec=0;
			setitimer(limit_timer_type,&itimerval,0);
		}
		error_probs[0]=.0125;
		error_probs[1]=.0075;
		error_probs[2]=.0050;
		error_probs[3]=.0100;
		error_probs[4]=.0025;
		error_probs[0]=0.02;
		estimate_freq(PEN_MODEL_EQUAL,error_probs);
		exit(0);
		if(!(n_bits=malloc(sizeof(int)*n_comp))) ABT_FUNC(MMsg);
		fx=count_bits(n_bits);
		if(fx) {
			if(!(fixed=malloc(sizeof(int)*ped_size))) ABT_FUNC(MMsg);
			for(j=0;j<ped_size;j++) fixed[j]=id_array[j].flag&24;
		}
		for(j=c=0;c<n_comp;c++) {
			if(n_bits[c]>max_bits) {
				fprintf(stderr,"Skipping component %d ",c+1);
				if(family_id) {
					fputc('(',stderr);
					(void)print_orig_family(stderr,comp_start[c]+1,0);
					fputc(')',stderr);
				}
				fprintf(stderr," - %d bits required, maximum bits set to %d\n",n_bits[c],max_bits);
			} else {
				j++;
/*				if(!(p=malloc(sizeof(double)*(1<<n_bits[c])))) ABT_FUNC(MMsg); */
				p=0;
				calc_pen(c,n_bits[c],0,PEN_MODEL_EQUAL,error_probs,p,fixed);
				free(p);
			}
		}
		if(!j) {
			fputs("No families analyzed\n",stderr);
			errcode=EXIT_FAILURE;
		}
		
	}
	return errcode?errcode:sig_caught;
}
