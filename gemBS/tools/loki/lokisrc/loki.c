/****************************************************************************
*                                                                          *
*     Loki - Programs for genetic analysis of complex traits using MCMC    *
*                                                                          *
*             Simon Heath - University of Washington                       *
*                                                                          *
*                       July 1997                                          *
*                                                                          *
* loki.c:                                                                  *
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
#include <errno.h>
#include <sys/types.h>
#include <sys/time.h>
#include <signal.h>

#include "lkgetopt.h"
#include "version.h"
#include "ranlib.h"
#include "utils.h"
#include "lk_malloc.h"
#include "libhdr.h"
#include "string_utils.h"
#include "loki.h"
#include "loki_peel.h"
#include "loki_utils.h"
#include "gen_elim.h"
#include "locus.h"
#include "ped_utils.h"
#include "loki_compress.h"
#include "loki_tlmoves.h"
#include "get_peelseq.h"
#include "pseudo_chrom.h"
#include "ibs_check.h"
#include "recomb.h"

static struct loki loki;
static int no_report=0;

/* Called when program receives signal */
static void int_handler(int i)
{
	char *sigs[]={"SIGTERM","SIGQUIT","SIGHUP","SIGINT","SIGALRM","SIGVTALRM","SIGPROF","UNKNOWN"};
	int j,signals[]={SIGTERM,SIGQUIT,SIGHUP,SIGINT,SIGALRM,SIGVTALRM,SIGPROF};
	
	for(j=0;j<7;j++) if(i==signals[j]) break;
	if(loki.sys.catch_sigs && !loki.sys.sig_caught) {
		loki.sys.sig_caught=i;
		if(j>=4 && j<7) (void)fprintf(stderr,"Received signal %s: %s closing down.\n",sigs[j],LOKI_NAME);
		else (void)fprintf(stderr,"Received signal %s: please wait while %s closes down.\nRepeat signal to abort immediately.\n",sigs[j],LOKI_NAME);
		return;
	}
	(void)fprintf(stderr,"Exiting on signal %s (%d)\n",sigs[j],i);
	exit(i);
}

void FreeStuff(void)
{
	int i;
	
	if(!no_report) print_end_time();
	if(from_abt) return;
	if(!no_report && (loki.params.ranseed_set&3)==2) {
		if(loki.names[LK_SEEDFILE]) (void)writeseed(loki.names[LK_SEEDFILE],1);
		else (void)writeseed("seedfile",1);
	}
	for(i=0;i<LK_NUM_NAMES;i++) if(loki.names[i]) free(loki.names[i]);
	FreeRemem(loki.sys.FirstRemBlock);
	free(loki.data);
	free(loki.models);
	free(loki.markers);
	free(loki.pedigree);
	free_string_lib();
	return;
}

static void ignore_handler(int i)
{
	i=i;
	/* Do nothing */
}

int main(int argc,char *argv[])
{
	FILE *fptr;
	int i=0,c,j,read_dump_flag=0,append_output_flag=0,no_pfile=0;
	static int ibscheck,no_unrel,recomb_chk;
	struct sigaction s_action;
	char *tfile=0,*fname;
	struct itimerval itimerval;
	double sec,frac;
	
	struct option longopts[]={
	{"ibs_check",no_argument,&ibscheck,1},
	{"ibscheck",no_argument,&ibscheck,1},
	{"no_between",no_argument,&no_unrel,1},
	{"nobetween",no_argument,&no_unrel,1},
	{"recomb_check",no_argument,&recomb_chk,1},
	{0,0,0,0}
	};
	
	/* Record when we start */
	loki.sys.lktime.start_time=time(0);
	
	init_utils(&loki);
	
	/* Set default output level */
	loki.params.verbose_level=OUTPUT_NORMAL;
	
	/* Handle commandline options */
	while((c=getopt_long(argc,argv,"arVvQqed:p:f:",longopts,0))!=-1) {
		switch(c) {
			case 'f':
				i=strlen(optarg);
				if(tfile) free(tfile);
					tfile=lk_malloc(i+1);
				(void)strncpy(tfile,optarg,i+1);
				read_dump_flag=1;
				break;
			case 'd':
				if((j=set_file_dir(optarg))) {
					fprintf(stderr,"Error setting default file directory: %s\n",j==UTL_BAD_STAT?strerror(errno):utl_error(j));
					exit(EXIT_FAILURE);
				}
				break;
			case 'p':
				if((j=set_file_prefix(optarg))) {
					fprintf(stderr,"Error setting default file prefix: %s\n",utl_error(j));
					exit(EXIT_FAILURE);
				}
				break;
			case 'r':
				read_dump_flag=1;
				break;
			case 'a':
				append_output_flag=1;
				break;
			case 'v':
				loki.params.verbose_level=OUTPUT_VERBOSE;
				break;
			case 'q':
				loki.params.verbose_level=OUTPUT_QUIET;
				break;
			case 'Q':
				loki.params.verbose_level=NO_OUTPUT;
				break;
			case 'e':
				loki.params.error_correct_type=1;
				break;
			case 'V':
				(void)printf("%s\n",LOKI_NAME);
				return 0;
		}
	}
	if(optind>=argc) {
		if(!ibscheck) message(WARN_MSG,"Warning - no parameter file specified\n");
		no_pfile=1;
	}
	if(read_dump_flag) tfile=make_file_name(".dump");
	if(ibscheck || no_pfile) no_report=1;
	
	/* Turn off some output routines if not writing to a terminal */
	if(!isatty(STDOUT_FILENO)) loki.params.verbose_level|=NON_INTERACTIVE;
	/* Allocate memory handlers (for memory that is difficult to explicitly free) */
	if(!(loki.sys.FirstRemBlock=malloc(sizeof(struct remember)))) ABT_FUNC(MMsg);
	loki.sys.RemBlock=loki.sys.FirstRemBlock;
	loki.sys.RemBlock->pos=0;
	loki.sys.RemBlock->next=0;
	
	/* Ignore SIGPIPE signals while reading input files (otherwise get problems if reading through filter) */
	s_action.sa_handler=ignore_handler;
	s_action.sa_flags=0;
	(void)sigaction(SIGPIPE,&s_action,0L);
	
	/* If program terminates normally or via a signal I want it to go through
		* FreeStuff() so that Memory can be cleared so I can check for unfree()d blocks */
	if(atexit(FreeStuff)) ABT_FUNC("Unable to register exit function FreeStuff()\n");
	
	/* Find available compression routines */
	loki.compress=init_compress();
	
	/* Allocate parameter storage structures */
	if(!(loki.models=calloc((size_t)1,sizeof(struct lk_model)))) ABT_FUNC(MMsg);
	if(!(loki.markers=calloc((size_t)1,sizeof(struct lk_markers)))) ABT_FUNC(MMsg);
	loki.markers->marker=0;
	loki.markers->linkage=0;
	if(!(loki.pedigree=calloc((size_t)1,sizeof(struct lk_ped)))) ABT_FUNC(MMsg);
	loki.pedigree->comp_size=0;
	if(!(loki.data=malloc(sizeof(struct lk_data)))) ABT_FUNC(MMsg);
	loki.data->id_variable=0;
	loki.data->nonid_variable=0;
	loki.pedigree->populations=0;
	/* Read files from prep */
	if(ReadXMLFile(&loki)) ABT_FUNC(AbMsg);
	/* Record start time to logfile */
	print_start_time(LOKI_NAME,"a",no_report?0:loki.names[LK_LOGFILE],&loki.sys.lktime);
	/* Read in parameter file */
	if(!no_pfile) {
		if((fptr=fopen(argv[optind],"r"))) i=ReadParam(fptr,argv[optind],&loki);
		else {
			(void)printf("Couldn't open '%s' for input as parameter file\nAborting...\n",argv[optind]);
			exit(EXIT_FAILURE);
		}
		(void)fclose(fptr);
		if(i) exit(EXIT_FAILURE);
	}
	if(loki.params.pseudo_flag) handle_pseudochrom(&loki);
	/* Catch all interesting signals */
	s_action.sa_handler=int_handler;
	s_action.sa_flags=0;
	(void)sigaction(SIGINT,&s_action,0L);
	(void)sigaction(SIGHUP,&s_action,0L);
	(void)sigaction(SIGQUIT,&s_action,0L);
	(void)sigaction(SIGTERM,&s_action,0L);
	(void)sigaction(SIGALRM,&s_action,0L);
	(void)sigaction(SIGVTALRM,&s_action,0L);
	(void)sigaction(SIGPROF,&s_action,0L);
	
	/* Set up timer if requested, to limit runtime */
	if(loki.params.limit_time>0.0) {
		frac=modf(loki.params.limit_time,&sec);
		itimerval.it_value.tv_sec=(time_t)sec;
		itimerval.it_value.tv_usec=(time_t)(frac*1000000.0);
		/* Set it to trigger another timer event after 5 minutes in case
			* we get stuck in a single iteration */
		itimerval.it_interval.tv_sec=300;
		itimerval.it_interval.tv_usec=0;
		setitimer(loki.params.limit_timer_type,&itimerval,0);
	}
	
	/* Allow over-riding of dumpfile name from commandline */
	if(tfile) {
		if(loki.names[LK_DUMPFILE]) free(loki.names[LK_DUMPFILE]);
		loki.names[LK_DUMPFILE]=tfile;
	}
	
	/* Initialize RNG */
	if(!(loki.params.ranseed_set&2))	{
		if(!loki.names[LK_SEEDFILE]) {
			if(getseed("seedfile")) init_ranf(135421);
		} else if(getseed(loki.names[LK_SEEDFILE])) init_ranf(135421);
		loki.params.ranseed_set|=2;
	}
	if(!no_report && loki.names[LK_LOGFILE])	{
		fname=add_file_dir(loki.names[LK_LOGFILE]);
		if(fname) {
			fptr=fopen(loki.names[LK_LOGFILE],"a");
			if(fptr) {
				(void)fputs("\n     Starting state of RNG (seedfile):\n\n",fptr);
				(void)dumpseed(fptr,1);
				(void)fputc('\n',fptr);
				(void)fclose(fptr);
			}
			free(fname);
		}
	}
	
	/* Initialize trait loci stuff */
	if(!(ibscheck||recomb_chk)) init_trait_loci(&loki);
	
	/* Check pedigree, do locus specific pruning, and initialize families */
	if(check_ped(&loki)) exit(EXIT_FAILURE);
	for(i=0;i<loki.models->n_models;i++) if(loki.models->models[i].polygenic_flag) {
		Calculate_NRM(&loki);
		break;
	}
		build_families(&loki);
	/*  get_loops(&loki); */
	prune_ped(&loki);
	
	/* Initialize Genotype model (set up frequencies, allele recoding, peeling sequence etc. */
	Init_Markers(no_pfile,&loki); /* Sets up marker allele frequencies and maps */
	if(ibscheck) ibs_check(no_unrel,&loki);
	if(recomb_chk) recomb_check(&loki);
	if(ibscheck || recomb_chk) return 0;
	/* Perform genotype elimination and allele recoding */
	if(Gen_Elim(&loki)) {
		fprintf(stderr,"Mendelian errors reported: check in file %s for details\n",loki.names[LK_GENERRFILE]);
		ABT_FUNC("Aborting\n");
	}
	Get_Peel_Seq(&loki);
	if(!no_pfile) {
		AllocEffects(&loki);
		InitValues(&loki);
		peel_alloc(&loki);
		/* Allocate space for trait loci */
		TL_Alloc(&loki);
		/* Sample */	
		SampleLoop(read_dump_flag,append_output_flag,&loki);
	}
	/* Exit */
	return loki.sys.sig_caught;
}
