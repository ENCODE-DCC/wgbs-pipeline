/****************************************************************************
 *                                                                          *
 *     Loki - Programs for genetic analysis of complex traits using MCMC    *
 *                                                                          *
 *                 Simon Heath - CNG, Paris                                 *
 *                                                                          *
 *                       August 2002                                        *
 *                                                                          *
 * check_het.c:                                                             *
 *                                                                          *
 * Small stand-alone program to check heterozygosities of markers           *
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
#include <stdio.h>
#include <signal.h>

#include "version.h"
#include "libhdr.h"
#include "ranlib.h"
#include "utils.h"
#include "scan.h"
#include "y.tab.h"
#include "check_het.h"

loki_time lt;

static char *LogFile;

void print_version_and_exit(void)
{
	(void)printf("%s\n",HET_NAME);
	exit(EXIT_SUCCESS);
}

int main(int argc,char *argv[])
{
	FILE *fptr;
	int err;

	init_stuff(&LogFile);
	lt.start_time=time(0);
	if(argc<2) abt(__FILE__,__LINE__,"No control file specified\n");
	if((fptr=fopen(argv[1],"r"))) err=ReadControl(fptr,argv[1],&LogFile);
	else {
		(void)printf("Couldn't open '%s' for input as control file\nAborting...\n",argv[1]);
		exit(EXIT_FAILURE);
	}
	(void)fclose(fptr);
	if(err) {
		LogFile=0;
		exit(EXIT_FAILURE);
	}
	print_start_time(HET_NAME,"w",LogFile,&lt);
	if(getseed("seedfile",0)) init_ranf(135421);
	if(!scan_error_n) ReadData(LogFile);          /* Read in the datafile(s) and recode (where necessary) */
	if(!pruned_ped_size) {
		(void)printf("Zero size pedigree\nAborting...\n");
		exit(EXIT_FAILURE);
	}
	if(!scan_error_n && n_markers) test_het(LogFile);
	(void)writeseed("seedfile",1);
	if(family) free(family);
	free_nodes();
	if(scan_error_n) (void)fprintf(stderr,"Errors: %d  ",scan_error_n);
	if(scan_warn_n) (void)fprintf(stderr,"Warnings: %d  ",scan_warn_n);
	if(scan_error_n || scan_warn_n) (void)fprintf(stderr,"\n");
	return sig_caught;
}
