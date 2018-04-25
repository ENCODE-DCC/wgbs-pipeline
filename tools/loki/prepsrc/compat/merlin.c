/****************************************************************************
 *                                                                          *
 *     Loki - Programs for genetic analysis of complex traits using MCMC    *
 *                                                                          *
 *                     Simon Heath - CNG, Evry                              *
 *                                                                          *
 *                         January 2003                                     *
 *                                                                          *
 * merlin.c:                                                                *
 *                                                                          *
 * MERLIN/QTDT compatability routines                                       *
 *                                                                          *
 * Copyright (C) Simon C. Heath 2003                                        *
 * This is free software.  You can distribute it and/or modify it           *
 * under the terms of the Modified BSD license, see the file COPYING        *
 *                                                                          *
 ****************************************************************************/

#include <config.h>
#include <stdlib.h>
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif
#include <string.h>
#ifdef USE_DMALLOC
#include <dmalloc.h>
#endif
#include <stdio.h>
#include <ctype.h>
#include <math.h>

#include "lkgetopt.h"
#include "version.h"
#include "utils.h"
#include "snprintf.h"
#include "scan.h"
#include "merlin.h"

#define BUFSIZE 512

int strip_vars;

static int check_if_header(char **p)
{
	int i;
	char **chk[3],**p1;
	static char *chk0[]={"CHR","CHROMOSOME",0},
	            *chk1[]={"MARKER","MARKER_NAME","MRK",0},
	            *chk2[]={"KOSAMBI","POS","POSITION","LOCATION","SEX_AVERAGED_POS","CM",0};

	chk[0]=chk0;
	chk[1]=chk1;
	chk[2]=chk2;
	i=0;
	while(p[i]) {
		p1=chk[i];
		while(*p1) {
			if(!(strcasecmp(p[i],*p1))) break;
			p1++;
		}
		if(!*p1) break;
		i++;
		if(i==3) break;
	}
	return i==3?1:0;
}

int process_qtdt(int argc,char *argv[],char **LogFile,int *error_check,loki_time *lt)
{
	static int noprune=0;
	size_t sz;
	int c,err,i,j,k,ncols=0,col_size=8,line,anon_var=0,nmark=0,chrom;
	FILE *fptr,*fin;
	char *datafile=0,*pedfile=0,*mapfile=0,*lfile=0,*missing=0,*p,*p1,*tmp_name;
	tokens *tok=0;
	static char buf[BUFSIZE];
	
	struct column {
		int c;
		char *name;
		double pos[3];
	} *cols=0;

	static struct option longopts[]={
		 {"noprune",no_argument,&noprune,1},
		 {"prune",no_argument,&noprune,0},
		 {"log",required_argument,0,'l'},
		 {0,0,0,0}
	};
	static char *fixed[]={
		  "Set skip_bad_ints 1\n",
		  "Set skip_bad_reals 1\n",
		  "Missing [\"RIF\"] \"x\"\n",
		  "Missing [\"RIF\"] \"X\"\n",
		  "Missing [\"F\"] \"0\"\n",
		  "Missing \"?\" sex_var\n",
		  "Pedigree family,id,father,mother\n",
		  "Sex sex_var \"1\",\"2\"\n",
		  "Sex sex_var \"M\",\"F\"\n",
		  "Sex sex_var \"m\",\"f\"\n",
		  0
	};
	
	while((c=getopt_long(argc,argv,"evm:d:p:x:l:X:",longopts,0))!=-1) switch(c) {
	 case 'e':
		error_check=0;
		break;
	 case 'v':
		print_version_and_exit();
		break; /* Previous command never returns */
	 case 'd':
		datafile=strdup(optarg);
		break;
	 case 'p':
		pedfile=strdup(optarg);
		break;
	 case 'm':
		mapfile=strdup(optarg);
		break;
	 case 'x':
		missing=strdup(optarg);
		break;
	 case 'l':
		lfile=strdup(optarg);
		break;
	 case 'X':
		fputs("-X option must occur as first argument\n",stderr);
		exit(EXIT_FAILURE);
	}
	if(!datafile) {
		fputs("No data file specified (use -d option)\n",stderr);
		exit(EXIT_FAILURE);
	}
	if(!pedfile) {
		fputs("No pedigree file specified (use -p option)\n",stderr);
		exit(EXIT_FAILURE);
	}
	/* Now we create a temporary command file which we can feed into the standard parse routines */
	if(!(fptr=tmpfile())) ABT_FUNC("Couldn't create temporary file\n"); 
	i=0;
	p=fixed[i];
	while(p) {
		(void)fputs(p,fptr);
		p=fixed[++i];
	}
	if(lfile) fprintf(fptr,"Log \"%s\"\n",lfile);
	else fputs("Log \"loki.log\"\n",fptr);
	if(missing) fprintf(fptr,"Missing [\"R\"] \"%s\"\n",missing);
	fprintf(fptr,"File [GS=' \\f\\r\\n\\t/'] \"%s\",family,id,father,mother,sex_var",pedfile);
	if(!(fin=fopen(datafile,"r"))) {
		perror("Couldn't open datafile for input");
		ABT_FUNC("Aborting\n");
	}
	if(!(cols=malloc(sizeof(struct column)*col_size))) ABT_FUNC(MMsg);
	line=0;
	while(fgets(buf,BUFSIZE,fin)) {
		line++;
		tok=tokenize(buf,0,tok);
		if(tok && tok->n_tok) {
			p=tok->toks[0];
			p1=tok->n_tok>1?tok->toks[1]:0;
			c=toupper((int)*p);
			if(i==1) fprintf(stderr,"Line %d: Item type %c has no name\n",line,c);
			(void)fputc(',',fptr);
			switch(c) {
			 case 'M':
				nmark++;
			 case 'A':
			 case 'T':
			 case 'C':
				if(ncols==col_size) {
					col_size+=2;
					if(!(cols=realloc(cols,sizeof(struct column)*col_size))) ABT_FUNC(MMsg);
				}
				if(p1) {
					sz=3+strlen(p1);
					if(!(cols[ncols].name=malloc(sz))) ABT_FUNC(MMsg);
					snprintf(cols[ncols].name,sz,"_%s_",p1);
				} else {
					i=1+(int)(log((double)++anon_var)/log(10.0));
					if(!(cols[ncols].name=malloc((size_t)(i+5)))) ABT_FUNC(MMsg);
					snprintf(cols[ncols].name,i+5,"var_%d",anon_var);
				}
				fputs(cols[ncols].name,fptr);
				cols[ncols++].c=c;
				break;
			 case 'Z':
			 case 'S':
			 case 'E':
				break;
			 default:
				fprintf(stderr,"Line %d: error in datafile - unknown item type (%c)\n",line,c);
				ABT_FUNC("Aborting\n");
			}
			if(c=='E') break;
		}
	}
	fclose(fin);
	fputc('\n',fptr);
	if(nmark) {
		fputs("Marker Loci",fptr);
		for(j=i=0;i<ncols;i++) if(cols[i].c=='M') {
			fputc(j++?',':' ',fptr);
			fputs(cols[i].name,fptr);
		}
		fputs("\nInteger",fptr);
		for(j=i=0;i<ncols;i++) if(cols[i].c=='M') {
			fputc(j++?',':' ',fptr);
			fputs(cols[i].name,fptr);
		}
		fputc('\n',fptr);
		if(mapfile) {
			chrom=0;
			if(!(fin=fopen(mapfile,"r"))) {
				perror("Couldn't open mapfile for input");
				ABT_FUNC("Aborting\n");
			}
			line=0;
			while(fgets(buf,BUFSIZE,fin)) {
				line++;
				tok=tokenize(buf,0,tok);
				if(tok) {
					i=tok->n_tok;
					if(!i) continue;
					if(i==3 || i==5) {
						if(check_if_header(tok->toks)) continue;
						c=toupper((int)*tok->toks[0]);
						if(c=='X') j=-1;
						else {
							j=strtol(tok->toks[0],&p,10);
							if(*p || j<1) {
								fprintf(stderr,"Line %d: error in mapfile - invalid chromosome '%s'\n",line,tok->toks[0]);
								continue;
							}
						}
						sz=3+strlen(tok->toks[1]);
						if(!(tmp_name=malloc(sz))) ABT_FUNC(MMsg);
						snprintf(tmp_name,sz,"_%s_",tok->toks[1]);
						for(k=0;k<ncols;k++) if(cols[k].c=='M') {
							if(!strcasecmp(tmp_name,cols[k].name)) break;
						}
						free(tmp_name);
						if(k==ncols) continue;
						cols[k].pos[0]=strtod(tok->toks[2],&p);
						if(*p) {
							fprintf(stderr,"Line %d: error in mapfile - invalid map position '%s'\n",line,tok->toks[2]);
							continue;
						}
						if(i==3) cols[k].c='m';
						else {
							cols[k].pos[2]=strtod(tok->toks[3],&p);
							if(*p) {
								fprintf(stderr,"Line %d: error in mapfile - invalid female map position '%s'\n",line,tok->toks[3]);
								continue;
							}
							cols[k].pos[1]=strtod(tok->toks[4],&p);
							if(*p) {
								fprintf(stderr,"Line %d: error in mapfile - invalid female map position '%s'\n",line,tok->toks[4]);
								continue;
							}
							cols[k].c='x';
						}
						if(j!=chrom) {
							if(j<0) fprintf(fptr,"Link [X] \'X\'");
							else fprintf(fptr,"Link \'%d\'",j);
							chrom=j;
						}
						cols[k].c='m';
						fputc(',',fptr);
						fputs(cols[k].name,fptr);
					} else {
						fprintf(stderr,"Line %d: error in mapfile - wrong number of columns (should be 3 or 5)\n",line);
						continue;
					}
				}
			}
			fputc('\n',fptr);
			fclose(fin);
			if(tok) free_tokens(tok);
		}
		for(k=j=i=0;i<ncols;i++) if(cols[i].c=='M') {
			if(!k) fputs("Link ",fptr);
			else fputc(',',fptr);
			fputs(cols[i].name,fptr);
		}
		fputc('\n',fptr);
		for(i=0;i<ncols;i++) if(cols[i].c=='m') {
			fprintf(fptr,"Position %s %g\n",cols[i].name,cols[i].pos[0]);
		}
		for(i=0;i<ncols;i++) if(cols[i].c=='x') {
			fprintf(fptr,"Position %s %g,%g,%g\n",cols[i].name,cols[i].pos[0],cols[i].pos[1],cols[i].pos[1]);
		}
	}
	for(i=0;i<ncols;i++) if(cols[i].c=='A') {
		p=cols[i].name;
		fprintf(fptr,"Discrete %s\n",p);
		fprintf(fptr,"Affected where (%s=='2' || %s=='D' || %s=='A' || %s=='Y')\n",p,p,p,p);
		fprintf(fptr,"Unaffected where (%s=='1' || %s=='U' || %s=='N')\n",p,p,p);
		break;
	}
	fseek(fptr,0,SEEK_SET);
	init_stuff(LogFile);
	strip_vars=1;
	err=ReadControl(fptr,"<TMP FILE>",LogFile);
	(void)fclose(fptr);
	if(!err) {
		print_start_time(PREP_NAME,"w",*LogFile,lt);
		if(!scan_error_n) ReadData(*LogFile);          /* Read in the datafile(s) and recode (where necessary) */
	}
	if(datafile) free(datafile);
	if(pedfile) free(pedfile);
	if(mapfile) free(mapfile);
	if(missing) free(missing);
	if(lfile) free(lfile);
	for(i=0;i<ncols;i++) free(cols[i].name);
	free(cols);
	return err;
}

