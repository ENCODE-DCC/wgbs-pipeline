/****************************************************************************
 *                                                                          *
 *     Loki - Programs for genetic analysis of complex traits using MCMC    *
 *                                                                          *
 *                 Simon Heath - CNG, Paris                                 *
 *                                                                          *
 *                       August 2002                                        *
 *                                                                          *
 * test_het.c:                                                              *
 *                                                                          *
 * Calculate heterozygosity, and test if different from expectation         *
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
#include <time.h>

#include "version.h"
#include "libhdr.h"
#include "ranlib.h"
#include "utils.h"
#include "scan.h"
#include "check_het.h"

static struct het_res *het_res;

static int cmp_p(const void *s1,const void *s2)
{
	double x1,x2;
	int i;
	
	i=*((const int *)s1);
	x1=het_res[i].p;
	i=*((const int *)s2);
	x2=het_res[i].p;
	if(x1<x2) return -1;
	if(x1>x2) return 1;
	return 0;
}

void test_het(char *Log)
{
	int locus,n_all,*counts,i,j,k,gt[2],hcount,hc1,*alist,n1;
	double nn,het,ehet,z,p1,p2;
	char *mname;
	time_t time1=0,time2;
	FILE *flog;
	
	fputs("Checking marker heterozygosity\n",stdout);
	if(!(het_res=malloc(sizeof(struct het_res)*n_markers))) ABT_FUNC(MMsg);
	for(i=locus=0;locus<n_markers;locus++) {
		n_all=markers[locus].element->n_levels;
		if(n_all>i) i=n_all;
	}
	if(i<2) {
		fputs("No segregating markers found\n",stdout);
		return;
	}
	if(!(counts=malloc(sizeof(int)*i))) ABT_FUNC(MMsg);
	for(k=i=0;i<ped_size;i++) if(id_array[i].haplo[0]) k+=2;
	if(!k) {
		fputs("No typed individuals found\n",stdout);
		free(counts);
		return;
	}
	if(!(alist=malloc(sizeof(int)*k))) ABT_FUNC(MMsg);
	sig_quiet=catch_sigs=1;
	for(locus=0;locus<n_markers;locus++) {
		sig_caught=0;
		mname=get_marker_name(locus);
		printf("Processing Marker %s:",mname);
		fflush(stdout);
		n_all=markers[locus].element->n_levels;
		if(n_all<2) {
			fputs("Skipped (< 2 alleles)\n",stdout);
			continue;
		}
		for(i=0;i<n_all;i++) counts[i]=0;
		for(nn=0.0,hcount=i=0;i<ped_size;i++) {
			if(id_array[i].haplo[0]) {
				for(k=0;k<2;k++) gt[k]=id_array[i].haplo[k][locus];
				if(gt[0] && gt[1]) {
					counts[gt[0]-1]++;
					counts[gt[1]-1]++;
					if(gt[0]==gt[1]) hcount++;
					nn+=2.0;
				}
			}
		}
		het=1.0-(double)(hcount*2)/nn;
		ehet=1.0;
		for(i=0;i<n_all;i++) {
			z=(double)counts[i]/nn;
			ehet-=z*z;
		}
		printf(" Het = %g, EHet = %g ",het,ehet);
		fflush(stdout);
		for(i=k=0;i<n_all;i++) for(j=0;j<counts[i];j++) alist[k++]=i;
		p1=p2=0.0;
		for(i=0;i<1000000;i++) {
			hc1=0;
			n1=(int)nn;
			while(n1) {
				j=(int)(safe_ranf()*n1);
				gt[0]=alist[j];
				alist[j]=alist[--n1];
				alist[n1]=gt[0];
				j=(int)(safe_ranf()*n1);
				gt[1]=alist[j];
				alist[j]=alist[--n1];
				alist[n1]=gt[1];
				if(gt[0]==gt[1]) hc1++;
			}
			if(hc1<hcount) p1++;
			else if(hc1>hcount) p2++;
			if(sig_caught || (i>=2000 && p1>=100 && p2>=100)) break;
		}
		if(sig_caught) fputs("\b\b",stdout);
		if(het>ehet) {
			fputs("p(Het>obs)=",stdout);
			z=p1;
			het_res[locus].flag=1;
		} else {
			fputs("p(Het<obs)=",stdout);
			z=p2;
			het_res[locus].flag=0;
		}
		if(z<100) {
			fputc('~',stdout);
			het_res[locus].flag|=2;
		}
		printf("%g\n",z/(double)i);
		het_res[locus].het=het;
		het_res[locus].ehet=ehet;
		het_res[locus].p=z/(double)i;
		het_res[locus].n=i;
		free(mname);
		if(sig_caught) {
			time2=time(0);
			z=difftime(time2,time1);
			if(z<.5) break;
			time1=time2;
		}
	}
	free(alist);
	k=locus;
	if(locus<n_markers) k++;
	if(k && Log && (flog=fopen(Log,"a"))) {
		(void)fprintf(flog,"\n***************** Heterozygosity Tests ******************\n\n");
	   if(!(alist=malloc(sizeof(int)*k))) ABT_FUNC(MMsg);
		for(i=0;i<k;i++) alist[i]=i;
		gnu_qsort(alist,(size_t)k,sizeof(int),cmp_p);
		(void)fputs("     Marker           Het      EHet      p-value    No. samples\n",flog);
		(void)fputs("     ------           ---      ----      -------    -----------\n",flog);
		for(i=0;i<k;i++) {
			locus=alist[i];
			mname=get_marker_name(locus);
			j=20-fprintf(flog,"     %s",mname);
			while((j--)>0) (void)fputc(' ',flog);
			fprintf(flog,"%.4f    %.4f    %c%-12g        %d\n",het_res[locus].het,het_res[locus].ehet,het_res[locus].flag&2?'~':' ',het_res[locus].p,het_res[locus].n);
			free(mname);
		}
		free(alist);
		fclose(flog);
	}
	free(counts);
	free(het_res);
}
