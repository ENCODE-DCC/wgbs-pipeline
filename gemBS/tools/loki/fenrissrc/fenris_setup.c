/****************************************************************************
 *                                                                          *
 *     Loki - Programs for genetic analysis of complex traits using MCMC    *
 *                                                                          *
 *             Simon Heath - Rockefeller University                         *
 *                                                                          *
 *                       October 1997                                       *
 *                                                                          *
 * loki_setup.c:                                                            *
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
#ifdef HAVE_LIMITS_H
#include <limits.h>
#endif
#include "ranlib.h"

#include "utils.h"
#include "loki.h"
#include "loki_peel.h"
#include "loki_output.h"
#include "mat_utils.h"
#include "sample_rand.h"
#include "fenris.h"

int censored_flag,censor_mode;
struct move_stats move_stats[N_MOVE_STATS];
static char *sexstr[2]={"female","male"};
#ifdef DEBUG
int *debug_level;
#endif

static double kosambi_to_haldane(double x)
{
	return x+50.0*log(cosh(.02*x));
}

static void marker_outside_error(int j,int k)
{
	int i;
	  
	i=marker[j].locus.link_group;
	(void)fputs("Marker",stderr);
	print_marker_name(stderr,j);
	if(sex_map) (void)fprintf(stderr," located outside of %s linkage group '%s'\n",sexstr[k],linkage[i].name);
	else (void)fprintf(stderr," located outside of linkage group '%s'\n",linkage[i].name);
}

void FenrisSetup(void)
{
	int i,j,k=-1,fx,k1,k1a,k2a,k2,er=0,comp,sire,dam,*temp_p,grp,*perm,*perm1,n_all;
	double p,min,max,*temp_dp,**temp_dpp,xx,x;
	struct Variable *group_var=0;
	struct Locus *loc,*loc1;
	
	if(sex_map) {
		for(j=i=0;i<ped_size;i++) {
			if(id_array[i].sex<1 || id_array[i].sex>2) {
				print_orig_id(stderr,i+1);
				(void)fputs(" has invalid sex information\n",stderr);
				er=1;
			}
		}
	}
	for(i=0;i<n_id_records;i++) if(id_variable[i].type&ST_GROUP) {
		group_var=id_variable+i;
		break;
	}
	for(i=j=0;i<ped_size;i++) {
		id_array[i].nkids=0;
		sire=id_array[i].sire;
		dam=id_array[i].dam;
		if(sire)	{
			id_array[sire-1].nkids++;
			j++;
		}
		if(dam) {
			id_array[dam-1].nkids++;
			j++;
		}
		id_array[i].flag=0;
	}
	if(j) {
		if(!(temp_p=malloc(sizeof(int)*j))) ABT_FUNC(MMsg);
		RemBlock=AddRemem(temp_p,RemBlock);
		for(i=0;i<ped_size;i++) {
			j=id_array[i].nkids;
			id_array[i].kids=temp_p;
			  if(j) {
				temp_p+=j;
				id_array[i].nkids=0;
			}
		}
		for(i=0;i<ped_size;i++) {
			sire=id_array[i].sire;
			dam=id_array[i].dam;
			if(sire) {
				j=id_array[sire-1].nkids++;
				id_array[sire-1].kids[j]=i;
			}
			if(dam) {
				j=id_array[dam-1].nkids++;
				id_array[dam-1].kids[j]=i;
			}
		}
	}
	for(k=0,comp=0;comp<n_comp;comp++)	{
		for(k2a=k1=0;k1<comp_size[comp];k1++) {
			k2=0;
			if(!id_array[k].sire) k2++;
			if(!id_array[k++].dam) k2++;
			comp_ngenes[comp]+=k2;
			if(k2<2) k2a++;
		}
		if(!k2a) {
			if(comp<n_comp-1) ABT_FUNC("Internal error - component has only singletons\n");
			singleton_flag=1;
			(void)printf("Last component consists of singletons\n");
		}
	}
	if(n_markers) {
		for(k=k2=i=0;i<n_markers;i++) {
			if(!marker[i].locus.n_alleles) continue;
			for(k1=grp=0;grp<n_genetic_groups;grp++) if(marker[i].count_flag[grp]) k1++;
			if(k1) {
				marker[i].counts=(void *)-1;
				k2+=k1*marker[i].locus.n_alleles;
				k++;
			} else marker[i].counts=0;
		}
		if(k) {
			if(!(temp_dpp=malloc(sizeof(void *)*k*n_genetic_groups))) ABT_FUNC(MMsg);
			RemBlock=AddRemem(temp_dpp,RemBlock);
			if(!(temp_dp=malloc(sizeof(double)*k2))) ABT_FUNC(MMsg);
			RemBlock=AddRemem(temp_dp,RemBlock);
			for(i=0;i<n_markers;i++) if(marker[i].counts) {
				marker[i].counts=temp_dpp;
				temp_dpp+=n_genetic_groups;
				for(grp=0;grp<n_genetic_groups;grp++) {
					if(marker[i].count_flag[grp]) {
						marker[i].counts[grp]=temp_dp;
						temp_dp+=marker[i].locus.n_alleles;
					} else marker[i].counts[grp]=0;
				}
			}
		}
	}
	for(i=0;i<n_markers;i++) if(marker[i].locus.n_alleles) {
		if(!extra_allele_flag) {
			j=marker[i].lumped;
			for(grp=0;grp<n_genetic_groups;grp++)
			  if(marker[i].freq_set[grp][j] && marker[i].locus.freq[grp][j]>0.0) break;
			if(grp==n_genetic_groups) {
				/* Lumped allele not used so remove */
				n_all=--marker[i].locus.n_alleles;
				for(comp=0;comp<n_comp;comp++) {
					if(marker[i].n_all1[comp]>n_all) marker[i].n_all1[comp]--;
				}
			}
		}
		for(grp=0;grp<n_genetic_groups;grp++) {
			if(marker[i].count_flag[grp]) {
				for(p=0.0,j=0;j<marker[i].locus.n_alleles;j++) {
					if(marker[i].freq_set[grp][j]) marker[i].counts[grp][j]=marker[i].locus.freq[grp][j];
					else marker[i].counts[grp][j]=1.0;
					p+=marker[i].counts[grp][j];
				}
				for(j=0;j<marker[i].locus.n_alleles;j++) {
					marker[i].locus.freq[grp][j]=marker[i].counts[grp][j]/p;
					marker[i].freq_set[grp][j]=0;
				}
			} else {
				for(p=0.0,j=0;j<marker[i].locus.n_alleles;j++) {
					if(!marker[i].freq_set[grp][j]) marker[i].locus.freq[grp][j]=0.1;
					p+=marker[i].locus.freq[grp][j];
				}
				if(fabs(p-1.0)>0.0001) {
					for(k=j=0;j<marker[i].locus.n_alleles;j++) if(!marker[i].freq_set[grp][j]) k++;
					if(k<j) {
						if(p<1.0) {
							for(j=0;j<marker[i].locus.n_alleles;j++) if(!marker[i].freq_set[grp][j]) marker[i].locus.freq[grp][j]=(1.0-p)/(double)k;
							p=1.0;
						} else {
							(void)fputs("Rescaling frequencies for marker ",stdout);
							print_marker_name(stdout,i);
							if(group_var) {
								(void)fputs(" in genetic group ",stdout);
								if(group_var->rec_flag==ST_STRING) (void)fputs(group_var->recode[grp].string,stdout);
								else (void)printf("%d",group_var->recode[grp].value);
							}
							for(j=0;j<marker[i].locus.n_alleles;j++) if(!marker[i].freq_set[grp][j]) marker[i].locus.freq[grp][j]=0.1;
							p+=.1*(double)k;
							(void)fputc('\n',stdout);
						}
					}
				}
				for(j=0;j<marker[i].locus.n_alleles;j++) marker[i].locus.freq[grp][j]/=p;
			}
		}
		if(!marker[i].pos_set) {
			er=1;
			(void)fputs("Position not set for marker ",stderr);
			print_marker_name(stderr,i);
			(void)fputc('\n',stderr);
		}
	}
	if(!er) {
		for(j=i=0;i<n_links;i++) if(linkage[i].n_markers>j) j=linkage[i].n_markers;
		if(j) {
			if(!(perm=malloc(sizeof(int)*2*j))) ABT_FUNC(MMsg);
			perm1=perm+j;
			for(i=0;i<n_links;i++) {
				get_locuslist(perm,i,&j,1);
				set_sort_sex(0);
				gnu_qsort(perm,(size_t)j,(size_t)sizeof(int),cmp_loci_pos);
				/* Check for zero recombination between markers */
				k1a=perm[0];
				for(k=1;k<j;k++) {
					k1=perm[k];
					for(k2=0;k2<=sex_map;k2++) {
						if(marker[k1].locus.pos[k2]==marker[k1a].locus.pos[k2]) {
							fputs("Zero ",stderr);
							if(sex_map) fputs(k2?"male":"female",stderr);
							fprintf(stderr,"recombination between markers ");
							print_marker_name(stderr,k1a);
							fputs(" and ",stderr);
							print_marker_name(stderr,k1);
							fputc('\n',stderr);
							er=1;
						}
					}
					k1a=k1;
				}
				if(er) continue;
				if(sex_map) {
					for(k=0;k<j;k++) perm1[k]=perm[k];
					set_sort_sex(1);
					gnu_qsort(perm1,(size_t)j,(size_t)sizeof(int),cmp_loci_pos);
					for(k=0;k<j;k++) if(perm[k]!=perm1[k]) {
						(void)fprintf(stderr,"Male and female marker maps for linkage group %s have different orders\n",linkage[i].name);
						er=1;
						break;
					}
					k=perm[0];
					p=marker[k].locus.pos[0];
					k1=(p>linkage[i].r1[0])?1:0;
					p=marker[k].locus.pos[1];
					k1^=(p>linkage[i].r1[1])?1:0;
					k=perm[j-1];
					p=marker[k].locus.pos[0];
					k1|=(p<linkage[i].r2[0])?2:0;
					p=marker[k].locus.pos[1];
					k1^=(p<linkage[i].r2[1])?2:0;
					if(!er && k1) {
						(void)fprintf(stderr,"Male and female marker maps for linkage group %s have different numbers of intervals\n",linkage[i].name);
						er=1;
					break;
					}
				}
				/* Convert from input Kosambi map to Haldane map (if necessary) */
				if(map_function==MAP_KOSAMBI) {
					for(k2=0;k2<=sex_map;k2++) {
						k1=perm[0];
						loc1=&marker[k1].locus;
						xx=loc1->pos[k2];
						if(linkage[i].range_set[k2]) {
							x=xx-linkage[i].r1[k2];
							if(x<0.0) {
								marker_outside_error(k1,k2);
								er=1;
							} else loc1->pos[k2]=linkage[i].r1[k2]+kosambi_to_haldane(x);
						}
						for(k=1;k<j;k++) {
							k1=perm[k];
							loc=&marker[k1].locus;
							x=kosambi_to_haldane(loc->pos[k2]-xx);
							x+=loc1->pos[k2];
							xx=loc->pos[k2];
							loc->pos[k2]=x;
							loc1=loc;
						}
						if(linkage[i].range_set[k2]) {
							x=linkage[i].r2[k2]-xx;
							if(x<0.0) {
								marker_outside_error(j,k2);
								er=1;
							} else linkage[i].r2[k2]=loc1->pos[k2]+kosambi_to_haldane(x);
						}
					}
					if(!sex_map) {
						for(k=0;k<j;k++) marker[perm[k]].locus.pos[1]=marker[perm[k]].locus.pos[0];
						linkage[i].r1[1]=linkage[i].r1[0];
						linkage[i].r2[1]=linkage[i].r2[0];
					}
				}
			}
			free(perm);
		}
	}
	if(!er && n_links) {
		for(i=0;i<n_links;i++) {
			linkage[i].sample_pos=0;
			fx=0;
			for(j=0;j<n_markers;j++) if(marker[j].locus.link_group==i) {
				if(marker[j].pos_set==2) {
					linkage[i].sample_pos=1;
					break;
				}
			}
			for(k=0;k<1+sex_map;k++) {
				if(!linkage[i].range_set[k]) {
					min=DBL_MAX;
					max=-DBL_MAX;
					for(j=0;j<n_markers;j++) if(marker[j].locus.link_group==i) {
						if(marker[j].pos_set==1) fx=1;
						if(marker[j].locus.pos[k]<min) min=marker[j].locus.pos[k];
						if(marker[j].locus.pos[k]>max) max=marker[j].locus.pos[k];
					}
					if(min==DBL_MAX) min=max=0.0;
					linkage[i].r1[k]=min;
					linkage[i].r2[k]=max;
					if(sex_map) (void)printf("Map range (%s) for linkage group '%s' set to %g-%gcM\n",sexstr[k],linkage[i].name,linkage[i].r1[k],linkage[i].r2[k]);
					else {
						(void)printf("Map range for linkage group '%s' set to %g-%gcM\n",linkage[i].name,linkage[i].r1[0],linkage[i].r2[0]);
						linkage[i].r1[1]=linkage[i].r1[0];
						linkage[i].r2[1]=linkage[i].r2[0];
					}
				} else {
					for(j=0;j<n_markers;j++) if(marker[j].locus.link_group==i) {
						if(marker[j].pos_set==1) fx=1;
						k1=0;
						if(marker[j].locus.pos[k]<linkage[i].r1[k]) k1=1;
						if(marker[j].locus.pos[k]>linkage[i].r2[k]) k1=1;
						if(k1) {
							marker_outside_error(j,k);
							er=1;
						}
					}
				}
			}
			if(!fx) { /* If no marker in linkage group has a fixed position, arbitrarily fix first marker */
				min=DBL_MAX;
				k1=-1;
				for(j=0;j<n_markers;j++) if(marker[j].locus.link_group==i) {
					if(marker[j].locus.pos[0]<min) {
						min=marker[j].locus.pos[0];
						k1=j;
					}
				}
				if(k1>=0) {
					marker[k1].pos_set=1;
					(void)fputs("Position for marker ",stdout);
					print_marker_name(stdout,k1);
					(void)printf(" fixed at %g",marker[k1].locus.pos[0]);
					if(sex_map) (void)printf(",%g",marker[k1].locus.pos[1]);
					(void)fputc('\n',stdout);
				}
			}
		}
		for(k=0;k<1+sex_map;k++) {
			p=0.0;
			for(i=0;i<n_links;i++) p+=linkage[i].r2[k]-linkage[i].r1[k];
			if(total_maplength[k]<0.0) {
				if(sex_map) (void)printf("No total %s map length set - no unlinked loci will be allowed\n",sexstr[k]);
				else (void)printf("No total map length set - no unlinked loci will be allowed\n");
				total_maplength[k]=p;
			} else if(p>=total_maplength[k]) {
				if(sex_map) (void)printf("Total %s map length <= sum of linkage group sizes - no unlinked loci will be allowed\n",sexstr[k]);
				else (void)printf("Total map length <= sum of linkage group sizes - no unlinked loci will be allowed\n");
				total_maplength[k]=p;
			}
			if(!sex_map) total_maplength[1]=total_maplength[0];
		}
	}
	if(er) exit(EXIT_FAILURE);
}

