#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <errno.h>
#include <sys/wait.h>
#include <assert.h>

#include "string_utils.h"
#include "lk_malloc.h"
#include "lkgetopt.h"
#include "read_output.h"
#include "libhdr.h"
#include "utils.h"
#include "loki_compress.h"

struct t_bin {
	int lg;
	int bin;
};

int main(int argc,char *argv[])
{
	int i,j,k,k1,k2,kk[2],rep,nqtl,lg,n_lg=0,**dist,*dist1,**nq2,*dist_lnk,*tdist,*nbins[2];
	int max_qtl,tot_bins,it,it_real,nqtl1,nlqtl,sw,sw1,c,start=1,stop=0,skip=0,fflag;
	double z,z1,z2,*bf,bin_width=1.0,map_length,*chr_length,*tp,**bin[2],pp[2];
	char *p,*p1,*p2,*outfile=0,*posfile=0;
	FILE *fptr;
	struct ro_data data;
	struct ro_link *link,**lnk;
	string *s=0;
	void *tbuf=0;
	struct t_bin *t_bin[2],*t_bin1;
	struct loki loki;
	while((c=getopt(argc,argv,"hf:p:b:i:"))!=-1) switch(c) {
	 case 'f':
		outfile=strdup(optarg);
		break;
	 case 'p':
		posfile=strdup(optarg);
		break;
	 case 'b':
		bin_width=atof(optarg);
		if(bin_width<=0.0) {
			fprintf(stderr,"Illegal value %g for bin_width\n",bin_width);
			exit(EXIT_FAILURE);
		}
		break;
	 case 'i':
		start=(int)strtol(optarg,&p,10);
		k=0;
		if(p!=optarg && start<1) k=1;
		else {
			if(p==optarg) start=1;
			while(*p && isspace((int)*p)) p++;
			if(*p) {
				if(*p==',' || *p=='-' || *p==':') {
					if(*(++p)) {
						stop=strtol(p,&p1,10);
						if(p==p1 || stop<0) k=1;
						else {
							while(*p1 && isspace((int)*p1)) p1++;
							if(*p1) k=1;
						}
					}
				} else k=1;
			}
		}
		if(k) {
			fprintf(stderr,"Illegal argument '%s' for -i option\n",optarg);
			exit(EXIT_FAILURE);
		}
		if(stop && stop<start) {
			k=stop;
			stop=start;
			start=k;
		}
		break;
	 case 'h':
	 case '?':
		fputs("usage: count -i [start iteration][,:-][stop iteration] -f outfile -p posfile -b binsize\n",stderr);
		exit(EXIT_FAILURE);
	}
	init_lib_utils();
	loki.params.verbose_level=OUTPUT_NORMAL;
	loki.compress=init_compress();
	if(!outfile) outfile="loki.out";
 	/*if(!(fptr=open_readfile_and_check(outfile,&fflag,loki.compress))) exit(EXIT_FAILURE);*/
        fptr=open_readfile_and_check(outfile,&fflag,loki.compress);
        if(!(fptr))
        {
            exit(EXIT_FAILURE);
        }

	i=read_header(fptr,&data,&s);
	fclose(fptr);
	if(fflag) {
		while(waitpid(-1,&i,WNOHANG)>0);
	}
	if(!i) {
		for(i=0;i<=data.sex_map;i++) {
			link=data.link;
			n_lg=1;
			while(link) {
				n_lg++;
				link=link->next;
			}
		}
		/* Allocate storage */
		if(!n_lg) return 0;
		max_qtl=data.n_qtl[1];
		nbins[0]=lk_malloc(sizeof(int)*n_lg*2);
		nbins[1]=nbins[0]+n_lg;
		lnk=lk_malloc(sizeof(void *)*n_lg);
		lnk[0]=lk_malloc(sizeof(struct ro_link));
		lnk[0]->name="Unlinked";
		bf=lk_malloc(sizeof(double)*n_lg*2);
		chr_length=bf+n_lg;
		for(i=0;i<n_lg*2;i++) bf[i]=0.0;
		tot_bins=0;
		map_length=0.0;
		/* Check map length */
		for(i=0;i<=data.sex_map;i++) {
			link=data.link;
			j=0;
			z=0.0;
			while(link) {
				nbins[i][j]=(int)(0.5+(link->map_r2[i]-link->map_r1[i])/bin_width);
				tot_bins+=nbins[i][j];
				j++;
				chr_length[j]+=link->map_r2[i]-link->map_r1[i];
				z+=chr_length[j];
				lnk[j]=link;
				link=link->next;
			}
			if(z>=data.total_map[i]) data.total_map[i]=z+1.0;
			pp[i]=log(1.0-(double)bin_width/(double)data.total_map[i]);
			map_length+=data.total_map[i];
		}
		if(map_length) {
			for(z=1.0,i=1;i<n_lg;i++) {
				chr_length[i]/=map_length;
				z-=chr_length[i];
				chr_length[i]=log(1.0-chr_length[i]);
			}
			chr_length[0]=log(1.0-z);
		}
		dist1=lk_malloc(n_lg*sizeof(int));
		dist_lnk=lk_calloc((size_t)(max_qtl+1)*2,sizeof(int));
		tdist=dist_lnk+max_qtl+1;
		dist=lk_malloc(sizeof(void *)*n_lg*2);
		nq2=dist+n_lg;
		dist[0]=lk_calloc((size_t)(max_qtl+1)*n_lg,sizeof(int));
		for(j=1;j<n_lg;j++) dist[j]=dist[j-1]+max_qtl+1;
		nq2[0]=lk_calloc((size_t)n_lg*n_lg,sizeof(int));
		for(j=1;j<n_lg;j++) nq2[j]=nq2[j-1]+n_lg;
		bin[0]=lk_malloc(sizeof(void *)*n_lg*2);
		bin[1]=bin[0]+n_lg;
		tp=lk_malloc(sizeof(double)*tot_bins);
		for(i=0;i<tot_bins;i++) tp[i]=0.0;
		t_bin[0]=lk_malloc(sizeof(struct t_bin)*2*max_qtl);
		t_bin[1]=t_bin[0]+max_qtl;
		for(i=0;i<=data.sex_map;i++) for(j=0;j<n_lg;j++) {
			bin[i][j]=tp;
			tp+=nbins[i][j];
		}
		it=it_real=sw=sw1=0;
		nqtl1=nlqtl=-1;
		if(!posfile) posfile="loki.pos";
                fptr=open_readfile_and_check(posfile,&fflag,loki.compress);
		if(!(fptr)) exit(EXIT_FAILURE);
		do {
			s=fget_string(fptr,s,&tbuf);
			p1=p=get_cstring(s);
			while(*p1 && *(p1++)!=':');
			if(*p1) {
				rep=atoi(p1);
				skip=0;
				k=rep;
				if(start>it_real+rep) skip=1;
				else if(start>it_real+1) k+=it_real+1-start;
				if(stop) {
					if(stop<=it) skip=2;
					else if(stop<=it_real+rep) k-=it_real+rep-stop;
				}
				it_real+=rep;
				rep=k;
				if(rep>0 && !skip) {
					it+=rep;
					memset(dist1,0,sizeof(int)*n_lg);
					nqtl=0;
					p2=p1;
					kk[0]=kk[1]=0;
					while(p<p2 && *p!=':') {
						lg=(int)strtol(p,&p1,10);
						j=(int)strtol(p1,&p,10);
						assert(j);
						nqtl+=j;
						dist[lg][j]+=rep;
						dist1[lg]=j;
						for(k=lg-1;k>0;k--) if(dist1[k]) nq2[lg][k]+=rep;
						if(lg) {
							for(k1=0;k1<j;k1++) {
								for(k=0;k<=data.sex_map;k++) {
									z=strtod(p,&p1);
									p=p1;
									t_bin[k][kk[k]].bin=(int)((z-lnk[lg]->map_r1[k])/bin_width);
									t_bin[k][kk[k]++].lg=lg-1;
								}
							}
						}
					}
					tdist[nqtl]+=rep;
					k=nqtl-dist1[0];
					dist_lnk[k]+=rep;
					if(nqtl!=nqtl1) {
						nqtl1=nqtl;
						sw++;
					}
					if(k!=nlqtl) {
						nlqtl=k;
						sw1++;
					}
					if(map_length && nqtl) {
						for(lg=0;lg<n_lg;lg++) if(dist1[lg]) {
							z=1.0-exp(chr_length[lg]*(double)nqtl);
							bf[lg]+=(double)rep/z;
						}
						for(k=0;k<=data.sex_map;k++) {
							z=(double)rep/(1.0-exp(pp[k]*(double)nqtl));
							t_bin1=t_bin[k];
							for(k1=0;k1<kk[k];k1++,t_bin1++) {
								bin[k][t_bin1->lg][t_bin1->bin]+=z; 
							}
						}
					}
				}
			}
			if(!s->len) break;
		} while(skip<2);
		fclose(fptr);
		if(s) {
			free_string(s);
			s=0;
		}
		printf("Iterations: %d\n",it);
		if(it) {
			z=1.0/(double)it;
			printf("Switching proportion for QTL no.: %.5f\n",sw*z);
			printf("Switching proportion for linked QTLs: %.5f\n\n",sw1*z);
			for(i=max_qtl;i>=0;i--) if(tdist[i]) break;
			max_qtl=i;
			if(map_length) {
				fputs("Linkage group      Count   Prop. linked    BF\n",stdout);
				fputs("------------------------------------------------\n",stdout);
				for(lg=0;lg<n_lg;lg++) {
					for(i=j=0;i<max_qtl;i++) j+=dist[lg][i+1];
					printf("%-15s %8d     %.5f     %.5f",lnk[lg]->name,j,z*(double)j,z*bf[lg]);
					if(lg) {
						for(k=0;k<=data.sex_map;k++) {
							z1=0.0;
							k2=0;
							for(k1=0;k1<nbins[k][lg-1];k1++) {
								z2=bin[k][lg-1][k1];
								if(z2>z1) {
									z1=z2;
									k2=k1;
								}
							}
							printf("  %8.2f %7.2fcM",z1*z,lnk[lg]->map_r1[k]+0.5+bin_width*(double)k2);
						}
					}
					fputc('\n',stdout);
				}
				fputc('\n',stdout);
				if(n_lg>2) {
					for(lg=2;lg<n_lg;lg++) {
						for(k=1;k<lg;k++) {
							s=strprintf(s,"%s-%s",lnk[lg]->name,lnk[k]->name);
							printf("%-31s    %.5f\n",get_cstring(s),z*(double)nq2[lg][k]);
						}
					}
					fputc('\n',stdout);
				}
			}
			fputs("QTL number ",stdout);
			for(i=0;i<=max_qtl;i++) printf(" %6d",i);
			fputs("\n---------------",stdout);
			for(i=0;i<=max_qtl;i++) printf("-------");
			fputs("\nOverall        ",stdout);
			for(i=0;i<=max_qtl;i++) printf(" %.4f",z*(double)tdist[i]);
			fputs("\nLinked         ",stdout);
			for(i=0;i<=max_qtl;i++) printf(" %.4f",z*(double)dist_lnk[i]);
			fputc('\n',stdout);
			for(lg=0;lg<n_lg;lg++) {
				printf("%-15s",lnk[lg]->name);
				k=it;
				for(i=1;i<=max_qtl;i++) k-=dist[lg][i];
				dist[lg][0]=k;
				for(i=0;i<=max_qtl;i++) printf(" %.4f",z*(double)dist[lg][i]);
				if(!lg) {
					fputs("\n---------------",stdout);
					for(i=0;i<=max_qtl;i++) printf("-------");
				}
				fputc('\n',stdout);
			}
		}
		if(fflag) {
			while(waitpid(-1,&i,WNOHANG)>0);
		}
	} else {
		assert(s);
		fprintf(stderr,"Error reading file '%s' - %s\n",outfile,get_cstring(s));
	}
	if(s) free_string(s);
	if(tbuf) free_fget_buffer(&tbuf);
	free_string_lib();
	return i;
}
