#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
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
	int i,j,k,k1,kk[2],rep,nqtl,lg,n_lg=0,*nbins[2],*lnk_flag;
	int max_qtl,tot_bins,it,it_real,nqtl1,nlqtl,sw,sw1,c,start=1,stop=0,skip=0,fflag;
	int nognuplot=0,colourps=0,chrom=-1,lopflag=0,logflag=0,*match=0;
	double z,z1,bin_width=1.0,map_length,*tp,**bin[2],pp[2],il10;
	char *p,*p1,*p2,*outfile=0,*posfile=0,*dname=0,*psfile=0;
	FILE *fptr,*fdat[2];
	struct ro_data data;
	struct ro_link *link,**lnk;
	struct ro_mark *mark;
	string *s=0,**tmpnames[2],*tmpcontrol=0;
	void *tbuf=0;
	struct t_bin *t_bin[2],*t_bin1;
	struct loki loki;
	
	while((c=getopt(argc,argv,"hc:d:b:i:p:f:o:CqlL"))!=-1) switch(c) {
	 case 'f':
		outfile=strdup(optarg);
		break;
	 case 'l':
		lopflag=1;
		break;
	 case 'L':
		logflag=1;
		break;
	 case 'p':
		posfile=strdup(optarg);
		break;
	 case 'd':
		dname=strdup(optarg);
		break;
	 case 'o':
		psfile=strdup(optarg);
		break;
	 case 'c':
		chrom=atoi(optarg);
		if(chrom<1) {
			fprintf(stderr,"Illegal chromosome (%d)\n",chrom);
			exit(EXIT_FAILURE);
		}
		break;
	 case 'b':
		bin_width=atof(optarg);
		if(bin_width<=0.0) {
			fprintf(stderr,"Illegal value %g for bin_width\n",bin_width);
			exit(EXIT_FAILURE);
		}
		break;
	 case 'q':
		nognuplot=1;
		break;
	 case 'C':
		colourps=1;
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
		fputs("usage: dist -qlC -i [start iteration][,:-][stop iteration] -f outfile -p posfile -b binsize\n       -o psfile -d datafile -c chromosome\n",stderr);
		exit(EXIT_FAILURE);
	}
	init_lib_utils();
	loki.params.verbose_level=OUTPUT_NORMAL;
	loki.compress=init_compress();
	il10=1.0/log(10.0);
	if(!outfile) outfile="loki.out";
        fptr=open_readfile_and_check(outfile,&fflag,loki.compress);
 	if(!(fptr)) exit(EXIT_FAILURE);
	i=read_header(fptr,&data,&s);
	fclose(fptr);
	if(fflag) while(waitpid(-1,&i,WNOHANG)>0);
	if(!i) {
		for(i=k=0;i<=data.sex_map;i++) {
			link=data.link;
			n_lg=1;
			while(link) {
				if(chrom<0 || chrom==n_lg) k++;
				n_lg++;
				link=link->next;
			}
		}
		if(!k) {
			fputs("No linked chromosomes found\n",stderr);
			return 0;
		}
		lnk=lk_malloc(sizeof(void *)*n_lg);
		lnk[0]=lk_malloc(sizeof(struct ro_link));
		lnk[0]->name="Unlinked";
		link=data.link;
		j=0;
		while(link) {
			lnk[++j]=link;
			link=link->next;
		}
		if(lopflag) {
			match=lk_calloc((size_t)n_lg,sizeof(int));
			if(chrom>0) {
				/* Look for matched chromosome */
				s=strprintf(s,"%s'",lnk[chrom]->name);
				p=get_cstring(s);
				for(i=1;i<n_lg;i++) if(i!=chrom) {
					if(!strcmp(p,lnk[i]->name)) {
						match[chrom]=i;
						match[i]=-1;
						break;
					}
				}
				/* Not found - try reverse LOP */
				if(i==n_lg) {
					j=(int)strlen(lnk[chrom]->name);
					if(j && lnk[chrom]->name[j-1]=='\'') {
						p=strdup(lnk[chrom]->name);
						p[j-1]=0;
						for(i=1;i<n_lg;i++) if(i!=chrom) {
							if(!strcmp(p,lnk[i]->name)) {
								match[chrom]=i;
								match[i]=-1;
								break;
							}
						}
						free(p);
						if(i==n_lg) ABT_FUNC("Matched chromsome for LOP not found\n");
					}
				}
			} else {
				for(lg=1;lg<n_lg;lg++) if(!match[lg]) {
					s=strprintf(s,"%s'",lnk[lg]->name);
					p=get_cstring(s);
					for(i=1;i<n_lg;i++) if(i!=lg) {
						if(!strcmp(p,lnk[i]->name)) {
							match[lg]=i;
							match[i]=-1;
							break;
						}
					}
				}
			}
		}
		/* Allocate storage */
		max_qtl=data.n_qtl[1];
		nbins[0]=lk_malloc(sizeof(int)*n_lg*3);
		nbins[1]=nbins[0]+n_lg;
		lnk_flag=nbins[1]+n_lg;
		tot_bins=0;
		map_length=0.0;
		/* Check map length */
		for(i=0;i<=data.sex_map;i++) {
			link=data.link;
			j=0;
			z=0.0;
			while(link) {
				if(chrom<0 || chrom==j+1 || (lopflag && match[j+1])) {
					nbins[i][j]=(int)(0.5+(link->map_r2[i]-link->map_r1[i])/bin_width);
					tot_bins+=nbins[i][j];
					lnk_flag[++j]=1;
				} else {
					lnk_flag[++j]=0;
					nbins[i][j]=0;
				}
				z+=link->map_r2[i]-link->map_r1[i];
				lnk[j]=link;
				link=link->next;
			}
			if(z>=data.total_map[i]) data.total_map[i]=z+1.0;
			pp[i]=log(1.0-(double)bin_width/(double)data.total_map[i]);
			map_length+=data.total_map[i];
		}
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
					nqtl=0;
					p2=p1;
					kk[0]=kk[1]=0;
					while(p<p2 && *p!=':') {
						lg=(int)strtol(p,&p1,10);
						j=(int)strtol(p1,&p,10);
						assert(j);
						nqtl+=j;
						if(lg) {
							for(k1=0;k1<j;k1++) {
								for(k=0;k<=data.sex_map;k++) {
									z=strtod(p,&p1);
									p=p1;
									if(lnk_flag[lg]) {
										t_bin[k][kk[k]].bin=(int)((z-lnk[lg]->map_r1[k])/bin_width);
										t_bin[k][kk[k]++].lg=lg-1;
									}
								}
							}
						}
					}
					if(map_length && nqtl) {
						for(k=0;k<=data.sex_map;k++) {
							z=lopflag?(double)rep:(double)rep/(1.0-exp(pp[k]*(double)nqtl));
							t_bin1=t_bin[k];
							for(k1=0;k1<kk[k];k1++,t_bin1++) bin[k][t_bin1->lg][t_bin1->bin]+=z;
						}
					}
				}
			}
			if(!s->len) break;
		} while(skip<2);
		fclose(fptr);
		if(it) {
			z=1.0/(double)it;
			tmpnames[0]=lk_malloc(sizeof(void *)*n_lg*2);
			tmpnames[1]=tmpnames[0]+n_lg;
			for(lg=0;lg<n_lg*2;lg++) tmpnames[0][lg]=0;
			if(dname) {
				tmpcontrol=strprintf(s,"%s",dname);
				s=0;
				if(!(fptr=fopen(get_cstring(tmpcontrol),"w"))) {
					fprintf(stderr,"Couldn't open file %s for output: %s\n",get_cstring(tmpcontrol),strerror(errno));
					exit(EXIT_FAILURE);
				}
			} else {
				tmpcontrol=strprintf(s,"/tmp/loki_distXXXXXX");
				s=0;
				j=mkstemp(get_cstring(tmpcontrol));
				if(j<0 || !(fptr=fdopen(j,"w"))) {
					fprintf(stderr,"Couldn't open temp. file %s for output: %s\n",get_cstring(tmpcontrol),strerror(errno));
					exit(EXIT_FAILURE);
				}
			}
			if(psfile) {
				fprintf(fptr,"set output \"%s\"\n",psfile);
				fputs("set term postscript \"Times-Roman\" 14",fptr);
				if(colourps) fputs("color solid",fptr);
			}
			if(lopflag) fputs("\nset ylabel \"LOP1\"\n",fptr);
			else if(logflag) fputs("\nset ylabel \"log_10(L-Score)\"\n",fptr);
			else fputs("\nset ylabel \"L-Score\"\n",fptr);
			fputs("set xlabel \"Position\" 0,-2.5\nset nokey\nset grid\n",fptr);
			for(lg=1;lg<n_lg;lg++) if(lnk_flag[lg]) {
				if(lopflag && match[lg]<1) continue;
				if(!data.sex_map) fprintf(fptr,"set title \"%s -- %s\"\n",data.model,lnk[lg]->name);
				if(dname) {
					if(data.sex_map) {
						tmpnames[0][lg]=strprintf(0,"%s_%s_m.dat",dname,lnk[lg]->name);
						tmpnames[1][lg]=strprintf(0,"%s_%s_f.dat",dname,lnk[lg]->name);
					} else tmpnames[0][lg]=strprintf(0,"%s_%s.dat",dname,lnk[lg]->name);
					for(i=0;i<=data.sex_map;i++) {
						if(!(fdat[i]=fopen(get_cstring(tmpnames[i][lg]),"w"))) {
							fprintf(stderr,"Couldn't open file %s for output: %s\n",get_cstring(tmpnames[i][lg]),strerror(errno));
							exit(EXIT_FAILURE);
						}
					}
				} else {
					for(i=0;i<=data.sex_map;i++) {
						if(data.sex_map) fprintf(fptr,"set title \"%s -- %s - %s map\"\n",data.model,lnk[lg]->name,i?"Female":"Male");
						tmpnames[i][lg]=strprintf(0,"/tmp/loki_distXXXXXX");
						j=mkstemp(get_cstring(tmpnames[i][lg]));
						if(j<0 || !(fdat[i]=fdopen(j,"w"))) {
							fprintf(stderr,"Couldn't open temp. file %s for output: %s\n",get_cstring(tmpnames[i][lg]),strerror(errno));
							exit(EXIT_FAILURE);
						}
					}
				}
				if(lopflag) k=match[lg];
				for(i=0;i<=data.sex_map;i++) {
					z1=lnk[lg]->map_r1[i]+.5*bin_width;
					if(lopflag) {
						assert(nbins[i][lg-1]==nbins[i][k-1]);
						for(j=0;j<nbins[i][lg-1];j++,z1+=bin_width) fprintf(fdat[i],"%g %g\n",z1,il10*log((1.0+bin[i][lg-1][j])/(1.0+bin[i][k-1][j])));
					} else if(logflag) for(j=0;j<nbins[i][lg-1];j++,z1+=bin_width) fprintf(fdat[i],"%g %g\n",z1,il10*log((1.0+bin[i][lg-1][j])*z));
					else for(j=0;j<nbins[i][lg-1];j++,z1+=bin_width) fprintf(fdat[i],"%g %g\n",z1,bin[i][lg-1][j]*z);
					fclose(fdat[i]);
					j=0;
					mark=lnk[lg]->markers;
					while(mark) {
						if(!j++) fputs("set xtics rotate (",fptr);
						else fputs (",\\\n",fptr);
						fprintf(fptr,"\"%s\" %g",mark->name,mark->pos[i]);
						mark=mark->next;
					}
					if(j) fputs(")\n",fptr);
					fprintf(fptr,"set xrange[%g:%g]\n",lnk[lg]->map_r1[i],lnk[lg]->map_r2[i]);
					fprintf(fptr,"plot \"%s\" w l\n\n",get_cstring(tmpnames[i][lg]));
					if(!psfile) fputs("pause -1\n",fptr);
				}
			}
			fclose(fptr);
			if(!nognuplot) {
				s=strprintf(s,"gnuplot %s\n",get_cstring(tmpcontrol));
				i=system(get_cstring(s));
			}
			if(!dname) unlink(get_cstring(tmpcontrol));
			free_string(tmpcontrol);
			for(i=0;i<=data.sex_map;i++) {
				for(lg=1;lg<n_lg;lg++) if(tmpnames[i][lg]) {
					if(!dname) unlink(get_cstring(tmpnames[i][lg]));
					free_string(tmpnames[i][lg]);
				}
			}
			free(tmpnames[0]);
		}
		if(fflag) while(waitpid(-1,&i,WNOHANG)>0);
		i=0;
	} else {
		assert(s);
		fprintf(stderr,"Error reading file '%s' - %s\n",outfile,get_cstring(s));
	}
	if(tbuf) free_fget_buffer(&tbuf);
	if(s) free_string(s);
	free_string_lib();
	return i;
}
