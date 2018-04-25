#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <fcntl.h>
#include <sys/types.h>
#include <sys/uio.h>
#include <sys/wait.h>
#include <ctype.h>
#include <string.h>
#include <unistd.h>
#include <errno.h>
#include <signal.h>

int isfpdigit(char c)
{
	if(c=='.' || c=='-' || (c>='0' && c<='9')) return 1;
	return 0;
}

int main(int argc,char *argv[])
{
	int i,j,c,k,l,idx=1,nbins=20,err=0,ofreq=1,min_set=0,max_set=0;
	unsigned int *bins;
	double min,max,t,nn,tx,x,*binps=0;
	FILE *fptr;
	char buf[4096],*p,*p1,*Outfile=0,*Infile=0,*PFile=0;
	
	min=max=0.0;
	t=tx=0.0;
	if(argc==1) err=2;
	else while ((c=getopt(argc,argv,"hcfo:r:b:i:p:"))!=-1) switch(c) {
	 case 'r':
		p=optarg;
		x=strtod(p,&p1);
		if(p1!=p) {
			min=x;
			min_set=1;
		}
		if(*p1) {
			if(*p1!=':' && *p1!=',' && *p1!='-') err=1;
			p1++;
			x=strtod(p1,&p);
			if(p1!=p)	{
				max=x;
				max_set=1;
			}
			if(*p) err=1;
		}
		break;
	 case 'o':
		l=strlen(optarg);
		if(Outfile) free(Outfile);
		if(!(Outfile=(char *)malloc(l+1))) {
			(void)fprintf(stderr,"hist: out of memory error\n");
			exit(EXIT_FAILURE);
		}
		(void)strncpy(Outfile,optarg,l+1);
		break;
	 case 'p':
		l=strlen(optarg);
		if(PFile) free(PFile);
		if(!(PFile=(char *)malloc(l+1)))	{
			(void)fprintf(stderr,"hist: out of memory error\n");
			exit(EXIT_FAILURE);
		}
		(void)strncpy(PFile,optarg,l+1);
		break;
	 case 'b':
		nbins=atoi(optarg);
		break;
	 case 'i':
		idx=atoi(optarg);
		break;
	 case 'c':
		ofreq=0;
		break;
	 case 'f':
		ofreq=1;
		break;
	 case 'h':
	 case '?':
		err=2;
	}
	if(err==1) {
		(void)fprintf(stderr,"Invalid range specified.\n");
		(void)fprintf(stderr,"Use -r min:max or -r min or -r :max to set range\n");
		exit(EXIT_FAILURE);
	}
	if(err) {
		(void)fprintf(stderr,"usage: hist [-cfh] [-i index] [-b bins] [-r min:max] [-o file] [-p file] [file]\n");
		(void)fprintf(stderr," -c          Output counts.\n");
		(void)fprintf(stderr," -f          Output frequencies (default).\n");
		(void)fprintf(stderr," -h          Display this help.\n");
		(void)fprintf(stderr," -i index    Specify which column to read.\n");
		(void)fprintf(stderr," -b bins     Specify the number of bins.\n");
		(void)fprintf(stderr," -r min:max  Set range.\n");
		(void)fprintf(stderr," -o file     Send output to file.\n");
		(void)fprintf(stderr," -p file     Read in points for bins from file.\n");
		exit(EXIT_FAILURE);
	}
	if(optind<argc) {
		l=strlen(argv[optind]);
		if(!(Infile=(char *)malloc(l+1))) {
			(void)fprintf(stderr,"hist: out of memory error\n");
			exit(EXIT_FAILURE);
		}
		(void)strncpy(Infile,argv[optind],l+1);
	}
	if(nbins<2) {
		(void)fprintf(stderr,"No. bins must be > 1\n");
		exit(EXIT_FAILURE);
	}
	if(idx<1) {
		(void)fprintf(stderr,"Column index must be >0\n");
		exit(EXIT_FAILURE);
	}
	if(PFile) {
		if(max_set || min_set) (void)fprintf(stderr,"Warning - range is ognored when -p option specified\n");
		if(!(fptr=fopen(PFile,"r"))) {
			(void)fprintf(stderr,"Can't open '%s' for input\n",PFile);
			exit(EXIT_FAILURE);
		}
		j=0;
		while(fgets(buf,4096,fptr)) if(sscanf(buf,"%lf",&x)==1) {
			if(j) {
				if(x<=tx) {
					(void)fprintf(stderr,"Points must be ordered in file '%s'\n",PFile);
					exit(EXIT_FAILURE);
				}
			}
			j++;
			tx=x;
		}
		if(j<3) {
			(void)fprintf(stderr,"Insufficient points read in from file '%s'\n",PFile);
			exit(EXIT_FAILURE);
		}
		nbins=j-1;
		if(!(binps=(double *)malloc(sizeof(double)*j)))	{
			(void)fprintf(stderr,"hist: out of memory error\n");
			exit(EXIT_FAILURE);
		}
		(void)fseek(fptr,0,0);
		k=0;
		while(fgets(buf,4096,fptr)) if(sscanf(buf,"%lf",&x)==1) {
			binps[k++]=x;
			if(k==j) break;
		}
		(void)fclose(fptr);
		free(PFile);
		min=binps[0];
		max=binps[j-1];
		min_set=max_set=1;
	}
	if(Infile) {
		if(!(fptr=fopen(Infile,"r"))) {
			(void)fprintf(stderr,"Can't open '%s' for input\n",Infile);
			exit(EXIT_FAILURE);
		}
	} else {
		fptr=stdin;
		if(!(min_set && max_set)) {
			(void)fprintf(stderr,"The range *must* be completely specified if input is not from a file!\n");
			exit(EXIT_FAILURE);
		}
	}
	if(!(max_set && min_set)) {
		j=0;
		while(fgets(buf,4096,fptr))	{
			p=buf;
			for(i=0;i<idx;i++) {
				while(*p && isspace((int)*p)) p++;
				if(!*p) break;
				if((i+1)==idx) t=atof(p);
				while(*p && !isspace((int)*p)) p++;
			}
			if(i==idx) {
				if(!min_set) {
					if(j)	{
						if(t<min) min=t;
					} else min=t;
				}
				if(!max_set) {
					if(j)	{
						if(t>max) max=t;
					} else max=t;
				}
				j=1;
			}
		}
		(void)fseek(fptr,0,0);
	}
	tx=(max-min)/(double)(nbins);
	if(!(bins=(unsigned int *)malloc(sizeof(unsigned int)*nbins)))	{
		(void)fprintf(stderr,"hist: out of memory error\n");
		exit(EXIT_FAILURE);
	}
	for(i=0;i<nbins;i++) bins[i]=0;
	nn=0.0;
	while(fgets(buf,4096,fptr)) {
		p=buf;
		for(i=0;i<idx;i++) {
			while(*p && isspace((int)*p)) p++;
			if(!*p) break;
			if((i+1)==idx) t=atof(p);
			while(*p && !isspace((int)*p)) p++;
		}
		if(i==idx) {
			if(PFile) {
				if(t<=max && t>=min) {
					for(j=0;j<nbins;j++) if(t<binps[j+1]) break;
					if(j==nbins) j--;
					bins[j]++;
					nn+=1.0;
				}
			} else {
				if(t==max) j=nbins-1;
				else j=(int)((t-min)/tx);
				if(j>=0 && j<nbins) {
					bins[j]++;
					nn+=1.0;
				}
			}
		}
	}
	if(Infile) {
		(void)fclose(fptr);
		free(Infile);
	}
	if(Outfile)	{
		if(!(fptr=fopen(Outfile,"w"))) {
			(void)fprintf(stderr,"Can't open '%s' for output\n",Outfile);
			exit(EXIT_FAILURE);
		}
	} else fptr=stdout;
	for(i=0;i<nbins;i++)	{
		if(PFile) {
			x=.5*(binps[i+1]+binps[i]);
			tx=binps[i+1]-binps[i];
			if(ofreq) (void)fprintf(fptr,"%g %g %g\n",x,(double)bins[i]/(nn*tx),tx);
			else (void)fprintf(fptr,"%g %d %g\n",x,bins[i],tx);
		} else {
			x=min+tx*(.5+(double)i);
			if(ofreq) (void)fprintf(fptr,"%g %g %g\n",x,(double)bins[i]/(nn*tx),tx);
			else (void)fprintf(fptr,"%g %d %g\n",x,bins[i],tx);
		}
	}
	if(Outfile)	{
		(void)fclose(fptr);
		free(Outfile);
	}
	free(bins);
	return 0;
}
