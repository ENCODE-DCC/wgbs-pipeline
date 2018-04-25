#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <ctype.h>
#include <string.h>
#include <float.h>

#define BLKSIZE 4096
static int bufsize,ncol,maxcol,**coldata,*ncols,*cols,blkptr,ldsize,rec,*nn;
static char *buf,*MMsg="Out of memory error\n";
static double *sum,*sum2,*min,*max,**linedata,*data;

static double addlog(double x1,double x2)
{
	double y,r;

	y=.5*(x1+x2);
	r=(log(cosh(x1-y))+y+log(2.0));
	return r;
}

static void abt(char *s) 
{
	(void)fputs(s,stderr);
	exit(EXIT_FAILURE);
}

static char *read_line(FILE *fptr) 
{
	char *p;
	int i=0;
	size_t l;
	
	for(;;) {
		p=fgets(buf+i,bufsize-i,fptr);
		if(p) {
			l=strlen(buf);
			if(buf[l-1]!='\n') {
				i=bufsize-1;
				bufsize*=2;
				if(!(buf=realloc(buf,(size_t)bufsize))) abt(MMsg);
			} else break;
		} else break;
	}
	return p;
}

static void add_data(double x,int col,int flag)
{
	double *data1;
	int *cols1;
	
	if(blkptr>=BLKSIZE) {
		if(col>=BLKSIZE) abt("No. columns > BLKSIZE\n");
		if(!(data1=malloc(sizeof(double)*BLKSIZE))) abt(MMsg);
		if(!(cols1=malloc(sizeof(int)*BLKSIZE))) abt(MMsg);
		(void)memcpy(data1,data+blkptr-col,col*sizeof(double));
		(void)memcpy(cols1,cols+blkptr-col,col*sizeof(int));
		linedata[rec]=data1;
		coldata[rec]=cols1;
		data=data1;
		cols=cols1;
		blkptr=col;
	}
	data[blkptr]=x;
	cols[blkptr++]=flag;
}

static void read_cols(FILE *fptr,int logopt,int start,int ppoint)
{
	int i,line=0,nc;
	char *p,*p1;
	double x;
	
	if(!buf) {
		/* Allocate initial buffers */
		bufsize=256;
		if(!(buf=malloc((size_t)bufsize))) abt(MMsg);
		maxcol=8;
		if(!(sum=malloc(sizeof(double)*maxcol))) abt(MMsg);
		if(!(sum2=malloc(sizeof(double)*maxcol))) abt(MMsg);
		if(!(nn=malloc(sizeof(int)*maxcol))) abt(MMsg);
		if(!(min=malloc(sizeof(double)*maxcol))) abt(MMsg);
		if(!(max=malloc(sizeof(double)*maxcol))) abt(MMsg);
		for(i=0;i<maxcol;i++) {
			sum[i]=sum2[i]=0.0;
			nn[i]=0;
			min[i]=DBL_MAX;
			max[i]=-DBL_MAX;
		}
		if(ppoint>=0.0) {
			if(!(data=malloc(sizeof(double)*BLKSIZE))) abt(MMsg);
			if(!(cols=malloc(sizeof(int)*BLKSIZE))) abt(MMsg);
			ldsize=32;
			if(!(linedata=malloc(sizeof(double *)*ldsize))) abt(MMsg);
			if(!(coldata=malloc(sizeof(int *)*ldsize))) abt(MMsg);
			if(!(ncols=malloc(sizeof(int *)*ldsize))) abt(MMsg);
			blkptr=0;
		}
	}
	while(read_line(fptr)) {
		line++;
		if(line<start) continue;
		nc=0;
		p=buf;
		if(ppoint>=0.0) {
			if(rec>=ldsize) {
				ldsize*=2;
				if(!(linedata=realloc(linedata,sizeof(double *)*ldsize))) abt(MMsg);
				if(!(coldata=realloc(coldata,sizeof(int *)*ldsize))) abt(MMsg);
				if(!(ncols=realloc(ncols,sizeof(int *)*ldsize))) abt(MMsg);
			}
			linedata[rec]=data+blkptr;
			coldata[rec]=cols+blkptr;
			ncols[rec]=0;
		}
		while(*p) {
			while(isspace((int)*p)) p++;
			if(*p) {
				if(nc>=maxcol) {
					i=maxcol;
					maxcol*=2;
					if(!(sum=realloc(sum,sizeof(double)*maxcol))) abt(MMsg);
					if(!(sum2=realloc(sum2,sizeof(double)*maxcol))) abt(MMsg);
					if(!(nn=realloc(nn,sizeof(int)*maxcol))) abt(MMsg);
					if(!(min=realloc(min,sizeof(double)*maxcol))) abt(MMsg);
					if(!(max=realloc(max,sizeof(double)*maxcol))) abt(MMsg);
					for(;i<maxcol;i++) {
						sum[i]=sum2[i]=0.0;
						nn[i]=0;
						min[i]=DBL_MAX;
						max[i]=-DBL_MAX;
					}
				}
				x=strtod(p,&p1);
				if(p!=p1) {
					if(!logopt) {
						sum[nc]+=x;
						sum2[nc]+=x*x;
					} else {
						if(nn[nc]>0.0) {
 							sum[nc]=addlog(sum[nc],x);
							sum2[nc]=addlog(sum[nc],2.0*x);
						} else {
							sum[nc]=x;
							sum2[nc]=2.0*x;
						}
					}
					if(x<min[nc]) min[nc]=x;
					if(x>max[nc]) max[nc]=x;
					nn[nc]++;
					if(ppoint>=0.0) add_data(x,nc,1);
				} else {
					if(ppoint>=0.0) add_data(0.0,nc,0);
				}
				nc++;
				p=p1;
				while(*p && !isspace((int)*p)) p++;
			}
		}
		if(ppoint>=0.0) ncols[rec]=nc;
		rec++;
		if(nc>ncol) ncol=nc;
	}
	return;
}

static int cmp_doubles(const void *s1,const void *s2)
{
	double x1,x2;

	x1=*((const double *)s1);
	x2=*((const double *)s2);
	if(x1>x2) return 1;
	if(x2>x1) return -1;
	return 0;
}

int main(int argc, char *argv[])
{
	char *p;
	FILE *fptr;
	int i,j,k,c,start=0,logopt=0;
	double ppoint=-1.0,x,v,*tx=0,z;

	while((c=getopt(argc,argv,"lr:p:"))!=-1) switch(c) {
	 case 'l':
		logopt=1;
		break;
	 case 'r':
		start=(int)strtol(optarg,&p,10);
		if(p==optarg) {
			(void)fprintf(stderr,"Invalid starting row number '%s' given to -r option\n",optarg);
			exit(EXIT_FAILURE);
		}
		break;
	 case 'p':
		ppoint=strtod(optarg,&p);
		if(p==optarg || ppoint<0.0 || ppoint>1.0) {
			(void)fprintf(stderr,"Invalid %% point '%s' given to -p option\n",optarg);
			exit(EXIT_FAILURE);
		}
		break;
	}
	if(optind>=argc) read_cols(stdin,logopt,start,ppoint);
	else for(i=optind;i<argc;i++) {
		if(!(fptr=fopen(argv[i],"r"))) {
			(void)fprintf(stderr,"Couldn't open file '%s' for reading: ",argv[i]);
			perror(0);
			exit(EXIT_FAILURE);
		}
		read_cols(fptr,logopt,start,ppoint);
		(void)fclose(fptr);
	}
	for(i=0;i<ncol;i++) if(!nn[i]) min[i]=max[i]=0.0;
	if(ppoint>=0 && rec) {
		if(!(tx=malloc(sizeof(double)*rec))) abt(MMsg);
	}
	for(i=0;i<ncol;i++) {
		if(!logopt) {
			x=nn[i]>0.0?sum[i]/(double)nn[i]:0.0;
			v=nn[i]>1.0?(sum2[i]-sum[i]*x)/(double)(nn[i]-1):0.0;
			(void)printf("Col %d - Sum = %g Mean = %g SD = %g n = %d Range %g -> %g\n",i+1,sum[i],x,sqrt(v),nn[i],min[i],max[i]);
		} else {
			x=nn[i]>0.0?sum[i]-log((double)nn[i]):0.0;
			(void)printf("Col %d - Sum = %g Mean = %g n = %d Range %g -> %g\n",i+1,sum[i],x,nn[i],min[i],max[i]);
		}
		if(ppoint>=0.0 && nn[i]) {
			for(z=0.0,k=j=0;j<rec;j++) if(i<ncols[j] && coldata[j][i]) {
				tx[k]=linedata[j][i];
				if(x>0.0) {
					if(tx[k]<=0.0) z++;
				} else {
					if(tx[k]>=0.0) z++;
				}
				k++;
			}
			if(k!=nn[i]) abt("Internal error - mismatch in record count\n");
			z/=(double)k;
			qsort(tx,(size_t)k,sizeof(double),cmp_doubles);
			(void)printf("     Q1 = %g, Q2 = %g, Q3 = %g, %g%% Lim = %g -> %g, p = %g\n",tx[(int)(nn[i]*.25)],
					 tx[(int)(nn[i]*.5)],tx[(int)(nn[i]*.75)],100.0*ppoint,tx[(int)(nn[i]*ppoint)],tx[(int)(nn[i]*(1.0-ppoint))],z);
		}
	}
	return 0;
}
