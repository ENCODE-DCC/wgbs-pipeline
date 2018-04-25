#include <config.h>
#include <stdlib.h>
#include <stdio.h>

#include "libhdr.h"
#include "ranlib.h"
#include "utils.h"

#define N              624                 
#define M              397                 
#define umask          0x80000000U
#define lmask          0x7FFFFFFFU
#define mixBits(u,v)   (((u)&umask)|((v)&lmask))
#define twist(u,v)     ((mixBits(u,v)>>1)^((v)&1U?aa:0U))
static unsigned int rgen[][3]={
	{0x9908B0DFU,0x9D2C5680U,0xEFC60000U}, /* Original MT19937 param */
	 {0x8c400000U,0xbb56ef00U,0xddd58000U}, /* Set of independent 2^19937 RNG parameters */
	 {0xb0b70001U,0x4e753a80U,0xf7d48000U},
	 {0x92db0002U,0xab55b780U,0x76538000U},
	 {0x852b0003U,0x9ddc6d80U,0xef558000U},
	 {0xd6db0004U,0x26eef680U,0xeec58000U},
	 {0xc36d0005U,0xdd4eec80U,0xefd38000U},
	 {0x938c0006U,0xcd757780U,0xee778000U},
	 {0xd6810007U,0xf2976c80U,0x6fcb8000U},
	 {0x872d0008U,0xbb54ee80U,0xddd58000U},
	 {0xf8950009U,0x65bd7380U,0xef4d8000U},
	 {0x8451000aU,0x38d4ed00U,0xddd78000U},
	 {0x97f4000bU,0x99e47580U,0xef458000U},
	 {0xd6d4000cU,0xbbd5be80U,0xf5d48000U},
	 {0xbea1000dU,0x39913a80U,0xf7d58000U},
	 {0xbcd0000eU,0x9dd4b500U,0xf7d58000U},
	 {0xae1d000fU,0x9d627680U,0xeff78000U},
	 {0x94df0010U,0x39573b80U,0xf7958000U},
	 {0xaf580011U,0xd9b77b80U,0xf7558000U},
	 {0xbaf2003cU,0xf735bf80U,0xba958000U},
	 {0xcd870013U,0x39b53a80U,0xf7d48000U},
	 {0x949f0014U,0x9b76bb80U,0xf7d18000U},
	 {0xaf1e0015U,0x9d727680U,0xeff78000U},
	 {0xf6100016U,0xec75bb80U,0xbb578000U},
	 {0x83d4003aU,0x6e557b80U,0xf7758000U},
	 {0xfe690018U,0xeed5bb80U,0x77558000U},
	 {0x88c00019U,0x6ea57b00U,0xf7f18000U},
	 {0xd751001aU,0x9dd47380U,0xef7d8000U},
	 {0xff4f001bU,0x9d8d7780U,0xed7d8000U},
	 {0xb8e0001cU,0xd674ea80U,0xefd48000U},
	 {0xf947001dU,0x2a753680U,0xf7d28000U},
	 {0x9b9d001eU,0xbb25bf80U,0xf6b58000U},
	 {0xad97001fU,0xb914ee80U,0xddd58000U},
	 {0x9e9d0020U,0xdf96ee80U,0xedd58000U},
	 {0xac830021U,0x93b4bb80U,0xf7568000U},
	 {0x88250022U,0xf3376f80U,0xdd558000U},
	 {0xed200023U,0xb7b57680U,0xddd58000U},
	 {0xf2300024U,0xb956ee00U,0xedd58000U},
	 {0xe08b0025U,0xab7d6d80U,0xdb578000U},
	 {0xde240026U,0xef55bb80U,0x77958000U},
	 {0xb0660027U,0x5b553b00U,0xf7d58000U},
	 {0xc4e00028U,0xb954ee80U,0xddd58000U},
	 {0x9ffc0029U,0x3b553b00U,0xf7d78000U},
	 {0x9fe5002aU,0x8f567600U,0xedd58000U},
	 {0xe52c002bU,0xf3d56d80U,0x6f5d8000U},
	 {0xd9b4002cU,0xbcb3b780U,0xf6578000U},
	 {0x958f002dU,0xf6f77e80U,0xdad58000U},
	 {0x8e91002fU,0x9d757680U,0xeef78000U},
	 {0xbc890031U,0xf5977680U,0xddf58000U},
	 {0xd8650032U,0xf5957780U,0xed558000U},
	 {0x87210033U,0xef55bb00U,0x77d58000U},
	 {0x9e000034U,0xcf557780U,0xed758000U},
	 {0xa1e60035U,0xf677bb80U,0xbb558000U},
	 {0xbe030036U,0x9b72bb80U,0xf7d78000U},
	 {0xca290037U,0x75b77780U,0xed558000U},
	 {0xd47d0038U,0x65ad7380U,0xef558000U},
	 {0x90110039U,0x9994bb80U,0xf7558000U},
	 {0,0,0}
};

static unsigned int state[N+1],set;     
static unsigned int *next,maskB,maskC,aa;
static int left=-1,shift0,gen_idx=0;
  
int set_mt_idx(int idx)
{
	int i,err=-1;
	
	if(idx>=0) {
		i=0;
		while(rgen[i][0]) i++;
		if(i>idx && rgen[idx][1]) {
			shift0=idx?12:11;
			aa=rgen[idx][0];
			maskB=rgen[idx][1];
			maskC=rgen[idx][2];
			if(set && idx!=gen_idx) err=1;
			else err=0;
			gen_idx=idx;
		}
	}
	set=1;
	return err;
}

void sgenrand(unsigned int seed)
{ 
	int j;

	if(!set) set_mt_idx(0);
	state[0]=seed&0xffffffffU;
	for(j=1;j<N;j++) {
		state[j]=(1812433253U*(state[j-1]^(state[j-1]>>30))+j);
		state[j]&=0xffffffffU;
	}
	left=0;
}

static void next_state(void)
{
	unsigned int *p=state;
	int j;

	left=N-1;
	next=state;
	for(j=N-M+1;--j;p++) *p=p[M]^twist(p[0],p[1]);
	for(j=M;--j;p++) *p=p[M-N]^twist(p[0],p[1]);
	*p=p[M-N]^twist(p[0],state[0]);
}

unsigned int genint(void)
{
	unsigned int y;
	
 	if(--left<0) next_state();
	y=*next++;
	y^=(y>>shift0);
	y^=(y<<7)&maskB;
	y^=(y<<15)&maskC;
	return (y^(y >> 18));
}

double genrand(void)
{
	unsigned int y;
	
 	if(--left<0) next_state();
	y=*next++;
	y^=(y>>shift0);
	y^=(y<<7)&maskB;
	y^=(y<<15)&maskC;
	return (y^(y >> 18))*(1.0/4294967295.0);
}

double safe_genrand(void)
{
	unsigned int y;
	
 	if(--left<0) next_state();
	y=*next++;
	y^=(y>>shift0);
	y^=(y<<7)&maskB;
	y^=(y<<15)&maskC;
	return ((double)(y^(y >> 18))+.5)*(1.0/4294967296.0);

}

int getseed(const char *fname)
{
	int flag=0,i,j,err=0;
	unsigned int sd;
	char buf[256],*fname1;
	FILE *fptr;
	
	if(!fname) {
		message(WARN_MSG,"getseed(): NULL filename\n");
		return -1;
	}
	fname1=add_file_dir(fname);
	if(!fname1) {
		message(WARN_MSG,"getseed(): Couldn't make filename\n");
		return -1;
	}
	fptr=fopen(fname1,"r");
	if(fptr)	{
		if(fgets(buf,256,fptr))	{
			free(fname1);
			flag=sscanf(buf,"mt19937b_seed = %u\n",&sd);
			if(flag!=1)	{
				if(sscanf(buf,"mt19937b_idx = %d\n",&left)==1) {
					if(fscanf(fptr,"mt19937b_gen = %d\n",&i)!=1) i=0;
					j=set_mt_idx(i);
					if(j==-1) {
						message(WARN_MSG,"getseed(): Bad random number generator (%d)\n",i);
						(void)fclose(fptr);
						set_mt_idx(0);
						return -1;
					}
					if(j==1) message(DEBUG_MSG,"getseed(): Changing random number generator to %d\n",i);
					if(left>=0 && left<=N) {
						for(i=0;i<N;i++) if(fscanf(fptr,"%u",state+i)!=1) break;
						next=state+N-left;
						if(i==N)	{
							message(DEBUG_MSG,"getseed(): Using stored state for mt19937b generator\n");
							return 0;
						}
					}
				}
			}
		}
		(void)fclose(fptr);
		if(flag!=1) {
			message(WARN_MSG,"getseed(): Couldn't read seed from file %s\n",fname1);
			err=-1;
		}
	} else {
		message(WARN_MSG,"getseed(): Couldn't open seedfile %s for reading\n",fname1);
		err=-1;
	}
	if(!err) {
		message(DEBUG_MSG,"getseed(): Using seed from file (%u)\n",sd);
		sgenrand(sd);
	}
	free(fname1);
	return err;
}

int dumpseed(FILE *fptr,const int flag)
{
	int i,err=0;
	unsigned int sd;
	
	if(flag)	{
		err=fprintf(fptr,"mt19937b_idx = %d\n",left);
		if(err>0 && gen_idx) err=fprintf(fptr,"mt19937b_gen = %d\n",gen_idx);
		if(err>0) for(i=0;i<N;i++)	{
			err=fprintf(fptr,"%u%c",state[i],(i+1)%16?' ':'\n');
			if(err<0) break;
		}
	} else {
		sd=genint();
		err=fprintf(fptr,"mt19937b_seed = %u\n",sd);
	}
	return err;
}

int bindumpseed(FILE *fptr)
{
	int i,err=0;
	
	if(gen_idx) {
		if(fprintf(fptr,"%x,%x,%x\n",N,left,gen_idx)<0) err=1;
	} else if(fprintf(fptr,"%x,%x\n",N,left)<0) err=1;
	for(i=0;!err && i<N;i++) if(fprintf(fptr,"%x\n",state[i])<0) err=1;
	return err;
}

int binreadseed(FILE *fptr,char *s)
{
	int err=0,k1;
	unsigned int i,j,k=0;
	
	if(sscanf(s,"%x,%x,%x\n",&i,&j,&k)!=3) {
		if(sscanf(s,"%x,%x\n",&i,&j)!=2) err=1;
	}
	k1=set_mt_idx((int)k);
	if(k1==-1) {
		set_mt_idx(0);
		err=1;
	}
	if(i!=N) err=1;
	if(!err) for(i=0;!err && i<N;i++) if(fscanf(fptr,"%x\n",state+i)!=1) err=1;
	left=(int)j;
	if(!err) next=state+N-left;
	return err;
}

int writeseed(const char *fname,const int flag)
{
	int err;
	FILE *fptr;
	char *fname1;
	
	if(!fname) {
		message(WARN_MSG,"writeseed(): NULL filename\n");
		return -1;
	}
	fname1=add_file_dir(fname);
	if(fname1) {
		(void)mkbackup(fname1,1);
		fptr=fopen(fname1,"w");
		free(fname1);
		if(fptr) {
			err=dumpseed(fptr,flag);
			(void)fclose(fptr);
		} else {
			message(WARN_MSG,"writeseed(): Couldn't open seedfile for writing\n");
			err=-1;
		}
	} else {
		message(WARN_MSG,"writeseed(): Couldn't create seedfile name\n");
		err=-1;
	}
	return err;
}
