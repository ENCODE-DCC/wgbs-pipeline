#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define NKIDS 12

int main(int argc,char *argv[])
{
	int i,j,k,l,n,fg[2][2]={{1,2},{3,4}};
	int kids[NKIDS]={0,0,0,1,1,1,2,2,2,3,3,3},s[NKIDS][2],g[2],n_alls=2;
	double z,z1,z2,z3;
	double p[6][4]={
		 {0.15,0.24,0.16,0.45},
		{0.17,0.12,0.18,.53},
		{0.05,0.38,0.39,0.18},
		{0.6,0.15,0.11,0.14},
		{0.21,0.36,0.36,0.07},
		{0.28,0.2,0.2,0.32}
	};
	double rf[6][4],ph[NKIDS][4];
	
	n=NKIDS;
	for(i=0;i<n;i++) s[i][0]=s[i][1]=0;
	/* Parental R-Functions */
	for(j=0;j<4;j++) {
		rf[4][j]=p[4][j];
		rf[5][j]=p[5][j];
	}
	/* Child R-Functions */
	for(i=0;i<n;i++) {
		j=kids[i];
		for(k=0;k<4;k++) ph[i][k]=p[j][k];
	}
	for(;;) {
		for(j=0;j<4;j++) for(k=0;k<4;k++) rf[j][k]=1.0;
		for(i=0;i<n;i++) {
			g[0]=fg[0][s[i][0]];
			g[1]=fg[1][s[i][1]];
			j=(s[i][0]<<1)|s[i][1];
			for(k=0;k<4;k++) rf[j][k]*=ph[i][k];
			printf("%d%d ",g[0],g[1]);
		}
		z=0.0;
		for(i=0;i<n_alls;i++) {
			for(j=0;j<n_alls;j++) {
				z1=rf[4][i*n_alls+j]; /* RF(1,2) */
				for(k=0;k<n_alls;k++) {
					z2=z1*rf[0][i*n_alls+k]* /* RF(1,3) */
					  rf[2][j*n_alls+k]; /* RF(2,3) */
					for(l=0;l<n_alls;l++) {
						z3=z2*rf[1][i*n_alls+l]* /* RF(1,4) */
						  rf[3][j*n_alls+l]* /* RF(2,4) */
						  rf[5][k*n_alls+l]; /* RF(3,4) */
						z+=z3;
/*						printf("%d %d %d %d %g %g %g %g\n",i,j,k,l,z1,z2,z3,z); */
					}
				}
			}
		}
		printf("%g\n",log(z));
/*		for(j=0;j<6;j++) {
			for(k=0;k<4;k++) {
				printf("%g ",rf[j][k]);
			}
			printf("\n");
		}
		break; */
		for(i=n-1;i>0;i--) {
			s[i][1]^=1;
			if(s[i][1]) break;
			s[i][0]^=1;
			if(s[i][0]) break;
		}
		if(i<=0) break;
	}
	return 0;
}
