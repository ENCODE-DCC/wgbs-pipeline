#include <config.h>
#include <stdlib.h>
#ifdef USE_DMALLOC
#include <dmalloc.h>
#endif
#include <math.h>
#include <stdio.h>
#include <float.h>

#include "ranlib.h"
#include "utils.h"
#include "loki.h"
#include "loki_peel.h"
#include "get_par_probs.h"
#include "loki_trait_simple_peel.h"
#ifndef DBL_MAX
#define DBL_MAX MAXDOUBLE
#endif

/* Similar to loki_simple_sample, but for a trait locus */
double loki_trait_simple_sample(const struct Simple_Element *element,const int locus,const int s_flag,double **freq,struct R_Func *rf,trait_pen_func *trait_pen,struct loki *loki)
{
	int ids,idd,is1,id1,i,j,k,m,n,n_off,*off,kid,*ix,link;
	int i1,i2,j3,j2,k1,k2,l1,l2,jj,n_all,n_idx;
	int ix1[]={0,3,12,15};
	int ix3[]={5,6,9,10};
	int par_type[]={0,1,1,0,2,3,3,2,2,3,3,2,0,1,1,0};
	double *tp,p1,p2,z,z1,prob=0.0,pp[4],*tmp,*tmp1,*tpp,*qval,*pval,*mval,*peel_famval;
	struct peel_mem *work;
	struct Id_Record *id_array;
	struct Locus *loc;
	
	ids=element->sire-1;
	idd=element->dam-1;
	work=&loki->peel->workspace;
	off=element->off;
	n_off=element->n_off;
	loc=&loki->models->tlocus[-1-locus];
	n_all=loc->n_alleles;
	n_idx=n_all*n_all;
	peel_famval=work->s2;
	qval=peel_famval+16; 
	pval=qval+n_idx;
	mval=pval+n_idx;
	id_array=loki->pedigree->id_array;
#ifdef TRACE_PEEL
	if(CHK_PEEL(TRACE_LEVEL_1)) (void)printf("In %s(%p,%d,%d)\n",__func__,(void *)element,locus,s_flag);
	if(CHK_PEEL(TRACE_LEVEL_2)) {
		if(family_id) {
			print_orig_family(stdout,off[0]+1,0);
			(void)fputc(' ',stdout);
		}
		print_orig_id1(stdout,ids+1);
		(void)fputc(',',stdout);
		print_orig_id1(stdout,idd+1);
		(void)fputc(' ',stdout);
		for(i=0;i<n_off;i++) {
			(void)fputc(i?',':'(',stdout);
			print_orig_id1(stdout,off[i]);
		}
		(void)fputs(")\n",stdout);
	}
#endif
	if(id_array[ids].flag&SAMPLED_MAT) {
#ifdef DEBUG
		if(!(id_array[ids].flag&SAMPLED_PAT)) ABT_FUNC("Internal error: 1 paternal allele sampled\n");
#endif
		tmp=pval;
		for(j=0;j<4;j++) *(tmp++)=0.0;
		j=id_array[ids].allele[X_MAT]-1;
		k=id_array[ids].allele[X_PAT]-1;
		is1=(k<<1)|j;
		pval[is1]=1.0;
	} else {
		is1=-1;
		p1=get_trait_par_probs(pval,ids,locus,trait_pen,freq,rf,loki);
		if(p1<=0.0)	{
			if(!(s_flag&(1|OP_SAMPLING))) return -DBL_MAX;
			ABT_FUNC("Zero probability in sampling operation\n");
		}
		prob+=log(p1);
	}
	if(id_array[idd].flag&SAMPLED_MAT) {
#ifdef DEBUG
		if((!id_array[idd].flag&SAMPLED_PAT)) ABT_FUNC("Internal error: 1 maternal allele sampled\n");
#endif
		tmp=mval;
		for(j=0;j<4;j++) *(tmp++)=0.0;
		j=id_array[idd].allele[X_MAT]-1;
		k=id_array[idd].allele[X_PAT]-1;
		id1=(k<<1)|j;
		mval[id1]=1.0;
	} else {
		id1=-1;
		p1=get_trait_par_probs(mval,idd,locus,trait_pen,freq,rf,loki);
		if(p1<=0.0)	{
			if(!(s_flag&(1|OP_SAMPLING))) return -DBL_MAX;
			ABT_FUNC("Zero probability in sampling operation\n");
		}
		prob+=log(p1);
	}
	link=loc->link_group;
	if(is1>=0 || id1>=0) {
		tmp=peel_famval;
		for(j=0;j<16;j++) *(tmp++)=0.0;
		if(is1>=0) {
			if(id1>=0) {
				peel_famval[(is1<<2)|id1]=mval[id1]*pval[is1];
			} else {
				tmp=peel_famval+(is1<<2);
				tmp1=mval;
				for(m=0;m<4;m++) *(tmp++)=*(tmp1++);
			}
		} else {
			tmp1=pval;
			for(m=id1;m<16;m+=4) peel_famval[m]=*(tmp1++);
		}
	} else {
		tmp=peel_famval;
		for(n=0;n<4;n++) {
			p1=pval[n];
			tmp1=mval;
			for(m=0;m<4;m++) *(tmp++)=*(tmp1++)*p1;
		}
	}
	for(m=0;m<n_off;m++) {
		kid=off[m]-1;
		if(id_array[kid].flag&SAMPLED_MAT) {
			tmp=qval;
			for(k=0;k<4;k++) *(tmp++)=0.0;
			j=id_array[kid].allele[X_MAT]-1;
			k=id_array[kid].allele[X_PAT]-1;
			qval[(k<<1)|j]=1.0;
		} else {
			tmp=qval;
			if((k=id_array[kid].rfp)>=0) { /* Insert Previously computed R_Func */
				tmp1=rf[k].p;
				for(j=0;j<n_idx;j++) *(tmp++)=*(tmp1++);
			} else for(j=0;j<n_idx;j++) *(tmp++)=1.0;
			if(id_array[kid].res[0]) trait_pen(qval,kid,loc,loki);
			tmp=qval;
			for(p1=0.0,j=0;j<4;j++) p1+=*(tmp++);
			if(p1<=0.0)	{
				if(!(s_flag&(1|OP_SAMPLING))) return -DBL_MAX;
				ABT_FUNC("Zero probability in sampling operation\n");
			}
			z=1.0/p1;
			tmp=qval;
			for(j=0;j<4;j++) *(tmp++)*=z;
			prob+=log(p1);
		}
		/* First do double homozygote configs */
		tmp=qval;
		ix=ix1;
		for(i=0;i<4;i++) peel_famval[*(ix++)]*=(*tmp++);
		if(link<0) {
			/* Then do pat_hom / mat_het configs */
			z1=.5*(qval[0]+qval[1]);
			peel_famval[1]*=z1;
			peel_famval[2]*=z1;
			z1=.5*(qval[2]+qval[3]);
			peel_famval[13]*=z1;
			peel_famval[14]*=z1;
			/* Then do pat_het / mat_hom configs */
			z1=.5*(qval[0]+qval[2]);
			peel_famval[4]*=z1;
			peel_famval[8]*=z1;
			z1=.5*(qval[1]+qval[3]);
			peel_famval[7]*=z1;
			peel_famval[11]*=z1;
			/* Then do double het configs */
			ix=ix3;
			for(i=0;i<4;i++) peel_famval[*(ix++)]*=.25;
		} else {
			/* Then do pat_hom / mat_het configs */
			tp=id_array[kid].tp;
			tpp=id_array[kid].tpp[X_MAT];
			peel_famval[1]*=(tpp[X_MAT]*qval[1]+tpp[X_PAT]*qval[0]);
			peel_famval[2]*=(tpp[X_MAT]*qval[0]+tpp[X_PAT]*qval[1]);
			peel_famval[13]*=(tpp[X_MAT]*qval[3]+tpp[X_PAT]*qval[2]);
			peel_famval[14]*=(tpp[X_MAT]*qval[2]+tpp[X_PAT]*qval[3]);
			/* Then do pat_het / mat_hom configs */
			tpp=id_array[kid].tpp[X_PAT];
			peel_famval[4]*=(tpp[X_MAT]*qval[2]+tpp[X_PAT]*qval[0]);
			peel_famval[8]*=(tpp[X_MAT]*qval[0]+tpp[X_PAT]*qval[2]);
			peel_famval[7]*=(tpp[X_MAT]*qval[3]+tpp[X_PAT]*qval[1]);
			peel_famval[11]*=(tpp[X_MAT]*qval[1]+tpp[X_PAT]*qval[3]);
			/* Then do double het configs */
			peel_famval[5]*=(tp[X_MM_PM]*qval[3]+tp[X_MP_PM]*qval[2]+tp[X_MM_PP]*qval[1]+tp[X_MP_PP]*qval[0]);
			peel_famval[6]*=(tp[X_MM_PM]*qval[2]+tp[X_MP_PM]*qval[3]+tp[X_MM_PP]*qval[0]+tp[X_MP_PP]*qval[1]);
			peel_famval[9]*=(tp[X_MM_PM]*qval[1]+tp[X_MP_PM]*qval[0]+tp[X_MM_PP]*qval[3]+tp[X_MP_PP]*qval[2]);
			peel_famval[10]*=(tp[X_MM_PM]*qval[0]+tp[X_MP_PM]*qval[1]+tp[X_MM_PP]*qval[2]+tp[X_MP_PP]*qval[3]);
		}
	}
	p1=0.0;
	tmp=peel_famval;
	for(n=0;n<16;n++) p1+=*(tmp++);
	if(p1<=0.0) {
		if(!(s_flag&(1|OP_SAMPLING))) return -DBL_MAX;
		ABT_FUNC("Zero probability in sampling operation\n");
	}
	prob+=log(p1);
	do {
		z=ranf()*p1;
		p2=0.0;
		tmp=peel_famval;
		for(i=0;i<16;i++) {
			p1=*(tmp++);
			if(p1>0.0) {
				p2+=p1;
				if(p2>=z) {
					id_array[idd].allele[X_MAT]=(i&1)?2:1;
					id_array[idd].allele[X_PAT]=(i&2)?2:1;
					id_array[ids].allele[X_MAT]=(i&4)?2:1;
					id_array[ids].allele[X_PAT]=(i&8)?2:1;
					id_array[ids].flag|=(SAMPLED_MAT|SAMPLED_PAT);
					id_array[idd].flag|=(SAMPLED_MAT|SAMPLED_PAT);
					break;
				}
			}
		}
	} while(i==16);
	jj=par_type[i];
	if(!jj) { /* Double Homozygotic parents */
		i2=id_array[ids].allele[X_MAT];
		k2=id_array[idd].allele[X_MAT];
		for(m=0;m<n_off;m++)	{
			kid=off[m]-1;
			if(!(id_array[kid].flag&SAMPLED_MAT)) {
				id_array[kid].allele[X_MAT]=k2;
				id_array[kid].allele[X_PAT]=i2;
 				id_array[kid].flag|=(SAMPLED_MAT|SAMPLED_PAT);
			} 
		}
	} else {
		i2=id_array[ids].allele[X_MAT];
		j2=id_array[ids].allele[X_PAT];
		k2=id_array[idd].allele[X_MAT];
		l2=id_array[idd].allele[X_PAT];
		i=(i2-1)<<1;
		j=(j2-1)<<1;
		i1=i+k2-1;
		j3=j+k2-1;
		k1=i+l2-1;
		l1=j+l2-1;
		for(m=0;m<n_off;m++)	{
			kid=off[m]-1;
			if(id_array[kid].flag&SAMPLED_MAT) continue;
			tmp=qval;
			if((k=id_array[kid].rfp)>=0) { /* Insert Previously computed R_Func */
				tmp1=rf[k].p;
				for(j=0;j<4;j++) *(tmp++)=*(tmp1++);
			} else for(j=0;j<4;j++) *(tmp++)=1.0;
			if(id_array[kid].res[0]) trait_pen(qval,kid,loc,loki);
			if(jj==1) {
				tpp=id_array[kid].tpp[X_MAT];
				pp[X_MAT]=tpp[X_MAT]*qval[i1];
				pp[X_PAT]=tpp[X_PAT]*qval[k1];
				p1=pp[X_MAT]+pp[X_PAT];
#ifdef DEBUG
				if(p1<=0.0) ABT_FUNC("Internal error - no offspring combination possible\n");
#endif
				z=ranf()*p1;
				id_array[kid].allele[X_MAT]=(z<=pp[X_MAT])?k2:l2;
				id_array[kid].allele[X_PAT]=i2;
			} else if(jj==2) {
				tpp=id_array[kid].tpp[X_PAT];
				pp[X_MAT]=tpp[X_MAT]*qval[i1];
				pp[X_PAT]=tpp[X_PAT]*qval[j3];
				p1=pp[X_MAT]+pp[X_PAT];
#ifdef DEBUG
				if(p1<=0.0) ABT_FUNC("Internal error - no offspring combination possible\n");
#endif
 				z=ranf()*p1;
				id_array[kid].allele[X_PAT]=(z<=pp[X_MAT])?i2:j2;
				id_array[kid].allele[X_MAT]=k2;
			} else {
				/* transmission probs */
				tp=id_array[kid].tp;
				p1=(pp[X_MM_PM]=tp[X_MM_PM]*qval[i1]);
				p1+=(pp[X_MM_PP]=tp[X_MM_PP]*qval[j3]);
				p1+=(pp[X_MP_PM]=tp[X_MP_PM]*qval[k1]);
				p1+=(pp[X_MP_PP]=tp[X_MP_PP]*qval[l1]);
#ifdef DEBUG
				if(p1<=0.0) ABT_FUNC("Internal error - no offspring combination possible\n");
#endif
				z=safe_ranf()*p1;
				p2=0.0;
				for(n=0;n<4;n++) {
					if(pp[n]>0.0) {
						p2+=pp[n];
						if(z<=p2) break;
					}
				}
				switch(n) {
				 case X_MM_PM:
					id_array[kid].allele[X_MAT]=k2;
					id_array[kid].allele[X_PAT]=i2;
					break;
				 case X_MM_PP:
					id_array[kid].allele[X_MAT]=k2;
					id_array[kid].allele[X_PAT]=j2;
					break;
				 case X_MP_PM:
					id_array[kid].allele[X_MAT]=l2;
					id_array[kid].allele[X_PAT]=i2;
					break;
				 case X_MP_PP:
					id_array[kid].allele[X_MAT]=l2;
					id_array[kid].allele[X_PAT]=j2;
					break;
#ifdef DEBUG
				 default:
					ABT_FUNC("Internal error - illegal sample\n");
#endif
				}
			}
			id_array[kid].flag|=(SAMPLED_MAT|SAMPLED_PAT);
		}
	}
	return prob;
}

