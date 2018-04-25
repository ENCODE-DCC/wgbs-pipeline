/****************************************************************************
 *                                                                          *
 *     Loki - Programs for genetic analysis of complex traits using MCMC    *
 *                                                                          *
 *                    Simon Heath - CNG, France                             *
 *                                                                          *
 *                         August 2002                                      *
 *                                                                          *
 * fenris_simple_peel.c:                                                    *
 *                                                                          *
 * Copyright (C) Simon C. Heath 2002                                        *
 * This is free software.  You can distribute it and/or modify it           *
 * under the terms of the Modified BSD license, see the file COPYING        *
 *                                                                          *
 ****************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#ifdef USE_DMALLOC
#include <dmalloc.h>
#endif
#include <math.h>

#include "utils.h"
#include "loki.h"
#include "loki_peel.h"
#include "fenris_peel.h"
#include "fenris.h"

#include "fenris_simple_macros.h"

static double *pval,*mval,*prf,*mrf,*val,*val1,**krf,**kpost,**simple_rf;
double **kval;
int *fc_state;

void fenris_peel_alloc(struct Peelseq_Head *pp)
{
	int i,j,k,max_off=0,max_inv=0,max_peel=0,n_piv=0,cf=0;
	struct Peelseq_Head *p;
	struct Fenris_Simple_Element *element;
	struct Complex_Element *c_element;
	  
	for(i=0;i<n_comp;i++) {
		p=pp+i;
		while(p->type) {
			switch(p->type) {
			 case FENRIS_PEEL_SIMPLE:
				element=p->ptr.fsimple;
				j=element->n_off;
				if(j>max_off) max_off=j;
				if(element->pivot) n_piv++;
				p=&element->next;
				break;
			 case PEEL_COMPLEX:
				c_element=p->ptr.complex;
				j=c_element->n_peel;
				k=c_element->n_involved;
				if(j>max_peel) max_peel=j;
				if(k>max_inv) max_inv=k;
				cf=1;
				p=&c_element->next;
				break;
			 default:
				ABT_FUNC("Peel type not handled\n");
			}
		}
	}
	/* Find maximum number of alleles required */
	for(i=j=0;i<n_markers;i++) {
		for(k=0;k<n_comp-singleton_flag;k++) {
			if(marker[i].n_all1[k]>j) j=marker[i].n_all1[k];
		}
	}
	/* And calulate maximum number of genotypes */
	i=j*(j+1)/2;
	/* Allocate space for simple R-Functions (those coming from simple peel ops) */
	if(n_piv) {
		if(!(simple_rf=malloc(sizeof(void *)*n_piv))) ABT_FUNC(MMsg);
		if(!(simple_rf[0]=malloc(sizeof(double)*n_piv*i))) ABT_FUNC(MMsg);
		for(k=1;k<n_piv;k++) simple_rf[k]=simple_rf[k-1]+i;
	}
	if(max_inv) {
		if(!(fc_state=malloc(sizeof(int)*max_inv*6))) ABT_FUNC(MMsg);
	}
	if(max_peel>max_off) max_off=max_peel;
	if(max_off) {
		if(!(kval=malloc(sizeof(void *)*max_off*3))) ABT_FUNC(MMsg);
		krf=kval+max_off;
		kpost=krf+max_off;
	}
	k=(max_off*3+6)*i;
	if(!(pval=malloc(sizeof(double)*k))) ABT_FUNC(MMsg); /* Storage for paternal function */
	mval=pval+i; /* maternal function */
	val=mval+i;
	val1=val+i;
	prf=val1+i;
	mrf=prf+i;
	kval[0]=mrf+i; /* child functions */
	for(j=1;j<max_off*3;j++) kval[j]=kval[j-1]+i;
	if(cf) setup_fenris_complex_peel();
}

void fenris_peel_free(void)
{
	if(pval) {
		free(pval);
		pval=0;
	}
	if(kval) {
		free(kval);
		kval=0;
	}
	if(simple_rf) {
		if(simple_rf[0]) free(simple_rf[0]);
		free(simple_rf);
		simple_rf=0;
	}
	if(fc_state) {
		free(fc_state);
		fc_state=0;
	}
	free_fenris_complex_peel();
}

double get_par_probs(double *p,int id,int locus,fenris_pen_func *pen,struct pen_par *ppar,double **freq)
{
	int i,j,k,n_all,n_gen,comp;
	double *fq,*p1,*p2,z,pp;
	
	comp=id_array[id].comp;
	n_all=marker[locus].n_all1[comp];
	n_gen=n_all*(n_all+1)/2;
	if(!founder_flag[locus][id]) { /* Not a founder */
		k=id_array[id].rfp;
		p1=p;
		if(k>=0) {
			p2=simple_rf[k];
			for(i=0;i<n_gen;i++) *p1++=*p2++;
		} else for(i=0;i<n_gen;i++) *p1++=1.0;
	} else {
		k=id_array[id].group-1;
		fq=freq[k];
		k=id_array[id].rfp;
		p1=p;
		if(k>=0) {
			p2=simple_rf[k];
			for(i=0;i<n_all;i++) {
				z=fq[i];
				for(j=0;j<i;j++) *p1++=2.0*z*fq[j]*(*p2++);
				*p1++=z*z*(*p2++);
			}
		} else {
			for(i=0;i<n_all;i++) {
				z=fq[i];
				for(j=0;j<i;j++) *p1++=2.0*z*fq[j];
				*p1++=z*z;
			}
		}
	}
	pp=1.0;
	if(marker[locus].haplo[id]) {
		z=pen(p,id,locus,ppar);
		if(z<RESCALE_LIMIT) {
			pp=z;
			z=1.0/z;
			for(i=0;i<n_gen;i++) (*p++)*=z;
		}
	}
	return pp;
}
  
double fenris_simple_peel(struct Fenris_Simple_Element *element,int locus,double **freq,fenris_pen_func *pen,struct pen_par *ppar)
{
	int i,j,k,l,kid,*off,n_off,comp,n_all,n_gen,pivot,ids,idd;
	int i1,j1,g1,g2,g3,g4,kk,pflag=1,ix1,ix2,*rfp;
	struct Marker *mark;
	double *p,*p1,*pa,*p1a,*p2,z,z1,za,prob=0.0,*pp1,*pp2;

	ids=element->sire-1;
	idd=element->dam-1;
	comp=id_array[ids].comp;
	/* Get information about marker */
	mark=marker+locus;
	n_all=mark->n_all1[comp];
	n_gen=n_all*(n_all+1)/2;
	pivot=element->pivot-1;
#ifdef DEBUG
	if(pivot<0) ABT_FUNC("Zero pivot in peeling op\n");
#endif
	rfp=element->rf;
	id_array[ids].rfp=rfp[0];
	id_array[idd].rfp=rfp[1];
	/* Create functions for parents */
	if(ids!=pivot) {
		z=get_par_probs(pval,ids,locus,pen,ppar,freq);
		prob+=log(z);
	} else {
		pflag=2;
		p=pval;
		k=id_array[ids].rfp;
		if(k>=0) {
			p1=simple_rf[k];
			for(j=0;j<n_gen;j++) *p++=*p1++;
		} else for(j=0;j<n_gen;j++) *p++=1.0;
	}
	if(idd!=pivot) {
		z=get_par_probs(mval,idd,locus,pen,ppar,freq);
		prob+=log(z);
	} else {
		pflag=3;
		p=mval;
		k=id_array[idd].rfp;
		if(k>=0) {
			p1=simple_rf[k];
			for(j=0;j<n_gen;j++) *p++=*p1++;
		} else for(j=0;j<n_gen;j++) *p++=1.0;
	}
	/* Create functions for offspring to be peeled */
	off=element->off;
	n_off=element->n_off;
	for(i=0;i<n_off;i++) {
		kid=off[i];
		id_array[kid].rfp=rfp[2+i];
	}
	if(pflag==1) n_off--; /* If peeling to an offspring, don't get R-func for that offspring (always the last one) */
	for(i=0;i<n_off;i++) {
		kid=off[i];
		p=kval[i];
		k=id_array[kid].rfp;
		if(k>=0) {
			p1=simple_rf[k];
			for(j=0;j<n_gen;j++) *p++=*p1++;
		} else for(j=0;j<n_gen;j++) *p++=1.0;
		if(mark->haplo[kid]) {
			z=pen(kval[i],kid,locus,ppar);
			if(z<RESCALE_LIMIT) {
				prob+=log(z);
				z=1.0/z;
				p=kval[i];
				for(j=0;j<n_gen;j++) (*p++)*=z;
			}
		}
	}
	/* Loop through parental haplotypes */
	/* i,j - paternal alleles
	 * k,l - maternal alleles
	 * 
	 * The loop structures are complicated to reduce tests within the loops.
	 * This way we know if l<k=j<i etc. without testing
	 */
	if(pflag==1) { /* Peeling to child */
		p=val;
		for(j=0;j<n_gen;j++) *p++=0.0;
		p=pval;
		pa=mval;
		for(i1=i=0;i<n_all;i++) {
			i1+=i;
			for(j1=j=0;j<i;j++) {
				j1+=j;
				z=*p++;
				za=*pa++;
				p1=mval;
				p1a=pval;
				g1=i1;
				g3=j1;
				for(k=0;k<=j;k++,g1++,g3++) {
					g2=i1;
					g4=j1;
					for(l=0;l<k;l++,g2++,g4++) {
						z1=.25*(z*(*p1++)+za*(*p1a++));
						/* l<k<=j<i */
						for(kk=0;kk<n_off;kk++) {
							p2=kval[kk];
							z1*=.25*(p2[g1]+p2[g2]+p2[g3]+p2[g4]);
						}
						val[g1]+=z1;
						val[g2]+=z1;
						val[g3]+=z1;
						val[g4]+=z1;
					}
					/* l=k<=j<i */
					z1=.5*(z*(*p1++)+za*(*p1a++));
					for(kk=0;kk<n_off;kk++) {
						p2=kval[kk];
						z1*=.5*(p2[g1]+p2[g3]);
					}
					val[g1]+=z1;
					val[g3]+=z1;
				}
				g3=k*(k-1)/2+j;
				for(;k<i;k++,g1++) {
					g3+=k;
					g2=i1;
					g4=j1;
					for(l=0;l<=j;l++,g2++,g4++) {
						z1=.25*(z*(*p1++)+za*(*p1a++));
						/* l<=j<k<i */
						for(kk=0;kk<n_off;kk++) {
							p2=kval[kk];
							z1*=.25*(p2[g1]+p2[g2]+p2[g3]+p2[g4]);
						}
						val[g1]+=z1;
						val[g2]+=z1;
						val[g3]+=z1;
						val[g4]+=z1;
					}
					g4=l*(l-1)/2+j;
					for(;l<k;l++,g2++) {
						g4+=l;
						z1=.25*(z*(*p1++)+za*(*p1a++));
						/* j<l<k<i */
						for(kk=0;kk<n_off;kk++) {
							p2=kval[kk];
							z1*=.25*(p2[g1]+p2[g2]+p2[g3]+p2[g4]);
						}
						val[g1]+=z1;
						val[g2]+=z1;
						val[g3]+=z1;
						val[g4]+=z1;
					}
					/* j<l=k<i */
					z1=.5*(z*(*p1++)+za*(*p1a++));
					for(kk=0;kk<n_off;kk++) {
						p2=kval[kk];
						z1*=.5*(p2[g1]+p2[g3]);
					}
					val[g1]+=z1;
					val[g3]+=z1;
				}
				g3+=k;
				g2=i1;
				g4=j1;
				for(l=0;l<j;l++,g2++,g4++) {
					z1=.25*(z*(*p1++)+za*(*p1a++));
					/* l<j<k=i */
					for(kk=0;kk<n_off;kk++) {
						p2=kval[kk];
						z1*=.25*(p2[g1]+p2[g2]+p2[g3]+p2[g4]);
					}
					val[g1]+=z1;
					val[g2]+=z1;
					val[g3]+=z1;
					val[g4]+=z1;
				}
				/* l=j<k=i */
				z1=z*(*p1++);
				for(kk=0;kk<n_off;kk++) {
					p2=kval[kk];
					z1*=.25*(p2[g1]+p2[g4])+.5*p2[g2];
				}
				val[g1]+=.25*z1;
				val[g2]+=.5*z1;
				val[g4]+=.25*z1;
			}
			z=*p++;
			za=*pa++;
			p1=mval;
			p1a=pval;
			g1=i1;
			for(k=0;k<i;k++,g1++) {
				g2=i1;
				for(l=0;l<k;l++,g2++) {
					z1=.5*(z*(*p1++)+za*(*p1a++));
					/* l<k<j=i */
					for(kk=0;kk<n_off;kk++) {
						p2=kval[kk];
						z1*=.5*(p2[g1]+p2[g2]);
					}
					val[g1]+=z1;
					val[g2]+=z1;
				}
				z1=z*(*p1++)+za*(*p1a++);
				/* l=k<j=i */
				for(kk=0;kk<n_off;kk++) {
					p2=kval[kk];
					z1*=p2[g1];
				}
				val[g1]+=z1;
			}
			g2=i1;
			for(l=0;l<k;l++,g2++) {
				z1=.5*(z*(*p1++)+za*(*p1a++));
				/* l<k=j=i */
				for(kk=0;kk<n_off;kk++) {
					p2=kval[kk];
					z1*=.5*(p2[g1]+p2[g2]);
				}
				val[g1]+=z1;
				val[g2]+=z1;
			}
			/* l=k=j=i */
			z1=z*(*p1++);
			for(kk=0;kk<n_off;kk++) {
				p2=kval[kk];
				z1*=p2[g1];
			}
			val[g1]+=z1;
		}
		p=val;
		z=0.0;
		/* Get existing R-Functions on pivot if present */
		k=id_array[pivot].rfp;
		if(k>=0) { /* Existing R-Function on this individual */
			p1=simple_rf[k];
			for(i=0;i<n_gen;i++) z+=((*p++)*=(*p1++));
		} else {
			for(i=0;i<n_gen;i++) z+=*p++;
		}
		prob+=log(z);
		z=1.0/z;
		k=element->out_index;
		p=simple_rf[k];
		p1=val;
		for(i=0;i<n_gen;i++) *p++=z*(*p1++);
		id_array[pivot].rfp=k;
	} else { /* Peeling to parent */
		p1=val;
		for(i=0;i<n_gen;i++) *p1++=0.0;
		if(pflag==2) { /* Father */
			pp1=pval;
			pp2=mval;
		} else { /* Mother */
			pp1=mval;
			pp2=pval;
		}
		p=pp1;
		pa=pp2;
		for(i1=ix1=i=0;i<n_all;i++) {
			i1+=i;
			for(j1=j=0;j<i;j++,ix1++) {
				j1+=j;
				z=*p++;
				za=*pa++;
				p1=pp2;
				p1a=pp1;
				g1=i1;
				g3=j1;
				for(ix2=k=0;k<=j;k++,g1++,g3++) {
					g2=i1;
					g4=j1;
					for(l=0;l<k;l++,g2++,g4++,ix2++) {
						/* l<k<=j<i */
						z1=1.0;
						for(kk=0;kk<n_off;kk++) {
							p2=kval[kk];
							z1*=.25*(p2[g1]+p2[g2]+p2[g3]+p2[g4]);
						}
						val[ix1]+=z1*z*(*p1++);
						val[ix2]+=z1*za*(*p1a++);
					}
					/* l=k<=j<i */
					z1=1.0;
					for(kk=0;kk<n_off;kk++) {
						p2=kval[kk];
						z1*=.5*(p2[g1]+p2[g3]);
					}
					val[ix1]+=z1*z*(*p1++);
					val[ix2++]+=z1*za*(*p1a++);
				}
				g3=k*(k-1)/2+j;
				for(;k<i;k++,g1++) {
					g3+=k;
					g2=i1;
					g4=j1;
					for(l=0;l<=j;l++,g2++,g4++,ix2++) {
						z1=1.0;
						/* l<=j<k<i */
						for(kk=0;kk<n_off;kk++) {
							p2=kval[kk];
							z1*=.25*(p2[g1]+p2[g2]+p2[g3]+p2[g4]);
						}
						val[ix1]+=z1*z*(*p1++);
						val[ix2]+=z1*za*(*p1a++);
					}
					g4=l*(l-1)/2+j;
					for(;l<k;l++,g2++,ix2++) {
						g4+=l;
						z1=1.0;
						/* j<l<k<i */
						for(kk=0;kk<n_off;kk++) {
							p2=kval[kk];
							z1*=.25*(p2[g1]+p2[g2]+p2[g3]+p2[g4]);
						}
						val[ix1]+=z1*z*(*p1++);
						val[ix2]+=z1*za*(*p1a++);
					}
					/* j<l=k<i */
					z1=1.0;
					for(kk=0;kk<n_off;kk++) {
						p2=kval[kk];
						z1*=.5*(p2[g1]+p2[g3]);
					}
					val[ix1]+=z1*z*(*p1++);
					val[ix2++]+=z1*za*(*p1a++);
				}
				g3+=k;
				g2=i1;
				g4=j1;
				for(l=0;l<j;l++,g2++,g4++,ix2++) {
					z1=1.0;
					/* l<j<k=i */
					for(kk=0;kk<n_off;kk++) {
						p2=kval[kk];
						z1*=.25*(p2[g1]+p2[g2]+p2[g3]+p2[g4]);
					}
					val[ix1]+=z1*z*(*p1++);
					val[ix2]+=z1*za*(*p1a++);
				}
				/* l=j<k=i */
				z1=z*(*p1++);
				for(kk=0;kk<n_off;kk++) {
					p2=kval[kk];
					z1*=.25*(p2[g1]+p2[g4])+.5*p2[g2];
				}
				val[ix1]+=z1;
			}
			z=*p++;
			za=*pa++;
			p1=pp2;
			p1a=pp1;
			g1=i1;
			for(ix2=k=0;k<i;k++,g1++) {
				g2=i1;
				for(l=0;l<k;l++,g2++,ix2++) {
					z1=1.0;
					/* l<k<j=i */
					for(kk=0;kk<n_off;kk++) {
						p2=kval[kk];
						z1*=.5*(p2[g1]+p2[g2]);
					}
					val[ix1]+=z1*z*(*p1++);
					val[ix2]+=z1*za*(*p1a++);
				}
				z1=1.0;
				/* l=k<j=i */
				for(kk=0;kk<n_off;kk++) {
					p2=kval[kk];
					z1*=p2[g1];
				}
				val[ix1]+=z1*z*(*p1++);
				val[ix2++]+=z1*za*(*p1a++);
			}
			g2=i1;
			for(l=0;l<k;l++,g2++,ix2++) {
				z1=1.0;
				/* l<k=j=i */
				for(kk=0;kk<n_off;kk++) {
					p2=kval[kk];
					z1*=.5*(p2[g1]+p2[g2]);
				}
				val[ix1]+=z1*z*(*p1++);
				val[ix2]+=z1*za*(*p1a++);
			}
			/* l=k=j=i */
			z1=z*(*p1++);
			for(kk=0;kk<n_off;kk++) {
				p2=kval[kk];
				z1*=p2[g1];
			}
			val[ix1++]+=z1;
		}
		p=val;
		for(z=0.0,i=0;i<n_gen;i++) z+=*p++;
		prob+=log(z);
		z=1.0/z;
		k=element->out_index;
		p=simple_rf[k];
		p1=val;
		for(i=0;i<n_gen;i++) *p++=z*(*p1++);
		id_array[pivot].rfp=k;
	}
	return prob;
}

/* Distribute evidence back to all nuclear families.
 * Collect information on expected distribution of founder alleles */
double fenris_simple_distribute(struct Fenris_Simple_Element *element,int locus,double **freq,double **nfreq,fenris_pen_func *pen,struct pen_par *ppar)
{
	int i,j,k,l,kid,*off,n_off,comp,n_all,n_gen,pivot,ids,idd,n_off1,peel_start;
	int i1,j1,g1,g2,g3,g4,kk,ix1,ix2,sire_rfp,dam_rfp,parflag,*rfp,piv_rfp;
	struct Marker *mark;
	double *p,*p1,*pa,*p1a,*p2,*p2a,z,z1,za,pp,prob=0.0,z2,z3,z4,z5;

	ids=element->sire-1;
	idd=element->dam-1;
	comp=id_array[ids].comp;
	/* Get information about marker */
	mark=marker+locus;
	n_all=mark->n_all1[comp];
	n_gen=n_all*(n_all+1)/2;
	pivot=element->pivot-1;
	rfp=element->rf;
	id_array[ids].rfp=rfp[0];
	id_array[idd].rfp=rfp[1];
	piv_rfp=element->out_index;
	parflag=(founder_flag[locus][ids]&&rfp[0]<0)?1:0;
	parflag|=(founder_flag[locus][idd]&&rfp[1]<0)?2:0;
	/* Create functions for parents */
	/* Get penetrance and allele frequency info. (if founder) only (i.e., no previous R-Functions) for sire */
	id_array[ids].rfp=-1;
	z=get_par_probs(pval,ids,locus,pen,ppar,freq);
	prob+=log(z);
	sire_rfp=rfp[0];
	/* Do we need to peel back to sire ? */
	if(sire_rfp>=0) {
		p=prf;
		p1=simple_rf[sire_rfp];
		for(i=0;i<n_gen;i++) {
			z=pval[i];
			z1=*p1++;
			*p++=z*z1;
			pval[i]*=z1;
		}
	}
	/* If sire is pivot, multiply pval by pivot R-Function */
	if(ids==pivot) {
		p=pval;
		p1=simple_rf[piv_rfp];
		for(i=0;i<n_gen;i++) (*p++)*=(*p1++);
	}
	p=pval;
	for(z=0.0,i=0;i<n_gen;i++) z+=*p++;
	if(z<RESCALE_LIMIT) {
		prob+=log(z);
		p=pval;
		z=1.0/z;
		for(i=0;i<n_gen;i++) (*p++)*=z;
	}
	/* Same for dam */
	id_array[idd].rfp=-1;
	z=get_par_probs(mval,idd,locus,pen,ppar,freq);
	prob+=log(z);
	dam_rfp=rfp[1];
	if(dam_rfp>=0) {
		p=mrf;
		p1=simple_rf[dam_rfp];
		for(i=0;i<n_gen;i++) {
			z=mval[i];
			z1=*p1++;
			*p++=z1*z;
			mval[i]*=z1;
		}
	}
	if(idd==pivot) {
		p=mval;
		p1=simple_rf[piv_rfp];
		for(i=0;i<n_gen;i++) (*p++)*=(*p1++);
	}
	p=mval;
	for(z=0.0,i=0;i<n_gen;i++) z+=*p++;
	if(z<RESCALE_LIMIT) {
		prob+=log(z);
		p=mval;
		z=1.0/z;
		for(i=0;i<n_gen;i++) (*p++)*=z;
	}
	/* Create functions for offspring to be peeled */
	off=element->off;
	n_off=n_off1=element->n_off;
	peel_start=-1;
	/* Loop through kids */
	for(i=0;i<n_off;i++) {
		kid=off[i];
		p=kval[i];
		/* Get penetrance probs. */
		for(j=0;j<n_gen;j++) *p++=1.0;
		if(mark->haplo[kid]) {
			z=pen(kval[i],kid,locus,ppar);
			if(z<RESCALE_LIMIT) {
				/* Normalize (if necessary) */
				prob+=log(z);
				z=1.0/z;
				p=kval[i];
				for(j=0;j<n_gen;j++) (*p++)*=z;
			}
		}
		/* Peeling back to kid? */
		k=rfp[2+i];
		if(k>=0) {
			if(peel_start<0) peel_start=i;
			p=krf[i];
			p1=simple_rf[k];
			p2=kval[i];
			for(j=0;j<n_gen;j++) {
				z=p2[j];
				z1=*p1++;
				*p++=z1*z;
				p2[j]*=z1;
			}
		}
		/* If old pivot, add in pivot R-Function */
		if(kid==pivot) {
			if(k<0) n_off1--;
			p=kval[i];
			p1=simple_rf[piv_rfp];
			for(j=0;j<n_gen;j++) (*p++)*=(*p1++);
		}
		p=kval[i];
		for(z=0.0,j=0;j<n_gen;j++) z+=*p++;
		if(z<RESCALE_LIMIT) {
			prob+=log(z);
			p=kval[i];
			z=1.0/z;
			for(j=0;j<n_gen;j++) (*p++)*=z;
		}
	}
	if(peel_start<0) peel_start=n_off1;
	/* Blank kid output R-Functions and posterior probs. */
	for(i=0;i<n_off;i++) {
		if(rfp[2+i]>=0) {
			p1=simple_rf[rfp[2+i]];
			for(j=0;j<n_gen;j++) *p1++=0.0;
		}
		p1=kpost[i];
		for(j=0;j<n_gen;j++) *p1++=0.0;
	}
	/* Blank parental output R-Functions */
	if(sire_rfp>=0) {
		p1=simple_rf[sire_rfp];
		for(i=0;i<n_gen;i++) *p1++=0.0;
	}
	if(dam_rfp>=0) {
		p1=simple_rf[dam_rfp];
		for(i=0;i<n_gen;i++) *p1++=0.0;
	}
	/* Blank parental posterior probs. */
	p1=val;
	p2=val1;
	for(i=0;i<n_gen;i++) *p1++=*p2++=0.0;
	pp=0.0; 	/* Overall prob. of family */
	/* Loop through parental haplotypes */
	/* i,j - paternal alleles
	 * k,l - maternal alleles
	 * 
	 * The loop structures are complicated to reduce tests within the loops.
	 * This way we know if l<k=j<i etc. without testing
	 */
	p=pval;
	pa=mval;
	for(i1=ix1=i=0;i<n_all;i++) {
		i1+=i;
		for(j1=j=0;j<i;j++,ix1++) {
			j1+=j;
			z=*p++;
			za=*pa++;
			p1=mval;
			p1a=pval;
			g1=i1;
			g3=j1;
			for(ix2=k=0;k<=j;k++,g1++,g3++) {
				g2=i1;
				g4=j1;
				for(l=0;l<k;l++,g2++,g4++,ix2++) {
					/* l<k<=j<i */
					ADD_INFO_KIDS4;
					HANDLE_SIRE_DAM;
					HANDLE_KIDS4;
				}
				/* l=k<=j<i */
				ADD_INFO_KIDS2(g1,g3);
				HANDLE_SIRE_DAM;
				HANDLE_KIDS2(g1,g3);
				ix2++;
			}
			g3=k*(k-1)/2+j;
			for(;k<i;k++,g1++) {
				g3+=k;
				g2=i1;
				g4=j1;
				for(l=0;l<=j;l++,g2++,g4++,ix2++) {
					/* l<=j<k<i */
					ADD_INFO_KIDS4;
					HANDLE_SIRE_DAM;
					HANDLE_KIDS4;
				}
				g4=l*(l-1)/2+j;
				for(;l<k;l++,g2++,ix2++) {
					g4+=l;
					/* j<l<k<i */
					ADD_INFO_KIDS4;
					HANDLE_SIRE_DAM;
					HANDLE_KIDS4;
				}
				/* j<l=k<i */
				ADD_INFO_KIDS2(g1,g3);
				HANDLE_SIRE_DAM;
				HANDLE_KIDS2(g1,g3);
				ix2++;
			}
			g3+=k;
			g2=i1;
			g4=j1;
			for(l=0;l<j;l++,g2++,g4++,ix2++) {
				/* l<j<k=i */
				ADD_INFO_KIDS4;
				HANDLE_SIRE_DAM;
				HANDLE_KIDS4;
			}
			/* l=j<k=i */
			ADD_INFO_KIDS3;
			HANDLE_SIRE_DAMa;
			HANDLE_KIDS3;
		}
		z=*p++;
		za=*pa++;
		p1=mval;
		p1a=pval;
		g1=i1;
		for(ix2=k=0;k<i;k++,g1++) {
			g2=i1;
			for(l=0;l<k;l++,g2++,ix2++) {
				/* l<k<j=i */
				ADD_INFO_KIDS2(g1,g2);
				HANDLE_SIRE_DAM;
				HANDLE_KIDS2(g1,g2);
			}
			/* l=k<j=i */
			ADD_INFO_KIDS1;
			HANDLE_SIRE_DAM;
			HANDLE_KIDS1;
			ix2++;
		}
		g2=i1;
		for(l=0;l<k;l++,g2++,ix2++) {
			/* l<k=j=i */
			ADD_INFO_KIDS2(g1,g2);
			HANDLE_SIRE_DAM;
			HANDLE_KIDS2(g1,g2);
		}
		/* l=k=j=i */
		ADD_INFO_KIDS1;
		HANDLE_SIRE_DAMa;
		HANDLE_KIDS1;
		ix1++;
	}
	/* If we have founders, add to allele distribution estimate */
	if(parflag&1) { /* Sire */
		p=val;
		for(z=0.0,i=0;i<n_gen;i++) z+=*p++;
		z=1.0/z;
		p=val;
		l=id_array[ids].group-1;
		p1=nfreq[l];
		for(k=i=0;i<n_all;i++) {
			for(j=0;j<i;j++) {
				z1=z*(*p++);
				p1[i]+=z1;
				p1[j]+=z1;
			}
			p1[i]+=2.0*z*(*p++);
		}
	}
	if(parflag&2) { /* Dam */
		p=val1;
		for(z=0.0,i=0;i<n_gen;i++) z+=*p++;
		z=1.0/z;
		p=val1;
		l=id_array[idd].group-1;
		p1=nfreq[l];
		for(k=i=0;i<n_all;i++) {
			for(j=0;j<i;j++) {
				z1=z*(*p++);
				p1[i]+=z1;
				p1[j]+=z1;
			}
			p1[i]+=2.0*z*(*p++);
		}
	}
	for(i=0;i<n_off+2;i++) if(rfp[i]>=0) {
		SCALE_RFUNC(rfp[i],p,z,j,n_gen);
	}
	prob+=log(pp);
	return prob;
}
