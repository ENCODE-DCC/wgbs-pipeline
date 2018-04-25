#include <config.h>
#include <stdlib.h>
#ifdef USE_DMALLOC
# include <dmalloc.h>
#endif
#include <math.h>
#include <stdio.h>
#include <assert.h>

#include "ranlib.h"
#include "utils.h"
#include "loki.h"
#include "loki_peel.h"
#include "get_par_probs.h"
#include "loki_simple_peel.h"

/* Performs simple (i.e., nuclear family based) peeling operation */
double loki_simple_peelop(const struct Simple_Element *element,const int locus,const int s_flag,pen_func pen,
								  lk_ulong **a_set,double **freq,struct R_Func *rf,struct loki *loki)
{
	int ids,idd,i,j,k,k1,l,l1,i1,j2,m,n,pivot,fsp=0,n_off,*off,kid,gt[4],of=0,nb1,nmc,no2=0;
	int comp,n_all,n_idx,n_bits,*id_set1,*id_set2;
	double prob=0.0,*tp,p1,z,*tmp,*tmp1,*tpp1,*tpp2;
	double *qval,*pval,*mval,*pivval,*id_set;
	lk_ulong a,b,a1,b1,cm[2],*tmp_idx,*tmp_idx1,*cmm[2],mask;
	lk_ulong *tt_all;
	struct fset *peel_fs,*t_fset;
	struct peel_mem *work;
	struct Id_Record *id_array;
	struct Marker *mark;
	
	pivot=element->pivot-1;
	if(pivot== -3) return peel_to_par(element,locus,pen,a_set,rf,loki);
	ids=abs(element->sire)-1;
	if(ids>=0 && s_flag && pivot== -1) return loki_simple_sample(element,locus,pen,a_set,freq,rf,loki);
	work=&loki->peel->workspace;
	id_array=loki->pedigree->id_array;
	idd=abs(element->dam)-1;
	off=element->off;
	n_off=element->n_off;
	kid=abs(off[0])-1;
	comp=id_array[kid].comp;
	mark=loki->markers->marker+locus;
	n_all=mark->n_all1[comp];
	n_bits=num_bits(n_all);
	n_idx=1<<(n_bits+n_bits);
	mask=(1<<n_bits)-1;
	k=n_all*n_all;
	id_set=work->s2;
	qval=id_set+k;
	if(ids<0) { /* Peeling singletons */
		for(m=0;m<n_off;m++)	{
			kid=off[m]-1;
			for(i=0;i<n_idx;i++) qval[i]=0.0;
			p1=get_par_probs(qval,kid,mark,pen,a_set,freq,rf,loki);
			prob+=log(p1);
			if(s_flag) {
				do {
					z=ranf();
					p1=0.0;
					for(i=0;i<n_idx;i++) if(qval[i]>0.0) {
						p1+=qval[i];
						if(z<=p1) break;
					}
				} while(i==n_idx);
				id_array[kid].allele[X_MAT]=1+(i&mask);
				id_array[kid].allele[X_PAT]=1+((i>>n_bits)&mask);
				id_array[kid].flag|=(SAMPLED_MAT|SAMPLED_PAT);
				
				assert(id_array[kid].allele[X_MAT]>0 && id_array[kid].allele[X_MAT]<=n_all);
				assert(id_array[kid].allele[X_PAT]>0 && id_array[kid].allele[X_PAT]<=n_all);
			}
		}
#ifdef TRACE_PEEL
		if(CHK_PEEL(TRACE_LEVEL_2)) {
			(void) printf("Returning from %s() with %g\n",__func__,prob);
		}
#endif
		return prob;
	}
	pval=qval+n_idx;
	mval=pval+n_idx;
	pivval=mval+n_idx;
	peel_fs=work->s0;
	id_set1=work->s1;
	id_set2=id_set1+k;
	tt_all=work->s3;
	cmm[0]=mark->req_set[0];
	cmm[1]=mark->req_set[1];
	nb1=1<<n_bits;
	if(idd!=pivot && pivot!= -2 && element->dam>0) {
		p1=get_par_probs(mval,idd,mark,pen,a_set,freq,rf,loki);
		prob+=log(p1);
	} else {
		a1=mark->temp[X_MAT][idd];
		j=0;
		if((k=id_array[idd].rfp)>=0) { /* Insert Previously computed R_Func */
			while(a1) {
				if(a1&1) {
					a=a_set[idd][j];
					l=j;
					while(a) {
						if(a&1) mval[l]=0.0;
						a>>=1;
						l+=nb1;
					}
				}
				a1>>=1;
				j++;
			}
			for(j=0;j<rf[k].n_terms;j++) mval[rf[k].index[j]]=rf[k].p[j];
		} else {
			while(a1) {
				if(a1&1) {
					a=a_set[idd][j];
					l=j;
					while(a) {
						if(a&1) mval[l]=1.0;
						a>>=1;
						l+=nb1;
					}
				}
				a1>>=1;
				j++;
			}
		}
	}
	if(ids!=pivot && pivot!= -2 && element->sire>0) {
		p1=get_par_probs(pval,ids,mark,pen,a_set,freq,rf,loki);
		prob+=log(p1);
	} else {
		a1=mark->temp[X_MAT][ids];
		j=0;
		if((k=id_array[ids].rfp)>=0) { /* Insert Previously computed R_Func */
			while(a1) {
				if(a1&1) {
					a=a_set[ids][j];
					l=j;
					while(a) {
						if(a&1) pval[l]=0.0;
						a>>=1;
						l+=nb1;
					}
				}
				a1>>=1;
				j++;
			}
			for(j=0;j<rf[k].n_terms;j++) pval[rf[k].index[j]]=rf[k].p[j];
		} else {
			while(a1) {
				if(a1&1) {
					a=a_set[ids][j];
					l=j;
					while(a) {
						if(a&1) pval[l]=1.0;
						a>>=1;
						l+=nb1;
					}
				}
				a1>>=1;
				j++;
			}
		}
	}
	/* Construct set of possible parental genotype combinations */
	nmc=0;
	k=0;
	b1=mark->temp[X_MAT][idd];
	while(b1) {
		if(b1&1) {
			b=a_set[idd][k];
			l=0;
			m=k;
			while(b) {
				if(b&1) {
					id_set1[nmc]=k;
					id_set2[nmc]=l;
					id_set[nmc++]=mval[m];
				}
				b>>=1;
				l++;
				m+=nb1;
			}
		}
		b1>>=1;
		k++;
	}
	tmp_idx=tt_all;
	k1=0;
	for(m=0;m<n_off;m++) {
		kid=abs(off[m])-1;
		if(off[m]<0) id_array[kid].flag|=(SAMPLED_MAT|SAMPLED_PAT);
		if(id_array[kid].flag&SAMPLED_MAT) continue;
		k1++;
		cm[X_MAT]=mark->req_set[X_MAT][kid];
		cm[X_PAT]=mark->req_set[X_PAT][kid];
		for(i=0;i<n_all;i++) tmp_idx[i]=0;
		for(i1=1,i=0;i<n_all;i++,i1<<=1) {
			a=a_set[kid][i];
			if(a) {
				if(a&cm[X_PAT]) a|=cm[X_PAT];
				if(i1&cm[X_MAT]) {
					b=cm[X_MAT];
					j=0;
					while(b) {
						if(b&1) tmp_idx[j]=a;
						j++;
						b>>=1;
					}
				} else tmp_idx[i]=a;
			}
		}
		tmp_idx1=tt_all;
		for(k=0;k<no2;k++) {
			for(i=0;i<n_all;i++) if(tmp_idx[i]!=tmp_idx1[i]) break;
			if(i==n_all) break;
			tmp_idx1+=n_all;
		} 
		if(k==no2) {
			no2++;
			tmp_idx+=n_all;
		}
	}
	t_fset=peel_fs;
	a1=mark->temp[X_MAT][ids];
	i=0;
	switch(no2) {
	 case 0:
		while(a1) {
			if(a1&1) {
				a=a_set[ids][i];
				j=0;
				while(a) {
					if(a&1) {
						p1=pval[(j<<n_bits)|i];
						for(k1=0;k1<nmc;k1++) {
							t_fset->pat_gene[X_MAT]=i;
							t_fset->pat_gene[X_PAT]=j;
							t_fset->mat_gene[X_MAT]=id_set1[k1];
							t_fset->mat_gene[X_PAT]=id_set2[k1];
							(t_fset++)->p=p1*id_set[k1];
							fsp++;
						} 
					}
					a>>=1;
					j++;
				}
			}
			a1>>=1;
			i++;
		}
		break;
	 case 1:
		while(a1) {
			if(a1&1) {
				a=a_set[ids][i];
				j=0;
				while(a) {
					if(a&1) {
						p1=pval[(j<<n_bits)|i];
						b=(1<<i)|(1<<j);
						for(k1=0;k1<nmc;k1++) {
							k=id_set1[k1];
							l=id_set2[k1];
							if((tt_all[k]&b)||(tt_all[l]&b)) {
								t_fset->pat_gene[X_MAT]=i;
								t_fset->pat_gene[X_PAT]=j;
								t_fset->mat_gene[X_MAT]=k;
								t_fset->mat_gene[X_PAT]=l;
								(t_fset++)->p=p1*id_set[k1];
								fsp++;
							} 
						}
					}
					a>>=1;
					j++;
				}
			}
			a1>>=1;
			i++;
		}
		break;
	 case 2:
		tmp_idx=tt_all+n_all;
		while(a1) {
			if(a1&1) {
				a=a_set[ids][i];
				j=0;
				while(a) {
					if(a&1) {
						p1=pval[(j<<n_bits)|i];
						b=(1<<i)|(1<<j);
						for(k1=0;k1<nmc;k1++) {
							k=id_set1[k1];
							l=id_set2[k1];
							if(((tt_all[k]&b)||(tt_all[l]&b))&&((tmp_idx[k]&b)||(tmp_idx[l]&b))) {
								t_fset->pat_gene[X_MAT]=i;
								t_fset->pat_gene[X_PAT]=j;
								t_fset->mat_gene[X_MAT]=k;
								t_fset->mat_gene[X_PAT]=l;
								(t_fset++)->p=p1*id_set[k1];
								fsp++;
							} 
						}
					}
					a>>=1;
					j++;
				}
			}
			a1>>=1;
			i++;
		}
		break;
	 default:
		while(a1) {
			if(a1&1) {
				a=a_set[ids][i];
				j=0;
				while(a) {
					if(a&1) {
						p1=pval[(j<<n_bits)|i];
						b=(1<<i)|(1<<j);
						for(k1=0;k1<nmc;k1++) {
							k=id_set1[k1];
							l=id_set2[k1];
							m=no2;
							tmp_idx=tt_all;
							while(m--) {
								if(!((tmp_idx[k]&b)||(tmp_idx[l]&b))) break;
								tmp_idx+=n_all;
							}
							if(m<0) {
								t_fset->pat_gene[X_MAT]=i;
								t_fset->pat_gene[X_PAT]=j;
								t_fset->mat_gene[X_MAT]=k;
								t_fset->mat_gene[X_PAT]=l;
								(t_fset++)->p=p1*id_set[k1];
								fsp++;
							} 
						}
					}
					a>>=1;
					j++;
				}
			}
			a1>>=1;
			i++;
		}
	}
	/* Add contributions from non-pivot offspring */
	for(m=0;m<n_off;m++) {
		kid=abs(off[m])-1;
		if(kid==pivot) continue;
		tp=id_array[kid].tp;
		cm[0]=cmm[0][kid];
		cm[1]=cmm[1][kid];
		t_fset=peel_fs;
		if(id_array[kid].flag&SAMPLED_MAT) { /* If kid is fixed */
			j=id_array[kid].allele[X_MAT]-1;
			k=id_array[kid].allele[X_PAT]-1;
			l1=(k<<n_bits)|j;
			if(!(cm[0] || cm[1])) {
				for(n=0;n<fsp;n++) {
					i=t_fset->pat_gene[X_MAT];
					j=t_fset->pat_gene[X_PAT];
					k=t_fset->mat_gene[X_MAT];
					l=t_fset->mat_gene[X_PAT];
					i1=i<<n_bits;
					j2=j<<n_bits;
					z=0.0;
					if((i1|k)==l1) z+=tp[X_MM_PM];
					if((j2|k)==l1) z+=tp[X_MM_PP];
					if((i1|l)==l1) z+=tp[X_MP_PM];
					if((j2|l)==l1) z+=tp[X_MP_PP];
					(t_fset++)->p*=z;
				}
			} else {
				l=(n_all-1)<<n_bits;
				j2=1<<n_bits;
				for(i=i1=0,j=1;i<n_all;i++,j<<=1,i1+=j2) {
					id_set1[i]=(cm[X_PAT]&j)?l:i1;
					id_set2[i]=(cm[X_MAT]&j)?n_all-1:i;
				}
				for(n=0;n<fsp;n++) {
					i1=id_set1[t_fset->pat_gene[X_MAT]];
					j2=id_set1[t_fset->pat_gene[X_PAT]];
					k=id_set2[t_fset->mat_gene[X_MAT]];
					l=id_set2[t_fset->mat_gene[X_PAT]];
					z=0.0;
					if((i1|k)==l1) z+=tp[X_MM_PM];
					if((j2|k)==l1) z+=tp[X_MM_PP];
					if((i1|l)==l1) z+=tp[X_MP_PM];
					if((j2|l)==l1) z+=tp[X_MP_PP];
					(t_fset++)->p*=z;
				}
			}
		} else { /* Kid not fixed */
			if(!pen) {
				if((k=id_array[kid].rfp)>=0) { /* Insert Previously computed R_Func */
					for(j=0;j<n_all;j++) {
						tmp=qval+(j<<n_bits);
						for(l=0;l<n_all;l++) *(tmp++)=0.0;
					}
					i1=rf[k].n_terms;
					tmp=rf[k].p;
					for(j=0;j<i1;j++) qval[rf[k].index[j]]=*(tmp++);
				} else {
					j=0;
					a1=mark->temp[X_MAT][kid];
					while(a1) {
						if(a1&1) {
							a=a_set[kid][j];
							tmp=qval+j;
							k=0;
							while(a) {
								*tmp=(a&1)?1.0:0.0;
								a>>=1;
								tmp+=nb1;
								k++;
							}
							for(;k<n_all;k++) {
								*tmp=0.0;
								tmp+=nb1;
							}
						} else {
							tmp=qval+j;
							for(k=0;k<n_all;k++) {
								*tmp=0.0;
								tmp+=nb1;
							}
						}
						a1>>=1;
						j++;
					}
					for(;j<n_all;j++) {
						tmp=qval+j;
						for(k=0;k<n_all;k++) {
							*tmp=0.0;
							tmp+=nb1;
						}
					}
				}
			} else {
				tmp=qval;
				for(j=0;j<n_idx;j++) *(tmp++)=0.0;
				if((k=id_array[kid].rfp)>=0) { /* Insert Previously computed R_Func */
					i1=rf[k].n_terms;
					tmp=rf[k].p;
					for(j=0;j<i1;j++) qval[rf[k].index[j]]=*(tmp++);
				} else {
					j=0;
					a1=mark->temp[X_MAT][kid];
					while(a1) {
						if(a1&1) {
							a=a_set[kid][j];
							tmp=qval+j;
							while(a) {
								if(a&1) *tmp=1.0;
								a>>=1;
								tmp+=nb1;
							}
						}
						a1>>=1;
						j++;
					}
				}
				pen(qval,kid,&mark->locus,n_all,n_bits,loki);
			}
			if(!(cm[0] || cm[1])) {
				for(n=0;n<fsp;n++) {
					i1=t_fset->pat_gene[X_MAT]<<n_bits;
					k=t_fset->mat_gene[X_MAT];
					j2=t_fset->pat_gene[X_PAT]<<n_bits;
					l=t_fset->mat_gene[X_PAT];
					z=tp[X_MM_PM]*qval[i1|k]+tp[X_MM_PP]*qval[j2|k]+tp[X_MP_PM]*qval[i1|l]+tp[X_MP_PP]*qval[j2|l];
					(t_fset++)->p*=z;
				}
			} else {
				tpp1=id_array[kid].tpp[X_PAT];
				tpp2=id_array[kid].tpp[X_MAT];
				l=(n_all-1)<<n_bits;
				j2=1<<n_bits;
				for(i=i1=0,j=1;i<n_all;i++,j<<=1,i1+=j2) {
					id_set1[i]=(cm[X_PAT]&j)?l:i1;
					id_set2[i]=(cm[X_MAT]&j)?n_all-1:i;
				}
				for(n=0;n<fsp;n++) {
					i1=id_set1[t_fset->pat_gene[X_MAT]];
					j2=id_set1[t_fset->pat_gene[X_PAT]];
					k=id_set2[t_fset->mat_gene[X_MAT]];
					l=id_set2[t_fset->mat_gene[X_PAT]];
					if(i1!=j2) {
						if(k!=l) z=tp[X_MM_PM]*qval[i1|k]+tp[X_MM_PP]*qval[j2|k]+tp[X_MP_PM]*qval[i1|l]+tp[X_MP_PP]*qval[j2|l];
						else z=tpp1[X_MAT]*qval[i1|k]+tpp1[X_PAT]*qval[j2|k];
					} else if(k!=l) z=tpp2[X_MAT]*qval[i1|k]+tpp2[X_PAT]*qval[i1|l];
					else z=qval[i1|k];
					(t_fset++)->p*=z;
				}
			}
		}
	}
	assert(fsp);
	if(pivot== -2) { /* Peeling to joint on both parents */
		p1=0.0;
		t_fset=peel_fs;
		for(n=0;n<fsp;n++) p1+=(t_fset++)->p;
		prob+=log(p1);
		k=element->out_index;
		rf[k].n_ind=4;
		rf[k].n_terms=n;
		get_rf_memory(rf+k,n,MRK_MBLOCK,loki);
		t_fset=peel_fs;
		for(n=0;n<fsp;n++) {
			gt[0]=t_fset->mat_gene[X_MAT]+1;
			gt[1]=t_fset->mat_gene[X_PAT]+1;
			gt[2]=t_fset->pat_gene[X_MAT]+1;
			gt[3]=t_fset->pat_gene[X_PAT]+1;
			rf[k].index[n]=get_index1(4,gt,n_bits);
			rf[k].p[n]=(t_fset++)->p/p1;
		}
#ifdef TRACE_PEEL
		if(CHK_PEEL(TRACE_LEVEL_2)) {
			(void) printf("Returning from %s() with %g\n",__func__,prob);
		}
#endif
		return prob;
	}
	/* If pivot is an offspring, bring in previous R-Function and zero out
	 * illegal genotypes */
	if(pivot>=0 && pivot!=ids && pivot!=idd) {
		tp=id_array[pivot].tp;
		cm[0]=cmm[0][pivot];
		cm[1]=cmm[1][pivot];
		j=0;
		if((k=id_array[pivot].rfp)>=0) { /* Insert Previously computed R_Func */
			for(j=0;j<n_all;j++) {
				a=a_set[pivot][j];
				l=j;
				while(a) {
					if(a&1) pivval[l]=0.0;
					a>>=1;
					l+=nb1;
				}
			}
			for(j=0;j<rf[k].n_terms;j++) pivval[rf[k].index[j]]=rf[k].p[j];
		} else {
			j=0;
			for(j=0;j<n_all;j++) {
				a=a_set[pivot][j];
				l=j;
				while(a) {
					if(a&1) pivval[l]=1.0;
					a>>=1;
					l+=nb1;
				}
			}
		}
		of=1;
	} else {
		tp=0;
		cm[0]=cm[1]=0;
	}
	/* Assemble output function in qval */
	t_fset=peel_fs;
	if(pivot<0) {
		p1=0.0;
		for(n=0;n<fsp;n++) p1+=(t_fset++)->p;
		prob+=log(p1);
	} else {
		for(j=0;j<n_all;j++) {
			tmp=qval+(j<<n_bits);
			for(l=0;l<n_all;l++) *(tmp++)=0.0;
		}
		if(!(idd==pivot || ids==pivot)) {
			l=(n_all-1)<<n_bits;
			j2=1<<n_bits;
			for(i=i1=0,j=1;i<n_all;i++,j<<=1,i1+=j2) {
				id_set1[i]=(cm[X_PAT]&j)?l:i1;
				id_set2[i]=(cm[X_MAT]&j)?n_all-1:i;
			}
			for(n=0;n<fsp;n++) {
				i1=id_set1[t_fset->pat_gene[X_MAT]];
				j2=id_set1[t_fset->pat_gene[X_PAT]];
				k=id_set2[t_fset->mat_gene[X_MAT]];
				l=id_set2[t_fset->mat_gene[X_PAT]];
				z=(t_fset++)->p;
				qval[i1|k]+=tp[X_MM_PM]*z;
				qval[j2|k]+=tp[X_MM_PP]*z;
				qval[i1|l]+=tp[X_MP_PM]*z;
				qval[j2|l]+=tp[X_MP_PP]*z;
			}
		} else  {
			if(idd==pivot) {
				for(n=0;n<fsp;n++) {
					k=t_fset->mat_gene[X_MAT];
					l=t_fset->mat_gene[X_PAT];
					z=(t_fset++)->p;
					qval[(l<<n_bits)|k]+=z;
				}
			} else {
				for(n=0;n<fsp;n++) {
					i=t_fset->pat_gene[X_MAT];
					j=t_fset->pat_gene[X_PAT];
					z=(t_fset++)->p;
					qval[(j<<n_bits)|i]+=z;
				}
			}
		}
		i=0;
		p1=0.0;
		if(of) {
			for(j=0;j<n_all;j++) {
				a=a_set[pivot][j];
				l=j;
				while(a) {
					if(a&1) {
						p1+=(qval[l]*=pivval[l]);
						i++;
					}
					l+=nb1;
					a>>=1;
				}
			}
		} else {
			for(j=0;j<n_all;j++) {
				a=a_set[pivot][j];
				l=j;
				while(a) {
					if(a&1) {
						p1+=qval[l];
						i++;
					}
					l+=nb1;
					a>>=1;
				}
			}
		}
#ifdef DEBUG
		if(p1<=0.0) {
			fprintf(stderr,"Prob. %g in peeling operation for locus %s",p1,mark->name);
			fputc('\n',stderr);
			ABT_FUNC("Aborting\n");
		}
#endif
		prob+=log(p1);
		k=element->out_index;
		id_array[pivot].rfp=k;
		rf[k].n_ind=2;
		rf[k].n_terms=i;
#ifdef DEBUG
		if(!i) ABT_FUNC("Internal error - zero possible combinations\n");
#endif
		get_rf_memory(rf+k,i,MRK_MBLOCK,loki);
		p1=1.0/p1;
		tmp1=rf[k].p;
		tmp_idx=rf[k].index;
		for(j=0;j<n_all;j++) {
			a=a_set[pivot][j];
			l=j;
			while(a) {
				if(a&1) {
					z=qval[l];
					*(tmp1++)=z*p1;
					*(tmp_idx++)=(lk_ulong)l;
				}
				a>>=1;
				l+=nb1;
			}
		}
#ifdef TRACE_PEEL
		if(CHK_PEEL(TRACE_LEVEL_3)) {
			for(j=0;j<i;j++) {
				l=(int)rf[k].index[j];
				printf("%d %d %g\n",(int)(1+(l&mask)),1+(l>>n_bits),rf[k].p[j]/p1);
			}
		}
#endif
	} 
#ifdef TRACE_PEEL
	if(CHK_PEEL(TRACE_LEVEL_2)) {
		(void) printf("Returning from %s() with %g\n",__func__,prob);
	}
#endif
	return prob;
}

/* Performs simple (i.e., nuclear family based) peeling operation on x-linked data*/
double loki_simple_peelop_x(const struct Simple_Element *element,const int locus,const int s_flag,pen_func pen,
									 lk_ulong **a_set,double **freq,struct R_Func *rf,struct loki *loki)
{
	int ids,idd,i,j,k,k1,l,l1,i1,j2,m,n,pivot,fsp=0,n_off,*off,kid,gt[4],of=0,nb1,nmc,no2=0;
	int comp,n_all,n_idx,n_bits,*id_set1,*id_set2,sex;
	double prob=0.0,*tp,p1,z,*tmp,*tmp1,*tpp1,*tpp2;
	double *qval,*pval,*mval,*pivval,*id_set;
	lk_ulong a,b,a1,b1,cm[2],*tmp_idx,*tmp_idx1,*cmm[2],mask;
	lk_ulong *tt_all;
	struct fset *peel_fs,*t_fset;
	struct peel_mem *work;
	struct Id_Record *id_array;
	struct Marker *mark;
	
	pivot=element->pivot-1;
	if(pivot== -3) return peel_to_par(element,locus,pen,a_set,rf,loki);
	ids=element->sire-1;
	if(ids>=0 && s_flag && pivot== -1) return loki_simple_sample(element,locus,pen,a_set,freq,rf,loki);
	work=&loki->peel->workspace;
	id_array=loki->pedigree->id_array;
	idd=element->dam-1;
	off=element->off;
	n_off=element->n_off;
	kid=off[0]-1;
	comp=id_array[kid].comp;
	mark=loki->markers->marker+locus;
	n_all=mark->n_all1[comp];
	n_bits=num_bits(n_all);
	n_idx=1<<(n_bits+n_bits);
	mask=(1<<n_bits)-1;
	k=n_all*n_all;
	id_set=work->s2;
	qval=id_set+k;
#ifdef TRACE_PEEL
	if(CHK_PEEL(TRACE_LEVEL_1)) (void)printf("In %s(%p,%d,%d,%p)\n",__func__,(void *)element,locus,s_flag,(void *)pen);
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
		(void)fputs(") -> ",stdout);
		print_orig_id1(stdout,pivot+1);
		(void)fputc('\n',stdout);
	}
#endif
	if(ids<0) { /* Peeling singletons */
		for(m=0;m<n_off;m++)	{
			kid=off[m]-1;
			sex=id_array[kid].sex;	
			p1=get_par_probs(qval,kid,mark,pen,a_set,freq,rf,loki);
			prob+=log(p1);
			if(s_flag) {
				if(sex==1) {
					do {
						z=ranf();
						p1=0.0;
						for(i=0;i<n_all;i++) if(qval[i]>0.0) {
							p1+=qval[i];
							if(z<=p1) break;
						}
					} while(i==n_all);
					id_array[kid].allele[X_MAT]=1+i;
					id_array[kid].allele[X_PAT]=1;
				} else {
					do {
						z=ranf();
						p1=0.0;
						for(i=0;i<n_idx;i++) if(qval[i]>0.0) {
							p1+=qval[i];
							if(z<=p1) break;
						}
					} while(i==n_idx);
					id_array[kid].allele[X_MAT]=1+(i&mask);
					id_array[kid].allele[X_PAT]=1+((i>>n_bits)&mask);
				}
				id_array[kid].flag|=(SAMPLED_MAT|SAMPLED_PAT);
			}
		}
#ifdef TRACE_PEEL
		if(CHK_PEEL(TRACE_LEVEL_2)) {
			(void) printf("Returning from %s() with %g\n",__func__,prob);
		}
#endif
		return prob;
	}
	pval=qval+n_idx;
	mval=pval+n_idx;
	pivval=mval+n_idx;
	peel_fs=work->s0;
	id_set1=work->s1;
	id_set2=id_set1+k;
	tt_all=work->s3;
	cmm[0]=mark->req_set[0];
	cmm[1]=mark->req_set[1];
	nb1=1<<n_bits;
	if(idd!=pivot && pivot!= -2) {
		p1=get_par_probs(mval,idd,mark,pen,a_set,freq,rf,loki);
		prob+=log(p1);
	} else {
		a1=mark->temp[X_MAT][idd];
		j=0;
		if((k=id_array[idd].rfp)>=0) { /* Insert Previously computed R_Func */
			while(a1) {
				if(a1&1) {
					a=a_set[idd][j];
					l=j;
					while(a) {
						if(a&1) mval[l]=0.0;
						a>>=1;
						l+=nb1;
					}
				}
				a1>>=1;
				j++;
			}
			for(j=0;j<rf[k].n_terms;j++) mval[rf[k].index[j]]=rf[k].p[j];
		} else {
			while(a1) {
				if(a1&1) {
					a=a_set[idd][j];
					l=j;
					while(a) {
						if(a&1) mval[l]=1.0;
						a>>=1;
						l+=nb1;
					}
				}
				a1>>=1;
				j++;
			}
		}
	}
	if(ids!=pivot && pivot!= -2) {
		p1=get_par_probs(pval,ids,mark,pen,a_set,freq,rf,loki);
		prob+=log(p1);
	} else {
		a1=mark->temp[X_MAT][ids];
		j=0;
		if((k=id_array[ids].rfp)>=0) { /* Insert Previously computed R_Func */
			while(a1) {
				if(a1&1) pval[j]=0.0;
				a1>>=1;
				j++;
			}
			for(j=0;j<rf[k].n_terms;j++) pval[rf[k].index[j]]=rf[k].p[j];
		} else {
			while(a1) {
				if(a1&1) pval[j]=1.0;
				a1>>=1;
				j++;
			}
		}
	}
	/* Construct set of possible parental genotype combinations */
	nmc=0;
	k=0;
	/* First, collect maternal combinations */
	b1=mark->temp[X_MAT][idd];
	while(b1) {
		if(b1&1) {
			b=a_set[idd][k];
			l=0;
			m=k;
			while(b) {
				if(b&1) {
					id_set1[nmc]=k;
					id_set2[nmc]=l;
					id_set[nmc++]=mval[m];
				}
				b>>=1;
				l++;
				m+=nb1;
			}
		}
		b1>>=1;
		k++;
	}
	/* Check requirements from each kid */
	tmp_idx=tt_all;
	k1=0;
	for(m=0;m<n_off;m++) {
		kid=off[m]-1;
		if(id_array[kid].flag&SAMPLED_MAT) continue;
		k1++;
		cm[X_MAT]=mark->req_set[X_MAT][kid];
		cm[X_PAT]=mark->req_set[X_PAT][kid];
		for(i=0;i<n_all;i++) tmp_idx[i]=0;
		for(i1=1,i=0;i<n_all;i++,i1<<=1) {
			a=a_set[kid][i];
			if(a) {
				if(a&cm[X_PAT]) a|=cm[X_PAT];
				if(i1&cm[X_MAT]) {
					b=cm[X_MAT];
					j=0;
					while(b) {
						if(b&1) tmp_idx[j]=a;
						j++;
						b>>=1;
					}
				} else tmp_idx[i]=a;
			}
		}
		tmp_idx1=tt_all;
		for(k=0;k<no2;k++) {
			for(i=0;i<n_all;i++) if(tmp_idx[i]!=tmp_idx1[i]) break;
			if(i==n_all) break;
			tmp_idx1+=n_all;
		} 
		if(k==no2) {
			no2++;
			tmp_idx+=n_all;
		}
	}
	t_fset=peel_fs;
	a1=mark->temp[X_MAT][ids];
	i=0;
	switch(no2) {
	 case 0:
		while(a1) {
			if(a1&1) {
				a=a_set[ids][i];
				j=0;
				while(a) {
					if(a&1) {
						p1=pval[(j<<n_bits)|i];
						for(k1=0;k1<nmc;k1++) {
							t_fset->pat_gene[X_MAT]=i;
							t_fset->pat_gene[X_PAT]=j;
							t_fset->mat_gene[X_MAT]=id_set1[k1];
							t_fset->mat_gene[X_PAT]=id_set2[k1];
							(t_fset++)->p=p1*id_set[k1];
							fsp++;
						} 
					}
					a>>=1;
					j++;
				}
			}
			a1>>=1;
			i++;
		}
		break;
	 case 1:
		while(a1) {
			if(a1&1) {
				a=a_set[ids][i];
				j=0;
				while(a) {
					if(a&1) {
						p1=pval[(j<<n_bits)|i];
						b=(1<<i)|(1<<j);
						for(k1=0;k1<nmc;k1++) {
							k=id_set1[k1];
							l=id_set2[k1];
							if((tt_all[k]&b)||(tt_all[l]&b)) {
								t_fset->pat_gene[X_MAT]=i;
								t_fset->pat_gene[X_PAT]=j;
								t_fset->mat_gene[X_MAT]=k;
								t_fset->mat_gene[X_PAT]=l;
								(t_fset++)->p=p1*id_set[k1];
								fsp++;
							} 
						}
					}
					a>>=1;
					j++;
				}
			}
			a1>>=1;
			i++;
		}
		break;
	 case 2:
		tmp_idx=tt_all+n_all;
		while(a1) {
			if(a1&1) {
				a=a_set[ids][i];
				j=0;
				while(a) {
					if(a&1) {
						p1=pval[(j<<n_bits)|i];
						b=(1<<i)|(1<<j);
						for(k1=0;k1<nmc;k1++) {
							k=id_set1[k1];
							l=id_set2[k1];
							if(((tt_all[k]&b)||(tt_all[l]&b))&&((tmp_idx[k]&b)||(tmp_idx[l]&b))) {
								t_fset->pat_gene[X_MAT]=i;
								t_fset->pat_gene[X_PAT]=j;
								t_fset->mat_gene[X_MAT]=k;
								t_fset->mat_gene[X_PAT]=l;
								(t_fset++)->p=p1*id_set[k1];
								fsp++;
							} 
						}
					}
					a>>=1;
					j++;
				}
			}
			a1>>=1;
			i++;
		}
		break;
	 default:
		while(a1) {
			if(a1&1) {
				p1=pval[i];
				b=(1<<i);
				for(k1=0;k1<nmc;k1++) {
					k=id_set1[k1];
					l=id_set2[k1];
					m=no2;
					tmp_idx=tt_all;
					while(m--) {
						if(!(tmp_idx[k]&tmp_idx[l]&b)) break;
						tmp_idx+=n_all;
					}
					if(m<0) {
						t_fset->pat_gene[X_MAT]=i;
						t_fset->mat_gene[X_MAT]=k;
						t_fset->mat_gene[X_PAT]=l;
						(t_fset++)->p=p1*id_set[k1];
						fsp++;
					}
				}
			}
			a1>>=1;
			i++;
		}
	}
	/* Add contributions from non-pivot offspring */
	for(m=0;m<n_off;m++) {
		kid=off[m]-1;
		if(kid==pivot) continue;
		tp=id_array[kid].tp;
		cm[0]=cmm[0][kid];
		cm[1]=cmm[1][kid];
		t_fset=peel_fs;
		if(id_array[kid].flag&SAMPLED_MAT) { /* If kid is fixed */
			j=id_array[kid].allele[X_MAT]-1;
			k=id_array[kid].allele[X_PAT]-1;
			l1=(k<<n_bits)|j;
			if(!(cm[0] || cm[1])) {
				for(n=0;n<fsp;n++) {
					i=t_fset->pat_gene[X_MAT];
					j=t_fset->pat_gene[X_PAT];
					k=t_fset->mat_gene[X_MAT];
					l=t_fset->mat_gene[X_PAT];
					i1=i<<n_bits;
					j2=j<<n_bits;
					z=0.0;
					if((i1|k)==l1) z+=tp[X_MM_PM];
					if((j2|k)==l1) z+=tp[X_MM_PP];
					if((i1|l)==l1) z+=tp[X_MP_PM];
					if((j2|l)==l1) z+=tp[X_MP_PP];
					(t_fset++)->p*=z;
				}
			} else {
				l=(n_all-1)<<n_bits;
				j2=1<<n_bits;
				for(i=i1=0,j=1;i<n_all;i++,j<<=1,i1+=j2) {
					id_set1[i]=(cm[X_PAT]&j)?l:i1;
					id_set2[i]=(cm[X_MAT]&j)?n_all-1:i;
				}
				for(n=0;n<fsp;n++) {
					i1=id_set1[t_fset->pat_gene[X_MAT]];
					j2=id_set1[t_fset->pat_gene[X_PAT]];
					k=id_set2[t_fset->mat_gene[X_MAT]];
					l=id_set2[t_fset->mat_gene[X_PAT]];
					z=0.0;
					if((i1|k)==l1) z+=tp[X_MM_PM];
					if((j2|k)==l1) z+=tp[X_MM_PP];
					if((i1|l)==l1) z+=tp[X_MP_PM];
					if((j2|l)==l1) z+=tp[X_MP_PP];
					(t_fset++)->p*=z;
				}
			}
		} else { /* Kid not fixed */
			if(!pen) {
				if((k=id_array[kid].rfp)>=0) { /* Insert Previously computed R_Func */
					for(j=0;j<n_all;j++) {
						tmp=qval+(j<<n_bits);
						for(l=0;l<n_all;l++) *(tmp++)=0.0;
					}
					i1=rf[k].n_terms;
					tmp=rf[k].p;
					for(j=0;j<i1;j++) qval[rf[k].index[j]]=*(tmp++);
				} else {
					j=0;
					a1=mark->temp[X_MAT][kid];
					while(a1) {
						if(a1&1) {
							a=a_set[kid][j];
							tmp=qval+j;
							k=0;
							while(a) {
								*tmp=(a&1)?1.0:0.0;
								a>>=1;
								tmp+=nb1;
								k++;
							}
							for(;k<n_all;k++) {
								*tmp=0.0;
								tmp+=nb1;
							}
						} else {
							tmp=qval+j;
							for(k=0;k<n_all;k++) {
								*tmp=0.0;
								tmp+=nb1;
							}
						}
						a1>>=1;
						j++;
					}
					for(;j<n_all;j++) {
						tmp=qval+j;
						for(k=0;k<n_all;k++) {
							*tmp=0.0;
							tmp+=nb1;
						}
					}
				}
			} else {
				tmp=qval;
				for(j=0;j<n_idx;j++) *(tmp++)=0.0;
				if((k=id_array[kid].rfp)>=0) { /* Insert Previously computed R_Func */
					i1=rf[k].n_terms;
					tmp=rf[k].p;
					for(j=0;j<i1;j++) qval[rf[k].index[j]]=*(tmp++);
				} else {
					j=0;
					a1=mark->temp[X_MAT][kid];
					while(a1) {
						if(a1&1) {
							a=a_set[kid][j];
							tmp=qval+j;
							while(a) {
								if(a&1) *tmp=1.0;
								a>>=1;
								tmp+=nb1;
							}
						}
						a1>>=1;
						j++;
					}
				}
				pen(qval,kid,&mark->locus,n_all,n_bits,loki);
			}
			if(!(cm[0] || cm[1])) {
				for(n=0;n<fsp;n++) {
					i1=t_fset->pat_gene[X_MAT]<<n_bits;
					k=t_fset->mat_gene[X_MAT];
					j2=t_fset->pat_gene[X_PAT]<<n_bits;
					l=t_fset->mat_gene[X_PAT];
					z=tp[X_MM_PM]*qval[i1|k]+tp[X_MM_PP]*qval[j2|k]+tp[X_MP_PM]*qval[i1|l]+tp[X_MP_PP]*qval[j2|l];
					(t_fset++)->p*=z;
				}
			} else {
				tpp1=id_array[kid].tpp[X_PAT];
				tpp2=id_array[kid].tpp[X_MAT];
				l=(n_all-1)<<n_bits;
				j2=1<<n_bits;
				for(i=i1=0,j=1;i<n_all;i++,j<<=1,i1+=j2) {
					id_set1[i]=(cm[X_PAT]&j)?l:i1;
					id_set2[i]=(cm[X_MAT]&j)?n_all:i;
				}
				for(n=0;n<fsp;n++) {
					i1=id_set1[t_fset->pat_gene[X_MAT]];
					j2=id_set1[t_fset->pat_gene[X_PAT]];
					k=id_set2[t_fset->mat_gene[X_MAT]];
					l=id_set2[t_fset->mat_gene[X_PAT]];
					if(i1!=j2) {
						if(k!=l) z=tp[X_MM_PM]*qval[i1|k]+tp[X_MM_PP]*qval[j2|k]+tp[X_MP_PM]*qval[i1|l]+tp[X_MP_PP]*qval[j2|l];
						else z=tpp1[X_MAT]*qval[i1|k]+tpp1[X_PAT]*qval[j2|k];
					} else if(k!=l) z=tpp2[X_MAT]*qval[i1|k]+tpp2[X_PAT]*qval[i1|l];
					else z=qval[i1|k];
					(t_fset++)->p*=z;
				}
			}
		}
	}
#ifdef debug
	if(!fsp) ABT_FUNC("No possible parental combinations\n");
#endif
	if(pivot== -2) { /* Peeling to joint on both parents */
		p1=0.0;
		t_fset=peel_fs;
		for(n=0;n<fsp;n++) p1+=(t_fset++)->p;
		prob+=log(p1);
		k=element->out_index;
		rf[k].n_ind=4;
		rf[k].n_terms=n;
		get_rf_memory(rf+k,n,MRK_MBLOCK,loki);
		t_fset=peel_fs;
		for(n=0;n<fsp;n++) {
			gt[0]=t_fset->mat_gene[X_MAT]+1;
			gt[1]=t_fset->mat_gene[X_PAT]+1;
			gt[2]=t_fset->pat_gene[X_MAT]+1;
			gt[3]=t_fset->pat_gene[X_PAT]+1;
			rf[k].index[n]=get_index1(4,gt,n_bits);
			rf[k].p[n]=(t_fset++)->p/p1;
		}
#ifdef TRACE_PEEL
		if(CHK_PEEL(TRACE_LEVEL_2)) {
			(void) printf("Returning from %s() with %g\n",__func__,prob);
		}
#endif
		return prob;
	}
	/* If pivot is an offspring, bring in previous R-Function and zero out
	 * illegal genotypes */
	if(pivot>=0 && pivot!=ids && pivot!=idd) {
		tp=id_array[pivot].tp;
		cm[0]=cmm[0][pivot];
		cm[1]=cmm[1][pivot];
		j=0;
		if((k=id_array[pivot].rfp)>=0) { /* Insert Previously computed R_Func */
			for(j=0;j<n_all;j++) {
				a=a_set[pivot][j];
				l=j;
				while(a) {
					if(a&1) pivval[l]=0.0;
					a>>=1;
					l+=nb1;
				}
			}
			for(j=0;j<rf[k].n_terms;j++) pivval[rf[k].index[j]]=rf[k].p[j];
		} else {
			j=0;
			for(j=0;j<n_all;j++) {
				a=a_set[pivot][j];
				l=j;
				while(a) {
					if(a&1) pivval[l]=1.0;
					a>>=1;
					l+=nb1;
				}
			}
		}
		of=1;
	} else {
		tp=0;
		cm[0]=cm[1]=0;
	}
	/* Assemble output function in qval */
	t_fset=peel_fs;
	if(pivot<0) {
		p1=0.0;
		for(n=0;n<fsp;n++) p1+=(t_fset++)->p;
		prob+=log(p1);
	} else {
		for(j=0;j<n_all;j++) {
			tmp=qval+(j<<n_bits);
			for(l=0;l<n_all;l++) *(tmp++)=0.0;
		}
		if(!(idd==pivot || ids==pivot)) {
			l=(n_all-1)<<n_bits;
			j2=1<<n_bits;
			for(i=i1=0,j=1;i<n_all;i++,j<<=1,i1+=j2) {
				id_set1[i]=(cm[X_PAT]&j)?l:i1;
				id_set2[i]=(cm[X_MAT]&j)?n_all-1:i;
			}
			for(n=0;n<fsp;n++) {
				i1=id_set1[t_fset->pat_gene[X_MAT]];
				j2=id_set1[t_fset->pat_gene[X_PAT]];
				k=id_set2[t_fset->mat_gene[X_MAT]];
				l=id_set2[t_fset->mat_gene[X_PAT]];
				z=(t_fset++)->p;
				qval[i1|k]+=tp[X_MM_PM]*z;
				qval[j2|k]+=tp[X_MM_PP]*z;
				qval[i1|l]+=tp[X_MP_PM]*z;
				qval[j2|l]+=tp[X_MP_PP]*z;
			}
		} else  {
			if(idd==pivot) {
				for(n=0;n<fsp;n++) {
					k=t_fset->mat_gene[X_MAT];
					l=t_fset->mat_gene[X_PAT];
					z=(t_fset++)->p;
					qval[(l<<n_bits)|k]+=z;
				}
			} else {
				for(n=0;n<fsp;n++) {
					i=t_fset->pat_gene[X_MAT];
					j=t_fset->pat_gene[X_PAT];
					z=(t_fset++)->p;
					qval[(j<<n_bits)|i]+=z;
				}
			}
		}
		i=0;
		p1=0.0;
		if(of) {
			for(j=0;j<n_all;j++) {
				a=a_set[pivot][j];
				l=j;
				while(a) {
					if(a&1) {
						p1+=(qval[l]*=pivval[l]);
						i++;
					}
					l+=nb1;
					a>>=1;
				}
			}
		} else {
			for(j=0;j<n_all;j++) {
				a=a_set[pivot][j];
				l=j;
				while(a) {
					if(a&1) {
						p1+=qval[l];
						i++;
					}
					l+=nb1;
					a>>=1;
				}
			}
		}
#ifdef DEBUG
		if(p1<=0) {
			fprintf(stderr,"Prob. %g in peeling operation for locus %s",p1,mark->name);
			fputc('\n',stderr);
			ABT_FUNC("Aborting\n");
		}
#endif
		prob+=log(p1);
		k=element->out_index;
		id_array[pivot].rfp=k;
		rf[k].n_ind=2;
		rf[k].n_terms=i;
#ifdef DEBUG
		if(!i) ABT_FUNC("Internal error - zero possible combinations\n");
#endif
		get_rf_memory(rf+k,i,MRK_MBLOCK,loki);
		p1=1.0/p1;
		tmp1=rf[k].p;
		tmp_idx=rf[k].index;
		for(j=0;j<n_all;j++) {
			a=a_set[pivot][j];
			l=j;
			while(a) {
				if(a&1) {
					z=qval[l];
					*(tmp1++)=z*p1;
					*(tmp_idx++)=(lk_ulong)l;
				}
				a>>=1;
				l+=nb1;
			}
		}
#ifdef TRACE_PEEL
		if(CHK_PEEL(TRACE_LEVEL_3)) {
			for(j=0;j<i;j++) {
				l=(int)rf[k].index[j];
				printf("%d %d %g\n",(int)(1+(l&mask)),1+(l>>n_bits),rf[k].p[j]/p1);
			}
		}
#endif
		
	} 
#ifdef TRACE_PEEL
	if(CHK_PEEL(TRACE_LEVEL_2)) {
		(void) printf("Returning from %s() with %g\n",__func__,prob);
	}
#endif
	return prob;
}

