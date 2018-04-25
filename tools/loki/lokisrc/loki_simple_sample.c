#include <config.h>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <assert.h>

#include "config.h"
#include "ranlib.h"
#include "utils.h"
#include "loki.h"
#include "loki_peel.h"
#include "get_par_probs.h"
#include "loki_simple_peel.h"

/* Sample a nuclear family */
double loki_simple_sample(const struct Simple_Element *element,const int locus,pen_func pen,
								  lk_ulong **a_set,double **freq,struct R_Func *rf,struct loki *loki)
{
	int ids,idd,i,j,i1,j2,k,l,l1,m,n,fsp=0,n_off,*off,kid,*off_flag,no=0,no1=0,nb1,nmc,no2=0;
	int n_all,n_idx,comp,n_bits,*peel_pt2,*id_set1,*id_set2;
	double *tp,p1,p2,z,z1,prob=0.0,pp[4],*tmp,*tpp1,*tpp2;
	double *qval,*pval,*mval,*id_set;
	lk_ulong a,a2,b,cm[2],*cm1[2],cmp,cmm,*t_cm,*t_all1,*t_all2;
	lk_ulong *tt_all,*peel_pt1,a1;
	struct fset *peel_fs,*t_fset;
	struct Id_Record *id_array;
	struct peel_mem *work;
	struct Marker *mark;
	
	id_array=loki->pedigree->id_array;
	work=&loki->peel->workspace;
	ids=abs(element->sire)-1;
	idd=abs(element->dam)-1;
	comp=id_array[ids].comp;
	mark=loki->markers->marker+locus;
	n_all=mark->n_all1[comp];
	n_bits=num_bits(n_all);
	n_idx=1<<(n_bits+n_bits);
	off=element->off;
	n_off=element->n_off;
	nb1=1<<n_bits;
	k=n_all*n_all;
	peel_fs=work->s0;
	id_set1=work->s1;
	id_set2=id_set1+k;
	id_set=work->s2;
	qval=id_set+k;
	pval=qval+n_idx;
	mval=pval+n_idx;
	tt_all=work->s3;
	peel_pt1=tt_all+loki->peel->max_peel_off*n_all;
	peel_pt2=work->s4;
	cm1[0]=peel_pt1;
	cm1[1]=cm1[0]+n_off;
	off_flag=peel_pt2;
	if(element->sire>0) {
		p1=get_par_probs(pval,ids,mark,pen,a_set,freq,rf,loki);
		prob+=log(p1);
	} else id_array[ids].flag|=(SAMPLED_MAT|SAMPLED_PAT);
	if(id_array[ids].flag&SAMPLED_MAT) {
#ifdef TRACE_PEEL
		if(CHK_PEEL(TRACE_LEVEL_3)) (void)fputs("Sire sampled\n",stdout);
#endif
	} else {
#ifdef TRACE_PEEL
		if(CHK_PEEL(TRACE_LEVEL_3)) (void)printf("Sire total (b): %g\n",p1);
#endif
		for(m=0;m<n_off;m++) {
			kid=abs(off[m])-1;
			if(off[m]<0) id_array[kid].flag|=(SAMPLED_MAT|SAMPLED_PAT);
			cmm=mark->req_set[X_PAT][kid];
			if(id_array[kid].flag&SAMPLED_MAT) cm1[X_PAT][m]=1L<<(id_array[kid].allele[X_PAT]-1);
			else cm1[X_PAT][m]=mark->temp[X_PAT][kid];
			if(cmm&cm1[X_PAT][m]) cm1[X_PAT][m]|=cmm;
			off_flag[m]=1;
		}
		if(n_off>1) {
			for(i=1;i<n_off;i++)	{	
				a=cm1[X_PAT][i];
				for(j=0;j<i;j++) {
					b=cm1[X_PAT][j];
					if((a&b)==a) off_flag[j]=0;
					else if((a&b)==b) off_flag[i]=0;
				}
			}
			no1=0;
			for(i=0;i<n_off;i++) if(off_flag[i]) cm1[X_PAT][no1++]=cm1[X_PAT][i];
		} else no1=n_off;
	}
	if(element->dam>0) {
		p1=get_par_probs(mval,idd,mark,pen,a_set,freq,rf,loki);
		prob+=log(p1);
	} else id_array[idd].flag|=(SAMPLED_MAT|SAMPLED_PAT);
	if(id_array[idd].flag&SAMPLED_MAT) {
#ifdef TRACE_PEEL
		if(CHK_PEEL(TRACE_LEVEL_3)) (void)fputs("Dam sampled\n",stdout);
#endif
	} else {
#ifdef TRACE_PEEL
		if(CHK_PEEL(TRACE_LEVEL_3)) (void)printf("Dam total (b): %g\n",p1);
#endif
		for(m=0;m<n_off;m++) {
			kid=abs(off[m])-1;
 			cmm=mark->req_set[X_MAT][kid];
			if(id_array[kid].flag&SAMPLED_MAT) cm1[X_MAT][m]=1L<<(id_array[kid].allele[X_MAT]-1);
			else cm1[X_MAT][m]=mark->temp[X_MAT][kid];
			if(cmm&cm1[X_MAT][m]) cm1[X_MAT][m]|=cmm;
			off_flag[m]=1;
		}
		if(n_off>1) {
			for(i=1;i<n_off;i++)	{	
				a=cm1[X_MAT][i];
				for(j=0;j<i;j++) {
					b=cm1[X_MAT][j];
					if((a&b)==a) off_flag[j]=0;
					else if((a&b)==b) off_flag[i]=0;
				}
			}
			no=0;
			for(i=0;i<n_off;i++) if(off_flag[i]) cm1[X_MAT][no++]=cm1[X_MAT][i];
		} else no=n_off;
	}
	if(!((id_array[ids].flag&SAMPLED_MAT)&&(id_array[idd].flag&SAMPLED_MAT))) {
		t_all1=tt_all;
		for(m=0;m<n_off;m++) {
			kid=abs(off[m])-1;
			if(id_array[kid].flag&SAMPLED_MAT) continue;
			cm[X_MAT]=mark->req_set[X_MAT][kid];
			cm[X_PAT]=mark->req_set[X_PAT][kid];
			for(i=0;i<n_all;i++) t_all1[i]=0;
			for(a1=1L,i=0;i<n_all;i++,a1<<=1) {
				a=a_set[kid][i];
				if(a) {
					if(a&cm[X_PAT]) a|=cm[X_PAT];
					if(a1&cm[X_MAT]) {
						b=cm[X_MAT];
						j=0;
						while(b) {
							if(b&1) t_all1[j]=a;
							j++;
							b>>=1;
						}
					} else t_all1[i]=a;
				}
			}
			t_all2=tt_all;
			for(k=0;k<no2;k++) {
				for(i=0;i<n_all;i++) if(t_all1[i]!=t_all2[i]) break;
				if(i==n_all) break;
				t_all2+=n_all;
			} 
			if(k==no2) {
				no2++;
				t_all1+=n_all;
			} 
		}
	} 
	t_fset=peel_fs;
	if(id_array[ids].flag&SAMPLED_MAT) {
		if(id_array[idd].flag&SAMPLED_MAT) {
			t_fset->pat_gene[X_MAT]=id_array[ids].allele[X_MAT]-1;
			t_fset->pat_gene[X_PAT]=id_array[ids].allele[X_PAT]-1;
			t_fset->mat_gene[X_MAT]=id_array[idd].allele[X_MAT]-1;
			t_fset->mat_gene[X_PAT]=id_array[idd].allele[X_PAT]-1;
			t_fset->p=1.0;
			fsp=1;
		} else {
			i=id_array[ids].allele[X_MAT]-1;
			j=id_array[ids].allele[X_PAT]-1;
			cmp=(1L<<i)|(1L<<j);
			for(k=0;k<n_all;k++) {
				b=a_set[idd][k];
				l=0;
				l1=k;
				while(b) {
					if(b&1) {
						cmm=(1L<<k)|(1L<<l);
						t_cm=cm1[X_MAT];
						for(m=0;m<no;m++) if(!(*(t_cm++)&cmm)) break;
						if(m==no) {
							t_all1=tt_all;
							for(m=0;m<no2;m++,t_all1+=n_all) if(!((t_all1[k]&cmp)||(t_all1[l]&cmp))) break;
							if(m==no2) {
								t_fset->pat_gene[X_MAT]=i;
								t_fset->pat_gene[X_PAT]=j;
								t_fset->mat_gene[X_MAT]=k;
								t_fset->mat_gene[X_PAT]=l;
								(t_fset++)->p=mval[l1];
								fsp++;
							}
						}
					}
					b>>=1;
					l++;
					l1+=nb1;
				}
			}
		}
	} else if(id_array[idd].flag&SAMPLED_MAT) {
		k=id_array[idd].allele[X_MAT]-1;
		l=id_array[idd].allele[X_PAT]-1;
		for(i=0;i<n_all;i++) {
			a=a_set[ids][i];
			j=0;
			l1=i;
			while(a) {
				if(a&1) {
					cmp=(1L<<i)|(1L<<j);
					t_cm=cm1[X_PAT];
					for(m=0;m<no1;m++) if(!(*(t_cm++)&cmp)) break;
					if(m==no1) {
						t_all1=tt_all;
						for(m=0;m<no2;m++,t_all1+=n_all) if(!((t_all1[k]&cmp)||(t_all1[l]&cmp))) break;
						if(m==no2) {
							t_fset->pat_gene[X_MAT]=i;
							t_fset->pat_gene[X_PAT]=j;
							t_fset->mat_gene[X_MAT]=k;
							t_fset->mat_gene[X_PAT]=l;
							(t_fset++)->p=pval[l1];
							fsp++;
						}
					}
				}
				a>>=1;
				j++;
				l1+=nb1;
			}
		}
	} else {
		for(nmc=k=0;k<n_all;k++) {
			b=a_set[idd][k];
			l=0;
			l1=k;
			while(b) {
				if(b&1) {
				cmm=(1L<<k)|(1L<<l);
					t_cm=cm1[X_MAT];
					for(m=0;m<no;m++) if(!(*(t_cm++)&cmm)) break;
					if(m==no) {
						id_set[nmc]=mval[l1];
						id_set1[nmc]=k;
						id_set2[nmc++]=l;
					}
				}
				b>>=1;
				l++;
				l1+=nb1;
			}
		}
		switch(no2) {
		 case 1:
			for(a1=1,i=0;i<n_all;i++,a1<<=1) {
				a=a_set[ids][i];
				j=0;
				a2=1;
				while(a) {
					if(a&1) {
						cmp=a1|a2;
						t_cm=cm1[X_PAT];
						for(m=0;m<no1;m++) if(!(*(t_cm++)&cmp)) break;
						if(m==no1) {
							z=pval[(j<<n_bits)|i];
							for(l1=0;l1<nmc;l1++) {
								k=id_set1[l1];
								l=id_set2[l1];
								if((tt_all[k]&cmp)||(tt_all[l]&cmp)) {
									t_fset->pat_gene[X_MAT]=i;
									t_fset->pat_gene[X_PAT]=j;
									t_fset->mat_gene[X_MAT]=k;
									t_fset->mat_gene[X_PAT]=l;
									(t_fset++)->p=z*id_set[l1];
									fsp++;
								} 
							}
						}
					}
					a>>=1;
					j++;
					a2<<=1;
				}
			}
			break;
		 case 2:
			t_all1=tt_all+n_all;
			for(a1=1,i=0;i<n_all;i++,a1<<=1) {
				a=a_set[ids][i];
				j=0;
				a2=1;
				while(a) {
					if(a&1) {
						cmp=a1|a2;
						t_cm=cm1[X_PAT];
						for(m=0;m<no1;m++) if(!(*(t_cm++)&cmp)) break;
						if(m==no1) {
							z=pval[(j<<n_bits)|i];
							for(l1=0;l1<nmc;l1++) {
								k=id_set1[l1];
								l=id_set2[l1];
								if(((tt_all[k]&cmp)||(tt_all[l]&cmp))&&((t_all1[k]&cmp)||(t_all1[l]&cmp))) {
									t_fset->pat_gene[X_MAT]=i;
									t_fset->pat_gene[X_PAT]=j;
									t_fset->mat_gene[X_MAT]=k;
									t_fset->mat_gene[X_PAT]=l;
									(t_fset++)->p=z*id_set[l1];
									fsp++;
								} 
							}
						}
					}
					a>>=1;
					j++;
					a2<<=1;
				}
			}
			break;
		 case 0:
			for(a1=1,i=0;i<n_all;i++,a1<<=1) {
				a=a_set[ids][i];
				j=0;
				a2=1;
				while(a) {
					if(a&1) {
						cmp=a1|a2;
						t_cm=cm1[X_PAT];
						for(m=0;m<no1;m++) if(!(*(t_cm++)&cmp)) break;
						if(m==no1) {
							z=pval[(j<<n_bits)|i];
							for(l1=0;l1<nmc;l1++) {
								t_fset->pat_gene[X_MAT]=i;
								t_fset->pat_gene[X_PAT]=j;
								t_fset->mat_gene[X_MAT]=id_set1[l1];
								t_fset->mat_gene[X_PAT]=id_set2[l1];
								(t_fset++)->p=z*id_set[l1];
								fsp++;
							}
						}
					}
					a>>=1;
					j++;
					a2<<=1;
				}
			}
			break;
		 default:
			for(a1=1,i=0;i<n_all;i++,a1<<=1) {
				a=a_set[ids][i];
				j=0;
				a2=1;
				while(a) {
					if(a&1) {
						cmp=a1|a2;
						t_cm=cm1[X_PAT];
						for(m=0;m<no1;m++) if(!(*(t_cm++)&cmp)) break;
						if(m==no1) {
							z=pval[(j<<n_bits)|i];
							for(l1=0;l1<nmc;l1++) {
								k=id_set1[l1];
								l=id_set2[l1];
								t_all1=tt_all;
								for(m=0;m<no2;m++,t_all1+=n_all) if(!((t_all1[k]&cmp)||(t_all1[l]&cmp))) break;
								if(m==no2) {
									t_fset->pat_gene[X_MAT]=i;
									t_fset->pat_gene[X_PAT]=j;
									t_fset->mat_gene[X_MAT]=k;
									t_fset->mat_gene[X_PAT]=l;
									(t_fset++)->p=z*id_set[l1];
									fsp++;
								} 
							}
						}
					}
					a>>=1;
					j++;
					a2<<=1;
				}
			}
		}
	}
	assert(fsp);
#ifndef NDEBUG
	for(k=0;k<fsp;k++) if(isnan(peel_fs[k].p)) break;
	if(k<fsp) {
		printf("OOOK!\n");
		for(k=0;k<fsp;k++) printf("%d %g\n",k,peel_fs[k].p);
		ABT_FUNC("Yo!\n");
	}
#endif
	for(m=0;m<n_off;m++) {
		kid=abs(off[m])-1;
		tp=id_array[kid].tp;
		cm[0]=mark->req_set[0][kid];
		cm[1]=mark->req_set[1][kid];
		t_fset=peel_fs;
		if(id_array[kid].flag&SAMPLED_MAT) {
			j=id_array[kid].allele[X_MAT]-1;
			k=id_array[kid].allele[X_PAT]-1;
			l1=(k<<n_bits)|j;
			if(!(cm[0] || cm[1])) for(n=0;n<fsp;n++) {
				i1=(t_fset->pat_gene[X_MAT]<<n_bits);
				j2=(t_fset->pat_gene[X_PAT]<<n_bits);
				k=peel_fs[n].mat_gene[X_MAT];
				l=peel_fs[n].mat_gene[X_PAT];
				z=0.0;
				if((i1|k)==l1) z+=tp[X_MM_PM];
				if((j2|k)==l1) z+=tp[X_MM_PP];
				if((i1|l)==l1) z+=tp[X_MP_PM];
				if((j2|l)==l1) z+=tp[X_MP_PP];
				assert(!isnan(z));
				(t_fset++)->p*=z;
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
					assert(!isnan(z));
					(t_fset++)->p*=z;
				}
			}
		} else {
			if(!pen) {
				if((k=id_array[kid].rfp)>=0) {/* Insert Previously computed R_Func */
					for(i1=0;i1<n_all;i1++) {
						tmp=qval+(i1<<n_bits);
						for(l=0;l<n_all;l++) *(tmp++)=0.0;
					}
					i1=rf[k].n_terms;
					tmp=rf[k].p;
					for(j=0;j<i1;j++) qval[rf[k].index[j]]=*(tmp++);
				} else {
					for(j=0;j<n_all;j++) {
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
					for(j=0;j<n_all;j++) {
						a=a_set[kid][j];
						tmp=qval+j;
						while(a) {
							if(a&1) *tmp=1.0;
							a>>=1;
							tmp+=nb1;
						}
					}
				}
				pen(qval,kid,&mark->locus,n_all,n_bits,loki);
			}
			if(!(cm[0] || cm[1])) for(n=0;n<fsp;n++) {
				i1=(t_fset->pat_gene[X_MAT])<<n_bits;
				j2=(t_fset->pat_gene[X_PAT])<<n_bits;
				k=t_fset->mat_gene[X_MAT];
				l=t_fset->mat_gene[X_PAT];
				z=tp[X_MM_PM]*qval[i1|k]+tp[X_MM_PP]*qval[j2|k]+tp[X_MP_PM]*qval[i1|l]+tp[X_MP_PP]*qval[j2|l];
				assert(!isnan(z));
				(t_fset++)->p*=z;
			} else {
				l=(n_all-1)<<n_bits;
				j2=1<<n_bits;
				for(i=i1=0,j=1;i<n_all;i++,j<<=1,i1+=j2) {
					id_set1[i]=(cm[X_PAT]&j)?l:i1;
					id_set2[i]=(cm[X_MAT]&j)?n_all-1:i;
				}
				tpp1=id_array[kid].tpp[X_PAT];
				tpp2=id_array[kid].tpp[X_MAT];
				for(n=0;n<fsp;n++) {
					i1=id_set1[t_fset->pat_gene[X_MAT]];
					j2=id_set1[t_fset->pat_gene[X_PAT]];
					k=id_set2[t_fset->mat_gene[X_MAT]];
					l=id_set2[t_fset->mat_gene[X_PAT]];
					z=tp[X_MM_PM]*qval[i1|k]+tp[X_MM_PP]*qval[j2|k]+tp[X_MP_PM]*qval[i1|l]+tp[X_MP_PP]*qval[j2|l];
					assert(!isnan(z));
					(t_fset++)->p*=z; 
				}
			}
		}
#ifndef NDEBUG
	for(k=0;k<fsp;k++) assert(!isnan(peel_fs[k].p));
#endif
	}
	p1=0.0;
	for(n=0;n<fsp;n++) p1+=peel_fs[n].p;
	assert(fsp);
	assert(p1>0.0 && !isnan(p1));
	prob+=log(p1);
	if(!((id_array[ids].flag&SAMPLED_MAT)&&(id_array[idd].flag&SAMPLED_MAT))) {
		do {
			z=ranf()*p1;
			p2=0.0;
			t_fset=peel_fs;
			for(n=0;n<fsp;n++,t_fset++) {
				p2+=t_fset->p;
				if(p2>=z) {
					id_array[ids].allele[X_MAT]=t_fset->pat_gene[X_MAT]+1;
					id_array[ids].allele[X_PAT]=t_fset->pat_gene[X_PAT]+1;
					id_array[idd].allele[X_MAT]=t_fset->mat_gene[X_MAT]+1;
					id_array[idd].allele[X_PAT]=t_fset->mat_gene[X_PAT]+1;
					id_array[ids].flag|=(SAMPLED_MAT|SAMPLED_PAT);
					id_array[idd].flag|=(SAMPLED_MAT|SAMPLED_PAT);
					assert(id_array[ids].allele[X_MAT]>0 && id_array[ids].allele[X_MAT]<=n_all);
					assert(id_array[ids].allele[X_PAT]>0 && id_array[ids].allele[X_PAT]<=n_all);
					assert(id_array[idd].allele[X_MAT]>0 && id_array[idd].allele[X_MAT]<=n_all);
					assert(id_array[idd].allele[X_PAT]>0 && id_array[idd].allele[X_PAT]<=n_all);
					break;
				}
			}
		} while(n==fsp);
	}
	for(m=0;m<n_off;m++) {
		kid=abs(off[m])-1;
		if(id_array[kid].flag&SAMPLED_MAT) continue;
		if(!pen) {
			if((k=id_array[kid].rfp)>=0) {/* Insert Previously computed R_Func */
				for(i1=0;i1<n_all;i1++) {
					tmp=qval+(i1<<n_bits);
					for(l=0;l<n_all;l++) *(tmp++)=0.0;
				}
				i1=rf[k].n_terms;
				tmp=rf[k].p;
				for(j=0;j<i1;j++) qval[rf[k].index[j]]=*(tmp++); 
			} else {
				for(j=0;j<n_all;j++) {
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
				}
			}
		} else {
			tmp=qval;
			for(j=0;j<n_idx;j++) *(tmp++)=0.0;
			if((k=id_array[kid].rfp)>=0) {/* Insert Previously computed R_Func */
				i1=rf[k].n_terms;
				tmp=rf[k].p;
				for(j=0;j<i1;j++) qval[rf[k].index[j]]=*(tmp++); 
			} else {
				for(j=0;j<n_all;j++) {
					a=a_set[kid][j];
					tmp=qval+j;
					while(a) {
						if(a&1) *tmp=1.0;
						a>>=1;
						tmp+=nb1;
					}
				}
			}
			pen(qval,kid,&mark->locus,n_all,n_bits,loki);
		}
		/* transmission probs */
		tp=id_array[kid].tp;
		cm[0]=mark->req_set[0][kid];
		cm[1]=mark->req_set[1][kid];
		i=id_array[ids].allele[X_MAT]-1;
		j=id_array[ids].allele[X_PAT]-1;
		k=id_array[idd].allele[X_MAT]-1;
		l=id_array[idd].allele[X_PAT]-1;
		if((1<<i)&cm[X_PAT]) i=n_all-1;
		if((1<<j)&cm[X_PAT]) j=n_all-1;
		if((1<<k)&cm[X_MAT]) k=n_all-1;
		if((1<<l)&cm[X_MAT]) l=n_all-1;
		i1=i<<n_bits;
		j2=j<<n_bits;
		p1=(pp[X_MM_PM]=tp[X_MM_PM]*qval[i1|k]);
		p1+=(pp[X_MM_PP]=tp[X_MM_PP]*qval[j2|k]);
		p1+=(pp[X_MP_PM]=tp[X_MP_PM]*qval[i1|l]);
		p1+=(pp[X_MP_PP]=tp[X_MP_PP]*qval[j2|l]);
		if(p1<=0.0) ABT_FUNC("Internal error - zero probability for child sample\n");
		do {
			z=ranf()*p1;
			p2=0.0;
			tmp=pp+3;
			n=4;
			while(n--) {
				z1=*(tmp--);
				if(z1>0.0) {
					p2+=z1;
					if(z<=p2) break;
				}
			}
		} while(n<0);
		switch(n) {
		 case X_MM_PM:
			id_array[kid].allele[X_MAT]=k+1;
			id_array[kid].allele[X_PAT]=i+1;
			break;
		 case X_MM_PP:
			id_array[kid].allele[X_MAT]=k+1;
			id_array[kid].allele[X_PAT]=j+1;
			break;
		 case X_MP_PM:
			id_array[kid].allele[X_MAT]=l+1;
			id_array[kid].allele[X_PAT]=i+1;
			break;
		 case X_MP_PP:
			id_array[kid].allele[X_MAT]=l+1;
			id_array[kid].allele[X_PAT]=j+1;
			break;
		 default:
			ABT_FUNC("Internal error - illegal sample\n");
		}
		id_array[kid].flag|=(SAMPLED_MAT|SAMPLED_PAT);
		assert(id_array[kid].allele[X_MAT]>0 && id_array[kid].allele[X_MAT]<=n_all);
		assert(id_array[kid].allele[X_PAT]>0 && id_array[kid].allele[X_PAT]<=n_all);
	}
	k=element->out_index;
#ifdef TRACE_PEEL
	if(CHK_PEEL(TRACE_LEVEL_2)) {
		(void) printf("Returning from %s() with %g\n",__func__,prob);
	}
#endif
	return prob;
}

