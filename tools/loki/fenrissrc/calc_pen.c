/****************************************************************************
 *                                                                          *
 *     Loki - Programs for genetic analysis of complex traits using MCMC    *
 *                                                                          *
 *                    Simon Heath - CNG, France                             *
 *                                                                          *
 *                         September 2002                                   *
 *                                                                          *
 * calc_pen.c:                                                              *
 *                                                                          *
 * Calculate penetrance vector for a component at a locus                   *
 *                                                                          *
 * Copyright (C) Simon C. Heath 2002                                        *
 * This is free software.  You can distribute it and/or modify it           *
 * under the terms of the Modified BSD license, see the file COPYING        *
 *                                                                          *
 ****************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#ifdef USE_DMALLOC
#include <dmalloc.h>
#endif
#include <math.h>

#include "utils.h"
#include "loki.h"
#include "loki_peel.h"
#include "loki_ibd.h"
#include "seg_pen.h"
#include "fenris_peel.h"
#include "fenris.h"
#include "calc_pen.h"

static int n_states,free_stack_ptr[ST_MAX_BITS],*free_stack[ST_MAX_BITS][ST_FREE_STACK_SIZE];
static struct cg_stack stack[256];
static struct int_store *store[ST_MAX_BITS],*root_store[ST_MAX_BITS];

static int *get_ints(int n) 
{
	int sz;
	int *p;
	struct int_store *pp,**pp1;

#ifdef DEBUG
	if(n<ST_MIN_BITS || n>=ST_MAX_BITS) ABT_FUNC("Illegal request size\n");
#endif
	n-=ST_MIN_BITS;
	if(free_stack_ptr[n]) {
		p=free_stack[n][--free_stack_ptr[n]];
		return p;
	}
	sz=(1<<(n+ST_MIN_BITS));
	pp1=store+n;
	pp=*pp1;
	if(pp) {
		if(pp->ptr>=ST_BLOCK_SIZE) {
			pp1=&pp->next;
			pp=*pp1;
			if(pp) pp->ptr=0;
		}
	}
	if(!pp) {
		if(!(pp=malloc(sizeof(struct int_store)))) ABT_FUNC(MMsg);
		*pp1=pp;
		if(!(pp->p=malloc(sizeof(int)*sz*ST_BLOCK_SIZE))) ABT_FUNC(MMsg);
		pp->ptr=0;
		pp->next=0;
		if(!root_store[n]) root_store[n]=store[n];
	}
	p=pp->p+sz*(pp->ptr++);
	store[n]=pp;
	return p;
}

static void free_ints(int n,int *p)
{
#ifdef DEBUG
	if(n<ST_MIN_BITS || n>=ST_MAX_BITS) ABT_FUNC("Illegal request size\n");
#endif
	n-=ST_MIN_BITS;
	if(free_stack_ptr[n]>=ST_FREE_STACK_SIZE) {
#ifdef DEBUG
		fprintf(stderr,"[%s:%d] %s(): Warning - stack limit reached, forgetting node (for this iteration only)\n",__FILE__,__LINE__,__func__);
#endif
	} else free_stack[n][free_stack_ptr[n]++]=p;
}

static int pass_genes(int *genes[2],int *seg[2],int *hap,int i,int par_flag)
{
	int par,fg=0;
	int g,j=0,k,k1,kid,stack_ptr=0;
	
	par=(par_flag==X_MAT)?id_array[i].dam-1:id_array[i].sire-1;
	g=genes[par_flag][i];
	genes[par_flag][i]=genes[seg[par_flag][i]][par];
	if(genes[par_flag][i]!=g) {
		if(hap[i]) fg=1;
		if(id_array[i].nkids) {
			for(;;) {
				k=id_array[i].sex==1?X_PAT:X_MAT;
				k1=0;
				for(;j<id_array[i].nkids;j++) {
					kid=id_array[i].kids[j];
					if(seg[k][kid]==par_flag) {
						g=genes[k][kid];
						genes[k][kid]=genes[par_flag][i];
						if(g!=genes[k][kid]) {
							if(hap[kid]) fg=1;
							if(id_array[kid].nkids) {
								stack[stack_ptr].id=i;
								stack[stack_ptr].par_flag=par_flag;
								stack[stack_ptr++].kid_ptr=j+1;
								i=kid;
								j=0;
								par_flag=k;
								k1=1;
								break;
							}
						}
					}
				}
				if(!k1) {
					if(stack_ptr) {
						j=stack[--stack_ptr].kid_ptr;
						par_flag=stack[stack_ptr].par_flag;
						i=stack[stack_ptr].id;
					} else break;
				}
			} 
		}
	}
	return fg;
}

void squash(int nobs,int nobs_fnd,int *obslist,int *genes[2],int *other,int *tt[2])
{
	int i,j,j1,k,k1,k2,*tt0,*tt1;
	
	tt0=tt[0];
	tt1=tt[1];
	k1=nobs_fnd*2;
	for(k=nobs_fnd;k<nobs;k++) {
		i=obslist[k];
		j=genes[0][i];
		j1=genes[1][i];
		if(tt0[j]<0 && tt0[j1]<0) {
			tt0[j]=tt1[j1]=k1++;
			tt1[j]=tt0[j1]=k1++;
			other[j]=j1;
			other[j1]=j;
		} else if(tt0[j]<0) {
			tt0[j]=tt1[j]=k1++;
			other[j]=j1;
			if(tt0[j1]!=tt1[j1]) {
				k2=(tt0[j1]>tt1[j1])?0:1;
				tt[k2][j1]=tt[k2^1][j1];
				tt[k2][other[j1]]=tt[k2^1][other[j1]];
			}
		} else if(tt0[j1]<0) {
			tt0[j1]=tt1[j1]=k1++;
			other[j1]=j1;
			if(tt0[j]!=tt1[j]) {
				k2=(tt0[j]>tt1[j])?0:1;
				tt[k2][j]=tt[k2^1][j];
				tt[k2][other[j]]=tt[k2^1][other[j]];
			}
		} else {
			if(other[j]!=j1) {
				if(tt0[j]!=tt1[j]) {
					k2=(tt0[j]>tt1[j])?0:1;
					tt[k2][j]=tt[k2^1][j];
					tt[k2][other[j]]=tt[k2^1][other[j]];
				}
				if(tt0[j1]!=tt1[j1]) {
					k2=(tt0[j1]>tt1[j1])?0:1;
					tt[k2][j1]=tt[k2^1][j1];
					tt[k2][other[j1]]=tt[k2^1][other[j1]];
				}
			}
		}
	}
}

static int find_node1(struct state_node1 *node,int *state,int i)
{
	int j,k,sz,*st,*p,*p1;
	
	st=node->state;
	k=state[i];
	sz=1<<node->size;
	for(j=0;j<node->n;j++) {
		if(st[j]==k) return st[j+sz];
	}
	node->n++;
	if(node->n>sz) {
		node->size++;
		p=get_ints(node->size+1);
		p1=p+sz*2;
		memcpy(p,st,sizeof(int)*j);
		memcpy(p1,st+sz,sizeof(int)*j);
		free_ints(node->size,st);
		node->state=p;
		p[j]=k;
		return p1[j]=n_states++;
	}
	st[j]=k;
	return st[j+sz]=n_states++;
}

static int find_node(struct state_node *node,int *state,int n,int i)
{
	int j,k,sz,*st;
	struct state_node1 *node1;
	
	st=node->state;
	k=state[i];
	sz=1<<node->size;
	for(j=0;j<node->n;j++) {
		if(st[j]==k) {
			if(i<n-2) return find_node(node->data.next+j,state,n,i+1);
			else return find_node1(node->data.last+j,state,i+1);
		}
	}
	node->n++;
	if(node->n>sz) {
		node->size++;
 		st=get_ints(node->size);
		memcpy(st,node->state,sizeof(int)*sz);
		free_ints(node->size-1,node->state);
		node->state=st;
		if(i<n-2) {
			if(!(node->data.next=realloc(node->data.next,sizeof(struct state_node)*sz*2))) ABT_FUNC(MMsg);
		} else {
			if(!(node->data.last=realloc(node->data.last,sizeof(struct state_node1)*sz*2))) ABT_FUNC(MMsg);
		}
	}
	st[j]=k;
	if(i<n-2) {
		node=node->data.next+j;
		i++;
		for(;i<n-2;i++) {
			node->size=INITIAL_NODE_BITS;
			node->state=get_ints(INITIAL_NODE_BITS);
			node->n=1;
			node->state[0]=state[i];
			if(!(node->data.next=malloc(sizeof(struct state_node)*INITIAL_NODE_SIZE))) ABT_FUNC(MMsg);
			node=node->data.next;
		}
		node->size=INITIAL_NODE_BITS;
		node->state=get_ints(INITIAL_NODE_BITS);
		node->n=1;
		node->state[0]=state[i];
		if(!(node->data.last=malloc(sizeof(struct state_node1)*INITIAL_NODE_SIZE))) ABT_FUNC(MMsg);
		j=0;
	}
	i++;
	node1=node->data.last+j;
	node1->state=get_ints(INITIAL_NODE_BITS+1);
	node1->n=1;
	node1->size=INITIAL_NODE_BITS;
	node1->state[0]=state[i];
	return node1->state[INITIAL_NODE_SIZE]=n_states++;
}

static int print_n(struct state_node *p,int n,int i)
{
	int j,k=0;
	struct state_node1 *p1;
	size_t t;
	
	t=i<n-2?sizeof(struct state_node):sizeof(struct state_node1);
	k+=(1<<p->size)*(sizeof(int)+t);
/*	k+=(p->n)*(sizeof(int)+t); */
	for(j=0;j<p->n;j++) {
		if(i<n-2) k+=print_n(p->data.next+j,n,i+1);
		else {
			p1=p->data.last+j;
			k+=sizeof(int)*(1<<(p1->size+1)); 
/*			k+=sizeof(int)*(2*p1->n); */
		}
	}
	return k;
}

void calc_pen(int comp,int n_bits,int locus,int pen_type,double *err_probs,double *p,int *fixed)
{
	int i,j,k,k1,k2,kk,*t,*map,nc,ct,iv,*genes[2],nn_all,n_all,*aflag,*a_trans,grp,*state,*pp;
	int *ff,cs,ng,nobs,nobs_fnd,nobs_fnd1,*obslist,ngen,*seg[2],*hap,*work,*work1[2];
	int *work2[2],g[2],nl,*n_longs,*inbr,*grps,*nstate[2],*alls[2];
	int chkstate[2][15]={{0,2,4,6,8,0,0,0,2,2,0,0,2,0,2},{1,3,5,7,9,2,2,2,4,6,7,8,2,7,7}};
	int chkstate1[2][15]={{0,2,4,6,8,0,0,0,2,2,0,0,2,0,2},{1,3,5,7,9,2,2,2,4,6,7,8,2,7,7}};
	int *obs_gps,*obs_genotypes,nobs1,nrepl;
	double **freq,*fq,z,**obspen,z1;
	unsigned long *tpl,*founders;
	struct pen_par ppar;
	struct state_node root;
	fenris_pen_func *pen;
	FILE *fptr;
	
	if(!(fptr=fopen("pen.dat","w"))) ABT_FUNC(MMsg);
	nc=1<<n_bits;
	printf("n_bits=%d, nc=0x%x\n",n_bits,nc);
	get_founder_params(&founders,&n_longs,0,&inbr);
	/* Set up allele frequencies */
	nn_all=marker[locus].locus.n_alleles;
	hap=marker[locus].haplo;
 	if(!(aflag=malloc(sizeof(int)*nn_all))) ABT_FUNC(MMsg);
	n_all=marker[locus].n_all1[comp];
	if(!(freq=malloc(sizeof(void *)*n_genetic_groups))) ABT_FUNC(MMsg);
	if(!(freq[0]=malloc(sizeof(double)*n_genetic_groups*n_all))) ABT_FUNC(MMsg);
	for(i=1;i<n_genetic_groups;i++) freq[i]=freq[i-1]+n_all;
	ngen=n_all*(n_all+1)/2;
	if(!(obs_genotypes=calloc((size_t)ngen,sizeof(int)))) ABT_FUNC(MMsg);
	z=1.0;
	a_trans=allele_trans[locus][comp];
	for(i=0;i<nn_all;i++) aflag[i]=0;
	for(j=0;j<n_all-1;j++) aflag[a_trans[j]]=1;
	for(grp=0;grp<n_genetic_groups;grp++) {
		fq=marker[locus].locus.freq[grp];
		for(j=0;j<n_all-1;j++) {
			k=a_trans[j];
			freq[grp][j]=fq[k];
			z-=fq[k];
		}
		freq[grp][j]=z;
	}
	free(aflag);
	/* Initialize node storage */
	for(i=0;i<ST_MAX_BITS;i++) {
		free_stack_ptr[i]=0;
		if((store[i]=root_store[i])) store[i]->ptr=0;
	}
	/* Initialize penetrance model */
	pen=get_pen_model(pen_type);
	ppar.e=err_probs;
	ppar.freq=freq;
	i=comp_start[comp];
	cs=comp_size[comp];
	/* Initialize founder_genes */
	if(!(genes[0]=malloc(sizeof(int)*cs*4))) ABT_FUNC(MMsg);
	genes[1]=genes[0]+cs;
	seg[0]=genes[1]+cs;
	seg[1]=seg[0]+cs;
	/* Set up penetrances for observed individuals */
	ff=founder_flag[locus];
	for(ng=nobs=nobs1=j=0;j<cs;j++,i++) {
		if(ff[i]) {
			genes[0][i]=ng++;
			if(id_array[i].nkids>1 || marker[locus].haplo[i]) genes[1][i]=ng++;
		} else {
			seg[0][j]=0;
			seg[1][j]=0;
			genes[X_MAT][j]=genes[0][id_array[i].dam-1];
			genes[X_PAT][j]=genes[0][id_array[i].sire-1];
		}
		if((k=marker[locus].haplo[i])) {
			k1=(k>>16)-1;
			k2=(k&65535)-1;
			k=k1*(k1+1)/2+k2;
			nobs++;
			if(!obs_genotypes[k]) {
				nobs1++;
				obs_genotypes[k]=-1;
			} else obs_genotypes[k]--;
		}
		id_array[i].flag=fixed?fixed[i]:0;
	}
	if(!(grps=malloc(sizeof(int)*ng))) ABT_FUNC(MMsg);
	setup_fseg_pen(ng,n_all);
	if(!nobs) ABT_FUNC("No observed individuals\n");
	if(!(obslist=malloc(sizeof(int)*nobs*7))) ABT_FUNC(MMsg);
	state=obslist+nobs;
	nstate[0]=state+nobs;
	nstate[1]=nstate[0]+nobs;
	obs_gps=nstate[1]+nobs;
	alls[0]=obs_gps+nobs;
	alls[1]=alls[0]+nobs;
	if(!(obspen=malloc(sizeof(void *)*nobs))) ABT_FUNC(MMsg);
	if(!(obspen[0]=malloc(sizeof(double)*nobs1*ngen))) ABT_FUNC(MMsg);
	for(i=0;i<nobs1*ngen;i++) obspen[0][i]=1.0;
/*	k1=-2;
	for(i=0;i<ngen;i++) if(obs_genotypes[i]<-1) {
		k=obs_genotypes[i];
		obs_genotypes[i]=k1;
		k1+=k;
	}
	for(i=0;i<ngen;i++) if(obs_genotypes[i]<0) {
		printf("%d %d\n",i,obs_genotypes[i]);
	}
	nrepl=-2-k1;
	i=comp_start[comp];
	if(nrepl) {
		for(j=0;j<cs;j++,i++) {
			if((k=marker[locus].haplo[i])) {
				k1=(k>>16)-1;
				k2=(k&65535)-1;
				k=k1*(k1+1)/2+k2;
				k1=obs_genotypes[k];
				if(k1<-1) {
					obs_gps[-2-k1]=j;
					obs_genotypes[k]--;
				}
			}
		}
	}
	for(i=0;i<nrepl;i++) {
		printf("%d\n",obs_gps[i]);
	}
	exit(0); */
	i=comp_start[comp];
	/* Add observed founders to observed list */
	for(nobs=nobs1=j=0;j<cs;j++,i++) {
		if(ff[i]) {
			if((k=marker[locus].haplo[i])) {
				k1=(k>>16)-1;
				k2=(k&65535)-1;
				alls[0][nobs]=k1;
				alls[1][nobs]=k2;
				k=k1*(k1+1)/2+k2;
				k1=obs_genotypes[k];
				if(k1<0) {
					obspen[nobs]=obspen[0]+(nobs1++)*ngen;
					obs_genotypes[k]=nobs;
					(void)pen(obspen[nobs],i,locus,&ppar);
				} else obspen[nobs]=obspen[k1];
				obslist[nobs++]=j;
			}
			k=id_array[i].kids[0];
			id_array[k].flag|=id_array[i].sex;
			k1=id_array[i].group-1;
			grps[genes[0][i]]=k1;
			if(id_array[i].nkids>1 || marker[locus].haplo[i]) grps[genes[1][i]]=k1;
		}
	}
	/* Add observed first children of founder couples to observed list */
	i=comp_start[comp];
	for(j=0;j<cs;j++,i++) {
		if(!ff[i] && (id_array[i].flag&3)==3 && marker[locus].haplo[i]) {
			k=marker[locus].haplo[i];
			k1=(k>>16)-1;
			k2=(k&65535)-1;
			alls[0][nobs]=k1;
			alls[1][nobs]=k2;
			k=k1*(k1+1)/2+k2;
			k1=obs_genotypes[k];
			if(k1<0) {
				obs_genotypes[k]=nobs;
				obspen[nobs]=obspen[0]+(nobs1++)*ngen;
				(void)pen(obspen[nobs],i,locus,&ppar);
			} else obspen[nobs]=obspen[k1];
			obslist[nobs++]=j;
			id_array[i].flag|=4;
		}
	}
	nobs_fnd=nobs;
	/* Add observed non-inbred individuals, not related to people already on observed list, to observed list */
	i=comp_start[comp];
	nl=n_longs[comp];
	tpl=founders+i*nl;
	for(j=0;j<cs;j++,i++) {
		if(!ff[i] && !(id_array[i].flag&4) && marker[locus].haplo[i]) {
			/* Check if inbred */
			if(!inbr[i]) {
				/* Check if related to people already on observed list */
				for(k1=0;k1<nobs;k1++) {
					k2=obslist[k1];
					for(k=0;k<nl;k++) if(tpl[j*nl+k]&tpl[k2*nl+k]) break;
					if(k<nl) break;
				}
				/* OK, add to list */
				if(k1==nobs) {
					k=marker[locus].haplo[i];
					k1=(k>>16)-1;
					k2=(k&65535)-1;
					alls[0][nobs]=k1;
					alls[1][nobs]=k2;
					k=k1*(k1+1)/2+k2;
					k1=obs_genotypes[k];
					if(k1<0) {
						obs_genotypes[k]=nobs;
						obspen[nobs]=obspen[0]+(nobs1++)*ngen;
						(void)pen(obspen[nobs],i,locus,&ppar);
					} else obspen[nobs]=obspen[k1];
					obslist[nobs++]=j;
					id_array[i].flag|=4;
				}
			}
		}
	}
	nobs_fnd1=nobs;
	/* Set up translation from bits to segregation indicators, and add remaining observed people 
	 * to observed list */
	if(!(map=malloc(sizeof(int)*n_bits))) ABT_FUNC(MMsg);
	i=comp_start[comp];
	kk=n_bits;
	for(j=0;j<cs;j++,i++) {
		if(!ff[i]) {
			if(!(id_array[i].flag&4) && marker[locus].haplo[i]) {
				k=marker[locus].haplo[i];
				k1=(k>>16)-1;
				k2=(k&65535)-1;
				alls[0][nobs]=k1;
				alls[1][nobs]=k2;
				k=k1*(k1+1)/2+k2;
				k1=obs_genotypes[k];
				if(k1<0) {
					obs_genotypes[k]=nobs;
					obspen[nobs]=obspen[0]+(nobs1++)*ngen;
					(void)pen(obspen[nobs],i,locus,&ppar);
				} else obspen[nobs]=obspen[k1];
				obslist[nobs++]=j;
			}
			k1=id_array[i].flag;
			if(!(k1&9)) map[--kk]=-1-j;
			if(!(k1&18)) map[--kk]=j+1;
		}
	}
	if(kk) ABT_FUNC("n_bits mismatch\n");
	printf("ng=%d, nobs=%d, nobs1=%d, nobs_fnd=%d, nobs_fnd1=%d\n",ng,nobs,nobs1,nobs_fnd,nobs_fnd1);
	if(!(work=malloc(sizeof(int)*ng*5))) ABT_FUNC(MMsg);
	work1[0]=work+ng;
	work1[1]=work1[0]+ng;
	work2[0]=work1[1]+ng;
	work2[1]=work2[0]+ng;
	root.n=0;
	root.size=INITIAL_NODE_BITS;
	root.state=get_ints(INITIAL_NODE_BITS);
	if(nobs-nobs_fnd>2) {
		if(!(root.data.next=malloc(sizeof(struct state_node)*INITIAL_NODE_SIZE))) ABT_FUNC(MMsg);
	} else if(!(root.data.last=malloc(sizeof(struct state_node1)*INITIAL_NODE_SIZE))) ABT_FUNC(MMsg);
	n_states=0;
	for(i=0;i<ng;i++) work2[0][i]=-1;
	squash(nobs_fnd,0,obslist,genes,work,work2);
	pp=work2[0];
	for(i=0;i<nobs_fnd;i++) {
		g[0]=genes[0][obslist[i]];
		g[1]=genes[1][obslist[i]];
		if(g[0]<g[1]) {
			nstate[0][i]=g[0];
			nstate[1][i]=g[1];
		} else {
			nstate[0][i]=g[1];
			nstate[1][i]=g[0];
		}
		g[0]=pp[g[0]];
		g[1]=pp[g[1]];
		state[i]=(g[0]>g[1]?g[0]*(g[0]+1)/2+g[1]:g[1]*(g[1]+1)/2+g[0]);
	}
/*	for(i=0;i<nobs;i++) {
		nstate[0][i]=chkstate[0][i];
		nstate[1][i]=chkstate[1][i];
	}
	z=fseg_pen(nstate,alls,obspen,grps,freq,nobs,ng,n_all);
	for(i=0;i<nobs;i++) printf("%d%d",nstate[0][i],nstate[1][i]);
	printf(" %.15f\n",z); 
   for(i=0;i<nobs;i++) {
		nstate[0][i]=chkstate1[0][i];
		nstate[1][i]=chkstate1[1][i];
	}
	z1=fseg_pen(nstate,alls,obspen,grps,freq,nobs,ng,n_all);
	for(i=0;i<nobs;i++) printf("%d%d",nstate[0][i],nstate[1][i]);
	printf(" %.15f %g\n",z1,z-z1);
	return; */
	/* Go through inheritance vectors in Gray code order */
	if(!(t=malloc(sizeof(int)*(n_bits+1)))) ABT_FUNC(MMsg);
	for(i=0;i<=n_bits;i++) t[i]=i;
	for(k=-1,j=iv=ct=0;ct<nc-1;ct++) {
		if(k) {
			memcpy(work1[0],work2[0],sizeof(int)*2*ng);
			squash(nobs,nobs_fnd,obslist,genes,work,work1);
			pp=work1[0];
			for(i=nobs_fnd1;i<nobs;i++) {
				g[0]=pp[genes[0][obslist[i]]];
				g[1]=pp[genes[1][obslist[i]]];
				state[i]=(g[0]>g[1]?g[0]*(g[0]+1)/2+g[1]:g[1]*(g[1]+1)/2+g[0]);
			}
			i=n_states;
			k=find_node(&root,state,nobs,nobs_fnd1);
			if(i!=n_states) {
				for(i=nobs_fnd;i<nobs;i++) {
					g[0]=genes[0][obslist[i]];
					g[1]=genes[1][obslist[i]];
					if(g[0]<g[1]) {
						nstate[0][i]=g[0];
						nstate[1][i]=g[1];
					} else {
						nstate[0][i]=g[1];
						nstate[1][i]=g[0];
					}
				}
/*				for(i=nobs_fnd;i<nobs;i++) printf("%d",state[i]);
				printf(" "); */
/*				for(i=0;i<nobs;i++) printf("%d%d",nstate[0][i],nstate[1][i]); */
				z=fseg_pen(nstate,alls,obspen,grps,freq,nobs,ng,n_all); 
				printf("%.10f\n",z); 
/*				fwrite(&z,sizeof(double),1,fptr); */
			}
		}
		i=t[0];
		j=map[i];
		if(j<0) {
			k=-1-j;
			k1=X_PAT;
		} else {
			k=j-1;
			k1=X_MAT;
		}
		seg[k1][k]^=1;
		k=pass_genes(genes,seg,hap,k,k1);
		t[0]=0;
		t[i]=t[i+1];
		t[i+1]=i+1;
		iv^=(1<<i);
	}
	if(k) {
		memcpy(work1[0],work2[0],sizeof(int)*2*ng);
		squash(nobs,nobs_fnd,obslist,genes,work,work1);
		pp=work1[0];
		for(i=nobs_fnd1;i<nobs;i++) {
			g[0]=pp[genes[X_MAT][obslist[i]]];
			g[1]=pp[genes[X_PAT][obslist[i]]];
			state[i]=(g[0]>g[1]?g[0]*(g[0]+1)/2+g[1]:g[1]*(g[1]+1)/2+g[0]);
		}
		i=n_states;
		k=find_node(&root,state,nobs,nobs_fnd1);
		if(i!=n_states) {
			for(i=nobs_fnd;i<nobs;i++) {
				g[0]=genes[0][obslist[i]];
				g[1]=genes[1][obslist[i]];
				if(g[0]<g[1]) {
					nstate[0][i]=g[0];
					nstate[1][i]=g[1];
				} else {
					nstate[0][i]=g[1];
					nstate[1][i]=g[0];
				}
			}
/*			for(i=nobs_fnd;i<nobs;i++) printf("%d",state[i]);
			printf(" "); 
			for(i=0;i<nobs;i++) printf("%d%d",nstate[0][i],nstate[1][i]); */
			z=fseg_pen(nstate,alls,obspen,grps,freq,nobs,ng,n_all);
			printf("%.10f\n",z); 
/*			fwrite(&z,sizeof(double),1,fptr); */
		}
	}
	for(i=0;i<ST_MAX_BITS;i++) {
		if(root_store[i]) {
			k=store[i]->ptr;
			store[i]=root_store[i];
			j=0;
			while(store[i]) {
				j++;
				store[i]=store[i]->next;
			}
			printf("%d %d %d\n",i,j,k);
		}
	}
	fclose(fptr);
	printf("n_states=%d\n",n_states);
	k=print_n(&root,nobs,nobs_fnd1);
	printf("size=%d\n",k);
	free(founders);
	free(grps);
	free(n_longs);
	free(work);
	free(genes[0]);
	free(map);
	free(obspen[0]);
	free(obspen);
	free(obslist);
	free(t);
	free(freq[0]);
	free(freq);
}

