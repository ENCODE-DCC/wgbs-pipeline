/****************************************************************************
 *                                                                          *
 *     Loki - Programs for genetic analysis of complex traits using MCMC    *
 *                                                                          *
 *                      Simon Heath - MSKCC                                 *
 *                                                                          *
 *                          August 2000                                     *
 *                                                                          *
 * handle_res.c:                                                            *
 *                                                                          *
 * Routines for handling residual variance                                  *
 *                                                                          *
 * Copyright (C) Simon C. Heath 1997, 2000, 2002                            *
 * This is free software.  You can distribute it and/or modify it           *
 * under the terms of the Modified BSD license, see the file COPYING        *
 *                                                                          *
 ****************************************************************************/

#include "config.h"
#include <math.h>
#include <stdio.h>
#include <float.h>
#ifndef DBL_MAX
#define DBL_MAX MAXDOUBLE
#endif

#include "utils.h"
#include "ranlib.h"
#include "loki.h"
#include "handle_res.h"

static void get_res_param(double *n,double *s,struct loki *loki)
{
	int i,j,nrec;
	double y;
	struct Id_Record *id_array;
	
	*s=*n=0.0;
	id_array=loki->pedigree->id_array;
	if(!loki->models->use_student_t) {
		for(i=0;i<loki->pedigree->ped_size;i++) if(id_array[i].res[0]) {
			nrec=id_array[i].n_rec;
			for(j=0;j<nrec;j++) {
				y=id_array[i].res[0][j];
				*s+=y*y;
			}
			*n+=(double)nrec;
		}
	} else {
		for(i=0;i<loki->pedigree->ped_size;i++) if(id_array[i].res[0]) {
			nrec=id_array[i].n_rec;
			for(j=0;j<nrec;j++) {
				y=id_array[i].res[0][j];
				*s+=y*y/id_array[i].vv[0][j];
			}
			*n+=(double)nrec;
		}
	}
}

double Calc_Res_Ratio(double v1,double v2,struct loki *loki)
{
	double s,nn,l1;

	get_res_param(&nn,&s,loki);
	l1=-.5*(nn*log(2.0*M_PI*v1)+s/v1);
	l1-=(-.5*(nn*log(2.0*M_PI*v2)+s/v2));
	return l1;
}

static double Calc_vprob(double v,double s,double x)
{
	v*=.5;
	return log(v)*v-lgamma(v)+log(s)*v-log(x)*(1.0+v)-v*s/x;
}

double Calc_Resprop(struct loki *loki)
{
	double s,v;

	get_res_param(&v,&s,loki);
	v+=RES_PRIOR_V0;
	s+=RES_PRIOR_V0*RES_PRIOR_S0;
	s/=v;
	return Calc_vprob(v,s,loki->models->residual_var[0]);
}

/* Sample residual variance */
double Sample_ResVar(struct loki *loki)
{
	double s,nn,y,v2;

	get_res_param(&nn,&s,loki);
	nn+=RES_PRIOR_V0;
	s+=RES_PRIOR_V0*RES_PRIOR_S0;
	y=loki->models->residual_var[0];
	loki->models->residual_var[0]=s/(sgamma(nn*0.5)*2.0);
	if(loki->models->residual_var[0]<loki->models->residual_var_limit[0]) loki->models->residual_var[0]=y;
	s/=nn;
	v2=nn*.5;
	y=Calc_vprob(v2,s,loki->models->residual_var[0]);
	return y;
}

double Calc_Var_Prior(double v,const struct loki *loki)
{
	return loki->sys.res_prior_konst-log(v)*(.5*RES_PRIOR_V0+1.0)-RES_PRIOR_S0*RES_PRIOR_V0/(2.0*v);
}

double Calc_ResLike(struct loki *loki)
{
	double s,nn;

	get_res_param(&nn,&s,loki);
	return -.5*(nn*log(2.0*M_PI*loki->models->residual_var[0])+s/loki->models->residual_var[0]);
}

double Calc_CensResLike(struct loki *loki)
{
	int i,j,type,idx;
	double s=0.0,nn=0.0,y,K=0.0,sd,z;
	struct Id_Record *id_array;
	
	id_array=loki->pedigree->id_array;
	sd=sqrt(loki->models->residual_var[0]*2.0);
	type=loki->models->models[0].var.type;
	idx=loki->models->models[0].var.var_index;
	if(type&ST_CONSTANT)	{
		for(i=0;i<loki->pedigree->ped_size;i++) if(id_array[i].res[0]) {
			y=id_array[i].res[0][0];
			if(!loki->models->censor_mode && (type&ST_CENSORED) && (id_array[i].data[idx].flag&2)) {
				z=.5*erfc((y-id_array[i].pseudo_qt[0][0])/sd);
				if(z>0.0) K+=log(z);
				else return -DBL_MAX;
			} else {
				s+=y*y;
				nn+=1.0;
			}
		}
	} else for(i=0;i<loki->pedigree->ped_size;i++) if(id_array[i].res[0])	{
		for(j=0;j<id_array[i].n_rec;j++)	{
			y=id_array[i].res[0][j];
			if(!loki->models->censor_mode && (id_array[i].data1[j][idx].flag&2))	{
				z=.5*erfc((y-id_array[i].pseudo_qt[0][j])/sd);
				if(z>0.0) K+=log(z);
				else return -DBL_MAX;
			} else {
				s+=y*y;
				nn+=1.0;
			}
		}
	}
	return K-.5*(nn*log(2.0*M_PI*loki->models->residual_var[0])+s/loki->models->residual_var[0]);
}

double Recalc_Res(int fg,struct loki *loki)
{
	int i,j,k,k1,k2,rec,nrec,type,mtype,idx,mod;
	struct id_data *data;
	double x,y,er=0.0;
	struct Id_Record *id_array;
	struct Locus *loc;
	struct Model *model;
	
	id_array=loki->pedigree->id_array;
	for(mod=0;mod<loki->models->n_models;mod++) {
		model=loki->models->models+mod;
		mtype=model->var.type;
		idx=model->var.var_index;
		for(i=0;i<loki->pedigree->ped_size;i++) if(id_array[i].res[mod]) {
			nrec=id_array[i].n_rec;
			for(rec=0;rec<nrec;rec++) {
				y=-loki->models->grand_mean[mod];
				if(mtype&ST_CONSTANT) {
					data=id_array[i].data+idx;
					if(data->flag&2) y+=id_array[i].pseudo_qt[mod][0];
				} else {
					data=id_array[i].data1[rec]+idx;
					if(data->flag&2) y+=id_array[i].pseudo_qt[mod][rec];
				}
				if(data->flag&ST_INTTYPE) y+=(double)data->data.value;
				else y+=data->data.rvalue;
				for(j=0;j<loki->params.n_tloci;j++) {
					loc=loki->models->tlocus+j;
					if(loc->flag && (loc->model_flag&(1<<mod))) {
						k2=loc->gt[i]-1;
						if(k2) y-=loc->eff[mod][k2-1];
					}
				}
				for(k=0;k<model->n_terms;k++) {
					type=model->term[k].vars[0].type;
					if(type&ST_TRAITLOCUS) continue;
					if(type&ST_ID) {
						y-=id_array[i].bv[mod];
						continue;
					}
					k1=model->term[k].vars[0].var_index;
					if(type&ST_MARKER) {
						k2=loki->markers->marker[k1].locus.gt[i]-1;
						if(k2) y-=model->term[k].eff[k2-1];
					} else {
						data=0;
						if(type&ST_CONSTANT) {
							if(id_array[i].data) data=id_array[i].data+k1;
						} else if(id_array[i].data1) data=id_array[i].data1[rec]+k1;
						if(!data) ABT_FUNC("Internal error - null data pointer\n");
						if(type&ST_FACTOR) {
							k2=(int)data->data.value-1;
							if(!(type&ST_RANDOM)) k2--;
							if(k2>=0) y-=model->term[k].eff[k2];
						} else {
							if(data->flag&ST_INTTYPE) x=(double)data->data.value;
							else x=data->data.rvalue;
							y-=model->term[k].eff[0]*x;
						}
					}
				}
				x=fabs(id_array[i].res[mod][rec]-y);
				er+=x;
				if(fg && x>0.0) (void)printf("%d %g %g\n",i,x,er);
				id_array[i].res[mod][rec]=y;
			}
		}
	}
	return er;
}
