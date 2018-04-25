/****************************************************************************
 *                                                                          *
 *     Loki - Programs for genetic analysis of complex traits using MCMC    *
 *                                                                          *
 *                 Simon Heath - CNG, Paris                                 *
 *                                                                          *
 *                       August 2002                                        *
 *                                                                          *
 * cleanup.c:                                                               *
 *                                                                          *
 * Clean up and free memory                                                 *
 *                                                                          *
 * Copyright (C) Simon C. Heath 1997, 2000, 2002                            *
 * This is free software.  You can distribute it and/or modify it           *
 * under the terms of the Modified BSD license, see the file COPYING        *
 *                                                                          *
 ****************************************************************************/

#include <config.h>
#include <stdlib.h>
#ifdef USE_DMALLOC
#include <dmalloc.h>
#endif
#include <stdio.h>
#include <signal.h>

#include "version.h"
#include "libhdr.h"
#include "ranlib.h"
#include "utils.h"
#include "scan.h"
#include "y.tab.h"

int catch_sigs,sig_caught,sig_quiet;
struct remember *RemBlock;
static struct remember *FirstRemBlock;
static char **LFile;

void free_infile(struct InFile *infile)
{
	struct DataBlock *db,*db1;
	
	db=infile->data;
	while(db) {
		if(db->records) free(db->records);
		if(db->type) free(db->type);
		db1=db->next;
		free(db);
		db=db1;
	}
	if(infile->format) {
		free(infile->format->f_atoms);
		free(infile->format);
	}
	if(infile->fformat) {
		if(infile->fformat->rs) free(infile->fformat->rs);
		if(infile->fformat->fs) free(infile->fformat->fs);
		if(infile->fformat->gs) free(infile->fformat->gs);
		free(infile->fformat);
	}
	free(infile->name);
	if(infile->element) free(infile->element);
	free(infile);
}

static void free_miss(struct Miss *m)
{
	if(m->Missing.type==ST_STRING) free(m->Missing.arg.string);
	if(m->element) free(m->element);
	if(m->scope) free(m->scope);
	free(m);
}

static void free_link(struct Link *l)
{
	if(l->name) free(l->name);
	if(l->element) free(l->element);
	free(l);
}

static void free_model(struct model *m)
{
	struct model_list *ml,*ml1;
	
	ml=m->model_list;
	while(ml) {
		ml1=ml->next;
		if(ml->element) free(ml->element);
		free(ml);
		ml=ml1;
	}
	free(m);
}

void free_op(struct operation *op)
{
	struct operation *op1;
	while(op) {
		if(op->type==STRING) free(op->arg.string);
		op1=op->next;
		free(op);
		op=op1;
	}
}

void free_restrict(struct Restrict *res)
{
	if(res->element) free(res->element);
	free_op(res->Op_List);
	free(res);
}

static void free_node(void *s)
{
	struct scan_data *sd;
	
	sd=s;
	if(sd->name) free(sd->name);
	if(sd->element) free(sd->element);
	free(s);
}

/*   This is called when the program finishes without calling abt(),
 *   i.e., normal termination or SIGINT etc.
 */
static void FreeStuff(void)
{
	int i;
	struct InFile *infile,*infile1;
	struct Miss *m,*m1;
	struct Link *l,*l1;
	struct model *model,*model1;
	struct Restrict *res,*res1;
	struct Censor *cen,*cen1;
	
	if(from_abt) return;
	if(id_array) {
		for(i=0;i<ped_size;i++) if(id_array[i].haplo[0]) free(id_array[i].haplo[0]);
		free(id_array[0].mkflag);
		free(id_array);
	}
	if(id_elements) free(id_elements);
	if(ped_recode) free(ped_recode);
	if(family_recode) free(family_recode);
	if(group_recode) free(group_recode);
	if(ped_recode1) free(ped_recode1);
	if(rec_tab) free(rec_tab);
	if(rec_tab1) free(rec_tab1);
	if(rsformat) free(rsformat);
	if(fsformat) free(fsformat);
	if(gsformat) free(gsformat);
	if(factor_recode) {
		for(i=0;i<(n_factors+n_markers);i++) if(factor_recode[i]) free(factor_recode[i]);
		free(factor_recode);
	}
	if(var_factors) free(var_factors);
	if(comp_size) free(comp_size);
	infile=Infiles;
	while(infile) {
		infile1=infile->next;
		free_infile(infile);
		infile=infile1;
	}
	m=Miss;
	while(m)	{
		m1=m->next;
		free_miss(m);
		m=m1;
	}
	l=links;
	while(l)	{
		l1=l->next;
		free_link(l);
		l=l1;
	}
	model=Models;
	while(model) {
		model1=model->next;
		free_model(model);
		model=model1;
	}
	res=Restrictions;
	while(res) {
		res1=res->next;
		free_restrict(res);
		res=res1;
	}
	cen=Censored;
	while(cen) {
		cen1=cen->next;
		free_op(cen->Op_List);
		free(cen);
		cen=cen1;
	}
	if(Affected) free_op(Affected);
	if(Unaffected) free_op(Unaffected);
	if(Proband) free_op(Proband);
	if(markers) {
		for(i=0;i<n_markers;i++) {
			if(markers[i].order)	{
				if(markers[i].order[0]) free(markers[i].order[0]);
				free(markers[i].order);
			}
			if(markers[i].allele_trans) {
				if(markers[i].allele_trans[0]) free(markers[i].allele_trans[0]);
				free(markers[i].allele_trans);
			}
			if(markers[i].o_size) free(markers[i].o_size);
			if(markers[i].sub_element) free(markers[i].sub_element);
		}
		free(markers);
	}
	if(traitlocus) {
		if(traitlocus->order) {	
			if(traitlocus->order[0]) free(traitlocus->order[0]);
			free(traitlocus->order);
		}
		if(traitlocus->o_size) free(traitlocus->o_size);
		free(traitlocus);
	}
	if(root_var) free_bin_tree(root_var,free_node);
	if(Filter) free(Filter);
	if(*LFile) free(*LFile);
	if(OutputFile) free(OutputFile);
	if(OutputLaurFile) free(OutputLaurFile);
	if(OutputRawFile) free(OutputRawFile);
	if(ErrorDir) free(ErrorDir);
	FreeRemem(FirstRemBlock);
}

/* Called when program receives signal */
static void int_handler(int i)
{
	static int oldsig=-1;
	char *sigs[]={"SIGTERM","SIGQUIT","SIGHUP","SIGINT","UNKNOWN"};
	int j,signals[]={SIGTERM,SIGQUIT,SIGHUP,SIGINT};
	
	for(j=0;j<4;j++) if(i==signals[j]) break;
	if(catch_sigs && oldsig!=i)	{
		sig_caught=oldsig=i;
		if(!sig_quiet) (void)printf("Caught sig %s\n",sigs[j]);
		return;
	}
	if(!sig_quiet) (void)fprintf(stderr,"Exiting on signal %s (%d)\n",sigs[j],i);
	exit(i);
}

void init_stuff(char **p)
{
	struct sigaction s_action;

	LFile=p;
	if(!(FirstRemBlock=malloc(sizeof(struct remember)))) ABT_FUNC(MMsg);
	RemBlock=FirstRemBlock;
	RemBlock->pos=0;
	RemBlock->next=0;
	/* If program terminates normally or via a signal I want it to go through
	 * FreeStuff() so that Memory can be cleared so I can check for unfreed blocks
	 */
	if(atexit(FreeStuff)) ABT_FUNC("Unable to register exit function FreeStuff()\n");
	if(atexit(print_end_time)) ABT_FUNC("Unable to register exit function print_end_time()\n");
	s_action.sa_handler=int_handler;
	s_action.sa_flags=0;
	(void)sigaction(SIGINT,&s_action,0L);
	(void)sigaction(SIGHUP,&s_action,0L);
	(void)sigaction(SIGQUIT,&s_action,0L);
	(void)sigaction(SIGTERM,&s_action,0L);
}
