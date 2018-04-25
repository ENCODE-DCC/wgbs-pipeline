#include <config.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>
#ifdef USE_DMALLOC
#include <dmalloc.h>
#endif
#include <stdio.h>
#include <math.h>
#include <ctype.h>
#include <assert.h>

#include "prep.h"
#include "parser.h"
#include "parse.tab.h"
#include "prep_input.h"

static parse_handler handle;
static struct loki *loki;
static char *keys[]={"USE","WHERE","LOCUS","CENSORED","START",0};
static char *special_name[N_SPECIAL_VAR]={"FS","RS","GS","SKIP"};
static struct special_vars special_vars;
static struct typeI_wait typeI_wait;

static void free_lk_assign(void *s)
{
	struct lk_assign *p;
	
	p=s;
	if(p->var) free_parse_term(p->var);
	if(p->arg) free_parse_term(p->arg);
	p->var=p->arg=0;
}

static void free_missing(void *s)
{
	struct Missing *m;
	
	m=s;
	if(m->miss) free_parse_term(m->miss);
	free_parse_var_list(m->vl);
	if(m->scope) free(m->scope);
}

static void free_files(void *s)
{
	struct prep_file *files;
	
	files=s;
	
	free(files->name);
	free_parse_var_list(files->varlist);
	if(files->fixed_start) free(files->fixed_start);
	if(files->fs) free_string(files->fs);
	if(files->rs) free_string(files->rs);
	if(files->gs) free_string(files->gs);
}

static void free_markers(void *s)
{
	struct lk_marker *mk;
	int i;
	
	mk=s;
	free_parse_term(mk->term);
	for(i=0;i<2;i++) if(mk->haps[i]) free_parse_term(mk->haps[i]);
}

static void free_sex_def(void *s)
{
	struct sex_define *sd;
	
	sd=s;
	free_parse_var_list(sd->sex_exp);
}

static void free_link(void *s)
{
	struct lk_link *lk;
	
	lk=s;
	if(lk->name) free_string(lk->name);
}

static void FreeReadControl(void)
{
	int i;
	struct lk_data *dat;
	
	for(i=0;i<LK_NUM_NAMES;i++) {
		if(loki->names[i]) {
			free(loki->names[i]);
			loki->names[i]=0;
		}
	}
	dat=loki->prep->data;
	free_list(dat->missing,free_missing);
	dat->missing=0;
	free_list(dat->infiles,free_files);
	dat->infiles=0;
	if(loki->prep->pedigree->ped_vars) {
		free_parse_var_list(loki->prep->pedigree->ped_vars);
		loki->prep->pedigree->ped_vars=0;
	}
	free_list(loki->prep->links,free_link);
	loki->prep->links=0;
	free_list(loki->prep->pedigree->sex_def,free_sex_def);
	loki->prep->pedigree->sex_def=0;
	free_list(loki->prep->markers,free_markers);
	loki->prep->markers=0;
	free_list(loki->prep->models,free_lk_assign);
	free_list(loki->prep->assignments,free_lk_assign);
	get_term_list(0,0);
}

static void model_com(struct parse_term *v,struct parse_term *ex)
{
	struct lk_assign *as;
	int type;
	
	assert(v && ex);
	type=ex->type&VT_TYPES;
	if(type!=VARIABLE && type!=EXPR_OP && type != FUNCTION) {
		parerror2("MODEL expression is constant\n");
	} else {
		fputs("MODEL ",stdout);
		print_parse_term_name(stdout,v);
		fputs(" = ",stdout);
		print_parse_exp(stdout,ex);
		fputc('\n',stdout);
		if(!(as=malloc(sizeof(struct lk_assign)))) ABT_FUNC(MMsg);
		as->next=loki->prep->models;
		loki->prep->models=as;
		as->var=copy_parse_term(v);
		as->arg=copy_parse_term(ex);
	}
}

static void assign(struct parse_term *v,struct parse_term *ex)
{
	struct lk_assign *as;
	int type;
	
	assert(v && ex);
	print_parse_term_name(stdout,v);
	fputs(" = ",stdout);
	print_parse_exp(stdout,ex);
	fputc('\n',stdout);
	type=ex->type&VT_TYPES;
	if(type!=VARIABLE && type!=EXPR_OP && type!=FUNCTION) return;
	if(!(as=malloc(sizeof(struct lk_assign)))) ABT_FUNC(MMsg);
	as->next=loki->prep->assignments;
	loki->prep->assignments=as;
	as->var=copy_parse_term(v);
	as->arg=copy_parse_term(ex);
}

static void bad_clause(char *com,struct parse_clause *c)
{
	string *s;
	
	if(c) {
		s=parse_clause_str(c);
		if(s) {
			parerror2("Unrecognized clause [%s] for %s command\n",get_cstring(s),com);
			free_string(s);
		}
	}
}

static void bad_arg(char *com,struct parse_term *t)
{
	string *s;
	
	if(t) {
 		s=parse_exp_str(t);
		if(s) {
			parerror2("Unused argument %s to %s command\n",get_cstring(s),com);
			free_string(s);
		}
	}
}

static void log_com(struct parse_clause *clause,struct parse_term *s,struct parse_var_list *vl)
{
	string *s1=0;
	
	if(clause) bad_clause("LOG",clause);
	if(vl) {
		s1=var_list_string(vl);
		parerror2("Unused variable list '%s' for LOG command\n",get_cstring(s1));
		free_string(s1);
		free_parse_var_list(vl);
	}
	if(s) s1=convert_to_string(s);
	if(!s1) parerror2("No logfile specified for LOG command\n");
	else {
		printf("LOG '%s'\n",get_cstring(s1));
		if(loki->names[LK_LOGFILE]) free(loki->names[LK_LOGFILE]);
		loki->names[LK_LOGFILE]=extract_cstring(s1);
	}
}

static void missing_com(struct parse_clause *clause,struct parse_term *s,struct parse_var_list *vl)
{
	int i;
	string *s1=0,*s2=0;
	char *p;
	struct Missing *m;
	
	if(clause && vl) parerror2("Can't have both explicit and implicit scope for MISSING command\n");
	else {
		if(!s) {
			parerror2("No argument to MISSING command\n");
			free_parse_var_list(vl);
		} else {
			if(!(m=malloc(sizeof(struct Missing)))) ABT_FUNC(MMsg);
			m->next=loki->prep->data->missing;
			loki->prep->data->missing=m;
			m->miss=copy_parse_term(s);
			m->vl=vl;
			m->scope=0;
			printf("MISSING ");
			if(clause) {
				do {
					if(clause->type==PC_TERM && clause->clause.term->type==STRING) s1=add_strings(s1,copy_string(clause->clause.term->elem.str));
					else {
						if(s2) s2=add_to_string(s2,';');
						if(clause->type==PC_TERM) s2=add_strings(s2,parse_exp_str(clause->clause.term));
						else s2=add_strings(s2,format_clause_str(clause->clause.format));
					}
					clause=clause->next;
				} while(clause);
				if(s1) {
					p=extract_cstring(s1);
					qstrip(p);
					m->scope=p;
					i=0;
					while(*p) {
						switch(toupper((int)*p)) {
						 case '!':
						 case 'F':
						 case 'G':
						 case 'P':
						 case 'C':
						 case 'R':
						 case 'I':
							break;
						 default: i=1;
						}
						if(i) break;
						p++;
					}
					i=0;
					if(*p) {
						i=1;
						parerror2("Illegal character '%c' in MISSING scope\n",*p);
					} else if(*(--p)=='!') {
						i=1;
						parerror2("MISSING scope can not end with a '!'\n");
					}
					if(!i) printf("[\"%s\"] ",m->scope);
					else {
						free(m->scope);
						m->scope=0;
					}
				}
			}
			print_parse_exp(stdout,s);
			if(vl) {
				fputc(' ',stdout);
				print_var_list(stdout,vl);
			}
			fputc('\n',stdout);
			if(s2) {
				parerror2("Unrecognized subclause(s) '%s' for MISSING statement\n",get_cstring(s2));
				free_string(s2);
			}
		}
	}
}

static void pedigree_com(struct parse_clause *clause,struct parse_term *s,struct parse_var_list *vl)
{
	if(clause) bad_clause("PEDIGREE",clause);
	if(s) bad_arg("PEDIGREE",s);
	if(!vl) parerror2("Missing variable list to PEDIGREE command\n");
	else {
		if(loki->prep->pedigree->ped_vars) {
			parerror2("Duplicate PEDIGREE command (ignored)\n");
			free_parse_var_list(vl);
		} else {
			loki->prep->pedigree->ped_vars=vl;
			printf("PEDIGREE ");
			print_var_list(stdout,vl);
			fputc('\n',stdout);
		}
	}
}

static int check_ped(void)
{
	struct parse_var_list *vl;
	struct parse_term *v;
	struct parse_var *vv;
	int err=0,i,k,ct,fam_flag=0;
	
	if((vl=loki->prep->pedigree->ped_vars)) {
		ct=0;
		while(vl) {
			v=&vl->term;
			k=1;
			if((v->type&VT_TYPES)==VARIABLE) {
				if(!(v->type&VT_ARRAY_ELEMENT)) {
					vv=v->elem.vv;
					if(vv->type&VT_ARRAY) {
						if(vv->type&VT_SIZE_SET) k=vv->size;
						else {
							fprintf(stderr,"Error in PEDIGREE command - size of array '%s' is undetermined\n",get_cstring(vv->name));
							err=1;
						}
					}
				}
			}
			ct+=k;
			vl=vl->next;
		}
		if(!err) {
			if(ct==3) fam_flag=0;
			else if(ct==4) fam_flag=1;
			else {
				fprintf(stderr,"Wrong number of variables (%d) specified for PEDIGREE command\n",ct);
				fputs("Only 3 (without family code) or 4 (with family code) variables allowed\n",stderr);
				err=1;
			}
		}
		if(!err) {
			vl=loki->prep->pedigree->ped_vars;
			ct=0;
			while(vl) {
				v=&vl->term;
				k=1;
				if((v->type&VT_TYPES)==VARIABLE) {
					if(v->type&VT_ARRAY_ELEMENT) {
						vv=get_array_var(v->elem.vv1->head,&v->elem.vv1->ex);
					} else {
						vv=v->elem.vv;
						if(vv->type&VT_ARRAY) k=vv->size;
					}
					for(i=0;i<k;i++) {
						switch(i+ct-fam_flag) {
						 case -1:
							vv[i].ltype|=(ST_FAMILY|ST_FACTOR|ST_CONSTANT|ST_REQUIRED);
							loki->prep->pedigree->fam_flag=1;
							break;
						 case 0:
							vv[i].ltype|=(ST_ID|ST_FACTOR|ST_CONSTANT|ST_REQUIRED);
							break;
						 case 1:
							vv[i].ltype|=(ST_SIRE|ST_FACTOR|ST_CONSTANT|ST_REQUIRED);
							break;
						 case 2:
							vv[i].ltype|=(ST_DAM|ST_FACTOR|ST_CONSTANT|ST_REQUIRED);
							break;
						}
					}
				} else {
					if(ct==fam_flag) {
						fputs("Error in PEDIGREE command - missing ID variable\n",stderr);
						err=1;
						break;
					}
				}
				ct+=k;
				vl=vl->next;
			}
		}
	} else {
		fputs("No PEDIGREE command found - can not construct pedigree\n",stderr);
		err=1;
	}
	return err;
}

static void frequency_com(int key _U_,struct parse_term *arg _U_,struct parse_var_list *vl _U_)
{
}

static void sex_com(int key,struct parse_term *arg,struct parse_var_list *vl)
{
	int i,fg=0;
	string *s;
	struct sex_define *sd;
	
	if(key>=0) {
		parerror2("Syntax error: can't use %s in front of SEX command\n",keys[key]);
	} else {
		if(!arg || !vl) parerror2("Error: missing arguments for SEX command\n");
		else {
			i=list_length(vl);
			if(i!=2) {
				s=var_list_string(vl);
				if(s) {
					parerror2("Error: SEX command has %d parameters [%s] (2 expected)\n",i,get_cstring(s));
					free_string(s);
				} else parerror2("Error: SEX command has %d parameters (2 expected)\n",i);
			} else {
				if(!(sd=malloc(sizeof(struct sex_define)))) ABT_FUNC(MMsg);
				sd->next=loki->prep->pedigree->sex_def;
				loki->prep->pedigree->sex_def=sd;
				sd->sex_exp=vl;
				fg=1;
				printf("SEX ");
				print_parse_exp(stdout,arg);
				fputc(' ',stdout);
				print_var_list(stdout,vl);
				fputc('\n',stdout);
			}
		}
	}
	if(!fg && vl) free_parse_var_list(vl);
}

static void marker_locus_com(int key,struct parse_var_list *vl,int com,struct lk_link *lk)
{
	int i,j,ix,type,er,ct=0;
	struct parse_var_list *vl1;
	struct lk_marker *mk;
	struct parse_clause *pc;
	struct parse_term *tm,*tl,*haps[2];
	struct parse_var *vv,*vv1,*vv2,*vhaps[2];
	string *s;
	char *p,*p1=0,*p2,*nms[]={"MARKER","TRAIT"};
	
	if(com<0 && com>1) ABT_FUNC("Internal error - invalid command\n");
	if(key!=2) parerror2("Syntax error: Unknown command %s %s\n",nms[com],keys[key]);
	else {
		if(lk) {
			if(lk->name) printf("LINK \"%s\"",get_cstring(lk->name));
			else fputs("LINK",stdout);
		} else printf("%s Locus",nms[com]);
		vl1=vl;
		while(vl1) {
			ct++;
			type=vl1->term.type&VT_TYPES;
			if(type==VARIABLE) {
				fputc(vl==vl1?' ':',',stdout);
				s=parse_exp_str(&vl1->term);
				fputs(get_cstring(s),stdout);
				free_string(s);
				ix=0;
				if((pc=vl1->clause)) {
					s=add_to_string(0,'[');
					s=add_strings(s,parse_clause_str(vl1->clause));
					s=add_to_string(s,']');
					fputs(get_cstring(s),stdout);
					free_string(s);
					while(pc) {
						if(!com && pc->type==PC_TERM) {
							tm=pc->clause.term;
							i=0;
							tl=get_term_list(tm,&i);
							for(j=0;j<i;j++) {
								type=tl[j].type&VT_TYPES;
								er=0;
								if(type==VARIABLE) {
									vv=find_parse_var(tl+j);
									if(!vv) er=-1;
									else if(vv->type&VT_IMMED) er=1;
									else if(ix==2) er=2;
									else {
										vhaps[ix]=vv;
										haps[ix++]=copy_parse_term(tl+j);
									}
								} else if(type!=EMPTY_VAR) er=3;
								if(er>0) {
									if(!p1) p1=parse_exp_cstr(&vl1->term);
									p=parse_exp_cstr(tl+j);
									switch(er) {
									 case 1:
										parerror2("Invalid use of non-vector variable %s for Marker %s\n",SAFE_PTR(p),SAFE_PTR(p1));
										break;
									 case 2:
										parerror2("Extra haplotype vector %s listed for Marker %s\n",SAFE_PTR(p),SAFE_PTR(p1));
										break;
									 case 3:
										parerror2("Expression %s inappropriate as a haplotype vector for Marker %s\n",SAFE_PTR(p),SAFE_PTR(p1));
										break;
									}
									if(p) free(p);
								}
							}
						} else {
							if(!p1) p1=parse_exp_cstr(&vl1->term);
							if(pc->type==PC_TERM) s=parse_exp_str(pc->clause.term);
							else s=format_clause_str(pc->clause.format);
							parerror2("Invalid clause [%s] for Locus %s\n",s?get_cstring(s):NULL_STR,SAFE_PTR(p1));
							if(s) free_string(s);
						}
						pc=pc->next;
					}
				}
				er=0;
				vv1=find_parse_var(&vl1->term);
				if(ix) {
					for(i=0;i<ix;i++) {
						if(vhaps[i]->ltype&ST_HAPLO) break;
						vhaps[i]->ltype|=ST_HAPLO;
					}
					if(i<ix) {
						mk=loki->prep->markers;
						while(mk) {
							for(i=0;i<2;i++) {
								tm=mk->haps[i];
								if(!tm) break;
								vv=find_parse_var(tm);
								for(j=0;j<ix;j++) {
									if(vv==vhaps[j]) {
										vv2=find_parse_var(mk->term);
										if(vv2!=vv1) break;
									}
								}
								if(j<ix) {
									if(!p1) p1=parse_exp_cstr(&vl1->term);
									p2=parse_exp_cstr(mk->term);
									p=parse_exp_cstr(tm);
									parerror2("Duplicate usage of haplotype vector %s for Markers %s and %s\n",SAFE_PTR(p),SAFE_PTR(p1),SAFE_PTR(p2));
									if(p) free(p);
									if(p2) free(p2);
									er=2;
								}
							}
							mk=mk->next;
						}
					} 
				}
				if(er) {
					for(i=0;i<ix;i++) free_parse_term(haps[i]);
				} else {
					if(vv1->ltype&(ST_MARKER|ST_LINKED)) {
						mk=loki->prep->markers;
						while(mk) {
							vv2=find_parse_var(mk->term);
							if(vv2==vv1) break;
							mk=mk->next;
						}
						if(!mk) ABT_FUNC("Internal error - marker not found\n");
						if(mk->type!=(com?LK_QTL:LK_MARKER)) {
 							if(!p1) p1=parse_exp_cstr(&vl1->term);
							parerror2("Locus %s used as both Trait and Marker locus\n",SAFE_PTR(p1));
							for(i=0;i<ix;i++) free_parse_term(haps[i]);
						} else {
							for(i=0;i<ix;i++) {
								for(j=0;j<2;j++) {
									if(!mk->haps[j]) break;
									vv=find_parse_var(mk->haps[j]);
									if(vv==vhaps[i]) break;
								}
								if(j==2) break;
								if(!mk->haps[j]) mk->haps[j]=haps[i];
								else free_parse_term(haps[i]);
							}
							if(i<ix) parerror2("Mismatch in haplotype vectors for Marker locus\n");
						}
						if(lk) {
							if(!mk->link) mk->link=lk;
							else if(mk->link!=lk) {
								if(!p1) p1=parse_exp_cstr(&vl1->term);
								parerror2("Marker locus %s appears in multiple linkage groups\n",SAFE_PTR(p1));
							}
						}
					} else {
						if(!(mk=malloc(sizeof(struct lk_marker)))) ABT_FUNC(MMsg);
						mk->next=loki->prep->markers;
						loki->prep->markers=mk;
						mk->term=copy_parse_term(&vl1->term);
						mk->link=lk;
						for(i=0;i<ix;i++) mk->haps[i]=haps[i];
						for(;i<2;i++) mk->haps[i]=0;
						mk->type=com?LK_QTL:LK_MARKER;
					}
				}
				vv1->ltype|=(lk?ST_LINKED:ST_MARKER);
			} else if(type!=EMPTY_VAR) {
				parerror2("Illegal use of immediate variable (no. %d in list) in LOCUS statement\n",ct);
			}
			if(p1) {
				free(p1);
				p1=0;
			}
			vl1=vl1->next;
		}
		fputc('\n',stdout);
	}
	free_parse_var_list(vl);
}

static void marker_com(struct parse_clause *clause,struct parse_term *s,struct parse_var_list *vl)
{
	if(clause) bad_clause("MARKER",clause);
	if(s) bad_arg("MARKER",s);
	marker_locus_com(2,vl,0,0);
}

static void file_com(struct parse_clause *pc,struct parse_term *s,struct parse_var_list *vl)
{
	int er=0;
	struct prep_file *file=0;
	
	if(!s) {
		parerror2("Missing file name for FILE command\n");
		er=1;
	}
	if(!er) {
		file=prep_file_com(pc,s,vl,&special_vars,loki);
	}
	if(vl) {
		if(er) free_parse_var_list(vl);
	}
	if(!vl) {
		typeI_wait.type=WAIT_FILE;
		typeI_wait.ptr.file=file;
		set_typeI_kludge();
	}
}

static void array_com(struct parse_clause *pc,struct parse_term *s,struct parse_var_list *vl)
{
	struct parse_var_list *vl1;
	struct parse_var *head;
	struct parse_clause *pc1;
	string *s1;
	char *p1;
	int type,er,size,ct=0;
	
	if(pc) bad_clause("ARRAY",pc);
	if(s) bad_arg("ARRAY",s);
	vl1=vl;
	while(vl1) {
		ct++;
		type=vl1->term.type&VT_TYPES;
		if(type==VARIABLE) {
			if((pc1=vl1->clause)) {
				p1=parse_exp_cstr(&vl1->term);
				if(pc1->type==PC_TERM) s1=parse_exp_str(pc1->clause.term);
				else s1=format_clause_str(pc1->clause.format);
				parerror2("Invalid clause [%s] in ARRAY command for variable %s\n",s1?get_cstring(s1):NULL_STR,SAFE_PTR(p1));
				if(p1) free(p1);
				if(s1) free_string(s1);
			} 
			/* Size specified ? */
			if(vl1->term.type&VT_ARRAY_ELEMENT) {
				head=vl1->term.elem.vv1->head;
				vl1->term.elem.vv1->ex=try_resolve_expr(vl1->term.elem.vv1->ex);
				size=convert_to_int(vl1->term.elem.vv1->ex,&er);
				if(er) ABT_FUNC("Internal error\n");
				if(size<head->size) parerror2("Illegal attempt to reduce array size for variable %s\n",get_cstring(head->name));
				else {
					if((head->type&VT_SIZE_SET) && size!=head->size) parerror2("Illegal attempt to change array size for variable %s\n",get_cstring(head->name));
					else head->type|=VT_SIZE_SET;
				}
			} else {
				head=vl1->term.elem.vv;
				if(head->type&VT_SCALAR) parerror2("Illegal attempt to change variable %s from a scalar to an array\n",get_cstring(head->name));
				else head->type|=VT_ARRAY;
			}
		} else if(type!=EMPTY_VAR) {
			parerror2("Illegal use of immediate variable (no. %d in list) in ARRAY statement\n",ct);
		}
		vl1=vl1->next;
	}
	free_parse_var_list(vl);
}

static void link_com(struct parse_clause *pc,struct parse_term *s,struct parse_var_list *vl)
{
	int er=0,type=0;
	string *s1=0;
	struct lk_link *lk;

	if(pc) {
		if(pc->next) er=1;
		else s1=parse_clause_str(pc);
		if(!s1) er=1;
		else {
			strip_quotes(s1);
			if(s1->len==1) {
				type=toupper((int)*(get_cstring(s1)));
				if(type!='X' && type!='Y') er=1;
			} else er=1;
		}
		if(er) {
			parerror2("Malformed clause [%s] for LINK command\n",s1?get_cstring(s1):NULL_STR);
			type=0;
		}
		if(s1) free_string(s1);
	}
	s1=0;
	if(s) s1=parse_exp_str(s);
	if(s1) strip_quotes(s1);
	/* See if we already have a link statement of the same name */
	lk=loki->prep->links;
	while(lk) {
		if(s1) {
			if(lk->name && !strcasecmp(get_cstring(s1),get_cstring(lk->name))) break;
		} else if(!lk->name) break;
		lk=lk->next;
	}
	if(!lk) {
		if(!(lk=malloc(sizeof(struct lk_link)))) ABT_FUNC(MMsg);
		if(s1) {
			lk->name=s1;
			s1=0;
		} else lk->name=0;
		lk->next=loki->prep->links;
		loki->prep->links=lk;
		lk->type=type;
	} else {
		if(lk->type!=type) parerror2("Linkage type mismatch with previous declaration\n");
	}
	if(s1) free_string(s1);
	if(!vl) {
		typeI_wait.type=WAIT_LINK;
		typeI_wait.ptr.lk=lk;
		set_typeI_kludge();
	} else {
		/* Set up marker structures */
		marker_locus_com(2,vl,0,lk);
	}
}

static void set_com(char *com,int flag,struct parse_clause *pc,struct parse_term *s,struct parse_var_list *vl)
{
	struct parse_term *v;
	struct parse_var *vv;
	
	if(pc) bad_clause(com,pc);
	if(s) bad_arg(com,s);
	while(vl) {
		v=&vl->term;
		if(v->type&VARIABLE) {
			vv=find_parse_var(v);
			vv->ltype|=flag;
		}
		vl=vl->next;
	}
}

static void set_int_com(struct parse_clause *pc,struct parse_term *s,struct parse_var_list *vl)
{
	set_com("INTEGER",ST_INTTYPE,pc,s,vl);
}

static void set_real_com(struct parse_clause *pc,struct parse_term *s,struct parse_var_list *vl)
{
	set_com("REAL",ST_REAL,pc,s,vl);
}

static void set_factor_com(struct parse_clause *pc,struct parse_term *s,struct parse_var_list *vl)
{
	set_com("FACTOR",ST_FACTOR,pc,s,vl);
}

static void set_constant_com(struct parse_clause *pc,struct parse_term *s,struct parse_var_list *vl)
{
	set_com("CONSTANT",ST_CONSTANT,pc,s,vl);
}

static void set_multiple_com(struct parse_clause *pc,struct parse_term *s,struct parse_var_list *vl)
{
	set_com("MULTIPLE",ST_MULTIPLE,pc,s,vl);
}

static void set_random_com(struct parse_clause *pc,struct parse_term *s,struct parse_var_list *vl)
{
	set_com("RANDOM",ST_RANDOM|ST_FACTOR,pc,s,vl);
}

static void ctypeI(char *com,struct parse_clause *clause,struct parse_term *s,struct parse_var_list *vl)
{
	int i;
	char *coms[]={"LOG","MISSING","PEDIGREE","MARKER","FILE","LINK","ARRAY",
		  "INTEGER","REAL","FACTOR","CONSTANT","MULTIPLE","RANDOM",0};
	string *s1;
	void (*funcs[])(struct parse_clause *,struct parse_term *,struct parse_var_list *)={
		log_com,missing_com,pedigree_com,marker_com,file_com,link_com,array_com,
		  set_int_com,set_real_com,set_factor_com,set_constant_com,set_multiple_com,set_random_com
	};
	
	if(!com) {
		switch(typeI_wait.type) {
		 case WAIT_LINK:
			/* Set up marker structures */
			marker_locus_com(2,vl,0,typeI_wait.ptr.lk);
			break;
		 case WAIT_FILE:
			if(typeI_wait.ptr.file) {
				typeI_wait.ptr.file->varlist=vl;
				print_var_list(stdout,vl);
				set_var_list_vector_flag(vl);
				fputc('\n',stdout);
			} else free_parse_var_list(vl);
			break;
		 default:
			printf("Null type I command received\n");
			break;
		}
		typeI_wait.type=0;
		return;
	}
	typeI_wait.type=0;
	i=-1;
	while(coms[++i]) if(!strcasecmp(coms[i],com)) break;
	if(coms[i]) funcs[i](clause,s,vl);
	else {
		printf("Received unknown type I command '%s'\n",com);
		if(clause) {
			printf("Clause: ");
			print_parse_clause(stdout,clause);
			fputc('\n',stdout);
		}
		if(s) {
			s1=parse_exp_str(s);
			if(s1) {
				printf("Arg: %s\n",get_cstring(s1));
				free_string(s1);
			}
		}
 		if(vl) {
			fputs("Var list:",stdout);
			print_var_list(stdout,vl);
			fputc('\n',stdout);
			free_parse_var_list(vl);
		}
	}
}

static void ctypeII(int key,char *com,struct parse_term *arg,struct parse_var_list *vl)
{
	int i;
	char *coms[]={"FREQUENCY","SEX",0};
	string *s1;
	void (*funcs[])(int key,struct parse_term *,struct parse_var_list *)={
		frequency_com,sex_com
	};
	
	i=-1;
	while(coms[++i]) if(!strcasecmp(coms[i],com)) break;
	if(coms[i]) funcs[i](key,arg,vl);
	else {
		printf("Received unknown type II command '%s'\n",com);
		if(key>=0) printf("With keyword %s\n",keys[key]);
		if(arg) {
			s1=parse_exp_str(arg);
			if(s1) {
				printf("Arg: %s\n",get_cstring(s1));
				free_string(s1);
			}
		}
		if(vl) {
			fputs("Var list:",stdout);
			print_var_list(stdout,vl);
			fputc('\n',stdout);
			free_parse_var_list(vl);
		}
	}
}

static void ctypeIII(char *com,struct parse_term *v,struct parse_term *arg)
{
	int i;
	char *coms[]={"MODEL","LET",0};
	string *s1;
	void (*funcs[])(struct parse_term *,struct parse_term *)={
		model_com,com_assign
	};
	
	i=-1;
	while(coms[++i]) if(!strcasecmp(coms[i],com)) break;
	if(coms[i]) funcs[i](v,arg);
	else {
		printf("Received unknown type III command '%s'\n",com);
		if(v) {
			s1=parse_exp_str(v);
			if(s1) {
				printf("Var: %s\n",get_cstring(s1));
				free_string(s1);
			}
			free_parse_term(v);
		}
		if(arg) {
			s1=parse_exp_str(arg);
			if(s1) {
				printf("Arg: %s\n",get_cstring(s1));
				free_string(s1);
			}
			free_parse_term(arg);
		}
	}
}

static void ctypeV(char *com,int key,struct parse_var_list *vl)
{
	int i;
	char *coms[]={"MARKER","TRAIT",0};
	
	i=-1;
	while(coms[++i]) if(!strcasecmp(coms[i],com)) break;
	if(coms[i]) marker_locus_com(key,vl,i,0);
	else {
		printf("Received unknown type V command '%s'\n",com);
		if(key>=0) printf("With keyword %s\n",keys[key]);
		if(vl) {
			fputs("Var list:",stdout);
			print_var_list(stdout,vl);
			fputc('\n',stdout);
			free_parse_var_list(vl);
		}
	}
}
static int check_markers(void)
{
	int err=0,k,k1;
	char *p;
	struct lk_marker *mk;
	struct parse_var *v,*elem;
	
	mk=loki->prep->markers;
	while(mk) {
		v=find_parse_var(mk->term);
		if((v->type&VT_ARRAY_ELEMENT) && (mk->term->elem.vv1->head->ltype&ST_MARKER))
		  v->ltype|=ST_MARKER;
		if(v->ltype&ST_LINKED) {
			if(!(v->ltype&ST_MARKER)) {
				k=k1=0;
				if(v->type&VT_ARRAY) {
					for(k=0;k<v->size;k++) {
						elem=v->arg.vv+k;
						if(elem->ltype&ST_MARKER) k1++;
					}
				}
				if(k) { 
					if(k==k1) {
						v->ltype|=(ST_MARKER|ST_REQUIRED);
						for(k=0;k<v->size;k++) v->arg.vv[k].ltype|=ST_REQUIRED;
					} else {
						p=parse_exp_cstr(mk->term);
						fprintf(stderr,"Array %s is declared in a LINK statement, but not all elements of %s were declared in a MARKER LOCUS statement\n",SAFE_PTR(p),SAFE_PTR(p));
						free(p);
						err=1;
					}
				} else {
					p=parse_exp_cstr(mk->term);
					fprintf(stderr,"Variable %s is declared in a LINK statement, but not in a MARKER LOCUS statement\n",SAFE_PTR(p));
					free(p);
					err=1;
				}
			} else {
				v->ltype|=ST_REQUIRED;
				if(v->type&VT_ARRAY) {
					for(k=0;k<v->size;k++) v->arg.vv[k].ltype|=ST_REQUIRED;
				}
			}
		}
		mk=mk->next;
	}
	return err;
}

static string *shell_com(string *s)
{
	char *p;
	int i,fg=0;
	
	if(s && s->len) {
		i=s->len;
		p=get_cstring(s);
		while(i && isspace((int)p[--i]));
		if(i>=0 && p[i]=='|') fg=1;
	}
	if(!fg) s=addn_to_string(s," |",2);
	return s;
}
static void set_parse_term_flag(struct parse_term *tm,int flag)
{
	int type;
	struct parse_var *vv;
	struct expr_op *eop;
	struct func_op *fop;
	
	type=tm->type&VT_TYPES;
	if(type==VARIABLE) {
		vv=find_parse_var(tm);
		if(!(vv->ltype&flag)) printf("Flagging '%s' with %x\n",get_cstring(vv->name),flag);
		vv->ltype|=flag;
	} else if(type==EXPR_OP) {
		eop=tm->elem.eop;
		if(eop->exp1) set_parse_term_flag(eop->exp1,flag);
		if(eop->exp2) set_parse_term_flag(eop->exp2,flag);
	} else if(type==FUNCTION) {
		fop=tm->elem.fop;
		set_parse_term_flag(fop->ex,flag);
	}
}

static int check_required(void) 
{
	struct lk_assign *models,*assign;
	struct parse_var *vv;
	int err=0;
	
	/* Flag variables in models */
	models=loki->prep->models;
	while(models) {
		vv=find_parse_var(models->var);
		vv->ltype|=ST_REQUIRED;
		printf("Flagging '%s' as required from MODEL statement (TRAIT)\n",get_cstring(vv->name));
		set_parse_term_flag(models->arg,ST_REQUIRED);
		models=models->next;
	}
	assign=loki->prep->assignments;
	while(assign) {
		vv=find_parse_var(assign->var);
		if(vv->ltype&ST_REQUIRED) {
			set_parse_term_flag(assign->arg,ST_REQUIRED);
		}
		assign=assign->next;
	}
	return err;
}

int NewReadControl(FILE *fptr,char *cname,struct loki *lk)
{
	int i,err=0;
	char *comms[]={"FILE","SEX","PEDIGREE","SET","MARKER","TRAIT","DISCRETE","MODEL","LOG",
		  "FILTER","MISSING","LINK","RANDOM","REAL","INTEGER","PRINT","MULTIPLE","CENSORED",
		  "GROUP","AFFECTED","UNAFFECTED","OUTPUT","INCLUDE","ERRORDIR","POSITION",
		  "FREQUENCY","ARRAY",0
	};
	func_def *funcs;
	
	message(INFO_MSG,"Reading in control file '%s'\n",cname);
	handle.fname=cname;
	handle.assign=assign;
	handle.ctypeI=ctypeI;
	handle.ctypeII=ctypeII;
	handle.ctypeIII=ctypeIII;
	handle.ctypeV=ctypeV;
	funcs=register_string_function(0,shell_com,"shell");
	/* Set up 'special' variables */
	for(i=0;i<N_SPECIAL_VAR;i++) {
		special_vars.string[i]=addn_to_string(0,special_name[i],strlen(special_name[i]));
		special_vars.var[i]=insert_prep_var(special_vars.string[i]);
		special_vars.type[i]=(i==3?INTEGER:STRING);
	}
	/* Register cleanup routine */
	loki=lk;
	if(atexit(FreeReadControl)) message(WARN_MSG,"Unable to register exit function FreeReadControl()\n");
	/* Call parser */
	Parse_File(fptr,keys,comms,funcs,&handle);
	/* Check structures for input errors */
	err=check_markers();
	err|=check_ped();
	err|=check_required();
	print_all_vars();
	free_functions(funcs);
	return err;
}
