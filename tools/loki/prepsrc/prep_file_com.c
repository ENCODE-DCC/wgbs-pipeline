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
#include "lk_malloc.h"

static int count_format_elem(struct format_clause *fc)
{
	int i=0;
	
	while(fc) {
		switch(fc->type) {
		 case FC_READ:
			i++;
			break;
		 case FC_SUBCLAUSE:
			i+=fc->multiplier*count_format_elem(fc->elem);
			break;
		}
		fc=fc->next;
	}
	return i;
}

static int collect_format_elem(struct format_clause *fc,int *start,int *end,int *idx,int col)
{
	int i;
	  
	while(fc) {
		switch(fc->type) {
		 case FC_READ:
			start[*idx]=col;
			end[(*idx)++]=col+fc->multiplier;
			/* This fall through is intentional! */
		 case FC_SKIP:
			col+=fc->multiplier;
			break;
		 case FC_SUBCLAUSE:
			for(i=0;i<fc->multiplier;i++) col=collect_format_elem(fc->elem,start,end,idx,col);
			break;
		}
		fc=fc->next;
	}
	return col;
}
										  
int count_input_cols(struct parse_var_list *vl) 
{
	int c=0;
	struct parse_term *v;
	struct parse_var *vv;
	
	while(vl) {
		v=&vl->term;
		if((v->type&VT_TYPES)==VARIABLE) {
			if(c<0) {
				parerror2("An array without a fixed size must come at the end of the parameter list for a FILE statement\n");
				c=0;
				break;
			}
			vv=v->elem.vv;
			if(vv->type&VT_ARRAY) {
				if(vv->type&VT_SIZE_SET) c+=vv->size;
				else c=-(c+1);
			} else c++;
		}
		vl=vl->next;
	}
	return c;
}

struct prep_file *prep_file_com(struct parse_clause *pc,struct parse_term *s,struct parse_var_list *vl,struct special_vars *special,struct loki *loki)
{
	int i,j,er,type,flag[N_SPECIAL_VAR];
	char *p;
	struct prep_file *file=0;
	struct format_clause *fc,*fc1;
	struct parse_term *tm,*tm1;
	struct parse_var *vv;
	string *s1;
	
	p=parse_exp_cstr(s);
	if(!p) {
		parerror2("Can't resolve filename for FILE command\n");
		if(vl) free_parse_var_list(vl);
	} else {
		for(i=0;i<N_SPECIAL_VAR;i++) flag[i]=0;
		fputs("FILE ",stdout);
		if(pc) {
			fputc('[',stdout);
			print_parse_clause(stdout,pc);
			fputs("] ",stdout);
		}
		file=lk_malloc(sizeof(struct prep_file));
		file->next=loki->prep->data->infiles;
		loki->prep->data->infiles=file;
		i=(int)strlen(p);
		assert(i);
		if(p[0]=='\"' && p[i-1]=='\"') {
			memcpy(p,p+1,i-2);
			p[i-2]=0;
		}
		file->name=p;
		file->varlist=vl;
		file->data=0;
		file->fs=file->gs=file->rs=0;
		file->skip=0;
		fc=fc1=0;
		/* Handle Clause */
		while(pc) {
			if(pc->type==PC_FORMAT_CLAUSE) {
				if(fc1) {
					fc1->next=pc->clause.format;
					fc1=fc1->next;
				} else fc=fc1=pc->clause.format;
			} else if((tm=pc->clause.term)) {
				if(tm->clause_assign) {
					tm1=tm->clause_assign->var;
					if(tm1->type&VT_ARRAY_ELEMENT) vv=0;
					else vv=tm1->elem.vv;
					if(vv) {
						for(i=0;i<N_SPECIAL_VAR;i++) if(vv==special->var[i]) break;
						if(i<N_SPECIAL_VAR) {
							tm1=tm->clause_assign->expr;
							type=tm1->type&VT_TYPES;
							if(type==STRING || type==INTEGER || type==REAL) {
								if(special->type[i]==STRING) {
									if(type!=STRING) tm1->type=convert_string(tm1);
									flag[i]=1;
									switch(i) {
									 case 0:
										if(file->fs) free_string(file->fs);
										file->fs=copy_string(tm1->elem.str);
										break;
									 case 1:
										if(file->rs) free_string(file->fs);
										file->rs=copy_string(tm1->elem.str);
										break;
									 case 2:
										if(file->gs) free_string(file->fs);
										file->gs=copy_string(tm1->elem.str);
										break;
									}
								} else {
									j=convert_to_int(tm1,&er);
									if(!er && j>=0) {
										flag[i]=1;
										file->skip=j;
									}
								}
							}
						}
					}
				}
			}
			pc=pc->next;
		}
		/* Take format values from global variables, if not already defined */
		for(i=0;i<N_SPECIAL_VAR;i++) if(!flag[i]) {
			vv=special->var[i];
			tm1=get_var_term(vv);
			tm1=try_resolve_expr(tm1);
			type=tm1->type&VT_TYPES;
			if(type==STRING || type==INTEGER || type==REAL) {
				if(special->type[i]==STRING) {
					flag[i]=1;
					if(type!=STRING) tm1->type=convert_string(tm1);
					switch(i) {
					 case 0:

						 break;
					 case 1:
						file->rs=copy_string(tm1->elem.str);
						break;
					 case 2:
						file->gs=copy_string(tm1->elem.str);
						break;
					}
				} else {
					j=convert_to_int(tm1,&er);
					if(!er && j>=0) {
						flag[i]=1;
						file->skip=j;
					}
				}
			}
			free_parse_term(tm1);
		}
		printf("[FS='%s', RS='%s', GS='%s', skip=%d] \"%s\"",file->fs?get_cstring(file->fs):NULL_STR,file->rs?get_cstring(file->rs):NULL_STR,file->gs?get_cstring(file->gs):NULL_STR,file->skip,p);
 		fputc(' ',stdout);
		if(vl) {
			print_var_list(stdout,vl);
			fputc('\n',stdout);
			set_var_list_vector_flag(vl);
		}
		i=0;
		file->ncol=0;
		file->fixed_start=file->fixed_end=0;
		if(fc) {
			i=count_format_elem(fc);
			if(i) {
				file->fixed_start=lk_malloc(sizeof(int)*2*i);
				file->fixed_end=file->fixed_start+i;
				i=0;
				(void)collect_format_elem(fc,file->fixed_start,file->fixed_end,&i,0);
				file->ncol=i;
			} else {
				s1=format_clause_str(fc);
				p=s1?get_cstring(s1):0;
				parerror2("Fixed format [%s] for FILE command has no data columns!\n",SAFE_PTR(p));
				if(s1) free_string(s1);
			}
		}
		/* Count variables in list */
		if(vl) i=count_input_cols(vl);
		if(file->ncol) {
			if(i!=file->ncol) {
				if(i>0) {
					parerror2("Wrong number of variables for fixed format read\n");
					file->ncol=0;
				} else if(i<0 && -i>file->ncol) {
					parerror2("Too many variables for fixed format read\n");
					file->ncol=0;
				}
			}
		} else file->ncol=i;
	}
	return file;
}

