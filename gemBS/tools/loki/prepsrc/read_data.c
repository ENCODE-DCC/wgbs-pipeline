/****************************************************************************
 *                                                                          *
 *     Loki - Programs for genetic analysis of complex traits using MCMC    *
 *                                                                          *
 *             Simon Heath - University of Washington                       *
 *                                                                          *
 *                       March 1997                                         *
 *                                                                          *
 * read_data.c:                                                             *
 *                                                                          *
 * Reads in data from all input files and recodes factorial data.           *
 * Also calls setup_pedigree() <setup_ped.c> to initialize pedigree         *
 * data structures.                                                         *
 *                                                                          *
 * Copyright (C) Simon C. Heath 1997, 2000, 2002                            *
 * This is free software.  You can distribute it and/or modify it           *
 * under the terms of the Modified BSD license, see the file COPYING        *
 *                                                                          *
 ****************************************************************************/

#include <config.h>
#include <stdlib.h>
#include <string.h>
#ifdef USE_DMALLOC
#include <dmalloc.h>
#endif
#include <math.h>
#include <stdio.h>
#include <ctype.h>
#include <errno.h>
#if HAVE_REGCOMP
#include <sys/types.h>
#include <regex.h>
#endif
#include <sys/wait.h>

#include "config.h"
#include "utils.h"
#include "y.tab.h"
#include "scan.h"

#define INIT_BLOCK_SIZE 128 /* Start allocating memory in blocks of INIT_BLOCK_SIZE records, doubling */
#define MAX_BLOCK_SIZE 512  /* the size if more space required up to MAX_BLOCK_SIZE */
#define S_BLOCK_SIZE 512   /* For string allocation: should be at least as big as BUFFER_SIZE */
#define BUFFER_SIZE 511   /* Maximum size of columns for free format reads */
#define LINE_COUNT 5000     /* How often to print 'At line' */

static char *StringData=0;
static char *default_rsformat="\n";
static size_t StringPos;
static struct miss_var_tag *miss_var;
static struct DataBlock *DataBlock;
static struct bin_node *node_strings,*node_ints,*node_ped_int,*node_fam_int,*node_ped_str,*node_fam_str;
static int num_nodes;
struct label_data **ped_recode,**family_recode,***factor_recode;
int *ped_recode1;
int ped_size;
  
struct miss_var_tag
{
	struct express **Missing;
	int nmiss;
};

static int get_scope(char *p)
{
	int i,j,scope;
	
	scope=j=0;
	while(*p) {
		i=toupper((int)*p);
		switch(i) {
		 case '!':
			j=1;
			break;
		 case 'P':
			scope|=j?~ST_PED:ST_PED;
			j=0;
			break;
		 case 'F':
			scope|=j?~ST_FACTOR:ST_FACTOR;
			j=0;
			break;
		 case 'G':
			scope|=j?~(ST_MARKER|ST_HAPLO):(ST_MARKER|ST_HAPLO);
			j=0;
			break;
		 case 'I':
			scope|=j?~ST_INTTYPE:ST_INTTYPE;
			break;
		 case 'C':
		 case 'R':
			scope|=j?~ST_REALTYPE:ST_REALTYPE;
			j=0;
			break;
		 default:
			ABT_FUNC("Illegal missing scope\n");
		}
		p++;
	}
	return scope;
}

static void process_missing(struct InFile *infile) 
{
	struct Miss *ms;
	int i,j,col=0,total_miss=0,scope=0;
	struct express **p;
	struct var_element *elem;
	
	if(!(miss_var=malloc(sizeof(struct miss_var_tag)*infile->ncol))) ABT_FUNC(MMsg);
	for(i=0;i<infile->ncol;i++) miss_var[i].nmiss=0;
	ms=Miss;
	while(ms) {
		if(ms->scope) scope=get_scope(ms->scope);
		for(col=j=0;j<infile->nvar;j++) if(infile->element[j]) {
			if(ms->element) {
				for(i=0;i<ms->nvar;i++) if(ms->element[i]==infile->element[j]) {
					miss_var[col].nmiss++;
					total_miss++;
				}
			} else if(ms->scope) {
				elem=infile->element[j];
				if(elem->type&scope) {
					miss_var[col].nmiss++;
					total_miss++;
				}
			} else {
				miss_var[col].nmiss++;
				total_miss++;
			}
			col++;
		}
		ms=ms->next;
	}
	if(total_miss) {
		if(!(p=malloc(sizeof(void *)*total_miss))) ABT_FUNC(MMsg);
	} else p=0;
	for(i=0;i<col;i++) {
		if(miss_var[i].nmiss) {
			miss_var[i].Missing=p;
			p+=miss_var[i].nmiss;
			miss_var[i].nmiss=0;
		} else miss_var[i].Missing=0;
	}
	ms=Miss;
	while(ms) {
		if(ms->scope) scope=get_scope(ms->scope);
		for(col=j=0;j<infile->nvar;j++) if(infile->element[j]) {
			if(ms->element) {
				for(i=0;i<ms->nvar;i++) if(ms->element[i]==infile->element[j])
				  miss_var[col].Missing[miss_var[col].nmiss++]=&ms->Missing;
			} else if(ms->scope) {
				elem=infile->element[j];
				if(elem->type&scope) 
					miss_var[col].Missing[miss_var[col].nmiss++]=&ms->Missing;
			} else miss_var[col].Missing[miss_var[col].nmiss++]=&ms->Missing;
			col++;
		}
		ms=ms->next;
	}
}

static struct bin_node *alloc_node(const void *s,int type)
{
	struct bin_node *node;
	struct label_data *data;
	size_t i;
	char *p,*p1;
	
	if(!(node=malloc(sizeof(struct bin_node)))) ABT_FUNC(MMsg);
	node->left=node->right=0;
	node->balance=0;
	if(!(data=malloc(sizeof(struct label_data)))) ABT_FUNC(MMsg);
	node->data=data;
	data->type=type;
	data->index=num_nodes++;
	switch(type) {
	 case STRING:
		p1=(char *)s;
		i=strlen(p1)+1;
		if(i>(S_BLOCK_SIZE-StringPos)) {
			if(i>S_BLOCK_SIZE) ABT_FUNC("String too long - increase S_BLOCK_SIZE\n");
			if(!(StringData=malloc(S_BLOCK_SIZE))) ABT_FUNC(MMsg);
			RemBlock=AddRemem(StringData,RemBlock);
			StringPos=0;
		}
		p=StringData+StringPos;
		data->data.string=p;
		if(!syst_var[IGNORE_CASE]) (void)strncpy(p,p1,i);
		else {
			p1=(char *)s;
			while(*p1) *p++=toupper((int)(*p1++));
			*p=0;
		}
		StringPos+=i;
		break;
	 case INTEGER:
		data->data.value= *(long *)s;
		break;
	 default:
		ABT_FUNC("Internal error - invalid type\n");
	}
	return node;
}

static struct label_data *find_inode(struct bin_node *node,const long s)
{
	long i;
	struct label_data *data;
	
	data=node->data;
	i=s-data->data.value;
	if(i<0) {
		if(node->left) data=find_inode(node->left,s);
		else data=0;
	} else if(i>0) {
		if(node->right) data=find_inode(node->right,s);
		else data=0;
	}
	return data;
}

static struct label_data *find_snode(struct bin_node *node,const char *s)
{
	long i;
	struct label_data *data;
	
	data=node->data;
	if(!syst_var[IGNORE_CASE]) i=strcmp(s,data->data.string);
	else i=strcasecmp(s,data->data.string);
	if(i<0) {
		if(node->left) data=find_snode(node->left,s);
		else data=0;
	} else if(i>0) {
		if(node->right) data=find_snode(node->right,s);
		else data=0;
	}
	return data;
}

static struct bin_node *insert_node(struct bin_node *node,const void *s,struct label_data **node1,int *balanced,int type)
{
	long i;
	struct label_data *data;
	
	data=node->data;
	if(type==INTEGER) i=*((long *)s)-data->data.value;
	else {
		if(!syst_var[IGNORE_CASE]) i=strcmp(s,data->data.string);
		else i=strcasecmp((char *)s,data->data.string);
	}
	if(i<0) {
		if(node->left) node->left=insert_node(node->left,s,node1,balanced,type);
		else {
			node->left=alloc_node(s,type);
			*node1=node->left->data;
			*balanced=0;
		}
		if(!(*balanced)) {
			switch(node->balance) {
			 case -1:
				node=rotate_left(node);
				*balanced=1;
				break;
			 case 0:
				node->balance=-1;
				break;
			 case 1:
				node->balance=0;
				*balanced=1;
			}
		}
	} else if(i>0) {
		if(node->right) node->right=insert_node(node->right,s,node1,balanced,type);
		else {
			node->right=alloc_node(s,type);
			*node1=node->right->data;
			*balanced=0;
		}
		if(!(*balanced)) {
			switch(node->balance) {
			 case -1:
				node->balance=0;
				*balanced=1;
				break;
			 case 0:
				node->balance=1;
				break;
			 case 1:
				node=rotate_right(node);
				*balanced=1;
			}
		}
	} else {
		*node1=node->data;
		*balanced=1;
	}
	return node;
}

struct label_data *find_node(const void *s,int type,int flag)
{
	struct bin_node *node;
	struct label_data *data=0;

	if(type==INTEGER) {
		if(flag) node=node_fam_int;
		else node=node_ped_int;
		if(node) data=find_inode(node,*((long *)s));
	} else {
		if(flag) node=node_fam_str;
		else node=node_ped_str;
		if(node) data=find_snode(node,(char *)s);
	}
	return data;
}

void free_nodes(void)
{
	if(node_ints) free_bin_tree(node_ints,free);
	if(node_strings) free_bin_tree(node_strings,free);
	if(node_ped_int) free_bin_tree(node_ped_int,free);
	if(node_ped_str) free_bin_tree(node_ped_str,free);
	if(node_fam_int) free_bin_tree(node_fam_int,free);
	if(node_fam_str) free_bin_tree(node_fam_str,free);
	node_ints=node_strings=node_ped_int=node_ped_str=node_fam_int=node_fam_str=0;
}

static void set_missing(int col,int ncol,int rec,struct DataBlock *db)
{
	unsigned int i,j;
	
	
	i=rec*ncol+col;
	db->records[i].value=0;
	j=7-(i&7);
	i>>=3;
	db->type[i]|=(unsigned char)(1<<j);
}

int check_missing(int col,int ncol,int rec,struct DataBlock *db)
{
	int i,j;
	i=rec*ncol+col;
	j=7-(i&7);
	i>>=3;
	return db->type[i]&(1<<j);
}

static int check_miss(char *string,struct express *miss,int elem_type)
{
	int match=0;
	double rval;
	long val;
	char *p;
	struct express ep;
	
	ep.type=0;
	if((elem_type&(ST_INTTYPE|ST_REALTYPE)) && (miss->type==ST_STRING)) {
		p=string;
		while(*p) {
			if(*p=='.') break;
			p++;
		}
		if(*p) {
			ep.arg.rvalue=strtod(miss->arg.string,&p);
			if(!*p) ep.type=ST_REAL;
		} else {
			ep.arg.value=strtol(miss->arg.string,&p,10);
			if(!*p) ep.type=ST_INTEGER;
		}
	} 
	if(!ep.type) {
		ep.type=miss->type;
		ep.arg=miss->arg;
	}
	switch(ep.type) {
	 case ST_STRING:
		match=!strcmp(string,ep.arg.string);
		break;
	 case ST_INTEGER:
		p=string;
		while(*p) {
			if(*p=='.') break;
			p++;
		}
		if(*p) {
			rval=strtod(string,&p);
			if(!*p && fabs(rval-(double)ep.arg.value)<1.0e-12) match=1;
		} else {
			val=strtol(string,&p,10);
			if(!*p && val==ep.arg.value) match=1;
		}
		break;
	 case ST_REAL:
		rval=strtod(string,&p);
		if(!*p && fabs(rval-ep.arg.rvalue)<1.0e-12) match=1;
		break;
	 case 0:
		break;
#ifdef DEBUG
	 default:
		ABT_FUNC("Internal error - bad missing type\n");
#endif
	}
	return match;
}

static int handle_string(char *string,struct InFile *infile,struct var_element *elem,int col,char *gs)
{
	int i,miss,ncol;
	char *sptr,ch=0,*p1,*p2,*gs1,flag=0;
	static char single1[2],single2[2];
	struct label_data *node;
	struct gt_data *gt;
	long value;

	ncol=infile->ncol;
	miss=0;
	if(miss_var) {
		for(i=0;i<miss_var[col].nmiss;i++) {
			if(check_miss(string,miss_var[col].Missing[i],elem->type)) {
				miss=1;
				break;  
			}
		}
	}
	if(!miss) {
		if(elem->type&(ST_ID|ST_SIRE|ST_DAM|ST_FAMILY|ST_FACTOR|ST_HAPLO|ST_MARKER)) {
			if(elem->type&ST_MARKER) {
				/* Genotype data - split into fields based on gsformat */
				p1=string;
				p2=0;
				if(gs) {
					if(!gs[0]) {
						single1[0]=*p1++;
						single2[0]=*p1++;
						if(*p1) print_scan_err("[%s:%d] %s(): Line %d column %d - >2 characters in simple genotype %s\n",__FILE__,__LINE__,__func__,lineno,col+1,string);
						p1=single1;
						p2=single2;
						flag=1;
					} else {
						while((ch=*(p1++))) {
							gs1=gs;
							while(*gs1) {
								if(ch==*gs1) {
									*(--p1)=0;
									break;
								}
								gs1++;
							}
							if(*gs1) break;
						}
					}
				} else {
					while((ch=*(p1++))) {
						if(isspace((int)ch)) {
							*(--p1)=0;
							break;
						}
					}
				}
				if(!flag) {
					if(ch) {
						p2=p1+1;
						qstrip(p2);
					}
					p1=string;
					qstrip(p1);
				}
				if(miss_var) {
					for(i=0;i<miss_var[col].nmiss;i++) {
						if(p1) if(check_miss(p1,miss_var[col].Missing[i],elem->type)) {
							p1=0;
						} 
						if(p2) if(check_miss(p2,miss_var[col].Missing[i],elem->type)) {
							p2=0;
						} 
						if(!(p1||p2)) break;
					}
				}
				if(p1||p2) {
					if(!(gt=malloc(sizeof(struct gt_data)))) ABT_FUNC(MMsg);
					gt->node1=gt->node2=0;	
					if(elem->type&(ST_INTTYPE)) {
						if(p1) {
							value=strtol(p1,&sptr,10);
							if(*sptr) {
								if(!syst_var[SKIP_BAD_INTS]) 
								  print_scan_err("[%s:%d] %s(): Line %d column %d - Malformed integer %s\n",__FILE__,__LINE__,__func__,lineno,col+1,p1);
								miss=1;
								free(gt);
								gt=0;
							} else {
								if(!node_ints) {
									node_ints=alloc_node(&value,INTEGER);
									node=node_ints->data;
								} else node_ints=insert_node(node_ints,&value,&node,&i,INTEGER);
								gt->node1=node;
							}
						}
						if(!miss && p2) {
							value=strtol(p2,&sptr,10);
							if(*sptr) {
								if(!syst_var[SKIP_BAD_INTS]) 
								  print_scan_err("[%s:%d] %s(): Line %d column %d - Malformed integer %s\n",__FILE__,__LINE__,__func__,lineno,col+1,p2);
								miss=1;
								free(gt);
								gt=0;
							} else {
								if(!node_ints) {
									node_ints=alloc_node(&value,INTEGER);
									node=node_ints->data;
								} else node_ints=insert_node(node_ints,&value,&node,&i,INTEGER);
								gt->node2=node;
							}
						}
						DataBlock->records[DataBlock->record_ptr*ncol+col].gt_data=gt;
					} else {
						if(p1) {
							if(!node_strings) {
								node_strings=alloc_node(p1,STRING);
								node=node_strings->data;
							} else node_strings=insert_node(node_strings,p1,&node,&i,STRING);
							gt->node1=node;
						}
						if(p2) {
							if(!node_strings) {
								node_strings=alloc_node(p2,STRING);
								node=node_strings->data;
							} else node_strings=insert_node(node_strings,p2,&node,&i,STRING);
							gt->node2=node;
						}
						DataBlock->records[DataBlock->record_ptr*ncol+col].gt_data=gt;
					}
				} else miss=1;
			} else {
				if(elem->type&ST_INTTYPE) {
					value=strtol(string,&sptr,10);
					if(*sptr) {
						if(!syst_var[SKIP_BAD_INTS]) 
						  print_scan_err("[%s:%d] %s(): Line %d column %d - Malformed integer %s\n",__FILE__,__LINE__,__func__,lineno,col+1,string);
						miss=1;
					} else {
						if(elem->type&(ST_ID|ST_SIRE|ST_DAM))	{
							if(!node_ped_int) {
								node_ped_int=alloc_node(&value,INTEGER);
								node=node_ped_int->data;
							} else node_ped_int=insert_node(node_ped_int,&value,&node,&i,INTEGER);
						} else if(elem->type&(ST_FAMILY)) {
							if(!node_fam_int) {
								node_fam_int=alloc_node(&value,INTEGER);
								node=node_fam_int->data;
							} else node_fam_int=insert_node(node_fam_int,&value,&node,&i,INTEGER);
						} else {
							if(!node_ints) {
								node_ints=alloc_node(&value,INTEGER);
								node=node_ints->data;
							} else node_ints=insert_node(node_ints,&value,&node,&i,INTEGER);
						}
						DataBlock->records[DataBlock->record_ptr*ncol+col].node=node;
					}
				} else {
					if(elem->type&(ST_ID|ST_SIRE|ST_DAM))	{
						if(!node_ped_str) {
							node_ped_str=alloc_node(string,STRING);
							node=node_ped_str->data;
						} else node_ped_str=insert_node(node_ped_str,string,&node,&i,STRING);
					} else if(elem->type&(ST_FAMILY)) {
						if(!node_fam_str) {
							node_fam_str=alloc_node(string,STRING);
							node=node_fam_str->data;
						} else node_fam_str=insert_node(node_fam_str,string,&node,&i,STRING);
					} else {
						if(!node_strings) {
							node_strings=alloc_node(string,STRING);
							node=node_strings->data;
						} else node_strings=insert_node(node_strings,string,&node,&i,STRING);
					}
					DataBlock->records[DataBlock->record_ptr*ncol+col].node=node;
				}
			}
		} else {
			if(elem->type&ST_INTTYPE) DataBlock->records[DataBlock->record_ptr*ncol+col].value=strtol(string,&sptr,10);
			else DataBlock->records[DataBlock->record_ptr*ncol+col].rvalue=strtod(string,&sptr);
			if(*sptr) {
				if(elem->type&ST_INTTYPE) {
					if(!syst_var[SKIP_BAD_INTS])
					  print_scan_err("[%s:%d] %s(): Line %d column %d - Malformed integer %s\n",__FILE__,__LINE__,__func__,lineno,col+1,string);
				} else {
					if(!syst_var[SKIP_BAD_REALS])
					  print_scan_err("[%s:%d] %s(): Line %d column %d - Malformed floating point number %s\n",__FILE__,__LINE__,__func__,lineno,col+1,string);
				}
				miss=1;
			}
		}
	}
	return miss;
}

static void update_datablock(int ncol)
{
	int i;
	
	DataBlock->record_ptr++;
	if(DataBlock->record_ptr==DataBlock->blocksize) {
		if(!(DataBlock->next=malloc(sizeof(struct DataBlock)))) ABT_FUNC(MMsg);
		DataBlock->next->blocksize=DataBlock->blocksize;
		DataBlock=DataBlock->next;
		if(DataBlock->blocksize<MAX_BLOCK_SIZE) DataBlock->blocksize*=2;
		i=DataBlock->blocksize*ncol;
		if(!(DataBlock->records=malloc(sizeof(union DataRec)*i))) ABT_FUNC(MMsg);
		if(!(DataBlock->type=calloc((size_t)(i>>3),sizeof(char)))) ABT_FUNC(MMsg);
		DataBlock->next=0;
		DataBlock->record_ptr=0;
	}
}

void ReadData(char *lfile)
{
	int i,j,ncol,num_records,col,realcol,eor,miss;
	struct InFile *infile;
	char *buffer,ch,*string,*string_start,oldchar,*rs,*fs,*gs,*gs1,*tname;
	int buffersize,fs_reg=0,skip,skip_this;
	FILE *fptr,*flog;
	size_t l,string_leng;
	struct recode_table_tag *recode_table=0;
	struct format *format;
	struct var_element *elem;
	struct sex_def *sd;
#ifdef HAVE_REGCOMP	
	regex_t preg;
	regmatch_t pmatch;
#endif

	infile=Infiles;
	if(!(StringData=malloc(S_BLOCK_SIZE))) ABT_FUNC(MMsg);
	RemBlock=AddRemem(StringData,RemBlock);
	StringPos=0;
	if(lfile && (tname=add_file_dir(lfile))) {
		flog=fopen(tname,"a");
		free(tname);
	} else flog=0;
	if(flog) i=fputs("\n******************** Reading in data ********************\n\n",flog);
	while(infile) {
		fs_reg=0;
		fs=0;
		rs=rsformat;
		gs=gsformat;
		skip=file_skip;
		if(infile->fformat) {
			fs=infile->fformat->fs;
			rs=infile->fformat->rs;
			gs=infile->fformat->gs;
			skip=infile->fformat->skip;
			if(fs) {
				if(strlen(fs)>1) {
#if HAVE_REGCOMP
					if((i=regcomp(&preg,fs,REG_EXTENDED))) fs=0;
					else fs_reg=1;
#endif				  
				} else if(fs[0]==' ') fs=0;
			}
		}
		if(!rs) rs=default_rsformat;
		if(rs[0]=='\0') rs=0;
		if(!scan_error_n) {
			ncol=infile->ncol;
			if(Miss) process_missing(infile);
			else miss_var=0;
			if(infile->shell_flag) {
				if(!(fptr=popen(infile->name,"r"))) {
					(void)fprintf(stderr,"Can't execute '%s'.\n",infile->name);
					perror("read_data()");
					scan_error_n++;
					break;
				}
				(void)printf("Reading in output from shell command '%s'\n",infile->name);
			} else {
				if(!(fptr=fopen(infile->name,"r"))) {
					(void)fprintf(stderr,"Can't open data file '%s' for reading.\n",infile->name);
					perror("read_data()");
					scan_error_n++;
					break;
				}
				(void)printf("Reading in data from file '%s'\n",infile->name);
			}
			num_records=0;
			lineno=1;
			if(!(DataBlock=malloc(sizeof(struct DataBlock)))) ABT_FUNC(MMsg);
			i=INIT_BLOCK_SIZE*ncol;
			if(!(DataBlock->records=malloc(sizeof(union DataRec)*i))) ABT_FUNC(MMsg);
			if(!(DataBlock->type=calloc((size_t)i>>3,sizeof(unsigned int)))) ABT_FUNC(MMsg);
			DataBlock->next=0;
			DataBlock->blocksize=INIT_BLOCK_SIZE;
			DataBlock->record_ptr=0;
			infile->ncol=ncol;
			infile->data=DataBlock;
			format=infile->format;
			if(format) buffersize=2+format->f_atoms[format->n_atoms-1].pos+format->f_atoms[format->n_atoms-1].size;
			else buffersize=BUFFER_SIZE+1;
			if(!(buffer=malloc((size_t)buffersize))) ABT_FUNC(MMsg);
			col=realcol=0;
			string=buffer;
			errno=0;
			for(;;) {
				if(format) {
					if(!fgets(buffer,buffersize,fptr)) {
						if(errno && errno!=ECHILD) {
							perror("ReadData()");
							abt(__FILE__,__LINE__,"Aborting\n");
						}
						errno=0;
						break;
					}
					l=strlen(buffer);
					if(lineno>skip) {
						for(i=0;i<(int)l;i++) if(!isspace((int)buffer[i])) break;
						if(i<(int)l) {
							for(col=realcol=0;realcol<infile->nvar;realcol++) {
								elem=infile->element[realcol];
								if(elem) {
									i=format->f_atoms[realcol].pos;
									if((size_t)i>=l) set_missing(col,ncol,DataBlock->record_ptr,DataBlock);
									else {
										j=format->f_atoms[realcol].size;
										string_start=buffer+i;
										oldchar=buffer[i+j];
										buffer[i+j]=0;
										if(string_start[0]) qstrip(string_start);
										miss=string_start[0]?handle_string(string_start,infile,elem,col,gs):1;
										if(miss) {
											set_missing(col,ncol,DataBlock->record_ptr,DataBlock);
											if(elem->type&ST_ID) print_scan_err("[%s:%d] %s(): Line %d column %d - Missing id variable\n",__FILE__,__LINE__,__func__,lineno,col+1);
											if(elem->type&ST_FAMILY) print_scan_err("[%s:%d] %s(): Line %d column %d - Missing family variable\n",__FILE__,__LINE__,__func__,lineno,col+1);
										}
										buffer[i+j]=oldchar;
									}
									col++;
								}
							}
							update_datablock(ncol);
							num_records++;
						}
					}
					lineno++;
					if(!(lineno%LINE_COUNT)) (void)printf("At %d\n",lineno);
					/* read to end of line */
					while(buffer[l-1]!='\n') {
						if(!fgets(buffer,buffersize,fptr)) {
							if(errno) {
								perror("ReadData()");
								abt(__FILE__,__LINE__,"Aborting\n");
							}
							break;
						}
						l=strlen(buffer);
					}
					continue;
				}
				l=BUFFER_SIZE;
				if(l<=(size_t)(string-buffer)) ABT_FUNC("column width exceeds buffer size\n");
				l-=(size_t)(string-buffer);
				errno=0;
				if(!(l=fread(string,1,l,fptr))) {
					if(errno) {
						perror("ReadData()");
						abt(__FILE__,__LINE__,"Aborting\n");
					}
					break;
				}
				*(string+l)=(char)0;
				string=buffer;
				for(;;) {
					eor=0;
					skip_this=(lineno<=skip)?1:0;
					if(fs) {
						string_start=string;
#if HAVE_REGCOMP						
						if(fs_reg) {
							if(regexec(&preg,string_start,1,&pmatch,0)) {
								while((ch= *(string++))) {
									if(ch=='\n') lineno++;
									if(rs) if(ch==rs[0])	{
										eor=1;
										break;
									}
								}
							} else {
								string=string_start+pmatch.rm_eo;
								for(i=0;i<(int)pmatch.rm_eo;i++)	{
									ch=string_start[i];
									if(ch=='\n') lineno++;
									if(rs) if(ch==rs[0])	{
										eor=1;
										string=string_start+i+1;
										break;
									}
								}
								if(!eor) string_start[pmatch.rm_so]='\0';
								ch=1;
							}
						} else
#endif						  
						{
							while((ch= *(string++))) {
								if(ch=='\n') lineno++;
								if(ch==fs[0]) break;
								if(rs) if(ch==rs[0])	{
									eor=1;
									break;
								}	
							}
						}
					} else {
						while((ch= *(string++))) {
							if(ch=='\n') lineno++;
							if(rs) if(ch==rs[0]) {
								eor=1;
								break;
							}
							if(!isspace((int)ch)) break;
						}
						if(!ch) {
							string=buffer;
							break;
						}
						string_start=string-1;
						elem=infile->element[realcol];
						if(elem && elem->type&ST_MARKER) {
							if(gs && *gs && (!isspace((int)*gs) || gs[1])) i=0;
							else i=1;
						} else i=0;
						if(!eor)	{
							while((ch= *(string++))) {
								if(rs) if(ch==rs[0]) {
									eor=1;
									break;
								}
								j=0;
								if(i && gs) {
									gs1=gs;
									while(*gs1) if(ch==*(gs1++)) {
										j=1; 
										break;
									}
								}
								if(j || isspace((int)ch)) {
									if(!i) break;
									i--;
									while((ch=*(string++))) {
										if(ch=='\n') lineno++;
										if(rs) if(ch==rs[0]) {
											eor=1;
											break;
										}
										if(!isspace((int)ch)) break;
									}
									if(!ch) break;
								}
							}
							if(ch=='\n') lineno++;
						}
					}
					*(string-1)=(char)0;
					string_leng=strlen(string_start);
					if(!ch) {
						if(buffer!=string_start) (void)memmove(buffer,string_start,string_leng);
						string=buffer+string_leng;
						break;
					}
					if(string_start[0]) qstrip(string_start);
					if(!skip_this) {
						elem=infile->element[realcol];
						if(elem && (string_start[0] || !eor)) {
							miss=string_start[0]?handle_string(string_start,infile,elem,col,gs):1;
							if(miss) {
								set_missing(col,ncol,DataBlock->record_ptr,DataBlock);
								if(elem->type&ST_ID) print_scan_err("[%s:%d] %s(): Line %d column %d - Missing id variable\n",__FILE__,__LINE__,__func__,lineno,col+1);
								if(elem->type&ST_FAMILY) print_scan_err("[%s:%d] %s(): Line %d column %d - Missing family variable\n",__FILE__,__LINE__,__func__,lineno,col+1);
							}
							col++;
						}
					}
					realcol++;
					if(realcol==infile->nvar || eor || (rs && col==ncol))	{
						if(col) {
							for(;col<ncol;col++) set_missing(col,ncol,DataBlock->record_ptr,DataBlock);
							num_records++;
							update_datablock(ncol);
						}
						realcol=col=0;
						if(rs) while(!eor) {
							while((ch= *(string++))) {
								if(ch=='\n') lineno++;
								if(ch==rs[0]) {
									eor=1;
									break;
								}
							}
							if(!ch) {
								string=buffer;
								if(!(l=fread(string,1,BUFFER_SIZE,fptr)))	{
									if(errno) {
										perror("ReadData()");
										abt(__FILE__,__LINE__,"Aborting\n");
									}
									eor=1;
									break;
								}
								*(string+l)=(char)0;
							}
						}
						if(!(lineno%LINE_COUNT)) (void)printf("At %d\n",lineno);
					}
				}
			}
			if(miss_var) {
				for(i=0;i<ncol;i++) if(miss_var[i].nmiss) {
					free(miss_var[i].Missing);
					break;
				}
				free(miss_var);
			}
			(void)fclose(fptr);
			if(infile->shell_flag) do i=wait(&j); while(i>0);
#if HAVE_REGCOMP			
			if(fs_reg) regfree(&preg);
#endif			
			infile->n_records=num_records;
			(void)printf("%d records of %d columns read in\n",num_records,ncol);
			if(flog) (void)fprintf(flog,"     Read %d records of %d columns from '%s'\n",num_records,ncol,infile->name);
			free(buffer);
		}
		infile=infile->next;
	}
	if(scan_error_n || !num_nodes) return;
	if(flog && flog!=stdout) (void)fclose(flog);
	if(!(recode_table=calloc((size_t)num_nodes,sizeof(struct recode_table_tag)))) ABT_FUNC(MMsg);
	for(i=j=0;i<n_markers;i++) {
		if(markers[i].element && !(markers[i].element->type&ST_NOT_REALLY_REQUIRED)) {
			memcpy(markers+j,markers+i,sizeof(struct Marker));
			j++;
		}
	}
	n_markers=j;
	if(!j) {
		free(markers);
		markers=0;
	}
	setup_pedigree(num_nodes,recode_table,lfile);
	if(scan_error_n) return;
	for(i=0;i<num_nodes;i++) recode_table[i].index=0;
	recode_factors(num_nodes,recode_table);
	if(scan_error_n) return;
	if(n_markers) recode_marker_data(num_nodes,recode_table);
	if(recode_table) free(recode_table);
	if(scan_error_n) return;
	match_records();
	if(Restrictions) restrict_data();
	if(Censored) censored_data();
	if(Affected) affected_data();
	if(Proband) proband_data();
	check_sex();
	handle_groups(lfile);
	cleanup_unused();
	sd=sex_def;
	while(sd) {
		for(i=0;i<2;i++) if(sd->sex_exp[i]) {
			if(sd->sex_exp[i]->type==ST_STRING) if(sd->sex_exp[i]->arg.string) free(sd->sex_exp[i]->arg.string);
			free(sd->sex_exp[i]);
		}
		sex_def=sd->next;
		free(sd);
		sd=sex_def;
	}
	prune_pedigree(lfile);
	check_sex2();
	count_components(lfile);
/*	count_relationships(); */
	if(OutputFile) Output_Data();
	if(OutputRawFile) Output_Raw_Data();
}

