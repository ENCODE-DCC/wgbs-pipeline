%{
/****************************************************************************
 *                                                                          *
 *     Loki - Programs for genetic analysis of complex traits using MCMC    *
 *                                                                          *
 *             Simon Heath - University of Washington                       *
 *                                                                          *
 *                       July 1997                                          *
 *                                                                          *
 * param_lex.l:                                                             *
 *                                                                          *
 * flex (NOT lex!) source for parameter file lexer.                         *
 *                                                                          *
 ****************************************************************************/

#include <config.h>
#include <stdlib.h>
#include <string.h>
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif
#ifdef USE_DMALLOC
#include <dmalloc.h>
#endif
#include <sys/stat.h>
#include "y.tab.h"
#include "loki_scan.h"
#include "utils.h"

#ifndef YY_STACK_USED
#define YY_STACK_USED 0
#endif
#ifndef YY_ALWAYS_INTERACTIVE
#define YY_ALWAYS_INTERACTIVE 0
#endif
#ifndef YY_NEVER_INTERACTIVE
#define YY_NEVER_INTERACTIVE 0
#endif
#ifndef YY_MAIN
#define YY_MAIN 0
#endif

#define LINEBUFSIZE 512

extern void yyerror(char *s),print_scan_err(char *fmt, ...);
extern int scan_error_n;
	
static int token,ii;
int tokenpos,lineno=1,lineno1=1;
char linebuf[LINEBUFSIZE],linebuf1[LINEBUFSIZE];
static int symbol_lookup(const char *p);
static char *get_string(const char *p);
static YY_BUFFER_STATE bufstate_list[MAX_INCLUDE];
char *fname_list[MAX_INCLUDE+1];
int list_ptr;
static int ln_list[MAX_INCLUDE],ln1_list[MAX_INCLUDE];
	
%}

nl (\r\n|\n|\r)
  
%x COMMENT
%s ISTATE

%%

<INITIAL>.*{nl} {
    (void)strncpy(linebuf,yytext,LINEBUFSIZE);
    (void)strncpy(linebuf1,linebuf,LINEBUFSIZE);
    yyless((int)0);
    BEGIN ISTATE;
    }

<<EOF>> {
	if(iflag) {
		iflag=0;
	} else {
		if(--list_ptr<0) yyterminate();
		else {
			free(fname_list[list_ptr+1]);
			yy_delete_buffer(YY_CURRENT_BUFFER);
			yy_switch_to_buffer(bufstate_list[list_ptr]);
			lineno=ln_list[list_ptr];
			lineno1=ln1_list[list_ptr];
			linebuf[0]=linebuf1[0]=0;
		}
	}
}

[\t ] {tokenpos+=yyleng;}

{nl}.* { 
    lineno++;
    tokenpos=0;
    (void)strncpy(linebuf,yytext+1,LINEBUFSIZE);
    yyless((int)1);
    }	

"/*" BEGIN(COMMENT);
<COMMENT>. {tokenpos+=yyleng;}
<COMMENT>{nl} {lineno++; tokenpos=0;}
<COMMENT>"*/" {BEGIN(ISTATE); tokenpos+=yyleng;}

#.*$ ;

\"[^\"\n\r]*({nl}|\") { tokenpos+=yyleng;
	          yylval.string=get_string(yytext+1);
	          ii=strlen(yylval.string);
				 if(ii) yylval.string[ii-1]='\0';
	          return STRING; }

\'[^\'\n\r]*({nl}|\') { tokenpos+=yyleng;
	          yylval.string=get_string(yytext+1);
	          ii=strlen(yylval.string);
				 if(ii) yylval.string[ii-1]='\0';
	          return STRING; }

[a-zA-Z_][a-zA-Z0-9_]* { token=symbol_lookup(yytext);
				    tokenpos+=yyleng;
	                   if(token==STRING) yylval.string=get_string(yytext);
				    return token; }

[+-]?([0-9]*\.[0-9]+)|(([0-9]+|([0-9]*\.[0-9]+))E[-+]?[0-9]+) {
         tokenpos+=yyleng; yylval.rvalue=atof(yytext); return REAL; }

[+-]?[0-9]+ { tokenpos+=yyleng; yylval.value=atoi(yytext); return INTEGER; }

. { tokenpos+=yyleng; return yytext[0]; }

%%

void yy_cleanup(void)
{
	yy_delete_buffer(yy_current_buffer);
}

#ifdef FUNC_NAME
#undef FUNC_NAME
#endif
#define FUNC_NAME "get_string"
static char *get_string(const char *s)
{
	char *s1;
	size_t size;
	
	size=strlen(s)+1;
	if(!(s1=malloc(size))) ABT_FUNC(MMsg);
	(void)memcpy(s1,s,size);
	return s1;
}

static int symbol_lookup(const char *p)
{
	static char *Coms[] = {"RESIDUAL","GENETIC","VARIANCE","POSITION","FREQUENCY",
		"START","ITERATIONS","SAMPLE","FROM","OUTPUT","MEAN","MAP","TOTAL",
		"SEED","FILE","SEEDFILE","TRAIT","LOCI","SET","ESTIMATE","IBD",
		"GROUP","ORDER","MALE","FEMALE","LIMIT","PHENOTYPE",
		"PHENOTYPES","PHENO","GENOTYPE","GENOTYPES","GENO","COUNTS","DUMP",
		"TYPE","ANALYZE","NORMAL","STUDENT_T","HAPLO","HAPLOTYPES","HAPLOTYPE","INCLUDE",
		  "FUNCTION","HALDANE","KOSAMBI",
		  "POLYGENIC","AFFECTED","TIME","VIRTUAL","MARKERS","GRID","COMPRESS","DIR",
		  "PSEUDO","PSEUDOCHROM","PSEUDOCHROMOSOME",(char *)0};
	
	static int Com_token[] = {RESIDUAL,GENETIC,VARIANCE,POSITION,FREQUENCY,
		START,ITERATIONS,SAMPLE,FROM,OUTPUT,MEAN,MAP,TOTAL,
		SEED,SFILE,SEEDFILE,TRAIT,LOCI,SET,ESTIMATE,IBD,
		GROUP,ORDER,MALE,FEMALE,LIMIT,PHENO,PHENO,PHENO,GENO,GENO,GENO,
		COUNTS,DUMP,TYPE,ANALYZE,NORMAL,STUDENT_T,HAPLO,HAPLO,HAPLO,INCLUDE,FUNCTION,
		  HALDANE,KOSAMBI,POLYGENIC,AFFECTED,
		  TIMECOM,VIRTUAL,MARKERS,GRID,COMPRESS,DIR,PSEUDO,PSEUDO,PSEUDO,SYSTEM_VAR,STRING};

	static char *Syst[] = {"NO_OVERDOMINANT","TAU_BETA","TAU_MODE","CENSOR_MODE",
		"DEBUG_LEVEL","LM_RATIO","PEEL_TRACE","BACKUPS",
		"SI_MODE","IBD_OUTPUT","RNG","GENV_OUT","MSCAN_PROBS","CAND_GENE_MODEL",(char *)0};

	int i=0,j;
	
	while(Coms[i])	{
		if(!strcasecmp(Coms[i],p)) break;
		i++;
	}
	if(Com_token[i]==SYSTEM_VAR) {
		i++;
		j=0;
		while(Syst[j]) {
			if(!strcasecmp(Syst[j],p))	{
				yylval.value=j;
				i--;
				break;
			}
			j++;
		}
	}
	return Com_token[i];
}

void include_param_file(char *fname)
{
	FILE *fptr;
	YY_BUFFER_STATE ystate;
	int i,j;
	struct stat st1,st2;
	  
	if(fname) {
		if(!fname[0]) yyerror("Null filename from INCLUDE command\n");
		if(stat(fname,&st1)) {
			(void)fprintf(stderr,"File: %s\n",fname);	
			yyerror("Include file could not be stat()'d\n");
			perror(0);
			exit(EXIT_FAILURE);
		}
		for(i=0;i<=list_ptr;i++) {
			if(stat(fname_list[i],&st2)) {
				(void)fprintf(stderr,"File: %s\n",fname_list[i]);	
				yyerror("Include file could not be stat()'d\n");
				perror(0);
				exit(EXIT_FAILURE);
			}
			if(st1.st_ino==st2.st_ino && st1.st_dev==st2.st_dev) {
				for(j=0;j<=list_ptr;j++) {
					if(j==i) (void)fprintf(stderr,"*%s* -> ",fname_list[j]);
					else (void)fprintf(stderr,"%s -> ",fname_list[j]);
				}
				(void)fprintf(stderr,"%s\n",fname);
				yyerror("Recursive include files!\n");
				exit(EXIT_FAILURE);
			}
		}
		if(list_ptr<MAX_INCLUDE) {
			if(!(fptr=fopen(fname,"r"))) {
				yyerror("Include file could not be opened");
				perror(0);
				exit(EXIT_FAILURE);
			}
			bufstate_list[list_ptr]=YY_CURRENT_BUFFER;
			ln_list[list_ptr]=lineno;
			ln1_list[list_ptr]=lineno1;
			lineno=lineno1=1;
			linebuf[0]=linebuf1[0]=0;
			fname_list[++list_ptr]=fname;
			ystate=yy_create_buffer(fptr,YY_BUF_SIZE);
			yy_switch_to_buffer(ystate);
		} else yyerror("Too many nested includes\n");
	}
}
