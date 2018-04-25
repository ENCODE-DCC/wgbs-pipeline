%{
/****************************************************************************
 *                                                                          *
 *     Loki - Programs for genetic analysis of complex traits using MCMC    *
 *                                                                          *
 *             Simon Heath - University of Washington                       *
 *                                                                          *
 *                       March 1997                                         *
 *                                                                          *
 * control_lex.l:                                                           *
 *                                                                          *
 * flex (NOT lex!) source for control file lexer.                           *
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
#include <ctype.h>
#include <sys/stat.h>
	
#include "y.tab.h"
#include "scan.h"
#include "scanner.h"
#include "utils.h"

#define YY_DECL int yy_orig_lex YY_PROTO(( void ))
#define MAX_STRING 512
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
	
int lineno=1,lineno1=1,at_file,at_model,at_use,begin_comm;
char linebuf[LINEBUFSIZE],linebuf1[LINEBUFSIZE];
static char *get_string(char *s,int size);
static YY_BUFFER_STATE bufstate_list[MAX_INCLUDE];
char *fname_list[MAX_INCLUDE+1];
int list_ptr;

static int tokenpos,stringbuf_ptr,ff_flag,ln_list[MAX_INCLUDE],ln1_list[MAX_INCLUDE];
static char stringbuf[MAX_STRING];
static void process_env_var(char *ptr);
static void append_to_sbuf(int *sptr,char c);

%}

nl (\r\n|\n|\r)
  
%x COMMENT
%x STR 
%s ATFORMAT
%s ISTATE

%%

%{
	int token;
%}

<<EOF>> {
	if(iflag) {
		iflag=0;
		return BREAK;
	}
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

{nl}?.*(\$[a-zA-Z][a-zA-Z0-9_]*)+.*$ { process_env_var(yytext); }

<INITIAL>.*{nl} {
    (void)strncpy(linebuf,yytext,LINEBUFSIZE);
    (void)strncpy(linebuf1,linebuf,LINEBUFSIZE);
    yyless((int)0);
    BEGIN ISTATE;
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

[\"\']   { stringbuf_ptr=0;
       BEGIN(STR);
 	  tokenpos+=yyleng; }

<STR>[\"\'] {
           BEGIN(ISTATE);
	      append_to_sbuf(&stringbuf_ptr,'\0');
		 yylval.string=get_string(stringbuf,stringbuf_ptr);
		 tokenpos+=yyleng;
		 return STRING; 
	   }

<STR>{nl} { BEGIN(ISTATE);
	     append_to_sbuf(&stringbuf_ptr,'\0');
		print_scan_warn("Line %d: Warning - Unterminanted string\n",lineno);
		tokenpos+=yyleng;
		yylval.string=get_string(stringbuf,stringbuf_ptr);
		return STRING; }

<STR>\\[0-7]{1,3} { unsigned int result;
                    (void) sscanf(yytext+1,"%uo",&result);
				tokenpos+=yyleng;
				if(result>0xff) print_scan_err("Line %d: Error - Octal constant too large\n",lineno);
				append_to_sbuf(&stringbuf_ptr,(char)result); }

<STR>\\x[a-fA-F0-9]{1,2} { unsigned int result;
                    (void) sscanf(yytext+1,"%ux",&result);
				tokenpos+=yyleng;
 				append_to_sbuf(&stringbuf_ptr,(char)result); }

<STR>\\[0-9]+ { tokenpos+=yyleng;
                print_scan_err("Line %d: Error - Syntax error in escape sequence\n",lineno); }

<STR>\\a { tokenpos+=yyleng;
	     append_to_sbuf(&stringbuf_ptr,'\a'); }

<STR>\\b { tokenpos+=yyleng;
 	     append_to_sbuf(&stringbuf_ptr,'\b'); }

<STR>\\f { tokenpos+=yyleng;
	     append_to_sbuf(&stringbuf_ptr,'\f'); }

<STR>\\n { tokenpos+=yyleng;
	     append_to_sbuf(&stringbuf_ptr,'\n'); }

<STR>\\r { tokenpos+=yyleng;
	     append_to_sbuf(&stringbuf_ptr,'\r'); }

<STR>\\t { tokenpos+=yyleng;
	     append_to_sbuf(&stringbuf_ptr,'\t'); }

<STR>\\v { tokenpos+=yyleng;
	     append_to_sbuf(&stringbuf_ptr,'\v'); }

<STR>\\(.|\n) { tokenpos++;
	     append_to_sbuf(&stringbuf_ptr,'\a'); }

<STR>[^\\\n\r\"\']+ { char *ss=yytext;
                  tokenpos+=yyleng;
			   while(*ss) append_to_sbuf(&stringbuf_ptr,*ss++); }

\[ { if(at_file) BEGIN ATFORMAT;
     tokenpos+=yyleng;
     return '['; }

<ATFORMAT>\] { if(at_file) BEGIN ISTATE;
              tokenpos+=yyleng;
              return ']'; }

<ATFORMAT>[xX] { tokenpos+=yyleng; return('x'); }

<ATFORMAT>[yY] { tokenpos+=yyleng; return('y'); }

"<=" { return LEQSYMBOL; }

">=" { return GEQSYMBOL; }

"!="|"<>" { return NEQSYMBOL; }

"==" { return '='; }

\||[Oo][Rr]|\|\| { return ORSYMBOL; }

&|[Aa][Nn][Dd]|&& { return ANDSYMBOL; }

!|[Nn][Oo][Tt] { return NOTSYMBOL; }

<ATFORMAT>[Ss][Kk][Ii][Pp][ \t]*= { return SKIPFORMAT; }
<ISTATE>[Ss][Kk][Ii][Pp][ \t]*= { if(at_file) return SKIPFORMAT;
								  if(!ff_flag) {yyless(0); token=BREAK;}			 
                          else token=SKIPFORMAT;
	                       ff_flag^=1;
	                       return token; }
<ATFORMAT>[Rr][Ss][ \t]*= { return RSFORMAT; }
<ISTATE>[Rr][Ss][ \t]*= { if(at_file) return RSFORMAT;
								  if(!ff_flag) {yyless(0); token=BREAK;}
                          else token=RSFORMAT;
	                       ff_flag^=1;
	                       return token; }
<ATFORMAT>[Gg][Ss][ \t]*= { return GSFORMAT; }
<ISTATE>[Gg][Ss][ \t]*= { if(at_file) return GSFORMAT;
								  if(!ff_flag) {yyless(0); token=BREAK;}
                          else token=GSFORMAT;
	                       ff_flag^=1;
	                       return token; }

<ATFORMAT>[Ff][Ss][ \t]*= { return FSFORMAT; }
<ISTATE>[Ff][Ss][ \t]*= { if(at_file) return FSFORMAT;
								  if(!ff_flag) {yyless(0); token=BREAK;}
                          else token=FSFORMAT;
	                       ff_flag^=1;
	                       return token; }

[a-zA-Z_][a-zA-Z0-9_]*/([\t ]*=)|([\t ]*\([^\n\"']*\)[\t ]*=) {
                               token=symbol_lookup(yytext,1);
				           if(token==BREAK) yyless(0);
				           else
						 {
						   tokenpos+=yyleng;
						   if(!at_model && !at_use && (token==VARIABLE || token==ARRAY_VAR)) token=ASSIGN;
						 }
						 if((!at_model && !at_use) || token!=BREAK) return token;
					    }

[_a-zA-Z][a-zA-Z0-9_]* {
	             token=symbol_lookup(yytext,0);
				    if(token==BREAK) yyless(0);
				    else tokenpos+=yyleng; 
				    return token; }

([0-9]*\.[0-9]+)|(([0-9]+|([0-9]*\.[0-9]+))E[-+]?[0-9]+) {
         tokenpos+=yyleng; yylval.rvalue=atof(yytext); return REAL; }

[0-9]+ { tokenpos+=yyleng; yylval.value=atoi(yytext); return INTEGER; }

. { tokenpos+=yyleng; return yytext[0]; }

%%
/*"*/

int yywrap(void) 
{
	return 1;
}
	
void yyerror(char *s)
{
     int i;
	
	if(loop_record== -1)	{
		print_scan_err("File %s, line %d: %s\n",fname_list[list_ptr],lineno,s);
		loop_record= -2;
	}
	else if(loop_record>=0)	{
		print_scan_err("File %s, line %d: %s\n%s\n",fname_list[list_ptr],lineno,s,linebuf);
		if(scan_error_n<=max_scan_errors && linebuf[0]) {
			for(i=1;i<tokenpos;i++) (void)putc('-',stderr);
			(void)fputs("^\n",stderr);
		}
	}
	scan_error|=SCAN_ERR;
}

void yyerror1(char *s)
{
	if(loop_record== -1)	{
		print_scan_err("File %s, line %d: %s\n",fname_list[list_ptr],lineno1,s);
		loop_record= -2;
	} else if(loop_record>=0) print_scan_err("File %s, line %d: %s\n%s\n",fname_list[list_ptr],lineno1,s,linebuf1);
	scan_error|=SCAN_ERR;
}

void include_control_file(char *fname)
{
	FILE *fptr;
	YY_BUFFER_STATE ystate;
	int i,j;
	struct stat st1,st2;
	  
	if(fname) {
		if(!fname[0]) yyerror1("Null filename from INCLUDE command\n");
		if(stat(fname,&st1)) {
			(void)fprintf(stderr,"File: %s\n",fname);	
			yyerror1("Include file could not be stat()'d\n");
			perror(0);
			exit(EXIT_FAILURE);
		}
		for(i=0;i<=list_ptr;i++) {
			if(stat(fname_list[i],&st2)) {
				(void)fprintf(stderr,"File: %s\n",fname_list[i]);	
				yyerror1("Include file could not be stat()'d\n");
				perror(0);
				exit(EXIT_FAILURE);
			}
			if(st1.st_ino==st2.st_ino && st1.st_dev==st2.st_dev) {
				for(j=0;j<=list_ptr;j++) {
					if(j==i) (void)fprintf(stderr,"*%s* -> ",fname_list[j]);
					else (void)fprintf(stderr,"%s -> ",fname_list[j]);
				}
				(void)fprintf(stderr,"%s\n",fname);
				yyerror1("Recursive include files!\n");
				exit(EXIT_FAILURE);
			}
		}
		if(list_ptr<MAX_INCLUDE) {
			if(!(fptr=fopen(fname,"r"))) {
				yyerror1("Include file could not be opened");
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
  
static void append_to_sbuf(int *sptr,char c)
{
	if(*sptr>=MAX_STRING) print_scan_err("Line %d: Maximum string size exceeded",lineno);
	stringbuf[(*sptr)++]=c;
}

#ifdef FUNC_NAME
#undef FUNC_NAME
#endif
#define FUNC_NAME "get_string"
static char *get_string(char *s,int size)
{
	char *s1;
	if(!(s1=malloc(size))) ABT_FUNC(MMsg);
	(void)memcpy(s1,s,size);
	return s1;
}

#ifdef FUNC_NAME
#undef FUNC_NAME
#endif
#define FUNC_NAME "yylex"
int yylex(void)
{
	int tok,i;
	
	if(loop_record==1) {
		tok=yy_orig_lex();
		if(!loop_stack) {
			if(!(loop_stack=malloc(sizeof(struct token_store)*loop_stack_size))) ABT_FUNC(MMsg);
			loop_main_ptr=0;
		} else if(loop_main_ptr==loop_stack_size)	{
			loop_stack_size*=2;
			if(!(loop_stack=realloc(loop_stack,sizeof(struct token_store)*loop_stack_size))) ABT_FUNC(MMsg);
		}
		loop_stack[loop_main_ptr].line=lineno;
		loop_stack[loop_main_ptr].line1=lineno1;
		loop_stack[loop_main_ptr].yylval=yylval;
		if(tok==STRING) yylval.string=strdup(yylval.string);
		loop_stack[loop_main_ptr++].token=tok;
		return tok;	   
	} else if(loop_record<0) {
		if(in_loopclause== -1) {
			in_loopclause= -2;
			return LOOP_CLAUSE_START;
		}
		yylval=loop_stack[loop_main_ptr].yylval;
		lineno=loop_stack[loop_main_ptr].line;
		lineno1=loop_stack[loop_main_ptr].line1;
		tok=loop_stack[loop_main_ptr++].token;
		if(tok==STRING) yylval.string=strdup(yylval.string);
		if(in_loopclause) {
			if(tok==',' && loop_stack[loop_main_ptr].token==BREAK) {
				loop_clause_element->arg.value+=loop_clause_step;
				if(loop_clause_element->arg.value==loop_clause_end) i=0;
				else {
					i=(loop_clause_element->arg.value>loop_clause_end);
					if(loop_clause_step<0) i=1-i;
				}
				if(i) {
					loop_record=loop_stat[--loop_level];
					if(!loop_record) loop_main_ptr=loop_level?loop_ptr[loop_level-1]:0;
					else loop_main_ptr=loop_clause_ptr;
					in_loopclause=0;
					return LOOP_CLAUSE_END;
				} else loop_main_ptr=loop_ptr[loop_level-1];
			}
		}
		return tok;
	}
	/* Dummy call to shut up compiler warnings about yy_flex_realloc() not being used */
	if(0) (void)yy_flex_realloc(0,0);
	return yy_orig_lex();
}

void yy_cleanup(void)
{
	yy_delete_buffer(yy_current_buffer);
	if(loop_stack) free(loop_stack);
}

#ifdef FUNC_NAME
#undef FUNC_NAME
#endif
#define FUNC_NAME "process_env_var"
static void process_env_var(char *ptr)
{
	int i=0,j,k,size=512,brack;
	char *s,old,*stringptr;
	
	if(!(stringptr=malloc(size))) ABT_FUNC(MMsg);
	for(;;) {
		while(*ptr && *ptr!='$') {
			if(i>=size) {
				size+=size;
				if(!(stringptr=realloc(stringptr,size))) ABT_FUNC(MMsg);
			}
			stringptr[i++] = *(ptr++);
		}
		if(!(*ptr)) break;
		ptr++;
		j=0;
		brack=0;
		if(*ptr=='(') {
			ptr++;
			while(ptr[j] && ptr[j]!=')') j++;
			if(ptr[j]==')') brack=1;
		} else {
			while(isalnum((int)ptr[j]) || ptr[j]=='_') j++;
		}
		old=ptr[j];
		ptr[j]=0;
		if((s=getenv(ptr))) {
			if((k=strlen(s)))	{
				if(k+i>=size) {
					while(k+i>size) size+=size;
					if(!(stringptr=realloc(stringptr,size))) ABT_FUNC(MMsg);
				}
				(void)memcpy(stringptr+i,s,k);
				i+=k;
			}
		}
		ptr[j]=old;
		ptr+=j+brack;
	}
	while(i) unput(stringptr[--i]);
	free(stringptr);
}
