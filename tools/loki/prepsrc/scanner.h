#ifndef _SCANNER_H_
#define _SCANNER_H_

struct token_store
{
	YYSTYPE yylval;
	int token,line,line1;
};

extern struct token_store *loop_stack;

extern int yylex(void);
extern int yyparse(void);
extern void include_control_file(char *);
extern char *fname_list[MAX_INCLUDE+1];
extern int list_ptr,iflag;

struct marker_info {
	struct marker_info *next;
	struct var_element *element;
	double pos[3];
	int pos_set[3];
};

#endif
