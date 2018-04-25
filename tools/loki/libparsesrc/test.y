%{
#undef YYLSP_NEEDED
#undef YYLEX_PARAM

%}

%union {
	int i;
	double x;
}

%token DIGIT
  
%%

clause: /* empty */ 
     | atom
     | clause ',' atom
	  | int '(' clause ')'
	  | '(' clause ')'
     ;

atom: int
     | 'x'
	  | int 'x'
	  ;
	  
int: DIGIT
     | int DIGIT
	  ;
	  
%%
