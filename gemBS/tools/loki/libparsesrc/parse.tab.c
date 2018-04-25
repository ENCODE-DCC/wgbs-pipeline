
/* A Bison parser, made by GNU Bison 2.4.1.  */

/* Skeleton implementation for Bison's Yacc-like parsers in C
   
      Copyright (C) 1984, 1989, 1990, 2000, 2001, 2002, 2003, 2004, 2005, 2006
   Free Software Foundation, Inc.
   
   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.
   
   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
   
   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.  */

/* As a special exception, you may create a larger work that contains
   part or all of the Bison parser skeleton and distribute that work
   under terms of your choice, so long as that work isn't itself a
   parser generator using the skeleton or a modified version thereof
   as a parser skeleton.  Alternatively, if you modify or redistribute
   the parser skeleton itself, you may (at your option) remove this
   special exception, which will cause the skeleton and the resulting
   Bison output files to be licensed under the GNU General Public
   License without this special exception.
   
   This special exception was added by the Free Software Foundation in
   version 2.2 of Bison.  */

/* C LALR(1) parser skeleton written by Richard Stallman, by
   simplifying the original so-called "semantic" parser.  */

/* All symbols defined below should begin with yy or YY, to avoid
   infringing on user name space.  This should be done even for local
   variables, as they might otherwise be expanded by user macros.
   There are some unavoidable exceptions within include files to
   define necessary library symbols; they are noted "INFRINGES ON
   USER NAME SPACE" below.  */

/* Identify Bison output.  */
#define YYBISON 1

/* Bison version.  */
#define YYBISON_VERSION "2.4.1"

/* Skeleton name.  */
#define YYSKELETON_NAME "yacc.c"

/* Pure parsers.  */
#define YYPURE 1

/* Push parsers.  */
#define YYPUSH 0

/* Pull parsers.  */
#define YYPULL 1

/* Using locations.  */
#define YYLSP_NEEDED 0

/* Substitute the variable and function names.  */
#define yyparse         parparse
#define yylex           parlex
#define yyerror         parerror
#define yylval          parlval
#define yychar          parchar
#define yydebug         pardebug
#define yynerrs         parnerrs


/* Copy the first part of user declarations.  */

/* Line 189 of yacc.c  */
#line 2 "parse.y"

/****************************************************************************
 *                                                                          *
 *     Loki - Programs for genetic analysis of complex traits using MCMC    *
 *                                                                          *
 *                     Simon Heath - CNG, Evry                              *
 *                                                                          *
 *                           June 2003                                      *
 *                                                                          *
 * parse.y:                                                                 *
 *                                                                          *
 * Basic control file parser                                                *
 *                                                                          *
 * Copyright (C) Simon C. Heath 2003                                        *
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
#include <stdio.h>
#include <stdarg.h>
#include <ctype.h>
#include <math.h>
#include <assert.h>
#include <sys/stat.h>
	
#include "string_utils.h"
#include "bin_tree.h"
#include "parser.h"
#include "parse.tab.h"
#include "lk_malloc.h"

#undef YYLSP_NEEDED
#undef YYLEX_PARAM

#define BUFSIZE 640
#define LOOKAHEAD 128

static int parerror(char *);
static int parlex(yystype *);
static char *key_err="Can not use keyword '%s' as a variable\n";
static parse_handler *handle;
static int in_clause;
	
#define ALPHA_BIT 1
#define DIGIT_BIT 2
#define E_BIT 4
#define USCORE_BIT 8
#define PERIOD_BIT 16
#define PM_BIT 32
#define MATH_BIT 64
#define REL_BIT 128
#define LBRACK_BIT 256
#define RBRACK_BIT 512
#define COMMA_BIT 1024
#define X_BIT 2048
#define TOKEN_BITS (PM_BIT-1)
#define FORMAT_BITS (DIGIT_BIT|X_BIT|LBRACK_BIT|RBRACK_BIT|COMMA_BIT)

struct lex_info {
	struct lex_info *next,*last;
	int line_no;
	int col_no,col_no1;
	int eof;
	char *name;
	int ptr,ptr1;
	int last_token;
	int flag;
	yystype last_lval;
	FILE *fptr;
	char buf[BUFSIZE];
};

static void include_file(string *);
static void copy_filepos(filepos *,filepos *);

static struct lex_info *lex_info;
static struct bin_node *root_var;
static int lex_tab[256];
static char **key_tab,**command_tab;
static func_def *func_list;
static struct parse_term *vtmp;
static int typeI_kludge;
	


/* Line 189 of yacc.c  */
#line 173 "y.tab.c"

/* Enabling traces.  */
#ifndef YYDEBUG
# define YYDEBUG 0
#endif

/* Enabling verbose error messages.  */
#ifdef YYERROR_VERBOSE
# undef YYERROR_VERBOSE
# define YYERROR_VERBOSE 1
#else
# define YYERROR_VERBOSE 0
#endif

/* Enabling the token table.  */
#ifndef YYTOKEN_TABLE
# define YYTOKEN_TABLE 0
#endif


/* Tokens.  */
#ifndef YYTOKENTYPE
# define YYTOKENTYPE
   /* Put the tokens into the symbol table, so that GDB and other debuggers
      know about them.  */
   enum yytokentype {
     ALIAS = 258,
     FOR = 259,
     IF = 260,
     ELSE = 261,
     DO = 262,
     WHILE = 263,
     INCLUDE = 264,
     NE_TOKEN = 265,
     LEQ_TOKEN = 266,
     GEQ_TOKEN = 267,
     CMP_TOKEN = 268,
     VARIABLE = 269,
     EXPR_OP = 270,
     INF = 271,
     XNAN = 272,
     DUMMY = 273,
     DOUBLE_DOT = 274,
     EMPTY_VAR = 275,
     EXPO = 276,
     CLAUSE_X = 277,
     EQ_TOKEN = 278,
     SUB_EXPR = 279,
     STRING = 280,
     INTEGER = 281,
     CLAUSE_INT = 282,
     PREFIX = 283,
     POSTFIX = 284,
     MATH_SHORT = 285,
     KEYWORD = 286,
     KEYWORD_B = 287,
     REAL = 288,
     TOKEN = 289,
     WORD_B = 290,
     COMMAND = 291,
     FUNCTION = 292,
     LOWER_THAN_ELSE = 293,
     UMINUS = 294
   };
#endif
/* Tokens.  */
#define ALIAS 258
#define FOR 259
#define IF 260
#define ELSE 261
#define DO 262
#define WHILE 263
#define INCLUDE 264
#define NE_TOKEN 265
#define LEQ_TOKEN 266
#define GEQ_TOKEN 267
#define CMP_TOKEN 268
#define VARIABLE 269
#define EXPR_OP 270
#define INF 271
#define XNAN 272
#define DUMMY 273
#define DOUBLE_DOT 274
#define EMPTY_VAR 275
#define EXPO 276
#define CLAUSE_X 277
#define EQ_TOKEN 278
#define SUB_EXPR 279
#define STRING 280
#define INTEGER 281
#define CLAUSE_INT 282
#define PREFIX 283
#define POSTFIX 284
#define MATH_SHORT 285
#define KEYWORD 286
#define KEYWORD_B 287
#define REAL 288
#define TOKEN 289
#define WORD_B 290
#define COMMAND 291
#define FUNCTION 292
#define LOWER_THAN_ELSE 293
#define UMINUS 294




#if ! defined YYSTYPE && ! defined YYSTYPE_IS_DECLARED
typedef union YYSTYPE
{

/* Line 214 of yacc.c  */
#line 93 "parse.y"

	int i;
	double x;
	string *str;
	struct parse_term *term;
	struct parse_var *vv;
	struct parse_var_list *vl;
	struct format_clause *fc;
	struct parse_clause *pc;
	func_def *f;



/* Line 214 of yacc.c  */
#line 301 "y.tab.c"
} YYSTYPE;
# define YYSTYPE_IS_TRIVIAL 1
# define yystype YYSTYPE /* obsolescent; will be withdrawn */
# define YYSTYPE_IS_DECLARED 1
#endif


/* Copy the second part of user declarations.  */


/* Line 264 of yacc.c  */
#line 313 "y.tab.c"

#ifdef short
# undef short
#endif

#ifdef YYTYPE_UINT8
typedef YYTYPE_UINT8 yytype_uint8;
#else
typedef unsigned char yytype_uint8;
#endif

#ifdef YYTYPE_INT8
typedef YYTYPE_INT8 yytype_int8;
#elif (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
typedef signed char yytype_int8;
#else
typedef short int yytype_int8;
#endif

#ifdef YYTYPE_UINT16
typedef YYTYPE_UINT16 yytype_uint16;
#else
typedef unsigned short int yytype_uint16;
#endif

#ifdef YYTYPE_INT16
typedef YYTYPE_INT16 yytype_int16;
#else
typedef short int yytype_int16;
#endif

#ifndef YYSIZE_T
# ifdef __SIZE_TYPE__
#  define YYSIZE_T __SIZE_TYPE__
# elif defined size_t
#  define YYSIZE_T size_t
# elif ! defined YYSIZE_T && (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
#  include <stddef.h> /* INFRINGES ON USER NAME SPACE */
#  define YYSIZE_T size_t
# else
#  define YYSIZE_T unsigned int
# endif
#endif

#define YYSIZE_MAXIMUM ((YYSIZE_T) -1)

#ifndef YY_
# if YYENABLE_NLS
#  if ENABLE_NLS
#   include <libintl.h> /* INFRINGES ON USER NAME SPACE */
#   define YY_(msgid) dgettext ("bison-runtime", msgid)
#  endif
# endif
# ifndef YY_
#  define YY_(msgid) msgid
# endif
#endif

/* Suppress unused-variable warnings by "using" E.  */
#if ! defined lint || defined __GNUC__
# define YYUSE(e) ((void) (e))
#else
# define YYUSE(e) /* empty */
#endif

/* Identity function, used to suppress warnings about constant conditions.  */
#ifndef lint
# define YYID(n) (n)
#else
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static int
YYID (int yyi)
#else
static int
YYID (yyi)
    int yyi;
#endif
{
  return yyi;
}
#endif

#if ! defined yyoverflow || YYERROR_VERBOSE

/* The parser invokes alloca or malloc; define the necessary symbols.  */

# ifdef YYSTACK_USE_ALLOCA
#  if YYSTACK_USE_ALLOCA
#   ifdef __GNUC__
#    define YYSTACK_ALLOC __builtin_alloca
#   elif defined __BUILTIN_VA_ARG_INCR
#    include <alloca.h> /* INFRINGES ON USER NAME SPACE */
#   elif defined _AIX
#    define YYSTACK_ALLOC __alloca
#   elif defined _MSC_VER
#    include <malloc.h> /* INFRINGES ON USER NAME SPACE */
#    define alloca _alloca
#   else
#    define YYSTACK_ALLOC alloca
#    if ! defined _ALLOCA_H && ! defined _STDLIB_H && (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
#     include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
#     ifndef _STDLIB_H
#      define _STDLIB_H 1
#     endif
#    endif
#   endif
#  endif
# endif

# ifdef YYSTACK_ALLOC
   /* Pacify GCC's `empty if-body' warning.  */
#  define YYSTACK_FREE(Ptr) do { /* empty */; } while (YYID (0))
#  ifndef YYSTACK_ALLOC_MAXIMUM
    /* The OS might guarantee only one guard page at the bottom of the stack,
       and a page size can be as small as 4096 bytes.  So we cannot safely
       invoke alloca (N) if N exceeds 4096.  Use a slightly smaller number
       to allow for a few compiler-allocated temporary stack slots.  */
#   define YYSTACK_ALLOC_MAXIMUM 4032 /* reasonable circa 2006 */
#  endif
# else
#  define YYSTACK_ALLOC YYMALLOC
#  define YYSTACK_FREE YYFREE
#  ifndef YYSTACK_ALLOC_MAXIMUM
#   define YYSTACK_ALLOC_MAXIMUM YYSIZE_MAXIMUM
#  endif
#  if (defined __cplusplus && ! defined _STDLIB_H \
       && ! ((defined YYMALLOC || defined malloc) \
	     && (defined YYFREE || defined free)))
#   include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
#   ifndef _STDLIB_H
#    define _STDLIB_H 1
#   endif
#  endif
#  ifndef YYMALLOC
#   define YYMALLOC malloc
#   if ! defined malloc && ! defined _STDLIB_H && (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
void *malloc (YYSIZE_T); /* INFRINGES ON USER NAME SPACE */
#   endif
#  endif
#  ifndef YYFREE
#   define YYFREE free
#   if ! defined free && ! defined _STDLIB_H && (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
void free (void *); /* INFRINGES ON USER NAME SPACE */
#   endif
#  endif
# endif
#endif /* ! defined yyoverflow || YYERROR_VERBOSE */


#if (! defined yyoverflow \
     && (! defined __cplusplus \
	 || (defined YYSTYPE_IS_TRIVIAL && YYSTYPE_IS_TRIVIAL)))

/* A type that is properly aligned for any stack member.  */
union yyalloc
{
  yytype_int16 yyss_alloc;
  YYSTYPE yyvs_alloc;
};

/* The size of the maximum gap between one aligned stack and the next.  */
# define YYSTACK_GAP_MAXIMUM (sizeof (union yyalloc) - 1)

/* The size of an array large to enough to hold all stacks, each with
   N elements.  */
# define YYSTACK_BYTES(N) \
     ((N) * (sizeof (yytype_int16) + sizeof (YYSTYPE)) \
      + YYSTACK_GAP_MAXIMUM)

/* Copy COUNT objects from FROM to TO.  The source and destination do
   not overlap.  */
# ifndef YYCOPY
#  if defined __GNUC__ && 1 < __GNUC__
#   define YYCOPY(To, From, Count) \
      __builtin_memcpy (To, From, (Count) * sizeof (*(From)))
#  else
#   define YYCOPY(To, From, Count)		\
      do					\
	{					\
	  YYSIZE_T yyi;				\
	  for (yyi = 0; yyi < (Count); yyi++)	\
	    (To)[yyi] = (From)[yyi];		\
	}					\
      while (YYID (0))
#  endif
# endif

/* Relocate STACK from its old location to the new one.  The
   local variables YYSIZE and YYSTACKSIZE give the old and new number of
   elements in the stack, and YYPTR gives the new location of the
   stack.  Advance YYPTR to a properly aligned location for the next
   stack.  */
# define YYSTACK_RELOCATE(Stack_alloc, Stack)				\
    do									\
      {									\
	YYSIZE_T yynewbytes;						\
	YYCOPY (&yyptr->Stack_alloc, Stack, yysize);			\
	Stack = &yyptr->Stack_alloc;					\
	yynewbytes = yystacksize * sizeof (*Stack) + YYSTACK_GAP_MAXIMUM; \
	yyptr += yynewbytes / sizeof (*yyptr);				\
      }									\
    while (YYID (0))

#endif

/* YYFINAL -- State number of the termination state.  */
#define YYFINAL  98
/* YYLAST -- Last index in YYTABLE.  */
#define YYLAST   1333

/* YYNTOKENS -- Number of terminals.  */
#define YYNTOKENS  59
/* YYNNTS -- Number of nonterminals.  */
#define YYNNTS  57
/* YYNRULES -- Number of rules.  */
#define YYNRULES  196
/* YYNRULES -- Number of states.  */
#define YYNSTATES  358

/* YYTRANSLATE(YYLEX) -- Bison symbol number corresponding to YYLEX.  */
#define YYUNDEFTOK  2
#define YYMAXUTOK   294

#define YYTRANSLATE(YYX)						\
  ((unsigned int) (YYX) <= YYMAXUTOK ? yytranslate[YYX] : YYUNDEFTOK)

/* YYTRANSLATE[YYLEX] -- Bison symbol number corresponding to YYLEX.  */
static const yytype_uint8 yytranslate[] =
{
       0,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,    50,     2,     2,     2,     2,    42,     2,
      56,    55,    48,    46,    39,    45,    47,    49,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,    54,
      43,    40,    44,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,    58,     2,    57,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,    52,    41,    53,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     1,     2,     3,     4,
       5,     6,     7,     8,     9,    10,    11,    12,    13,    14,
      15,    16,    17,    18,    19,    20,    21,    22,    23,    24,
      25,    26,    27,    28,    29,    30,    31,    32,    33,    34,
      35,    36,    37,    38,    51
};

#if YYDEBUG
/* YYPRHS[YYN] -- Index of the first RHS symbol of rule number YYN in
   YYRHS.  */
static const yytype_uint16 yyprhs[] =
{
       0,     0,     3,     5,     7,    10,    12,    14,    18,    20,
      23,    25,    27,    29,    31,    33,    35,    37,    39,    41,
      43,    45,    48,    53,    56,    60,    66,    70,    73,    77,
      82,    88,    90,    94,    96,   100,   104,   109,   113,   118,
     120,   124,   127,   130,   134,   139,   146,   151,   157,   161,
     162,   167,   168,   175,   176,   182,   183,   191,   194,   198,
     200,   202,   204,   206,   208,   210,   212,   214,   216,   218,
     220,   222,   224,   226,   228,   229,   233,   238,   242,   244,
     246,   248,   252,   254,   256,   258,   262,   266,   269,   271,
     281,   287,   295,   302,   308,   309,   311,   315,   317,   321,
     323,   325,   329,   332,   335,   337,   339,   343,   345,   349,
     351,   355,   356,   361,   370,   375,   383,   387,   389,   393,
     397,   399,   405,   409,   411,   414,   416,   418,   421,   424,
     426,   429,   432,   436,   440,   442,   446,   450,   454,   458,
     462,   466,   470,   474,   478,   482,   486,   490,   494,   498,
     501,   505,   510,   512,   514,   516,   520,   524,   526,   531,
     533,   535,   537,   541,   543,   547,   551,   555,   559,   563,
     567,   571,   575,   579,   583,   587,   591,   595,   599,   603,
     607,   611,   614,   617,   621,   623,   625,   627,   629,   631,
     634,   637,   642,   644,   646,   648,   650
};

/* YYRHS -- A `-1'-separated list of the rules' RHS.  */
static const yytype_int8 yyrhs[] =
{
      60,     0,    -1,    62,    -1,    61,    -1,    60,    62,    -1,
      90,    -1,    64,    -1,    52,    63,    53,    -1,    64,    -1,
      63,    64,    -1,    65,    -1,    66,    -1,    68,    -1,    69,
      -1,    70,    -1,    71,    -1,    85,    -1,    86,    -1,    87,
      -1,    88,    -1,    54,    -1,    36,   105,    -1,    36,   105,
      39,   104,    -1,    36,   104,    -1,    36,    83,   105,    -1,
      36,    83,   105,    39,   104,    -1,    36,    83,   104,    -1,
      36,    18,    -1,    36,    18,    83,    -1,    36,    18,    39,
     104,    -1,    36,    18,    83,    39,   104,    -1,    98,    -1,
      98,    39,   104,    -1,    97,    -1,    97,    39,   104,    -1,
      36,    95,    91,    -1,    31,    36,    95,    91,    -1,    95,
      40,   112,    -1,    36,    95,    40,   112,    -1,    67,    -1,
      67,    39,   110,    -1,    95,    29,    -1,    28,    95,    -1,
      95,    30,   112,    -1,    36,    32,   112,    55,    -1,    31,
     104,    31,    56,   112,    55,    -1,    31,    32,   112,    55,
      -1,    32,   112,    55,    31,   104,    -1,    36,    31,   104,
      -1,    -1,    31,    40,    72,   112,    -1,    -1,    32,   112,
      55,    40,    73,   112,    -1,    -1,    36,    31,    40,    74,
     112,    -1,    -1,    36,    32,   112,    55,    40,    75,   112,
      -1,    34,     1,    -1,    79,    76,   112,    -1,    46,    -1,
      48,    -1,    49,    -1,    47,    -1,    45,    -1,    43,    -1,
      44,    -1,    11,    -1,    12,    -1,    10,    -1,    42,    -1,
      41,    -1,    21,    -1,    27,    -1,    22,    -1,    -1,    27,
      22,    78,    -1,    27,    56,    80,    55,    -1,    56,    80,
      55,    -1,    34,    -1,    36,    -1,    77,    -1,    80,    39,
      77,    -1,   111,    -1,    80,    -1,    81,    -1,    82,    54,
      81,    -1,    84,    82,    57,    -1,    84,    57,    -1,    58,
      -1,     4,    56,    89,    54,    89,    54,    89,    55,    62,
      -1,     5,    56,   112,    55,    62,    -1,     5,    56,   112,
      55,    62,     6,    62,    -1,     7,    62,     8,    56,   112,
      55,    -1,     8,    56,   112,    55,    62,    -1,    -1,   112,
      -1,     9,   107,    18,    -1,    92,    -1,    91,    39,    92,
      -1,   114,    -1,   112,    -1,   112,    19,   112,    -1,    19,
     112,    -1,   112,    19,    -1,    27,    -1,    93,    -1,    94,
      39,    93,    -1,    79,    -1,    35,    94,    55,    -1,    79,
      -1,    35,    94,    55,    -1,    -1,    35,    18,    94,    55,
      -1,    56,    18,   100,    39,    79,    40,   111,    55,    -1,
      56,    18,   100,    55,    -1,    56,   100,    39,    79,    40,
     111,    55,    -1,    56,   100,    55,    -1,    96,    -1,   100,
      39,    96,    -1,    32,   112,    55,    -1,    31,    -1,   100,
      39,    32,   112,    55,    -1,   100,    39,    31,    -1,    95,
      -1,    95,    83,    -1,   101,    -1,    99,    -1,    39,    39,
      -1,   103,    39,    -1,   102,    -1,    39,   102,    -1,   103,
     102,    -1,   104,    39,   102,    -1,   104,   103,   102,    -1,
     114,    -1,   115,    46,   112,    -1,   115,    47,   112,    -1,
     115,    45,   112,    -1,   115,    49,   112,    -1,   115,    48,
     112,    -1,   115,    43,   112,    -1,   115,    44,   112,    -1,
     115,    42,   112,    -1,   115,    41,   112,    -1,   115,    23,
     112,    -1,   115,    10,   112,    -1,   115,    11,   112,    -1,
     115,    12,   112,    -1,   115,    21,   112,    -1,    50,   112,
      -1,    56,   105,    55,    -1,    37,    56,   105,    55,    -1,
      46,    -1,    47,    -1,   108,    -1,   107,   106,   109,    -1,
      95,   106,   109,    -1,    25,    -1,    37,    56,   107,    55,
      -1,   108,    -1,    95,    -1,   112,    -1,   110,    39,   112,
      -1,   112,    -1,   111,    39,   112,    -1,   112,    46,   112,
      -1,   112,    45,   112,    -1,   112,    47,   112,    -1,   112,
      48,   112,    -1,   112,    49,   112,    -1,   112,    40,   112,
      -1,   112,    23,   112,    -1,   112,    43,   112,    -1,   112,
      44,   112,    -1,   112,    41,   112,    -1,   112,    42,   112,
      -1,   112,    21,   112,    -1,   112,    11,   112,    -1,   112,
      12,   112,    -1,   112,    10,   112,    -1,   112,    30,   112,
      -1,    45,   112,    -1,    50,   112,    -1,    56,   110,    55,
      -1,   113,    -1,   114,    -1,    31,    -1,    32,    -1,    95,
      -1,    95,    29,    -1,    28,    95,    -1,    37,    56,   112,
      55,    -1,    26,    -1,    33,    -1,    25,    -1,   114,    -1,
      95,    -1
};

/* YYRLINE[YYN] -- source line where rule number YYN was defined.  */
static const yytype_uint16 yyrline[] =
{
       0,   139,   139,   140,   141,   144,   147,   148,   151,   152,
     155,   156,   157,   158,   159,   160,   161,   162,   163,   164,
     165,   169,   170,   171,   172,   173,   174,   175,   176,   177,
     178,   179,   180,   181,   182,   186,   187,   190,   194,   195,
     196,   197,   198,   199,   204,   205,   206,   207,   211,   215,
     215,   216,   216,   217,   217,   218,   218,   219,   220,   223,
     224,   225,   226,   227,   228,   229,   230,   231,   232,   233,
     234,   235,   238,   239,   240,   240,   241,   242,   245,   246,
     249,   250,   253,   254,   257,   258,   261,   262,   265,   268,
     271,   272,   275,   278,   281,   282,   285,   288,   289,   292,
     295,   296,   297,   298,   299,   302,   303,   306,   307,   310,
     311,   312,   315,   318,   319,   322,   323,   326,   327,   328,
     329,   330,   331,   334,   335,   338,   339,   342,   343,   346,
     347,   348,   349,   350,   353,   354,   355,   356,   357,   358,
     359,   360,   361,   362,   363,   364,   365,   366,   367,   368,
     369,   370,   373,   374,   377,   378,   379,   382,   383,   386,
     387,   390,   391,   394,   395,   398,   399,   400,   401,   402,
     403,   404,   405,   406,   407,   408,   409,   410,   411,   412,
     413,   421,   422,   423,   424,   427,   428,   429,   430,   431,
     432,   433,   436,   437,   438,   441,   442
};
#endif

#if YYDEBUG || YYERROR_VERBOSE || YYTOKEN_TABLE
/* YYTNAME[SYMBOL-NUM] -- String name of the symbol SYMBOL-NUM.
   First, the terminals, then, starting at YYNTOKENS, nonterminals.  */
static const char *const yytname[] =
{
  "$end", "error", "$undefined", "ALIAS", "FOR", "IF", "ELSE", "DO",
  "WHILE", "INCLUDE", "NE_TOKEN", "LEQ_TOKEN", "GEQ_TOKEN", "CMP_TOKEN",
  "VARIABLE", "EXPR_OP", "INF", "XNAN", "DUMMY", "DOUBLE_DOT", "EMPTY_VAR",
  "EXPO", "CLAUSE_X", "EQ_TOKEN", "SUB_EXPR", "STRING", "INTEGER",
  "CLAUSE_INT", "PREFIX", "POSTFIX", "MATH_SHORT", "KEYWORD", "KEYWORD_B",
  "REAL", "TOKEN", "WORD_B", "COMMAND", "FUNCTION", "LOWER_THAN_ELSE",
  "','", "'='", "'|'", "'&'", "'<'", "'>'", "'-'", "'+'", "'.'", "'*'",
  "'/'", "'!'", "UMINUS", "'{'", "'}'", "';'", "')'", "'('", "']'", "'['",
  "$accept", "control_file", "directive", "statement", "statement_block",
  "ind_statement", "ctypeI", "ctypeII", "assignment", "ctypeIII",
  "ctypeIV", "ctypeV", "ctype_bad", "$@1", "$@2", "$@3", "$@4", "expr_op",
  "format_atom", "$@5", "word", "format_clause", "clause_elem",
  "clause_content", "clause", "begin_clause", "for_statement",
  "if_statement", "do_statement", "while_statement", "opt_expr",
  "include_com", "tok_list", "tok", "array_index_elem", "array_index",
  "var1", "var2", "var3", "loop_clause1", "loop_clause", "array_var_list",
  "var", "var_list1", "poly_commas", "var_list", "imm_expr", "str_op",
  "str_expr", "str_elem", "str_elem1", "gen_expr_list", "gen_expr_list1",
  "gen_expr", "gen_element", "imm_element", "imm_element1", 0
};
#endif

# ifdef YYPRINT
/* YYTOKNUM[YYLEX-NUM] -- Internal token number corresponding to
   token YYLEX-NUM.  */
static const yytype_uint16 yytoknum[] =
{
       0,   256,   257,   258,   259,   260,   261,   262,   263,   264,
     265,   266,   267,   268,   269,   270,   271,   272,   273,   274,
     275,   276,   277,   278,   279,   280,   281,   282,   283,   284,
     285,   286,   287,   288,   289,   290,   291,   292,   293,    44,
      61,   124,    38,    60,    62,    45,    43,    46,    42,    47,
      33,   294,   123,   125,    59,    41,    40,    93,    91
};
# endif

/* YYR1[YYN] -- Symbol number of symbol that rule YYN derives.  */
static const yytype_uint8 yyr1[] =
{
       0,    59,    60,    60,    60,    61,    62,    62,    63,    63,
      64,    64,    64,    64,    64,    64,    64,    64,    64,    64,
      64,    65,    65,    65,    65,    65,    65,    65,    65,    65,
      65,    65,    65,    65,    65,    66,    66,    67,    68,    68,
      68,    68,    68,    68,    69,    69,    69,    69,    70,    72,
      71,    73,    71,    74,    71,    75,    71,    71,    71,    76,
      76,    76,    76,    76,    76,    76,    76,    76,    76,    76,
      76,    76,    77,    77,    78,    77,    77,    77,    79,    79,
      80,    80,    81,    81,    82,    82,    83,    83,    84,    85,
      86,    86,    87,    88,    89,    89,    90,    91,    91,    92,
      93,    93,    93,    93,    93,    94,    94,    95,    95,    96,
      96,    96,    97,    98,    98,    99,    99,   100,   100,   100,
     100,   100,   100,   101,   101,   102,   102,   103,   103,   104,
     104,   104,   104,   104,   105,   105,   105,   105,   105,   105,
     105,   105,   105,   105,   105,   105,   105,   105,   105,   105,
     105,   105,   106,   106,   107,   107,   107,   108,   108,   109,
     109,   110,   110,   111,   111,   112,   112,   112,   112,   112,
     112,   112,   112,   112,   112,   112,   112,   112,   112,   112,
     112,   112,   112,   112,   112,   113,   113,   113,   113,   113,
     113,   113,   114,   114,   114,   115,   115
};

/* YYR2[YYN] -- Number of symbols composing right hand side of rule YYN.  */
static const yytype_uint8 yyr2[] =
{
       0,     2,     1,     1,     2,     1,     1,     3,     1,     2,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     2,     4,     2,     3,     5,     3,     2,     3,     4,
       5,     1,     3,     1,     3,     3,     4,     3,     4,     1,
       3,     2,     2,     3,     4,     6,     4,     5,     3,     0,
       4,     0,     6,     0,     5,     0,     7,     2,     3,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     0,     3,     4,     3,     1,     1,
       1,     3,     1,     1,     1,     3,     3,     2,     1,     9,
       5,     7,     6,     5,     0,     1,     3,     1,     3,     1,
       1,     3,     2,     2,     1,     1,     3,     1,     3,     1,
       3,     0,     4,     8,     4,     7,     3,     1,     3,     3,
       1,     5,     3,     1,     2,     1,     1,     2,     2,     1,
       2,     2,     3,     3,     1,     3,     3,     3,     3,     3,
       3,     3,     3,     3,     3,     3,     3,     3,     3,     2,
       3,     4,     1,     1,     1,     3,     3,     1,     4,     1,
       1,     1,     3,     1,     3,     3,     3,     3,     3,     3,
       3,     3,     3,     3,     3,     3,     3,     3,     3,     3,
       3,     2,     2,     3,     1,     1,     1,     1,     1,     2,
       2,     4,     1,     1,     1,     1,     1
};

/* YYDEFACT[STATE-NAME] -- Default rule to reduce with in state
   STATE-NUM when YYTABLE doesn't specify something else to do.  Zero
   means the default is an error.  */
static const yytype_uint8 yydefact[] =
{
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,    79,     0,    20,     0,     0,     3,     2,     6,    10,
      11,    39,    12,    13,    14,    15,   107,    16,    17,    18,
      19,     5,     0,    33,    31,    94,     0,     0,     0,   157,
      78,     0,    79,     0,   107,     0,     0,   154,    42,     0,
      79,     0,    49,   111,   123,   126,   125,   129,     0,     0,
     194,   192,     0,   186,   187,   193,     0,     0,     0,     0,
     188,     0,   184,   185,    57,     0,     0,   104,   105,     0,
     100,    27,     0,     0,     0,     0,   111,    88,     0,     0,
     123,    23,    21,   134,     0,     0,     8,   111,     1,     4,
       0,    68,    66,    67,    71,    70,    69,    64,    65,    63,
      59,    62,    60,    61,     0,    41,     0,     0,     0,     0,
       0,    95,     0,     0,     0,     0,   152,   153,     0,    96,
       0,     0,     0,   127,   130,     0,   120,     0,     0,   109,
     117,     0,   124,   128,   131,     0,     0,     0,   190,     0,
     181,   182,     0,   161,   189,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,   102,     0,   108,   103,     0,    28,    53,
      48,     0,     0,   149,     0,     0,   107,   196,     0,   123,
      26,    24,    73,    72,     0,    87,    80,    83,    84,     0,
      82,   163,     0,    35,    97,    99,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     7,     9,     0,    40,    58,    43,    37,    34,    32,
      94,     0,     0,     0,     0,   160,   159,   156,   155,    46,
      36,    50,     0,     0,   111,   116,     0,   132,   133,     0,
       0,   183,   179,   177,   178,   176,   171,   180,   170,   174,
     175,   172,   173,   166,   165,   167,   168,   169,     0,    51,
     112,   106,   101,    29,     0,     0,    44,     0,     0,   150,
       0,    74,     0,     0,     0,     0,    86,     0,    38,     0,
      22,   145,   146,   147,   148,   144,   143,   142,   140,   141,
     137,   135,   136,   139,   138,   111,   114,     0,    90,     0,
      93,   158,   119,   110,   122,     0,   109,   118,     0,   191,
     162,    47,     0,    30,    54,    55,   151,   108,    25,    75,
       0,     0,    77,    81,    85,   164,    98,   109,    94,     0,
      92,     0,     0,    45,    52,     0,    76,     0,     0,    91,
     121,     0,    56,     0,     0,   115,   113,    89
};

/* YYDEFGOTO[NTERM-NUM].  */
static const yytype_int16 yydefgoto[] =
{
      -1,    15,    16,    17,    95,    18,    19,    20,    21,    22,
      23,    24,    25,   135,   322,   275,   345,   114,   196,   329,
      44,   197,   198,   199,   142,    89,    27,    28,    29,    30,
     120,    31,   203,   204,    78,    79,    70,   140,    33,    34,
      55,   141,    56,    57,    58,    59,   188,   130,    46,    47,
     237,   152,   200,    80,    72,    73,    94
};

/* YYPACT[STATE-NUM] -- Index in YYTABLE of the portion describing
   STATE-NUM.  */
#define YYPACT_NINF -309
static const yytype_int16 yypact[] =
{
     376,   -26,   -18,   389,    46,   276,   155,   130,  1277,   516,
    1171,   960,    36,  -309,    59,   103,  -309,  -309,  -309,  -309,
    -309,    44,  -309,  -309,  -309,  -309,   208,  -309,  -309,  -309,
    -309,  -309,   232,    67,    70,  1277,  1277,   104,  1277,  -309,
    -309,  1198,  -309,    88,  -309,    10,    27,  -309,  -309,  1277,
     155,   166,  -309,   314,    83,  -309,  -309,  -309,   413,    -4,
    -309,  -309,   155,  -309,  -309,  -309,    89,  1277,  1277,  1277,
     119,   568,  -309,  -309,  -309,  1198,  1277,  -309,  -309,    -7,
     990,   -16,   302,  1277,    95,  1277,   253,  -309,   497,  1224,
     467,    84,   115,  1044,  1084,   431,  -309,   314,  -309,  -309,
    1277,  -309,  -309,  -309,  -309,  -309,  -309,  -309,  -309,  -309,
    -309,  -309,  -309,  -309,  1277,  -309,  1277,  1277,   419,   419,
     102,  1030,   608,   105,   648,   276,  -309,  -309,   276,  -309,
     276,   688,   214,  -309,  -309,  1277,  -309,  1277,  1198,  -309,
    -309,    -6,  -309,  -309,  -309,   107,   166,   413,  -309,  1277,
     137,   137,     8,  1030,  -309,  1277,  1277,  1277,  1277,  1277,
    1277,  1277,  1277,  1277,  1277,  1277,  1277,  1277,  1277,  1277,
    1277,     6,    20,  1030,  1198,  -309,  1277,   419,   129,  -309,
      84,   728,  1144,  1030,  1198,  1144,    21,  -309,   116,   528,
      84,   133,  -309,    -1,  1251,  -309,  -309,   142,  -309,    12,
     145,  1030,  1277,   148,  -309,  -309,   419,  1277,  1277,  1277,
    1277,  1277,  1277,  1277,  1277,  1277,  1277,  1277,  1277,  1277,
    1277,  -309,  -309,    45,   153,  1030,  1030,  1030,    84,    84,
    1277,   389,  1277,   389,   128,  -309,  -309,  -309,  -309,  -309,
     148,  1030,   768,    65,   367,  -309,  1277,  -309,  -309,   808,
    1277,  -309,  1099,   167,   167,   137,  1099,  1030,  1030,  1092,
    1139,   167,   167,     4,     4,     4,   137,   137,   419,  -309,
    -309,  -309,  1030,    84,   419,  1277,   136,   149,    74,  -309,
     419,  -309,     9,    77,     9,  1251,  -309,  1277,  1030,   214,
      84,  1030,  1030,  1030,  1030,  1030,  1030,  1030,  1030,  1030,
    1030,  1030,  1030,  1030,  1030,   367,  -309,   163,   197,   848,
    -309,  -309,  -309,  -309,  -309,  1277,   181,  -309,   888,  -309,
    1030,    84,  1277,    84,  1030,  -309,  -309,    85,    84,  -309,
       9,    87,  -309,  -309,  -309,  1030,  -309,   183,  1277,   389,
    -309,   928,  1277,  -309,  1030,  1277,  -309,  1277,   152,  -309,
    -309,    94,  1030,    97,   389,  -309,  -309,  -309
};

/* YYPGOTO[NTERM-NUM].  */
static const yytype_int16 yypgoto[] =
{
    -309,  -309,  -309,     5,  -309,    -2,  -309,  -309,  -309,  -309,
    -309,  -309,  -309,  -309,  -309,  -309,  -309,  -309,   -60,  -309,
       1,  -177,   -55,  -309,    13,  -309,  -309,  -309,  -309,  -309,
    -221,  -309,    93,   -52,    72,   -57,     0,  -216,  -309,  -309,
    -309,   144,  -309,   -32,   -30,     3,    -9,   213,   138,   -50,
     134,   165,  -308,   159,  -309,    11,  -309
};

/* YYTABLE[YYPACT[STATE-NUM]].  What to do in state STATE-NUM.  If
   positive, shift that token.  If negative, reduce the rule which
   number is the opposite.  If zero, do what YYDEFACT says.
   If YYTABLE_NINF, syntax error.  */
#define YYTABLE_NINF -197
static const yytype_int16 yytable[] =
{
      32,    26,    92,    32,    26,    45,    48,    54,    37,   307,
      96,    90,    32,    26,    91,    32,    26,   283,   172,   134,
      99,   281,    93,   177,    88,   158,   144,   145,   317,   147,
      35,   192,   174,   244,   351,   146,   193,   268,    36,   353,
       1,     2,    87,     3,     4,   129,   269,   250,   175,   245,
     132,    54,   169,   170,   139,   282,   126,   127,    54,   174,
    -109,   147,   148,   251,     6,   330,   285,     7,     8,   286,
       9,    10,    11,   126,   127,   270,  -109,    97,   236,   191,
     236,   243,    54,   100,   305,   180,   187,   186,   189,   317,
      13,   190,    14,   222,   178,    32,    26,    93,   139,    93,
     306,   205,    38,    98,   174,   331,   118,     1,     2,   119,
       3,     4,   123,   174,   247,   248,   284,   348,    54,    54,
     313,   228,   229,   146,  -110,    45,   284,   278,   235,   327,
     235,     6,   332,   287,     7,     8,   287,     9,    10,    11,
    -110,    87,   346,   205,   125,   149,    54,    54,   154,   355,
     147,   182,   356,   283,   206,    12,   230,    13,   158,    14,
     147,   232,    49,   246,    40,    41,    50,    71,   274,    51,
      52,   279,   280,   277,   126,   127,   325,    54,  -197,  -197,
     273,   284,   187,   311,   287,   187,    53,   289,   158,    40,
      41,    42,   250,    93,   121,   122,    93,   124,   147,   147,
      40,    41,    42,   339,   326,   133,    54,   354,   131,   290,
    -197,  -197,   166,   167,   168,   169,   170,   338,   101,   102,
     103,   342,    53,   347,   333,   240,   150,   151,   153,   104,
     334,    32,    26,    32,    26,   173,   308,   336,   310,    60,
      61,   223,   181,   147,   183,   316,   271,    65,   201,   105,
     106,   107,   108,   109,   110,   111,   112,   113,   128,   153,
     147,   115,   116,   234,   238,   224,     0,     0,    54,     0,
       0,   321,   117,   225,    54,   226,   227,   323,    60,    61,
      54,     0,     0,   328,   136,   137,    65,    40,   184,    42,
      84,   147,     0,   147,   241,     0,   242,     0,   147,     0,
     205,    39,     0,    85,     0,     0,   337,     0,   249,   185,
      40,    41,    42,    43,   252,   253,   254,   255,   256,   257,
     258,   259,   260,   261,   262,   263,   264,   265,   266,   267,
       0,     0,     0,     0,     0,   272,    40,    41,    42,    32,
      26,    51,   179,     0,   349,   136,   137,     0,    40,   138,
      42,     0,     0,   153,    32,    26,     0,     0,    53,   357,
       0,   288,     0,     0,     0,     0,   291,   292,   293,   294,
     295,   296,   297,   298,   299,   300,   301,   302,   303,   304,
       1,     2,     0,     3,     4,     5,     0,     0,     0,   121,
       0,   309,     0,     1,     2,     0,     3,     4,   314,   315,
       0,    40,   138,    42,     6,   318,     0,     7,     8,   320,
       9,    10,    11,     0,     0,     0,     0,     6,     0,     0,
       7,     8,     0,     9,    10,    11,     0,     0,    12,     0,
      13,     0,    14,     0,   324,     1,     2,     0,     3,     4,
       0,    12,     0,    13,   201,    14,   335,    40,    41,    42,
       0,     0,   143,    40,    41,    42,     0,     0,    51,     6,
       0,     0,     7,     8,     0,     9,    10,    11,     0,    53,
       0,     0,     0,     0,   341,    53,     0,  -196,  -196,  -196,
       0,   344,     0,     0,   221,    13,     0,    14,  -196,     0,
    -196,     0,    60,    61,     0,     0,     0,   121,     0,     0,
      65,   201,     0,     0,   352,     0,   201,   202,  -196,  -196,
    -196,  -196,  -196,  -196,  -196,  -196,  -196,    74,     0,     0,
       0,     0,    60,    61,     0,    87,   -78,   -78,   -78,     0,
      65,    40,    41,    42,    84,     0,    51,   -78,  -196,  -196,
    -196,     0,     0,     0,     0,   -78,   -78,    85,     0,  -196,
       0,  -196,     0,    86,     0,     0,   -78,   -78,   -78,   -78,
     -78,   -78,   -78,   -78,   -78,   -78,     0,     0,     0,  -196,
    -196,  -196,  -196,  -196,  -196,  -196,  -196,  -196,   155,   156,
     157,     0,     0,     0,     0,     0,    87,     0,     0,   158,
       0,   159,     0,     0,     0,     0,     0,     0,   160,     0,
       0,     0,     0,     0,     0,     0,     0,     0,   161,   162,
     163,   164,   165,   166,   167,   168,   169,   170,   155,   156,
     157,     0,     0,   171,     0,     0,     0,     0,     0,   158,
       0,   159,     0,     0,     0,     0,     0,     0,   160,     0,
       0,     0,     0,     0,     0,     0,     0,     0,   161,   162,
     163,   164,   165,   166,   167,   168,   169,   170,   155,   156,
     157,     0,     0,   231,     0,     0,     0,     0,     0,   158,
       0,   159,     0,     0,     0,     0,     0,     0,   160,     0,
       0,     0,     0,     0,     0,     0,     0,     0,   161,   162,
     163,   164,   165,   166,   167,   168,   169,   170,   155,   156,
     157,     0,     0,   233,     0,     0,     0,     0,     0,   158,
       0,   159,     0,     0,     0,     0,     0,     0,   160,     0,
       0,     0,     0,     0,     0,     0,     0,     0,   161,   162,
     163,   164,   165,   166,   167,   168,   169,   170,   155,   156,
     157,     0,     0,   239,     0,     0,     0,     0,     0,   158,
       0,   159,     0,     0,     0,     0,     0,     0,   160,     0,
       0,     0,     0,     0,     0,     0,     0,     0,   161,   162,
     163,   164,   165,   166,   167,   168,   169,   170,   155,   156,
     157,     0,     0,   276,     0,     0,     0,     0,     0,   158,
       0,   159,     0,     0,     0,     0,     0,     0,   160,     0,
       0,     0,     0,     0,     0,     0,     0,     0,   161,   162,
     163,   164,   165,   166,   167,   168,   169,   170,   155,   156,
     157,     0,     0,   312,     0,     0,     0,     0,     0,   158,
       0,   159,     0,     0,     0,     0,     0,     0,   160,     0,
       0,     0,     0,     0,     0,     0,     0,     0,   161,   162,
     163,   164,   165,   166,   167,   168,   169,   170,   155,   156,
     157,     0,     0,   319,     0,     0,     0,     0,     0,   158,
       0,   159,     0,     0,     0,     0,     0,     0,   160,     0,
       0,     0,     0,     0,     0,     0,     0,     0,   161,   162,
     163,   164,   165,   166,   167,   168,   169,   170,   155,   156,
     157,     0,     0,   340,     0,     0,     0,     0,     0,   158,
       0,   159,     0,     0,     0,     0,     0,     0,   160,     0,
       0,     0,     0,     0,     0,     0,     0,     0,   161,   162,
     163,   164,   165,   166,   167,   168,   169,   170,   155,   156,
     157,     0,     0,   343,     0,     0,     0,     0,     0,   158,
       0,   159,     0,     0,     0,     0,     0,     0,   160,     0,
       0,     0,     0,     0,     0,     0,     0,     0,   161,   162,
     163,   164,   165,   166,   167,   168,   169,   170,    81,     0,
       0,     0,     0,   350,     0,    60,    61,     0,     0,     0,
       0,    82,    83,    65,    40,    41,    42,    84,     0,    51,
     155,   156,   157,     0,     0,     0,     0,     0,     0,   176,
      85,   158,     0,   159,     0,     0,    86,     0,    87,     0,
     160,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     161,   162,   163,   164,   165,   166,   167,   168,   169,   170,
     155,   156,   157,     0,     0,     0,     0,     0,     0,     0,
       0,   158,     0,   159,  -195,  -195,  -195,     0,     0,     0,
     160,     0,     0,     0,     0,  -195,     0,  -195,     0,     0,
     161,   162,   163,   164,   165,   166,   167,   168,   169,   170,
       0,     0,     0,     0,     0,  -195,  -195,  -195,  -195,  -195,
    -195,  -195,  -195,  -195,   207,   208,   209,     0,     0,     0,
       0,     0,   155,   156,   157,   210,     0,   211,     0,  -197,
     156,   157,     0,   158,     0,   159,     0,     0,     0,     0,
     158,     0,  -197,     0,     0,   212,   213,   214,   215,   216,
     217,   218,   219,   220,   163,   164,   165,   166,   167,   168,
     169,   170,   164,   165,   166,   167,   168,   169,   170,   155,
     156,   157,     0,     0,     0,     0,     0,     0,     0,     0,
     158,     0,   159,     0,     0,     0,     0,     0,     0,    60,
      61,     0,     0,     0,     0,     0,     0,    65,    40,    41,
      42,    84,   164,   165,   166,   167,   168,   169,   170,    75,
      76,     0,     0,     0,    85,     0,    60,    61,    77,    62,
     185,     0,    63,    64,    65,    40,    41,    42,    66,     0,
       0,     0,     0,     0,     0,     0,    67,    76,     0,     0,
       0,    68,     0,    60,    61,    77,    62,    69,     0,    63,
      64,    65,    40,    41,    42,    66,     0,     0,     0,     0,
       0,     0,     0,    67,     0,     0,   192,     0,    68,    60,
      61,   193,    62,     0,    69,    63,    64,    65,    40,    41,
      42,    66,     0,     0,     0,     0,     0,     0,     0,    67,
       0,     0,     0,   192,    68,     0,    60,    61,   193,    62,
     194,   195,    63,    64,    65,    40,    41,    42,    66,     0,
       0,     0,     0,     0,     0,     0,    67,     0,     0,     0,
       0,    68,    60,    61,     0,    62,     0,   194,    63,    64,
      65,    40,    41,    42,    66,     0,     0,     0,     0,     0,
       0,     0,    67,     0,     0,     0,     0,    68,     0,     0,
       0,     0,     0,    69
};

static const yytype_int16 yycheck[] =
{
       0,     0,    11,     3,     3,     5,     6,     7,     3,   230,
      12,    11,    12,    12,    11,    15,    15,   194,    75,    51,
      15,    22,    11,    39,    11,    21,    58,    31,   244,    59,
      56,    22,    39,    39,   342,    39,    27,    31,    56,   347,
       4,     5,    58,     7,     8,    18,    40,    39,    55,    55,
      50,    51,    48,    49,    53,    56,    46,    47,    58,    39,
      39,    91,    62,    55,    28,    56,    54,    31,    32,    57,
      34,    35,    36,    46,    47,    55,    55,    18,   128,    88,
     130,   138,    82,    39,    39,    82,    86,    86,    88,   305,
      54,    88,    56,    95,    81,    95,    95,    86,    97,    88,
      55,    90,    56,     0,    39,   282,    39,     4,     5,    39,
       7,     8,     8,    39,   146,   147,    39,   338,   118,   119,
      55,   118,   119,    39,    39,   125,    39,   184,   128,    55,
     130,    28,    55,    39,    31,    32,    39,    34,    35,    36,
      55,    58,    55,   132,    56,    56,   146,   147,    29,    55,
     180,    56,    55,   330,    39,    52,    54,    54,    21,    56,
     190,    56,    32,    56,    34,    35,    36,     8,    39,    39,
      40,    55,    39,   182,    46,    47,    40,   177,    11,    12,
     177,    39,   182,    55,    39,   185,    56,    39,    21,    34,
      35,    36,    39,   182,    35,    36,   185,    38,   228,   229,
      34,    35,    36,     6,    55,    39,   206,    55,    49,   206,
      43,    44,    45,    46,    47,    48,    49,    54,    10,    11,
      12,    40,    56,    40,   284,   132,    67,    68,    69,    21,
     285,   231,   231,   233,   233,    76,   231,   289,   233,    25,
      26,    97,    83,   273,    85,   244,   174,    33,    89,    41,
      42,    43,    44,    45,    46,    47,    48,    49,    45,   100,
     290,    29,    30,   125,   130,   100,    -1,    -1,   268,    -1,
      -1,   268,    40,   114,   274,   116,   117,   274,    25,    26,
     280,    -1,    -1,   280,    31,    32,    33,    34,    35,    36,
      37,   321,    -1,   323,   135,    -1,   137,    -1,   328,    -1,
     289,    25,    -1,    50,    -1,    -1,   305,    -1,   149,    56,
      34,    35,    36,    37,   155,   156,   157,   158,   159,   160,
     161,   162,   163,   164,   165,   166,   167,   168,   169,   170,
      -1,    -1,    -1,    -1,    -1,   176,    34,    35,    36,   339,
     339,    39,    40,    -1,   339,    31,    32,    -1,    34,    35,
      36,    -1,    -1,   194,   354,   354,    -1,    -1,    56,   354,
      -1,   202,    -1,    -1,    -1,    -1,   207,   208,   209,   210,
     211,   212,   213,   214,   215,   216,   217,   218,   219,   220,
       4,     5,    -1,     7,     8,     9,    -1,    -1,    -1,   230,
      -1,   232,    -1,     4,     5,    -1,     7,     8,    31,    32,
      -1,    34,    35,    36,    28,   246,    -1,    31,    32,   250,
      34,    35,    36,    -1,    -1,    -1,    -1,    28,    -1,    -1,
      31,    32,    -1,    34,    35,    36,    -1,    -1,    52,    -1,
      54,    -1,    56,    -1,   275,     4,     5,    -1,     7,     8,
      -1,    52,    -1,    54,   285,    56,   287,    34,    35,    36,
      -1,    -1,    39,    34,    35,    36,    -1,    -1,    39,    28,
      -1,    -1,    31,    32,    -1,    34,    35,    36,    -1,    56,
      -1,    -1,    -1,    -1,   315,    56,    -1,    10,    11,    12,
      -1,   322,    -1,    -1,    53,    54,    -1,    56,    21,    -1,
      23,    -1,    25,    26,    -1,    -1,    -1,   338,    -1,    -1,
      33,   342,    -1,    -1,   345,    -1,   347,    40,    41,    42,
      43,    44,    45,    46,    47,    48,    49,     1,    -1,    -1,
      -1,    -1,    25,    26,    -1,    58,    10,    11,    12,    -1,
      33,    34,    35,    36,    37,    -1,    39,    21,    10,    11,
      12,    -1,    -1,    -1,    -1,    29,    30,    50,    -1,    21,
      -1,    23,    -1,    56,    -1,    -1,    40,    41,    42,    43,
      44,    45,    46,    47,    48,    49,    -1,    -1,    -1,    41,
      42,    43,    44,    45,    46,    47,    48,    49,    10,    11,
      12,    -1,    -1,    -1,    -1,    -1,    58,    -1,    -1,    21,
      -1,    23,    -1,    -1,    -1,    -1,    -1,    -1,    30,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    40,    41,
      42,    43,    44,    45,    46,    47,    48,    49,    10,    11,
      12,    -1,    -1,    55,    -1,    -1,    -1,    -1,    -1,    21,
      -1,    23,    -1,    -1,    -1,    -1,    -1,    -1,    30,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    40,    41,
      42,    43,    44,    45,    46,    47,    48,    49,    10,    11,
      12,    -1,    -1,    55,    -1,    -1,    -1,    -1,    -1,    21,
      -1,    23,    -1,    -1,    -1,    -1,    -1,    -1,    30,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    40,    41,
      42,    43,    44,    45,    46,    47,    48,    49,    10,    11,
      12,    -1,    -1,    55,    -1,    -1,    -1,    -1,    -1,    21,
      -1,    23,    -1,    -1,    -1,    -1,    -1,    -1,    30,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    40,    41,
      42,    43,    44,    45,    46,    47,    48,    49,    10,    11,
      12,    -1,    -1,    55,    -1,    -1,    -1,    -1,    -1,    21,
      -1,    23,    -1,    -1,    -1,    -1,    -1,    -1,    30,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    40,    41,
      42,    43,    44,    45,    46,    47,    48,    49,    10,    11,
      12,    -1,    -1,    55,    -1,    -1,    -1,    -1,    -1,    21,
      -1,    23,    -1,    -1,    -1,    -1,    -1,    -1,    30,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    40,    41,
      42,    43,    44,    45,    46,    47,    48,    49,    10,    11,
      12,    -1,    -1,    55,    -1,    -1,    -1,    -1,    -1,    21,
      -1,    23,    -1,    -1,    -1,    -1,    -1,    -1,    30,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    40,    41,
      42,    43,    44,    45,    46,    47,    48,    49,    10,    11,
      12,    -1,    -1,    55,    -1,    -1,    -1,    -1,    -1,    21,
      -1,    23,    -1,    -1,    -1,    -1,    -1,    -1,    30,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    40,    41,
      42,    43,    44,    45,    46,    47,    48,    49,    10,    11,
      12,    -1,    -1,    55,    -1,    -1,    -1,    -1,    -1,    21,
      -1,    23,    -1,    -1,    -1,    -1,    -1,    -1,    30,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    40,    41,
      42,    43,    44,    45,    46,    47,    48,    49,    10,    11,
      12,    -1,    -1,    55,    -1,    -1,    -1,    -1,    -1,    21,
      -1,    23,    -1,    -1,    -1,    -1,    -1,    -1,    30,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    40,    41,
      42,    43,    44,    45,    46,    47,    48,    49,    18,    -1,
      -1,    -1,    -1,    55,    -1,    25,    26,    -1,    -1,    -1,
      -1,    31,    32,    33,    34,    35,    36,    37,    -1,    39,
      10,    11,    12,    -1,    -1,    -1,    -1,    -1,    -1,    19,
      50,    21,    -1,    23,    -1,    -1,    56,    -1,    58,    -1,
      30,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      40,    41,    42,    43,    44,    45,    46,    47,    48,    49,
      10,    11,    12,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    21,    -1,    23,    10,    11,    12,    -1,    -1,    -1,
      30,    -1,    -1,    -1,    -1,    21,    -1,    23,    -1,    -1,
      40,    41,    42,    43,    44,    45,    46,    47,    48,    49,
      -1,    -1,    -1,    -1,    -1,    41,    42,    43,    44,    45,
      46,    47,    48,    49,    10,    11,    12,    -1,    -1,    -1,
      -1,    -1,    10,    11,    12,    21,    -1,    23,    -1,    10,
      11,    12,    -1,    21,    -1,    23,    -1,    -1,    -1,    -1,
      21,    -1,    23,    -1,    -1,    41,    42,    43,    44,    45,
      46,    47,    48,    49,    42,    43,    44,    45,    46,    47,
      48,    49,    43,    44,    45,    46,    47,    48,    49,    10,
      11,    12,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      21,    -1,    23,    -1,    -1,    -1,    -1,    -1,    -1,    25,
      26,    -1,    -1,    -1,    -1,    -1,    -1,    33,    34,    35,
      36,    37,    43,    44,    45,    46,    47,    48,    49,    18,
      19,    -1,    -1,    -1,    50,    -1,    25,    26,    27,    28,
      56,    -1,    31,    32,    33,    34,    35,    36,    37,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    45,    19,    -1,    -1,
      -1,    50,    -1,    25,    26,    27,    28,    56,    -1,    31,
      32,    33,    34,    35,    36,    37,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    45,    -1,    -1,    22,    -1,    50,    25,
      26,    27,    28,    -1,    56,    31,    32,    33,    34,    35,
      36,    37,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    45,
      -1,    -1,    -1,    22,    50,    -1,    25,    26,    27,    28,
      56,    57,    31,    32,    33,    34,    35,    36,    37,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    45,    -1,    -1,    -1,
      -1,    50,    25,    26,    -1,    28,    -1,    56,    31,    32,
      33,    34,    35,    36,    37,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    45,    -1,    -1,    -1,    -1,    50,    -1,    -1,
      -1,    -1,    -1,    56
};

/* YYSTOS[STATE-NUM] -- The (internal number of the) accessing
   symbol of state STATE-NUM.  */
static const yytype_uint8 yystos[] =
{
       0,     4,     5,     7,     8,     9,    28,    31,    32,    34,
      35,    36,    52,    54,    56,    60,    61,    62,    64,    65,
      66,    67,    68,    69,    70,    71,    79,    85,    86,    87,
      88,    90,    95,    97,    98,    56,    56,    62,    56,    25,
      34,    35,    36,    37,    79,    95,   107,   108,    95,    32,
      36,    39,    40,    56,    95,    99,   101,   102,   103,   104,
      25,    26,    28,    31,    32,    33,    37,    45,    50,    56,
      95,   112,   113,   114,     1,    18,    19,    27,    93,    94,
     112,    18,    31,    32,    37,    50,    56,    58,    83,    84,
      95,   104,   105,   114,   115,    63,    64,    18,     0,    62,
      39,    10,    11,    12,    21,    41,    42,    43,    44,    45,
      46,    47,    48,    49,    76,    29,    30,    40,    39,    39,
      89,   112,   112,     8,   112,    56,    46,    47,   106,    18,
     106,   112,    95,    39,   102,    72,    31,    32,    35,    79,
      96,   100,    83,    39,   102,    31,    39,   103,    95,    56,
     112,   112,   110,   112,    29,    10,    11,    12,    21,    23,
      30,    40,    41,    42,    43,    44,    45,    46,    47,    48,
      49,    55,    94,   112,    39,    55,    19,    39,    83,    40,
     104,   112,    56,   112,    35,    56,    79,    95,   105,    95,
     104,   105,    22,    27,    56,    57,    77,    80,    81,    82,
     111,   112,    40,    91,    92,   114,    39,    10,    11,    12,
      21,    23,    41,    42,    43,    44,    45,    46,    47,    48,
      49,    53,    64,   100,   110,   112,   112,   112,   104,   104,
      54,    55,    56,    55,   107,    95,   108,   109,   109,    55,
      91,   112,   112,    94,    39,    55,    56,   102,   102,   112,
      39,    55,   112,   112,   112,   112,   112,   112,   112,   112,
     112,   112,   112,   112,   112,   112,   112,   112,    31,    40,
      55,    93,   112,   104,    39,    74,    55,   105,    94,    55,
      39,    22,    56,    80,    39,    54,    57,    39,   112,    39,
     104,   112,   112,   112,   112,   112,   112,   112,   112,   112,
     112,   112,   112,   112,   112,    39,    55,    89,    62,   112,
      62,    55,    55,    55,    31,    32,    79,    96,   112,    55,
     112,   104,    73,   104,   112,    40,    55,    55,   104,    78,
      56,    80,    55,    77,    81,   112,    92,    79,    54,     6,
      55,   112,    40,    55,   112,    75,    55,    40,    89,    62,
      55,   111,   112,   111,    55,    55,    55,    62
};

#define yyerrok		(yyerrstatus = 0)
#define yyclearin	(yychar = YYEMPTY)
#define YYEMPTY		(-2)
#define YYEOF		0

#define YYACCEPT	goto yyacceptlab
#define YYABORT		goto yyabortlab
#define YYERROR		goto yyerrorlab


/* Like YYERROR except do call yyerror.  This remains here temporarily
   to ease the transition to the new meaning of YYERROR, for GCC.
   Once GCC version 2 has supplanted version 1, this can go.  */

#define YYFAIL		goto yyerrlab

#define YYRECOVERING()  (!!yyerrstatus)

#define YYBACKUP(Token, Value)					\
do								\
  if (yychar == YYEMPTY && yylen == 1)				\
    {								\
      yychar = (Token);						\
      yylval = (Value);						\
      yytoken = YYTRANSLATE (yychar);				\
      YYPOPSTACK (1);						\
      goto yybackup;						\
    }								\
  else								\
    {								\
      yyerror (YY_("syntax error: cannot back up")); \
      YYERROR;							\
    }								\
while (YYID (0))


#define YYTERROR	1
#define YYERRCODE	256


/* YYLLOC_DEFAULT -- Set CURRENT to span from RHS[1] to RHS[N].
   If N is 0, then set CURRENT to the empty location which ends
   the previous symbol: RHS[0] (always defined).  */

#define YYRHSLOC(Rhs, K) ((Rhs)[K])
#ifndef YYLLOC_DEFAULT
# define YYLLOC_DEFAULT(Current, Rhs, N)				\
    do									\
      if (YYID (N))                                                    \
	{								\
	  (Current).first_line   = YYRHSLOC (Rhs, 1).first_line;	\
	  (Current).first_column = YYRHSLOC (Rhs, 1).first_column;	\
	  (Current).last_line    = YYRHSLOC (Rhs, N).last_line;		\
	  (Current).last_column  = YYRHSLOC (Rhs, N).last_column;	\
	}								\
      else								\
	{								\
	  (Current).first_line   = (Current).last_line   =		\
	    YYRHSLOC (Rhs, 0).last_line;				\
	  (Current).first_column = (Current).last_column =		\
	    YYRHSLOC (Rhs, 0).last_column;				\
	}								\
    while (YYID (0))
#endif


/* YY_LOCATION_PRINT -- Print the location on the stream.
   This macro was not mandated originally: define only if we know
   we won't break user code: when these are the locations we know.  */

#ifndef YY_LOCATION_PRINT
# if YYLTYPE_IS_TRIVIAL
#  define YY_LOCATION_PRINT(File, Loc)			\
     fprintf (File, "%d.%d-%d.%d",			\
	      (Loc).first_line, (Loc).first_column,	\
	      (Loc).last_line,  (Loc).last_column)
# else
#  define YY_LOCATION_PRINT(File, Loc) ((void) 0)
# endif
#endif


/* YYLEX -- calling `yylex' with the right arguments.  */

#ifdef YYLEX_PARAM
# define YYLEX yylex (&yylval, YYLEX_PARAM)
#else
# define YYLEX yylex (&yylval)
#endif

/* Enable debugging if requested.  */
#if YYDEBUG

# ifndef YYFPRINTF
#  include <stdio.h> /* INFRINGES ON USER NAME SPACE */
#  define YYFPRINTF fprintf
# endif

# define YYDPRINTF(Args)			\
do {						\
  if (yydebug)					\
    YYFPRINTF Args;				\
} while (YYID (0))

# define YY_SYMBOL_PRINT(Title, Type, Value, Location)			  \
do {									  \
  if (yydebug)								  \
    {									  \
      YYFPRINTF (stderr, "%s ", Title);					  \
      yy_symbol_print (stderr,						  \
		  Type, Value); \
      YYFPRINTF (stderr, "\n");						  \
    }									  \
} while (YYID (0))


/*--------------------------------.
| Print this symbol on YYOUTPUT.  |
`--------------------------------*/

/*ARGSUSED*/
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yy_symbol_value_print (FILE *yyoutput, int yytype, YYSTYPE const * const yyvaluep)
#else
static void
yy_symbol_value_print (yyoutput, yytype, yyvaluep)
    FILE *yyoutput;
    int yytype;
    YYSTYPE const * const yyvaluep;
#endif
{
  if (!yyvaluep)
    return;
# ifdef YYPRINT
  if (yytype < YYNTOKENS)
    YYPRINT (yyoutput, yytoknum[yytype], *yyvaluep);
# else
  YYUSE (yyoutput);
# endif
  switch (yytype)
    {
      default:
	break;
    }
}


/*--------------------------------.
| Print this symbol on YYOUTPUT.  |
`--------------------------------*/

#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yy_symbol_print (FILE *yyoutput, int yytype, YYSTYPE const * const yyvaluep)
#else
static void
yy_symbol_print (yyoutput, yytype, yyvaluep)
    FILE *yyoutput;
    int yytype;
    YYSTYPE const * const yyvaluep;
#endif
{
  if (yytype < YYNTOKENS)
    YYFPRINTF (yyoutput, "token %s (", yytname[yytype]);
  else
    YYFPRINTF (yyoutput, "nterm %s (", yytname[yytype]);

  yy_symbol_value_print (yyoutput, yytype, yyvaluep);
  YYFPRINTF (yyoutput, ")");
}

/*------------------------------------------------------------------.
| yy_stack_print -- Print the state stack from its BOTTOM up to its |
| TOP (included).                                                   |
`------------------------------------------------------------------*/

#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yy_stack_print (yytype_int16 *yybottom, yytype_int16 *yytop)
#else
static void
yy_stack_print (yybottom, yytop)
    yytype_int16 *yybottom;
    yytype_int16 *yytop;
#endif
{
  YYFPRINTF (stderr, "Stack now");
  for (; yybottom <= yytop; yybottom++)
    {
      int yybot = *yybottom;
      YYFPRINTF (stderr, " %d", yybot);
    }
  YYFPRINTF (stderr, "\n");
}

# define YY_STACK_PRINT(Bottom, Top)				\
do {								\
  if (yydebug)							\
    yy_stack_print ((Bottom), (Top));				\
} while (YYID (0))


/*------------------------------------------------.
| Report that the YYRULE is going to be reduced.  |
`------------------------------------------------*/

#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yy_reduce_print (YYSTYPE *yyvsp, int yyrule)
#else
static void
yy_reduce_print (yyvsp, yyrule)
    YYSTYPE *yyvsp;
    int yyrule;
#endif
{
  int yynrhs = yyr2[yyrule];
  int yyi;
  unsigned long int yylno = yyrline[yyrule];
  YYFPRINTF (stderr, "Reducing stack by rule %d (line %lu):\n",
	     yyrule - 1, yylno);
  /* The symbols being reduced.  */
  for (yyi = 0; yyi < yynrhs; yyi++)
    {
      YYFPRINTF (stderr, "   $%d = ", yyi + 1);
      yy_symbol_print (stderr, yyrhs[yyprhs[yyrule] + yyi],
		       &(yyvsp[(yyi + 1) - (yynrhs)])
		       		       );
      YYFPRINTF (stderr, "\n");
    }
}

# define YY_REDUCE_PRINT(Rule)		\
do {					\
  if (yydebug)				\
    yy_reduce_print (yyvsp, Rule); \
} while (YYID (0))

/* Nonzero means print parse trace.  It is left uninitialized so that
   multiple parsers can coexist.  */
int yydebug;
#else /* !YYDEBUG */
# define YYDPRINTF(Args)
# define YY_SYMBOL_PRINT(Title, Type, Value, Location)
# define YY_STACK_PRINT(Bottom, Top)
# define YY_REDUCE_PRINT(Rule)
#endif /* !YYDEBUG */


/* YYINITDEPTH -- initial size of the parser's stacks.  */
#ifndef	YYINITDEPTH
# define YYINITDEPTH 200
#endif

/* YYMAXDEPTH -- maximum size the stacks can grow to (effective only
   if the built-in stack extension method is used).

   Do not make this value too large; the results are undefined if
   YYSTACK_ALLOC_MAXIMUM < YYSTACK_BYTES (YYMAXDEPTH)
   evaluated with infinite-precision integer arithmetic.  */

#ifndef YYMAXDEPTH
# define YYMAXDEPTH 10000
#endif



#if YYERROR_VERBOSE

# ifndef yystrlen
#  if defined __GLIBC__ && defined _STRING_H
#   define yystrlen strlen
#  else
/* Return the length of YYSTR.  */
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static YYSIZE_T
yystrlen (const char *yystr)
#else
static YYSIZE_T
yystrlen (yystr)
    const char *yystr;
#endif
{
  YYSIZE_T yylen;
  for (yylen = 0; yystr[yylen]; yylen++)
    continue;
  return yylen;
}
#  endif
# endif

# ifndef yystpcpy
#  if defined __GLIBC__ && defined _STRING_H && defined _GNU_SOURCE
#   define yystpcpy stpcpy
#  else
/* Copy YYSRC to YYDEST, returning the address of the terminating '\0' in
   YYDEST.  */
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static char *
yystpcpy (char *yydest, const char *yysrc)
#else
static char *
yystpcpy (yydest, yysrc)
    char *yydest;
    const char *yysrc;
#endif
{
  char *yyd = yydest;
  const char *yys = yysrc;

  while ((*yyd++ = *yys++) != '\0')
    continue;

  return yyd - 1;
}
#  endif
# endif

# ifndef yytnamerr
/* Copy to YYRES the contents of YYSTR after stripping away unnecessary
   quotes and backslashes, so that it's suitable for yyerror.  The
   heuristic is that double-quoting is unnecessary unless the string
   contains an apostrophe, a comma, or backslash (other than
   backslash-backslash).  YYSTR is taken from yytname.  If YYRES is
   null, do not copy; instead, return the length of what the result
   would have been.  */
static YYSIZE_T
yytnamerr (char *yyres, const char *yystr)
{
  if (*yystr == '"')
    {
      YYSIZE_T yyn = 0;
      char const *yyp = yystr;

      for (;;)
	switch (*++yyp)
	  {
	  case '\'':
	  case ',':
	    goto do_not_strip_quotes;

	  case '\\':
	    if (*++yyp != '\\')
	      goto do_not_strip_quotes;
	    /* Fall through.  */
	  default:
	    if (yyres)
	      yyres[yyn] = *yyp;
	    yyn++;
	    break;

	  case '"':
	    if (yyres)
	      yyres[yyn] = '\0';
	    return yyn;
	  }
    do_not_strip_quotes: ;
    }

  if (! yyres)
    return yystrlen (yystr);

  return yystpcpy (yyres, yystr) - yyres;
}
# endif

/* Copy into YYRESULT an error message about the unexpected token
   YYCHAR while in state YYSTATE.  Return the number of bytes copied,
   including the terminating null byte.  If YYRESULT is null, do not
   copy anything; just return the number of bytes that would be
   copied.  As a special case, return 0 if an ordinary "syntax error"
   message will do.  Return YYSIZE_MAXIMUM if overflow occurs during
   size calculation.  */
static YYSIZE_T
yysyntax_error (char *yyresult, int yystate, int yychar)
{
  int yyn = yypact[yystate];

  if (! (YYPACT_NINF < yyn && yyn <= YYLAST))
    return 0;
  else
    {
      int yytype = YYTRANSLATE (yychar);
      YYSIZE_T yysize0 = yytnamerr (0, yytname[yytype]);
      YYSIZE_T yysize = yysize0;
      YYSIZE_T yysize1;
      int yysize_overflow = 0;
      enum { YYERROR_VERBOSE_ARGS_MAXIMUM = 5 };
      char const *yyarg[YYERROR_VERBOSE_ARGS_MAXIMUM];
      int yyx;

# if 0
      /* This is so xgettext sees the translatable formats that are
	 constructed on the fly.  */
      YY_("syntax error, unexpected %s");
      YY_("syntax error, unexpected %s, expecting %s");
      YY_("syntax error, unexpected %s, expecting %s or %s");
      YY_("syntax error, unexpected %s, expecting %s or %s or %s");
      YY_("syntax error, unexpected %s, expecting %s or %s or %s or %s");
# endif
      char *yyfmt;
      char const *yyf;
      static char const yyunexpected[] = "syntax error, unexpected %s";
      static char const yyexpecting[] = ", expecting %s";
      static char const yyor[] = " or %s";
      char yyformat[sizeof yyunexpected
		    + sizeof yyexpecting - 1
		    + ((YYERROR_VERBOSE_ARGS_MAXIMUM - 2)
		       * (sizeof yyor - 1))];
      char const *yyprefix = yyexpecting;

      /* Start YYX at -YYN if negative to avoid negative indexes in
	 YYCHECK.  */
      int yyxbegin = yyn < 0 ? -yyn : 0;

      /* Stay within bounds of both yycheck and yytname.  */
      int yychecklim = YYLAST - yyn + 1;
      int yyxend = yychecklim < YYNTOKENS ? yychecklim : YYNTOKENS;
      int yycount = 1;

      yyarg[0] = yytname[yytype];
      yyfmt = yystpcpy (yyformat, yyunexpected);

      for (yyx = yyxbegin; yyx < yyxend; ++yyx)
	if (yycheck[yyx + yyn] == yyx && yyx != YYTERROR)
	  {
	    if (yycount == YYERROR_VERBOSE_ARGS_MAXIMUM)
	      {
		yycount = 1;
		yysize = yysize0;
		yyformat[sizeof yyunexpected - 1] = '\0';
		break;
	      }
	    yyarg[yycount++] = yytname[yyx];
	    yysize1 = yysize + yytnamerr (0, yytname[yyx]);
	    yysize_overflow |= (yysize1 < yysize);
	    yysize = yysize1;
	    yyfmt = yystpcpy (yyfmt, yyprefix);
	    yyprefix = yyor;
	  }

      yyf = YY_(yyformat);
      yysize1 = yysize + yystrlen (yyf);
      yysize_overflow |= (yysize1 < yysize);
      yysize = yysize1;

      if (yysize_overflow)
	return YYSIZE_MAXIMUM;

      if (yyresult)
	{
	  /* Avoid sprintf, as that infringes on the user's name space.
	     Don't have undefined behavior even if the translation
	     produced a string with the wrong number of "%s"s.  */
	  char *yyp = yyresult;
	  int yyi = 0;
	  while ((*yyp = *yyf) != '\0')
	    {
	      if (*yyp == '%' && yyf[1] == 's' && yyi < yycount)
		{
		  yyp += yytnamerr (yyp, yyarg[yyi++]);
		  yyf += 2;
		}
	      else
		{
		  yyp++;
		  yyf++;
		}
	    }
	}
      return yysize;
    }
}
#endif /* YYERROR_VERBOSE */


/*-----------------------------------------------.
| Release the memory associated to this symbol.  |
`-----------------------------------------------*/

/*ARGSUSED*/
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yydestruct (const char *yymsg, int yytype, YYSTYPE *yyvaluep)
#else
static void
yydestruct (yymsg, yytype, yyvaluep)
    const char *yymsg;
    int yytype;
    YYSTYPE *yyvaluep;
#endif
{
  YYUSE (yyvaluep);

  if (!yymsg)
    yymsg = "Deleting";
  YY_SYMBOL_PRINT (yymsg, yytype, yyvaluep, yylocationp);

  switch (yytype)
    {

      default:
	break;
    }
}

/* Prevent warnings from -Wmissing-prototypes.  */
#ifdef YYPARSE_PARAM
#if defined __STDC__ || defined __cplusplus
int yyparse (void *YYPARSE_PARAM);
#else
int yyparse ();
#endif
#else /* ! YYPARSE_PARAM */
#if defined __STDC__ || defined __cplusplus
int yyparse (void);
#else
int yyparse ();
#endif
#endif /* ! YYPARSE_PARAM */





/*-------------------------.
| yyparse or yypush_parse.  |
`-------------------------*/

#ifdef YYPARSE_PARAM
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
int
yyparse (void *YYPARSE_PARAM)
#else
int
yyparse (YYPARSE_PARAM)
    void *YYPARSE_PARAM;
#endif
#else /* ! YYPARSE_PARAM */
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
int
yyparse (void)
#else
int
yyparse ()

#endif
#endif
{
/* The lookahead symbol.  */
int yychar;

/* The semantic value of the lookahead symbol.  */
YYSTYPE yylval;

    /* Number of syntax errors so far.  */
    int yynerrs;

    int yystate;
    /* Number of tokens to shift before error messages enabled.  */
    int yyerrstatus;

    /* The stacks and their tools:
       `yyss': related to states.
       `yyvs': related to semantic values.

       Refer to the stacks thru separate pointers, to allow yyoverflow
       to reallocate them elsewhere.  */

    /* The state stack.  */
    yytype_int16 yyssa[YYINITDEPTH];
    yytype_int16 *yyss;
    yytype_int16 *yyssp;

    /* The semantic value stack.  */
    YYSTYPE yyvsa[YYINITDEPTH];
    YYSTYPE *yyvs;
    YYSTYPE *yyvsp;

    YYSIZE_T yystacksize;

  int yyn;
  int yyresult;
  /* Lookahead token as an internal (translated) token number.  */
  int yytoken;
  /* The variables used to return semantic value and location from the
     action routines.  */
  YYSTYPE yyval;

#if YYERROR_VERBOSE
  /* Buffer for error messages, and its allocated size.  */
  char yymsgbuf[128];
  char *yymsg = yymsgbuf;
  YYSIZE_T yymsg_alloc = sizeof yymsgbuf;
#endif

#define YYPOPSTACK(N)   (yyvsp -= (N), yyssp -= (N))

  /* The number of symbols on the RHS of the reduced rule.
     Keep to zero when no symbol should be popped.  */
  int yylen = 0;

  yytoken = 0;
  yyss = yyssa;
  yyvs = yyvsa;
  yystacksize = YYINITDEPTH;

  YYDPRINTF ((stderr, "Starting parse\n"));

  yystate = 0;
  yyerrstatus = 0;
  yynerrs = 0;
  yychar = YYEMPTY; /* Cause a token to be read.  */

  /* Initialize stack pointers.
     Waste one element of value and location stack
     so that they stay on the same level as the state stack.
     The wasted elements are never initialized.  */
  yyssp = yyss;
  yyvsp = yyvs;

  goto yysetstate;

/*------------------------------------------------------------.
| yynewstate -- Push a new state, which is found in yystate.  |
`------------------------------------------------------------*/
 yynewstate:
  /* In all cases, when you get here, the value and location stacks
     have just been pushed.  So pushing a state here evens the stacks.  */
  yyssp++;

 yysetstate:
  *yyssp = yystate;

  if (yyss + yystacksize - 1 <= yyssp)
    {
      /* Get the current used size of the three stacks, in elements.  */
      YYSIZE_T yysize = yyssp - yyss + 1;

#ifdef yyoverflow
      {
	/* Give user a chance to reallocate the stack.  Use copies of
	   these so that the &'s don't force the real ones into
	   memory.  */
	YYSTYPE *yyvs1 = yyvs;
	yytype_int16 *yyss1 = yyss;

	/* Each stack pointer address is followed by the size of the
	   data in use in that stack, in bytes.  This used to be a
	   conditional around just the two extra args, but that might
	   be undefined if yyoverflow is a macro.  */
	yyoverflow (YY_("memory exhausted"),
		    &yyss1, yysize * sizeof (*yyssp),
		    &yyvs1, yysize * sizeof (*yyvsp),
		    &yystacksize);

	yyss = yyss1;
	yyvs = yyvs1;
      }
#else /* no yyoverflow */
# ifndef YYSTACK_RELOCATE
      goto yyexhaustedlab;
# else
      /* Extend the stack our own way.  */
      if (YYMAXDEPTH <= yystacksize)
	goto yyexhaustedlab;
      yystacksize *= 2;
      if (YYMAXDEPTH < yystacksize)
	yystacksize = YYMAXDEPTH;

      {
	yytype_int16 *yyss1 = yyss;
	union yyalloc *yyptr =
	  (union yyalloc *) YYSTACK_ALLOC (YYSTACK_BYTES (yystacksize));
	if (! yyptr)
	  goto yyexhaustedlab;
	YYSTACK_RELOCATE (yyss_alloc, yyss);
	YYSTACK_RELOCATE (yyvs_alloc, yyvs);
#  undef YYSTACK_RELOCATE
	if (yyss1 != yyssa)
	  YYSTACK_FREE (yyss1);
      }
# endif
#endif /* no yyoverflow */

      yyssp = yyss + yysize - 1;
      yyvsp = yyvs + yysize - 1;

      YYDPRINTF ((stderr, "Stack size increased to %lu\n",
		  (unsigned long int) yystacksize));

      if (yyss + yystacksize - 1 <= yyssp)
	YYABORT;
    }

  YYDPRINTF ((stderr, "Entering state %d\n", yystate));

  if (yystate == YYFINAL)
    YYACCEPT;

  goto yybackup;

/*-----------.
| yybackup.  |
`-----------*/
yybackup:

  /* Do appropriate processing given the current state.  Read a
     lookahead token if we need one and don't already have one.  */

  /* First try to decide what to do without reference to lookahead token.  */
  yyn = yypact[yystate];
  if (yyn == YYPACT_NINF)
    goto yydefault;

  /* Not known => get a lookahead token if don't already have one.  */

  /* YYCHAR is either YYEMPTY or YYEOF or a valid lookahead symbol.  */
  if (yychar == YYEMPTY)
    {
      YYDPRINTF ((stderr, "Reading a token: "));
      yychar = YYLEX;
    }

  if (yychar <= YYEOF)
    {
      yychar = yytoken = YYEOF;
      YYDPRINTF ((stderr, "Now at end of input.\n"));
    }
  else
    {
      yytoken = YYTRANSLATE (yychar);
      YY_SYMBOL_PRINT ("Next token is", yytoken, &yylval, &yylloc);
    }

  /* If the proper action on seeing token YYTOKEN is to reduce or to
     detect an error, take that action.  */
  yyn += yytoken;
  if (yyn < 0 || YYLAST < yyn || yycheck[yyn] != yytoken)
    goto yydefault;
  yyn = yytable[yyn];
  if (yyn <= 0)
    {
      if (yyn == 0 || yyn == YYTABLE_NINF)
	goto yyerrlab;
      yyn = -yyn;
      goto yyreduce;
    }

  /* Count tokens shifted since error; after three, turn off error
     status.  */
  if (yyerrstatus)
    yyerrstatus--;

  /* Shift the lookahead token.  */
  YY_SYMBOL_PRINT ("Shifting", yytoken, &yylval, &yylloc);

  /* Discard the shifted token.  */
  yychar = YYEMPTY;

  yystate = yyn;
  *++yyvsp = yylval;

  goto yynewstate;


/*-----------------------------------------------------------.
| yydefault -- do the default action for the current state.  |
`-----------------------------------------------------------*/
yydefault:
  yyn = yydefact[yystate];
  if (yyn == 0)
    goto yyerrlab;
  goto yyreduce;


/*-----------------------------.
| yyreduce -- Do a reduction.  |
`-----------------------------*/
yyreduce:
  /* yyn is the number of a rule to reduce with.  */
  yylen = yyr2[yyn];

  /* If YYLEN is nonzero, implement the default value of the action:
     `$$ = $1'.

     Otherwise, the following line sets YYVAL to garbage.
     This behavior is undocumented and Bison
     users should not rely upon it.  Assigning to YYVAL
     unconditionally makes the parser a bit smaller, and it avoids a
     GCC warning that YYVAL may be used uninitialized.  */
  yyval = yyvsp[1-yylen];


  YY_REDUCE_PRINT (yyn);
  switch (yyn)
    {
        case 6:

/* Line 1455 of yacc.c  */
#line 147 "parse.y"
    {copy_filepos(&handle->start_pos,&handle->begin_tok);}
    break;

  case 8:

/* Line 1455 of yacc.c  */
#line 151 "parse.y"
    {copy_filepos(&handle->start_pos,&handle->begin_tok);}
    break;

  case 9:

/* Line 1455 of yacc.c  */
#line 152 "parse.y"
    {copy_filepos(&handle->start_pos,&handle->begin_tok);}
    break;

  case 21:

/* Line 1455 of yacc.c  */
#line 169 "parse.y"
    {com_ctypeI((yyvsp[(1) - (2)].vv),0,(yyvsp[(2) - (2)].term),0);}
    break;

  case 22:

/* Line 1455 of yacc.c  */
#line 170 "parse.y"
    {com_ctypeI((yyvsp[(1) - (4)].vv),0,(yyvsp[(2) - (4)].term),(yyvsp[(4) - (4)].vl));}
    break;

  case 23:

/* Line 1455 of yacc.c  */
#line 171 "parse.y"
    {com_ctypeI((yyvsp[(1) - (2)].vv),0,0,(yyvsp[(2) - (2)].vl));}
    break;

  case 24:

/* Line 1455 of yacc.c  */
#line 172 "parse.y"
    {com_ctypeI((yyvsp[(1) - (3)].vv),(yyvsp[(2) - (3)].pc),(yyvsp[(3) - (3)].term),0);}
    break;

  case 25:

/* Line 1455 of yacc.c  */
#line 173 "parse.y"
    {com_ctypeI((yyvsp[(1) - (5)].vv),(yyvsp[(2) - (5)].pc),(yyvsp[(3) - (5)].term),(yyvsp[(5) - (5)].vl));}
    break;

  case 26:

/* Line 1455 of yacc.c  */
#line 174 "parse.y"
    {com_ctypeI((yyvsp[(1) - (3)].vv),(yyvsp[(2) - (3)].pc),0,(yyvsp[(3) - (3)].vl));}
    break;

  case 27:

/* Line 1455 of yacc.c  */
#line 175 "parse.y"
    {com_ctypeI(0,0,0,make_var_list(get_var_term((yyvsp[(1) - (2)].vv)),0));}
    break;

  case 28:

/* Line 1455 of yacc.c  */
#line 176 "parse.y"
    {com_ctypeI(0,0,0,make_var_list(get_var_term((yyvsp[(1) - (3)].vv)),(yyvsp[(3) - (3)].pc)));}
    break;

  case 29:

/* Line 1455 of yacc.c  */
#line 177 "parse.y"
    {com_ctypeI(0,0,0,merge_var_varlist((yyvsp[(1) - (4)].vv),0,(yyvsp[(4) - (4)].vl)));}
    break;

  case 30:

/* Line 1455 of yacc.c  */
#line 178 "parse.y"
    {com_ctypeI(0,0,0,merge_var_varlist((yyvsp[(1) - (5)].vv),(yyvsp[(3) - (5)].pc),(yyvsp[(5) - (5)].vl)));}
    break;

  case 31:

/* Line 1455 of yacc.c  */
#line 179 "parse.y"
    {com_ctypeI(0,0,0,(yyvsp[(1) - (1)].vl));}
    break;

  case 32:

/* Line 1455 of yacc.c  */
#line 180 "parse.y"
    {com_ctypeI(0,0,0,addto_var_list((yyvsp[(1) - (3)].vl),(yyvsp[(3) - (3)].vl)));}
    break;

  case 33:

/* Line 1455 of yacc.c  */
#line 181 "parse.y"
    {com_ctypeI(0,0,0,(yyvsp[(1) - (1)].vl));}
    break;

  case 34:

/* Line 1455 of yacc.c  */
#line 182 "parse.y"
    {com_ctypeI(0,0,0,addto_var_list((yyvsp[(1) - (3)].vl),(yyvsp[(3) - (3)].vl)));}
    break;

  case 35:

/* Line 1455 of yacc.c  */
#line 186 "parse.y"
    {com_ctypeII(-1,(yyvsp[(1) - (3)].vv),(yyvsp[(2) - (3)].term),(yyvsp[(3) - (3)].vl));}
    break;

  case 36:

/* Line 1455 of yacc.c  */
#line 187 "parse.y"
    {com_ctypeII((yyvsp[(1) - (4)].i),(yyvsp[(2) - (4)].vv),(yyvsp[(3) - (4)].term),(yyvsp[(4) - (4)].vl));}
    break;

  case 37:

/* Line 1455 of yacc.c  */
#line 190 "parse.y"
    {if((yyvsp[(1) - (3)].term) && (yyvsp[(3) - (3)].term)) com_assign((yyvsp[(1) - (3)].term),(yyvsp[(3) - (3)].term));}
    break;

  case 38:

/* Line 1455 of yacc.c  */
#line 194 "parse.y"
    {com_ctypeIII((yyvsp[(1) - (4)].vv),(yyvsp[(2) - (4)].term),(yyvsp[(4) - (4)].term));}
    break;

  case 40:

/* Line 1455 of yacc.c  */
#line 196 "parse.y"
    {free_parse_term((yyvsp[(3) - (3)].term));}
    break;

  case 41:

/* Line 1455 of yacc.c  */
#line 197 "parse.y"
    {free_parse_term(parse_incr_decr((yyvsp[(1) - (2)].term),PREFLAG|((yyvsp[(2) - (2)].i))));}
    break;

  case 42:

/* Line 1455 of yacc.c  */
#line 198 "parse.y"
    {free_parse_term(parse_incr_decr((yyvsp[(2) - (2)].term),PREFLAG|((yyvsp[(1) - (2)].i))));}
    break;

  case 43:

/* Line 1455 of yacc.c  */
#line 199 "parse.y"
    {vtmp=do_op(copy_parse_term((yyvsp[(1) - (3)].term)),(yyvsp[(3) - (3)].term),(yyvsp[(2) - (3)].i));
		                           if((yyvsp[(1) - (3)].term) && vtmp) com_assign((yyvsp[(1) - (3)].term),vtmp);}
    break;

  case 44:

/* Line 1455 of yacc.c  */
#line 204 "parse.y"
    {}
    break;

  case 45:

/* Line 1455 of yacc.c  */
#line 205 "parse.y"
    {}
    break;

  case 46:

/* Line 1455 of yacc.c  */
#line 206 "parse.y"
    {}
    break;

  case 47:

/* Line 1455 of yacc.c  */
#line 207 "parse.y"
    {}
    break;

  case 48:

/* Line 1455 of yacc.c  */
#line 211 "parse.y"
    {com_ctypeV((yyvsp[(1) - (3)].vv),(yyvsp[(2) - (3)].i),(yyvsp[(3) - (3)].vl));}
    break;

  case 49:

/* Line 1455 of yacc.c  */
#line 215 "parse.y"
    {parerror1(key_err,key_tab[(yyvsp[(1) - (2)].i)]);}
    break;

  case 51:

/* Line 1455 of yacc.c  */
#line 216 "parse.y"
    {parerror1(key_err,key_tab[(yyvsp[(1) - (4)].i)]);}
    break;

  case 53:

/* Line 1455 of yacc.c  */
#line 217 "parse.y"
    {parerror1(key_err,key_tab[(yyvsp[(2) - (3)].i)]);}
    break;

  case 55:

/* Line 1455 of yacc.c  */
#line 218 "parse.y"
    {parerror1(key_err,key_tab[(yyvsp[(2) - (5)].i)]);}
    break;

  case 57:

/* Line 1455 of yacc.c  */
#line 219 "parse.y"
    {}
    break;

  case 58:

/* Line 1455 of yacc.c  */
#line 220 "parse.y"
    {parerror2("Syntax error\n");}
    break;

  case 72:

/* Line 1455 of yacc.c  */
#line 238 "parse.y"
    {(yyval.fc)=make_format_clause((yyvsp[(1) - (1)].i),0);}
    break;

  case 73:

/* Line 1455 of yacc.c  */
#line 239 "parse.y"
    {(yyval.fc)=make_format_clause(-1,0);}
    break;

  case 74:

/* Line 1455 of yacc.c  */
#line 240 "parse.y"
    {}
    break;

  case 75:

/* Line 1455 of yacc.c  */
#line 240 "parse.y"
    {(yyval.fc)=make_format_clause(-(yyvsp[(1) - (3)].i),0);}
    break;

  case 76:

/* Line 1455 of yacc.c  */
#line 241 "parse.y"
    {(yyval.fc)=make_format_clause((yyvsp[(1) - (4)].i),reverse_list((yyvsp[(3) - (4)].fc)));}
    break;

  case 77:

/* Line 1455 of yacc.c  */
#line 242 "parse.y"
    {(yyval.fc)=(yyvsp[(2) - (3)].fc);}
    break;

  case 81:

/* Line 1455 of yacc.c  */
#line 250 "parse.y"
    {(yyvsp[(3) - (3)].fc)->next=(yyvsp[(1) - (3)].fc); (yyval.fc)=(yyvsp[(3) - (3)].fc);}
    break;

  case 82:

/* Line 1455 of yacc.c  */
#line 253 "parse.y"
    {(yyval.pc)=make_clause((yyvsp[(1) - (1)].term));}
    break;

  case 83:

/* Line 1455 of yacc.c  */
#line 254 "parse.y"
    {(yyval.pc)=convert_format_clause(reverse_list((yyvsp[(1) - (1)].fc)));}
    break;

  case 85:

/* Line 1455 of yacc.c  */
#line 258 "parse.y"
    {(yyvsp[(3) - (3)].pc)->next=(yyvsp[(1) - (3)].pc); (yyval.pc)=(yyvsp[(3) - (3)].pc);}
    break;

  case 86:

/* Line 1455 of yacc.c  */
#line 261 "parse.y"
    {(yyval.pc)=reverse_list((yyvsp[(2) - (3)].pc)); in_clause=0;}
    break;

  case 87:

/* Line 1455 of yacc.c  */
#line 262 "parse.y"
    {(yyval.pc)=0; in_clause=0;}
    break;

  case 88:

/* Line 1455 of yacc.c  */
#line 265 "parse.y"
    {in_clause=1;}
    break;

  case 95:

/* Line 1455 of yacc.c  */
#line 282 "parse.y"
    {}
    break;

  case 96:

/* Line 1455 of yacc.c  */
#line 285 "parse.y"
    {include_file((yyvsp[(2) - (3)].str));}
    break;

  case 98:

/* Line 1455 of yacc.c  */
#line 289 "parse.y"
    {(yyval.vl)=addto_var_list((yyvsp[(1) - (3)].vl),(yyvsp[(3) - (3)].vl));}
    break;

  case 99:

/* Line 1455 of yacc.c  */
#line 292 "parse.y"
    {(yyval.vl)=make_var_list((yyvsp[(1) - (1)].term),0); copy_filepos(&((yyval.vl)->pos),&handle->begin_tok);}
    break;

  case 102:

/* Line 1455 of yacc.c  */
#line 297 "parse.y"
    {(yyval.term)=(yyvsp[(2) - (2)].term);}
    break;

  case 104:

/* Line 1455 of yacc.c  */
#line 299 "parse.y"
    {(yyval.term)=parse_make_texp((yyvsp[(1) - (1)].i),0,0,INTEGER);}
    break;

  case 107:

/* Line 1455 of yacc.c  */
#line 306 "parse.y"
    {(yyval.term)=get_var_term((yyvsp[(1) - (1)].vv));}
    break;

  case 108:

/* Line 1455 of yacc.c  */
#line 307 "parse.y"
    {(yyval.term)=get_array_term((yyvsp[(1) - (3)].vv),(yyvsp[(2) - (3)].term));}
    break;

  case 109:

/* Line 1455 of yacc.c  */
#line 310 "parse.y"
    {(yyval.vl)=make_var_list(get_var_term((yyvsp[(1) - (1)].vv)),0);}
    break;

  case 110:

/* Line 1455 of yacc.c  */
#line 311 "parse.y"
    {(yyval.vl)=get_array_var_list((yyvsp[(1) - (3)].vv),(yyvsp[(2) - (3)].term));}
    break;

  case 111:

/* Line 1455 of yacc.c  */
#line 312 "parse.y"
    {(yyval.vl)=make_var_list(0,0);}
    break;

  case 112:

/* Line 1455 of yacc.c  */
#line 315 "parse.y"
    {(yyval.vl)=get_array_var_list((yyvsp[(1) - (4)].vv),(yyvsp[(3) - (4)].term));}
    break;

  case 113:

/* Line 1455 of yacc.c  */
#line 318 "parse.y"
    {(yyval.vl)=make_loop_clause((yyvsp[(3) - (8)].vl),(yyvsp[(5) - (8)].vv),(yyvsp[(7) - (8)].term));}
    break;

  case 114:

/* Line 1455 of yacc.c  */
#line 319 "parse.y"
    {(yyval.vl)=make_loop_clause((yyvsp[(3) - (4)].vl),0,0);}
    break;

  case 115:

/* Line 1455 of yacc.c  */
#line 322 "parse.y"
    {(yyval.vl)=make_loop_clause((yyvsp[(2) - (7)].vl),(yyvsp[(4) - (7)].vv),(yyvsp[(6) - (7)].term));}
    break;

  case 116:

/* Line 1455 of yacc.c  */
#line 323 "parse.y"
    {(yyval.vl)=make_loop_clause((yyvsp[(2) - (3)].vl),0,0);}
    break;

  case 118:

/* Line 1455 of yacc.c  */
#line 327 "parse.y"
    {(yyval.vl)=addto_var_list((yyvsp[(1) - (3)].vl),(yyvsp[(3) - (3)].vl));}
    break;

  case 119:

/* Line 1455 of yacc.c  */
#line 328 "parse.y"
    {parerror1(key_err,key_tab[(yyvsp[(1) - (3)].i)]); (yyval.vl)=0;}
    break;

  case 120:

/* Line 1455 of yacc.c  */
#line 329 "parse.y"
    {parerror1(key_err,key_tab[(yyvsp[(1) - (1)].i)]); (yyval.vl)=0;}
    break;

  case 121:

/* Line 1455 of yacc.c  */
#line 330 "parse.y"
    {parerror1(key_err,key_tab[(yyvsp[(3) - (5)].i)]); (yyval.vl)=0;}
    break;

  case 122:

/* Line 1455 of yacc.c  */
#line 331 "parse.y"
    {parerror1(key_err,key_tab[(yyvsp[(3) - (3)].i)]); (yyval.vl)=0;}
    break;

  case 123:

/* Line 1455 of yacc.c  */
#line 334 "parse.y"
    {(yyval.vl)=make_var_list((yyvsp[(1) - (1)].term),0); copy_filepos(&((yyval.vl)->pos),&handle->begin_tok);}
    break;

  case 124:

/* Line 1455 of yacc.c  */
#line 335 "parse.y"
    {(yyval.vl)=make_var_list((yyvsp[(1) - (2)].term),(yyvsp[(2) - (2)].pc));}
    break;

  case 127:

/* Line 1455 of yacc.c  */
#line 342 "parse.y"
    {(yyval.vl)=make_var_list(0,0);}
    break;

  case 128:

/* Line 1455 of yacc.c  */
#line 343 "parse.y"
    {(yyval.vl)=addto_var_list((yyvsp[(1) - (2)].vl),make_var_list(0,0));}
    break;

  case 130:

/* Line 1455 of yacc.c  */
#line 347 "parse.y"
    {(yyval.vl)=addto_var_list(make_var_list(0,0),(yyvsp[(2) - (2)].vl));}
    break;

  case 131:

/* Line 1455 of yacc.c  */
#line 348 "parse.y"
    {(yyval.vl)=addto_var_list(addto_var_list((yyvsp[(1) - (2)].vl),make_var_list(0,0)),(yyvsp[(2) - (2)].vl));}
    break;

  case 132:

/* Line 1455 of yacc.c  */
#line 349 "parse.y"
    {(yyval.vl)=addto_var_list((yyvsp[(1) - (3)].vl),(yyvsp[(3) - (3)].vl));}
    break;

  case 133:

/* Line 1455 of yacc.c  */
#line 350 "parse.y"
    {(yyval.vl)=addto_var_list(addto_var_list((yyvsp[(1) - (3)].vl),(yyvsp[(2) - (3)].vl)),(yyvsp[(3) - (3)].vl));}
    break;

  case 135:

/* Line 1455 of yacc.c  */
#line 354 "parse.y"
    {(yyval.term)=do_op1((yyvsp[(1) - (3)].term),(yyvsp[(3) - (3)].term),'+');}
    break;

  case 136:

/* Line 1455 of yacc.c  */
#line 355 "parse.y"
    {(yyval.term)=do_op1((yyvsp[(1) - (3)].term),(yyvsp[(3) - (3)].term),'.');}
    break;

  case 137:

/* Line 1455 of yacc.c  */
#line 356 "parse.y"
    {(yyval.term)=do_op1((yyvsp[(1) - (3)].term),(yyvsp[(3) - (3)].term),'-');}
    break;

  case 138:

/* Line 1455 of yacc.c  */
#line 357 "parse.y"
    {(yyval.term)=do_op1((yyvsp[(1) - (3)].term),(yyvsp[(3) - (3)].term),'/');}
    break;

  case 139:

/* Line 1455 of yacc.c  */
#line 358 "parse.y"
    {(yyval.term)=do_op1((yyvsp[(1) - (3)].term),(yyvsp[(3) - (3)].term),'*');}
    break;

  case 140:

/* Line 1455 of yacc.c  */
#line 359 "parse.y"
    {(yyval.term)=do_op1((yyvsp[(1) - (3)].term),(yyvsp[(3) - (3)].term),'<');}
    break;

  case 141:

/* Line 1455 of yacc.c  */
#line 360 "parse.y"
    {(yyval.term)=do_op1((yyvsp[(1) - (3)].term),(yyvsp[(3) - (3)].term),'>');}
    break;

  case 142:

/* Line 1455 of yacc.c  */
#line 361 "parse.y"
    {(yyval.term)=do_op1((yyvsp[(1) - (3)].term),(yyvsp[(3) - (3)].term),'&');}
    break;

  case 143:

/* Line 1455 of yacc.c  */
#line 362 "parse.y"
    {(yyval.term)=do_op1((yyvsp[(1) - (3)].term),(yyvsp[(3) - (3)].term),'|');}
    break;

  case 144:

/* Line 1455 of yacc.c  */
#line 363 "parse.y"
    {(yyval.term)=do_op1((yyvsp[(1) - (3)].term),(yyvsp[(3) - (3)].term),EQ_TOKEN);}
    break;

  case 145:

/* Line 1455 of yacc.c  */
#line 364 "parse.y"
    {(yyval.term)=do_op1((yyvsp[(1) - (3)].term),(yyvsp[(3) - (3)].term),NE_TOKEN);}
    break;

  case 146:

/* Line 1455 of yacc.c  */
#line 365 "parse.y"
    {(yyval.term)=do_op1((yyvsp[(1) - (3)].term),(yyvsp[(3) - (3)].term),LEQ_TOKEN);}
    break;

  case 147:

/* Line 1455 of yacc.c  */
#line 366 "parse.y"
    {(yyval.term)=do_op1((yyvsp[(1) - (3)].term),(yyvsp[(3) - (3)].term),GEQ_TOKEN);}
    break;

  case 148:

/* Line 1455 of yacc.c  */
#line 367 "parse.y"
    {(yyval.term)=do_op1((yyvsp[(1) - (3)].term),(yyvsp[(3) - (3)].term),EXPO);}
    break;

  case 149:

/* Line 1455 of yacc.c  */
#line 368 "parse.y"
    {(yyval.term)=do_op1((yyvsp[(2) - (2)].term),0,'!');}
    break;

  case 150:

/* Line 1455 of yacc.c  */
#line 369 "parse.y"
    {(yyval.term)=check_immediate((yyvsp[(2) - (3)].term));}
    break;

  case 151:

/* Line 1455 of yacc.c  */
#line 370 "parse.y"
    {(yyval.term)=do_func((yyvsp[(1) - (4)].f),(yyvsp[(3) - (4)].term));}
    break;

  case 155:

/* Line 1455 of yacc.c  */
#line 378 "parse.y"
    {(yyval.str)=add_strings((yyvsp[(1) - (3)].str),(yyvsp[(3) - (3)].str));}
    break;

  case 156:

/* Line 1455 of yacc.c  */
#line 379 "parse.y"
    {(yyval.str)=add_strings(check_string((yyvsp[(1) - (3)].term)),(yyvsp[(3) - (3)].str));}
    break;

  case 158:

/* Line 1455 of yacc.c  */
#line 383 "parse.y"
    {(yyval.str)=(yyvsp[(1) - (4)].f)->f.f_string((yyvsp[(3) - (4)].str));}
    break;

  case 160:

/* Line 1455 of yacc.c  */
#line 387 "parse.y"
    {(yyval.str)=check_string((yyvsp[(1) - (1)].term));}
    break;

  case 161:

/* Line 1455 of yacc.c  */
#line 390 "parse.y"
    {if((yyvsp[(1) - (1)].term)->defer) do_deferred((yyvsp[(1) - (1)].term)); (yyval.term)=(yyvsp[(1) - (1)].term); }
    break;

  case 162:

/* Line 1455 of yacc.c  */
#line 391 "parse.y"
    {if((yyvsp[(3) - (3)].term)->defer) do_deferred((yyvsp[(3) - (3)].term)); (yyval.term)=do_op((yyvsp[(1) - (3)].term),(yyvsp[(3) - (3)].term),',');}
    break;

  case 164:

/* Line 1455 of yacc.c  */
#line 395 "parse.y"
    {(yyval.term)=do_op2((yyvsp[(1) - (3)].term),(yyvsp[(3) - (3)].term),SUB_EXPR,in_clause);}
    break;

  case 165:

/* Line 1455 of yacc.c  */
#line 398 "parse.y"
    {(yyval.term)=do_op((yyvsp[(1) - (3)].term),(yyvsp[(3) - (3)].term),'+');}
    break;

  case 166:

/* Line 1455 of yacc.c  */
#line 399 "parse.y"
    {(yyval.term)=do_op((yyvsp[(1) - (3)].term),(yyvsp[(3) - (3)].term),'-');}
    break;

  case 167:

/* Line 1455 of yacc.c  */
#line 400 "parse.y"
    {(yyval.term)=do_op((yyvsp[(1) - (3)].term),(yyvsp[(3) - (3)].term),'.');}
    break;

  case 168:

/* Line 1455 of yacc.c  */
#line 401 "parse.y"
    {(yyval.term)=do_op((yyvsp[(1) - (3)].term),(yyvsp[(3) - (3)].term),'*');}
    break;

  case 169:

/* Line 1455 of yacc.c  */
#line 402 "parse.y"
    {(yyval.term)=do_op((yyvsp[(1) - (3)].term),(yyvsp[(3) - (3)].term),'/');}
    break;

  case 170:

/* Line 1455 of yacc.c  */
#line 403 "parse.y"
    {(yyval.term)=do_op2((yyvsp[(1) - (3)].term),(yyvsp[(3) - (3)].term),'=',in_clause);}
    break;

  case 171:

/* Line 1455 of yacc.c  */
#line 404 "parse.y"
    {(yyval.term)=do_op((yyvsp[(1) - (3)].term),(yyvsp[(3) - (3)].term),EQ_TOKEN);}
    break;

  case 172:

/* Line 1455 of yacc.c  */
#line 405 "parse.y"
    {(yyval.term)=do_op((yyvsp[(1) - (3)].term),(yyvsp[(3) - (3)].term),'<');}
    break;

  case 173:

/* Line 1455 of yacc.c  */
#line 406 "parse.y"
    {(yyval.term)=do_op((yyvsp[(1) - (3)].term),(yyvsp[(3) - (3)].term),'>');}
    break;

  case 174:

/* Line 1455 of yacc.c  */
#line 407 "parse.y"
    {(yyval.term)=do_op((yyvsp[(1) - (3)].term),(yyvsp[(3) - (3)].term),'|');}
    break;

  case 175:

/* Line 1455 of yacc.c  */
#line 408 "parse.y"
    {(yyval.term)=do_op((yyvsp[(1) - (3)].term),(yyvsp[(3) - (3)].term),'&');}
    break;

  case 176:

/* Line 1455 of yacc.c  */
#line 409 "parse.y"
    {(yyval.term)=do_op((yyvsp[(1) - (3)].term),(yyvsp[(3) - (3)].term),EXPO);}
    break;

  case 177:

/* Line 1455 of yacc.c  */
#line 410 "parse.y"
    {(yyval.term)=do_op((yyvsp[(1) - (3)].term),(yyvsp[(3) - (3)].term),LEQ_TOKEN);}
    break;

  case 178:

/* Line 1455 of yacc.c  */
#line 411 "parse.y"
    {(yyval.term)=do_op((yyvsp[(1) - (3)].term),(yyvsp[(3) - (3)].term),GEQ_TOKEN);}
    break;

  case 179:

/* Line 1455 of yacc.c  */
#line 412 "parse.y"
    {(yyval.term)=do_op((yyvsp[(1) - (3)].term),(yyvsp[(3) - (3)].term),NE_TOKEN);}
    break;

  case 180:

/* Line 1455 of yacc.c  */
#line 413 "parse.y"
    {
		  if(IS_VAR((yyvsp[(1) - (3)].term))) {
			  vtmp=do_op(copy_parse_term((yyvsp[(1) - (3)].term)),(yyvsp[(3) - (3)].term),(yyvsp[(2) - (3)].i)); if((yyvsp[(1) - (3)].term) && vtmp) com_assign(copy_parse_term((yyvsp[(1) - (3)].term)),vtmp); (yyval.term)=(yyvsp[(1) - (3)].term);
		  } else {
			  parerror1("Can't assign to constant: '=' used instead of '==' perhaps?\n");
			  free_parse_term((yyvsp[(3) - (3)].term));
			  (yyval.term)=(yyvsp[(1) - (3)].term);
		  } }
    break;

  case 181:

/* Line 1455 of yacc.c  */
#line 421 "parse.y"
    {(yyval.term)=do_op((yyvsp[(2) - (2)].term),0,'-');}
    break;

  case 182:

/* Line 1455 of yacc.c  */
#line 422 "parse.y"
    {(yyval.term)=do_op((yyvsp[(2) - (2)].term),0,'!');}
    break;

  case 183:

/* Line 1455 of yacc.c  */
#line 423 "parse.y"
    {(yyval.term)=(yyvsp[(2) - (3)].term);}
    break;

  case 186:

/* Line 1455 of yacc.c  */
#line 428 "parse.y"
    {parerror1(key_err,key_tab[(yyvsp[(1) - (1)].i)]);(yyval.term)=0;}
    break;

  case 187:

/* Line 1455 of yacc.c  */
#line 429 "parse.y"
    {parerror1(key_err,key_tab[(yyvsp[(1) - (1)].i)]);(yyval.term)=0;}
    break;

  case 189:

/* Line 1455 of yacc.c  */
#line 431 "parse.y"
    {(yyval.term)=parse_incr_decr((yyvsp[(1) - (2)].term),POSTFLAG|((yyvsp[(2) - (2)].i)));}
    break;

  case 190:

/* Line 1455 of yacc.c  */
#line 432 "parse.y"
    {(yyval.term)=parse_incr_decr((yyvsp[(2) - (2)].term),PREFLAG|((yyvsp[(1) - (2)].i)));}
    break;

  case 191:

/* Line 1455 of yacc.c  */
#line 433 "parse.y"
    {(yyval.term)=do_func((yyvsp[(1) - (4)].f),(yyvsp[(3) - (4)].term));}
    break;

  case 192:

/* Line 1455 of yacc.c  */
#line 436 "parse.y"
    {(yyval.term)=parse_make_texp((yyvsp[(1) - (1)].i),0,0,INTEGER);}
    break;

  case 193:

/* Line 1455 of yacc.c  */
#line 437 "parse.y"
    {(yyval.term)=parse_make_texp(0,(yyvsp[(1) - (1)].x),0,REAL);}
    break;

  case 194:

/* Line 1455 of yacc.c  */
#line 438 "parse.y"
    {(yyval.term)=parse_make_texp(0,0,(yyvsp[(1) - (1)].str),STRING);}
    break;

  case 196:

/* Line 1455 of yacc.c  */
#line 442 "parse.y"
    {(yyval.term)=check_immediate((yyvsp[(1) - (1)].term));}
    break;



/* Line 1455 of yacc.c  */
#line 2955 "y.tab.c"
      default: break;
    }
  YY_SYMBOL_PRINT ("-> $$ =", yyr1[yyn], &yyval, &yyloc);

  YYPOPSTACK (yylen);
  yylen = 0;
  YY_STACK_PRINT (yyss, yyssp);

  *++yyvsp = yyval;

  /* Now `shift' the result of the reduction.  Determine what state
     that goes to, based on the state we popped back to and the rule
     number reduced by.  */

  yyn = yyr1[yyn];

  yystate = yypgoto[yyn - YYNTOKENS] + *yyssp;
  if (0 <= yystate && yystate <= YYLAST && yycheck[yystate] == *yyssp)
    yystate = yytable[yystate];
  else
    yystate = yydefgoto[yyn - YYNTOKENS];

  goto yynewstate;


/*------------------------------------.
| yyerrlab -- here on detecting error |
`------------------------------------*/
yyerrlab:
  /* If not already recovering from an error, report this error.  */
  if (!yyerrstatus)
    {
      ++yynerrs;
#if ! YYERROR_VERBOSE
      yyerror (YY_("syntax error"));
#else
      {
	YYSIZE_T yysize = yysyntax_error (0, yystate, yychar);
	if (yymsg_alloc < yysize && yymsg_alloc < YYSTACK_ALLOC_MAXIMUM)
	  {
	    YYSIZE_T yyalloc = 2 * yysize;
	    if (! (yysize <= yyalloc && yyalloc <= YYSTACK_ALLOC_MAXIMUM))
	      yyalloc = YYSTACK_ALLOC_MAXIMUM;
	    if (yymsg != yymsgbuf)
	      YYSTACK_FREE (yymsg);
	    yymsg = (char *) YYSTACK_ALLOC (yyalloc);
	    if (yymsg)
	      yymsg_alloc = yyalloc;
	    else
	      {
		yymsg = yymsgbuf;
		yymsg_alloc = sizeof yymsgbuf;
	      }
	  }

	if (0 < yysize && yysize <= yymsg_alloc)
	  {
	    (void) yysyntax_error (yymsg, yystate, yychar);
	    yyerror (yymsg);
	  }
	else
	  {
	    yyerror (YY_("syntax error"));
	    if (yysize != 0)
	      goto yyexhaustedlab;
	  }
      }
#endif
    }



  if (yyerrstatus == 3)
    {
      /* If just tried and failed to reuse lookahead token after an
	 error, discard it.  */

      if (yychar <= YYEOF)
	{
	  /* Return failure if at end of input.  */
	  if (yychar == YYEOF)
	    YYABORT;
	}
      else
	{
	  yydestruct ("Error: discarding",
		      yytoken, &yylval);
	  yychar = YYEMPTY;
	}
    }

  /* Else will try to reuse lookahead token after shifting the error
     token.  */
  goto yyerrlab1;


/*---------------------------------------------------.
| yyerrorlab -- error raised explicitly by YYERROR.  |
`---------------------------------------------------*/
yyerrorlab:

  /* Pacify compilers like GCC when the user code never invokes
     YYERROR and the label yyerrorlab therefore never appears in user
     code.  */
  if (/*CONSTCOND*/ 0)
     goto yyerrorlab;

  /* Do not reclaim the symbols of the rule which action triggered
     this YYERROR.  */
  YYPOPSTACK (yylen);
  yylen = 0;
  YY_STACK_PRINT (yyss, yyssp);
  yystate = *yyssp;
  goto yyerrlab1;


/*-------------------------------------------------------------.
| yyerrlab1 -- common code for both syntax error and YYERROR.  |
`-------------------------------------------------------------*/
yyerrlab1:
  yyerrstatus = 3;	/* Each real token shifted decrements this.  */

  for (;;)
    {
      yyn = yypact[yystate];
      if (yyn != YYPACT_NINF)
	{
	  yyn += YYTERROR;
	  if (0 <= yyn && yyn <= YYLAST && yycheck[yyn] == YYTERROR)
	    {
	      yyn = yytable[yyn];
	      if (0 < yyn)
		break;
	    }
	}

      /* Pop the current state because it cannot handle the error token.  */
      if (yyssp == yyss)
	YYABORT;


      yydestruct ("Error: popping",
		  yystos[yystate], yyvsp);
      YYPOPSTACK (1);
      yystate = *yyssp;
      YY_STACK_PRINT (yyss, yyssp);
    }

  *++yyvsp = yylval;


  /* Shift the error token.  */
  YY_SYMBOL_PRINT ("Shifting", yystos[yyn], yyvsp, yylsp);

  yystate = yyn;
  goto yynewstate;


/*-------------------------------------.
| yyacceptlab -- YYACCEPT comes here.  |
`-------------------------------------*/
yyacceptlab:
  yyresult = 0;
  goto yyreturn;

/*-----------------------------------.
| yyabortlab -- YYABORT comes here.  |
`-----------------------------------*/
yyabortlab:
  yyresult = 1;
  goto yyreturn;

#if !defined(yyoverflow) || YYERROR_VERBOSE
/*-------------------------------------------------.
| yyexhaustedlab -- memory exhaustion comes here.  |
`-------------------------------------------------*/
yyexhaustedlab:
  yyerror (YY_("memory exhausted"));
  yyresult = 2;
  /* Fall through.  */
#endif

yyreturn:
  if (yychar != YYEMPTY)
     yydestruct ("Cleanup: discarding lookahead",
		 yytoken, &yylval);
  /* Do not reclaim the symbols of the rule which action triggered
     this YYABORT or YYACCEPT.  */
  YYPOPSTACK (yylen);
  YY_STACK_PRINT (yyss, yyssp);
  while (yyssp != yyss)
    {
      yydestruct ("Cleanup: popping",
		  yystos[*yyssp], yyvsp);
      YYPOPSTACK (1);
    }
#ifndef yyoverflow
  if (yyss != yyssa)
    YYSTACK_FREE (yyss);
#endif
#if YYERROR_VERBOSE
  if (yymsg != yymsgbuf)
    YYSTACK_FREE (yymsg);
#endif
  /* Make sure YYID is used.  */
  return YYID (yyresult);
}



/* Line 1675 of yacc.c  */
#line 444 "parse.y"


void set_typeI_kludge(void)
{
	typeI_kludge=1;
}

static void set_filepos(filepos *fp)
{
	fp->fname=lex_info->name;
	fp->line=lex_info->line_no;
	fp->col=lex_info->col_no;
	if(!fp->col) fp->col++;
}

static void copy_filepos(filepos *fp1,filepos *fp2)
{
	memcpy(fp1,fp2,sizeof(filepos));
}

static void include_file(string *s)
{
	struct lex_info *li;
	struct stat st1,st2;
	
	if(s) {
		/* Check that we will not end up in a loop... */
		if(stat(get_cstring(s),&st1)) {
			fprintf(stderr,"File Error.  Couldn't stat include file '%s' for input: ",get_cstring(s));
			perror(0);
			ABT_FUNC(AbMsg);
		}
		li=lex_info;
		while(li) {
			if(stat(li->name,&st2)) {
				fprintf(stderr,"File Error.  Couldn't stat include file '%s': ",li->name);
				perror(0);
				ABT_FUNC(AbMsg);
			}
			if(st1.st_ino==st2.st_ino && st1.st_dev==st2.st_dev) break;
			li=li->last;
		}
		if(li) {
			fprintf(stderr,"Recursive include files - '%s'\n",get_cstring(s));
			li=lex_info;
			while(li) {
				if(!stat(li->name,&st2) && (st1.st_ino==st2.st_ino && st1.st_dev==st2.st_dev))
				  fprintf(stderr,"         -> included from '%s' <-\n",li->name);
				else fprintf(stderr,"            included from '%s'\n",li->name);
				li=li->last;
			}
			ABT_FUNC(AbMsg);
		}
		li=lex_info->next;
		if(!li) {
			li=lk_malloc(sizeof(struct lex_info));
			li->next=0;
			li->last=lex_info;
			lex_info->next=li;
		}
		li->name=extract_cstring(s);
		li->ptr=li->ptr1=li->col_no=li->col_no1=li->line_no=li->eof=li->last_token=li->flag=0;
		if(!(li->fptr=fopen(li->name,"r"))) abt(__FILE__,__LINE__,"%s(): File Error.  Couldn't open '%s' for input\n",__func__,li->name);
		printf("-->Moving to include file '%s'\n",li->name);
		lex_info=li;
		handle->fname=lex_info->name;
		handle->line=&lex_info->line_no;
	}
}
  
static int read_buf(void) 
{
	int s;
	
	if(lex_info->eof) return 0;
	s=(int)fread(lex_info->buf+lex_info->ptr1,1,BUFSIZE-lex_info->ptr1,lex_info->fptr);
	lex_info->ptr1+=s;
	lex_info->ptr=0;
	if(feof(lex_info->fptr)) lex_info->eof=1;
	return s;
}

static int get_char(int flag) 
{
	int i;
	
	if(lex_info->flag) {
		lex_info->flag=0;
		while(lex_info) {
			fclose(lex_info->fptr); /* Yes, close file and move back */
			free(lex_info->name);
			lex_info=lex_info->last;
			printf("--> Moving back to include file '%s'\n",lex_info->name);
			handle->fname=lex_info->name;
			handle->line=&lex_info->line_no;
			if(lex_info->last_token) break;
		}
	}
	/* Check if we have enough characters in the buffer */
	while(((i=lex_info->ptr1-lex_info->ptr))<(flag?1:LOOKAHEAD)) {
		assert(i>=0);
		if(flag) return -1;
		if(!(lex_info->eof)) {
			/* Shift unread part of buffer to beginning */
			if(lex_info->ptr && i) {
				memmove(lex_info->buf,lex_info->buf+lex_info->ptr,i);
				lex_info->ptr1=i;
				lex_info->ptr=0;
			}
			/* Read in more characters (if possible) */
			if(read_buf()) break; /* OK, refilled buffer */
		}
		/* Can't read in more, are we at the end of the buffer ? */
		if(!i) {
			/* Yes, is there a parent file to return to? */
			if(lex_info->last) {
				i=-1;
				lex_info->flag=1;
				break;
			} else {
				lex_info->buf[lex_info->ptr1]=0; 
				break;
			}
		} else break;
	}
	if(i>=0) {
		/* Get next character */
		if(lex_info->ptr<lex_info->ptr1) {
			i=lex_info->buf[lex_info->ptr++];
			lex_info->col_no++;
		} else i=EOF;
		/* Handle DOS or Mac line ending sequences */
		if(i=='\r') {
			if(lex_info->ptr<lex_info->ptr1 && lex_info->buf[lex_info->ptr]=='\n') 
			  lex_info->ptr++;
			i='\n';
		}
		if(i=='\n') {
			lex_info->line_no++;
			lex_info->col_no1=lex_info->col_no;
			lex_info->col_no=0;
		}
	} else i=0;
	return i;
}

static struct bin_node *alloc_var(string *s)
{
	struct bin_node *node;
	struct parse_var *vv;
	
	node=lk_malloc(sizeof(struct bin_node));
	node->left=node->right=0;
	node->balance=0;
	vv=lk_malloc(sizeof(struct parse_var));
	node->data=vv;
	vv->name=s;
	vv->type=VARIABLE;
	vv->size=0;
	return node;
}

static struct bin_node *insert_var(string *s,struct bin_node *node,struct bin_node **node1,int *balanced)
{
	int i;
	struct parse_var *vv;
	
	vv=node->data;
	if((i=string_cmp(s,vv->name))) {
		if(i<0) {
			if(node->left) {
				node->left=insert_var(s,node->left,node1,balanced);
			} else {
				*node1=node->left=alloc_var(s);
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
		} else {
			if(node->right) {
				node->right=insert_var(s,node->right,node1,balanced);
			} else {
				*node1=node->right=alloc_var(s);
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
		}
	} else {
		*node1=node;
		*balanced=1;
		free_string(s);
	}
	return node;
}

struct parse_var *insert_prep_var(string *s)
{
	int i;
	struct bin_node *node;
	
	if(!root_var) node=root_var=alloc_var(s);
	else root_var=insert_var(s,root_var,&node,&i);
	return node->data;
}

static int check_token(string *s,yystype *lval,char c)
{
	int i,j;
	static char *builtin[]={"FOR","IF","ELSE","WHILE","DO","INCLUDE","ALIAS","AND","OR","NOT","EQ","NE",0};
	int builtin_tok[]={FOR,IF,ELSE,WHILE,DO,INCLUDE,ALIAS,'&','|','!',EQ_TOKEN,NE_TOKEN,0};
	struct bin_node *node;
	func_def *fl;
	
	i=-1;
	while(builtin[++i]) if(!strcasecmp(get_cstring(s),builtin[i])) break;
	j=builtin_tok[i];;
	if(!j && key_tab) {
		i=-1;
		while(key_tab[++i]) if(!strcasecmp(get_cstring(s),key_tab[i])) break;
		if(key_tab[i]) {
			j=KEYWORD;
			lval->i=i;
			free_string(s);
		}
	}
	if(!j && func_list && c=='(') {
		fl=func_list;
		i=0;
		while(fl) {
			if(!strcasecmp(get_cstring(s),fl->com)) break;
			fl=fl->next;
			i++;
		}
		if(fl) {
			j=FUNCTION;
			lval->f=fl;
			free_string(s);
		}
	}
	if(!j) {
		i=-1;
		if(command_tab) {
			while(command_tab[++i]) if(!strcasecmp(get_cstring(s),command_tab[i])) break;
		}
		if(!root_var) node=root_var=alloc_var(s);
		else root_var=insert_var(s,root_var,&node,&i);
		if(i>=0 && command_tab[i]) j=COMMAND;
		else j=TOKEN;
		lval->vv=node->data;
	}
	return j;
}

static void copy_lval(yystype *s1,yystype *s2,int tok)
{
	switch(tok) {
	 case INTEGER:
		s1->i=s2->i;
		break;
	 case REAL:
		s1->x=s2->x;
		break;
	 case TOKEN:
		s1->vv=s2->vv;
		break;
	 case STRING:
		s1->str=copy_string(s2->str);
		break;
	 default:
		ABT_FUNC("OOOK!\n");
		break;
	}
}

static int parlex(yystype *lval)
{
	string *s=0;
	int c,c1,c2,state,res,pt1=0;
	static int inc_state,postfix_flag,in_clause;
	struct lex_info *li;
	filepos fp;
	
	if(typeI_kludge) {
		typeI_kludge=0;
		return DUMMY;
	}
	for(;;) {
		c=get_char(0);
		if(c==EOF) {
			c=0;
			break;
		}
		if(isspace((int)c)) {
			postfix_flag=0;
			continue;
		}
		/* Remove comments */
		c1=lex_info->buf[lex_info->ptr];
		if(c=='#' || (c=='/' && c1=='/')) {
			do c=get_char(0); while(c && c!='\n');
		} else if(c=='/' && c1=='*') {
			do c=get_char(0); while (c!='*' || lex_info->buf[lex_info->ptr]!='/');
			(void)get_char(0);
		} else break;
	}
	res=0;
	set_filepos(&fp);
	/* Strings */
	if(c=='\"' || c=='\'') {
		for(;;) {
			c1=get_char(0);
			if(c1==EOF) c1=0;
			if(c1==c || !c1 || c1=='\n') break;
			if(c1=='\\') { /* Check for escape sequences */
				c1=get_char(0);
				switch(c1) {
				 case 'a':
					c1='\a';
					break;
				 case 'b':
					c1='\b';
					break;
				 case 'f':
					c1='\f';
					break;
				 case 'n':
					c1='\n';
					break;
				 case 'r':
					c1='\r';
					break;
				 case 't':
					c1='\t';
					break;
				 case 'v':
					c1='\v';
					break;
				}
			}
			s=add_to_string(s,c1);
		}
		if(c1!=c) parerror1("Unterminated string\n");
		c=res=STRING;
		lval->str=s;
	} else if(c=='[' && !in_clause) { /* Format Clauses */
		in_clause=1;
	} else if(c==']' && in_clause) {
		in_clause=0;
	}
	c1=lex_tab[c];
	if(!res) {
		if(in_clause) {
			if(in_clause==1 && (c1&FORMAT_BITS)) {
				if(c1&DIGIT_BIT) {
					pt1=lex_info->ptr-1;
					do {
						c=get_char(1);
						assert(c>=0);
						c1=lex_tab[c];
					} while(c1&DIGIT_BIT);
					res=CLAUSE_INT;
					lval->i=atoi(lex_info->buf+pt1);
					while(c && isspace(c)) c=get_char(0);
					c1=lex_tab[c];
					if(c1&(MATH_BIT|REL_BIT)) res=INTEGER;
				} else if(c1&X_BIT) {
					c2=lex_info->buf[lex_info->ptr];
					if(!(lex_tab[c2]&(TOKEN_BITS|REL_BIT|MATH_BIT|LBRACK_BIT))) {
						res=CLAUSE_X;
						c=get_char(0);
					}
				}
			} else if(c==';') in_clause=1;
		}
		if(!res && (c1&TOKEN_BITS)) {
			if(c1&(ALPHA_BIT|USCORE_BIT)) state=1;
			else if(c1&DIGIT_BIT) state=2;
			else if((c1&PERIOD_BIT)&&!postfix_flag) state=7;
			else state=0;
			/* Starting position of token */
			pt1=lex_info->ptr-1;
			while(state>0) {
				c=get_char(1);
				assert(c>=0);
				c1=lex_tab[c];
				switch(state) {
				 case 1: /* TOKEN */
					if(!(c1&(ALPHA_BIT|USCORE_BIT|DIGIT_BIT))) {
						res=TOKEN;
						state=0;
					}
					break;
				 case 2: /* INTEGER */
					if(c1&E_BIT) {
						res=INTEGER;
						state=4;
					} else if(c1&PERIOD_BIT) state=3;
					else if(!(c1&DIGIT_BIT)) {
						res=INTEGER;
						state=0;
					}
					break;
				 case 3: /* REAL */
					if(c1&E_BIT) {
						res=REAL;
						state=4;
					} else if(!(c1&DIGIT_BIT)) {
						res=REAL;
						state=0;
					}
					break;
				 case 4: /* REAL (exponential) */
					if(c1&DIGIT_BIT) state=5;
					else if(c1&PM_BIT) state=6;
					else {
						lex_info->ptr--;
						if(lex_info->col_no) lex_info->col_no--;
						else {
							lex_info->line_no--;
							lex_info->col_no=lex_info->col_no1;
						}
						state=0;
					}
					break;
				 case 5: /* REAL (exp. after digit) */
					if(!(c1&DIGIT_BIT)) {
						res=REAL;
						state=0;
					}
					break;
				 case 6: /* REAL (exp. after +-) */
					if(c1&DIGIT_BIT) {
						res=REAL;
						state=4;
					} else {
						lex_info->ptr-=2;
						state=0;
					}
					break;
				 case 7: /* After initial . */
					if(c1&DIGIT_BIT) {
						res=REAL;
						state=3;
					} else {
						state=0;
						res='.';
					}
				}
			}
		}
		if(in_clause==1 && res && res!=CLAUSE_INT && res!=CLAUSE_X) in_clause=2;
	}
	if(res) {
		copy_filepos(&handle->begin_tok,&fp);
	}
	if(res && res!=STRING) {
		lex_info->ptr--;
		if(lex_info->col_no) lex_info->col_no--;
		else {
			lex_info->line_no--;
			lex_info->col_no=lex_info->col_no1;
		}
		switch(res) {
		 case TOKEN:
			s=addn_to_string(s,lex_info->buf+pt1,lex_info->ptr-pt1);
			for(;;) {
				c1=lex_info->buf[lex_info->ptr];
				if(c1 && isspace(c1)) (void)get_char(0);
				else break;
			}
			res=check_token(s,lval,c1);
			if(res==TOKEN || res==KEYWORD || res==COMMAND) {
				if(c1=='(') {
					res=(res==KEYWORD?KEYWORD_B:WORD_B);
					lex_info->col_no++;
					lex_info->ptr++;
				}
			}
			break;
		 case INTEGER:
			lval->i=atoi(lex_info->buf+pt1);
			break;
		 case REAL:
			lval->x=atof(lex_info->buf+pt1);
			break;
		}
		c=res;
	} else if(c1&REL_BIT) {
		switch(c) {
		 case '<':
			c1=lex_info->buf[lex_info->ptr];
			if(c1=='=') {
				c=LEQ_TOKEN;
				lex_info->ptr++;
				lex_info->col_no++;
			} else if(c1=='>') {
				c=NE_TOKEN;
				lex_info->ptr++;
				lex_info->col_no++;
			}
			break;
		 case '>':
			if(lex_info->buf[lex_info->ptr]=='=') {
				c=GEQ_TOKEN;
				lex_info->ptr++;
				lex_info->col_no++;
			}
			break;
		 case '=':
			if(lex_info->buf[lex_info->ptr]=='=') {
				lex_info->ptr++;
				lex_info->col_no++;
				c=EQ_TOKEN;
			}
			break;
			/* Handle doubled symbols (no special meaning...) */
		 case '|':
			if(lex_info->buf[lex_info->ptr]=='|') {
				lex_info->ptr++;
				lex_info->col_no++;
			}
			break;
		 case '&':
			if(lex_info->buf[lex_info->ptr]=='&') {
				lex_info->ptr++;
				lex_info->col_no++;
			}
			break;
		 case '!':
			if(lex_info->buf[lex_info->ptr]=='=') {
				c=NE_TOKEN;
				lex_info->ptr++;
				lex_info->col_no++;
			}
		}
	} else if(c1&MATH_BIT) {
		c1=lex_info->buf[lex_info->ptr];
		if(c1==c) {
			if(c=='+' || c=='-') {
				lex_info->ptr++;
				lex_info->col_no++;
				lval->i=c;
				c=(postfix_flag?POSTFIX:PREFIX);
			} else if(c=='.') {
				c=DOUBLE_DOT;
				lex_info->ptr++;
				lex_info->col_no++;
			} else if(c=='*') {
				c=EXPO;
				lex_info->ptr++;
				lex_info->col_no++;
			}
		} else if(c1=='=') {
			lex_info->ptr++;
			lex_info->col_no++;
			lval->i=c;
			c=MATH_SHORT;
		}
	}
	if(!c && lex_info->flag) {
		li=lex_info->last;
		while(li) {
			c=li->last_token;
			if(c) {
				copy_lval(lval,&li->last_lval,c);
				break;
			}
			li=li->last;
		}
	}
	/* Kludge to avoid sending next token if we've finished an include statement 
	 * As this is a pure parse we can not use yyclearin to discard lookup token... */
	if(inc_state==1) {
		if(c==STRING || c==TOKEN || c==COMMAND) inc_state=2;
		else inc_state=3;
	} else if(inc_state==2) {
		if(c=='+') inc_state=1;
		else inc_state=3;
	} 
	if(inc_state==3) {
		if(!lex_info->flag) {
			/* Store last token, to be used when we return to this file */
			printf("storing token '%d' for file '%s'\n",c,lex_info->name);
			lex_info->last_token=c;
			copy_lval(&lex_info->last_lval,lval,c);
		} else lex_info->last_token=0;
		c=DUMMY;
		inc_state=0;
	}
	postfix_flag=(c==TOKEN || c==COMMAND || c==')')?1:0;
	return c;
}

static int parerror(char *s)
{
	int c;
	
	c=lex_info->col_no;
	if(!c) c++;
	if(s) fprintf(stderr,"%s, line %d column %d: %s\n",lex_info->name,lex_info->line_no+1,c,s);
	else fprintf(stderr,"%s, line %d column %d: ",lex_info->name,lex_info->line_no+1,c);
	return 0;
}

/* Report file position given by fp */
void parerror_fp(const filepos *fp)
{
	if(stdout) (void)fflush(stdout);
	if(fp->fname) fprintf(stderr,"%s, line %d column %d: ",fp->fname,fp->line+1,fp->col);
	else fputs("<Unknown position>: ",stderr);
}

/* Report file position at beginning of last token*/
void parerror1(const char *fmt, ...)
{
	va_list args;

	parerror_fp(&handle->begin_tok);
	va_start(args,fmt);
	(void)vfprintf(stderr,fmt,args);
	va_end(args);
}

/* Report file position at beginning of last command */
void parerror2(const char *fmt, ...)
{
	va_list args;

	parerror_fp(&handle->start_pos);
	va_start(args,fmt);
	(void)vfprintf(stderr,fmt,args);
	va_end(args);
}

void parerror3(filepos *fp,const char *fmt, ...)
{
	va_list args;

	parerror_fp(fp);
	va_start(args,fmt);
	(void)vfprintf(stderr,fmt,args);
	va_end(args);
}

/* Set up lexer table */
static void init_lex_tab(void) {
	int i,j;
	
	for(i=0;i<256;i++) {
		j=0;
		if(isalpha(i)) {
			j=ALPHA_BIT;
			if(i=='e' || i=='E') j|=E_BIT;
			else if(i=='x' || i=='X') j|=X_BIT;
		} else if(isdigit(i)) j=DIGIT_BIT;
		else switch(i) {
		 case '+':
		 case '-':
			j=MATH_BIT|PM_BIT;
			break;
		 case '_':
			j=USCORE_BIT;
			break;
		 case '.':
			j=MATH_BIT|PERIOD_BIT;
			break;
		 case ',':
			j=COMMA_BIT;
			break;
		 case '(':
			j=LBRACK_BIT;
			break;
		 case ')':
			j=RBRACK_BIT;
			break;
		 case '<':
		 case '>':
		 case '=':
		 case '!':
		 case '|':
		 case '&':
			j=REL_BIT;
			break;
		 case '/':
		 case '*':
			j=MATH_BIT;
			break;
		}
		lex_tab[i]=j;
	}
}

static void init_func_list(func_def *funcs)
{
	funcs=register_real_function(funcs,log,"log");
	funcs=register_real_function(funcs,log10,"log10");
	funcs=register_real_function(funcs,exp,"exp");
	funcs=register_real_function(funcs,sqrt,"sqrt");
	funcs=register_real_function(funcs,fabs,"abs");
	funcs=register_real_function(funcs,cbrt,"cbrt");
	funcs=register_real_function(funcs,sin,"sin");
	funcs=register_real_function(funcs,cos,"cos");
	funcs=register_real_function(funcs,tan,"tan");
	funcs=register_real_function(funcs,sinh,"sinh");
	funcs=register_real_function(funcs,cosh,"cosh");
	funcs=register_real_function(funcs,tanh,"tanh");
	funcs=register_real_function(funcs,asin,"asin");
	funcs=register_real_function(funcs,acos,"acos");
	funcs=register_real_function(funcs,atan,"atan");
	func_list=funcs;
}

static void free_vars(struct bin_node *node)
{
	struct parse_var *v;
	
	if(node->left) free_vars(node->left);
	if(node->right) free_vars(node->right);
	v=node->data;
	free_parse_var(v);
	free(node);
}

static void free_parse_stuff(void)
{
	if(root_var) free_vars(root_var);
	free(lex_info);
}

static void print_vars(struct bin_node *node)
{
	struct parse_var *v;
	
	if(node->left) print_vars(node->left);
	if(node->right) print_vars(node->right);
	v=node->data;
	fputs(v->name?get_cstring(v->name):"<NULL>",stdout);
	printf(" %x %x\n",v->type,v->ltype);
}

void print_all_vars(void)
{
	if(root_var) print_vars(root_var);
}

int Parse_File(FILE *fptr,char **keywords,char **commands,func_def *funcs,parse_handler *hand)
{
	int err;
	
#if YYDEBUG
	pardebug=1;
#endif
	
	set_strcmpfunc(strcasecmp);
	handle=hand;
	init_parse_utils(hand);
	lex_info=lk_calloc((size_t)1,sizeof(struct lex_info));
	lex_info->next=lex_info->last=0;
	lex_info->fptr=fptr;
	handle->line=&lex_info->line_no;
	lex_info->name=handle->fname;
	handle->start_pos.fname=handle->fname;
	handle->start_pos.col=1;
	handle->start_pos.line=1;
	init_lex_tab();
	key_tab=keywords;
	command_tab=commands;
	init_func_list(funcs);
	atexit(free_parse_stuff);
	err=parparse();
	return err;
}

