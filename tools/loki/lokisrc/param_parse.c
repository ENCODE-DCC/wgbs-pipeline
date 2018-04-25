
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
#define YYPURE 0

/* Push parsers.  */
#define YYPUSH 0

/* Pull parsers.  */
#define YYPULL 1

/* Using locations.  */
#define YYLSP_NEEDED 0



/* Copy the first part of user declarations.  */

/* Line 189 of yacc.c  */
#line 1 "param_parse.y"

/****************************************************************************
 *                                                                          *
 *     Loki - Programs for genetic analysis of complex traits using MCMC    *
 *                                                                          *
 *             Simon Heath - University of Washington                       *
 *                                                                          *
 *                       July 1997                                          *
 *                                                                          *
 * param_parse.y:                                                           *
 *                                                                          *
 * yacc source for parameter file parser.                                   *
 *                                                                          *
 ****************************************************************************/

#include <config.h>
#include <stdlib.h>
#include <string.h>
#ifdef USE_DMALLOC
#include <dmalloc.h>
#endif
#include <stdarg.h>
#include <stdio.h>
#include <math.h>
#include <sys/time.h>
	
#include "utils.h"
#include "loki.h"
#include "loki_scan.h"
#include "loki_ibd.h"
#include "loki_utils.h"
#include "lk_malloc.h"
#include "shared_peel.h"
#include "mat_utils.h"
#include "meiosis_scan.h"
#include "ranlib.h"
#include "snprintf.h"

static struct loki *loki;
static struct Marker *freq_marker;
static int check_variance(double,int);
static void print_scan_warn(char *, ...);
static void set_position(struct lk_variable *, double, double);
static struct lk_variable *find_var(char *, int, int);
static struct Marker *check_marker(struct lk_variable *);
static int find_allele(char *, struct Marker *);
static int find_group(char *,int, int);
static int find_trait(struct lk_variable *);
static void set_output_gen(char *,char *);
static struct IBD_List *add_ibd_list(double,struct IBD_List *);
static void set_freq(struct Marker *, double,int);
static void set_map_range(char *,double,double, int);
static void set_tloci(int,int);	
static void set_ibd_list(char *,struct IBD_List *,int);
static void set_ibd_markers(char *);
static void set_output(struct lk_variable *);
static void set_group_order(int);
static void set_analyze(char *);
static void set_ibd_mode(char *);
static void set_pseudo_chrome(struct string_list *,struct clause_atom *);
static struct num_array *make_num_array(int);
static void free_num_array(struct num_array *);
static void add_real_to_num_array(struct num_array *,double);
static void add_int_to_num_array(struct num_array *,int);
static void set_syst_var(int,struct num_array *);	
static int group_ptr,*group_order,group_counter,freq_allele,c_flag;
static struct clause_atom *new_clause_atom(char *,int);
static struct string_list *new_string_list(char *);

extern void yyerror(char *s),print_scan_err(char *fmt, ...);
extern int yyparse(void),yylex(void),lineno,lineno1,tokenpos;
extern char *yytext,linebuf[];

extern FILE *yyin;
int scan_error_n,iflag;
static int max_scan_errors=30,n_gen_grps,start_flag;
static int scan_warn_n,max_scan_warnings=30;


/* Line 189 of yacc.c  */
#line 153 "y.tab.c"

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
     RESIDUAL = 258,
     GENETIC = 259,
     VARIANCE = 260,
     POSITION = 261,
     FREQUENCY = 262,
     VIRTUAL = 263,
     START = 264,
     MEAN = 265,
     ITERATIONS = 266,
     SAMPLE = 267,
     FROM = 268,
     OUTPUT = 269,
     MAP = 270,
     TOTAL = 271,
     SEED = 272,
     SFILE = 273,
     SEEDFILE = 274,
     TRAIT = 275,
     LOCI = 276,
     SET = 277,
     SYSTEM_VAR = 278,
     TIMECOM = 279,
     ESTIMATE = 280,
     IBD = 281,
     GROUP = 282,
     ORDER = 283,
     MALE = 284,
     FEMALE = 285,
     LIMIT = 286,
     AFFECTED = 287,
     PHENO = 288,
     GENO = 289,
     COUNTS = 290,
     DUMP = 291,
     TYPE = 292,
     ANALYZE = 293,
     NORMAL = 294,
     STUDENT_T = 295,
     HAPLO = 296,
     INCLUDE = 297,
     FUNCTION = 298,
     HALDANE = 299,
     KOSAMBI = 300,
     POLYGENIC = 301,
     MARKERS = 302,
     GRID = 303,
     COMPRESS = 304,
     DIR = 305,
     PSEUDO = 306,
     STRING = 307,
     INTEGER = 308,
     REAL = 309
   };
#endif
/* Tokens.  */
#define RESIDUAL 258
#define GENETIC 259
#define VARIANCE 260
#define POSITION 261
#define FREQUENCY 262
#define VIRTUAL 263
#define START 264
#define MEAN 265
#define ITERATIONS 266
#define SAMPLE 267
#define FROM 268
#define OUTPUT 269
#define MAP 270
#define TOTAL 271
#define SEED 272
#define SFILE 273
#define SEEDFILE 274
#define TRAIT 275
#define LOCI 276
#define SET 277
#define SYSTEM_VAR 278
#define TIMECOM 279
#define ESTIMATE 280
#define IBD 281
#define GROUP 282
#define ORDER 283
#define MALE 284
#define FEMALE 285
#define LIMIT 286
#define AFFECTED 287
#define PHENO 288
#define GENO 289
#define COUNTS 290
#define DUMP 291
#define TYPE 292
#define ANALYZE 293
#define NORMAL 294
#define STUDENT_T 295
#define HAPLO 296
#define INCLUDE 297
#define FUNCTION 298
#define HALDANE 299
#define KOSAMBI 300
#define POLYGENIC 301
#define MARKERS 302
#define GRID 303
#define COMPRESS 304
#define DIR 305
#define PSEUDO 306
#define STRING 307
#define INTEGER 308
#define REAL 309




#if ! defined YYSTYPE && ! defined YYSTYPE_IS_DECLARED
typedef union YYSTYPE
{

/* Line 214 of yacc.c  */
#line 80 "param_parse.y"

	char *string;
	int value;
	double rvalue;
	struct IBD_List *rlist;
	struct lk_variable *lk_var;
	struct num_array *num_array;
	struct clause_atom *clause;
	struct string_list *string_list;



/* Line 214 of yacc.c  */
#line 310 "y.tab.c"
} YYSTYPE;
# define YYSTYPE_IS_TRIVIAL 1
# define yystype YYSTYPE /* obsolescent; will be withdrawn */
# define YYSTYPE_IS_DECLARED 1
#endif


/* Copy the second part of user declarations.  */


/* Line 264 of yacc.c  */
#line 322 "y.tab.c"

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
#define YYFINAL  4
/* YYLAST -- Last index in YYTABLE.  */
#define YYLAST   401

/* YYNTOKENS -- Number of terminals.  */
#define YYNTOKENS  63
/* YYNNTS -- Number of nonterminals.  */
#define YYNNTS  57
/* YYNRULES -- Number of rules.  */
#define YYNRULES  173
/* YYNRULES -- Number of states.  */
#define YYNSTATES  318

/* YYTRANSLATE(YYLEX) -- Bison symbol number corresponding to YYLEX.  */
#define YYUNDEFTOK  2
#define YYMAXUTOK   309

#define YYTRANSLATE(YYX)						\
  ((unsigned int) (YYX) <= YYMAXUTOK ? yytranslate[YYX] : YYUNDEFTOK)

/* YYTRANSLATE[YYLEX] -- Bison symbol number corresponding to YYLEX.  */
static const yytype_uint8 yytranslate[] =
{
       0,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
      60,    61,    62,     2,    55,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,    57,
       2,    56,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,    58,     2,    59,     2,     2,     2,     2,     2,     2,
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
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     1,     2,     3,     4,
       5,     6,     7,     8,     9,    10,    11,    12,    13,    14,
      15,    16,    17,    18,    19,    20,    21,    22,    23,    24,
      25,    26,    27,    28,    29,    30,    31,    32,    33,    34,
      35,    36,    37,    38,    39,    40,    41,    42,    43,    44,
      45,    46,    47,    48,    49,    50,    51,    52,    53,    54
};

#if YYDEBUG
/* YYPRHS[YYN] -- Index of the first RHS symbol of rule number YYN in
   YYRHS.  */
static const yytype_uint16 yyprhs[] =
{
       0,     0,     3,     4,     7,     9,    10,    14,    17,    19,
      21,    22,    26,    28,    30,    32,    34,    36,    38,    40,
      42,    44,    46,    48,    50,    52,    54,    56,    58,    60,
      62,    64,    66,    68,    70,    74,    80,    84,    90,    94,
      98,    99,   103,   104,   106,   108,   110,   114,   118,   122,
     126,   128,   132,   134,   137,   142,   148,   154,   158,   163,
     167,   174,   179,   183,   187,   191,   193,   195,   197,   199,
     203,   204,   208,   212,   216,   221,   226,   229,   234,   238,
     244,   247,   253,   260,   267,   271,   277,   282,   287,   291,
     297,   302,   307,   311,   315,   319,   324,   330,   335,   338,
     342,   348,   352,   356,   362,   366,   370,   376,   379,   383,
     387,   392,   397,   402,   406,   410,   414,   417,   421,   425,
     431,   435,   440,   445,   451,   456,   462,   466,   471,   476,
     482,   487,   493,   496,   500,   504,   510,   511,   516,   517,
     523,   524,   529,   531,   533,   538,   540,   544,   546,   547,
     552,   554,   556,   558,   560,   564,   566,   568,   570,   571,
     576,   577,   583,   585,   587,   591,   595,   597,   601,   603,
     605,   609,   613,   615
};

/* YYRHS -- A `-1'-separated list of the rules' RHS.  */
static const yytype_int8 yyrhs[] =
{
      64,     0,    -1,    -1,    65,    67,    -1,     1,    -1,    -1,
      64,    66,    67,    -1,    64,     1,    -1,    69,    -1,    70,
      -1,    -1,     9,    68,    69,    -1,    94,    -1,    96,    -1,
      99,    -1,   100,    -1,    98,    -1,    92,    -1,    71,    -1,
      93,    -1,    90,    -1,    89,    -1,    91,    -1,    72,    -1,
      81,    -1,   108,    -1,    95,    -1,    97,    -1,    86,    -1,
      73,    -1,    83,    -1,    88,    -1,    82,    -1,    80,    -1,
      12,    13,    53,    -1,    12,    13,    53,    75,    53,    -1,
       9,    14,    53,    -1,     9,    14,    53,    75,    53,    -1,
      22,    23,   118,    -1,    22,    52,   118,    -1,    -1,    42,
      74,    52,    -1,    -1,    55,    -1,    52,    -1,    76,    -1,
      77,    55,    76,    -1,    52,    56,    53,    -1,     9,    56,
      53,    -1,     7,    56,    53,    -1,    78,    -1,    79,    57,
      78,    -1,    51,    -1,    51,    77,    -1,    51,    58,    79,
      59,    -1,    51,    58,    79,    59,    77,    -1,    25,    26,
      52,    75,   117,    -1,    25,    26,   117,    -1,    25,    26,
      47,    52,    -1,    25,    26,    47,    -1,    25,    26,    48,
      52,    75,   117,    -1,    25,    26,    48,   117,    -1,    49,
      26,    14,    -1,    49,    14,    26,    -1,    25,    32,     7,
      -1,    52,    -1,    32,    -1,    26,    -1,    84,    -1,    85,
      55,    84,    -1,    -1,    38,    87,    85,    -1,    24,    31,
     119,    -1,    31,    24,   119,    -1,    31,     8,    24,   119,
      -1,     8,    24,    31,   119,    -1,    19,    52,    -1,    19,
      52,    55,    53,    -1,    17,    18,    52,    -1,    17,    18,
      52,    55,    53,    -1,    17,    53,    -1,    15,    52,   119,
      55,   119,    -1,    29,    15,    52,   119,    55,   119,    -1,
      30,    15,    52,   119,    55,   119,    -1,    16,    15,   119,
      -1,    16,    15,   119,    55,   119,    -1,    16,    29,    15,
     119,    -1,    16,    30,    15,   119,    -1,    15,    16,   119,
      -1,    15,    16,   119,    55,   119,    -1,    29,    15,    16,
     119,    -1,    30,    15,    16,   119,    -1,    15,    43,    44,
      -1,    15,    43,    45,    -1,    20,    21,    53,    -1,     9,
      20,    21,    53,    -1,    20,    21,    53,    55,    53,    -1,
      20,    21,    10,   119,    -1,    11,    53,    -1,    14,     7,
      53,    -1,    14,     7,    53,    55,    53,    -1,    14,     7,
      52,    -1,    12,     7,    53,    -1,    12,     7,    53,    55,
      53,    -1,    14,    33,    52,    -1,    14,    34,    52,    -1,
      14,    34,    52,    55,    52,    -1,    14,   106,    -1,    14,
      37,    53,    -1,    14,    18,    52,    -1,    14,     6,    18,
      52,    -1,    14,    26,    18,    52,    -1,    14,    26,    50,
      52,    -1,    36,    18,    52,    -1,    36,     7,    53,    -1,
      14,    41,    52,    -1,    14,    41,    -1,    14,    46,    52,
      -1,    14,    26,    52,    -1,    14,    26,    52,    55,    52,
      -1,     3,     5,   119,    -1,     3,     5,   104,   119,    -1,
       3,     5,    31,   119,    -1,     3,     5,    31,   104,   119,
      -1,    31,     3,     5,   119,    -1,    31,     3,     5,   104,
     119,    -1,     4,     5,   119,    -1,     4,     5,   104,   119,
      -1,     4,     5,    31,   119,    -1,     4,     5,    31,   104,
     119,    -1,    31,     4,     5,   119,    -1,    31,     4,     5,
     104,   119,    -1,    10,   119,    -1,    10,   104,   119,    -1,
       6,   105,   119,    -1,     6,   105,   119,    55,   119,    -1,
      -1,     7,   101,   107,   113,    -1,    -1,     7,   102,    35,
     107,   113,    -1,    -1,    35,   103,   107,   113,    -1,   105,
      -1,    52,    -1,    52,    60,    53,    61,    -1,   105,    -1,
     106,    55,   105,    -1,   105,    -1,    -1,    27,    28,   109,
     111,    -1,    53,    -1,    54,    -1,    52,    -1,   110,    -1,
     111,    55,   110,    -1,    53,    -1,    54,    -1,    52,    -1,
      -1,   112,    55,   114,   116,    -1,    -1,   113,   112,    55,
     115,   116,    -1,   119,    -1,    62,    -1,   116,    55,   119,
      -1,   116,    55,    62,    -1,   119,    -1,   117,    55,   119,
      -1,    54,    -1,    53,    -1,   118,    55,    54,    -1,   118,
      55,    53,    -1,    54,    -1,    53,    -1
};

/* YYRLINE[YYN] -- source line where rule number YYN was defined.  */
static const yytype_uint16 yyrline[] =
{
       0,   112,   112,   112,   113,   114,   114,   115,   118,   119,
     120,   120,   123,   124,   125,   126,   127,   130,   131,   132,
     133,   134,   135,   136,   137,   138,   139,   140,   141,   142,
     143,   144,   145,   146,   149,   150,   151,   152,   155,   156,
     159,   159,   162,   163,   166,   169,   170,   173,   174,   175,
     178,   179,   182,   183,   184,   185,   188,   189,   190,   191,
     192,   193,   196,   197,   200,   203,   204,   205,   208,   209,
     212,   212,   215,   216,   217,   218,   221,   222,   223,   224,
     225,   233,   234,   235,   236,   237,   238,   239,   240,   241,
     242,   243,   244,   245,   248,   249,   250,   251,   254,   256,
     257,   258,   259,   260,   261,   262,   263,   264,   265,   266,
     267,   268,   269,   270,   271,   272,   273,   274,   275,   276,
     279,   280,   283,   284,   285,   286,   289,   290,   293,   294,
     295,   296,   299,   302,   306,   307,   310,   310,   311,   311,
     312,   312,   315,   317,   318,   321,   322,   325,   327,   327,
     330,   331,   332,   335,   336,   339,   340,   341,   344,   344,
     345,   345,   348,   349,   350,   351,   354,   355,   358,   359,
     360,   361,   364,   365
};
#endif

#if YYDEBUG || YYERROR_VERBOSE || YYTOKEN_TABLE
/* YYTNAME[SYMBOL-NUM] -- String name of the symbol SYMBOL-NUM.
   First, the terminals, then, starting at YYNTOKENS, nonterminals.  */
static const char *const yytname[] =
{
  "$end", "error", "$undefined", "RESIDUAL", "GENETIC", "VARIANCE",
  "POSITION", "FREQUENCY", "VIRTUAL", "START", "MEAN", "ITERATIONS",
  "SAMPLE", "FROM", "OUTPUT", "MAP", "TOTAL", "SEED", "SFILE", "SEEDFILE",
  "TRAIT", "LOCI", "SET", "SYSTEM_VAR", "TIMECOM", "ESTIMATE", "IBD",
  "GROUP", "ORDER", "MALE", "FEMALE", "LIMIT", "AFFECTED", "PHENO", "GENO",
  "COUNTS", "DUMP", "TYPE", "ANALYZE", "NORMAL", "STUDENT_T", "HAPLO",
  "INCLUDE", "FUNCTION", "HALDANE", "KOSAMBI", "POLYGENIC", "MARKERS",
  "GRID", "COMPRESS", "DIR", "PSEUDO", "STRING", "INTEGER", "REAL", "','",
  "'='", "';'", "'['", "']'", "'('", "')'", "'*'", "$accept", "parmfile",
  "$@1", "$@2", "command1", "$@3", "command", "command_a", "samplecommand",
  "setcommand", "includecommand", "$@4", "opt_comma", "pchrom_list_atom",
  "pchrom_list", "clause_atom", "clause", "pseudocommand", "ibdcommand",
  "compresscommand", "aff_freqcommand", "analyzecom", "analyzelist",
  "analyzecommand", "$@5", "limit_timecommand", "seedcommand",
  "mapcommand", "tlocicommand", "itercommand", "outputcommand",
  "resvarcommand", "limitresvarcommand", "addvarcommand",
  "limitaddvarcommand", "meancommand", "positioncommand",
  "frequencycommand", "$@6", "$@7", "$@8", "trait_var", "lkvar",
  "lkvarlist", "lkmarker", "groupcommand", "$@9", "group", "grouplist",
  "allele", "freqlist", "$@10", "$@11", "freqlist1", "ibdlist", "array",
  "rnum", 0
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
     285,   286,   287,   288,   289,   290,   291,   292,   293,   294,
     295,   296,   297,   298,   299,   300,   301,   302,   303,   304,
     305,   306,   307,   308,   309,    44,    61,    59,    91,    93,
      40,    41,    42
};
# endif

/* YYR1[YYN] -- Symbol number of symbol that rule YYN derives.  */
static const yytype_uint8 yyr1[] =
{
       0,    63,    65,    64,    64,    66,    64,    64,    67,    67,
      68,    67,    69,    69,    69,    69,    69,    70,    70,    70,
      70,    70,    70,    70,    70,    70,    70,    70,    70,    70,
      70,    70,    70,    70,    71,    71,    71,    71,    72,    72,
      74,    73,    75,    75,    76,    77,    77,    78,    78,    78,
      79,    79,    80,    80,    80,    80,    81,    81,    81,    81,
      81,    81,    82,    82,    83,    84,    84,    84,    85,    85,
      87,    86,    88,    88,    88,    88,    89,    89,    89,    89,
      89,    90,    90,    90,    90,    90,    90,    90,    90,    90,
      90,    90,    90,    90,    91,    91,    91,    91,    92,    93,
      93,    93,    93,    93,    93,    93,    93,    93,    93,    93,
      93,    93,    93,    93,    93,    93,    93,    93,    93,    93,
      94,    94,    95,    95,    95,    95,    96,    96,    97,    97,
      97,    97,    98,    98,    99,    99,   101,   100,   102,   100,
     103,   100,   104,   105,   105,   106,   106,   107,   109,   108,
     110,   110,   110,   111,   111,   112,   112,   112,   114,   113,
     115,   113,   116,   116,   116,   116,   117,   117,   118,   118,
     118,   118,   119,   119
};

/* YYR2[YYN] -- Number of symbols composing right hand side of rule YYN.  */
static const yytype_uint8 yyr2[] =
{
       0,     2,     0,     2,     1,     0,     3,     2,     1,     1,
       0,     3,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     3,     5,     3,     5,     3,     3,
       0,     3,     0,     1,     1,     1,     3,     3,     3,     3,
       1,     3,     1,     2,     4,     5,     5,     3,     4,     3,
       6,     4,     3,     3,     3,     1,     1,     1,     1,     3,
       0,     3,     3,     3,     4,     4,     2,     4,     3,     5,
       2,     5,     6,     6,     3,     5,     4,     4,     3,     5,
       4,     4,     3,     3,     3,     4,     5,     4,     2,     3,
       5,     3,     3,     5,     3,     3,     5,     2,     3,     3,
       4,     4,     4,     3,     3,     3,     2,     3,     3,     5,
       3,     4,     4,     5,     4,     5,     3,     4,     4,     5,
       4,     5,     2,     3,     3,     5,     0,     4,     0,     5,
       0,     4,     1,     1,     4,     1,     3,     1,     0,     4,
       1,     1,     1,     1,     3,     1,     1,     1,     0,     4,
       0,     5,     1,     1,     3,     3,     1,     3,     1,     1,
       3,     3,     1,     1
};

/* YYDEFACT[STATE-NAME] -- Default rule to reduce with in state
   STATE-NUM when YYTABLE doesn't specify something else to do.  Zero
   means the default is an error.  */
static const yytype_uint8 yydefact[] =
{
       0,     4,     0,     0,     1,     7,     0,     0,     0,     0,
     136,     0,    10,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,   140,
       0,    70,    40,     0,    52,     3,     8,     9,    18,    23,
      29,    33,    24,    32,    30,    28,    31,    21,    20,    22,
      17,    19,    12,    26,    13,    27,    16,    14,    15,    25,
       6,     0,     0,   143,     0,     0,     0,     0,     0,     0,
       0,   173,   172,     0,   142,   132,    98,     0,     0,     0,
       0,     0,     0,     0,     0,     0,   116,     0,   145,   107,
       0,     0,     0,     0,     0,     0,     0,    80,    76,     0,
       0,     0,     0,     0,     0,   148,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,    44,
       0,    45,    53,     0,     0,   120,     0,     0,   126,     0,
     134,   147,     0,     0,     0,    36,     0,     0,     0,    11,
     133,   102,    34,     0,   101,    99,   109,     0,     0,   118,
     104,   105,   108,   115,   117,     0,    88,    92,    93,     0,
      84,     0,     0,    78,     0,     0,    94,   169,   168,    38,
      39,    72,    59,     0,    42,    57,   166,    64,     0,     0,
       0,     0,     0,     0,     0,     0,    73,     0,   114,   113,
      67,    66,    65,    68,    71,    41,    63,    62,     0,     0,
       0,    50,     0,     0,     0,   122,   121,     0,   128,   127,
       0,     0,   157,   155,   156,     0,   137,     0,    75,    43,
       0,    95,     0,     0,     0,     0,   110,     0,   111,   112,
       0,     0,   146,     0,     0,     0,    86,    87,     0,    77,
      97,     0,     0,    58,    42,    61,     0,     0,   152,   150,
     151,   153,   149,    90,     0,    91,     0,     0,   124,     0,
     130,    74,   141,     0,     0,     0,     0,     0,    54,    46,
     123,   129,   144,   135,   158,     0,   139,    37,   103,    35,
     100,   119,   106,    89,    81,    85,    79,    96,   171,   170,
       0,    56,   167,     0,     0,     0,   125,   131,    69,    49,
      48,    47,    51,    55,     0,   160,    60,   154,    82,    83,
     163,   159,   162,     0,     0,   161,   165,   164
};

/* YYDEFGOTO[NTERM-NUM].  */
static const yytype_int16 yydefgoto[] =
{
      -1,     2,     3,     6,    35,    70,    36,    37,    38,    39,
      40,   116,   220,   121,   122,   201,   202,    41,    42,    43,
      44,   193,   194,    45,   115,    46,    47,    48,    49,    50,
      51,    52,    53,    54,    55,    56,    57,    58,    65,    66,
     112,   124,    74,    89,   132,    59,   178,   251,   252,   215,
     216,   304,   313,   311,   175,   169,   176
};

/* YYPACT[STATE-NUM] -- Index in YYTABLE of the portion describing
   STATE-NUM.  */
#define YYPACT_NINF -213
static const yytype_int16 yypact[] =
{
     301,  -213,   248,   350,  -213,  -213,   350,    17,    25,   -38,
      10,    35,    49,    53,    23,    59,     1,     9,    52,   -12,
      21,    41,     6,    55,    94,    73,    68,   110,    20,  -213,
       8,  -213,  -213,    64,    72,  -213,  -213,  -213,  -213,  -213,
    -213,  -213,  -213,  -213,  -213,  -213,  -213,  -213,  -213,  -213,
    -213,  -213,  -213,  -213,  -213,  -213,  -213,  -213,  -213,  -213,
    -213,    40,    43,    57,    78,   -38,    74,   113,   105,   138,
      81,  -213,  -213,    78,  -213,  -213,  -213,   107,   108,   144,
      98,   111,    18,   112,   121,   124,   122,   126,  -213,   125,
      78,   102,    78,    78,   150,   164,   129,  -213,   127,    -7,
     100,   100,    78,    75,   176,  -213,     4,     5,   179,   180,
     162,    78,   -38,   134,   136,   -15,   137,   166,   181,  -213,
       3,  -213,   135,    53,    78,  -213,    53,    78,  -213,   140,
     141,  -213,    82,   -38,    78,    90,   146,   192,   195,  -213,
    -213,   147,    90,   149,  -213,   148,  -213,   152,   153,   151,
    -213,   156,  -213,  -213,  -213,   -38,   157,  -213,  -213,   158,
     159,    78,    78,   161,   154,    78,   163,  -213,  -213,   168,
     168,  -213,   165,    85,   169,   170,  -213,  -213,    88,    78,
      78,    78,    78,    53,    53,    78,  -213,    82,  -213,  -213,
    -213,  -213,  -213,  -213,   171,  -213,  -213,  -213,   172,   173,
     174,  -213,    45,   167,    78,  -213,  -213,    78,  -213,  -213,
     175,    78,  -213,  -213,  -213,   177,    82,    82,  -213,  -213,
     155,  -213,    53,    53,   178,   182,  -213,   184,  -213,  -213,
     186,   187,  -213,    78,    78,    78,  -213,  -213,   188,  -213,
    -213,   189,   103,  -213,   169,   170,    78,    78,  -213,  -213,
    -213,  -213,   185,  -213,   190,  -213,   198,    78,  -213,    78,
    -213,  -213,    82,   -15,   194,   197,   208,     3,   167,  -213,
    -213,  -213,  -213,  -213,  -213,   211,    82,  -213,  -213,  -213,
    -213,  -213,  -213,  -213,  -213,  -213,  -213,  -213,  -213,  -213,
      78,   170,  -213,    88,    78,    78,  -213,  -213,  -213,  -213,
    -213,  -213,  -213,   135,   -22,  -213,   170,  -213,  -213,  -213,
    -213,   214,  -213,   -22,    46,   214,  -213,  -213
};

/* YYPGOTO[NTERM-NUM].  */
static const yytype_int16 yypgoto[] =
{
    -213,  -213,  -213,  -213,   221,  -213,   201,  -213,  -213,  -213,
    -213,  -213,  -141,    30,   -25,     7,  -213,  -213,  -213,  -213,
    -213,    13,  -213,  -213,  -213,  -213,  -213,  -213,  -213,  -213,
    -213,  -213,  -213,  -213,  -213,  -213,  -213,  -213,  -213,  -213,
    -213,    -8,     0,  -213,   -94,  -213,  -213,    -6,  -213,  -212,
    -174,  -213,  -213,   -33,  -171,   191,   -13
};

/* YYTABLE[YYPACT[STATE-NUM]].  What to do in state STATE-NUM.  If
   positive, shift that token.  If negative, reduce the rule which
   number is the opposite.  If zero, do what YYDEFACT says.
   If YYTABLE_NINF, syntax error.  */
#define YYTABLE_NINF -139
static const yytype_int16 yytable[] =
{
      75,   225,   245,   165,   275,    73,    96,    79,    80,    64,
     198,   190,   199,   262,    63,   113,    88,   191,   187,    81,
     179,   181,    61,   108,   109,    90,   114,    82,   110,   100,
      62,    71,    72,   246,    83,    84,   147,   192,    85,   217,
     310,    97,    86,   276,   111,  -138,   166,    87,   125,   128,
     275,   130,    91,    63,   127,   200,   180,   182,   101,    67,
     140,    92,    99,    68,   275,   131,    77,    93,   148,    69,
     149,   123,    78,    98,   126,   291,    76,   156,   117,   159,
     160,    94,    95,   106,   137,   138,   102,     9,    10,   171,
     118,    13,    63,    71,    72,    63,    71,    72,   186,    71,
      72,   105,   267,   290,   268,    63,    71,    72,   316,   133,
     205,   206,   131,   208,   209,   204,    29,   129,   207,   306,
     103,   218,   172,   173,   119,   107,   104,   174,    71,    72,
     120,    71,    72,   131,   212,   213,   214,   244,    71,    72,
     248,   249,   250,   -42,   134,   219,   157,   158,   236,   237,
     144,   145,   240,   167,   168,   232,   288,   289,   135,   136,
     141,   142,   143,   146,   150,   161,   253,   254,   255,   256,
     258,   260,   261,   151,   153,   257,   259,   152,   154,   162,
     155,   163,   164,   177,   183,   184,   185,   188,   189,   195,
     203,   270,   196,   210,   271,   197,   211,   222,   273,   221,
     223,   226,   224,   227,   228,   229,   230,   239,   277,   125,
     128,   231,   233,   234,   235,   127,   238,   243,   241,   119,
     283,   284,   285,   242,   219,   247,   263,    60,   264,   265,
     266,   278,   274,   269,   292,   279,   272,   280,   281,   282,
     293,   286,   287,   303,   296,   294,   297,   299,     4,     5,
     300,    -5,    -5,   295,    -5,    -5,    -5,    -5,    -5,    -5,
      -5,   301,    -5,    -5,    -5,    -5,   305,    -5,    -5,   314,
      -5,   139,    -5,    -5,   302,    -5,   298,    -5,    -5,    -5,
     315,   308,   309,    -5,    -5,     0,    -5,   307,     0,     0,
      -5,   312,   170,     0,     0,     0,     0,    -5,     0,    -5,
     312,   317,     1,     0,    -2,    -2,     0,    -2,    -2,    -2,
      -2,    -2,    -2,    -2,     0,    -2,    -2,    -2,    -2,     0,
      -2,    -2,     0,    -2,     0,    -2,    -2,     0,    -2,     0,
      -2,    -2,    -2,     0,     0,     0,    -2,    -2,     0,    -2,
       0,     0,     0,    -2,     0,     0,     0,     0,     0,     0,
      -2,     0,    -2,     7,     8,     0,     9,    10,    11,    12,
      13,    14,    15,     0,    16,    17,    18,    19,     0,    20,
      21,     0,    22,     0,    23,    24,     0,    25,     0,    26,
      27,    28,     0,     0,     0,    29,    30,     0,    31,     0,
       0,     0,    32,     0,     0,     0,     0,     0,     0,    33,
       0,    34
};

static const yytype_int16 yycheck[] =
{
      13,   142,   173,    10,   216,    13,    18,     6,     7,     9,
       7,    26,     9,   187,    52,     7,    16,    32,   112,    18,
      16,    16,     5,     3,     4,    16,    18,    26,     8,    23,
       5,    53,    54,   174,    33,    34,    18,    52,    37,   133,
      62,    53,    41,   217,    24,    35,    53,    46,    61,    62,
     262,    64,    43,    52,    62,    52,    52,    52,    52,    24,
      73,    52,    21,    14,   276,    65,     7,    15,    50,    20,
      52,    31,    13,    52,    31,   246,    53,    90,    14,    92,
      93,    29,    30,    15,     3,     4,    31,     6,     7,   102,
      26,    10,    52,    53,    54,    52,    53,    54,   111,    53,
      54,    28,    57,   244,    59,    52,    53,    54,    62,    35,
     123,   124,   112,   126,   127,   123,    35,    60,   126,   290,
      26,   134,    47,    48,    52,    15,    32,    52,    53,    54,
      58,    53,    54,   133,    52,    53,    54,    52,    53,    54,
      52,    53,    54,    53,    31,    55,    44,    45,   161,   162,
      52,    53,   165,    53,    54,   155,    53,    54,    53,    21,
      53,    53,    18,    52,    52,    15,   179,   180,   181,   182,
     183,   184,   185,    52,    52,   183,   184,    53,    52,    15,
      55,    52,    55,     7,     5,     5,    24,    53,    52,    52,
      55,   204,    26,    53,   207,    14,    55,     5,   211,    53,
       5,    52,    55,    55,    52,    52,    55,    53,    53,   222,
     223,    55,    55,    55,    55,   223,    55,    52,    55,    52,
     233,   234,   235,    55,    55,    55,    55,     6,    56,    56,
      56,    53,    55,   203,   247,    53,    61,    53,    52,    52,
      55,    53,    53,   268,   257,    55,   259,    53,     0,     1,
      53,     3,     4,    55,     6,     7,     8,     9,    10,    11,
      12,    53,    14,    15,    16,    17,    55,    19,    20,    55,
      22,    70,    24,    25,   267,    27,   263,    29,    30,    31,
     313,   294,   295,    35,    36,    -1,    38,   293,    -1,    -1,
      42,   304,   101,    -1,    -1,    -1,    -1,    49,    -1,    51,
     313,   314,     1,    -1,     3,     4,    -1,     6,     7,     8,
       9,    10,    11,    12,    -1,    14,    15,    16,    17,    -1,
      19,    20,    -1,    22,    -1,    24,    25,    -1,    27,    -1,
      29,    30,    31,    -1,    -1,    -1,    35,    36,    -1,    38,
      -1,    -1,    -1,    42,    -1,    -1,    -1,    -1,    -1,    -1,
      49,    -1,    51,     3,     4,    -1,     6,     7,     8,     9,
      10,    11,    12,    -1,    14,    15,    16,    17,    -1,    19,
      20,    -1,    22,    -1,    24,    25,    -1,    27,    -1,    29,
      30,    31,    -1,    -1,    -1,    35,    36,    -1,    38,    -1,
      -1,    -1,    42,    -1,    -1,    -1,    -1,    -1,    -1,    49,
      -1,    51
};

/* YYSTOS[STATE-NUM] -- The (internal number of the) accessing
   symbol of state STATE-NUM.  */
static const yytype_uint8 yystos[] =
{
       0,     1,    64,    65,     0,     1,    66,     3,     4,     6,
       7,     8,     9,    10,    11,    12,    14,    15,    16,    17,
      19,    20,    22,    24,    25,    27,    29,    30,    31,    35,
      36,    38,    42,    49,    51,    67,    69,    70,    71,    72,
      73,    80,    81,    82,    83,    86,    88,    89,    90,    91,
      92,    93,    94,    95,    96,    97,    98,    99,   100,   108,
      67,     5,     5,    52,   105,   101,   102,    24,    14,    20,
      68,    53,    54,   104,   105,   119,    53,     7,    13,     6,
       7,    18,    26,    33,    34,    37,    41,    46,   105,   106,
      16,    43,    52,    15,    29,    30,    18,    53,    52,    21,
      23,    52,    31,    26,    32,    28,    15,    15,     3,     4,
       8,    24,   103,     7,    18,    87,    74,    14,    26,    52,
      58,    76,    77,    31,   104,   119,    31,   104,   119,    60,
     119,   105,   107,    35,    31,    53,    21,     3,     4,    69,
     119,    53,    53,    18,    52,    53,    52,    18,    50,    52,
      52,    52,    53,    52,    52,    55,   119,    44,    45,   119,
     119,    15,    15,    52,    55,    10,    53,    53,    54,   118,
     118,   119,    47,    48,    52,   117,   119,     7,   109,    16,
      52,    16,    52,     5,     5,    24,   119,   107,    53,    52,
      26,    32,    52,    84,    85,    52,    26,    14,     7,     9,
      52,    78,    79,    55,   104,   119,   119,   104,   119,   119,
      53,    55,    52,    53,    54,   112,   113,   107,   119,    55,
      75,    53,     5,     5,    55,    75,    52,    55,    52,    52,
      55,    55,   105,    55,    55,    55,   119,   119,    55,    53,
     119,    55,    55,    52,    52,   117,    75,    55,    52,    53,
      54,   110,   111,   119,   119,   119,   119,   104,   119,   104,
     119,   119,   113,    55,    56,    56,    56,    57,    59,    76,
     119,   119,    61,   119,    55,   112,   113,    53,    53,    53,
      53,    52,    52,   119,   119,   119,    53,    53,    53,    54,
      75,   117,   119,    55,    55,    55,   119,   119,    84,    53,
      53,    53,    78,    77,   114,    55,   117,   110,   119,   119,
      62,   116,   119,   115,    55,   116,    62,   119
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
# define YYLEX yylex (YYLEX_PARAM)
#else
# define YYLEX yylex ()
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


/* The lookahead symbol.  */
int yychar;

/* The semantic value of the lookahead symbol.  */
YYSTYPE yylval;

/* Number of syntax errors so far.  */
int yynerrs;



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
        case 2:

/* Line 1455 of yacc.c  */
#line 112 "param_parse.y"
    {lineno1=lineno;}
    break;

  case 3:

/* Line 1455 of yacc.c  */
#line 112 "param_parse.y"
    {iflag=0;}
    break;

  case 4:

/* Line 1455 of yacc.c  */
#line 113 "param_parse.y"
    {iflag=0;}
    break;

  case 5:

/* Line 1455 of yacc.c  */
#line 114 "param_parse.y"
    {lineno1=lineno;}
    break;

  case 6:

/* Line 1455 of yacc.c  */
#line 114 "param_parse.y"
    {iflag=0;}
    break;

  case 7:

/* Line 1455 of yacc.c  */
#line 115 "param_parse.y"
    {iflag=0;}
    break;

  case 10:

/* Line 1455 of yacc.c  */
#line 120 "param_parse.y"
    {start_flag=2;}
    break;

  case 11:

/* Line 1455 of yacc.c  */
#line 120 "param_parse.y"
    {start_flag=1;}
    break;

  case 34:

/* Line 1455 of yacc.c  */
#line 149 "param_parse.y"
    {loki->params.sample_from[0]=loki->params.sample_from[1]=(yyvsp[(3) - (3)].value);}
    break;

  case 35:

/* Line 1455 of yacc.c  */
#line 150 "param_parse.y"
    {loki->params.sample_from[0]=(yyvsp[(3) - (5)].value); loki->params.sample_from[1]=(yyvsp[(5) - (5)].value);}
    break;

  case 36:

/* Line 1455 of yacc.c  */
#line 151 "param_parse.y"
    {loki->params.sample_from[0]=loki->params.sample_from[1]=(yyvsp[(3) - (3)].value);}
    break;

  case 37:

/* Line 1455 of yacc.c  */
#line 152 "param_parse.y"
    {loki->params.sample_from[0]=(yyvsp[(3) - (5)].value); loki->params.sample_from[1]=(yyvsp[(5) - (5)].value);}
    break;

  case 38:

/* Line 1455 of yacc.c  */
#line 155 "param_parse.y"
    {set_syst_var((yyvsp[(2) - (3)].value),(yyvsp[(3) - (3)].num_array)); free_num_array((yyvsp[(3) - (3)].num_array));}
    break;

  case 39:

/* Line 1455 of yacc.c  */
#line 156 "param_parse.y"
    { free((yyvsp[(2) - (3)].string)); free_num_array((yyvsp[(3) - (3)].num_array)); yyerror("Unrecognized system variable"); }
    break;

  case 40:

/* Line 1455 of yacc.c  */
#line 159 "param_parse.y"
    {iflag=1;}
    break;

  case 41:

/* Line 1455 of yacc.c  */
#line 159 "param_parse.y"
    {include_param_file((yyvsp[(3) - (3)].string));}
    break;

  case 44:

/* Line 1455 of yacc.c  */
#line 166 "param_parse.y"
    {(yyval.string_list)=new_string_list((yyvsp[(1) - (1)].string));}
    break;

  case 46:

/* Line 1455 of yacc.c  */
#line 170 "param_parse.y"
    { (yyvsp[(3) - (3)].string_list)->next=(yyvsp[(1) - (3)].string_list); (yyval.string_list)=(yyvsp[(3) - (3)].string_list); }
    break;

  case 47:

/* Line 1455 of yacc.c  */
#line 173 "param_parse.y"
    {(yyval.clause)=new_clause_atom((yyvsp[(1) - (3)].string),(yyvsp[(3) - (3)].value));}
    break;

  case 48:

/* Line 1455 of yacc.c  */
#line 174 "param_parse.y"
    {(yyval.clause)=new_clause_atom(strdup("INIT"),(yyvsp[(3) - (3)].value));}
    break;

  case 49:

/* Line 1455 of yacc.c  */
#line 175 "param_parse.y"
    {(yyval.clause)=new_clause_atom(strdup("FREQ"),(yyvsp[(3) - (3)].value));}
    break;

  case 51:

/* Line 1455 of yacc.c  */
#line 179 "param_parse.y"
    { (yyvsp[(3) - (3)].clause)->next=(yyvsp[(1) - (3)].clause); (yyval.clause)=(yyvsp[(3) - (3)].clause); }
    break;

  case 52:

/* Line 1455 of yacc.c  */
#line 182 "param_parse.y"
    {set_pseudo_chrome(0,0);}
    break;

  case 53:

/* Line 1455 of yacc.c  */
#line 183 "param_parse.y"
    {set_pseudo_chrome((yyvsp[(2) - (2)].string_list),0);}
    break;

  case 54:

/* Line 1455 of yacc.c  */
#line 184 "param_parse.y"
    {set_pseudo_chrome(0,(yyvsp[(3) - (4)].clause));}
    break;

  case 55:

/* Line 1455 of yacc.c  */
#line 185 "param_parse.y"
    {set_pseudo_chrome((yyvsp[(5) - (5)].string_list),(yyvsp[(3) - (5)].clause));}
    break;

  case 56:

/* Line 1455 of yacc.c  */
#line 188 "param_parse.y"
    { set_ibd_list((yyvsp[(3) - (5)].string),(yyvsp[(5) - (5)].rlist),IBD_EST_DISCRETE); free((yyvsp[(3) - (5)].string));}
    break;

  case 57:

/* Line 1455 of yacc.c  */
#line 189 "param_parse.y"
    { set_ibd_list(0,(yyvsp[(3) - (3)].rlist),IBD_EST_DISCRETE); }
    break;

  case 58:

/* Line 1455 of yacc.c  */
#line 190 "param_parse.y"
    { set_ibd_markers((yyvsp[(4) - (4)].string)); free((yyvsp[(4) - (4)].string));}
    break;

  case 59:

/* Line 1455 of yacc.c  */
#line 191 "param_parse.y"
    { set_ibd_markers(0); }
    break;

  case 60:

/* Line 1455 of yacc.c  */
#line 192 "param_parse.y"
    {set_ibd_list((yyvsp[(4) - (6)].string),(yyvsp[(6) - (6)].rlist),IBD_EST_GRID); free((yyvsp[(4) - (6)].string));}
    break;

  case 61:

/* Line 1455 of yacc.c  */
#line 193 "param_parse.y"
    { set_ibd_list(0,(yyvsp[(4) - (4)].rlist),IBD_EST_GRID); }
    break;

  case 62:

/* Line 1455 of yacc.c  */
#line 196 "param_parse.y"
    {loki->params.compress_ibd=1;}
    break;

  case 63:

/* Line 1455 of yacc.c  */
#line 197 "param_parse.y"
    {loki->params.compress_ibd=1;}
    break;

  case 64:

/* Line 1455 of yacc.c  */
#line 200 "param_parse.y"
    {loki->params.est_aff_freq=1;}
    break;

  case 65:

/* Line 1455 of yacc.c  */
#line 203 "param_parse.y"
    {set_analyze((yyvsp[(1) - (1)].string)); free((yyvsp[(1) - (1)].string));}
    break;

  case 66:

/* Line 1455 of yacc.c  */
#line 204 "param_parse.y"
    {set_analyze("AFFECTED");}
    break;

  case 67:

/* Line 1455 of yacc.c  */
#line 205 "param_parse.y"
    {set_analyze("IBD");}
    break;

  case 70:

/* Line 1455 of yacc.c  */
#line 212 "param_parse.y"
    {loki->params.analysis=0;}
    break;

  case 72:

/* Line 1455 of yacc.c  */
#line 215 "param_parse.y"
    {loki->params.limit_time=(yyvsp[(3) - (3)].rvalue),loki->params.limit_timer_type=ITIMER_REAL;}
    break;

  case 73:

/* Line 1455 of yacc.c  */
#line 216 "param_parse.y"
    {loki->params.limit_time=(yyvsp[(3) - (3)].rvalue),loki->params.limit_timer_type=ITIMER_REAL;}
    break;

  case 74:

/* Line 1455 of yacc.c  */
#line 217 "param_parse.y"
    {loki->params.limit_time=(yyvsp[(4) - (4)].rvalue),loki->params.limit_timer_type=ITIMER_VIRTUAL;}
    break;

  case 75:

/* Line 1455 of yacc.c  */
#line 218 "param_parse.y"
    {loki->params.limit_time=(yyvsp[(4) - (4)].rvalue),loki->params.limit_timer_type=ITIMER_VIRTUAL;}
    break;

  case 76:

/* Line 1455 of yacc.c  */
#line 221 "param_parse.y"
    {if(loki->names[LK_SEEDFILE]) free(loki->names[LK_SEEDFILE]); loki->names[LK_SEEDFILE]=(yyvsp[(2) - (2)].string);}
    break;

  case 77:

/* Line 1455 of yacc.c  */
#line 222 "param_parse.y"
    {if(loki->names[LK_SEEDFILE]) free(loki->names[LK_SEEDFILE]); loki->names[LK_SEEDFILE]=(yyvsp[(2) - (4)].string); if((yyvsp[(4) - (4)].value)) loki->params.ranseed_set|=1; else loki->params.ranseed_set&=~1;}
    break;

  case 78:

/* Line 1455 of yacc.c  */
#line 223 "param_parse.y"
    {if(loki->names[LK_SEEDFILE]) free(loki->names[LK_SEEDFILE]); loki->names[LK_SEEDFILE]=(yyvsp[(3) - (3)].string);}
    break;

  case 79:

/* Line 1455 of yacc.c  */
#line 224 "param_parse.y"
    {if(loki->names[LK_SEEDFILE]) free(loki->names[LK_SEEDFILE]); loki->names[LK_SEEDFILE]=(yyvsp[(3) - (5)].string); if((yyvsp[(5) - (5)].value)) loki->params.ranseed_set|=1; else loki->params.ranseed_set&=~1;}
    break;

  case 80:

/* Line 1455 of yacc.c  */
#line 225 "param_parse.y"
    {
			 if((yyvsp[(2) - (2)].value)<0) yyerror("Seedvalue out of range");
			 else {
				 init_ranf((yyvsp[(2) - (2)].value));
				 loki->params.ranseed_set|=2;
			 }
		 }
    break;

  case 81:

/* Line 1455 of yacc.c  */
#line 233 "param_parse.y"
    {set_map_range((yyvsp[(2) - (5)].string),(yyvsp[(3) - (5)].rvalue),(yyvsp[(5) - (5)].rvalue),-1); free((yyvsp[(2) - (5)].string));}
    break;

  case 82:

/* Line 1455 of yacc.c  */
#line 234 "param_parse.y"
    {set_map_range((yyvsp[(3) - (6)].string),(yyvsp[(4) - (6)].rvalue),(yyvsp[(6) - (6)].rvalue),X_PAT); free((yyvsp[(3) - (6)].string)); }
    break;

  case 83:

/* Line 1455 of yacc.c  */
#line 235 "param_parse.y"
    {set_map_range((yyvsp[(3) - (6)].string),(yyvsp[(4) - (6)].rvalue),(yyvsp[(6) - (6)].rvalue),X_MAT); free((yyvsp[(3) - (6)].string)); }
    break;

  case 84:

/* Line 1455 of yacc.c  */
#line 236 "param_parse.y"
    {set_map_range(0,(yyvsp[(3) - (3)].rvalue),(yyvsp[(3) - (3)].rvalue),-1);}
    break;

  case 85:

/* Line 1455 of yacc.c  */
#line 237 "param_parse.y"
    {set_map_range(0,(yyvsp[(3) - (5)].rvalue),(yyvsp[(5) - (5)].rvalue),-2);}
    break;

  case 86:

/* Line 1455 of yacc.c  */
#line 238 "param_parse.y"
    {set_map_range(0,(yyvsp[(4) - (4)].rvalue),(yyvsp[(4) - (4)].rvalue),X_PAT);}
    break;

  case 87:

/* Line 1455 of yacc.c  */
#line 239 "param_parse.y"
    {set_map_range(0,(yyvsp[(4) - (4)].rvalue),(yyvsp[(4) - (4)].rvalue),X_MAT);}
    break;

  case 88:

/* Line 1455 of yacc.c  */
#line 240 "param_parse.y"
    {set_map_range(0,(yyvsp[(3) - (3)].rvalue),(yyvsp[(3) - (3)].rvalue),-1);}
    break;

  case 89:

/* Line 1455 of yacc.c  */
#line 241 "param_parse.y"
    {set_map_range(0,(yyvsp[(3) - (5)].rvalue),(yyvsp[(5) - (5)].rvalue),-2);}
    break;

  case 90:

/* Line 1455 of yacc.c  */
#line 242 "param_parse.y"
    {set_map_range(0,(yyvsp[(4) - (4)].rvalue),(yyvsp[(4) - (4)].rvalue),X_PAT);}
    break;

  case 91:

/* Line 1455 of yacc.c  */
#line 243 "param_parse.y"
    {set_map_range(0,(yyvsp[(4) - (4)].rvalue),(yyvsp[(4) - (4)].rvalue),X_MAT);}
    break;

  case 92:

/* Line 1455 of yacc.c  */
#line 244 "param_parse.y"
    {loki->params.map_function=MAP_HALDANE;}
    break;

  case 93:

/* Line 1455 of yacc.c  */
#line 245 "param_parse.y"
    {loki->params.map_function=MAP_KOSAMBI;}
    break;

  case 94:

/* Line 1455 of yacc.c  */
#line 248 "param_parse.y"
    {set_tloci(-1,(yyvsp[(3) - (3)].value));}
    break;

  case 95:

/* Line 1455 of yacc.c  */
#line 249 "param_parse.y"
    {set_tloci(-2,(yyvsp[(4) - (4)].value));}
    break;

  case 96:

/* Line 1455 of yacc.c  */
#line 250 "param_parse.y"
    {set_tloci((yyvsp[(3) - (5)].value),(yyvsp[(5) - (5)].value));}
    break;

  case 97:

/* Line 1455 of yacc.c  */
#line 251 "param_parse.y"
    {loki->models->tloci_mean=(yyvsp[(4) - (4)].rvalue); loki->models->tloci_mean_set=1;}
    break;

  case 98:

/* Line 1455 of yacc.c  */
#line 254 "param_parse.y"
    { loki->params.num_iter=(yyvsp[(2) - (2)].value); }
    break;

  case 99:

/* Line 1455 of yacc.c  */
#line 256 "param_parse.y"
    {loki->params.sample_freq[0]=loki->params.sample_freq[1]=(yyvsp[(3) - (3)].value);}
    break;

  case 100:

/* Line 1455 of yacc.c  */
#line 257 "param_parse.y"
    {loki->params.sample_freq[0]=(yyvsp[(3) - (5)].value); loki->params.sample_freq[1]=(yyvsp[(5) - (5)].value);}
    break;

  case 101:

/* Line 1455 of yacc.c  */
#line 258 "param_parse.y"
    {if(loki->names[LK_FREQFILE]) free(loki->names[LK_FREQFILE]); loki->names[LK_FREQFILE]=(yyvsp[(3) - (3)].string);}
    break;

  case 102:

/* Line 1455 of yacc.c  */
#line 259 "param_parse.y"
    {loki->params.sample_freq[0]=loki->params.sample_freq[1]=(yyvsp[(3) - (3)].value);}
    break;

  case 103:

/* Line 1455 of yacc.c  */
#line 260 "param_parse.y"
    {loki->params.sample_freq[0]=(yyvsp[(3) - (5)].value); loki->params.sample_freq[1]=(yyvsp[(5) - (5)].value);}
    break;

  case 104:

/* Line 1455 of yacc.c  */
#line 261 "param_parse.y"
    {loki->names[LK_PHENFILE]=(yyvsp[(3) - (3)].string); }
    break;

  case 105:

/* Line 1455 of yacc.c  */
#line 262 "param_parse.y"
    {set_output_gen((yyvsp[(3) - (3)].string),0);}
    break;

  case 106:

/* Line 1455 of yacc.c  */
#line 263 "param_parse.y"
    {set_output_gen((yyvsp[(3) - (5)].string),(yyvsp[(5) - (5)].string)); free((yyvsp[(5) - (5)].string)); }
    break;

  case 108:

/* Line 1455 of yacc.c  */
#line 265 "param_parse.y"
    {loki->params.output_type=(yyvsp[(3) - (3)].value);}
    break;

  case 109:

/* Line 1455 of yacc.c  */
#line 266 "param_parse.y"
    {if(loki->names[LK_OUTPUTFILE]) free(loki->names[LK_OUTPUTFILE]); loki->names[LK_OUTPUTFILE]=(yyvsp[(3) - (3)].string);}
    break;

  case 110:

/* Line 1455 of yacc.c  */
#line 267 "param_parse.y"
    {if(loki->names[LK_POSFILE]) free(loki->names[LK_POSFILE]); loki->names[LK_POSFILE]=(yyvsp[(4) - (4)].string);}
    break;

  case 111:

/* Line 1455 of yacc.c  */
#line 268 "param_parse.y"
    {if(loki->names[LK_IBDFILE]) free(loki->names[LK_IBDFILE]); loki->names[LK_IBDFILE]=(yyvsp[(4) - (4)].string);}
    break;

  case 112:

/* Line 1455 of yacc.c  */
#line 269 "param_parse.y"
    {if(loki->names[LK_IBDDIR]) free(loki->names[LK_IBDDIR]); loki->names[LK_IBDDIR]=(yyvsp[(4) - (4)].string);}
    break;

  case 113:

/* Line 1455 of yacc.c  */
#line 270 "param_parse.y"
    {if(loki->names[LK_DUMPFILE]) free(loki->names[LK_DUMPFILE]); loki->names[LK_DUMPFILE]=(yyvsp[(3) - (3)].string);}
    break;

  case 114:

/* Line 1455 of yacc.c  */
#line 271 "param_parse.y"
    {loki->params.dump_freq=(yyvsp[(3) - (3)].value);}
    break;

  case 115:

/* Line 1455 of yacc.c  */
#line 272 "param_parse.y"
    {loki->params.output_haplo=1;if(loki->names[LK_HAPLOFILE]) free(loki->names[LK_HAPLOFILE]); loki->names[LK_HAPLOFILE]=(yyvsp[(3) - (3)].string);}
    break;

  case 116:

/* Line 1455 of yacc.c  */
#line 273 "param_parse.y"
    {loki->params.output_haplo=1;}
    break;

  case 117:

/* Line 1455 of yacc.c  */
#line 274 "param_parse.y"
    {if(loki->names[LK_POLYFILE]) free(loki->names[LK_POLYFILE]); loki->names[LK_POLYFILE]=(yyvsp[(3) - (3)].string);}
    break;

  case 118:

/* Line 1455 of yacc.c  */
#line 275 "param_parse.y"
    {set_ibd_mode((yyvsp[(3) - (3)].string)); free((yyvsp[(3) - (3)].string));}
    break;

  case 119:

/* Line 1455 of yacc.c  */
#line 276 "param_parse.y"
    {set_ibd_mode((yyvsp[(3) - (5)].string)); free((yyvsp[(3) - (5)].string)); set_ibd_mode((yyvsp[(5) - (5)].string)); free((yyvsp[(5) - (5)].string));}
    break;

  case 120:

/* Line 1455 of yacc.c  */
#line 279 "param_parse.y"
    { if(!check_variance((yyvsp[(3) - (3)].rvalue),0)) {loki->models->res_var_set[0]=start_flag; loki->models->residual_var[0]=(yyvsp[(3) - (3)].rvalue);} }
    break;

  case 121:

/* Line 1455 of yacc.c  */
#line 280 "param_parse.y"
    { if((yyvsp[(3) - (4)].value)>=0 && !check_variance((yyvsp[(4) - (4)].rvalue),1)) {loki->models->res_var_set[(yyvsp[(3) - (4)].value)]=start_flag; BB(loki->models->residual_var,(yyvsp[(3) - (4)].value),(yyvsp[(3) - (4)].value))=(yyvsp[(4) - (4)].rvalue);} }
    break;

  case 122:

/* Line 1455 of yacc.c  */
#line 283 "param_parse.y"
    { if(!check_variance((yyvsp[(4) - (4)].rvalue),0)) loki->models->residual_var_limit[0]=(yyvsp[(4) - (4)].rvalue); }
    break;

  case 123:

/* Line 1455 of yacc.c  */
#line 284 "param_parse.y"
    { if((yyvsp[(4) - (5)].value)>=0 && !check_variance((yyvsp[(5) - (5)].rvalue),0)) loki->models->residual_var_limit[(yyvsp[(4) - (5)].value)]=(yyvsp[(5) - (5)].rvalue); }
    break;

  case 124:

/* Line 1455 of yacc.c  */
#line 285 "param_parse.y"
    { if(!check_variance((yyvsp[(4) - (4)].rvalue),0)) loki->models->residual_var_limit[0]=(yyvsp[(4) - (4)].rvalue); }
    break;

  case 125:

/* Line 1455 of yacc.c  */
#line 286 "param_parse.y"
    { if((yyvsp[(4) - (5)].value)>=0 && !check_variance((yyvsp[(5) - (5)].rvalue),0)) loki->models->residual_var_limit[(yyvsp[(4) - (5)].value)]=(yyvsp[(5) - (5)].rvalue); }
    break;

  case 126:

/* Line 1455 of yacc.c  */
#line 289 "param_parse.y"
    { if(!check_variance((yyvsp[(3) - (3)].rvalue),0)) {loki->models->add_var_set[0]=start_flag; loki->models->additive_var[0]=(yyvsp[(3) - (3)].rvalue);} }
    break;

  case 127:

/* Line 1455 of yacc.c  */
#line 290 "param_parse.y"
    { if((yyvsp[(3) - (4)].value)>=0 && !check_variance((yyvsp[(4) - (4)].rvalue),0)) {loki->models->add_var_set[(yyvsp[(3) - (4)].value)]=start_flag; BB(loki->models->additive_var,(yyvsp[(3) - (4)].value),(yyvsp[(3) - (4)].value))=(yyvsp[(4) - (4)].rvalue);} }
    break;

  case 128:

/* Line 1455 of yacc.c  */
#line 293 "param_parse.y"
    { if(!check_variance((yyvsp[(4) - (4)].rvalue),0)) loki->models->additive_var_limit[0]=(yyvsp[(4) - (4)].rvalue); }
    break;

  case 129:

/* Line 1455 of yacc.c  */
#line 294 "param_parse.y"
    { if((yyvsp[(4) - (5)].value)>=0 && !check_variance((yyvsp[(5) - (5)].rvalue),0)) loki->models->additive_var_limit[(yyvsp[(4) - (5)].value)]=(yyvsp[(5) - (5)].rvalue); }
    break;

  case 130:

/* Line 1455 of yacc.c  */
#line 295 "param_parse.y"
    { if(!check_variance((yyvsp[(4) - (4)].rvalue),0)) loki->models->additive_var_limit[0]=(yyvsp[(4) - (4)].rvalue); }
    break;

  case 131:

/* Line 1455 of yacc.c  */
#line 296 "param_parse.y"
    { if((yyvsp[(4) - (5)].value)>=0 && !check_variance((yyvsp[(5) - (5)].rvalue),0)) loki->models->additive_var_limit[(yyvsp[(4) - (5)].value)]=(yyvsp[(5) - (5)].rvalue); }
    break;

  case 132:

/* Line 1455 of yacc.c  */
#line 299 "param_parse.y"
    { if(loki->models->n_models>1) yyerror("Model must be specified when multiple models are present");
	                      else if(!loki->models->n_models) print_scan_warn("No model present - MEAN command ignored\n");
	                      else {loki->models->grand_mean_set[0]=start_flag; loki->models->grand_mean[0]=(yyvsp[(2) - (2)].rvalue); } }
    break;

  case 133:

/* Line 1455 of yacc.c  */
#line 302 "param_parse.y"
    { if(!loki->models->n_models) print_scan_warn("No model present - MEAN command ignored\n");
				                      else if((yyvsp[(2) - (3)].value)>=0) {loki->models->grand_mean_set[(yyvsp[(2) - (3)].value)]=start_flag; loki->models->grand_mean[(yyvsp[(2) - (3)].value)]=(yyvsp[(3) - (3)].rvalue);} }
    break;

  case 134:

/* Line 1455 of yacc.c  */
#line 306 "param_parse.y"
    { set_position((yyvsp[(2) - (3)].lk_var),(yyvsp[(3) - (3)].rvalue),(yyvsp[(3) - (3)].rvalue)); }
    break;

  case 135:

/* Line 1455 of yacc.c  */
#line 307 "param_parse.y"
    { loki->markers->sex_map=1; set_position((yyvsp[(2) - (5)].lk_var),(yyvsp[(3) - (5)].rvalue),(yyvsp[(5) - (5)].rvalue)); }
    break;

  case 136:

/* Line 1455 of yacc.c  */
#line 310 "param_parse.y"
    {c_flag=0;}
    break;

  case 138:

/* Line 1455 of yacc.c  */
#line 311 "param_parse.y"
    {c_flag=1;}
    break;

  case 140:

/* Line 1455 of yacc.c  */
#line 312 "param_parse.y"
    {c_flag=1;}
    break;

  case 142:

/* Line 1455 of yacc.c  */
#line 315 "param_parse.y"
    {(yyval.value)=find_trait((yyvsp[(1) - (1)].lk_var));}
    break;

  case 143:

/* Line 1455 of yacc.c  */
#line 317 "param_parse.y"
    { (yyval.lk_var)=find_var((yyvsp[(1) - (1)].string),0,0); }
    break;

  case 144:

/* Line 1455 of yacc.c  */
#line 318 "param_parse.y"
    { (yyval.lk_var)=find_var((yyvsp[(1) - (4)].string),(yyvsp[(3) - (4)].value),1); }
    break;

  case 145:

/* Line 1455 of yacc.c  */
#line 321 "param_parse.y"
    { if((yyvsp[(1) - (1)].lk_var)) {set_output((yyvsp[(1) - (1)].lk_var)); free((yyvsp[(1) - (1)].lk_var));} }
    break;

  case 146:

/* Line 1455 of yacc.c  */
#line 322 "param_parse.y"
    { if((yyvsp[(3) - (3)].lk_var)) {set_output((yyvsp[(3) - (3)].lk_var)); free((yyvsp[(3) - (3)].lk_var));} }
    break;

  case 147:

/* Line 1455 of yacc.c  */
#line 325 "param_parse.y"
    { freq_marker=check_marker((yyvsp[(1) - (1)].lk_var)); }
    break;

  case 148:

/* Line 1455 of yacc.c  */
#line 327 "param_parse.y"
    {group_ptr=0;}
    break;

  case 149:

/* Line 1455 of yacc.c  */
#line 327 "param_parse.y"
    {if(group_ptr<n_gen_grps) print_scan_err("Line %d: Too few groups in order statement\n",lineno1);}
    break;

  case 150:

/* Line 1455 of yacc.c  */
#line 330 "param_parse.y"
    { (yyval.value)=find_group(yytext,(yyvsp[(1) - (1)].value),1); }
    break;

  case 151:

/* Line 1455 of yacc.c  */
#line 331 "param_parse.y"
    { (yyval.value)=find_group(yytext,0,0); }
    break;

  case 152:

/* Line 1455 of yacc.c  */
#line 332 "param_parse.y"
    { (yyval.value)=find_group((yyvsp[(1) - (1)].string),0,0); free((yyvsp[(1) - (1)].string));}
    break;

  case 153:

/* Line 1455 of yacc.c  */
#line 335 "param_parse.y"
    {set_group_order((yyvsp[(1) - (1)].value));}
    break;

  case 154:

/* Line 1455 of yacc.c  */
#line 336 "param_parse.y"
    {set_group_order((yyvsp[(3) - (3)].value));}
    break;

  case 155:

/* Line 1455 of yacc.c  */
#line 339 "param_parse.y"
    { (yyval.value)=find_allele(yytext,freq_marker); }
    break;

  case 156:

/* Line 1455 of yacc.c  */
#line 340 "param_parse.y"
    { (yyval.value)=find_allele(yytext,freq_marker); }
    break;

  case 157:

/* Line 1455 of yacc.c  */
#line 341 "param_parse.y"
    { (yyval.value)=find_allele((yyvsp[(1) - (1)].string),freq_marker); free((yyvsp[(1) - (1)].string)); }
    break;

  case 158:

/* Line 1455 of yacc.c  */
#line 344 "param_parse.y"
    {group_counter=0; freq_allele=(yyvsp[(1) - (2)].value);}
    break;

  case 159:

/* Line 1455 of yacc.c  */
#line 344 "param_parse.y"
    {if(freq_marker && group_counter<n_gen_grps) print_scan_err("Line %d: Too few frequencies specified\n",lineno);}
    break;

  case 160:

/* Line 1455 of yacc.c  */
#line 345 "param_parse.y"
    {group_counter=0; freq_allele=(yyvsp[(2) - (3)].value);}
    break;

  case 161:

/* Line 1455 of yacc.c  */
#line 345 "param_parse.y"
    {if(freq_marker && group_counter<n_gen_grps) print_scan_err("Line %d: Too few frequencies specified\n",lineno);}
    break;

  case 162:

/* Line 1455 of yacc.c  */
#line 348 "param_parse.y"
    {set_freq(freq_marker,(yyvsp[(1) - (1)].rvalue),freq_allele);}
    break;

  case 163:

/* Line 1455 of yacc.c  */
#line 349 "param_parse.y"
    {group_counter++;}
    break;

  case 164:

/* Line 1455 of yacc.c  */
#line 350 "param_parse.y"
    {set_freq(freq_marker,(yyvsp[(3) - (3)].rvalue),freq_allele);}
    break;

  case 165:

/* Line 1455 of yacc.c  */
#line 351 "param_parse.y"
    {group_counter++;}
    break;

  case 166:

/* Line 1455 of yacc.c  */
#line 354 "param_parse.y"
    { (yyval.rlist)=add_ibd_list((yyvsp[(1) - (1)].rvalue),0); }
    break;

  case 167:

/* Line 1455 of yacc.c  */
#line 355 "param_parse.y"
    { (yyval.rlist)=add_ibd_list((yyvsp[(3) - (3)].rvalue),(yyvsp[(1) - (3)].rlist)); }
    break;

  case 168:

/* Line 1455 of yacc.c  */
#line 358 "param_parse.y"
    {add_real_to_num_array(((yyval.num_array)=make_num_array(8)),(yyvsp[(1) - (1)].rvalue));}
    break;

  case 169:

/* Line 1455 of yacc.c  */
#line 359 "param_parse.y"
    {add_int_to_num_array(((yyval.num_array)=make_num_array(8)),(yyvsp[(1) - (1)].value));}
    break;

  case 170:

/* Line 1455 of yacc.c  */
#line 360 "param_parse.y"
    { add_real_to_num_array((yyvsp[(1) - (3)].num_array),(yyvsp[(3) - (3)].rvalue)); }
    break;

  case 171:

/* Line 1455 of yacc.c  */
#line 361 "param_parse.y"
    { add_int_to_num_array((yyvsp[(1) - (3)].num_array),(yyvsp[(3) - (3)].value)); }
    break;

  case 173:

/* Line 1455 of yacc.c  */
#line 365 "param_parse.y"
    { (yyval.rvalue)=(double)(yyvsp[(1) - (1)].value); }
    break;



/* Line 1455 of yacc.c  */
#line 2802 "y.tab.c"
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
#line 368 "param_parse.y"


static void set_syst_var(int com,struct num_array *r)
{
	int i;
	double x[4];
	
	if(!r) ABT_FUNC("passed zero pointer\n");
	if(com==SYST_MSCAN_PROBS) {
		if(r->ptr!=4) yyerror("System variable MSCAN_PROBS requires 4 parameters\n");
		else {
			for(i=0;i<4;i++) x[i]=(r->x[i].flag==ST_INTEGER?(double)r->x[i].data.value:r->x[i].data.rvalue);
			i=set_mscan_probs(x);
			if(i) yyerror("Illegal parameters specified for MSCAN_PROBS\n");
		}
	} else {
		if(r->ptr!=1) yyerror("Multiple values not supported for this system variable\n");
		else loki->sys.syst_var[com]=r->x[0];
	}
}

static struct num_array *make_num_array(int n)
{
	struct num_array *r;
	
	if(n<1) n=8;
	if(!(r=malloc(sizeof(struct num_array)))) ABT_FUNC(MMsg);
	r->ptr=0;
	r->size=n;
	if(!(r->x=malloc(sizeof(struct id_data)*n))) ABT_FUNC(MMsg);
	return r;
}

static void free_num_array(struct num_array *r)
{
	if(!r) ABT_FUNC("passed zero pointer\n");
	if(r->x) free(r->x);
	free(r);
}

static void add_real_to_num_array(struct num_array *r,double x)
{
	if(r->ptr>=r->size) {
		if(r->ptr>r->size) ABT_FUNC("num_array corrupted?\n");
		r->size<<=1;
		if(!(r->x=realloc(r->x,sizeof(struct id_data)*r->size))) ABT_FUNC(MMsg);
	}
	r->x[r->ptr].data.rvalue=x;
	r->x[r->ptr++].flag=ST_REAL;
}

static void add_int_to_num_array(struct num_array *r,int i)
{
	if(r->ptr>=r->size) {
		if(r->ptr>r->size) ABT_FUNC("num_array corrupted?\n");
		r->size<<=1;
		if(!(r->x=realloc(r->x,sizeof(struct id_data)*r->size))) ABT_FUNC(MMsg);
	}
	r->x[r->ptr].data.value=i;
	r->x[r->ptr++].flag=ST_INTEGER;
}

static int find_trait(struct lk_variable *lkv)
{
	int i,type,mod;
	struct Variable *var;
	
	if(!loki->models->n_models || (lkv->type!=LK_TYPE_IDVAR && lkv->type!=LK_TYPE_NONIDVAR)) {
		yyerror("Not a trait variable");
		return -1;
	}
	for(mod=0;mod<loki->models->n_models;mod++) {
		type=loki->models->models[mod].var.type;
		i=loki->models->models[mod].var.var_index;
		var=(type&ST_CONSTANT)?loki->data->id_variable+i:loki->data->nonid_variable+i;
		if(var==lkv->var.var) break;
	}
	if(mod==loki->models->n_models) {
		yyerror("Not a trait variable");
		mod= -1;
	}
	return mod;
}

static void set_analyze(char *p)
{
	int i;
	char *com[]={"AFFECTED","NULL","IBD",0};
		
	if(p) {
		i=0;
		while(com[i]) {
			if(!strcasecmp(com[i],p)) break;
			i++;
		}
		if(com[i]) loki->params.analysis|=(1<<i);
		else yyerror("Invalid parameter to analyze statement");
	}
}

static void set_ibd_mode(char *p)
{
	int i;
	char *com[]={"LOKI","MERLIN","SOLAR","SINGLEPOINT","SINGLE","SINGLE_POINT",0};
		
	if(p) {
		i=0;
		while(com[i]) {
			if(!strcasecmp(com[i],p)) break;
			i++;
		}
		if(com[i]) {
			if(i<3) loki->params.ibd_mode=i;
			else loki->params.ibd_mode |=4;
		} else yyerror("Unknown ibd mode");
	}
}

static void set_group_order(int gp)
{
	int i;
	
	if(gp<0) return;
	for(i=0;i<group_ptr;i++) if(group_order[i]==gp)	{
		yyerror("Group repeated in order statement");
		return;
	}
	if(group_ptr>=n_gen_grps) yyerror("Too many groups - internal error?");
	else group_order[group_ptr++]=gp;
}

static struct IBD_List *add_ibd_list(double x,struct IBD_List *p)
{
	if(!p) {
		if(!(p=malloc(sizeof(struct IBD_List)))) ABT_FUNC(MMsg);
		p->idx=0;
		p->size=32;
		if(!(p->pos=malloc(sizeof(double)*p->size))) ABT_FUNC(MMsg);
	}
	if(p->size==p->idx) {
		p->size*=2;
		if(!(p->pos=realloc(p->pos,sizeof(double)*p->size))) ABT_FUNC(MMsg);
	}
	p->pos[p->idx++]=x;
	return p;
}

static int check_link_name(char *name)
{
	int i=0;
	
	if(name)	{
		for(i=0;i<loki->markers->n_links;i++) if(loki->markers->linkage[i].name && !strcasecmp(name,loki->markers->linkage[i].name)) break;
	} else for(i=0;i<loki->markers->n_links;i++) if(!loki->markers->linkage[i].name) break;
	if(i==loki->markers->n_links) {
		(void)fprintf(stderr,"Warning: linkage group %s not found; command ignored\n",name);
		i= -1;
	}
	return i;
}

static struct clause_atom *new_clause_atom(char *v,int i)
{
	struct clause_atom *c;
	
	c=lk_malloc(sizeof(struct clause_atom));
	c->next=0;
	c->v=v;
	c->i=i;
	return c;
}

static struct string_list *new_string_list(char *v)
{
	struct string_list *s;
	
	s=lk_malloc(sizeof(struct string_list));
	s->next=0;
	s->v=v;
	return s;
}

static void set_pseudo_chrome(struct string_list *s,struct clause_atom *c)
{
	int i;
	static char *opt[]={"INIT","FREQ",0};

	if(c) {
		while(c) {
			i=0;
			while(opt[i]) {
				if(!strcasecmp(c->v,opt[i])) break;
				i++;
			}
			switch(i) {
			 case 0:
				loki->params.pseudo_start=c->i;
				break;
			 case 1:
				loki->params.pseudo_freq=c->i;
				break;
			 default:
				(void)fprintf(stderr,"Warning: Unknown option %s in PSEUDOCHROMOSOME command; option ignored\n",c->v);
				break;
			}
			free(c->v);
			c=c->next;
		}
		free_list(c,0);
	}
	if(s) {
		while(s) {
			i=check_link_name(s->v);
			if(i>=0) {
				loki->markers->linkage[i].type|=LINK_MIRRORED;
				loki->params.pseudo_flag=1;
			}
			free(s->v);
			s=s->next;
		}
		free_list(s,0);
	} else {
		for(i=0;i<loki->markers->n_links;i++) loki->markers->linkage[i].type|=LINK_MIRRORED;
		loki->params.pseudo_flag=1;
	}
}
  
static void set_output_gen(char *file,char *link)
{
	int i;
	struct output_gen *p;

	if(link)	{
		i=check_link_name(link);
		if(i<0) return;
		i++;
	} else i=0;
	p=loki->params.Output_Gen;
	if(!(loki->params.Output_Gen=malloc(sizeof(struct output_gen)))) ABT_FUNC(MMsg);
	loki->params.Output_Gen->next=p;
	loki->params.Output_Gen->file=file;
	loki->params.Output_Gen->link_group=i;
}

static void check_previous_list(int i)
{
	if(loki->markers->linkage[i].ibd_est_type) {
		(void)fprintf(stderr,"Warning: overwriting previous IBD settings for linkage group %s\n",loki->markers->linkage[i].name);
		if(loki->markers->linkage[i].ibd_list) {
			free(loki->markers->linkage[i].ibd_list->pos);
			free(loki->markers->linkage[i].ibd_list);
		}
		loki->markers->linkage[i].ibd_list=0;
		loki->markers->linkage[i].ibd_est_type=0;
	}
}

static void set_ibd_list(char *name,struct IBD_List *p,int type)
{
	int i=0,k;
	
	if(type==IBD_EST_GRID) {
		i=-1;
		if(p->idx<3) (void)fprintf(stderr,"Warning: too few parameters (%d) for IBD Grid (3 required); IBD request ignored\n",p->idx);
		else if(p->idx>3) (void)fprintf(stderr,"Warning: too many parameters (%d) for IBD Grid (3 required); IBD request ignored\n",p->idx);
		else if(fabs(p->pos[2])<IBD_MIN_GRID_STEP) (void)fprintf(stderr,"Warning: step size (%g) for IBD Grid < IBD_MIN_GRID_STEP (%g) in loki_ibd.h ; IBD request ignored\n",p->pos[2],IBD_MIN_GRID_STEP);
		else {
			k=1+(int)(.5+(p->pos[1]-p->pos[0])/p->pos[2]);
			if(k>IBD_MAX_GRID) (void)fprintf(stderr,"Warning: grid evaluations requested (%d) for IBD Grid > IBD_MAX_GRID (%d) in loki_ibd.h ; IBD request ignored\n",k,IBD_MAX_GRID);
			else i=0;
		}
	}
	if(!i) i=check_link_name(name);
	if(i<0) {
		free(p->pos);
		free(p);
	} else {
		check_previous_list(i);
		loki->markers->linkage[i].ibd_list=p;
		loki->markers->linkage[i].ibd_est_type=type;
	}
}

static void set_ibd_markers(char *name) 
{
	int i;
	
	i=check_link_name(name);
	if(i>=0) {
		for(;i<loki->markers->n_links;i++) {
			check_previous_list(i);
			loki->markers->linkage[i].ibd_est_type=IBD_EST_MARKERS;
			if(name) break;
		}
	}
}

static void set_tloci(int a,int b)
{
	if(a<0) {
		if(a== -1) loki->params.min_tloci=loki->params.max_tloci=b;
		else loki->params.start_tloci=b;
	} else if(a<b) {
		loki->params.min_tloci=a;
		loki->params.max_tloci=b;
	} else {
		loki->params.min_tloci=b;
		loki->params.max_tloci=a;
	}
}

static struct Marker *check_marker(struct lk_variable *lkvar)
{
	struct Marker *mk=0;
	
	if(!lkvar) return 0;
	if(lkvar->type!=LK_TYPE_MARKER)
		yyerror("Attempting to set frequency of a non-marker");
	else mk=lkvar->var.marker;
	free(lkvar);
	return mk;
}

static void set_output(struct lk_variable *lkvar)
{
	int i,j,type,mod;
	struct Variable *var;
	
	if(!lkvar) return;
	for(mod=0;mod<loki->models->n_models;mod++) {
		for(i=0;i<loki->models->models[mod].n_terms;i++) {
			type=loki->models->models[mod].term[i].vars[0].type;
			j=loki->models->models[mod].term[i].vars[0].var_index;
			if(type&ST_MARKER) {
				if(lkvar->type==LK_TYPE_MARKER && lkvar->var.marker==loki->markers->marker+j) {
					loki->models->models[mod].term[i].out_flag=1;
				}
			} else if(type&(ST_TRAITLOCUS|ST_ID|ST_SIRE|ST_DAM)) continue;
			else {
				if(type&ST_CONSTANT) var=loki->data->id_variable+j;
				else var=loki->data->nonid_variable+j;
				if((lkvar->type==LK_TYPE_IDVAR || lkvar->type==LK_TYPE_NONIDVAR) && lkvar->var.var==var) {
					loki->models->models[mod].term[i].out_flag=1;
				}
			}
		}
	}
}

static void set_map_range(char *name,double r1,double r2,int flag)
{
	int i;
	double t;
	static char *sexstr[2]={"female","male"};
	
	if(flag!= -1) loki->markers->sex_map=1;
	if(name)	{
		i=check_link_name(name);
		if(i>=0) {
			if(r2<r1) {
				t=r2;
				r2=r1;
				r1=t;
			}
			if(flag== -1) {
				loki->markers->linkage[i].r1[0]=loki->markers->linkage[i].r1[1]=r1;
				loki->markers->linkage[i].r2[0]=loki->markers->linkage[i].r2[1]=r2;
				loki->markers->linkage[i].range_set[0]=loki->markers->linkage[i].range_set[1]=1;
				message(INFO_MSG,"Map range for linkage group '%s' set to %g-%gcM\n",name,r1,r2);
			} else {
				loki->markers->linkage[i].r1[flag]=r1;
				loki->markers->linkage[i].r2[flag]=r2;
				loki->markers->linkage[i].range_set[flag]=1;
				message(INFO_MSG,"Map range (%s) for linkage group '%s' set to %g-%gcM\n",sexstr[flag],name,r1,r2);
			}
		}
	} else {
		if(flag<0) {
			loki->markers->total_maplength[X_PAT]=r1;
			loki->markers->total_maplength[X_MAT]=r2;
			message(INFO_MSG,"Total (genome) map length set to (%g,%g)cM\n",r1,r2);
		} else {
			loki->markers->total_maplength[flag]=r1;
			message(INFO_MSG,"Total (genome) %s map length set to %gcM\n",sexstr[flag],r1);
		}
	}
}

static int find_group(char *p,int gp, int flag)
{
	int i;
	char *s;
	struct Id_Recode *rec;
	
	rec=&loki->pedigree->group_recode;
	if(!rec->recode) return -1;
	if(rec->flag==ST_STRING) {
		for(i=0;i<n_gen_grps;i++) if(!(strcasecmp(p,rec->recode[i].string))) return i;
	} else {
		if(!flag) {
			gp=strtol(p,&s,10);
			if(!(*s)) for(i=0;i<n_gen_grps;i++) if(gp==rec->recode[i].value) return i;
		} else for(i=0;i<n_gen_grps;i++) if(gp==rec->recode[i].value) return i;
	}
	yyerror("Group not found\n");
	return -1;
}

static int find_allele( char *p, struct Marker *mk)
{
	int i,j;
	
	if(!mk) return -1;
	j=mk->locus.n_alleles-1;
	for(i=0;i<j;i++) if(!(strcmp(p,mk->recode[i]))) return i;
	return -1;
}

static struct lk_variable *find_var(char *p, int idx, int flag)
{
	int i,j=0;
	struct lk_variable *lkv=0;
	char *p1;
	
	if(flag) {
		if(idx<1) {
			yyerror("Illegal array index");
			free(p);
			return 0;
		}
		i=(int)strlen(p)+(int)(log((double)idx)/log(10.0)+4.00001);
		if(!(p1=malloc((size_t)i))) ABT_FUNC(MMsg);
		snprintf(p1,i,"%s(%d)",p,idx);
		free(p);
		p=p1;
	}
	for(i=0;i<loki->markers->n_markers;i++) if(!strcasecmp(loki->markers->marker[i].name,p)) {
		j=LK_TYPE_MARKER;
		break;
	}
	if(!j) for(i=0;i<loki->data->n_id_records;i++) if(!strcasecmp(loki->data->id_variable[i].name,p)) {
		j=LK_TYPE_IDVAR;
		break;
	}
	if(!j) for(i=0;i<loki->data->n_nonid_records;i++) if(!strcasecmp(loki->data->nonid_variable[i].name,p)) {
		j=LK_TYPE_NONIDVAR;
		break;
	}
	if(!j) for(i=0;i<loki->markers->n_links;i++) if(loki->markers->linkage[i].name && !strcasecmp(loki->markers->linkage[i].name,p)) {
		j=LK_TYPE_LINK;
		break;
	}
	if(j) {
		if(!(lkv=malloc(sizeof(struct lk_variable)))) ABT_FUNC(MMsg);
		lkv->type=j;
		switch(j) {
		 case LK_TYPE_MARKER:
			lkv->var.marker=loki->markers->marker+i;
			break;
		 case LK_TYPE_LINK:
			lkv->var.link=loki->markers->linkage+i;
			break;
		 case LK_TYPE_IDVAR:
			lkv->var.var=loki->data->id_variable+i;
			break;
		 case LK_TYPE_NONIDVAR:
			lkv->var.var=loki->data->nonid_variable+i;
			break;
		}
	}
	free(p);
	return lkv;
}

static void set_freq(struct Marker *mk, double freq,int allele)
{
	static int fg,fg1;
	int i=0;
	
	group_counter++;
	if(mk) {
		if(n_gen_grps>1) {
			if(!group_ptr)	{
				if(!fg) yyerror("Genetic group order not set");
				fg=1;
				return;
			}
			if(group_counter>n_gen_grps) {
				if(!fg1) yyerror("Too many frequencies specified (only 1 per genetic group)");
				fg1=1;
				return;
			}
			i=group_order[group_counter-1];
			if(i<0) return;
		}
		if(freq<0.0) {
			yyerror("Invalid (negative) frequency");
			return;
		}
		if(allele>=0 && freq==0.0)	{
			yyerror("Can not set frequency of observed allele to zero\n");
			return;
		}
		if(allele>=0) {
			mk->locus.freq[i][allele]=freq;
		} else {
			allele=mk->locus.n_alleles-1;
			mk->locus.freq[i][allele]+=freq;
		}
		mk->freq_set[i][allele]=start_flag;
		mk->count_flag[i]=c_flag;
	}
	fg1=0;
}

static void set_position(struct lk_variable *lkvar, double pos1, double pos2)
{
	if(lkvar) {
		if(lkvar->type!=LK_TYPE_MARKER) {
			yyerror("Attempting to set position of a non-marker");
			return;
		}
		lkvar->var.marker->locus.pos[X_PAT]=pos1;
		lkvar->var.marker->locus.pos[X_MAT]=pos2;
		lkvar->var.marker->pos_set=start_flag;
		free(lkvar);
	}
}

static int check_variance(const double v,const int fg)
{
	if(v<=0.0) {
		yyerror("Variance must be positive");
		return 1;
	}
	if(loki->models->n_models>1 && !fg) {
		yyerror("Must specify which model when multiple models are present");
		return 1;
	}
	if(!loki->models->n_models) {
		print_scan_warn("No model present - VARIANCE command ignored\n");
		return 1;
	}
	return 0;
}

void print_scan_err(char *fmt, ...)
{
	va_list args;
	
	va_start(args,fmt);
	(void)vfprintf(stderr,fmt,args);
	va_end(args);
	if((++scan_error_n)>=max_scan_errors) abt(__FILE__,__LINE__,"Too many errors - aborting\n");
}

void print_scan_warn(char *fmt, ...)
{
	va_list args;
	
	if(scan_warn_n<max_scan_warnings)
	{
		va_start(args,fmt);
		(void)vfprintf(stderr,fmt,args);
		va_end(args);
	}
	scan_warn_n++;
}

void yyerror(char *s)
{
     int i;
	
	print_scan_err("Line %d: %s\n%s\n",lineno,s,linebuf);
	if(scan_error_n<=max_scan_errors)
	{
		for(i=1;i<tokenpos;i++) (void)putc('-',stderr);
		(void)fputs("^\n",stderr);
	}
}

int yywrap(void)
{
	return 1;
}

int ReadParam(FILE *fptr,char *cname,struct loki *lk_struct)
{
	int i,j;
	void yy_cleanup(void);

#if YYDEBUG
	yydebug=1;
#endif
	loki=lk_struct;
	loki->params.start_tloci=-1;
	start_flag=1;
	loki->params.max_tloci=DEFAULT_MAX_TLOCI;
	loki->params.prune_option=PRUNE_LOCUS_SPECIFIC;
	loki->params.analysis=0;
	loki->params.output_type=DEFAULT_OUTPUT_TYPE;
	for(i=0;i<2;i++) loki->params.sample_freq[i]=1;
	fname_list[0]=cname;
	list_ptr=0;
	for(i=0;i<NUM_SYSTEM_VAR;i++) loki->sys.syst_var[i].flag=0;
	n_gen_grps=loki->pedigree->n_genetic_groups;
	if(!(group_order=malloc(sizeof(int)*n_gen_grps))) ABT_FUNC(MMsg);
	yyin=fptr;
	if((i=yyparse())) {
	  (void)fprintf(stderr,"Error: yyparse returned error %d\n",i);
	  scan_error_n++;
	}
	yy_cleanup();
	if(group_order) free(group_order);
	if(loki->params.start_tloci<0) loki->params.start_tloci=loki->params.min_tloci;
	else if(loki->params.start_tloci<loki->params.min_tloci || loki->params.start_tloci>loki->params.max_tloci) {
		(void)fprintf(stderr,"ReadParam(): Starting no. trait loci (%d) is outside set range (%d-%d)\n",loki->params.start_tloci,loki->params.min_tloci,loki->params.max_tloci);
		scan_error_n++;
	}
	if(loki->models->n_models>1 && loki->params.output_type<2) {
		(void)fprintf(stderr,"ReadParam(): Ouput type %d not supported with multilpe trait loci\n",loki->params.output_type);
		scan_error_n++;
	}
	if(!loki->sys.syst_var[SYST_IBD_OUTPUT].flag) {
		j=loki->params.ibd_mode;
		if(loki->params.compress_ibd) j|=COMPRESS_IBD;
		if(j) {
			loki->sys.syst_var[SYST_IBD_OUTPUT].flag=ST_INTEGER;
			loki->sys.syst_var[SYST_IBD_OUTPUT].data.value=j;
		}
	}
	if(scan_error_n) return 1;
	return 0;
}

