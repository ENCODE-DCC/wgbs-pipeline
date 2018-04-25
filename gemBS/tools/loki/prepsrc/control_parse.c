
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
#line 1 "control_parse.y"

/****************************************************************************
 *                                                                          *
 *     Loki - Programs for genetic analysis of complex traits using MCMC    *
 *                                                                          *
 *             Simon Heath - University of Washington                       *
 *                                                                          *
 *                       March 1997                                         *
 *                                                                          *
 * control_parse.y:                                                         *
 *                                                                          *
 * yacc source for control file parser.                                     *
 *                                                                          *
 * Copyright (C) Simon C. Heath 1997, 2000, 2002                            *
 * This is free software.  You can distribute it and/or modify it           *
 * under the terms of the Modified BSD license, see the file COPYING        *
 *                                                                          *
 ****************************************************************************/

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

#include "config.h"
#include "utils.h"
#include "lk_malloc.h"
	
#ifndef YYDEBUG
#define YYDEBUG 0
#endif


/* Line 189 of yacc.c  */
#line 113 "y.tab.c"

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
     FILEC = 258,
     MARKER = 259,
     LOCUS = 260,
     TRAIT = 261,
     RANDOM = 262,
     PEDIGREE = 263,
     LOG = 264,
     SUPER = 265,
     MODEL = 266,
     FILTER = 267,
     LINK = 268,
     MISSING = 269,
     FACTOR = 270,
     BREAK = 271,
     DOLOOP = 272,
     WHILE = 273,
     USE = 274,
     WHERE = 275,
     ORSYMBOL = 276,
     ANDSYMBOL = 277,
     NEQSYMBOL = 278,
     LEQSYMBOL = 279,
     GEQSYMBOL = 280,
     NOTSYMBOL = 281,
     LOGICAL = 282,
     SHELL = 283,
     ARRAY = 284,
     PRINTEXP = 285,
     INCLUDE = 286,
     RAWOUTPUT = 287,
     LOOP_CLAUSE_START = 288,
     LOOP_CLAUSE_END = 289,
     CONSTANT = 290,
     MULTIPLE = 291,
     RSFORMAT = 292,
     FSFORMAT = 293,
     SKIPFORMAT = 294,
     GSFORMAT = 295,
     CENSORED = 296,
     GROUP = 297,
     SET = 298,
     GENDER = 299,
     AFFECTED = 300,
     OUTPUT = 301,
     ERRORDIR = 302,
     LAUROUTPUT = 303,
     UNAFFECTED = 304,
     POSITION = 305,
     FREQUENCY = 306,
     PROBAND = 307,
     STRING = 308,
     VARIABLE = 309,
     ASSIGN = 310,
     ARRAY_VAR = 311,
     INTEGER = 312,
     SYSTEM_VAR = 313,
     REAL = 314,
     UMINUS = 315
   };
#endif
/* Tokens.  */
#define FILEC 258
#define MARKER 259
#define LOCUS 260
#define TRAIT 261
#define RANDOM 262
#define PEDIGREE 263
#define LOG 264
#define SUPER 265
#define MODEL 266
#define FILTER 267
#define LINK 268
#define MISSING 269
#define FACTOR 270
#define BREAK 271
#define DOLOOP 272
#define WHILE 273
#define USE 274
#define WHERE 275
#define ORSYMBOL 276
#define ANDSYMBOL 277
#define NEQSYMBOL 278
#define LEQSYMBOL 279
#define GEQSYMBOL 280
#define NOTSYMBOL 281
#define LOGICAL 282
#define SHELL 283
#define ARRAY 284
#define PRINTEXP 285
#define INCLUDE 286
#define RAWOUTPUT 287
#define LOOP_CLAUSE_START 288
#define LOOP_CLAUSE_END 289
#define CONSTANT 290
#define MULTIPLE 291
#define RSFORMAT 292
#define FSFORMAT 293
#define SKIPFORMAT 294
#define GSFORMAT 295
#define CENSORED 296
#define GROUP 297
#define SET 298
#define GENDER 299
#define AFFECTED 300
#define OUTPUT 301
#define ERRORDIR 302
#define LAUROUTPUT 303
#define UNAFFECTED 304
#define POSITION 305
#define FREQUENCY 306
#define PROBAND 307
#define STRING 308
#define VARIABLE 309
#define ASSIGN 310
#define ARRAY_VAR 311
#define INTEGER 312
#define SYSTEM_VAR 313
#define REAL 314
#define UMINUS 315




#if ! defined YYSTYPE && ! defined YYSTYPE_IS_DECLARED
typedef union YYSTYPE
{

/* Line 214 of yacc.c  */
#line 41 "control_parse.y"

	char *string;
	struct bin_node *var;
	int value;
	double rvalue;
	struct format_clause *format_clause;
	struct fformat *fformat;
	struct format_atom *f_atom;
	struct model_list *model_list;
	struct var_list *var_list;
	struct var_element *element;
	struct express *express;



/* Line 214 of yacc.c  */
#line 285 "y.tab.c"
} YYSTYPE;
# define YYSTYPE_IS_TRIVIAL 1
# define yystype YYSTYPE /* obsolescent; will be withdrawn */
# define YYSTYPE_IS_DECLARED 1
#endif


/* Copy the second part of user declarations.  */

/* Line 264 of yacc.c  */
#line 87 "control_parse.y"

#include "scan.h"
#include "scanner.h"

static struct format_atom *make_f_atom(int,int);
static struct format_clause *add_f_atom(struct format_clause *,struct format_atom *);
static struct format_clause *add_f_list(struct format_clause *,struct format_clause *,int);
static struct format *setup_format(struct format_clause *);
static struct bin_node *create_var(char *);
static struct model_list *add_to_model(struct model_list *,struct var_list *);
static struct var_list *add_to_var_list(struct var_list *,struct bin_node *,struct express *);
static struct var_list *add_var_lists(struct var_list *,struct var_list *);
static struct var_element *get_element(struct bin_node *,struct express *);
static struct var_element *assign_var(struct bin_node *,struct express *,struct express *);
static struct express *alloc_express(void);
static struct express *do_express_op(struct express *,struct express *,int);
static struct express *do_logical_op(struct express *,struct express *,int);
static struct fformat *add_fformat(struct fformat *,struct fformat *);
static struct fformat *create_fformat(void *,int);
static void begin_looping(struct var_element *, struct express *, struct express *);
static void free_vlist(struct var_list *);
static void do_ped_com(struct var_list *);
static void set_position(struct var_element *,struct express *,struct express *,struct express *);
static void add_restriction(struct var_list *);
static void add_censored(struct var_element *,const int);
static void do_file_com(char *,struct format_clause *,struct fformat *,struct var_list *);
static void set_locus_array(struct bin_node *);
static void set_locus_element(struct var_element *);
static void set_slocus_element(struct var_element *,struct var_list *);
static void set_haplo_element(struct var_element *,struct var_element *);
static void do_link_com(char *,int, struct var_list *);
static void do_missing_com(struct express *,struct var_list *,char *);
static void change_type(int,struct var_list *);
static void do_model_com(struct model_list *,struct bin_node *,struct express *);
static void add_operation(void *,int,int);
static void set_array_var(struct scan_data *,struct express *);
static void check_element_add_op(struct var_element *);
static void enter_loop(void);
static void start_loopclause(void);
static void do_while_com(struct express *);
static void print_exp(struct express *);
static void new_command(void);
static void set_sex(struct var_element *,struct express *,struct express *);
static void set_group(struct var_element *);
static int count_var_list(struct var_list *),shell_flag;
static struct format_atom *f_atom_list;
static struct var_element *pedlist[4];
static int f_atom_n,f_atom_size=32,pedflag;
struct operation *Affected,*Unaffected,*Proband;
struct bin_node *root_var;
struct InFile *Infiles;
struct Link *links;
struct Miss *Miss;
struct Restrict *Restrictions;
struct Censor *Censored;
struct model *Models;
static struct operation *Op_List;
struct Marker *markers,*traitlocus;
static struct var_element *hap_list[2];
struct express *sex_exp[2];
struct sex_def *sex_def;
struct var_element *group_elem;
static char *string_copy(char *s1,char *s2);
static char *LogFile;
static struct marker_info *m_info;
	
int scan_error,scan_error_n,scan_warn_n;
int max_scan_errors=30,max_scan_warnings=30,n_markers,iflag,file_skip;
char *Filter,*ErrorDir,*rsformat,*fsformat,*gsformat,*OutputFile,*OutputRawFile,*OutputLaurFile;
int loop_level,loop_ptr[MAX_LOOP],loop_stat[MAX_LOOP],loop_record,loop_stack_size=256;
int loop_main_ptr,in_loopclause,loop_clause_end,loop_clause_step,loop_clause_ptr;
int syst_var[NUM_PREP_SYSTEM_VAR];
int family_id;
struct var_element *loop_clause_element;
struct token_store *loop_stack;
	  


/* Line 264 of yacc.c  */
#line 376 "y.tab.c"

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
#define YYFINAL  140
/* YYLAST -- Last index in YYTABLE.  */
#define YYLAST   910

/* YYNTOKENS -- Number of terminals.  */
#define YYNTOKENS  77
/* YYNNTS -- Number of nonterminals.  */
#define YYNNTS  61
/* YYNRULES -- Number of rules.  */
#define YYNRULES  208
/* YYNRULES -- Number of states.  */
#define YYNSTATES  417

/* YYTRANSLATE(YYLEX) -- Bison symbol number corresponding to YYLEX.  */
#define YYUNDEFTOK  2
#define YYMAXUTOK   315

#define YYTRANSLATE(YYX)						\
  ((unsigned int) (YYX) <= YYMAXUTOK ? yytranslate[YYX] : YYUNDEFTOK)

/* YYTRANSLATE[YYLEX] -- Bison symbol number corresponding to YYLEX.  */
static const yytype_uint8 yytranslate[] =
{
       0,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
      69,    70,    65,    63,    71,    64,    66,    67,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,    74,
      61,    60,    62,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,    72,     2,    73,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
      75,    76,     2,     2,     2,     2,     2,     2,     2,     2,
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
      45,    46,    47,    48,    49,    50,    51,    52,    53,    54,
      55,    56,    57,    58,    59,    68
};

#if YYDEBUG
/* YYPRHS[YYN] -- Index of the first RHS symbol of rule number YYN in
   YYRHS.  */
static const yytype_uint16 yyprhs[] =
{
       0,     0,     3,     5,     7,    10,    14,    18,    20,    22,
      24,    26,    28,    30,    32,    34,    36,    38,    40,    42,
      44,    46,    48,    50,    52,    54,    56,    58,    60,    62,
      64,    66,    68,    70,    72,    76,    83,    85,    87,    89,
      91,    95,    99,   103,   107,   110,   114,   118,   122,   124,
     125,   129,   136,   142,   148,   154,   160,   165,   169,   175,
     183,   188,   194,   196,   200,   204,   208,   212,   216,   220,
     224,   228,   231,   234,   236,   240,   243,   246,   249,   252,
     255,   260,   267,   274,   280,   287,   292,   294,   296,   298,
     300,   304,   308,   312,   316,   320,   323,   327,   331,   335,
     339,   343,   347,   351,   355,   358,   363,   371,   379,   381,
     386,   388,   392,   397,   404,   406,   410,   414,   417,   420,
     423,   426,   428,   431,   433,   436,   438,   441,   444,   447,
     450,   453,   456,   460,   465,   471,   474,   478,   482,   486,
     489,   492,   495,   498,   501,   504,   506,   511,   516,   519,
     524,   527,   530,   535,   543,   546,   550,   554,   556,   558,
     562,   566,   570,   574,   576,   578,   580,   582,   587,   589,
     594,   596,   600,   603,   605,   607,   611,   617,   621,   623,
     627,   632,   634,   639,   641,   645,   647,   656,   667,   672,
     674,   675,   677,   679,   684,   686,   688,   692,   696,   698,
     700,   704,   708,   710,   714,   718,   720,   722,   726
};

/* YYRHS -- A `-1'-separated list of the rules' RHS.  */
static const yytype_int16 yyrhs[] =
{
      78,     0,    -1,    79,    -1,     1,    -1,    78,    16,    -1,
      78,    16,    79,    -1,    78,    16,     1,    -1,   101,    -1,
     116,    -1,    85,    -1,   107,    -1,   111,    -1,   112,    -1,
     113,    -1,   115,    -1,   117,    -1,   110,    -1,    98,    -1,
      99,    -1,    80,    -1,    96,    -1,    93,    -1,    84,    -1,
      90,    -1,    87,    -1,    88,    -1,   118,    -1,    83,    -1,
      95,    -1,    89,    -1,   108,    -1,   109,    -1,    91,    -1,
      81,    -1,    55,    60,    82,    -1,    55,    69,    82,    70,
      60,    82,    -1,    59,    -1,    57,    -1,    53,    -1,   123,
      -1,    82,    63,    82,    -1,    82,    64,    82,    -1,    82,
      65,    82,    -1,    82,    67,    82,    -1,    64,    82,    -1,
      69,    82,    70,    -1,    43,    58,    57,    -1,    43,   121,
      57,    -1,    17,    -1,    -1,    31,    86,   137,    -1,    41,
     123,    20,    69,   100,    70,    -1,    45,    20,    69,   100,
      70,    -1,    49,    20,    69,   100,    70,    -1,    52,    20,
      69,   100,    70,    -1,    44,   123,    82,    71,    82,    -1,
      18,    69,    92,    70,    -1,    50,   123,    82,    -1,    50,
     123,    82,    71,    82,    -1,    50,   123,    82,    71,    82,
      71,    82,    -1,    50,   123,    71,    82,    -1,    50,   123,
      71,    71,    82,    -1,    82,    -1,    92,    60,    92,    -1,
      92,    23,    92,    -1,    92,    24,    92,    -1,    92,    25,
      92,    -1,    92,    61,    92,    -1,    92,    62,    92,    -1,
      92,    21,    92,    -1,    92,    22,    92,    -1,    26,    92,
      -1,    30,    94,    -1,    82,    -1,    94,    71,    82,    -1,
      37,   137,    -1,    38,   137,    -1,    40,   137,    -1,    39,
      57,    -1,    29,    97,    -1,    54,    69,    82,    70,    -1,
      97,    71,    54,    69,    82,    70,    -1,    19,   135,    20,
      69,   100,    70,    -1,    19,    20,    69,   100,    70,    -1,
      20,    69,   100,    70,    19,   135,    -1,    20,    69,   100,
      70,    -1,    57,    -1,    59,    -1,    53,    -1,   123,    -1,
     100,    63,   100,    -1,   100,    64,   100,    -1,   100,    65,
     100,    -1,   100,    67,   100,    -1,    69,   100,    70,    -1,
      64,   100,    -1,   100,    60,   100,    -1,   100,    23,   100,
      -1,   100,    24,   100,    -1,   100,    25,   100,    -1,   100,
      61,   100,    -1,   100,    62,   100,    -1,   100,    21,   100,
      -1,   100,    22,   100,    -1,    26,   100,    -1,     3,   102,
      71,   134,    -1,     3,    72,   103,    73,   102,    71,   134,
      -1,     3,    72,   104,    73,   102,    71,   134,    -1,   137,
      -1,    28,    69,   137,    70,    -1,   106,    -1,   103,    71,
     106,    -1,    57,    69,   103,    70,    -1,   103,    71,    57,
      69,   103,    70,    -1,   105,    -1,   104,    71,   105,    -1,
     104,    74,   105,    -1,    38,   137,    -1,    37,   137,    -1,
      40,   137,    -1,    39,    57,    -1,    57,    -1,    57,    75,
      -1,    75,    -1,    57,     1,    -1,     1,    -1,     9,   137,
      -1,    46,   137,    -1,    48,   137,    -1,    32,   137,    -1,
      47,   137,    -1,    14,    82,    -1,    14,    82,   135,    -1,
      14,    82,    71,   135,    -1,    14,    72,   137,    73,    82,
      -1,     8,   135,    -1,     4,     5,   124,    -1,     6,     5,
     135,    -1,    10,     5,   127,    -1,    35,   135,    -1,    36,
     135,    -1,     7,   135,    -1,    15,   135,    -1,    59,   135,
      -1,    57,   135,    -1,    13,    -1,    13,    72,    75,    73,
      -1,    13,    72,    76,    73,    -1,   114,   135,    -1,   114,
     136,    71,   135,    -1,   114,   136,    -1,    12,   137,    -1,
      11,   121,    60,   119,    -1,    11,    56,    69,    82,    70,
      60,   119,    -1,    42,   123,    -1,   119,    63,   122,    -1,
     119,    63,   120,    -1,   120,    -1,   122,    -1,   122,    65,
     122,    -1,   122,    66,   122,    -1,   120,    65,   122,    -1,
     120,    66,   122,    -1,    54,    -1,    44,    -1,   121,    -1,
      56,    -1,    56,    69,    82,    70,    -1,   121,    -1,    56,
      69,    82,    70,    -1,   125,    -1,   124,    71,   125,    -1,
     123,   126,    -1,   123,    -1,    56,    -1,    72,   123,    73,
      -1,    72,   123,    71,   123,    73,    -1,    72,     1,    73,
      -1,   128,    -1,   127,    71,   125,    -1,   123,    72,   135,
      73,    -1,   123,    -1,   123,    72,     1,    73,    -1,   122,
      -1,   129,    71,   122,    -1,    69,    -1,   130,   129,    71,
      16,    81,    71,    82,    70,    -1,   130,   129,    71,    16,
      81,    71,    82,    71,    82,    70,    -1,   131,    33,   129,
      34,    -1,   131,    -1,    -1,   121,    -1,    56,    -1,    56,
      69,    82,    70,    -1,   133,    -1,   132,    -1,   134,    71,
     133,    -1,   134,    71,   132,    -1,   122,    -1,   132,    -1,
     135,    71,   122,    -1,   135,    71,   132,    -1,    53,    -1,
     136,    63,    53,    -1,   123,    63,   136,    -1,    53,    -1,
     123,    -1,   137,    63,    53,    -1,   137,    63,   123,    -1
};

/* YYRLINE[YYN] -- source line where rule number YYN was defined.  */
static const yytype_uint16 yyrline[] =
{
       0,   167,   167,   168,   169,   170,   171,   174,   175,   176,
     177,   178,   179,   180,   181,   182,   183,   184,   185,   186,
     187,   188,   189,   190,   191,   192,   193,   194,   195,   196,
     197,   198,   199,   202,   205,   206,   209,   210,   211,   212,
     214,   215,   216,   217,   218,   219,   222,   223,   226,   229,
     229,   232,   235,   236,   237,   240,   243,   246,   247,   248,
     249,   250,   253,   254,   255,   256,   257,   258,   259,   260,
     261,   262,   265,   268,   269,   272,   273,   274,   275,   278,
     280,   281,   284,   285,   288,   289,   292,   293,   294,   295,
     296,   297,   298,   299,   300,   301,   302,   303,   304,   305,
     306,   307,   308,   309,   310,   313,   314,   315,   318,   319,
     322,   323,   324,   325,   328,   329,   330,   333,   334,   335,
     336,   339,   340,   341,   342,   343,   346,   350,   351,   352,
     355,   359,   360,   361,   362,   365,   368,   369,   370,   373,
     374,   375,   376,   377,   378,   381,   382,   383,   386,   387,
     388,   391,   399,   400,   403,   406,   407,   408,   409,   412,
     413,   414,   415,   418,   419,   422,   423,   424,   427,   428,
     431,   431,   433,   434,   435,   438,   439,   440,   443,   443,
     445,   446,   447,   450,   451,   454,   457,   458,   461,   462,
     465,   466,   467,   468,   471,   472,   473,   474,   477,   478,
     479,   480,   483,   484,   485,   488,   489,   490,   491
};
#endif

#if YYDEBUG || YYERROR_VERBOSE || YYTOKEN_TABLE
/* YYTNAME[SYMBOL-NUM] -- String name of the symbol SYMBOL-NUM.
   First, the terminals, then, starting at YYNTOKENS, nonterminals.  */
static const char *const yytname[] =
{
  "$end", "error", "$undefined", "FILEC", "MARKER", "LOCUS", "TRAIT",
  "RANDOM", "PEDIGREE", "LOG", "SUPER", "MODEL", "FILTER", "LINK",
  "MISSING", "FACTOR", "BREAK", "DOLOOP", "WHILE", "USE", "WHERE",
  "ORSYMBOL", "ANDSYMBOL", "NEQSYMBOL", "LEQSYMBOL", "GEQSYMBOL",
  "NOTSYMBOL", "LOGICAL", "SHELL", "ARRAY", "PRINTEXP", "INCLUDE",
  "RAWOUTPUT", "LOOP_CLAUSE_START", "LOOP_CLAUSE_END", "CONSTANT",
  "MULTIPLE", "RSFORMAT", "FSFORMAT", "SKIPFORMAT", "GSFORMAT", "CENSORED",
  "GROUP", "SET", "GENDER", "AFFECTED", "OUTPUT", "ERRORDIR", "LAUROUTPUT",
  "UNAFFECTED", "POSITION", "FREQUENCY", "PROBAND", "STRING", "VARIABLE",
  "ASSIGN", "ARRAY_VAR", "INTEGER", "SYSTEM_VAR", "REAL", "'='", "'<'",
  "'>'", "'+'", "'-'", "'*'", "'.'", "'/'", "UMINUS", "'('", "')'", "','",
  "'['", "']'", "';'", "'x'", "'y'", "$accept", "comfile", "command",
  "assigncommand", "assignment", "expression", "setcommand", "docommand",
  "includecommand", "$@1", "censorcommand", "affectedcommand",
  "sexcommand", "whilecommand", "positioncommand", "condition",
  "printcommand", "printlist", "defformatcommand", "arraycommand",
  "arraylist", "usecommand", "wherecommand", "res_condition",
  "filecommand", "filename_string", "formatlist", "fformatlist", "fformat",
  "format", "logcommand", "outputcommand", "errordircommand",
  "missingcommand", "pedcommand", "locicommand", "changetypecommand",
  "linkcom1", "linkcommand", "filtercommand", "modelcommand",
  "groupcommand", "modellist", "interactionlist", "variable",
  "single_vlist", "single_element", "locuslist", "locus", "lociclause",
  "slocuslist", "slocus", "simple_varlist", "open_bracket", "loop_clause1",
  "loop_clause", "fsingle_vlist", "filevarlist", "varlist",
  "complex_string1", "complex_string", 0
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
     305,   306,   307,   308,   309,   310,   311,   312,   313,   314,
      61,    60,    62,    43,    45,    42,    46,    47,   315,    40,
      41,    44,    91,    93,    59,   120,   121
};
# endif

/* YYR1[YYN] -- Symbol number of symbol that rule YYN derives.  */
static const yytype_uint8 yyr1[] =
{
       0,    77,    78,    78,    78,    78,    78,    79,    79,    79,
      79,    79,    79,    79,    79,    79,    79,    79,    79,    79,
      79,    79,    79,    79,    79,    79,    79,    79,    79,    79,
      79,    79,    79,    80,    81,    81,    82,    82,    82,    82,
      82,    82,    82,    82,    82,    82,    83,    83,    84,    86,
      85,    87,    88,    88,    88,    89,    90,    91,    91,    91,
      91,    91,    92,    92,    92,    92,    92,    92,    92,    92,
      92,    92,    93,    94,    94,    95,    95,    95,    95,    96,
      97,    97,    98,    98,    99,    99,   100,   100,   100,   100,
     100,   100,   100,   100,   100,   100,   100,   100,   100,   100,
     100,   100,   100,   100,   100,   101,   101,   101,   102,   102,
     103,   103,   103,   103,   104,   104,   104,   105,   105,   105,
     105,   106,   106,   106,   106,   106,   107,   108,   108,   108,
     109,   110,   110,   110,   110,   111,   112,   112,   112,   113,
     113,   113,   113,   113,   113,   114,   114,   114,   115,   115,
     115,   116,   117,   117,   118,   119,   119,   119,   119,   120,
     120,   120,   120,   121,   121,   122,   122,   122,   123,   123,
     124,   124,   125,   125,   125,   126,   126,   126,   127,   127,
     128,   128,   128,   129,   129,   130,   131,   131,   132,   132,
     133,   133,   133,   133,   134,   134,   134,   134,   135,   135,
     135,   135,   136,   136,   136,   137,   137,   137,   137
};

/* YYR2[YYN] -- Number of symbols composing right hand side of rule YYN.  */
static const yytype_uint8 yyr2[] =
{
       0,     2,     1,     1,     2,     3,     3,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     3,     6,     1,     1,     1,     1,
       3,     3,     3,     3,     2,     3,     3,     3,     1,     0,
       3,     6,     5,     5,     5,     5,     4,     3,     5,     7,
       4,     5,     1,     3,     3,     3,     3,     3,     3,     3,
       3,     2,     2,     1,     3,     2,     2,     2,     2,     2,
       4,     6,     6,     5,     6,     4,     1,     1,     1,     1,
       3,     3,     3,     3,     3,     2,     3,     3,     3,     3,
       3,     3,     3,     3,     2,     4,     7,     7,     1,     4,
       1,     3,     4,     6,     1,     3,     3,     2,     2,     2,
       2,     1,     2,     1,     2,     1,     2,     2,     2,     2,
       2,     2,     3,     4,     5,     2,     3,     3,     3,     2,
       2,     2,     2,     2,     2,     1,     4,     4,     2,     4,
       2,     2,     4,     7,     2,     3,     3,     1,     1,     3,
       3,     3,     3,     1,     1,     1,     1,     4,     1,     4,
       1,     3,     2,     1,     1,     3,     5,     3,     1,     3,
       4,     1,     4,     1,     3,     1,     8,    10,     4,     1,
       0,     1,     1,     4,     1,     1,     3,     3,     1,     1,
       3,     3,     1,     3,     3,     1,     1,     3,     3
};

/* YYDEFACT[STATE-NAME] -- Default rule to reduce with in state
   STATE-NUM when YYTABLE doesn't specify something else to do.  Zero
   means the default is an error.  */
static const yytype_uint8 yydefact[] =
{
       0,     3,     0,     0,     0,     0,     0,     0,     0,     0,
       0,   145,     0,     0,    48,     0,     0,     0,     0,     0,
      49,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     2,    19,    33,    27,    22,     9,    24,
      25,    29,    23,    32,    21,    28,    20,    17,    18,     7,
      10,    30,    31,    16,    11,    12,    13,     0,    14,     8,
      15,    26,     0,   164,   205,   163,     0,     0,     0,   168,
     206,   108,     0,     0,   166,   185,   165,   198,     0,   189,
     199,   141,   135,   126,     0,     0,     0,   151,     0,    38,
      37,    36,     0,     0,     0,   131,    39,   142,     0,     0,
       0,     0,     0,    79,    73,    72,     0,   129,   139,   140,
      75,    76,    78,    77,     0,   154,     0,     0,     0,     0,
     127,   130,   128,     0,     0,     0,     0,     0,   144,   143,
       1,     0,   202,   166,   165,     0,   148,   150,     0,     0,
     125,     0,     0,     0,     0,     0,   123,     0,     0,   114,
     110,   190,     0,   174,   173,   136,   170,   137,     0,   183,
       0,     0,     0,   181,   138,   178,     0,     0,     0,     0,
      44,     0,     0,     0,     0,     0,     0,     0,   132,     0,
      62,     0,     0,     0,     0,    88,    86,    87,     0,     0,
       0,    89,     0,     0,     0,    50,     0,    46,    47,     0,
       0,     0,     0,    57,     0,    34,     0,     6,     5,     0,
       0,     0,     0,     0,     0,   118,   117,   120,   119,   124,
       0,   122,     0,     0,     0,     0,     0,   192,   191,   195,
     194,   105,   207,   208,     0,   172,     0,     0,     0,     0,
     200,   201,     0,     0,     0,   152,   157,   158,   146,   147,
      45,     0,    40,    41,    42,    43,   133,    71,     0,     0,
       0,     0,     0,     0,     0,     0,    56,     0,     0,   104,
      95,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,    85,     0,     0,    74,     0,     0,
       0,     0,     0,    60,     0,     0,     0,     0,   204,   203,
     149,   109,   169,     0,     0,   111,     0,   115,     0,   116,
       0,   190,     0,     0,   171,   167,     0,   184,   188,     0,
       0,     0,   179,     0,     0,     0,     0,     0,     0,   134,
      69,    70,    64,    65,    66,    63,    67,    68,    83,     0,
      94,   102,   103,    97,    98,    99,    96,   100,   101,    90,
      91,    92,    93,     0,    80,     0,     0,    55,    52,    53,
      61,    58,    54,     0,   167,   112,     0,   190,   190,     0,
     197,   196,   177,     0,   175,     0,   182,   180,     0,   156,
     155,   161,   162,   159,   160,    82,    84,     0,    51,     0,
      35,     0,   106,   107,   193,     0,     0,   153,    81,    59,
     113,   176,     0,   186,     0,     0,   187
};

/* YYDEFGOTO[NTERM-NUM].  */
static const yytype_int16 yydefgoto[] =
{
      -1,    42,    43,    44,    45,   190,    46,    47,    48,   116,
      49,    50,    51,    52,    53,   191,    54,   115,    55,    56,
     113,    57,    58,   200,    59,    78,   157,   158,   159,   160,
      60,    61,    62,    63,    64,    65,    66,    67,    68,    69,
      70,    71,   255,   256,    79,    87,   106,   165,   166,   245,
     174,   175,   170,    88,    89,    90,   240,   241,    91,   147,
      81
};

/* YYPACT[STATE-NUM] -- Index in YYTABLE of the portion describing
   STATE-NUM.  */
#define YYPACT_NINF -230
static const yytype_int16 yypact[] =
{
     561,  -230,   204,     6,    11,   345,   345,    -2,    25,   168,
      -2,    -3,   669,   345,  -230,    -9,   114,    14,    35,   780,
    -230,    -2,   345,   345,    -2,    -2,    16,    -2,   339,   339,
     -21,   339,    76,    -2,    -2,    -2,    80,   339,    85,   -10,
     345,   345,    19,  -230,  -230,  -230,  -230,  -230,  -230,  -230,
    -230,  -230,  -230,  -230,  -230,  -230,  -230,  -230,  -230,  -230,
    -230,  -230,  -230,  -230,  -230,  -230,  -230,   140,  -230,  -230,
    -230,  -230,    45,  -230,  -230,  -230,    66,   102,    73,  -230,
    -230,    94,   412,   345,    92,  -230,  -230,  -230,   587,   133,
    -230,    98,    98,    94,   339,   110,   125,    94,    78,  -230,
    -230,  -230,   780,   780,    -2,   680,  -230,    98,   267,   126,
     -18,   728,   129,   130,   597,   165,    -2,    94,    98,    98,
      94,    94,  -230,    94,   230,  -230,   181,   185,   780,   184,
      94,    94,    94,   192,   742,   194,   780,   780,    98,    98,
    -230,   502,  -230,   208,   212,   215,    98,   -14,    -2,   780,
    -230,    -2,    -2,   223,    -2,    24,  -230,   -42,    72,  -230,
    -230,   386,   218,    66,   217,   220,  -230,    98,   780,  -230,
     224,   587,   345,   227,   225,  -230,   780,   587,   221,   229,
    -230,   265,   -35,   780,   780,   780,   780,   345,    98,   267,
     597,   127,   728,   234,   728,  -230,  -230,  -230,   728,   728,
     400,  -230,   780,   250,   780,    94,   236,  -230,  -230,   712,
     728,   728,   761,   745,   728,   597,   348,  -230,  -230,   780,
     586,   256,   345,    -7,   460,    94,    94,  -230,    94,  -230,
       5,  -230,    31,   375,    81,   375,    81,   243,  -230,  -230,
    -230,   245,  -230,  -230,    23,  -230,   412,   792,    -8,   -30,
    -230,  -230,    54,   412,   800,   264,   145,   149,  -230,  -230,
    -230,   780,     3,     3,  -230,  -230,    98,  -230,   267,   267,
     267,   267,   267,   267,   267,   267,  -230,   411,   728,  -230,
    -230,   562,   728,   728,   728,   728,   728,   728,   728,   728,
     728,   728,   728,   728,   300,   808,   268,   597,   728,   780,
     612,   623,   780,   597,   780,   634,   274,   816,  -230,  -230,
      98,  -230,  -230,    -6,    42,  -230,   277,  -230,   284,  -230,
     780,   386,   283,    65,  -230,  -230,   304,  -230,  -230,   587,
     291,    89,  -230,   303,   587,   587,   587,   587,   587,   597,
     222,   290,   193,  -230,  -230,   193,  -230,  -230,  -230,   645,
    -230,   384,   706,   282,   726,   726,   282,   726,   726,   115,
     115,  -230,  -230,   345,  -230,   780,   695,   597,  -230,  -230,
     597,   764,  -230,   780,   321,  -230,     5,   386,   386,   824,
    -230,  -230,  -230,   339,  -230,   295,  -230,  -230,   587,   145,
     149,  -230,  -230,  -230,  -230,  -230,    98,   832,  -230,   780,
     597,   158,   245,   245,  -230,   313,   780,   264,  -230,   597,
    -230,  -230,   783,  -230,   780,   840,  -230
};

/* YYPGOTO[NTERM-NUM].  */
static const yytype_int16 yypgoto[] =
{
    -230,  -230,   247,  -230,    64,   -12,  -230,  -230,  -230,  -230,
    -230,  -230,  -230,  -230,  -230,  -142,  -230,  -230,  -230,  -230,
    -230,  -230,  -230,    87,  -230,   -29,  -229,  -230,   -11,   159,
    -230,  -230,  -230,  -230,  -230,  -230,  -230,  -230,  -230,  -230,
    -230,  -230,     8,    63,     4,   -70,   206,  -230,  -172,  -230,
    -230,  -230,   239,  -230,  -230,  -158,   106,  -108,    -1,   233,
      51
};

/* YYTABLE[YYPACT[STATE-NUM]].  What to do in state STATE-NUM.  If
   positive, shift that token.  If negative, reduce the rule which
   number is the opposite.  If zero, do what YYDEFACT says.
   If YYTABLE_NINF, syntax error.  */
#define YYTABLE_NINF -170
static const yytype_int16 yytable[] =
{
     105,   313,   193,   239,   328,    92,   150,   114,   326,    86,
      86,    82,   107,    96,   251,   110,    83,    86,   169,   140,
      86,   118,   119,    73,   322,   229,    86,    86,   162,   232,
      94,   233,   150,    75,   127,   141,    73,   126,   261,   138,
     139,   329,    73,   229,    86,    86,    75,   267,    84,   221,
     136,    74,    75,   172,    76,   330,   162,   222,    93,   137,
     108,    97,   155,   311,   375,   232,   146,    73,   185,    98,
     186,   144,   117,   122,   324,   120,   121,    75,   123,    76,
     156,   332,   167,   111,   130,   131,   132,    86,   314,   112,
     180,   181,    86,   230,  -121,  -121,   129,  -121,    73,   231,
     133,   169,   250,   150,   188,   135,   156,   257,    75,    86,
      84,   376,  -121,  -121,   148,  -121,   209,   231,   151,   152,
     153,   154,   213,    85,   215,   216,   340,   341,   342,   343,
     344,   345,   346,   347,   109,   149,   383,   224,   384,   151,
     152,   153,   154,   234,   161,   235,   236,   401,   268,   269,
     270,   271,   272,   178,   179,   182,   247,   162,    73,   155,
     172,   168,   387,   380,   254,   238,   171,   205,    75,   172,
      84,   262,   263,   264,   265,    86,    86,   156,   327,   176,
     292,    86,   293,    85,    73,   177,   266,   273,   274,   275,
     295,    86,   297,   142,    75,   192,   143,   276,   202,   223,
     303,   203,   225,   226,   316,   228,   318,   307,    80,    85,
     335,   336,    73,    80,   337,   338,    80,   271,   272,   239,
     239,   310,    75,   317,    95,   319,    86,    80,   410,   232,
      80,    80,    72,    80,   124,   125,   204,   128,   207,    80,
      80,    80,   208,   134,   269,   270,   271,   272,    73,   339,
     206,   331,    86,   210,   274,   275,    86,    74,    75,   327,
      76,   211,    73,   214,   390,   391,   392,   393,   394,   402,
     403,   242,    75,   145,    76,  -168,    77,   219,   220,   277,
     227,   279,   273,   274,   275,   280,   281,   367,   164,   244,
     370,   246,   371,   189,   258,   248,   253,   300,   301,   252,
     173,   305,   259,   278,   296,   298,   285,   286,   379,   309,
      80,    73,   320,   270,   271,   272,   321,   201,   257,   363,
      99,    75,    80,    76,   100,   238,   101,   334,   183,   184,
     185,   102,   186,    86,   373,   260,   103,   365,    86,    86,
      86,    86,    86,   288,   289,   290,   291,   292,   377,   293,
     273,   274,   275,   397,    80,   378,   382,    80,    80,    39,
      80,   400,   396,   388,   386,   349,   406,    86,   243,   351,
     352,   353,   354,   355,   356,   357,   358,   359,   360,   361,
     362,   238,   238,    73,  -169,   366,   411,   409,   218,    73,
     385,   315,    86,    75,   412,    76,   407,   389,   201,    75,
     201,    84,   415,    72,   201,   201,   283,   284,   285,   286,
     249,   183,   184,   185,    85,   186,   201,   201,   306,    73,
     201,   282,   283,   284,   285,   286,   145,   381,    74,    75,
      73,    76,   282,   283,   284,   285,   286,     0,     0,    80,
      75,    80,   237,     0,   287,   288,   289,   290,   291,   292,
     323,   293,   164,   308,     0,    85,    73,     0,     0,   164,
     287,   288,   289,   290,   291,   292,    75,   293,   163,     0,
     294,   287,   288,   289,   290,   291,   292,     0,   293,     0,
       0,   348,     0,     0,   201,     0,     0,     0,   201,   201,
     201,   201,   201,   201,   201,   201,   201,   201,   201,   201,
       0,     0,    -4,   217,   201,     2,     3,     0,     4,     5,
       6,     7,     8,     9,    10,    11,    12,    13,    -4,    14,
      15,    16,    17,   183,   184,   185,     0,   186,     0,     0,
     312,    18,    19,    20,    21,     0,     0,    22,    23,    24,
      25,    26,    27,    28,    29,    30,    31,    32,    33,    34,
      35,    36,    37,     0,    38,     0,     0,    39,     0,    40,
       0,    41,     1,     0,     2,     3,     0,     4,     5,     6,
       7,     8,     9,    10,    11,    12,    13,     0,    14,    15,
      16,    17,     0,   282,   283,   284,   285,   286,     0,   405,
      18,    19,    20,    21,     0,     0,    22,    23,    24,    25,
      26,    27,    28,    29,    30,    31,    32,    33,    34,    35,
      36,    37,     0,    38,     0,     0,    39,     0,    40,     0,
      41,     0,   287,   288,   289,   290,   291,   292,     0,   293,
      73,    73,   350,   282,   283,   284,   285,   286,     0,   142,
      75,    75,    76,    84,   282,   283,   284,   285,   286,     0,
       0,     0,     0,     0,     0,   282,   283,   284,   285,   286,
     183,   184,   185,     0,   186,     0,   282,   283,   284,   285,
     286,     0,   287,   288,   289,   290,   291,   292,     0,   293,
       0,     0,   368,   287,   288,   289,   290,   291,   292,     0,
     293,     0,     0,   369,   287,   288,   289,   290,   291,   292,
       0,   293,     0,     0,   372,   287,   288,   289,   290,   291,
     292,     0,   293,    73,     0,   395,   282,   283,   284,   285,
     286,     0,    99,    75,    73,    76,   100,     0,   101,   284,
     285,   286,     0,   102,    75,     0,    84,     0,   103,     0,
       0,   104,     0,   183,   184,   185,     0,   186,     0,    85,
       0,   187,     0,     0,   194,   287,   288,   289,   290,   291,
     292,     0,   293,     0,     0,   398,   287,   288,   289,   290,
     291,   292,    73,   293,     0,   183,   184,   185,     0,   186,
       0,   195,    75,   299,    76,   196,    73,   197,     0,   290,
     291,   292,   198,   293,     0,    99,    75,   199,    76,   100,
       0,   101,     0,     0,     0,    73,   102,     0,   183,   184,
     185,   103,   186,   212,    99,    75,   304,    76,   100,     0,
     101,     0,     0,     0,    73,   102,     0,   183,   184,   185,
     103,   186,   302,    99,    75,   399,    76,   100,     0,   101,
       0,     0,     0,     0,   102,     0,   183,   184,   185,   103,
     186,     0,     0,   413,   414,   183,   184,   185,     0,   186,
       0,     0,   325,   183,   184,   185,     0,   186,     0,     0,
     333,   183,   184,   185,     0,   186,     0,     0,   364,   183,
     184,   185,     0,   186,     0,     0,   374,   183,   184,   185,
       0,   186,     0,     0,   404,   183,   184,   185,     0,   186,
       0,     0,   408,   183,   184,   185,     0,   186,     0,     0,
     416
};

static const yytype_int16 yycheck[] =
{
      12,   230,    20,   161,    34,     6,     1,    19,    16,     5,
       6,     5,    13,     9,   172,    16,     5,    13,    88,     0,
      16,    22,    23,    44,     1,     1,    22,    23,    63,    71,
       5,    73,     1,    54,    30,    16,    44,    58,    73,    40,
      41,    71,    44,     1,    40,    41,    54,   189,    56,    63,
      60,    53,    54,    71,    56,     1,    63,    71,     7,    69,
      69,    10,    57,    70,    70,    71,    67,    44,    65,    72,
      67,    67,    21,    57,   246,    24,    25,    54,    27,    56,
      75,   253,    83,    69,    33,    34,    35,    83,    57,    54,
     102,   103,    88,    69,    70,    71,    20,    73,    44,    75,
      20,   171,   172,     1,   105,    20,    75,   177,    54,   105,
      56,    69,    70,    71,    69,    73,   128,    75,    37,    38,
      39,    40,   134,    69,   136,   137,   268,   269,   270,   271,
     272,   273,   274,   275,    20,    69,    71,   149,    73,    37,
      38,    39,    40,    71,    71,    73,    74,   376,    21,    22,
      23,    24,    25,    75,    76,   104,   168,    63,    44,    57,
      71,    69,    73,   321,   176,   161,    33,   116,    54,    71,
      56,   183,   184,   185,   186,   171,   172,    75,   248,    69,
      65,   177,    67,    69,    44,    60,   187,    60,    61,    62,
     202,   187,   204,    53,    54,    69,    56,    70,    69,   148,
     212,    71,   151,   152,   233,   154,   235,   219,     2,    69,
      65,    66,    44,     7,    65,    66,    10,    24,    25,   377,
     378,   222,    54,   234,    56,   236,   222,    21,    70,    71,
      24,    25,    28,    27,    28,    29,    71,    31,    57,    33,
      34,    35,    57,    37,    22,    23,    24,    25,    44,   261,
      20,   252,   248,    69,    61,    62,   252,    53,    54,   329,
      56,    69,    44,    69,   334,   335,   336,   337,   338,   377,
     378,    53,    54,    67,    56,    63,    72,    69,    63,   192,
      57,   194,    60,    61,    62,   198,   199,   299,    82,    72,
     302,    71,   304,    26,    73,    71,    71,   210,   211,    72,
      94,   214,    73,    69,    54,    69,    24,    25,   320,    53,
     104,    44,    69,    23,    24,    25,    71,   111,   388,    19,
      53,    54,   116,    56,    57,   321,    59,    63,    63,    64,
      65,    64,    67,   329,    60,    70,    69,    69,   334,   335,
     336,   337,   338,    61,    62,    63,    64,    65,    71,    67,
      60,    61,    62,   365,   148,    71,    73,   151,   152,    55,
     154,   373,   363,    60,    73,   278,    71,   363,   162,   282,
     283,   284,   285,   286,   287,   288,   289,   290,   291,   292,
     293,   377,   378,    44,    63,   298,    73,   399,   141,    44,
     326,   232,   388,    54,   406,    56,   388,   334,   192,    54,
     194,    56,   414,    28,   198,   199,    22,    23,    24,    25,
     171,    63,    64,    65,    69,    67,   210,   211,    70,    44,
     214,    21,    22,    23,    24,    25,   220,   321,    53,    54,
      44,    56,    21,    22,    23,    24,    25,    -1,    -1,   233,
      54,   235,    56,    -1,    60,    61,    62,    63,    64,    65,
     244,    67,   246,   220,    -1,    69,    44,    -1,    -1,   253,
      60,    61,    62,    63,    64,    65,    54,    67,    56,    -1,
      70,    60,    61,    62,    63,    64,    65,    -1,    67,    -1,
      -1,    70,    -1,    -1,   278,    -1,    -1,    -1,   282,   283,
     284,   285,   286,   287,   288,   289,   290,   291,   292,   293,
      -1,    -1,     0,     1,   298,     3,     4,    -1,     6,     7,
       8,     9,    10,    11,    12,    13,    14,    15,    16,    17,
      18,    19,    20,    63,    64,    65,    -1,    67,    -1,    -1,
      70,    29,    30,    31,    32,    -1,    -1,    35,    36,    37,
      38,    39,    40,    41,    42,    43,    44,    45,    46,    47,
      48,    49,    50,    -1,    52,    -1,    -1,    55,    -1,    57,
      -1,    59,     1,    -1,     3,     4,    -1,     6,     7,     8,
       9,    10,    11,    12,    13,    14,    15,    -1,    17,    18,
      19,    20,    -1,    21,    22,    23,    24,    25,    -1,   383,
      29,    30,    31,    32,    -1,    -1,    35,    36,    37,    38,
      39,    40,    41,    42,    43,    44,    45,    46,    47,    48,
      49,    50,    -1,    52,    -1,    -1,    55,    -1,    57,    -1,
      59,    -1,    60,    61,    62,    63,    64,    65,    -1,    67,
      44,    44,    70,    21,    22,    23,    24,    25,    -1,    53,
      54,    54,    56,    56,    21,    22,    23,    24,    25,    -1,
      -1,    -1,    -1,    -1,    -1,    21,    22,    23,    24,    25,
      63,    64,    65,    -1,    67,    -1,    21,    22,    23,    24,
      25,    -1,    60,    61,    62,    63,    64,    65,    -1,    67,
      -1,    -1,    70,    60,    61,    62,    63,    64,    65,    -1,
      67,    -1,    -1,    70,    60,    61,    62,    63,    64,    65,
      -1,    67,    -1,    -1,    70,    60,    61,    62,    63,    64,
      65,    -1,    67,    44,    -1,    70,    21,    22,    23,    24,
      25,    -1,    53,    54,    44,    56,    57,    -1,    59,    23,
      24,    25,    -1,    64,    54,    -1,    56,    -1,    69,    -1,
      -1,    72,    -1,    63,    64,    65,    -1,    67,    -1,    69,
      -1,    71,    -1,    -1,    26,    60,    61,    62,    63,    64,
      65,    -1,    67,    -1,    -1,    70,    60,    61,    62,    63,
      64,    65,    44,    67,    -1,    63,    64,    65,    -1,    67,
      -1,    53,    54,    71,    56,    57,    44,    59,    -1,    63,
      64,    65,    64,    67,    -1,    53,    54,    69,    56,    57,
      -1,    59,    -1,    -1,    -1,    44,    64,    -1,    63,    64,
      65,    69,    67,    71,    53,    54,    71,    56,    57,    -1,
      59,    -1,    -1,    -1,    44,    64,    -1,    63,    64,    65,
      69,    67,    71,    53,    54,    71,    56,    57,    -1,    59,
      -1,    -1,    -1,    -1,    64,    -1,    63,    64,    65,    69,
      67,    -1,    -1,    70,    71,    63,    64,    65,    -1,    67,
      -1,    -1,    70,    63,    64,    65,    -1,    67,    -1,    -1,
      70,    63,    64,    65,    -1,    67,    -1,    -1,    70,    63,
      64,    65,    -1,    67,    -1,    -1,    70,    63,    64,    65,
      -1,    67,    -1,    -1,    70,    63,    64,    65,    -1,    67,
      -1,    -1,    70,    63,    64,    65,    -1,    67,    -1,    -1,
      70
};

/* YYSTOS[STATE-NUM] -- The (internal number of the) accessing
   symbol of state STATE-NUM.  */
static const yytype_uint8 yystos[] =
{
       0,     1,     3,     4,     6,     7,     8,     9,    10,    11,
      12,    13,    14,    15,    17,    18,    19,    20,    29,    30,
      31,    32,    35,    36,    37,    38,    39,    40,    41,    42,
      43,    44,    45,    46,    47,    48,    49,    50,    52,    55,
      57,    59,    78,    79,    80,    81,    83,    84,    85,    87,
      88,    89,    90,    91,    93,    95,    96,    98,    99,   101,
     107,   108,   109,   110,   111,   112,   113,   114,   115,   116,
     117,   118,    28,    44,    53,    54,    56,    72,   102,   121,
     123,   137,     5,     5,    56,    69,   121,   122,   130,   131,
     132,   135,   135,   137,     5,    56,   121,   137,    72,    53,
      57,    59,    64,    69,    72,    82,   123,   135,    69,    20,
     135,    69,    54,    97,    82,    94,    86,   137,   135,   135,
     137,   137,    57,   137,   123,   123,    58,   121,   123,    20,
     137,   137,   137,    20,   123,    20,    60,    69,   135,   135,
       0,    16,    53,    56,   121,   123,   135,   136,    69,    69,
       1,    37,    38,    39,    40,    57,    75,   103,   104,   105,
     106,    71,    63,    56,   123,   124,   125,   135,    69,   122,
     129,    33,    71,   123,   127,   128,    69,    60,    75,    76,
      82,    82,   137,    63,    64,    65,    67,    71,   135,    26,
      82,    92,    69,    20,    26,    53,    57,    59,    64,    69,
     100,   123,    69,    71,    71,   137,    20,    57,    57,    82,
      69,    69,    71,    82,    69,    82,    82,     1,    79,    69,
      63,    63,    71,   137,    82,   137,   137,    57,   137,     1,
      69,    75,    71,    73,    71,    73,    74,    56,   121,   132,
     133,   134,    53,   123,    72,   126,    71,    82,    71,   129,
     122,   132,    72,    71,    82,   119,   120,   122,    73,    73,
      70,    73,    82,    82,    82,    82,   135,    92,    21,    22,
      23,    24,    25,    60,    61,    62,    70,   100,    69,   100,
     100,   100,    21,    22,    23,    24,    25,    60,    61,    62,
      63,    64,    65,    67,    70,    82,    54,    82,    69,    71,
     100,   100,    71,    82,    71,   100,    70,    82,   136,    53,
     135,    70,    70,   103,    57,   106,   102,   105,   102,   105,
      69,    71,     1,   123,   125,    70,    16,   122,    34,    71,
       1,   135,   125,    70,    63,    65,    66,    65,    66,    82,
      92,    92,    92,    92,    92,    92,    92,    92,    70,   100,
      70,   100,   100,   100,   100,   100,   100,   100,   100,   100,
     100,   100,   100,    19,    70,    69,   100,    82,    70,    70,
      82,    82,    70,    60,    70,    70,    69,    71,    71,    82,
     132,   133,    73,    71,    73,    81,    73,    73,    60,   120,
     122,   122,   122,   122,   122,    70,   135,    82,    70,    71,
      82,   103,   134,   134,    70,   123,    71,   119,    70,    82,
      70,    73,    82,    70,    71,    82,    70
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
#line 167 "control_parse.y"
    { new_command(); }
    break;

  case 3:

/* Line 1455 of yacc.c  */
#line 168 "control_parse.y"
    { new_command(); }
    break;

  case 5:

/* Line 1455 of yacc.c  */
#line 170 "control_parse.y"
    { new_command(); }
    break;

  case 6:

/* Line 1455 of yacc.c  */
#line 171 "control_parse.y"
    { new_command(); }
    break;

  case 33:

/* Line 1455 of yacc.c  */
#line 202 "control_parse.y"
    {}
    break;

  case 34:

/* Line 1455 of yacc.c  */
#line 205 "control_parse.y"
    { (yyval.element)=assign_var((yyvsp[(1) - (3)].var),0,(yyvsp[(3) - (3)].express)); if((yyvsp[(3) - (3)].express)) free((yyvsp[(3) - (3)].express));}
    break;

  case 35:

/* Line 1455 of yacc.c  */
#line 206 "control_parse.y"
    { (yyval.element)=assign_var((yyvsp[(1) - (6)].var),(yyvsp[(3) - (6)].express),(yyvsp[(6) - (6)].express)); if((yyvsp[(3) - (6)].express)) free((yyvsp[(3) - (6)].express)); if((yyvsp[(6) - (6)].express)) free((yyvsp[(6) - (6)].express));}
    break;

  case 36:

/* Line 1455 of yacc.c  */
#line 209 "control_parse.y"
    { (yyval.express)=alloc_express(); (yyval.express)->type=ST_REAL; (yyval.express)->arg.rvalue=(yyvsp[(1) - (1)].rvalue); }
    break;

  case 37:

/* Line 1455 of yacc.c  */
#line 210 "control_parse.y"
    { (yyval.express)=alloc_express(); (yyval.express)->type=ST_INTEGER; (yyval.express)->arg.value=(yyvsp[(1) - (1)].value); }
    break;

  case 38:

/* Line 1455 of yacc.c  */
#line 211 "control_parse.y"
    { (yyval.express)=alloc_express(); (yyval.express)->type=ST_STRING; (yyval.express)->arg.string=(yyvsp[(1) - (1)].string); }
    break;

  case 39:

/* Line 1455 of yacc.c  */
#line 212 "control_parse.y"
    { (yyval.express)=alloc_express(); (yyval.express)->type=(yyvsp[(1) - (1)].element)->type; if((yyvsp[(1) - (1)].element)->type==ST_STRING) (yyval.express)->arg.string=string_copy(0,(yyvsp[(1) - (1)].element)->arg.string);
		  else (yyval.express)->arg=(yyvsp[(1) - (1)].element)->arg; }
    break;

  case 40:

/* Line 1455 of yacc.c  */
#line 214 "control_parse.y"
    { (yyval.express)=do_express_op((yyvsp[(1) - (3)].express),(yyvsp[(3) - (3)].express),'+'); }
    break;

  case 41:

/* Line 1455 of yacc.c  */
#line 215 "control_parse.y"
    { (yyval.express)=do_express_op((yyvsp[(1) - (3)].express),(yyvsp[(3) - (3)].express),'-'); }
    break;

  case 42:

/* Line 1455 of yacc.c  */
#line 216 "control_parse.y"
    { (yyval.express)=do_express_op((yyvsp[(1) - (3)].express),(yyvsp[(3) - (3)].express),'*'); }
    break;

  case 43:

/* Line 1455 of yacc.c  */
#line 217 "control_parse.y"
    { (yyval.express)=do_express_op((yyvsp[(1) - (3)].express),(yyvsp[(3) - (3)].express),'/'); }
    break;

  case 44:

/* Line 1455 of yacc.c  */
#line 218 "control_parse.y"
    { (yyval.express)=do_express_op((yyvsp[(2) - (2)].express),0,'-'); }
    break;

  case 45:

/* Line 1455 of yacc.c  */
#line 219 "control_parse.y"
    { (yyval.express)=(yyvsp[(2) - (3)].express); }
    break;

  case 46:

/* Line 1455 of yacc.c  */
#line 222 "control_parse.y"
    { syst_var[(yyvsp[(2) - (3)].value)]=(yyvsp[(3) - (3)].value); }
    break;

  case 47:

/* Line 1455 of yacc.c  */
#line 223 "control_parse.y"
    { yyerror("Unrecognized system variable"); }
    break;

  case 48:

/* Line 1455 of yacc.c  */
#line 226 "control_parse.y"
    { enter_loop(); }
    break;

  case 49:

/* Line 1455 of yacc.c  */
#line 229 "control_parse.y"
    {iflag=1;}
    break;

  case 50:

/* Line 1455 of yacc.c  */
#line 229 "control_parse.y"
    {include_control_file((yyvsp[(3) - (3)].string));}
    break;

  case 51:

/* Line 1455 of yacc.c  */
#line 232 "control_parse.y"
    { add_censored((yyvsp[(2) - (6)].element),1); at_use=0;}
    break;

  case 52:

/* Line 1455 of yacc.c  */
#line 235 "control_parse.y"
    { add_censored(0,0); at_use=0;}
    break;

  case 53:

/* Line 1455 of yacc.c  */
#line 236 "control_parse.y"
    { add_censored(0,2); at_use=0;}
    break;

  case 54:

/* Line 1455 of yacc.c  */
#line 237 "control_parse.y"
    { add_censored(0,3); at_use=0;}
    break;

  case 55:

/* Line 1455 of yacc.c  */
#line 240 "control_parse.y"
    { set_sex((yyvsp[(2) - (5)].element),(yyvsp[(3) - (5)].express),(yyvsp[(5) - (5)].express)); }
    break;

  case 56:

/* Line 1455 of yacc.c  */
#line 243 "control_parse.y"
    { do_while_com((yyvsp[(3) - (4)].express)); if((yyvsp[(3) - (4)].express)) free((yyvsp[(3) - (4)].express));}
    break;

  case 57:

/* Line 1455 of yacc.c  */
#line 246 "control_parse.y"
    { set_position((yyvsp[(2) - (3)].element),(yyvsp[(3) - (3)].express),0,0); }
    break;

  case 58:

/* Line 1455 of yacc.c  */
#line 247 "control_parse.y"
    { set_position((yyvsp[(2) - (5)].element),0,(yyvsp[(3) - (5)].express),(yyvsp[(5) - (5)].express)); }
    break;

  case 59:

/* Line 1455 of yacc.c  */
#line 248 "control_parse.y"
    { set_position((yyvsp[(2) - (7)].element),(yyvsp[(3) - (7)].express),(yyvsp[(5) - (7)].express),(yyvsp[(7) - (7)].express)); }
    break;

  case 60:

/* Line 1455 of yacc.c  */
#line 249 "control_parse.y"
    { set_position((yyvsp[(2) - (4)].element),0,0,(yyvsp[(4) - (4)].express)); }
    break;

  case 61:

/* Line 1455 of yacc.c  */
#line 250 "control_parse.y"
    { set_position((yyvsp[(2) - (5)].element),0,0,(yyvsp[(5) - (5)].express)); }
    break;

  case 63:

/* Line 1455 of yacc.c  */
#line 254 "control_parse.y"
    {(yyval.express)=do_logical_op((yyvsp[(1) - (3)].express),(yyvsp[(3) - (3)].express),'=');}
    break;

  case 64:

/* Line 1455 of yacc.c  */
#line 255 "control_parse.y"
    {(yyval.express)=do_logical_op((yyvsp[(1) - (3)].express),(yyvsp[(3) - (3)].express),NEQSYMBOL);}
    break;

  case 65:

/* Line 1455 of yacc.c  */
#line 256 "control_parse.y"
    {(yyval.express)=do_logical_op((yyvsp[(1) - (3)].express),(yyvsp[(3) - (3)].express),LEQSYMBOL);}
    break;

  case 66:

/* Line 1455 of yacc.c  */
#line 257 "control_parse.y"
    {(yyval.express)=do_logical_op((yyvsp[(1) - (3)].express),(yyvsp[(3) - (3)].express),GEQSYMBOL);}
    break;

  case 67:

/* Line 1455 of yacc.c  */
#line 258 "control_parse.y"
    {(yyval.express)=do_logical_op((yyvsp[(1) - (3)].express),(yyvsp[(3) - (3)].express),'<');}
    break;

  case 68:

/* Line 1455 of yacc.c  */
#line 259 "control_parse.y"
    {(yyval.express)=do_logical_op((yyvsp[(1) - (3)].express),(yyvsp[(3) - (3)].express),'>');}
    break;

  case 69:

/* Line 1455 of yacc.c  */
#line 260 "control_parse.y"
    {(yyval.express)=do_logical_op((yyvsp[(1) - (3)].express),(yyvsp[(3) - (3)].express),ORSYMBOL);}
    break;

  case 70:

/* Line 1455 of yacc.c  */
#line 261 "control_parse.y"
    {(yyval.express)=do_logical_op((yyvsp[(1) - (3)].express),(yyvsp[(3) - (3)].express),ANDSYMBOL);}
    break;

  case 71:

/* Line 1455 of yacc.c  */
#line 262 "control_parse.y"
    {(yyval.express)=do_logical_op((yyvsp[(2) - (2)].express),0,NOTSYMBOL);}
    break;

  case 72:

/* Line 1455 of yacc.c  */
#line 265 "control_parse.y"
    { (void)fputc('\n',stdout); }
    break;

  case 73:

/* Line 1455 of yacc.c  */
#line 268 "control_parse.y"
    { print_exp((yyvsp[(1) - (1)].express)); if((yyvsp[(1) - (1)].express)) free((yyvsp[(1) - (1)].express)); }
    break;

  case 74:

/* Line 1455 of yacc.c  */
#line 269 "control_parse.y"
    { print_exp((yyvsp[(3) - (3)].express)); if((yyvsp[(3) - (3)].express)) free((yyvsp[(3) - (3)].express)); }
    break;

  case 75:

/* Line 1455 of yacc.c  */
#line 272 "control_parse.y"
    {if(rsformat) free(rsformat); rsformat=(yyvsp[(2) - (2)].string);}
    break;

  case 76:

/* Line 1455 of yacc.c  */
#line 273 "control_parse.y"
    {if(fsformat) free(fsformat); fsformat=(yyvsp[(2) - (2)].string);}
    break;

  case 77:

/* Line 1455 of yacc.c  */
#line 274 "control_parse.y"
    {if(gsformat) free(gsformat); gsformat=(yyvsp[(2) - (2)].string);}
    break;

  case 78:

/* Line 1455 of yacc.c  */
#line 275 "control_parse.y"
    {file_skip=(yyvsp[(2) - (2)].value);}
    break;

  case 80:

/* Line 1455 of yacc.c  */
#line 280 "control_parse.y"
    {set_array_var((yyvsp[(1) - (4)].var)->data,(yyvsp[(3) - (4)].express)); if((yyvsp[(3) - (4)].express)) free((yyvsp[(3) - (4)].express)); }
    break;

  case 81:

/* Line 1455 of yacc.c  */
#line 281 "control_parse.y"
    {set_array_var((yyvsp[(3) - (6)].var)->data,(yyvsp[(5) - (6)].express)); if((yyvsp[(5) - (6)].express)) free((yyvsp[(5) - (6)].express)); }
    break;

  case 82:

/* Line 1455 of yacc.c  */
#line 284 "control_parse.y"
    {add_restriction((yyvsp[(2) - (6)].var_list));at_use=0;}
    break;

  case 83:

/* Line 1455 of yacc.c  */
#line 285 "control_parse.y"
    {add_restriction(0);at_use=0;}
    break;

  case 84:

/* Line 1455 of yacc.c  */
#line 288 "control_parse.y"
    {add_restriction((yyvsp[(6) - (6)].var_list));at_use=0;}
    break;

  case 85:

/* Line 1455 of yacc.c  */
#line 289 "control_parse.y"
    {add_restriction(0);at_use=0;}
    break;

  case 86:

/* Line 1455 of yacc.c  */
#line 292 "control_parse.y"
    {add_operation(&((yyvsp[(1) - (1)].value)),INTEGER,0);}
    break;

  case 87:

/* Line 1455 of yacc.c  */
#line 293 "control_parse.y"
    {add_operation(&((yyvsp[(1) - (1)].rvalue)),REAL,0);}
    break;

  case 88:

/* Line 1455 of yacc.c  */
#line 294 "control_parse.y"
    {add_operation((yyvsp[(1) - (1)].string),STRING,0);}
    break;

  case 89:

/* Line 1455 of yacc.c  */
#line 295 "control_parse.y"
    {if((yyvsp[(1) - (1)].element)) check_element_add_op((yyvsp[(1) - (1)].element));}
    break;

  case 90:

/* Line 1455 of yacc.c  */
#line 296 "control_parse.y"
    {add_operation(0,0,'+');}
    break;

  case 91:

/* Line 1455 of yacc.c  */
#line 297 "control_parse.y"
    {add_operation(0,0,'-');}
    break;

  case 92:

/* Line 1455 of yacc.c  */
#line 298 "control_parse.y"
    {add_operation(0,0,'*');}
    break;

  case 93:

/* Line 1455 of yacc.c  */
#line 299 "control_parse.y"
    {add_operation(0,0,'/');}
    break;

  case 95:

/* Line 1455 of yacc.c  */
#line 301 "control_parse.y"
    {add_operation(0,0,UMINUS);}
    break;

  case 96:

/* Line 1455 of yacc.c  */
#line 302 "control_parse.y"
    {add_operation(0,0,'=');}
    break;

  case 97:

/* Line 1455 of yacc.c  */
#line 303 "control_parse.y"
    {add_operation(0,0,NEQSYMBOL);}
    break;

  case 98:

/* Line 1455 of yacc.c  */
#line 304 "control_parse.y"
    {add_operation(0,0,LEQSYMBOL);}
    break;

  case 99:

/* Line 1455 of yacc.c  */
#line 305 "control_parse.y"
    {add_operation(0,0,GEQSYMBOL);}
    break;

  case 100:

/* Line 1455 of yacc.c  */
#line 306 "control_parse.y"
    {add_operation(0,0,'<');}
    break;

  case 101:

/* Line 1455 of yacc.c  */
#line 307 "control_parse.y"
    {add_operation(0,0,'>');}
    break;

  case 102:

/* Line 1455 of yacc.c  */
#line 308 "control_parse.y"
    {add_operation(0,0,ORSYMBOL);}
    break;

  case 103:

/* Line 1455 of yacc.c  */
#line 309 "control_parse.y"
    {add_operation(0,0,ANDSYMBOL);}
    break;

  case 104:

/* Line 1455 of yacc.c  */
#line 310 "control_parse.y"
    {add_operation(0,0,NOTSYMBOL);}
    break;

  case 105:

/* Line 1455 of yacc.c  */
#line 313 "control_parse.y"
    { do_file_com((yyvsp[(2) - (4)].string),0,0,(yyvsp[(4) - (4)].var_list)); }
    break;

  case 106:

/* Line 1455 of yacc.c  */
#line 314 "control_parse.y"
    { do_file_com((yyvsp[(5) - (7)].string),(yyvsp[(3) - (7)].format_clause),0,(yyvsp[(7) - (7)].var_list)); }
    break;

  case 107:

/* Line 1455 of yacc.c  */
#line 315 "control_parse.y"
    { do_file_com((yyvsp[(5) - (7)].string),0,(yyvsp[(3) - (7)].fformat),(yyvsp[(7) - (7)].var_list)); }
    break;

  case 108:

/* Line 1455 of yacc.c  */
#line 318 "control_parse.y"
    { (yyval.string)=(yyvsp[(1) - (1)].string); }
    break;

  case 109:

/* Line 1455 of yacc.c  */
#line 319 "control_parse.y"
    { (yyval.string)=(yyvsp[(3) - (4)].string); shell_flag=1; }
    break;

  case 110:

/* Line 1455 of yacc.c  */
#line 322 "control_parse.y"
    {(yyval.format_clause)=add_f_atom(0,(yyvsp[(1) - (1)].f_atom)); }
    break;

  case 111:

/* Line 1455 of yacc.c  */
#line 323 "control_parse.y"
    {(yyval.format_clause)=add_f_atom((yyvsp[(1) - (3)].format_clause),(yyvsp[(3) - (3)].f_atom)); }
    break;

  case 112:

/* Line 1455 of yacc.c  */
#line 324 "control_parse.y"
    {(yyval.format_clause)=add_f_list(0,(yyvsp[(3) - (4)].format_clause),(yyvsp[(1) - (4)].value)); }
    break;

  case 113:

/* Line 1455 of yacc.c  */
#line 325 "control_parse.y"
    {(yyval.format_clause)=add_f_list((yyvsp[(1) - (6)].format_clause),(yyvsp[(5) - (6)].format_clause),(yyvsp[(3) - (6)].value)); }
    break;

  case 115:

/* Line 1455 of yacc.c  */
#line 329 "control_parse.y"
    {(yyval.fformat)=add_fformat((yyvsp[(1) - (3)].fformat),(yyvsp[(3) - (3)].fformat)); }
    break;

  case 116:

/* Line 1455 of yacc.c  */
#line 330 "control_parse.y"
    {(yyval.fformat)=add_fformat((yyvsp[(1) - (3)].fformat),(yyvsp[(3) - (3)].fformat)); }
    break;

  case 117:

/* Line 1455 of yacc.c  */
#line 333 "control_parse.y"
    {(yyval.fformat)=create_fformat((yyvsp[(2) - (2)].string),2); }
    break;

  case 118:

/* Line 1455 of yacc.c  */
#line 334 "control_parse.y"
    {(yyval.fformat)=create_fformat((yyvsp[(2) - (2)].string),1); }
    break;

  case 119:

/* Line 1455 of yacc.c  */
#line 335 "control_parse.y"
    {(yyval.fformat)=create_fformat((yyvsp[(2) - (2)].string),4); }
    break;

  case 120:

/* Line 1455 of yacc.c  */
#line 336 "control_parse.y"
    {(yyval.fformat)=create_fformat(&(yyvsp[(2) - (2)].value),3); }
    break;

  case 121:

/* Line 1455 of yacc.c  */
#line 339 "control_parse.y"
    {(yyval.f_atom)=make_f_atom((yyvsp[(1) - (1)].value),0);}
    break;

  case 122:

/* Line 1455 of yacc.c  */
#line 340 "control_parse.y"
    {(yyval.f_atom)=make_f_atom((yyvsp[(1) - (2)].value),1);}
    break;

  case 123:

/* Line 1455 of yacc.c  */
#line 341 "control_parse.y"
    {(yyval.f_atom)=make_f_atom(1,1);}
    break;

  case 124:

/* Line 1455 of yacc.c  */
#line 342 "control_parse.y"
    {(yyval.f_atom)=make_f_atom(0,1); scan_error|=FORMAT_ERR; }
    break;

  case 125:

/* Line 1455 of yacc.c  */
#line 343 "control_parse.y"
    {(yyval.f_atom)=make_f_atom(0,1); scan_error|=FORMAT_ERR; }
    break;

  case 126:

/* Line 1455 of yacc.c  */
#line 347 "control_parse.y"
    { if(LogFile) free(LogFile); LogFile=(yyvsp[(2) - (2)].string); }
    break;

  case 127:

/* Line 1455 of yacc.c  */
#line 350 "control_parse.y"
    { if(OutputFile) free(OutputFile); OutputFile=(yyvsp[(2) - (2)].string); }
    break;

  case 128:

/* Line 1455 of yacc.c  */
#line 351 "control_parse.y"
    { if(OutputLaurFile) free(OutputLaurFile); OutputLaurFile=(yyvsp[(2) - (2)].string); }
    break;

  case 129:

/* Line 1455 of yacc.c  */
#line 352 "control_parse.y"
    { if(OutputRawFile) free(OutputRawFile); OutputRawFile=(yyvsp[(2) - (2)].string); }
    break;

  case 130:

/* Line 1455 of yacc.c  */
#line 356 "control_parse.y"
    { if(ErrorDir) free(ErrorDir); ErrorDir=(yyvsp[(2) - (2)].string);}
    break;

  case 131:

/* Line 1455 of yacc.c  */
#line 359 "control_parse.y"
    { do_missing_com((yyvsp[(2) - (2)].express),0,0); free((yyvsp[(2) - (2)].express)); }
    break;

  case 132:

/* Line 1455 of yacc.c  */
#line 360 "control_parse.y"
    { do_missing_com((yyvsp[(2) - (3)].express),(yyvsp[(3) - (3)].var_list),0); free((yyvsp[(2) - (3)].express)); }
    break;

  case 133:

/* Line 1455 of yacc.c  */
#line 361 "control_parse.y"
    { do_missing_com((yyvsp[(2) - (4)].express),(yyvsp[(4) - (4)].var_list),0); free((yyvsp[(2) - (4)].express)); }
    break;

  case 134:

/* Line 1455 of yacc.c  */
#line 362 "control_parse.y"
    { do_missing_com((yyvsp[(5) - (5)].express),0,(yyvsp[(3) - (5)].string)); free((yyvsp[(5) - (5)].express)); }
    break;

  case 135:

/* Line 1455 of yacc.c  */
#line 365 "control_parse.y"
    { do_ped_com((yyvsp[(2) - (2)].var_list)); }
    break;

  case 137:

/* Line 1455 of yacc.c  */
#line 369 "control_parse.y"
    { change_type(ST_TRAITLOCUS,(yyvsp[(3) - (3)].var_list)); free_vlist((yyvsp[(3) - (3)].var_list));}
    break;

  case 139:

/* Line 1455 of yacc.c  */
#line 373 "control_parse.y"
    { change_type(ST_CONSTANT,(yyvsp[(2) - (2)].var_list)); free_vlist((yyvsp[(2) - (2)].var_list));}
    break;

  case 140:

/* Line 1455 of yacc.c  */
#line 374 "control_parse.y"
    { change_type(ST_MULTIPLE,(yyvsp[(2) - (2)].var_list)); free_vlist((yyvsp[(2) - (2)].var_list));}
    break;

  case 141:

/* Line 1455 of yacc.c  */
#line 375 "control_parse.y"
    { change_type(ST_RANDOM|ST_FACTOR,(yyvsp[(2) - (2)].var_list)); free_vlist((yyvsp[(2) - (2)].var_list));}
    break;

  case 142:

/* Line 1455 of yacc.c  */
#line 376 "control_parse.y"
    { change_type(ST_FACTOR,(yyvsp[(2) - (2)].var_list)); free_vlist((yyvsp[(2) - (2)].var_list));}
    break;

  case 143:

/* Line 1455 of yacc.c  */
#line 377 "control_parse.y"
    {change_type(ST_REALTYPE,(yyvsp[(2) - (2)].var_list)); free_vlist((yyvsp[(2) - (2)].var_list));}
    break;

  case 144:

/* Line 1455 of yacc.c  */
#line 378 "control_parse.y"
    {change_type(ST_INTTYPE,(yyvsp[(2) - (2)].var_list)); free_vlist((yyvsp[(2) - (2)].var_list));}
    break;

  case 145:

/* Line 1455 of yacc.c  */
#line 381 "control_parse.y"
    { (yyval.value)=LINK_AUTO; }
    break;

  case 146:

/* Line 1455 of yacc.c  */
#line 382 "control_parse.y"
    { (yyval.value)=LINK_X; }
    break;

  case 147:

/* Line 1455 of yacc.c  */
#line 383 "control_parse.y"
    { (yyval.value)=LINK_Y; }
    break;

  case 148:

/* Line 1455 of yacc.c  */
#line 386 "control_parse.y"
    { do_link_com(0,(yyvsp[(1) - (2)].value),(yyvsp[(2) - (2)].var_list)); }
    break;

  case 149:

/* Line 1455 of yacc.c  */
#line 387 "control_parse.y"
    { do_link_com((yyvsp[(2) - (4)].string),(yyvsp[(1) - (4)].value),(yyvsp[(4) - (4)].var_list)); }
    break;

  case 150:

/* Line 1455 of yacc.c  */
#line 388 "control_parse.y"
    { do_link_com((yyvsp[(2) - (2)].string),(yyvsp[(1) - (2)].value),0); }
    break;

  case 151:

/* Line 1455 of yacc.c  */
#line 391 "control_parse.y"
    {
	if(Filter) {
	  print_scan_warn("Line %d: Warning - Filter defined twice\n",lineno);
	  free(Filter);
   }
	Filter=(yyvsp[(2) - (2)].string); }
    break;

  case 152:

/* Line 1455 of yacc.c  */
#line 399 "control_parse.y"
    {do_model_com((yyvsp[(4) - (4)].model_list),(yyvsp[(2) - (4)].var),0);}
    break;

  case 153:

/* Line 1455 of yacc.c  */
#line 400 "control_parse.y"
    {do_model_com((yyvsp[(7) - (7)].model_list),(yyvsp[(2) - (7)].var),(yyvsp[(4) - (7)].express)); if((yyvsp[(4) - (7)].express)) free((yyvsp[(4) - (7)].express)); }
    break;

  case 154:

/* Line 1455 of yacc.c  */
#line 403 "control_parse.y"
    {set_group((yyvsp[(2) - (2)].element));}
    break;

  case 155:

/* Line 1455 of yacc.c  */
#line 406 "control_parse.y"
    { (yyval.model_list)=add_to_model((yyvsp[(1) - (3)].model_list),(yyvsp[(3) - (3)].var_list)); }
    break;

  case 156:

/* Line 1455 of yacc.c  */
#line 407 "control_parse.y"
    { (yyval.model_list)=add_to_model((yyvsp[(1) - (3)].model_list),(yyvsp[(3) - (3)].var_list)); }
    break;

  case 157:

/* Line 1455 of yacc.c  */
#line 408 "control_parse.y"
    { (yyval.model_list)=add_to_model(0,(yyvsp[(1) - (1)].var_list)); }
    break;

  case 158:

/* Line 1455 of yacc.c  */
#line 409 "control_parse.y"
    { (yyval.model_list)=add_to_model(0,(yyvsp[(1) - (1)].var_list)); }
    break;

  case 159:

/* Line 1455 of yacc.c  */
#line 412 "control_parse.y"
    { (yyval.var_list)=add_var_lists((yyvsp[(1) - (3)].var_list),(yyvsp[(3) - (3)].var_list)); }
    break;

  case 160:

/* Line 1455 of yacc.c  */
#line 413 "control_parse.y"
    { (yyval.var_list)=add_var_lists((yyvsp[(1) - (3)].var_list),(yyvsp[(3) - (3)].var_list)); }
    break;

  case 161:

/* Line 1455 of yacc.c  */
#line 414 "control_parse.y"
    { (yyval.var_list)=add_var_lists((yyvsp[(1) - (3)].var_list),(yyvsp[(3) - (3)].var_list)); }
    break;

  case 162:

/* Line 1455 of yacc.c  */
#line 415 "control_parse.y"
    { (yyval.var_list)=add_var_lists((yyvsp[(1) - (3)].var_list),(yyvsp[(3) - (3)].var_list)); }
    break;

  case 164:

/* Line 1455 of yacc.c  */
#line 419 "control_parse.y"
    { (yyval.var)=create_var("SEX"); }
    break;

  case 165:

/* Line 1455 of yacc.c  */
#line 422 "control_parse.y"
    { (yyval.var_list)=add_to_var_list(0,(yyvsp[(1) - (1)].var),0); }
    break;

  case 166:

/* Line 1455 of yacc.c  */
#line 423 "control_parse.y"
    { (yyval.var_list)=add_to_var_list(0,(yyvsp[(1) - (1)].var),0); }
    break;

  case 167:

/* Line 1455 of yacc.c  */
#line 424 "control_parse.y"
    { (yyval.var_list)=add_to_var_list(0,(yyvsp[(1) - (4)].var),(yyvsp[(3) - (4)].express)); if((yyvsp[(3) - (4)].express)) free((yyvsp[(3) - (4)].express)); }
    break;

  case 168:

/* Line 1455 of yacc.c  */
#line 427 "control_parse.y"
    { (yyval.element)=get_element((yyvsp[(1) - (1)].var),0); }
    break;

  case 169:

/* Line 1455 of yacc.c  */
#line 428 "control_parse.y"
    { (yyval.element)=get_element((yyvsp[(1) - (4)].var),(yyvsp[(3) - (4)].express)); if((yyvsp[(3) - (4)].express)) free((yyvsp[(3) - (4)].express)); }
    break;

  case 172:

/* Line 1455 of yacc.c  */
#line 433 "control_parse.y"
    { if((yyvsp[(1) - (2)].element)) set_locus_element((yyvsp[(1) - (2)].element)); }
    break;

  case 173:

/* Line 1455 of yacc.c  */
#line 434 "control_parse.y"
    { if((yyvsp[(1) - (1)].element)) set_locus_element((yyvsp[(1) - (1)].element)); }
    break;

  case 174:

/* Line 1455 of yacc.c  */
#line 435 "control_parse.y"
    { if((yyvsp[(1) - (1)].var)) set_locus_array((yyvsp[(1) - (1)].var)); }
    break;

  case 175:

/* Line 1455 of yacc.c  */
#line 438 "control_parse.y"
    { if((yyvsp[(2) - (3)].element)) set_haplo_element((yyvsp[(2) - (3)].element),0); }
    break;

  case 176:

/* Line 1455 of yacc.c  */
#line 439 "control_parse.y"
    { if((yyvsp[(2) - (5)].element)) set_haplo_element((yyvsp[(2) - (5)].element),(yyvsp[(4) - (5)].element)); }
    break;

  case 180:

/* Line 1455 of yacc.c  */
#line 445 "control_parse.y"
    { if((yyvsp[(1) - (4)].element)) set_slocus_element((yyvsp[(1) - (4)].element),(yyvsp[(3) - (4)].var_list)); }
    break;

  case 181:

/* Line 1455 of yacc.c  */
#line 446 "control_parse.y"
    { yyerror ("No marker list for Super Locus\n"); }
    break;

  case 184:

/* Line 1455 of yacc.c  */
#line 451 "control_parse.y"
    { (yyval.var_list)=add_var_lists((yyvsp[(1) - (3)].var_list),(yyvsp[(3) - (3)].var_list)); }
    break;

  case 185:

/* Line 1455 of yacc.c  */
#line 454 "control_parse.y"
    { start_loopclause(); }
    break;

  case 186:

/* Line 1455 of yacc.c  */
#line 457 "control_parse.y"
    { free_vlist((yyvsp[(2) - (8)].var_list)); begin_looping((yyvsp[(5) - (8)].element),(yyvsp[(7) - (8)].express),0); }
    break;

  case 187:

/* Line 1455 of yacc.c  */
#line 458 "control_parse.y"
    { free_vlist((yyvsp[(2) - (10)].var_list)); begin_looping((yyvsp[(5) - (10)].element),(yyvsp[(7) - (10)].express),(yyvsp[(9) - (10)].express)); }
    break;

  case 188:

/* Line 1455 of yacc.c  */
#line 461 "control_parse.y"
    { (yyval.var_list)=(yyvsp[(3) - (4)].var_list); in_loopclause=0; }
    break;

  case 189:

/* Line 1455 of yacc.c  */
#line 462 "control_parse.y"
    { (yyval.var_list)=0; in_loopclause=0; }
    break;

  case 190:

/* Line 1455 of yacc.c  */
#line 465 "control_parse.y"
    { (yyval.var_list)=add_to_var_list(0,0,0); }
    break;

  case 191:

/* Line 1455 of yacc.c  */
#line 466 "control_parse.y"
    { (yyval.var_list)=add_to_var_list(0,(yyvsp[(1) - (1)].var),0); }
    break;

  case 192:

/* Line 1455 of yacc.c  */
#line 467 "control_parse.y"
    { (yyval.var_list)=add_to_var_list(0,(yyvsp[(1) - (1)].var),0); }
    break;

  case 193:

/* Line 1455 of yacc.c  */
#line 468 "control_parse.y"
    { (yyval.var_list)=add_to_var_list(0,(yyvsp[(1) - (4)].var),(yyvsp[(3) - (4)].express)); if((yyvsp[(3) - (4)].express)) free((yyvsp[(3) - (4)].express)); }
    break;

  case 196:

/* Line 1455 of yacc.c  */
#line 473 "control_parse.y"
    { (yyval.var_list)=add_var_lists((yyvsp[(1) - (3)].var_list),(yyvsp[(3) - (3)].var_list)); }
    break;

  case 197:

/* Line 1455 of yacc.c  */
#line 474 "control_parse.y"
    { (yyval.var_list)=add_var_lists((yyvsp[(1) - (3)].var_list),(yyvsp[(3) - (3)].var_list)); }
    break;

  case 200:

/* Line 1455 of yacc.c  */
#line 479 "control_parse.y"
    { (yyval.var_list)=add_var_lists((yyvsp[(1) - (3)].var_list),(yyvsp[(3) - (3)].var_list)); }
    break;

  case 201:

/* Line 1455 of yacc.c  */
#line 480 "control_parse.y"
    { (yyval.var_list)=add_var_lists((yyvsp[(1) - (3)].var_list),(yyvsp[(3) - (3)].var_list)); }
    break;

  case 202:

/* Line 1455 of yacc.c  */
#line 483 "control_parse.y"
    { (yyval.string) = (yyvsp[(1) - (1)].string); }
    break;

  case 203:

/* Line 1455 of yacc.c  */
#line 484 "control_parse.y"
    { (yyval.string) = string_copy((yyvsp[(1) - (3)].string),(yyvsp[(3) - (3)].string)); free((yyvsp[(3) - (3)].string)); }
    break;

  case 204:

/* Line 1455 of yacc.c  */
#line 485 "control_parse.y"
    { if((yyvsp[(1) - (3)].element) && ((yyvsp[(1) - (3)].element)->type&ST_STRING)) (yyval.string) = string_copy((yyvsp[(3) - (3)].string),(yyvsp[(1) - (3)].element)->arg.string); else (yyval.string)=(yyvsp[(3) - (3)].string); }
    break;

  case 205:

/* Line 1455 of yacc.c  */
#line 488 "control_parse.y"
    { (yyval.string) = (yyvsp[(1) - (1)].string); }
    break;

  case 206:

/* Line 1455 of yacc.c  */
#line 489 "control_parse.y"
    { if((yyvsp[(1) - (1)].element) && ((yyvsp[(1) - (1)].element)->type&ST_STRING)) (yyval.string) = string_copy(0,(yyvsp[(1) - (1)].element)->arg.string); else (yyval.string)=0; }
    break;

  case 207:

/* Line 1455 of yacc.c  */
#line 490 "control_parse.y"
    { (yyval.string) = string_copy((yyvsp[(1) - (3)].string),(yyvsp[(3) - (3)].string)); free((yyvsp[(3) - (3)].string)); }
    break;

  case 208:

/* Line 1455 of yacc.c  */
#line 491 "control_parse.y"
    { if((yyvsp[(3) - (3)].element) && ((yyvsp[(3) - (3)].element)->type&ST_STRING)) (yyval.string) = string_copy((yyvsp[(1) - (3)].string),(yyvsp[(3) - (3)].element)->arg.string); else (yyval.string)=(yyvsp[(1) - (3)].string); }
    break;



/* Line 1455 of yacc.c  */
#line 3199 "y.tab.c"
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
#line 493 "control_parse.y"


static void enter_loop(void)
{	
	if(loop_level<MAX_LOOP)	{
		loop_stat[loop_level]=loop_record;
		loop_ptr[loop_level]=loop_main_ptr;
		loop_level++;
		if(!loop_record) loop_record=1;
	} else yyerror("Too many nested loops\n");
}

static void start_loopclause(void)
{	
	if(loop_level<MAX_LOOP)	{
		loop_stat[loop_level]=loop_record;
		loop_ptr[loop_level]=loop_main_ptr;
		in_loopclause=1;
		loop_level++;
		if(!loop_record) loop_record=1;
	} else yyerror("Too many nested loops\n");
}

static void begin_looping(struct var_element *element, struct express *exp1, struct express *exp_2)
{	
	int er=0,i;
	
	if(element && exp1) {
		if(exp1->type==ST_INTEGER) loop_clause_end=(int)exp1->arg.value;
		else er=1;
		if(!er && exp_2) {
			if(exp_2->type==ST_INTEGER) loop_clause_step=(int)exp_2->arg.value;
			else er=1;
		} else loop_clause_step=1;
		if(element->type&ST_INTEGER) {
			i=(int)element->arg.value;
			if(loop_clause_step<0) {
				if(i<loop_clause_end) er= -1;
			} else if(i>loop_clause_end) er= -1;
		} else er=1;
	} else er=2;
	if(er) {
		switch(er) {
		 case 1:
			yyerror("Loop variable not integer type\n");
			break;
		 case 2:
			yyerror("Syntax error\n");
			break;
		}
		loop_record=loop_stat[--loop_level];
		if(!loop_record) loop_main_ptr=loop_level?loop_ptr[loop_level-1]:0;
		in_loopclause=0;
	} else {
		in_loopclause= -1;
		loop_record= -1;
		loop_clause_ptr=loop_main_ptr;
		loop_clause_element=element;
		loop_main_ptr=loop_ptr[loop_level-1];
	}
	if(exp1) free(exp1);
	if(exp_2) free(exp_2);
}

static int check_pos(struct express *e,double *x)
{
	int er=0;
	
	if(e->type==ST_INTEGER) *x=(double)e->arg.value;
	else if(e->type==ST_REAL) *x=e->arg.rvalue;
	else {
		yyerror1("Non-numeric arguments to position command\n");
		er=1;
	}
	return er;
}

static void set_position(struct var_element *elem,struct express *sex_avg,struct express *male,struct express *female)
{
	int i,set[3],er=0;
	double x[3];
	struct marker_info *mi;
	
	/* First check if positions are valid */
	if(sex_avg) {
		if(check_pos(sex_avg,x+2)) er=1;
		else set[2]=1;
	} else set[2]=0;
	if(!er && female) {
		if(check_pos(male,x)) er=1;
		else set[0]=1;
	} else set[0]=0;
	if(!er && male) {
		if(check_pos(female,x+1)) er=1;
		else set[1]=1;
	} else set[1]=0;
	/* If OK, then check if this element has already got an info entry */
	if(!er && (set[0] || set[1] || set[2])) {
		mi=m_info;
		while(mi) {
			if(mi->element==elem) break;
			mi=mi->next;
		}
		if(mi) {
			if((mi->pos_set[0] && set[0]) || (mi->pos_set[1] && set[1]) || (mi->pos_set[2] && set[2])) {
				yyerror1("Error: duplicate positions set for marker\n");
				er=1;
			} else {
				for(i=0;i<3;i++) mi->pos_set[i]|=set[i];
			}
		} else {
			mi=lk_malloc(sizeof(struct marker_info));
			mi->next=m_info;
			m_info=mi;
			for(i=0;i<3;i++) mi->pos_set[i]=set[i];
			mi->element=elem;
		}
		if(!er) for(i=0;i<3;i++) if(set[i]) mi->pos[i]=x[i];
	}
}

static int if_true(struct express *express)
{
	int l=0;
	
	if(express)	{
		switch(express->type) {
		 case ST_INTEGER:
			l=(express->arg.value!=0);
			break;
		 case ST_REAL:
			l=(express->arg.rvalue!=0.0);
			break;
		 case ST_STRING:
			if(express->arg.string && express->arg.string[0]) l=1;
		}
	}
	return l;
}

static void do_while_com(struct express *express)
{
	int i;
	
	if(loop_level)	{
		if(!scan_error_n && if_true(express)) {
			loop_record= -1;
			loop_main_ptr=loop_ptr[loop_level-1];
		} else {
			loop_record=loop_stat[--loop_level];
			if(!loop_record) {
				i=loop_main_ptr-1;
				loop_main_ptr=loop_level?loop_ptr[loop_level-1]:0;
				for(;i>=loop_main_ptr;i--)	{
					if(loop_stack[i].token==STRING) free(loop_stack[i].yylval.string);
				}
			}
		}
	} else yyerror("WHILE outside of do loop\n");
}

static void print_exp(struct express *express)
{
	if(!express) return;
	switch(express->type) {
	 case ST_STRING:
		(void)fputs(express->arg.string,stdout);
		free(express->arg.string);
		break;
	 case ST_INTEGER:
		(void)printf("%ld",express->arg.value);
		break;
	 case ST_REAL:
		(void)printf("%g",express->arg.rvalue);
		break;
	}
}

static void set_sex(struct var_element *elem,struct express *exp_1,struct express *exp_2)
{
	struct sex_def *se;
	
	if(!exp_1 || !exp_2) yyerror1("Null arguments to sex command\n");
	else {
		if(exp_1->type != exp_2->type) yyerror1("Arguments to sex command of different type\n");
		else if(exp_1->type!=ST_INTEGER && exp_1->type!=ST_STRING) yyerror1("Arguments to sex command of invalid type\n");
		else {
			se=lk_malloc(sizeof(struct sex_def));
			se->sex_exp[0]=exp_1;
			se->sex_exp[1]=exp_2;
			se->sex_elem=elem;
			elem->type|=(ST_SEX|ST_FACTOR|ST_CONSTANT|ST_DATA);
			if(exp_1->type==ST_INTEGER) elem->type|=ST_INTTYPE;
			se->next=sex_def;
			sex_def=se;
		}
	}
}

static void set_group(struct var_element *elem)
{
	if(group_elem) yyerror1("Error: Multiple group commands");
	else {
		group_elem=elem;
		group_elem->type|=(ST_GROUP|ST_FACTOR|ST_CONSTANT);
	}
}

static void do_ped_com(struct var_list *vlist)
{
	int i,j,n;
	struct var_list *vlist1;
	struct scan_data *sd;
	
	if(pedflag) {
		yyerror1("Error: Multiple pedigree commands");
		scan_error|=PED_ERR;
	}
	pedflag=1;
	n=0;
	while(vlist) {
		sd=vlist->var->data;
		if(sd->vtype&ST_ARRAY) {
			if(vlist->index) {
				if(n<4) pedlist[n]=sd->element+vlist->index-1;
				n++;
			} else for(i=0;i<sd->n_elements;i++)	{
				if(n<4) pedlist[n]=sd->element+i;
				n++;
			}
		} else {
			if(n<4) pedlist[n]=sd->element;
			n++;
		}
		vlist1=vlist->next;
		free(vlist);
		vlist=vlist1;
	}
	if(n!=3 && n!=4) {
		yyerror1("Error: Wrong no. variables for pedigree command (3 or 4 required)");
		scan_error|=PED_ERR;
	} else {
		for(i=1;i<n;i++) {
			for(j=0;j<i;j++) if(pedlist[i]==pedlist[j]) {
				yyerror1("Error: Repeated variables in pedigree command");
				scan_error|=PED_ERR;
			}
		}
		if(n==4) {
			pedlist[0]->type|=(ST_FAMILY|ST_FACTOR|ST_CONSTANT);
			family_id=1;
			i=1;
		} else i=0;
		pedlist[i++]->type|=(ST_ID|ST_FACTOR|ST_CONSTANT);
		pedlist[i++]->type|=(ST_SIRE|ST_FACTOR|ST_CONSTANT);
		pedlist[i++]->type|=(ST_DAM|ST_FACTOR|ST_CONSTANT);
	}
}

static struct express *alloc_express(void)
{
	struct express *e;
	
	e=lk_malloc(sizeof(struct express));
	return e;
}

static struct express *do_logical_op(struct express *ex1,struct express *ex2,int op)
{
	int i=0,l=0;
	double rv1,rv2;
	char *s1,*s2;
	
	s1=s2=0;
	if(ex1->type&ST_STRING) {
		s1=ex1->arg.string;
		i++;
	}
	if(ex2 && ex2->type&ST_STRING) {
		s2=ex2->arg.string;
		i++;
	}
	if(i==2)	{
		switch(op) {
		 case '=':
			l=mystrcmp(s1,s2)?0:1;
			break;
		 case NEQSYMBOL:
			l=mystrcmp(s1,s2)?1:0;
			break;
		 case '<':
			l=mystrcmp(s1,s2)<0?1:0;
			break;
		 case '>':
			l=mystrcmp(s1,s2)>0?1:0;
			break;
		 case LEQSYMBOL:
			l=mystrcmp(s1,s2)>=0?1:0;
			break;
		 case GEQSYMBOL:
			l=mystrcmp(s1,s2)>=0?1:0;
			break;
		 case ORSYMBOL:
			l=(s1 || s2);
			break;
		 case ANDSYMBOL:
			l=(s1 && s2);
			break;
		 default:
			ABT_FUNC("Internal error - invalid string op\n");
		}
	} else if(i && !ex2)	{
		assert(op==NOTSYMBOL);
		l=s1?0:1;
	} else if(i) {
		switch(op) {
		 case ORSYMBOL:
			if(ex1->type&ST_STRING) l=(s1 || ex2->arg.value);
			else l=(s2 || ex1->arg.value);
			break;
		 case ANDSYMBOL:
			if(ex1->type&ST_STRING) l=(s1 && ex2->arg.value);
			else l=(s2 && ex1->arg.value);
			break;
		 default:
			ABT_FUNC("Internal error - invalid string op\n");
		}
	} else {
		rv1=rv2=0.0;
		if(ex1->type&ST_INTEGER) rv1=(double)ex1->arg.value;
		else if(ex1->type&ST_REAL) rv1=ex1->arg.rvalue;
		if(ex2) {
			if(ex2->type&ST_INTEGER) rv2=(double)ex2->arg.value;
			else if(ex2->type&ST_REAL) rv2=ex2->arg.rvalue;
		}
		switch(op) {
		 case '=':
			l=(rv1==rv2);
			break;
		 case NEQSYMBOL:
			l=(rv1!=rv2);
			break;
		 case '<':
			l=(rv1<rv2);
			break;
		 case '>':
			l=(rv1>rv2);
			break;
		 case LEQSYMBOL:
			l=(rv1<=rv2);
			break;
		 case GEQSYMBOL:
			l=(rv1>=rv2);
			break;
		 case ORSYMBOL:
			l=(rv1 || rv2);
			break;
		 case ANDSYMBOL:
			l=(rv1 && rv2);
			break;
		 case NOTSYMBOL:
			l=(rv1==0.0);
			break;
		 default:
			ABT_FUNC("Internal error - invalid op\n");
		}
	}
	if(ex2) free(ex2);
	ex1->type=ST_INTEGER;
	ex1->arg.value=l;
	return ex1;
}

static struct express *do_express_op(struct express *ex1,struct express *ex2,int op)
{
	double rv1,rv2;
	int i;
	
	if(ex1->type&ST_STRING)	{
		if(ex2 && ex2->type&ST_STRING) {
			if(op!='+') yyerror("Illegal string operation\n");
			else {
				ex1->arg.string=string_copy(ex1->arg.string,ex2->arg.string);
				free(ex2->arg.string);
			}
		} else yyerror("Can't mix numeric and string expressions\n");
	} else if(ex2 && ex2->type&ST_STRING) yyerror("Can't mix numeric and string expressions\n");
	else {
		rv1=rv2=0.0;
		if(ex1->type&ST_INTEGER) rv1=(double)ex1->arg.value;
		else if(ex1->type&ST_REAL) rv1=ex1->arg.rvalue;
		if(ex2) {
			if(ex2->type&ST_INTEGER) rv2=(double)ex2->arg.value;
			else if(ex2->type&ST_REAL) rv2=ex2->arg.rvalue;
		}
		switch(op) {
		 case '+': 
			rv1+=rv2;
			break;
		 case '-':
			if(ex2) rv1-=rv2;
			else rv1= -rv1;
			break;
		 case '*':
			rv1*=rv2;
			break;
		 case '/':
			if(rv2==0.0) {
				yyerror("Divide by zero error\n");
				rv1=0.0;
			} else rv1/=rv2;
			break;
		}
		i=(int)rv1;
		if((double)i==rv1) {
			ex1->type=ST_INTEGER;
			ex1->arg.value=i;
		} else {
			ex1->type=ST_REAL;
			ex1->arg.rvalue=rv1;
		}
	}
	if(ex2) free(ex2);
	return ex1;
}

static void check_element_add_op(struct var_element *element)
{
	switch(element->type&(ST_REAL|ST_INTEGER|ST_STRING)) {
	 case ST_STRING:
		add_operation(string_copy(0,element->arg.string),STRING,0);
		break;
	 case ST_INTEGER:
		add_operation(&element->arg,INTEGER,0);
		break;
	 case ST_REAL:
		add_operation(&element->arg,REAL,0);
		break;
	 case 0:
		add_operation(element,VARIABLE,0);
		break;
	 default:
		ABT_FUNC("Internal error - illegal element type\n");
	}
}

static int check_index(struct scan_data *sd,struct express *express)
{
	int i;

	if(express->type!=ST_INTEGER) {
		if(in_loopclause<=0) yyerror("Non-integral expression for array index");
	} else if(sd->vtype&ST_ARRAY) {
		i=(int)express->arg.value;
		if(i<1 || i>sd->n_elements) {
			if(in_loopclause<=0) yyerror("Array index out of bounds");
		} else return i;
	} else yyerror("Not an array");
	return 0;
}

static struct var_element *get_element(struct bin_node *node,struct express *express)
{
	int i;
	struct scan_data *sd;
	
	sd=node->data;
	if(express) {
		if(!(i=check_index(sd,express))) return 0;
		return sd->element+i-1;
	} else {
		if(sd->vtype&ST_ARRAY)	{
			yyerror("Illegal reference to array");
			return 0;
		}
	}
	return sd->element;
}

static void set_array_var(struct scan_data *sd,struct express *express)
{
	int i;

	if(express->type!=ST_INTEGER) yyerror("Non-integral expression for array size");
	else if((i=(int)express->arg.value)<1) yyerror("Illegal array size");
	else if(sd->vtype) yyerror("Can't redefine variable");
	else {
		sd->vtype|=ST_ARRAY;
		sd->n_elements=i;
		free(sd->element);
		sd->element=lk_calloc((size_t)sd->n_elements,sizeof(struct var_element));
	}
}

static int count_var_list(struct var_list *vlist)
{
	int i=0;
	struct scan_data *sd;
	
	while(vlist) {
		sd=vlist->var?vlist->var->data:0;
		if(sd && (sd->vtype&ST_ARRAY) && !vlist->index) i+=sd->n_elements;
		else i++;
		vlist=vlist->next;
	}
	return i;
}

static struct var_element *assign_var(struct bin_node *node,struct express *ix,struct express *express)
{
	struct var_element *element;
	struct scan_data *sd;
	
	if(!express) return 0;
	if(!(element=get_element(node,ix))) return 0;
	switch(express->type) {
	 case ST_STRING:
	     element->arg.string=express->arg.string;
		RemBlock=AddRemem(element->arg.string,RemBlock);
		break;
	 case ST_REAL:
	 case ST_INTEGER:
	     element->arg=express->arg;
		break;
	 case 0:
	     yyerror1("Undefined assignment\n");
		element->type=0;
		element->arg.string=0;
		break;
	 default:
		ABT_FUNC(IntErr);
	}
	if(!ix) {
		sd=node->data;
		sd->vtype|=ST_SCALAR;
	}
	element->type=express->type;
	return element;
}

void check_vars(struct bin_node *node,int *i,void check_func(struct bin_node *,int *))
{
	if(node->left) {
		check_vars(node->left,i,check_func);
	}
	check_func(node,i);
	if(node->right) {
		check_vars(node->right,i,check_func);
	}
}

static void check_vars_1(struct bin_node *node,void check_func(struct bin_node *))
{
	if(node->left) {
		check_vars_1(node->left,check_func);
	}
	check_func(node);
	if(node->right) {
		check_vars_1(node->right,check_func);
	}
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
	
	if(scan_warn_n<max_scan_warnings) {
		va_start(args,fmt);
		(void)vfprintf(stderr,fmt,args);
		va_end(args);
	}
	scan_warn_n++;
}

static void add_operation(void *arg,int type, int op)
{
	struct operation *o;
	
	assert(arg || !type);
	o=lk_malloc(sizeof(struct operation));
	o->next=Op_List;
	Op_List=o;
	o->type=type;
	o->op=op;
	switch(type) {
	 case VARIABLE:
		o->arg.element= (struct var_element *)arg;
		break;
	 case INTEGER:
		o->arg.value= *(int *)arg;
		break;
	 case REAL:
		o->arg.rvalue= *(double *)arg;
		break;
	 case STRING:
		o->arg.string= (char *)arg;
		break;
	}
}

static void new_command(void)
{
	shell_flag=in_loopclause=0;
	Op_List=0;
	iflag=0;
}

static struct model_list *add_to_model(struct model_list *model,struct var_list *vlist)
{
	struct var_list *vlist1;
	struct scan_data *sd;
	struct model_list *m1;
	int i;
	
	m1=lk_malloc(sizeof(struct model_list));
	if(vlist) {
		i=count_var_list(vlist);
		m1->element=lk_malloc(sizeof(void *)*i);
		i=0;
		while(vlist) {
			sd=vlist->var->data;
			if(sd->vtype&ST_ARRAY)	{
				if(vlist->index) {
					sd->element[vlist->index-1].type|=ST_MODEL;
					sd->element[vlist->index-1].index=vlist->index;
					sd->element[vlist->index-1].oindex=vlist->index;
					m1->element[i++]=sd->element+vlist->index-1;
				} else yyerror("Error - Can't use whole arrays as model parameters");
			} else {
				sd->element[0].type|=ST_MODEL;
				sd->element[0].index=0;
				m1->element[i++]=sd->element;
			}
			vlist1=vlist->next;
			free(vlist);
			vlist=vlist1;
		}
		m1->nvar=i;
	} else ABT_FUNC("Nothing to add...\n");
	m1->next=model;
	return m1;
}

static char *string_copy(char *s1,char *s2)
{
   size_t sz;	
	if(s1) { 
	   sz=strlen(s2);
		s1=lk_realloc(s1,strlen(s1)+sz+1);
		(void)strncat(s1,s2,sz);
	} else {
	  sz=strlen(s2)+1;
		s1=lk_malloc(sz);
		(void)strncpy(s1,s2,sz);
	}
	return s1;
}

static struct format *setup_format(struct format_clause *fc)
{
	int i,n=0,pp=0;
	struct format_atom **fa;
	struct format *format;
	
	fa=fc->f_atoms;
	for(i=0;i<fc->n_atoms;i++) if(!fa[i]->pos) n++;
	if(!n) {
		if(!(scan_error&FORMAT_ERR)) yyerror("Error - Empty format clause");
		free(fa);
		free(fc);
		scan_error|=FORMAT_ERR;
		scan_error_n++;
		return 0;
	}
	format=lk_malloc(sizeof(struct format));
	format->line=lineno;
	format->f_atoms=lk_malloc(sizeof(struct format_atom)*n);
	for(i=n=0;i<fc->n_atoms;i++) {
		if(!fa[i]->pos) {
			format->f_atoms[n].size=fa[i]->size;
			format->f_atoms[n++].pos=pp;
		}
		pp+=fa[i]->size;
	}
	free(fa);
	format->n_atoms=n;
	f_atom_n=0;
	free(fc);
	return format;
}

static struct format_atom *make_f_atom(int n,int flag)
{
	if(f_atom_n>=f_atom_size) {
		f_atom_size*=2;
		f_atom_list=lk_realloc(f_atom_list,sizeof(struct format_atom)*f_atom_size);
	}
	f_atom_list[f_atom_n].size=n;
	f_atom_list[f_atom_n].pos=flag;
	return &f_atom_list[f_atom_n++];
}

static struct fformat *add_fformat(struct fformat *f1,struct fformat *f2)
{
	if(f2->rs) {
		if(f1->rs) free(f1->rs);
		f1->rs=f2->rs;
	}
	if(f2->fs) {
		if(f1->fs) free(f1->fs);
		f1->fs=f2->fs;
	}
	if(f2->gs) {
		if(f1->gs) free(f1->gs);
		f1->gs=f2->gs;
	}	
	if(f2->skip) f1->skip=f2->skip;
	free(f2);
	return f1;
}

static struct fformat *create_fformat(void *p,int fg)
{
	struct fformat *ff;
	int *i;
	
	ff=lk_malloc(sizeof(struct fformat));
	ff->rs=ff->fs=ff->gs=0;
	ff->skip=0;
	switch(fg) {
	 case 1:
		ff->rs=p;
		break;
	 case 2:
		ff->fs=p;
		break;
	 case 3:
		i=p;
		ff->skip=*i;
		break;
	 case 4:
		ff->gs=p;
		break;
	 default:
		ABT_FUNC("Internal error - incorrect flag\n");
	}
	return ff;
}

static struct format_clause *add_f_atom(struct format_clause *fc,struct format_atom *fa)
{
	if(!fc) {
		fc=lk_malloc(sizeof(struct format_clause));
		fc->fc_size=16;
		fc->n_atoms=0;
		fc->f_atoms=lk_malloc(sizeof(struct format_atom *)*fc->fc_size);
	}
	if(fc->n_atoms>=fc->fc_size) {
		fc->fc_size*=2;
		fc->f_atoms=lk_realloc(fc->f_atoms,sizeof(struct format_atom *)*fc->fc_size);
	}
	fc->f_atoms[fc->n_atoms++]=fa;
	return fc;
}

static struct format_clause *add_f_list(struct format_clause *fc,struct format_clause *fc1,int n)
{
	int sz,i,j;
	
	sz=fc1->n_atoms*n;
	if(!fc) {
		fc=lk_malloc(sizeof(struct format_clause));
		fc->fc_size=16;
		if(sz>16) fc->fc_size=sz;
		fc->n_atoms=0;
		fc->f_atoms=lk_malloc(sizeof(struct format_atom *)*fc->fc_size);
	} else {
		if(sz>(fc->fc_size-fc->n_atoms)) {
			fc->fc_size=sz+fc->n_atoms;
			fc->f_atoms=lk_realloc(fc->f_atoms,sizeof(struct format_atom *)*fc->fc_size);
		}
	}
	for(i=0;i<n;i++) for(j=0;j<fc1->n_atoms;j++)
	  fc->f_atoms[fc->n_atoms++]=fc1->f_atoms[j];
	free(fc1->f_atoms);
	free(fc1);
	return fc;
}

static struct bin_node *alloc_var(char *p)
{
	struct bin_node *node;
	struct scan_data *sd;
	int i;
	
	node=lk_malloc(sizeof(struct bin_node));
	node->left=node->right=0;
	node->balance=0;
	sd=lk_malloc(sizeof(struct scan_data));
	node->data=sd;
	sd->vtype=0;
	i=(int)strlen(p);
	sd->name=lk_malloc((size_t)i+1);
	sd->name[i--]=0;
	for(;i>=0;i--) sd->name[i]=toupper((int)p[i]);
	sd->n_elements=1;
	sd->element=lk_calloc(1,sizeof(struct var_element));
	sd->element->arg.element=0;
	return node;
}

static struct bin_node *find_var(char *p,struct bin_node *node,struct bin_node **node1,int *balanced)
{
	int i;
	struct scan_data *sd;
	
	sd=node->data;
	if((i=strcasecmp(p,sd->name))) {
		if(i<0) {
			if(node->left) {
				node->left=find_var(p,node->left,node1,balanced);
			} else {
				*node1=node->left=alloc_var(p);
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
				node->right=find_var(p,node->right,node1,balanced);
			} else {
				*node1=node->right=alloc_var(p);
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
	}
	return node;
}

static void Check_var(struct bin_node *node)
{
	int i;
	struct var_element *element;
	struct scan_data *sd;
	char *nbuf;
	
	if(node->left) Check_var(node->left);
	sd=node->data;
	i=strlen(sd->name)+4+log((double)(sd->n_elements+1))/log(10.0);
	nbuf=lk_malloc((size_t)i);
	for(i=0;i<sd->n_elements;i++) {
		if(sd->vtype&ST_ARRAY) (void)snprintf(nbuf,i,"%s(%d)",sd->name,i+1);
		else (void)strncpy(nbuf,sd->name,(size_t)i);
		element=sd->element+i;
		if(!(element->type&(ST_DATA|ST_TRAITLOCUS|ST_LINKED|ST_LUMPED))) {
			if(element->type&(ST_MODEL|ST_SEX|ST_ID|ST_FAMILY|ST_SIRE|ST_DAM|ST_TRAIT|ST_GROUP))
			  print_scan_err("Error: No data for variable %s\n",nbuf);
		}
		if((element->type&ST_DATA) && (element->type&ST_TRAITLOCUS))
		  print_scan_err("Error: Variable %s can not have data\n",nbuf);
		else if((element->type&ST_LINKED) && !(element->type&(ST_TRAITLOCUS|ST_MARKER)))
		  print_scan_err("Error: Variable %s is not a locus and so can not be linked\n",nbuf);
		else if((element->type&ST_TRAIT) && (element->type&(ST_GROUP|ST_MARKER|ST_TRAITLOCUS|ST_ID|ST_FAMILY|ST_SIRE|ST_DAM|ST_HAPLO|ST_LUMPED|ST_LINKED|ST_STRING|ST_REAL|ST_INTEGER)))
		  print_scan_err("Error: Variable %s inappropriate type for trait\n",nbuf);
		else if((element->type&ST_TRAITLOCUS) && (element->type&(ST_SEX|ST_GROUP|ST_CENSORED|ST_RANDOM|ST_MARKER|ST_ID|ST_FAMILY|ST_SIRE|ST_DAM|ST_HAPLO|ST_STRING|ST_REAL|ST_INTEGER|ST_REALTYPE|ST_INTTYPE)))
		  print_scan_err("Error: Variable %s inappropriate type for trait locus\n",nbuf);
		else if((element->type&ST_MARKER) && (element->type&(ST_SEX|ST_GROUP|ST_CENSORED|ST_REAL|ST_INTEGER|ST_RANDOM|ST_ID|ST_FAMILY|ST_SIRE|ST_DAM|ST_HAPLO|ST_STRING|ST_INTEGER)))
		  print_scan_err("Error: Variable %s inappropriate type for marker\n",nbuf);
		else if((element->type&ST_HAPLO) && (element->type&(ST_SEX|ST_GROUP|ST_CENSORED|ST_REAL|ST_RANDOM|ST_ID|ST_FAMILY|ST_SIRE|ST_DAM|ST_STRING|ST_REAL|ST_INTEGER)))
		  print_scan_err("Error: Variable %s inappropriate type for haplotype\n",nbuf);
		else if((element->type&ST_RANDOM) && (element->type&(ST_SEX|ST_GROUP|ST_STRING|ST_REAL|ST_INTEGER|ST_REAL)))
		  print_scan_err("Error: Variable %s inappropriate type to be random\n",nbuf);
		else if((element->type&(ST_INTTYPE|ST_REALTYPE)) && (element->type&(ST_STRING|ST_REAL|ST_INTEGER)))
		  print_scan_err("Error: Type collision for variable %s\n",nbuf);
		else if((element->type&ST_INTTYPE) && (element->type&ST_REALTYPE))
		  print_scan_err("Error: Real variable %s can not also be integer type\n",nbuf);
		else if((element->type&(ST_STRING|ST_REAL)) && (element->type&(ST_SEX|ST_ID|ST_FAMILY|ST_SIRE|ST_DAM)))
		  print_scan_err("Error: Variable %s can not be a pedigree or sex variable\n",nbuf);
		else if((element->type&ST_REAL) && (element->type&(ST_ID|ST_FAMILY|ST_SIRE|ST_DAM|ST_SEX)))
		  print_scan_err("Error: Real variable %s can not be a pedigree or sex variable\n",nbuf);
		else if((element->type&ST_FACTOR) && (element->type&ST_REAL))
		  print_scan_err("Error: Real variable %s can not be a factor\n",nbuf);
		else if((element->type&ST_CONSTANT)&&(element->type&ST_MULTIPLE))
		  print_scan_err("Error: Variable %s can not be in multiple records and be constant\n",nbuf);
		else if(element->type&(ST_SEX|ST_ID|ST_FAMILY|ST_SIRE|ST_DAM|ST_TRAITLOCUS|ST_GROUP|ST_LINKED|ST_MODEL|ST_TRAIT))
		  element->type|=ST_REQUIRED;
		else if(element->type&ST_HAPLO) {
			if(element->arg.element && element->arg.element->type&ST_LINKED) {
				element->type|=ST_REQUIRED;
				if(!(element->type&ST_DATA)) print_scan_err("Error: No data for variable %s\n",nbuf);
			}
		}
		if(element->type&ST_MARKER) n_markers++;
		if(!(element->type&(ST_CONSTANT|ST_MULTIPLE))) 
		  element->type|=syst_var[MULTIPLE_RECORDS]?ST_MULTIPLE:ST_CONSTANT;
		if(element->type&(ST_MARKER|ST_REQUIRED|ST_RESTRICT))	{
			if(!(element->type&(ST_HAPLO|ST_LUMPED))) element->arg.var=node;
		} else element->type=0;
	}
	free(nbuf);
	if(node->right) Check_var(node->right);
}

static struct bin_node *create_var(char *p)
{
	int k;
	struct bin_node *node;
	
	if(!root_var) node=root_var=alloc_var(p);
	else {
		root_var=find_var(p,root_var,&node,&k);
	}
	return node;
}

int symbol_lookup(char *p,int fg)
{
	static char *Coms[] = {"FILE","LOCUS","LOCI","MARKER","DISCRETE","MODEL","PEDIGREE","LOG",
		"FILTER","MISSING","MODEL","LINK","RANDOM","TRAIT","WHERE","USE",
		"REAL","INTEGER","SHELL","ARRAY","PRINT","DO","WHILE","CONSTANT",
		"MULTIPLE","CENSORED","GROUP","SET","SEX","AFFECTED","UNAFFECTED","PROBAND","OUTPUT","INCLUDE","ERRORDIR",
		"LAUROUTPUT","RAWOUTPUT","POSITION","FREQUENCY","SUPER",(char *)0};
	static int Com_token[] = {FILEC,LOCUS,LOCUS,MARKER,FACTOR,MODEL,PEDIGREE,LOG,
		FILTER,MISSING,MODEL,LINK,RANDOM,TRAIT,WHERE,USE,
		REAL,INTEGER,SHELL,ARRAY,PRINTEXP,DOLOOP,WHILE,CONSTANT,
		MULTIPLE,CENSORED,GROUP,SET,GENDER,AFFECTED,UNAFFECTED,PROBAND,OUTPUT,INCLUDE,ERRORDIR,
		LAUROUTPUT,RAWOUTPUT,POSITION,FREQUENCY,SUPER,SYSTEM_VAR,VARIABLE,ARRAY_VAR};
	static char *Syst[] = {"PRUNE_OPTION","RECODE_OPTION","NO_EXTRA_ALLELE",
		"PEEL_OPTION","TRACE_RESTRICT","TRACE_CENSORED","TRACE_AFFECTED",
		"CORRECT_ERRORS","TRACE_PEEL","MULTIPLE_RECORDS","MULTIVARIATE_TEST",
		"ERROR_CHECK","NO_DEFAULT_MISSING","SKIP_BAD_REALS","SKIP_BAD_INTS","IGNORE_CASE",(char *)0};
	int i=0,j=0;
	static struct scan_data *sd;
	
	while(Coms[i])	{
		if(!strcasecmp(Coms[i],p)) break;
		i++;
	}
	at_file=0;
	if(Com_token[i]==FILEC || Com_token[i]==LINK) at_file=1;
	if(Com_token[i]==SYSTEM_VAR) {
		i++;
		while(Syst[j])	{
			if(!strcasecmp(Syst[j],p))	{
				yylval.value=j;
				i--;
				break;
			}
			j++;
		}
	}
	if(Com_token[i]==VARIABLE) {
		if(fg==1 && begin_comm) {
			begin_comm=0;
			return BREAK;
		}
		yylval.var=create_var(p);
		sd=yylval.var->data;
		if(sd->vtype&ST_ARRAY) i++;
		if(fg==1) {
			begin_comm=1;
			(void)strncpy(linebuf1,linebuf,LINEBUFSIZE);
			lineno1=lineno;
		}
	} else if(begin_comm && Com_token[i]!=SYSTEM_VAR && Com_token[i]!=LOCUS && Com_token[i]!=SHELL
				 && !(at_use==1 && Com_token[i]==WHERE) && !(at_use==2 && Com_token[i]==USE))	{
		begin_comm=0;
		at_use=0;
		return BREAK;
	} else {
		begin_comm=1;
		(void)strncpy(linebuf1,linebuf,LINEBUFSIZE);
		lineno1=lineno;
		if(Com_token[i]==MODEL) at_model=1;
		else at_model=0;
		if(Com_token[i]==USE || Com_token[i]==CENSORED || Com_token[i]==AFFECTED || Com_token[i]==UNAFFECTED || Com_token[i]==PROBAND) at_use|=1;
		else if(Com_token[i]==WHERE) at_use|=2;
		else at_use=0;
	}
	return Com_token[i];
}

static struct var_list *add_to_var_list(struct var_list *vlist,struct bin_node *node,struct express *express)
{
	struct var_list *vlist1,*vlist2;
	struct scan_data *sd=0;
	int i;

	if(node)	sd=node->data;
	if(express) i=check_index(sd,express);
	else {
		i=0;
		if(sd && !(sd->vtype&ST_ARRAY)) sd->vtype|=ST_SCALAR;
	}
	vlist1=lk_malloc(sizeof(struct var_list));
	vlist1->next=0;
	vlist1->var=node;
	vlist1->index=i;
	vlist2=vlist;
	if(vlist2) {
		while(vlist2->next) vlist2=vlist2->next;
		vlist2->next=vlist1;
	} else vlist=vlist1;
	return vlist;
}

struct var_list *add_var_lists(struct var_list *vlist,struct var_list *vlist1)
{
	struct var_list *vlist2;
	
	vlist2=vlist;
	if(vlist2) {
		while(vlist2->next) vlist2=vlist2->next;
		vlist2->next=vlist1;
	} else vlist=vlist1;
	return vlist;
}

static void set_locus_array(struct bin_node *node)
{
	struct scan_data *sd;
	int i;
	
	sd=node->data;
	if(sd->vtype&ST_ARRAY) {
		for(i=0;i<sd->n_elements;i++) {
			set_locus_element(sd->element+i);
		}
	} else yyerror("Not an array");
}

static void set_locus_element(struct var_element *element)
{
	element->type|=(ST_MARKER|ST_FACTOR|ST_CONSTANT);
	if(hap_list[0]) {
		if(hap_list[0]->arg.element && hap_list[0]->arg.element!=element)	{
			yyerror1("Haplotype vector (left) used twice");
			hap_list[0]->arg.element=0;
		} else hap_list[0]->arg.element=element;
	}
	if(hap_list[1]) {
		if(hap_list[1]->arg.element && hap_list[1]->arg.element!=element)	{
			yyerror1("Haplotype vector (right) used twice");
			hap_list[1]->arg.element=0;
		} else hap_list[1]->arg.element=element;
	}
	hap_list[0]=hap_list[1]=0;
}

static void set_slocus_element(struct var_element *element,struct var_list *vlist)
{
	int j;
	struct scan_data *sd;
	
	element->type|=(ST_MARKER|ST_FACTOR|ST_CONSTANT);
	while(vlist) {
		sd=vlist->var->data;
		assert(sd);
		if(sd->vtype&ST_ARRAY) {
			if(vlist->index) {
				sd->element[vlist->index-1].type|=ST_LUMPED;
				sd->element[vlist->index-1].arg.element=element;
			} else for(j=0;j<sd->n_elements;j++) {
				sd->element[j].type|=ST_LUMPED;
				sd->element[j].arg.element=element;
			}
		} else {
			sd->element[0].type|=ST_LUMPED;
			sd->element[0].arg.element=element;
		}
		vlist=vlist->next;
	}
}

static void set_haplo_element(struct var_element *element,struct var_element *element1)
{
	element->type|=(ST_HAPLO|ST_FACTOR|ST_CONSTANT);
	if(element1) element1->type|=(ST_HAPLO|ST_FACTOR|ST_CONSTANT);
	hap_list[0]=element;
	hap_list[1]=element1;
}

static void do_file_com(char *fname,struct format_clause *fc,struct fformat *ff,struct var_list *vlist)
{
	int i,j;
	struct InFile *file;
	struct format *format;
	struct var_list *vlist1;
	struct var_element *element;
	struct scan_data *sd;
	
	if(!vlist) {
		yyerror1("No variables listed for FILE command\n");
		return;
	} else if(!fname) {
		free_vlist(vlist);
		return;
	} else if(!fname[0]) {
		yyerror1("Zero length filename for FILE command\n");
		free_vlist(vlist);
		return;
	}
	file=Infiles;
	Infiles=lk_calloc(1,sizeof(struct InFile));
	Infiles->next=file;
	Infiles->nvar=count_var_list(vlist);
	Infiles->element=lk_malloc(sizeof(void *)*Infiles->nvar);
	i=0;
	while(vlist) {
		if(vlist->var) {
			sd=vlist->var->data;
			if(sd->vtype&ST_ARRAY) {
				if(vlist->index) {
					element=sd->element+vlist->index-1;
					element->type|=ST_DATA;
					Infiles->element[i++]=element;
				} else {
					for(j=0;j<sd->n_elements;j++) {
						element=sd->element+j;
						element->type|=ST_DATA;
						Infiles->element[i++]=element;
					}
				}
			} else {
				element=sd->element;
				element->type|=ST_DATA;
				Infiles->element[i++]=element;
			}
		} else Infiles->element[i++]=0;
		vlist1=vlist->next;
		free(vlist);
		vlist=vlist1;
	}
	if(fc) {
		format=setup_format(fc);
		Infiles->format=format;
		if(!(scan_error&FORMAT_ERR)) {
			if(format->n_atoms<i) {
				(void)printf("format->n_atoms = %d\n",format->n_atoms);
				(void)printf("i = %d\n",i);
				print_scan_err("Line %d: Error - Too many variables for format clause\n",format->line);
				scan_error|=FORMAT_ERR;
			} else if(format->n_atoms>i)
				  print_scan_warn("Line %d: Warning - Too few variables for format clause\n",format->line);
		}
	} else if(ff) Infiles->fformat=ff;
	Infiles->name=fname;
	Infiles->shell_flag=shell_flag;
}

static void change_type(int type,struct var_list *vlist)
{
	int j;
	struct scan_data *sd;
	
	
	while(vlist) {
		sd=vlist->var->data;
		if(sd->vtype&ST_ARRAY) {
			if(vlist->index) sd->element[vlist->index-1].type|=type;
			else for(j=0;j<sd->n_elements;j++)
			  sd->element[j].type|=type;
		} else sd->element[0].type|=type;
		vlist=vlist->next;
	}
}

static void free_vlist(struct var_list *vlist)
{
	struct var_list *vlist1;
	
	while(vlist) {
		vlist1=vlist->next;
		free(vlist);
		vlist=vlist1;
	}
}

static void do_link_com(char *s,int type,struct var_list *vlist)
{
	struct Link *l,*l1,**ll;
	struct var_list *vlist1;
	struct var_element *element;
	struct scan_data *sd=0;
	int i,j,k;
	
	if(vlist) sd=vlist->var->data;
	if(!s && sd) {
		if(sd->vtype&ST_ARRAY && vlist->index) {
			element=sd->element+vlist->index-1;
			if(element->type&ST_STRING) {
				s=element->arg.string;
				vlist1=vlist->next;
				free(vlist);
				vlist=vlist1;
			}
		} else {
			element=sd->element;
			if(element->type&ST_STRING) {
				s=element->arg.string;
				vlist1=vlist->next;
				free(vlist);
				vlist=vlist1;
			}
		}
	}
	ll=&links;
	while(*ll) {
		l=*ll;
		if(s) {
			if(l->name) {
				if(!strcasecmp(s,l->name)) break;
			}
		} else if(!l->name) break;
		ll=&l->next;
	}
	if(*ll) l1=*ll;
	else {
		l1=lk_malloc(sizeof(struct Link));
		l1->next=0;
		l1->name=s;
		l1->n_loci=0;
		l1->element=0;
		l1->type=-1;
		*ll=l1;
	}
	i=count_var_list(vlist);
	if(l1->type>=0 && l1->type!=type) print_scan_err("Error: Linkage group has inconsistent linkage type\n");
	l1->type=type;
	if(i) {
		k=i+l1->n_loci;
		if(l1->element) {
			l1->element=lk_realloc(l1->element,sizeof(void *)*k);
		} else l1->element=lk_malloc(sizeof(void *)*k);
		while(vlist) {
			sd=vlist->var->data;
			if(sd->vtype&ST_ARRAY) {
				if(vlist->index) {
					element=sd->element+vlist->index-1;
					for(i=0;i<l1->n_loci;i++) {
						if(l1->element[i]==element) break;
					}
					if(i==l1->n_loci) {
						if(element->type&ST_LINKED) {
							print_scan_err("Error: %s(%d) appears in multiple linkage groups\n",sd->name,vlist->index);
							scan_error|=LINK_ERR;
						} else {
							element->type|=ST_LINKED;
							l1->element[l1->n_loci++]=element;
						}
					}
				} else {
					for(j=0;j<sd->n_elements;j++) {
						element=sd->element+j;
						for(i=0;i<l1->n_loci;i++) {
							if(l1->element[i]==element) break;
						}
						if(i==l1->n_loci) {
							if(element->type&ST_LINKED) {
								print_scan_err("Error: %s(%d) appears in multiple linkage groups\n",sd->name,vlist->index);
								scan_error|=LINK_ERR;
							} else {
								element->type|=ST_LINKED;
								l1->element[l1->n_loci++]=element;
							}
						}
						sd->vtype|=ST_LINKED;
					}
				}
			} else {
				element=sd->element;
				for(i=0;i<l1->n_loci;i++) {
					if(l1->element[i]==element) break;
				}
				if(i==l1->n_loci) {
					if(element->type&ST_LINKED) {
						print_scan_err("Error: %s appears in multiple linkage groups\n",sd->name);
						scan_error|=LINK_ERR;
					} else {
						element->type|=ST_LINKED;
						l1->element[l1->n_loci++]=element;
					}
				}
			}
			vlist1=vlist->next;
			free(vlist);
			vlist=vlist1;
		}
		if(l1->n_loci<k) l1->element=lk_realloc(l1->element,sizeof(void *)*l1->n_loci);
	}
}

static void do_missing_com(struct express *expr,struct var_list *vlist,char *s1)
{
	struct var_list *vlist1;
	struct scan_data *sd;
	struct var_element **elem;
	struct Miss *m;
	int i,j;
	char *p;
	
	if(s1) {
		assert(!vlist);
		if(s1[0]==0) {
			print_scan_err("Empty scope - MISSING directive ignored\n");
			if(expr->type==ST_STRING) free(expr->arg.string);
			free(s1);
			return;
		}
		qstrip(s1);
		p=s1;
		i=j=0;
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
		if(*p) {
			j=1;
			print_scan_err("Illegal character '%c' in MISSING scope\n",*p);
		} else if(*(--p)=='!') {
			j=1;
			print_scan_err("MISSING scope can not end with a '!'\n",*p);
		}
		if(j) {
			free(s1);
			if(expr->type==ST_STRING) free(expr->arg.string);
			return;
		}
	}
	m=Miss;
	Miss=lk_malloc(sizeof(struct Miss));
	Miss->Missing.arg=expr->arg;
	Miss->Missing.type=expr->type;
	Miss->next=m;
	Miss->element=0;
	Miss->scope=0;
	if((i=count_var_list(vlist))) {
		elem=lk_malloc(sizeof(void *)*i);
		i=0;
		while(vlist) {
			sd=vlist->var->data;
			if(sd->vtype&ST_ARRAY) {
				if(vlist->index) elem[i++]=sd->element+vlist->index-1;
				else for(j=0;j<sd->n_elements;j++) elem[i++]=sd->element+j;
			} else elem[i++]=sd->element;
			Miss->element=elem;
			vlist1=vlist->next;
			free(vlist);
			vlist=vlist1;
		}
	} else if(s1) Miss->scope=s1;
	Miss->nvar=i;
}

static void do_model_com(struct model_list *mlist,struct bin_node *node,struct express *express)
{
	struct model *model,*model1;
	struct var_element *element;
	struct scan_data *sd;
	
	sd=node->data;
	model=lk_malloc(sizeof(struct model));
	model->next=0;
   if(Models) {
		model1=Models;
		while(model1->next) model1=model1->next;
		model1->next=model;
	} else Models=model;
	model->trait=sd;
	if(!express) {
		model->index=0;
		sd->element[0].type|=ST_TRAIT;
	} else {
		element=get_element(node,express);
		if(element) {
			model->index=(int)express->arg.value;
			element->index=element->oindex=model->index;
			element->type|=ST_TRAIT;
		} else model->trait=0;
	}
	model->model_list=mlist;
}

static void add_censored(struct var_element *element,const int fg)
{
	struct operation *ops;
	struct Censor *cen;
	
	if(fg==1 && !element) {
		print_scan_err("Error: Nothing to censor!\n");
		return;
	}
	/* Reverse list order (really return list to original order! */
	ops=Op_List=reverse_list(Op_List);
	switch(fg) {
	 case 1:
		cen=lk_malloc(sizeof(struct Censor));
		cen->next=Censored;
		Censored=cen;
		cen->Op_List=ops;
		cen->element=element;
		element->type|=ST_CENSORED;
		break;
	 case 0:
		if(Affected) {
			print_scan_warn("Warning - new affected statement overrules previous statement\n");
			free_op(Affected);
		}
		Affected=ops;
		break;
	 case 2:
		if(Unaffected) {
			print_scan_warn("Warning - new unaffected statement overrules previous statement\n");
			free_op(Unaffected);
		}
		Unaffected=ops;
		break;
	 case 3:
		if(Proband) {
			print_scan_warn("Warning - new proband statement overrules previous statement\n");
			free_op(Proband);
		}
		Proband=ops;
		break;
	}
	ops=Op_List;
	while(ops) {
		if(ops->type==VARIABLE) ops->arg.element->type|=ST_RESTRICT;
		ops=ops->next;
	}
}

static void add_restriction(struct var_list *vlist)
{
	struct operation *ops;
	struct Restrict *res;
	struct var_list *vlist1;
	struct scan_data *sd;
	int i,j;
	
	/* Reverse list order (really return list to original order! */
	Op_List=ops=reverse_list(Op_List);
	res=lk_malloc(sizeof(struct Restrict));
	res->next=Restrictions;
	Restrictions=res;
	res->Op_List=ops;
	if((res->nvar=count_var_list(vlist))) {
		res->element=lk_malloc(sizeof(void *)*res->nvar);
		i=0;
		while(vlist) {
			sd=vlist->var->data;
			if(sd->vtype&ST_ARRAY) {
				if(vlist->index) res->element[i++]=sd->element+vlist->index-1;
				else for(j=0;j<sd->n_elements;j++) res->element[i++]=sd->element+j;
			} else res->element[i++]=sd->element;
			vlist1=vlist->next;
			free(vlist);
			vlist=vlist1;
		}	
	} else res->element=0;
	while(ops) {
		if(ops->type==VARIABLE) ops->arg.element->type|=ST_RESTRICT;
		ops=ops->next;
	}
}

static void find_markers(struct bin_node *node,int *i)
{
	int j;
	struct scan_data *sd;
	
	sd=node->data;
	for(j=0;j<sd->n_elements;j++) if(sd->element[j].type&ST_MARKER) {
		if(sd->element[j].type&ST_REQUIRED) {
			markers[*i].element=sd->element+j;
			markers[*i].var=sd;
			markers[(*i)++].index=sd->n_elements>1?j+1:0;
		}
	}
}

static void find_trait_loci(struct bin_node *node,int *i)
{
	int j;
	struct scan_data *sd;
	
	sd=node->data;
	for(j=0;j<sd->n_elements;j++) if(sd->element[j].type&ST_TRAITLOCUS) {
		if(traitlocus) {
			traitlocus[*i].element=sd->element+j;
			traitlocus[*i].var=sd;
			traitlocus[(*i)].index=j+1;
		}
		(*i)++;
	}
}

static void handle_super_loci(struct bin_node *node)
{
	int k;
	struct scan_data *sd;
	
	sd=node->data;
	for(k=0;k<sd->n_elements;k++) if(sd->element[k].type&ST_LUMPED) {
		if(sd->element[k].arg.element->type&ST_LINKED) {
			if(!(sd->element[k].type&ST_MARKER)) {
				if(sd->n_elements>1) print_scan_err("Error: Variable '%s(%d)' is not a marker and so can not be part of a super locus\n",sd->name,k+1);
				else print_scan_err("Error: Variable '%s' is not a marker and so can not be part of a super locus\n",sd->name);
			} else sd->element[k].type|=ST_LINKED;
		}
	}
}

static void find_haplo(struct bin_node *node)
{
	int j,k;
	struct scan_data *sd;
	
	sd=node->data;
	for(k=0;k<sd->n_elements;k++) if(sd->element[k].type&ST_HAPLO) {
		for(j=0;j<n_markers;j++) if(sd->element[k].arg.element==markers[j].element) {
			if(!markers[j].hap_element[0]) markers[j].hap_element[0]=sd->element+k;
			else if(!markers[j].hap_element[1])	markers[j].hap_element[1]=sd->element+k;
			else {
				if(markers[j].index) print_scan_err("Error: marker %s(%d) has >2 haplotype vectors associated with it\n",markers[j].var->name,markers[j].index);
				else print_scan_err("Error: marker %s has >2 haplotype vectors associated with it\n",markers[j].var->name);
			}
			break;
		}
		assert(j<n_markers);
	}
}

static void find_lumped1(struct bin_node *node)
{
	int j,k;
	struct scan_data *sd;
	
	sd=node->data;
	for(k=0;k<sd->n_elements;k++) if(sd->element[k].type&ST_LUMPED) {
		if(sd->element[k].type&ST_LINKED) {
			for(j=0;j<n_markers;j++) if(sd->element[k].arg.element==markers[j].element) {
				if(markers[j].hap_element[0]) {
					if(markers[j].index) print_scan_err("Error: super locus %s(%d) has haplotype vectors associated with it\n",markers[j].var->name,markers[j].index);
					else print_scan_err("Error: super locus %s has haplotype vectors associated with it\n",markers[j].var->name);
				}
				markers[j].n_sub_elements++;
				break;
			}
			assert(j<n_markers);
		}
	}
}

static void find_lumped2(struct bin_node *node)
{
	int j,k;
	struct scan_data *sd;
	
	sd=node->data;
	for(k=0;k<sd->n_elements;k++) if(sd->element[k].type&ST_LUMPED) {
		if(sd->element[k].type&ST_LINKED) {
			for(j=0;j<n_markers;j++) if(sd->element[k].arg.element==markers[j].element) {
				markers[j].sub_element[markers[j].n_sub_elements++]=sd->element+k;
				break;
			}
		}
	}
}

static void correct_linkage(struct bin_node *node)
{
	int j,k,k1;
	struct scan_data *sd;
	char *p1,*p2;
	
	sd=node->data;
	for(k=0;k<sd->n_elements;k++) if(sd->element[k].type&ST_LUMPED) {
		if(sd->element[k].type&ST_LINKED) {
			for(j=0;j<n_markers;j++) if(sd->element+k==markers[j].element) break;
			assert(j<n_markers);
			for(k1=0;k1<n_markers;k1++) if(sd->element[k].arg.element==markers[k1].element) {
				break;
			}
			assert(k1<n_markers);
			if(markers[j].link) {
				p1=get_marker_name(j);
				p2=get_marker_name(k1);
				print_scan_err("Error: Marker '%s' is a part of superlocus '%s' so can not appear directly in a LINK statement\n",p1,p2);
				free(p1);
				free(p2);
			}
			markers[j].link=markers[k1].link;
		}
	}
}

static void strip_names(struct bin_node *node)
{
	char *p;
	int i;
	struct scan_data *sd;
	
	sd=node->data;
	if((p=sd->name)) {
		i=strlen(p);
		if(i>2) {
			if(p[i-1]=='_' && p[0]=='_') {
				p[i-1]=0;
				memmove(p,p+1,i-1);
			}
		}
	}
}

int ReadControl(FILE *fptr,char *cname,char **lfile)
{
	int i,j,k;
	void yy_cleanup(void);
	struct InFile *infile,**infile_p;
	struct Restrict *res,*res1,**res_p;
	struct Censor *cen,**cen_p;
	struct var_element *elem;
	struct Link *linkp;
	struct operation *ops;
	struct express tmp_expr;
	struct marker_info *mi;
	
	yyin=fptr;
	fname_list[0]=cname;
	list_ptr=0;
	f_atom_list=lk_malloc(sizeof(struct format_atom)*f_atom_size);
	for(i=0;i<NUM_SYSTEM_VAR;i++) syst_var[i]=0;
	syst_var[PRUNE_OPTION]=syst_var[RECODE_OPTION]=2;
	syst_var[ERROR_CHECK]=1;
	if((i=yyparse())) print_scan_err("Error: yyparse returned error %d\n",i);
	yy_cleanup();
	if(strip_vars) check_vars_1(root_var,strip_names);
	/* Sanity check! */
	if(!scan_error_n)	{
		if(!Infiles) print_scan_err("Error: No input files specified\n");
		if(!pedflag) print_scan_err("Error: No pedigree variables specified\n");
		else {
			for(i=0;i<3;i++) if(pedlist[i+family_id]->type&ST_INTTYPE) break;
			if(i<3) for(i=0;i<3;i++) pedlist[i+family_id]->type|=ST_INTTYPE;
		}
		if(root_var) {
			check_vars_1(root_var,handle_super_loci);
			Check_var(root_var);
		}
		/* Flag variables used as the operands to a restriction statement *whose result is used* as ST_REQUIRED */
		res=0;
		while(res!=Restrictions) {	
			res1=Restrictions;
			while(res1->next!=res) res1=res1->next;
			for(i=j=0;i<res1->nvar;i++) if(res1->element[i]->type&ST_REQUIRED) {
				j=1;
				break;
			}
			if(!res1->nvar || j) {
				ops=res1->Op_List;
				while(ops) {
					if(ops->type==VARIABLE) ops->arg.element->type|=ST_REQUIRED;
					ops=ops->next;
				}
			}
			res=res1;
		}
		/* Delete restrict structures that are not used */
		res=Restrictions;
		res_p= &Restrictions;
		while(res) {
			for(i=j=0;i<res->nvar;i++) if(res->element[i]->type&ST_REQUIRED) {
				j=1;
				break;
			}
			if(res->nvar && !j) {
				*res_p=res->next;
				free_restrict(res);
				res= *res_p;
			} else {
				res_p= &res->next;
				res=res->next;
			}
		}
		if(Unaffected && !Affected) print_scan_err("Error: Unaffected definition without affected definition\n");
		/* Flag variables used in censored statements as required.  Delete unused censored statements */
		cen=Censored;
		cen_p= &Censored;
		while(cen) {
			if(cen->element->type&ST_TRAIT) {
				ops=cen->Op_List;
				while(ops) {
					if(ops->type==VARIABLE) ops->arg.element->type|=ST_REQUIRED;
					ops=ops->next;
				}
				cen_p= &cen->next;
				cen=cen->next;
			} else {
				*cen_p=cen->next;
				free_op(cen->Op_List);
				free(cen);
				cen= *cen_p;
			}
		}
		if((ops=Affected)) {
			while(ops) {
				if(ops->type==VARIABLE) ops->arg.element->type|=ST_REQUIRED;
				ops=ops->next;
			}
		}
		if((ops=Unaffected)) {
			while(ops) {
				if(ops->type==VARIABLE) ops->arg.element->type|=ST_REQUIRED;
				ops=ops->next;
			}
		}
		if((ops=Proband)) {
			while(ops) {
				if(ops->type==VARIABLE) ops->arg.element->type|=ST_REQUIRED;
				ops=ops->next;
			}
		}
		/* Check file structures - remove ones that aren't needed */
		infile=Infiles;
		infile_p= &Infiles;
		while(infile) {
			infile->ncol=0;
			for(i=0;i<infile->nvar;i++) {
				elem=infile->element[i];
			   if(elem) elem->type&=~ST_FLAG;
			}
			for(k=infile->nvar-1;k>=0;k--) {
				elem=infile->element[k];
				if(elem && (elem->type&ST_REQUIRED)) break;
			}
			for(i=j=0;i<infile->nvar;i++)	{
				elem=infile->element[i];
				if(elem) {
					if((elem->type&ST_MARKER) && !(elem->type&ST_REQUIRED)) {
						if(j<k) elem->type|=(ST_REQUIRED|ST_NOT_REALLY_REQUIRED);
						else elem->type=0;
					}
					if(elem->type&ST_REQUIRED) {
						if(elem->type&ST_FLAG) {
							print_scan_err("Error: Duplicate variables for file %s\n",infile->name);
							break;
						}
						elem->type|=ST_FLAG;
						if(elem->type&ST_ID)	{
							j|=1;
							infile->id_col=infile->ncol;
						} else if(elem->type&ST_FAMILY) {
							j|=2;
							infile->family_col=infile->ncol;
						}
						infile->ncol++;
					} else infile->element[i]=0;
				}
			}
			for(i=0;i<infile->nvar;i++) if(infile->element[i]) infile->element[i]->type&=~ST_FLAG;
			if(!(j&1)) print_scan_err("Error: No id column for file %s\n",infile->name);
			else if(family_id && j!=3) print_scan_err("Error: No family column for file %s\n",infile->name);
			if(infile->ncol==1) {
				*infile_p=infile->next;
				free_infile(infile);
				infile= *infile_p;
			} else {
				infile_p= &infile->next;
				infile=infile->next;
			}
		}
		if(!Infiles) print_scan_err("Error: No input files with data\n");
		free(f_atom_list);
		/* Count markers and link up with haplotype vectors */
		if(n_markers) {
			markers=lk_calloc((size_t)n_markers,sizeof(struct Marker));
			for(i=0;i<n_markers;i++) {
				markers[i].allele_trans=0;
			   markers[i].order=0;
				markers[i].o_size=0;
				markers[i].pos_set[0]=markers[i].pos_set[1]=0;
				markers[i].hap_element[0]=markers[i].hap_element[1];
				markers[i].sub_element=0;
				markers[i].n_sub_elements=0;
			}
			i=0;
			if(root_var) {
				check_vars(root_var,&i,find_markers);
				check_vars_1(root_var,find_haplo);
				check_vars_1(root_var,find_lumped1);
				for(i=0;i<n_markers;i++) if(markers[i].n_sub_elements) {
					markers[i].sub_element=lk_malloc(sizeof(void *)*markers[i].n_sub_elements);
					markers[i].n_sub_elements=0;
				}
				check_vars_1(root_var,find_lumped2);
			}
			n_markers=i;
			while(m_info) {
				mi=m_info->next;
				for(i=0;i<n_markers;i++) if(markers[i].element==m_info->element) {
					for(k=0;k<3;k++) {
						markers[i].pos_set[k]=m_info->pos_set[k];
						markers[i].pos[k]=m_info->pos[k];
					}
					break;
				}
				free(m_info);
				m_info=mi;
			}
			for(i=0;i<n_markers;i++) {
				if(!markers[i].element || markers[i].element->type&ST_NOT_REALLY_REQUIRED) continue;
				linkp=0;
				j=0;
				linkp=links;
				while(linkp) {
					j++;
					for(k=0;k<linkp->n_loci;k++) {
						if(linkp->element[k]==markers[i].element) {
							markers[i].link=j;
							break;
						}
					}
					if(k<linkp->n_loci) break;
					linkp=linkp->next;
				}
				if(!linkp && !(markers[i].element->type&ST_LUMPED)) {
					if(markers[i].var->vtype&ST_ARRAY)
					  abt(__FILE__,__LINE__,"%s(): No linkage group specified for candidate gene %s(%d)\n",__func__,markers[i].var->name,markers[i].index);
					else abt(__FILE__,__LINE__,"%s(): No linkage group specified for candidate gene %s\n",__func__,markers[i].var->name);
				}
				if(markers[i].n_sub_elements) {
					assert(!markers[i].hap_element[0]);
					if(markers[i].element->type&ST_DATA) {
						if(markers[i].var->vtype&ST_ARRAY)
						  print_scan_err("Error: super locus %s(%d) can not have its own data\n",markers[i].var->name,markers[i].index);
						else
						  print_scan_err("Error: super locus %s can not have its own data\n",markers[i].var->name);
					}
				} else if(markers[i].hap_element[0]) {
					if(markers[i].element->type&ST_DATA) {
						if(markers[i].var->vtype&ST_ARRAY)
						  print_scan_err("Error: marker variable %s(%d) can not have both genotype and haplotype data\n",markers[i].var->name,markers[i].index);
						else
						  print_scan_err("Error: marker variable %s can not have both genotype and haplotype data\n",markers[i].var->name);
					}
					if(markers[i].hap_element[0]->type&ST_INTTYPE) markers[i].hap_element[1]->type|=ST_INTTYPE;
					if(markers[i].hap_element[1] && markers[i].hap_element[1]->type&ST_INTTYPE) markers[i].hap_element[0]->type|=ST_INTTYPE;
				} else {
					if(!(markers[i].element->type&ST_DATA)) {
						if(markers[i].var->vtype&ST_ARRAY)
						  print_scan_err("Error: marker variable %s(%d) has no data\n",markers[i].var->name,markers[i].index);
						else
						  print_scan_err("Error: marker variable %s has no data\n",markers[i].var->name);
					}
				}
			}
		}
		i=0;
		if(root_var) {
			check_vars_1(root_var,correct_linkage);
			check_vars(root_var,&i,find_trait_loci);
		}
		if(i) {
			if(i>1) print_scan_err("Error: multiple trait loci indicated\n");
			else {
				traitlocus=lk_calloc(1,sizeof(struct Marker));
				traitlocus->sub_element=0;
				traitlocus->n_sub_elements=0;
				traitlocus->order=0;
				traitlocus->o_size=0;
				i=0;
				check_vars(root_var,&i,find_trait_loci);
			}
		}
		if(Models && Models->next && !syst_var[MULTIVARIATE_TEST]) {
			print_scan_err("Error: Multiple models not currently supported\n");
		}
	}
	*lfile=LogFile;
	if(!scan_error_n && !Miss && !syst_var[NO_DEFAULT_MISSING]) {
		tmp_expr.arg.string=strdup("0");
		tmp_expr.type=ST_STRING;
		do_missing_com(&tmp_expr,0,strdup("PF"));
	}
	return scan_error_n;
}

