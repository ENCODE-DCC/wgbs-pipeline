
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
#define yyparse         xmlparse
#define yylex           xmllex
#define yyerror         xmlerror
#define yylval          xmllval
#define yychar          xmlchar
#define yydebug         xmldebug
#define yynerrs         xmlnerrs


/* Copy the first part of user declarations.  */

/* Line 189 of yacc.c  */
#line 2 "xml_parse.y"

  /****************************************************************************
   *                                                                          *
   *     Loki - Programs for genetic analysis of complex traits using MCMC    *
   *                                                                          *
   *                     Simon Heath - CNG, Evry                              *
   *                                                                          *
   *                         February 2003                                    *
   *                                                                          *
   * xml_parse.y:                                                             *
   *                                                                          *
   * Basic XML parser                                                         *
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
#include <ctype.h>

#include "string_utils.h"
#include "xml.h"
#include "lk_malloc.h"
#include "y.tab.h"

#undef YYLSP_NEEDED
#undef YYLEX_PARAM

  static int xmllex(yystype *);
  static int xmlerror(char *);
  static FILE *xmlin;
  static int *line_no,col_pos,abt_flag;
  static args *make_arg(string *,string *);
  static void free_att_list(args *);
  static void check_start(string *,args *);
  static void check_end(string *,int);
  static void check_declaration(args *);
  static int check_ws(string *);
  static int check_dtd_tok(string *,char *);
  static XML_handler *call;
  
  typedef struct state {
    struct state *next;
    struct state *prev;
    string *name;
  } state;

  static state *curr_state,*free_state_list;
  static int pi_state,in_comment,element_state,attlist_state;

  

/* Line 189 of yacc.c  */
#line 142 "y.tab.c"

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
     END = 258,
     XMLTOK = 259,
     QCLOSE = 260,
     QOPEN = 261,
     DOCTYPE = 262,
     COMMENTSTART = 263,
     ELEMENT = 264,
     ENTITY = 265,
     ATTLIST = 266,
     NOTATION = 267,
     CDATASTART = 268,
     DOUBLEDASH = 269,
     EMPTY = 270,
     ANY = 271,
     PCDATA = 272,
     REQUIRED = 273,
     IMPLIED = 274,
     FIXED = 275,
     CDATA = 276,
     ID = 277,
     IDREF = 278,
     IDREFS = 279,
     ENTITIES = 280,
     NMTOKEN = 281,
     NMTOKENS = 282,
     LETTER = 283,
     DIGIT = 284,
     SQ = 285,
     DQ = 286,
     DASH = 287,
     DOT = 288,
     USCORE = 289,
     COLON = 290,
     CHAR = 291,
     OPEN = 292,
     CLOSE = 293,
     SLASH = 294,
     AMPERSAND = 295,
     EQ = 296,
     LSBRACK = 297,
     RSBRACK = 298,
     PERCENT = 299,
     SEMICOLON = 300,
     LBRACK = 301,
     RBRACK = 302,
     STAR = 303,
     PLUS = 304,
     COMMA = 305,
     QUERY = 306,
     BAR = 307,
     SPACE = 308
   };
#endif
/* Tokens.  */
#define END 258
#define XMLTOK 259
#define QCLOSE 260
#define QOPEN 261
#define DOCTYPE 262
#define COMMENTSTART 263
#define ELEMENT 264
#define ENTITY 265
#define ATTLIST 266
#define NOTATION 267
#define CDATASTART 268
#define DOUBLEDASH 269
#define EMPTY 270
#define ANY 271
#define PCDATA 272
#define REQUIRED 273
#define IMPLIED 274
#define FIXED 275
#define CDATA 276
#define ID 277
#define IDREF 278
#define IDREFS 279
#define ENTITIES 280
#define NMTOKEN 281
#define NMTOKENS 282
#define LETTER 283
#define DIGIT 284
#define SQ 285
#define DQ 286
#define DASH 287
#define DOT 288
#define USCORE 289
#define COLON 290
#define CHAR 291
#define OPEN 292
#define CLOSE 293
#define SLASH 294
#define AMPERSAND 295
#define EQ 296
#define LSBRACK 297
#define RSBRACK 298
#define PERCENT 299
#define SEMICOLON 300
#define LBRACK 301
#define RBRACK 302
#define STAR 303
#define PLUS 304
#define COMMA 305
#define QUERY 306
#define BAR 307
#define SPACE 308




#if ! defined YYSTYPE && ! defined YYSTYPE_IS_DECLARED
typedef union YYSTYPE
{

/* Line 214 of yacc.c  */
#line 62 "xml_parse.y"

  string *str;
  args *arg;
  char c;
  int i;



/* Line 214 of yacc.c  */
#line 293 "y.tab.c"
} YYSTYPE;
# define YYSTYPE_IS_TRIVIAL 1
# define yystype YYSTYPE /* obsolescent; will be withdrawn */
# define YYSTYPE_IS_DECLARED 1
#endif


/* Copy the second part of user declarations.  */


/* Line 264 of yacc.c  */
#line 305 "y.tab.c"

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
#define YYFINAL  14
/* YYLAST -- Last index in YYTABLE.  */
#define YYLAST   661

/* YYNTOKENS -- Number of terminals.  */
#define YYNTOKENS  54
/* YYNNTS -- Number of nonterminals.  */
#define YYNNTS  63
/* YYNRULES -- Number of rules.  */
#define YYNRULES  161
/* YYNRULES -- Number of states.  */
#define YYNSTATES  277

/* YYTRANSLATE(YYLEX) -- Bison symbol number corresponding to YYLEX.  */
#define YYUNDEFTOK  2
#define YYMAXUTOK   308

#define YYTRANSLATE(YYX)						\
  ((unsigned int) (YYX) <= YYMAXUTOK ? yytranslate[YYX] : YYUNDEFTOK)

/* YYTRANSLATE[YYLEX] -- Bison symbol number corresponding to YYLEX.  */
static const yytype_uint8 yytranslate[] =
{
       0,     2,     2,     2,     2,     2,     2,     2,     2,     2,
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
      45,    46,    47,    48,    49,    50,    51,    52,    53
};

#if YYDEBUG
/* YYPRHS[YYN] -- Index of the first RHS symbol of rule number YYN in
   YYRHS.  */
static const yytype_uint16 yyprhs[] =
{
       0,     0,     3,     7,    10,    15,    17,    18,    20,    22,
      25,    27,    30,    34,    37,    40,    46,    53,    56,    63,
      72,    76,    78,    81,    83,    85,    87,    89,    93,    94,
     103,   110,   116,   118,   120,   121,   122,   131,   132,   133,
     143,   145,   147,   149,   151,   153,   155,   157,   159,   161,
     168,   174,   179,   185,   187,   189,   193,   195,   196,   198,
     200,   202,   205,   212,   220,   225,   231,   232,   234,   236,
     238,   245,   252,   258,   263,   269,   274,   280,   283,   286,
     288,   294,   300,   305,   307,   309,   312,   315,   317,   320,
     322,   324,   326,   328,   330,   332,   334,   336,   338,   340,
     342,   344,   346,   348,   350,   352,   354,   356,   358,   360,
     362,   364,   366,   368,   370,   372,   374,   377,   379,   381,
     383,   385,   387,   390,   393,   394,   400,   402,   405,   409,
     412,   416,   420,   422,   425,   430,   434,   438,   440,   442,
     445,   448,   450,   452,   455,   458,   459,   461,   466,   471,
     473,   475,   478,   481,   483,   485,   488,   490,   493,   495,
     498,   500
};

/* YYRHS -- A `-1'-separated list of the rules' RHS.  */
static const yytype_int8 yyrhs[] =
{
      55,     0,    -1,    56,   103,    99,    -1,    87,    57,    -1,
      87,    57,    60,    57,    -1,    57,    -1,    -1,    58,    -1,
      99,    -1,    58,    99,    -1,    41,    -1,    53,    41,    -1,
      53,    41,    53,    -1,    41,    53,    -1,    61,    38,    -1,
      61,    42,    43,   111,    38,    -1,    61,    42,    63,    43,
     111,    38,    -1,    62,   111,    -1,    62,    53,    90,    53,
     108,   111,    -1,    62,    53,    90,    53,   108,    53,   108,
     111,    -1,     7,    53,    98,    -1,    64,    -1,    63,    64,
      -1,    67,    -1,    65,    -1,    66,    -1,    53,    -1,    44,
      98,    45,    -1,    -1,     9,    53,    98,    53,    68,    79,
     111,    38,    -1,    11,    53,    98,    69,   111,    38,    -1,
      11,    53,    98,   111,    38,    -1,    88,    -1,   101,    -1,
      -1,    -1,    53,    98,    53,    70,    74,    53,    71,    77,
      -1,    -1,    -1,    69,    53,    98,    53,    72,    74,    53,
      73,    77,    -1,    21,    -1,    22,    -1,    23,    -1,    24,
      -1,    10,    -1,    25,    -1,    26,    -1,    27,    -1,    75,
      -1,    46,   111,    98,    76,   111,    47,    -1,    46,   111,
      98,   111,    47,    -1,   111,    52,   111,    98,    -1,    76,
     111,    52,   111,    98,    -1,    18,    -1,    19,    -1,    20,
      53,   108,    -1,   108,    -1,    -1,    48,    -1,    15,    -1,
      16,    -1,    82,    81,    -1,    46,   111,    17,   111,    47,
      78,    -1,    46,   111,    17,    80,   111,    47,    48,    -1,
     111,    52,   111,    98,    -1,    80,   111,    52,   111,    98,
      -1,    -1,    49,    -1,    48,    -1,    51,    -1,    46,   111,
      85,    83,   111,    47,    -1,    46,   111,    85,    84,   111,
      47,    -1,    46,   111,    85,   111,    47,    -1,   111,    52,
     111,    85,    -1,    83,   111,    52,   111,    85,    -1,   111,
      50,   111,    85,    -1,    84,   111,    50,   111,    85,    -1,
      98,    81,    -1,    82,    81,    -1,     6,    -1,    86,     4,
     106,   111,     5,    -1,    86,    98,    53,    89,     5,    -1,
      86,    98,   111,     5,    -1,    91,    -1,    53,    -1,    89,
      91,    -1,    89,    53,    -1,    28,    -1,    90,    28,    -1,
      92,    -1,    37,    -1,    40,    -1,    93,    -1,    30,    -1,
      31,    -1,    96,    -1,    36,    -1,    42,    -1,    43,    -1,
      45,    -1,    44,    -1,    38,    -1,    41,    -1,    39,    -1,
      93,    -1,    31,    -1,    93,    -1,    30,    -1,    97,    -1,
      32,    -1,    33,    -1,    28,    -1,    29,    -1,    34,    -1,
      35,    -1,    97,    -1,    98,    96,    -1,   101,    -1,    88,
      -1,    53,    -1,    91,    -1,    53,    -1,   100,    91,    -1,
     100,    53,    -1,    -1,     8,   102,   100,    14,    38,    -1,
     113,    -1,   105,   112,    -1,   105,   115,   112,    -1,    37,
      98,    -1,   104,   106,    38,    -1,   104,   111,    38,    -1,
     107,    -1,   106,   107,    -1,    53,    98,    59,   108,    -1,
      30,   109,    30,    -1,    31,   110,    31,    -1,    94,    -1,
      53,    -1,   109,    94,    -1,   109,    53,    -1,    95,    -1,
      53,    -1,   110,    95,    -1,   110,    53,    -1,    -1,    53,
      -1,     3,    98,   111,    38,    -1,   104,   106,    39,    38,
      -1,    92,    -1,    53,    -1,   114,    92,    -1,   114,    53,
      -1,   114,    -1,   116,    -1,   115,   116,    -1,   103,    -1,
     103,   114,    -1,   101,    -1,   101,   114,    -1,    88,    -1,
      88,   114,    -1
};

/* YYRLINE[YYN] -- source line where rule number YYN was defined.  */
static const yytype_uint16 yyrline[] =
{
       0,    81,    81,    83,    84,    85,    87,    88,    90,    91,
      93,    94,    95,    96,    98,    99,   100,   102,   103,   104,
     106,   108,   109,   111,   112,   114,   115,   117,   119,   119,
     120,   121,   122,   123,   125,   125,   125,   126,   126,   126,
     128,   129,   130,   131,   132,   133,   134,   135,   136,   138,
     139,   141,   142,   144,   145,   146,   147,   149,   150,   152,
     153,   154,   155,   156,   158,   159,   161,   162,   163,   164,
     166,   167,   168,   170,   171,   173,   174,   176,   177,   179,
     181,   183,   184,   186,   187,   188,   189,   191,   192,   194,
     195,   196,   198,   199,   200,   202,   203,   204,   205,   206,
     207,   208,   209,   210,   212,   213,   215,   216,   218,   219,
     220,   222,   223,   224,   225,   227,   228,   230,   231,   232,
     234,   235,   236,   237,   239,   239,   241,   242,   243,   245,
     247,   248,   250,   251,   253,   255,   256,   258,   259,   260,
     261,   263,   264,   265,   266,   268,   269,   271,   273,   275,
     276,   277,   278,   280,   281,   282,   284,   285,   286,   287,
     288,   289
};
#endif

#if YYDEBUG || YYERROR_VERBOSE || YYTOKEN_TABLE
/* YYTNAME[SYMBOL-NUM] -- String name of the symbol SYMBOL-NUM.
   First, the terminals, then, starting at YYNTOKENS, nonterminals.  */
static const char *const yytname[] =
{
  "$end", "error", "$undefined", "END", "XMLTOK", "QCLOSE", "QOPEN",
  "DOCTYPE", "COMMENTSTART", "ELEMENT", "ENTITY", "ATTLIST", "NOTATION",
  "CDATASTART", "DOUBLEDASH", "EMPTY", "ANY", "PCDATA", "REQUIRED",
  "IMPLIED", "FIXED", "CDATA", "ID", "IDREF", "IDREFS", "ENTITIES",
  "NMTOKEN", "NMTOKENS", "LETTER", "DIGIT", "SQ", "DQ", "DASH", "DOT",
  "USCORE", "COLON", "CHAR", "OPEN", "CLOSE", "SLASH", "AMPERSAND", "EQ",
  "LSBRACK", "RSBRACK", "PERCENT", "SEMICOLON", "LBRACK", "RBRACK", "STAR",
  "PLUS", "COMMA", "QUERY", "BAR", "SPACE", "$accept", "document",
  "prolog", "misc1", "misc2", "eq", "doctypedecl", "dtd_start1",
  "dtd_start", "int_subset", "int_subset1", "declsep", "pereference",
  "markupdecl", "$@1", "attdef", "$@2", "$@3", "$@4", "$@5", "atttype",
  "enum_type", "enums", "defaultdecl", "opt_star", "contentspec", "mixed",
  "opt_modif", "choice_seq", "choice", "seq", "cp", "pi_start", "xmldecl",
  "pi", "pi_data", "dtd_tok", "chardata1", "chardata", "attchar",
  "attchar1", "attchar2", "name_char", "name_char1", "name", "misc",
  "comment_content", "comment", "$@6", "element", "start1", "start_tag",
  "att_list", "attribute", "att_value", "att_val1", "att_val2",
  "opt_space", "end_tag", "empty_elem", "cdata", "content", "content1", 0
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
     305,   306,   307,   308
};
# endif

/* YYR1[YYN] -- Symbol number of symbol that rule YYN derives.  */
static const yytype_uint8 yyr1[] =
{
       0,    54,    55,    56,    56,    56,    57,    57,    58,    58,
      59,    59,    59,    59,    60,    60,    60,    61,    61,    61,
      62,    63,    63,    64,    64,    65,    65,    66,    68,    67,
      67,    67,    67,    67,    70,    71,    69,    72,    73,    69,
      74,    74,    74,    74,    74,    74,    74,    74,    74,    75,
      75,    76,    76,    77,    77,    77,    77,    78,    78,    79,
      79,    79,    79,    79,    80,    80,    81,    81,    81,    81,
      82,    82,    82,    83,    83,    84,    84,    85,    85,    86,
      87,    88,    88,    89,    89,    89,    89,    90,    90,    91,
      91,    91,    92,    92,    92,    93,    93,    93,    93,    93,
      93,    93,    93,    93,    94,    94,    95,    95,    96,    96,
      96,    97,    97,    97,    97,    98,    98,    99,    99,    99,
     100,   100,   100,   100,   102,   101,   103,   103,   103,   104,
     105,   105,   106,   106,   107,   108,   108,   109,   109,   109,
     109,   110,   110,   110,   110,   111,   111,   112,   113,   114,
     114,   114,   114,   115,   115,   115,   116,   116,   116,   116,
     116,   116
};

/* YYR2[YYN] -- Number of symbols composing right hand side of rule YYN.  */
static const yytype_uint8 yyr2[] =
{
       0,     2,     3,     2,     4,     1,     0,     1,     1,     2,
       1,     2,     3,     2,     2,     5,     6,     2,     6,     8,
       3,     1,     2,     1,     1,     1,     1,     3,     0,     8,
       6,     5,     1,     1,     0,     0,     8,     0,     0,     9,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     6,
       5,     4,     5,     1,     1,     3,     1,     0,     1,     1,
       1,     2,     6,     7,     4,     5,     0,     1,     1,     1,
       6,     6,     5,     4,     5,     4,     5,     2,     2,     1,
       5,     5,     4,     1,     1,     2,     2,     1,     2,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     2,     1,     1,     1,
       1,     1,     2,     2,     0,     5,     1,     2,     3,     2,
       3,     3,     1,     2,     4,     3,     3,     1,     1,     2,
       2,     1,     1,     2,     2,     0,     1,     4,     4,     1,
       1,     2,     2,     1,     1,     2,     1,     2,     1,     2,
       1,     2
};

/* YYDEFACT[STATE-NAME] -- Default rule to reduce with in state
   STATE-NUM when YYTABLE doesn't specify something else to do.  Zero
   means the default is an error.  */
static const yytype_uint8 yydefact[] =
{
       6,    79,   124,   119,     0,     0,     5,     7,     0,     6,
     118,     8,   117,     0,     1,     0,     0,   145,     0,   126,
       0,     9,     0,   111,   112,   113,   114,   115,   145,     3,
      93,    94,   109,   110,    96,    90,   101,   103,    91,   102,
      97,    98,   100,    99,   121,   120,    89,    92,    95,   108,
       0,   129,     2,   146,     0,   132,     0,     0,   150,   160,
     149,   158,   156,   127,   153,     0,   154,     0,   145,   146,
     116,     0,     0,     6,     0,   145,     0,   123,   122,     0,
     130,     0,   133,   131,   145,   161,   159,   157,   152,   151,
     128,   155,     0,    84,     0,    83,    82,     0,     4,    14,
       0,   146,    17,   125,    10,     0,     0,   148,   146,     0,
      80,    81,    86,    85,    20,     0,     0,   145,     0,    26,
       0,    21,    24,    25,    23,    32,    33,    87,     0,    13,
      11,     0,     0,   134,   147,     0,     0,     0,     0,   145,
      22,    88,     0,    12,   105,   138,   104,   137,     0,   107,
     142,   106,   141,     0,     0,   145,    15,    27,     0,   145,
     135,   140,   139,   136,   144,   143,    28,   146,   145,     0,
      16,   146,    18,     0,     0,   146,     0,    31,   145,    59,
      60,   145,   145,    66,    34,     0,    30,    19,     0,     0,
      68,    67,    69,    61,     0,    37,   145,   145,    66,   145,
      66,    29,    44,    40,    41,    42,    43,    45,    46,    47,
     145,     0,    48,     0,   145,     0,     0,    78,   145,   145,
       0,    77,     0,    35,     0,     0,    57,   145,     0,     0,
      72,   145,   145,   145,     0,    38,     0,   145,    58,    62,
       0,    70,   145,    71,   145,     0,     0,   145,     0,    53,
      54,     0,    36,    56,     0,    63,     0,    64,     0,     0,
      75,    73,     0,    50,   145,     0,    39,    65,    74,    76,
      49,   145,     0,    55,     0,    51,    52
};

/* YYDEFGOTO[NTERM-NUM].  */
static const yytype_int16 yydefgoto[] =
{
      -1,     4,     5,     6,     7,   106,    73,    74,    75,   120,
     121,   122,   123,   124,   173,   168,   194,   234,   213,   254,
     211,   212,   247,   252,   239,   182,   214,   193,   198,   218,
     219,   199,    20,     9,    10,    94,   128,    45,    46,    47,
     147,   152,    48,    49,   200,    11,    50,    12,    13,    62,
      17,    18,    54,    55,   253,   148,   153,    56,    63,    19,
      64,    65,    66
};

/* YYPACT[STATE-NUM] -- Index in YYTABLE of the portion describing
   STATE-NUM.  */
#define YYPACT_NINF -159
static const yytype_int16 yypact[] =
{
      25,  -159,  -159,  -159,    43,   -13,  -159,    25,     6,    25,
    -159,  -159,  -159,   383,  -159,   118,    25,   -28,   158,  -159,
     118,  -159,    37,  -159,  -159,  -159,  -159,  -159,   110,    29,
    -159,  -159,  -159,  -159,  -159,  -159,  -159,  -159,  -159,  -159,
    -159,  -159,  -159,  -159,  -159,  -159,  -159,  -159,  -159,  -159,
     357,   347,  -159,   118,    65,  -159,    17,   118,  -159,   435,
    -159,   435,   435,  -159,   461,    36,  -159,   118,   -28,   409,
    -159,    94,    58,    25,   -36,    59,    95,  -159,  -159,    88,
    -159,    96,  -159,  -159,   193,   461,   461,   461,  -159,  -159,
    -159,  -159,   135,  -159,   325,  -159,  -159,   118,  -159,  -159,
       9,   120,  -159,  -159,   101,   116,   -17,  -159,  -159,   134,
    -159,  -159,  -159,  -159,   347,   105,   122,   130,   118,  -159,
      21,  -159,  -159,  -159,  -159,  -159,  -159,  -159,    -5,  -159,
     131,   536,   558,  -159,  -159,   118,   118,   147,   303,   130,
    -159,  -159,   -17,  -159,  -159,  -159,  -159,  -159,   487,  -159,
    -159,  -159,  -159,   513,   269,   580,  -159,  -159,   166,   152,
    -159,  -159,  -159,  -159,  -159,  -159,  -159,   118,   156,   179,
    -159,   -17,  -159,    31,   588,   118,   181,  -159,   130,  -159,
    -159,   130,   130,   167,  -159,   596,  -159,  -159,    79,   191,
    -159,  -159,  -159,  -159,   260,  -159,   130,   130,   167,   130,
     610,  -159,  -159,  -159,  -159,  -159,  -159,  -159,  -159,  -159,
     130,   177,  -159,   260,   130,    24,   178,  -159,   130,   130,
      32,  -159,   118,  -159,   180,    83,   183,   130,   103,    10,
    -159,   130,   130,   193,   151,  -159,   187,   130,  -159,  -159,
     118,  -159,   130,  -159,   130,   178,   178,   130,   126,  -159,
    -159,   186,  -159,  -159,   151,  -159,   118,   347,   178,   178,
    -159,  -159,   127,  -159,   130,   -17,  -159,   347,  -159,  -159,
    -159,   130,   118,  -159,   118,   347,   347
};

/* YYPGOTO[NTERM-NUM].  */
static const yytype_int16 yypgoto[] =
{
    -159,  -159,  -159,     2,  -159,  -159,  -159,  -159,  -159,  -159,
     123,  -159,  -159,  -159,  -159,  -159,  -159,  -159,  -159,  -159,
      34,  -159,  -159,   -12,  -159,  -159,  -159,   -74,    72,  -159,
    -159,  -158,   254,  -159,   -15,  -159,  -159,   -31,     8,   -51,
     107,   112,   -23,    -8,     1,    76,  -159,   -14,  -159,   251,
    -159,  -159,   245,   -46,  -105,  -159,  -159,    81,   203,  -159,
     199,  -159,   204
};

/* YYTABLE[YYPACT[STATE-NUM]].  What to do in state STATE-NUM.  If
   positive, shift that token.  If negative, reduce the rule which
   number is the opposite.  If zero, do what YYDEFACT says.
   If YYTABLE_NINF, syntax error.  */
#define YYTABLE_NINF -1
static const yytype_uint16 yytable[] =
{
      27,   133,    99,    59,    61,    70,   100,    27,    82,    28,
      22,    29,    27,   131,   132,     1,    51,     2,   115,    78,
     116,    28,    82,   141,    15,    53,    60,     1,    70,     2,
     115,     1,   116,     2,    23,    24,    72,   159,    95,    57,
      25,    26,     1,    14,     2,    27,   179,   180,   142,    27,
      59,    61,   117,   118,    79,    83,    70,   243,    84,    27,
     244,    70,   119,   113,   139,   118,   178,    60,    79,    60,
      60,   226,    89,    15,   119,    98,   227,   181,     3,   230,
     146,   151,   231,    21,   232,   125,   126,   260,   261,    27,
      67,    70,    52,    89,    89,    89,   196,   146,   114,    96,
     268,   269,   151,    80,    81,   125,   126,    23,    24,    71,
      27,    97,   101,    25,    26,    70,    23,    24,    67,   138,
      32,    33,    25,    26,   217,   197,   221,    27,    27,   104,
     236,    70,    70,   103,   107,   237,   154,   155,    23,    24,
     110,   105,    32,    33,    25,    26,    23,    24,   127,    92,
     241,    70,    25,    26,   129,   242,   102,   130,   135,    27,
     273,    57,    70,    69,     1,   109,     2,    27,   174,   249,
     250,   251,   134,   263,   270,   136,   185,    70,   264,   271,
      27,   131,   132,   108,   143,   156,    23,    24,    30,    31,
      32,    33,    25,    26,    34,    15,    36,    37,   137,    39,
      40,    41,    42,    43,   170,   171,    23,    24,    27,   175,
      70,    58,    25,    26,    27,   190,   191,   177,   192,   186,
     158,    23,    24,   233,   197,    32,    33,    25,    26,   201,
     223,   238,    27,   235,    70,   255,   169,    27,    27,   265,
     172,   257,   266,   140,    70,   183,   108,   224,    27,   176,
      27,    27,    70,    70,     8,   162,    16,   267,    85,   187,
      86,    87,   188,   189,    27,   165,    27,    68,    90,    91,
     202,     0,     0,   275,     0,   276,     0,   215,   216,     0,
     220,   203,   204,   205,   206,   207,   208,   209,     0,     0,
       0,   222,     0,     0,     0,   225,     0,    23,    24,   228,
     229,    32,    33,    25,    26,     0,   210,     0,   240,     0,
       0,     0,   245,   246,   248,     0,     0,     0,   256,     0,
       0,     0,   166,   258,     0,   259,     0,     0,   262,     0,
     111,    23,    24,     0,     0,    32,    33,    25,    26,     0,
       0,     0,     0,     0,     0,   272,     0,     0,   157,     0,
       0,     0,   274,    23,    24,    30,    31,    32,    33,    25,
      26,    34,    35,    36,    37,    38,    39,    40,    41,    42,
      43,    76,     0,     0,     0,    23,    24,     0,   112,    32,
      33,    25,    26,     0,     0,    23,    24,    30,    31,    32,
      33,    25,    26,    34,    35,    36,    37,    38,    39,    40,
      41,    42,    43,     0,     0,     0,     0,     0,     0,     0,
      77,    23,    24,    30,    31,    32,    33,    25,    26,    34,
      35,    36,    37,    38,    39,    40,    41,    42,    43,     0,
       0,     0,     0,     0,     0,     0,    44,    23,    24,    30,
      31,    32,    33,    25,    26,    34,    35,    36,    37,    38,
      39,    40,    41,    42,    43,     0,     0,     0,     0,     0,
       0,     0,    93,    23,    24,    30,    31,    32,    33,    25,
      26,    34,     0,    36,    37,     0,    39,    40,    41,    42,
      43,     0,     0,     0,     0,     0,     0,     0,    58,    23,
      24,    30,    31,    32,    33,    25,    26,    34,     0,    36,
      37,     0,    39,    40,    41,    42,    43,     0,     0,     0,
       0,     0,     0,     0,    88,    23,    24,   160,   144,    32,
      33,    25,    26,    34,     0,    36,    37,     0,    39,    40,
      41,    42,    43,     0,     0,     0,     0,     0,     0,     0,
     161,    23,    24,   149,   163,    32,    33,    25,    26,    34,
       0,    36,    37,     0,    39,    40,    41,    42,    43,     0,
       0,     0,     0,     0,    23,    24,   164,   144,    32,    33,
      25,    26,    34,     0,    36,    37,     0,    39,    40,    41,
      42,    43,     0,     0,     0,     0,    23,    24,   149,   145,
      32,    33,    25,    26,    34,     0,    36,    37,     0,    39,
      40,    41,    42,    43,     0,     0,     0,     0,    23,    24,
       0,   150,    32,    33,    25,    26,    23,    24,     0,     0,
      32,    33,    25,    26,    23,    24,     0,     0,    32,    33,
      25,    26,     0,   167,     0,     0,     0,     0,    23,    24,
       0,   184,    32,    33,    25,    26,     0,     0,     0,   195,
       0,     0,     0,     0,     0,     0,     0,     0,   190,   191,
       0,   192
};

static const yytype_int16 yycheck[] =
{
       8,   106,    38,    18,    18,    28,    42,    15,    54,     8,
       4,     9,    20,    30,    31,     6,    15,     8,     9,    50,
      11,    20,    68,    28,    37,    53,    18,     6,    51,     8,
       9,     6,    11,     8,    28,    29,     7,   142,    69,     3,
      34,    35,     6,     0,     8,    53,    15,    16,    53,    57,
      65,    65,    43,    44,    53,    38,    79,    47,    57,    67,
      50,    84,    53,    94,    43,    44,   171,    59,    67,    61,
      62,    47,    64,    37,    53,    73,    52,    46,    53,    47,
     131,   132,    50,     7,    52,   100,   100,   245,   246,    97,
      53,   114,    16,    85,    86,    87,    17,   148,    97,     5,
     258,   259,   153,    38,    39,   120,   120,    28,    29,    28,
     118,    53,    53,    34,    35,   138,    28,    29,    53,   118,
      32,    33,    34,    35,   198,    46,   200,   135,   136,    41,
      47,   154,   155,    38,    38,    52,   135,   136,    28,    29,
       5,    53,    32,    33,    34,    35,    28,    29,    28,    68,
      47,   174,    34,    35,    53,    52,    75,    41,    53,   167,
     265,     3,   185,    53,     6,    84,     8,   175,   167,    18,
      19,    20,    38,    47,    47,    53,   175,   200,    52,    52,
     188,    30,    31,    53,    53,    38,    28,    29,    30,    31,
      32,    33,    34,    35,    36,    37,    38,    39,   117,    41,
      42,    43,    44,    45,    38,    53,    28,    29,   216,    53,
     233,    53,    34,    35,   222,    48,    49,    38,    51,    38,
     139,    28,    29,   222,    46,    32,    33,    34,    35,    38,
      53,    48,   240,    53,   257,    48,   155,   245,   246,    53,
     159,   240,   254,   120,   267,   173,    53,   213,   256,   168,
     258,   259,   275,   276,     0,   148,     5,   256,    59,   178,
      61,    62,   181,   182,   272,   153,   274,    22,    65,    65,
      10,    -1,    -1,   272,    -1,   274,    -1,   196,   197,    -1,
     199,    21,    22,    23,    24,    25,    26,    27,    -1,    -1,
      -1,   210,    -1,    -1,    -1,   214,    -1,    28,    29,   218,
     219,    32,    33,    34,    35,    -1,    46,    -1,   227,    -1,
      -1,    -1,   231,   232,   233,    -1,    -1,    -1,   237,    -1,
      -1,    -1,    53,   242,    -1,   244,    -1,    -1,   247,    -1,
       5,    28,    29,    -1,    -1,    32,    33,    34,    35,    -1,
      -1,    -1,    -1,    -1,    -1,   264,    -1,    -1,    45,    -1,
      -1,    -1,   271,    28,    29,    30,    31,    32,    33,    34,
      35,    36,    37,    38,    39,    40,    41,    42,    43,    44,
      45,    14,    -1,    -1,    -1,    28,    29,    -1,    53,    32,
      33,    34,    35,    -1,    -1,    28,    29,    30,    31,    32,
      33,    34,    35,    36,    37,    38,    39,    40,    41,    42,
      43,    44,    45,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      53,    28,    29,    30,    31,    32,    33,    34,    35,    36,
      37,    38,    39,    40,    41,    42,    43,    44,    45,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    53,    28,    29,    30,
      31,    32,    33,    34,    35,    36,    37,    38,    39,    40,
      41,    42,    43,    44,    45,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    53,    28,    29,    30,    31,    32,    33,    34,
      35,    36,    -1,    38,    39,    -1,    41,    42,    43,    44,
      45,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    53,    28,
      29,    30,    31,    32,    33,    34,    35,    36,    -1,    38,
      39,    -1,    41,    42,    43,    44,    45,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    53,    28,    29,    30,    31,    32,
      33,    34,    35,    36,    -1,    38,    39,    -1,    41,    42,
      43,    44,    45,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      53,    28,    29,    30,    31,    32,    33,    34,    35,    36,
      -1,    38,    39,    -1,    41,    42,    43,    44,    45,    -1,
      -1,    -1,    -1,    -1,    28,    29,    53,    31,    32,    33,
      34,    35,    36,    -1,    38,    39,    -1,    41,    42,    43,
      44,    45,    -1,    -1,    -1,    -1,    28,    29,    30,    53,
      32,    33,    34,    35,    36,    -1,    38,    39,    -1,    41,
      42,    43,    44,    45,    -1,    -1,    -1,    -1,    28,    29,
      -1,    53,    32,    33,    34,    35,    28,    29,    -1,    -1,
      32,    33,    34,    35,    28,    29,    -1,    -1,    32,    33,
      34,    35,    -1,    53,    -1,    -1,    -1,    -1,    28,    29,
      -1,    53,    32,    33,    34,    35,    -1,    -1,    -1,    53,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    48,    49,
      -1,    51
};

/* YYSTOS[STATE-NUM] -- The (internal number of the) accessing
   symbol of state STATE-NUM.  */
static const yytype_uint8 yystos[] =
{
       0,     6,     8,    53,    55,    56,    57,    58,    86,    87,
      88,    99,   101,   102,     0,    37,   103,   104,   105,   113,
      86,    99,     4,    28,    29,    34,    35,    97,    98,    57,
      30,    31,    32,    33,    36,    37,    38,    39,    40,    41,
      42,    43,    44,    45,    53,    91,    92,    93,    96,    97,
     100,    98,    99,    53,   106,   107,   111,     3,    53,    88,
      92,   101,   103,   112,   114,   115,   116,    53,   106,    53,
      96,   111,     7,    60,    61,    62,    14,    53,    91,    98,
      38,    39,   107,    38,    98,   114,   114,   114,    53,    92,
     112,   116,   111,    53,    89,    91,     5,    53,    57,    38,
      42,    53,   111,    38,    41,    53,    59,    38,    53,   111,
       5,     5,    53,    91,    98,     9,    11,    43,    44,    53,
      63,    64,    65,    66,    67,    88,   101,    28,    90,    53,
      41,    30,    31,   108,    38,    53,    53,   111,    98,    43,
      64,    28,    53,    53,    31,    53,    93,    94,   109,    30,
      53,    93,    95,   110,    98,    98,    38,    45,   111,   108,
      30,    53,    94,    31,    53,    95,    53,    53,    69,   111,
      38,    53,   111,    68,    98,    53,   111,    38,   108,    15,
      16,    46,    79,    82,    53,    98,    38,   111,   111,   111,
      48,    49,    51,    81,    70,    53,    17,    46,    82,    85,
      98,    38,    10,    21,    22,    23,    24,    25,    26,    27,
      46,    74,    75,    72,    80,   111,   111,    81,    83,    84,
     111,    81,   111,    53,    74,   111,    47,    52,   111,   111,
      47,    50,    52,    98,    71,    53,    47,    52,    48,    78,
     111,    47,    52,    47,    50,   111,   111,    76,   111,    18,
      19,    20,    77,   108,    73,    48,   111,    98,   111,   111,
      85,    85,   111,    47,    52,    53,    77,    98,    85,    85,
      47,    52,   111,   108,   111,    98,    98
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
        case 10:

/* Line 1455 of yacc.c  */
#line 93 "xml_parse.y"
    {}
    break;

  case 11:

/* Line 1455 of yacc.c  */
#line 94 "xml_parse.y"
    {free_string((yyvsp[(1) - (2)].str));}
    break;

  case 12:

/* Line 1455 of yacc.c  */
#line 95 "xml_parse.y"
    {free_string((yyvsp[(1) - (3)].str)); free_string((yyvsp[(3) - (3)].str));}
    break;

  case 13:

/* Line 1455 of yacc.c  */
#line 96 "xml_parse.y"
    {free_string((yyvsp[(2) - (2)].str));}
    break;

  case 14:

/* Line 1455 of yacc.c  */
#line 98 "xml_parse.y"
    {}
    break;

  case 17:

/* Line 1455 of yacc.c  */
#line 102 "xml_parse.y"
    {}
    break;

  case 18:

/* Line 1455 of yacc.c  */
#line 103 "xml_parse.y"
    {if(!check_dtd_tok((yyvsp[(3) - (6)].str),"SYSTEM")) printf("System '%s'\n",get_cstring((yyvsp[(5) - (6)].str)));}
    break;

  case 19:

/* Line 1455 of yacc.c  */
#line 104 "xml_parse.y"
    {if(!check_dtd_tok((yyvsp[(3) - (8)].str),"PUBLIC")) printf("Public '%s' '%s'\n",get_cstring((yyvsp[(5) - (8)].str)),get_cstring((yyvsp[(7) - (8)].str)));}
    break;

  case 20:

/* Line 1455 of yacc.c  */
#line 106 "xml_parse.y"
    {printf("DTD name %s\n",get_cstring((yyvsp[(3) - (3)].str)));}
    break;

  case 21:

/* Line 1455 of yacc.c  */
#line 108 "xml_parse.y"
    {}
    break;

  case 22:

/* Line 1455 of yacc.c  */
#line 109 "xml_parse.y"
    {}
    break;

  case 23:

/* Line 1455 of yacc.c  */
#line 111 "xml_parse.y"
    {}
    break;

  case 24:

/* Line 1455 of yacc.c  */
#line 112 "xml_parse.y"
    {}
    break;

  case 25:

/* Line 1455 of yacc.c  */
#line 114 "xml_parse.y"
    {}
    break;

  case 26:

/* Line 1455 of yacc.c  */
#line 115 "xml_parse.y"
    {}
    break;

  case 27:

/* Line 1455 of yacc.c  */
#line 117 "xml_parse.y"
    {}
    break;

  case 28:

/* Line 1455 of yacc.c  */
#line 119 "xml_parse.y"
    {element_state=1;}
    break;

  case 29:

/* Line 1455 of yacc.c  */
#line 119 "xml_parse.y"
    {element_state=0;}
    break;

  case 32:

/* Line 1455 of yacc.c  */
#line 122 "xml_parse.y"
    {}
    break;

  case 33:

/* Line 1455 of yacc.c  */
#line 123 "xml_parse.y"
    {}
    break;

  case 34:

/* Line 1455 of yacc.c  */
#line 125 "xml_parse.y"
    {attlist_state=1;}
    break;

  case 35:

/* Line 1455 of yacc.c  */
#line 125 "xml_parse.y"
    {attlist_state=2;}
    break;

  case 36:

/* Line 1455 of yacc.c  */
#line 125 "xml_parse.y"
    {attlist_state=0;}
    break;

  case 37:

/* Line 1455 of yacc.c  */
#line 126 "xml_parse.y"
    {attlist_state=1;}
    break;

  case 38:

/* Line 1455 of yacc.c  */
#line 126 "xml_parse.y"
    {attlist_state=2;}
    break;

  case 39:

/* Line 1455 of yacc.c  */
#line 126 "xml_parse.y"
    {attlist_state=0;}
    break;

  case 40:

/* Line 1455 of yacc.c  */
#line 128 "xml_parse.y"
    {}
    break;

  case 41:

/* Line 1455 of yacc.c  */
#line 129 "xml_parse.y"
    {}
    break;

  case 42:

/* Line 1455 of yacc.c  */
#line 130 "xml_parse.y"
    {}
    break;

  case 43:

/* Line 1455 of yacc.c  */
#line 131 "xml_parse.y"
    {}
    break;

  case 44:

/* Line 1455 of yacc.c  */
#line 132 "xml_parse.y"
    {}
    break;

  case 45:

/* Line 1455 of yacc.c  */
#line 133 "xml_parse.y"
    {}
    break;

  case 46:

/* Line 1455 of yacc.c  */
#line 134 "xml_parse.y"
    {}
    break;

  case 47:

/* Line 1455 of yacc.c  */
#line 135 "xml_parse.y"
    {}
    break;

  case 48:

/* Line 1455 of yacc.c  */
#line 136 "xml_parse.y"
    {}
    break;

  case 49:

/* Line 1455 of yacc.c  */
#line 138 "xml_parse.y"
    {}
    break;

  case 50:

/* Line 1455 of yacc.c  */
#line 139 "xml_parse.y"
    {}
    break;

  case 51:

/* Line 1455 of yacc.c  */
#line 141 "xml_parse.y"
    {}
    break;

  case 52:

/* Line 1455 of yacc.c  */
#line 142 "xml_parse.y"
    {}
    break;

  case 53:

/* Line 1455 of yacc.c  */
#line 144 "xml_parse.y"
    {}
    break;

  case 54:

/* Line 1455 of yacc.c  */
#line 145 "xml_parse.y"
    {}
    break;

  case 55:

/* Line 1455 of yacc.c  */
#line 146 "xml_parse.y"
    {}
    break;

  case 56:

/* Line 1455 of yacc.c  */
#line 147 "xml_parse.y"
    {}
    break;

  case 58:

/* Line 1455 of yacc.c  */
#line 150 "xml_parse.y"
    {}
    break;

  case 59:

/* Line 1455 of yacc.c  */
#line 152 "xml_parse.y"
    {}
    break;

  case 60:

/* Line 1455 of yacc.c  */
#line 153 "xml_parse.y"
    {}
    break;

  case 61:

/* Line 1455 of yacc.c  */
#line 154 "xml_parse.y"
    {}
    break;

  case 62:

/* Line 1455 of yacc.c  */
#line 155 "xml_parse.y"
    {}
    break;

  case 64:

/* Line 1455 of yacc.c  */
#line 158 "xml_parse.y"
    {}
    break;

  case 65:

/* Line 1455 of yacc.c  */
#line 159 "xml_parse.y"
    {}
    break;

  case 67:

/* Line 1455 of yacc.c  */
#line 162 "xml_parse.y"
    {}
    break;

  case 68:

/* Line 1455 of yacc.c  */
#line 163 "xml_parse.y"
    {}
    break;

  case 69:

/* Line 1455 of yacc.c  */
#line 164 "xml_parse.y"
    {}
    break;

  case 70:

/* Line 1455 of yacc.c  */
#line 166 "xml_parse.y"
    {}
    break;

  case 71:

/* Line 1455 of yacc.c  */
#line 167 "xml_parse.y"
    {}
    break;

  case 72:

/* Line 1455 of yacc.c  */
#line 168 "xml_parse.y"
    {}
    break;

  case 73:

/* Line 1455 of yacc.c  */
#line 170 "xml_parse.y"
    {}
    break;

  case 74:

/* Line 1455 of yacc.c  */
#line 171 "xml_parse.y"
    {}
    break;

  case 75:

/* Line 1455 of yacc.c  */
#line 173 "xml_parse.y"
    {}
    break;

  case 76:

/* Line 1455 of yacc.c  */
#line 174 "xml_parse.y"
    {}
    break;

  case 77:

/* Line 1455 of yacc.c  */
#line 176 "xml_parse.y"
    {}
    break;

  case 78:

/* Line 1455 of yacc.c  */
#line 177 "xml_parse.y"
    {}
    break;

  case 79:

/* Line 1455 of yacc.c  */
#line 179 "xml_parse.y"
    {pi_state=1;}
    break;

  case 80:

/* Line 1455 of yacc.c  */
#line 181 "xml_parse.y"
    {pi_state=0; check_declaration((yyvsp[(3) - (5)].arg));}
    break;

  case 81:

/* Line 1455 of yacc.c  */
#line 183 "xml_parse.y"
    {pi_state=0; free_string((yyvsp[(3) - (5)].str)); if(call->pi) call->pi((yyvsp[(2) - (5)].str),(yyvsp[(4) - (5)].str)); free_string((yyvsp[(2) - (5)].str)); free_string((yyvsp[(4) - (5)].str)); }
    break;

  case 82:

/* Line 1455 of yacc.c  */
#line 184 "xml_parse.y"
    {pi_state=0; if(call->pi) call->pi((yyvsp[(2) - (4)].str),0); free_string((yyvsp[(2) - (4)].str)); }
    break;

  case 83:

/* Line 1455 of yacc.c  */
#line 186 "xml_parse.y"
    {(yyval.str)=add_to_string(0,(yyvsp[(1) - (1)].c));}
    break;

  case 85:

/* Line 1455 of yacc.c  */
#line 188 "xml_parse.y"
    {(yyval.str)=add_to_string((yyvsp[(1) - (2)].str),(yyvsp[(2) - (2)].c));}
    break;

  case 86:

/* Line 1455 of yacc.c  */
#line 189 "xml_parse.y"
    {(yyval.str)=add_strings((yyvsp[(1) - (2)].str),(yyvsp[(2) - (2)].str));}
    break;

  case 87:

/* Line 1455 of yacc.c  */
#line 191 "xml_parse.y"
    {(yyval.str)=add_to_string(0,(yyvsp[(1) - (1)].c));}
    break;

  case 88:

/* Line 1455 of yacc.c  */
#line 192 "xml_parse.y"
    {(yyval.str)=add_to_string((yyvsp[(1) - (2)].str),(yyvsp[(2) - (2)].c));}
    break;

  case 115:

/* Line 1455 of yacc.c  */
#line 227 "xml_parse.y"
    {(yyval.str)=add_to_string(0,(yyvsp[(1) - (1)].c));}
    break;

  case 116:

/* Line 1455 of yacc.c  */
#line 228 "xml_parse.y"
    {(yyval.str)=add_to_string((yyvsp[(1) - (2)].str),(yyvsp[(2) - (2)].c));}
    break;

  case 117:

/* Line 1455 of yacc.c  */
#line 230 "xml_parse.y"
    {}
    break;

  case 118:

/* Line 1455 of yacc.c  */
#line 231 "xml_parse.y"
    {}
    break;

  case 119:

/* Line 1455 of yacc.c  */
#line 232 "xml_parse.y"
    {free_string((yyvsp[(1) - (1)].str));}
    break;

  case 120:

/* Line 1455 of yacc.c  */
#line 234 "xml_parse.y"
    {(yyval.str)=add_to_string(0,(yyvsp[(1) - (1)].c));}
    break;

  case 122:

/* Line 1455 of yacc.c  */
#line 236 "xml_parse.y"
    {(yyval.str)=add_to_string((yyvsp[(1) - (2)].str),(yyvsp[(2) - (2)].c));}
    break;

  case 123:

/* Line 1455 of yacc.c  */
#line 237 "xml_parse.y"
    {(yyval.str)=add_strings((yyvsp[(1) - (2)].str),(yyvsp[(2) - (2)].str));}
    break;

  case 124:

/* Line 1455 of yacc.c  */
#line 239 "xml_parse.y"
    {in_comment=1;}
    break;

  case 125:

/* Line 1455 of yacc.c  */
#line 239 "xml_parse.y"
    {in_comment=0; if(call->comment) call->comment((yyvsp[(3) - (5)].str)); free_string((yyvsp[(3) - (5)].str));}
    break;

  case 129:

/* Line 1455 of yacc.c  */
#line 245 "xml_parse.y"
    {(yyval.str)=(yyvsp[(2) - (2)].str);}
    break;

  case 130:

/* Line 1455 of yacc.c  */
#line 247 "xml_parse.y"
    {check_start((yyvsp[(1) - (3)].str),(yyvsp[(2) - (3)].arg)); }
    break;

  case 131:

/* Line 1455 of yacc.c  */
#line 248 "xml_parse.y"
    {check_start((yyvsp[(1) - (3)].str),0); }
    break;

  case 133:

/* Line 1455 of yacc.c  */
#line 251 "xml_parse.y"
    {(yyvsp[(2) - (2)].arg)->next=(yyvsp[(1) - (2)].arg); (yyval.arg)=(yyvsp[(2) - (2)].arg);}
    break;

  case 134:

/* Line 1455 of yacc.c  */
#line 253 "xml_parse.y"
    {free_string((yyvsp[(1) - (4)].str)); (yyval.arg)=make_arg((yyvsp[(2) - (4)].str),(yyvsp[(4) - (4)].str));}
    break;

  case 135:

/* Line 1455 of yacc.c  */
#line 255 "xml_parse.y"
    {(yyval.str)=(yyvsp[(2) - (3)].str);}
    break;

  case 136:

/* Line 1455 of yacc.c  */
#line 256 "xml_parse.y"
    {(yyval.str)=(yyvsp[(2) - (3)].str);}
    break;

  case 137:

/* Line 1455 of yacc.c  */
#line 258 "xml_parse.y"
    {(yyval.str)=add_to_string(0,(yyvsp[(1) - (1)].c));}
    break;

  case 139:

/* Line 1455 of yacc.c  */
#line 260 "xml_parse.y"
    {(yyval.str)=add_to_string((yyvsp[(1) - (2)].str),(yyvsp[(2) - (2)].c));}
    break;

  case 140:

/* Line 1455 of yacc.c  */
#line 261 "xml_parse.y"
    {(yyval.str)=add_strings((yyvsp[(1) - (2)].str),(yyvsp[(2) - (2)].str));}
    break;

  case 141:

/* Line 1455 of yacc.c  */
#line 263 "xml_parse.y"
    {(yyval.str)=add_to_string(0,(yyvsp[(1) - (1)].c));}
    break;

  case 143:

/* Line 1455 of yacc.c  */
#line 265 "xml_parse.y"
    {(yyval.str)=add_to_string((yyvsp[(1) - (2)].str),(yyvsp[(2) - (2)].c));}
    break;

  case 144:

/* Line 1455 of yacc.c  */
#line 266 "xml_parse.y"
    {(yyval.str)=add_strings((yyvsp[(1) - (2)].str),(yyvsp[(2) - (2)].str));}
    break;

  case 146:

/* Line 1455 of yacc.c  */
#line 269 "xml_parse.y"
    {free_string((yyvsp[(1) - (1)].str));}
    break;

  case 147:

/* Line 1455 of yacc.c  */
#line 271 "xml_parse.y"
    {check_end((yyvsp[(2) - (4)].str),0);}
    break;

  case 148:

/* Line 1455 of yacc.c  */
#line 273 "xml_parse.y"
    {check_start((yyvsp[(1) - (4)].str),(yyvsp[(2) - (4)].arg)); check_end((yyvsp[(1) - (4)].str),1);}
    break;

  case 149:

/* Line 1455 of yacc.c  */
#line 275 "xml_parse.y"
    {(yyval.str)=add_to_string(0,(yyvsp[(1) - (1)].c));}
    break;

  case 151:

/* Line 1455 of yacc.c  */
#line 277 "xml_parse.y"
    {(yyval.str)=add_to_string((yyvsp[(1) - (2)].str),(yyvsp[(2) - (2)].c));}
    break;

  case 152:

/* Line 1455 of yacc.c  */
#line 278 "xml_parse.y"
    {(yyval.str)=add_strings((yyvsp[(1) - (2)].str),(yyvsp[(2) - (2)].str));}
    break;

  case 153:

/* Line 1455 of yacc.c  */
#line 280 "xml_parse.y"
    { if(check_ws((yyvsp[(1) - (1)].str)) && call->content) call->content((yyvsp[(1) - (1)].str)); free_string((yyvsp[(1) - (1)].str)); }
    break;

  case 154:

/* Line 1455 of yacc.c  */
#line 281 "xml_parse.y"
    { if(check_ws((yyvsp[(1) - (1)].str)) && call->content) call->content((yyvsp[(1) - (1)].str)); free_string((yyvsp[(1) - (1)].str)); }
    break;

  case 155:

/* Line 1455 of yacc.c  */
#line 282 "xml_parse.y"
    { if(check_ws((yyvsp[(2) - (2)].str)) && call->content) call->content((yyvsp[(2) - (2)].str)); free_string((yyvsp[(2) - (2)].str)); }
    break;

  case 156:

/* Line 1455 of yacc.c  */
#line 284 "xml_parse.y"
    {(yyval.str)=0;}
    break;

  case 157:

/* Line 1455 of yacc.c  */
#line 285 "xml_parse.y"
    {(yyval.str)=(yyvsp[(2) - (2)].str);}
    break;

  case 158:

/* Line 1455 of yacc.c  */
#line 286 "xml_parse.y"
    {(yyval.str)=0;}
    break;

  case 159:

/* Line 1455 of yacc.c  */
#line 287 "xml_parse.y"
    {(yyval.str)=(yyvsp[(2) - (2)].str);}
    break;

  case 160:

/* Line 1455 of yacc.c  */
#line 288 "xml_parse.y"
    {(yyval.str)=0;}
    break;

  case 161:

/* Line 1455 of yacc.c  */
#line 289 "xml_parse.y"
    {(yyval.str)=(yyvsp[(2) - (2)].str);}
    break;



/* Line 1455 of yacc.c  */
#line 2609 "y.tab.c"
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
#line 291 "xml_parse.y"


static int check_ws(string *s)
{
  int i=0;
  char *p;

  if(s) {
    p=get_cstring(s);
    while(*p) {
      if(!isspace((int)*p)) break;
      p++;
    }
    if(*p) i=1;
  }
  return i;
}

static int check_dtd_tok(string *s,char *p) 
{
  int i;

  if(s && !strcmp(get_cstring(s),p)) i=0;
  else {
    i=1;
    if(call->error) call->error(*line_no+1,"Illegal DTD token");
    abt_flag=1;
  }
  return i;
}

static void check_declaration(args *attr)
{
  int st=0,err=0,i;
  char *atts[]={"version","encoding","standalone"};
  args *attr1;
	
  attr1=attr;
  while(attr && !err) {
    for(i=0;i<3;i++) if(!strcmp(attr->name,atts[i])) break;
    if(i==3) {
      if(call->warning) call->warning(*line_no+1,"Illegal XML declaration attribute '%s'\n",attr->name);
      err=1;
      break;
    }
    switch(i) {
    case 0:
      if(st<3) st=3;
      else err=16;
      break;
    case 1:
      if(st<2) st=2;
      else err=16;
      break;
    case 2:
      if(st<1) {
	if(strcmp(attr->att,"yes") && strcmp(attr->att,"no")) {
	  err=17;
	  if(call->error) call->error(*line_no+1,"Illegal XML standalone declaration value '%s'\n",attr->att);
	} else st=1;
      } else err=16;
      break;
    }
    attr=attr->next;
  }
  if(err&16) {
    if(!(err&15) && call->error) call->error(*line_no+1,"Illegal XML declaration attribute order\n");
    abt_flag=1;
  } else if(call->declaration) call->declaration(attr1);
  free_att_list(attr1);
}

static void check_start(string *s,args *attr)
{
  state *st;
	
  if((st=free_state_list)) free_state_list=st->next;
  else st=malloc(sizeof(state));
  st->next=0;
  st->prev=curr_state;
  st->name=s;
  if(curr_state) curr_state->next=st;
  curr_state=st;
  if(call->start_element) call->start_element(s,attr);
  free_att_list(attr);
}

static void check_end(string *s,int fg)
{
  state *st;

  if(!(curr_state)) {
    if(call->error) call->error(*line_no+1,"spurious close of '%s' element\n",get_cstring(s));
    abt_flag=1;
  } else if(strcmp(get_cstring(curr_state->name),get_cstring(s))) {
    if(call->error) call->error(*line_no+1,"mismatched close of '%s' element by '%s'\n",get_cstring(curr_state->name),get_cstring(s));
    abt_flag=1;
  } else {
    st=curr_state->prev;
    if(!fg) free_string(curr_state->name); 
    curr_state->name=0;
    curr_state->next=free_state_list;
    free_state_list=curr_state;
    curr_state=st;
    if(call->end_element) call->end_element(s);
  }
  free_string(s);
}

static args *make_arg(string *s1,string *s2)
{
  args *arg;
	
  arg=malloc(sizeof(args));
  arg->next=0;
  arg->name=extract_cstring(s1);
  arg->att=extract_cstring(s2);
  return arg;
}

static void free_att_list(args *arg)
{
  struct args *arg1;
	
  while(arg) {
    arg1=arg->next;
    free(arg->name);
    free(arg->att);
    free(arg);
    arg=arg1;
  }
}

static int xmlerror(char *s)
{
  if(!(abt_flag)) fprintf(stderr,"line %d, col %d: %s\n",*line_no+1,col_pos,s);
  return 0;
}

static int lex_tab[]={
  CHAR,CHAR,CHAR,CHAR,CHAR,CHAR,CHAR,CHAR,SPACE,SPACE,-1,CHAR,CHAR,-2,CHAR,CHAR,
  CHAR,CHAR,CHAR,CHAR,CHAR,CHAR,CHAR,CHAR,CHAR,CHAR,CHAR,CHAR,CHAR,CHAR,CHAR,CHAR,
  SPACE,CHAR,DQ,CHAR,CHAR,PERCENT,AMPERSAND,SQ,CHAR,CHAR,CHAR,CHAR,CHAR,DASH,DOT,SLASH,
  DIGIT,DIGIT,DIGIT,DIGIT,DIGIT,DIGIT,DIGIT,DIGIT,DIGIT,DIGIT,COLON,SEMICOLON,OPEN,EQ,CLOSE,CHAR,
  CHAR,LETTER,LETTER,LETTER,LETTER,LETTER,LETTER,LETTER,LETTER,LETTER,LETTER,LETTER,LETTER,LETTER,LETTER,LETTER,
  LETTER,LETTER,LETTER,LETTER,LETTER,LETTER,LETTER,LETTER,LETTER,LETTER,LETTER,LSBRACK,CHAR,RSBRACK,CHAR,USCORE,
  CHAR,LETTER,LETTER,LETTER,LETTER,LETTER,LETTER,LETTER,LETTER,LETTER,LETTER,LETTER,LETTER,LETTER,LETTER,LETTER,
  LETTER,LETTER,LETTER,LETTER,LETTER,LETTER,LETTER,LETTER,LETTER,LETTER,LETTER,CHAR,CHAR,CHAR,CHAR,CHAR
};

static int lex_tab1[]={
  CHAR,CHAR,CHAR,CHAR,CHAR,CHAR,CHAR,CHAR,SPACE,SPACE,-1,CHAR,CHAR,-2,CHAR,CHAR,
  CHAR,CHAR,CHAR,CHAR,CHAR,CHAR,CHAR,CHAR,CHAR,CHAR,CHAR,CHAR,CHAR,CHAR,CHAR,CHAR,
  SPACE,CHAR,DQ,CHAR,CHAR,PERCENT,AMPERSAND,SQ,LBRACK,RBRACK,STAR,PLUS,COMMA,DASH,DOT,SLASH,
  DIGIT,DIGIT,DIGIT,DIGIT,DIGIT,DIGIT,DIGIT,DIGIT,DIGIT,DIGIT,COLON,SEMICOLON,OPEN,EQ,CLOSE,QUERY,
  CHAR,LETTER,LETTER,LETTER,LETTER,LETTER,LETTER,LETTER,LETTER,LETTER,LETTER,LETTER,LETTER,LETTER,LETTER,LETTER,
  LETTER,LETTER,LETTER,LETTER,LETTER,LETTER,LETTER,LETTER,LETTER,LETTER,LETTER,LSBRACK,CHAR,RSBRACK,CHAR,USCORE,
  CHAR,LETTER,LETTER,LETTER,LETTER,LETTER,LETTER,LETTER,LETTER,LETTER,LETTER,LETTER,LETTER,LETTER,LETTER,LETTER,
  LETTER,LETTER,LETTER,LETTER,LETTER,LETTER,LETTER,LETTER,LETTER,LETTER,LETTER,CHAR,BAR,CHAR,CHAR,CHAR
};
	  
static int xmllex(yystype *lval)
{
  static int c,fg;
  char *str[]={"OCTYPE","EMENT","TITY","TTLIST","OTATION","CDATA[","MPTY","NY","PCDATA",
	       "EQUIRED","MPLIED","IXED"};
  char *str1[]={"CDATA","IDREFS","IDREF","ID","ENTITY","ENTITIES","NMTOKENS","NMTOKEN",0};
  int atttype[]={CDATA,IDREFS,IDREF,ID,ENTITY,ENTITIES,NMTOKENS,NMTOKEN};
  char *p,buf[16];

  string *s;
  int i=-1,j,c1,err=0,buflen=16;
  
  if(abt_flag) return 0;
  if(!fg) c=fgetc(xmlin);
  else fg=0;
  col_pos++;
  if(c==EOF) return 0;
  if(pi_state) {
    i=0;
    if(pi_state==1 && c=='x') {
      c1=fgetc(xmlin);
      if(c1=='m') {
	c1=fgetc(xmlin);
	if(c1=='l') {
	  i=XMLTOK;
	  col_pos+=2;
	  fg=0;
	} else ungetc(c1,xmlin);
      }
      if(!i) {
	i=lex_tab[c];
	lval->c=(char)c;
	c=c1;
	fg=1;
      }
    } else if(c=='?') {
      c1=fgetc(xmlin);
      if(c1=='>') {
	i=QCLOSE;
	col_pos++;
      } else {
	i=lex_tab[c];
	lval->c=(char)c;
	c=c1;
	fg=1;
      }
    } else {
      i=lex_tab[c];
      lval->c=(char)c;
    }
    pi_state=2;
  } else if(in_comment && c=='-') {
    c1=fgetc(xmlin);
    if(c1=='-') {
      i=DOUBLEDASH;
      col_pos++;
    } else {
      i=lex_tab[c];
      lval->c=(char)c;
      c=c1;
      fg=1;
    }
  } else if(element_state) {
    j=0;
    if(element_state==1) {
      if(c=='E') {
	p=str[6];
	while(*p) {
	  c1=fgetc(xmlin);
	  if(c1!=*p) break;
	  p++;
	}
	if(*p) err=1;
	else i=EMPTY;
	j=1;
      } else if (c=='A') {
	p=str[7];
	while(*p) {
	  c1=fgetc(xmlin);
	  if(c1!=*p) break;
	  p++;
	}
	if(*p) err=1;
	else i=ANY;
	j=1;
      }
    } else if(element_state==2 && c=='#') {
      p=str[8];
      while(*p) {
	c1=fgetc(xmlin);
	if(c1!=*p) break;
	p++;
      }
      if(*p) err=1;
      else i=PCDATA;
      j=1;
    }
    element_state=2;
    if(!j) {
      i=lex_tab1[c];
      lval->c=(char)c;
    }
  } else if(attlist_state==1) {
    i=lex_tab1[c];
    lval->c=(char)c;
    if(i==LETTER) {
      buf[0]=c;
      for(j=1;j<buflen-1;j++) {
	c1=fgetc(xmlin);
	if(c1==EOF || isspace(c1)) {
	  c=c1;
	  fg=1;
	  break;
	}
	buf[j]=c1;
      }
      buf[j]=0;
      i=0;
      p=str1[i];
      while((p=str1[i])) {
	if(!strcmp(p,buf)) break;
	i++;
      }
      if(p) i=atttype[i];
      else err=1;
    }
    attlist_state=2;
  } else if(attlist_state==2) {
    if(c=='#') {
      c1=fgetc(xmlin);
      switch(c1) {
      case 'R':
	p=str[9];
	i=REQUIRED;
	break;
      case 'I':
	p=str[10];
	i=IMPLIED;
	break;
      case 'F':
	p=str[11];
	i=FIXED;
	break;
      default:
	err=1;
      }
      if(!err) {
	while(*p) {
	  c1=fgetc(xmlin);
	  if(c1!=*p) break;
	  p++;
	}
	if(*p) err=1;
      }
    } else {
      i=lex_tab1[c];
      lval->c=(char)c;
    }
  } else if(!in_comment && c=='<') {
    c1=fgetc(xmlin);
    col_pos++;
    err=0;
    if(c1=='?') i=QOPEN;
    else if(c1=='!') {
      col_pos++;
      c1=fgetc(xmlin);
      switch(c1) {
      case '-':
	c1=fgetc(xmlin);
	if(c1=='-') i=COMMENTSTART;
	else err=1;
	break;
      case '[':
	p=str[5];
	while(*p) {
	  c1=fgetc(xmlin);
	  if(c1!=*p) break;
	  p++;
	}
	if(*p) err=1;
	else i=CDATASTART;
	break;
      case 'D':
	p=str[0];
	while(*p) {
	  c1=fgetc(xmlin);
	  if(c1!=*p) break;
	  p++;
	}
	if(*p) err=1;
	else i=DOCTYPE;
	break;
      case 'E':
	c1=fgetc(xmlin);
	if(c1=='L') j=1;
	else if(c1=='N') j=2;
	else err=1;
	if(!err) {
	  p=str[j];
	  while(*p) {
	    c1=fgetc(xmlin);
	    if(c1!=*p) break;
	    p++;
	  }
	  if(*p) err=1;
	  else if(j==1) {
	    i=ELEMENT;
	  } else i=ENTITY;
	}
	break;
      case 'A':
 	p=str[3];
	while(*p) {
	  c1=fgetc(xmlin);
	  if(c1!=*p) break;
	  p++;
	}
	if(*p) err=1;
	else i=ATTLIST;
	break;
      case 'N':
 	p=str[4];
	while(*p) {
	  c1=fgetc(xmlin);
	  if(c1!=*p) break;
	  p++;
	}
	if(*p) err=1;
	else i=NOTATION;
	break;
      default:
	err=1;
	break;
      }
    } else if(c1=='/') i=END;
    else {
      col_pos--;
      i=lex_tab[c];
      lval->c=(char)c;
      c=c1;
      fg=1;
    }
  } else {
    i=lex_tab[c];
    lval->c=(char)c;
  }
  if(err) {
    if(call->error) call->error(*line_no+1,"Illegal construction\n");
    abt_flag=1;
  }
  if(i<0) {
    (*line_no)++;
    col_pos=0;
    if(i==-2) {
      c=fgetc(xmlin);
      if(c!='\n') fg=1;
      else col_pos=1;
      lval->c='\n';
    }
    i=SPACE;
  }
  if(i==SPACE) {
    s=0;
    c1=lval->c;
    do {
      s=add_to_string(s,c1);
      if(!fg) c=fgetc(xmlin);
      else fg=0;
      col_pos++;
      if(c==EOF) break;
      c1=c;
      i=lex_tab[c];
      if(i<0) {
	(*line_no)++;
	col_pos=0;
	if(i==-2) {
	  c=fgetc(xmlin);
	  if(c!='\n') fg=1;
	  else col_pos++;
	  c1='\n';
	}
	i=SPACE;
      }
    } while(i==SPACE);
    fg=1;
    col_pos--;
    lval->str=s;
    return SPACE;
  }
  return i;
}

int Read_XML(FILE *fptr,XML_handler *handler)
{
  int err,line_bk=0;
  state *st;
	
#if YYDEBUG
  xmldebug=1;
#endif
  curr_state=free_state_list=0;
  call=handler;
  xmlin=fptr;
  line_no=handler->line;
  if(!line_no) line_no=&line_bk;
  err=xmlparse();
  if(!err && abt_flag) err=-1;
  if(!err && curr_state) {
    if(call->error) call->error(0,"The following elements were not closed:");
    while(curr_state) {
      st=curr_state->prev;
      if(call->error) call->error(0," %s",get_cstring(curr_state->name));
      free_string(curr_state->name);
      free(curr_state);
      curr_state=st;
    }
    err=-2;
  }
  if(free_state_list) {
    while(free_state_list) {
      st=free_state_list->next;
      free(free_state_list);
      free_state_list=st;
    }
  }
  return err;
}

