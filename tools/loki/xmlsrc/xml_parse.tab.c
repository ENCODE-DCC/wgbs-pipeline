/* A Bison parser, made by GNU Bison 2.0.  */

/* Skeleton parser for Yacc-like parsing with Bison,
   Copyright (C) 1984, 1989, 1990, 2000, 2001, 2002, 2003, 2004 Free Software Foundation, Inc.

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2, or (at your option)
   any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place - Suite 330,
   Boston, MA 02111-1307, USA.  */

/* As a special exception, when this file is copied by Bison into a
   Bison output file, you may use that output file without restriction.
   This special exception was added by the Free Software Foundation
   in version 1.24 of Bison.  */

/* Written by Richard Stallman by simplifying the original so called
   ``semantic'' parser.  */

/* All symbols defined below should begin with yy or YY, to avoid
   infringing on user name space.  This should be done even for local
   variables, as they might otherwise be expanded by user macros.
   There are some unavoidable exceptions within include files to
   define necessary library symbols; they are noted "INFRINGES ON
   USER NAME SPACE" below.  */

/* Identify Bison output.  */
#define YYBISON 1

/* Skeleton name.  */
#define YYSKELETON_NAME "yacc.c"

/* Pure parsers.  */
#define YYPURE 1

/* Using locations.  */
#define YYLSP_NEEDED 0



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




/* Copy the first part of user declarations.  */
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

#if ! defined (YYSTYPE) && ! defined (YYSTYPE_IS_DECLARED)
#line 62 "xml_parse.y"
typedef union YYSTYPE {
  string *str;
  args *arg;
  char c;
  int i;
} YYSTYPE;
/* Line 190 of yacc.c.  */
#line 249 "xml_parse.tab.c"
# define yystype YYSTYPE /* obsolescent; will be withdrawn */
# define YYSTYPE_IS_DECLARED 1
# define YYSTYPE_IS_TRIVIAL 1
#endif



/* Copy the second part of user declarations.  */


/* Line 213 of yacc.c.  */
#line 261 "xml_parse.tab.c"

#if ! defined (yyoverflow) || YYERROR_VERBOSE

# ifndef YYFREE
#  define YYFREE free
# endif
# ifndef YYMALLOC
#  define YYMALLOC malloc
# endif

/* The parser invokes alloca or malloc; define the necessary symbols.  */

# ifdef YYSTACK_USE_ALLOCA
#  if YYSTACK_USE_ALLOCA
#   ifdef __GNUC__
#    define YYSTACK_ALLOC __builtin_alloca
#   else
#    define YYSTACK_ALLOC alloca
#   endif
#  endif
# endif

# ifdef YYSTACK_ALLOC
   /* Pacify GCC's `empty if-body' warning. */
#  define YYSTACK_FREE(Ptr) do { /* empty */; } while (0)
# else
#  if defined (__STDC__) || defined (__cplusplus)
#   include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
#   define YYSIZE_T size_t
#  endif
#  define YYSTACK_ALLOC YYMALLOC
#  define YYSTACK_FREE YYFREE
# endif
#endif /* ! defined (yyoverflow) || YYERROR_VERBOSE */


#if (! defined (yyoverflow) \
     && (! defined (__cplusplus) \
	 || (defined (YYSTYPE_IS_TRIVIAL) && YYSTYPE_IS_TRIVIAL)))

/* A type that is properly aligned for any stack member.  */
union yyalloc
{
  short int yyss;
  YYSTYPE yyvs;
  };

/* The size of the maximum gap between one aligned stack and the next.  */
# define YYSTACK_GAP_MAXIMUM (sizeof (union yyalloc) - 1)

/* The size of an array large to enough to hold all stacks, each with
   N elements.  */
# define YYSTACK_BYTES(N) \
     ((N) * (sizeof (short int) + sizeof (YYSTYPE))			\
      + YYSTACK_GAP_MAXIMUM)

/* Copy COUNT objects from FROM to TO.  The source and destination do
   not overlap.  */
# ifndef YYCOPY
#  if defined (__GNUC__) && 1 < __GNUC__
#   define YYCOPY(To, From, Count) \
      __builtin_memcpy (To, From, (Count) * sizeof (*(From)))
#  else
#   define YYCOPY(To, From, Count)		\
      do					\
	{					\
	  register YYSIZE_T yyi;		\
	  for (yyi = 0; yyi < (Count); yyi++)	\
	    (To)[yyi] = (From)[yyi];		\
	}					\
      while (0)
#  endif
# endif

/* Relocate STACK from its old location to the new one.  The
   local variables YYSIZE and YYSTACKSIZE give the old and new number of
   elements in the stack, and YYPTR gives the new location of the
   stack.  Advance YYPTR to a properly aligned location for the next
   stack.  */
# define YYSTACK_RELOCATE(Stack)					\
    do									\
      {									\
	YYSIZE_T yynewbytes;						\
	YYCOPY (&yyptr->Stack, Stack, yysize);				\
	Stack = &yyptr->Stack;						\
	yynewbytes = yystacksize * sizeof (*Stack) + YYSTACK_GAP_MAXIMUM; \
	yyptr += yynewbytes / sizeof (*yyptr);				\
      }									\
    while (0)

#endif

#if defined (__STDC__) || defined (__cplusplus)
   typedef signed char yysigned_char;
#else
   typedef short int yysigned_char;
#endif

/* YYFINAL -- State number of the termination state. */
#define YYFINAL  14
/* YYLAST -- Last index in YYTABLE.  */
#define YYLAST   670

/* YYNTOKENS -- Number of terminals. */
#define YYNTOKENS  54
/* YYNNTS -- Number of nonterminals. */
#define YYNNTS  63
/* YYNRULES -- Number of rules. */
#define YYNRULES  161
/* YYNRULES -- Number of states. */
#define YYNSTATES  277

/* YYTRANSLATE(YYLEX) -- Bison symbol number corresponding to YYLEX.  */
#define YYUNDEFTOK  2
#define YYMAXUTOK   308

#define YYTRANSLATE(YYX) 						\
  ((unsigned int) (YYX) <= YYMAXUTOK ? yytranslate[YYX] : YYUNDEFTOK)

/* YYTRANSLATE[YYLEX] -- Bison symbol number corresponding to YYLEX.  */
static const unsigned char yytranslate[] =
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
static const unsigned short int yyprhs[] =
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

/* YYRHS -- A `-1'-separated list of the rules' RHS. */
static const yysigned_char yyrhs[] =
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
      -1,    46,   111,    90,    76,   111,    47,    -1,    46,   111,
      90,   111,    47,    -1,   111,    52,   111,    90,    -1,    76,
     111,    52,   111,    90,    -1,    18,    -1,    19,    -1,    20,
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
static const unsigned short int yyrline[] =
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

#if YYDEBUG || YYERROR_VERBOSE
/* YYTNME[SYMBOL-NUM] -- String name of the symbol SYMBOL-NUM.
   First, the terminals, then, starting at YYNTOKENS, nonterminals. */
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
  "markupdecl", "@1", "attdef", "@2", "@3", "@4", "@5", "atttype",
  "enum_type", "enums", "defaultdecl", "opt_star", "contentspec", "mixed",
  "opt_modif", "choice_seq", "choice", "seq", "cp", "pi_start", "xmldecl",
  "pi", "pi_data", "dtd_tok", "chardata1", "chardata", "attchar",
  "attchar1", "attchar2", "name_char", "name_char1", "name", "misc",
  "comment_content", "comment", "@6", "element", "start1", "start_tag",
  "att_list", "attribute", "att_value", "att_val1", "att_val2",
  "opt_space", "end_tag", "empty_elem", "cdata", "content", "content1", 0
};
#endif

# ifdef YYPRINT
/* YYTOKNUM[YYLEX-NUM] -- Internal token number corresponding to
   token YYLEX-NUM.  */
static const unsigned short int yytoknum[] =
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
static const unsigned char yyr1[] =
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
static const unsigned char yyr2[] =
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
static const unsigned char yydefact[] =
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

/* YYDEFGOTO[NTERM-NUM]. */
static const short int yydefgoto[] =
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
#define YYPACT_NINF -224
static const short int yypact[] =
{
       0,  -224,  -224,  -224,    18,    27,  -224,     0,    48,     0,
    -224,  -224,  -224,   384,  -224,   121,     0,    20,   291,  -224,
     121,  -224,    28,  -224,  -224,  -224,  -224,  -224,   183,    32,
    -224,  -224,  -224,  -224,  -224,  -224,  -224,  -224,  -224,  -224,
    -224,  -224,  -224,  -224,  -224,  -224,  -224,  -224,  -224,  -224,
     358,   313,  -224,   121,    69,  -224,    33,   121,  -224,   436,
    -224,   436,   436,  -224,   462,    11,  -224,   121,    20,   410,
    -224,   104,    73,     0,    42,    87,   101,  -224,  -224,   113,
    -224,   115,  -224,  -224,   234,   462,   462,   462,  -224,  -224,
    -224,  -224,   138,  -224,   326,  -224,  -224,   121,  -224,  -224,
      35,   133,  -224,  -224,   111,   124,    93,  -224,  -224,   143,
    -224,  -224,  -224,  -224,   313,   129,   130,   139,   121,  -224,
      91,  -224,  -224,  -224,  -224,  -224,  -224,  -224,   -26,  -224,
     140,   537,   559,  -224,  -224,   121,   121,   153,   273,   139,
    -224,  -224,    93,  -224,  -224,  -224,  -224,  -224,   488,  -224,
    -224,  -224,  -224,   514,   581,   589,  -224,  -224,   159,   151,
    -224,  -224,  -224,  -224,  -224,  -224,  -224,   121,   156,   172,
    -224,    93,  -224,     9,   597,   121,   175,  -224,   139,  -224,
    -224,   139,   139,   -18,  -224,   611,  -224,  -224,   161,   181,
    -224,  -224,  -224,  -224,   148,  -224,   139,   139,   -18,   139,
     619,  -224,  -224,  -224,  -224,  -224,  -224,  -224,  -224,  -224,
     139,   167,  -224,   148,   139,    -5,    83,  -224,   139,   139,
      40,  -224,   133,  -224,   169,    49,   158,   139,    86,    66,
    -224,   139,   139,   -13,   168,  -224,   173,   139,  -224,  -224,
     121,  -224,   139,  -224,   139,    83,    83,   139,   105,  -224,
    -224,   170,  -224,  -224,   168,  -224,   121,   313,    83,    83,
    -224,  -224,   132,  -224,   139,    93,  -224,   313,  -224,  -224,
    -224,   139,   133,  -224,   133,   197,   197
};

/* YYPGOTO[NTERM-NUM].  */
static const short int yypgoto[] =
{
    -224,  -224,  -224,     2,  -224,  -224,  -224,  -224,  -224,  -224,
     106,  -224,  -224,  -224,  -224,  -224,  -224,  -224,  -224,  -224,
      16,  -224,  -224,   -24,  -224,  -224,  -224,  -135,    58,  -224,
    -224,  -223,   233,  -224,   -15,  -224,  -212,   -37,     8,   -28,
      92,    82,   -23,    -8,     1,    22,  -224,   -14,  -224,   241,
    -224,  -224,   217,   -34,  -105,  -224,  -224,    46,   182,  -224,
     141,  -224,   184
};

/* YYTABLE[YYPACT[STATE-NUM]].  What to do in state STATE-NUM.  If
   positive, shift that token.  If negative, reduce the rule which
   number is the opposite.  If zero, do what YYDEFACT says.
   If YYTABLE_NINF, syntax error.  */
#define YYTABLE_NINF -1
static const unsigned short int yytable[] =
{
      27,   133,   141,    59,    61,    70,     1,    27,     2,    28,
     233,    29,    27,    78,    57,   141,    51,     1,    14,     2,
      82,    28,   260,   261,   179,   180,    60,   142,    70,    21,
     190,   191,    95,   192,    82,   268,   269,   159,    52,    72,
     108,     1,   226,     2,   115,    27,   116,   227,    15,    27,
      59,    61,    22,     3,    79,   181,    70,   113,    84,    27,
     275,    70,   276,   217,    15,   221,   178,    60,    79,    60,
      60,    83,    89,    53,    71,    98,    23,    24,   117,   118,
      99,    67,    25,    26,   100,   125,   126,   230,   119,    27,
     231,    70,   232,    89,    89,    89,   236,     1,   114,     2,
     115,   237,   116,   146,   151,   125,   126,    80,    81,    96,
      27,    23,    24,   243,    92,    70,   244,    25,    26,   138,
     146,   102,    67,   131,   132,   151,    97,    27,    27,   197,
     109,    70,    70,   241,   139,   118,   154,   155,   242,   103,
     101,    23,    24,   110,   119,    32,    33,    25,    26,    23,
      24,    70,   263,   107,   104,    25,    26,   264,   202,    27,
     273,   127,    70,   137,   129,   130,   105,    27,   174,   203,
     204,   205,   206,   207,   208,   209,   185,    70,   196,   270,
      27,   134,   135,   136,   271,   158,   249,   250,   251,    23,
      24,   156,   108,   143,   210,    25,    26,   170,   131,   132,
      85,   169,    86,    87,   171,   172,   238,   197,    27,   175,
     177,    23,    24,   186,   176,    32,    33,    25,    26,   201,
     223,   255,   235,   265,   187,   141,   140,   188,   189,   224,
     266,   183,    27,     8,    70,   165,    69,    27,    27,    68,
     162,   257,   215,   216,    70,   220,    16,    90,    27,    91,
      27,    27,     0,     0,     0,     0,   222,   267,     0,     0,
     225,     0,    23,    24,   228,   229,    32,    33,    25,    26,
       0,     0,     0,   240,     0,     0,     0,   245,   246,   248,
       0,     0,     0,   256,     0,     0,     0,   108,   258,     0,
     259,     0,     0,   262,    57,     0,     0,     1,     0,     2,
       0,    23,    24,     0,     0,    32,    33,    25,    26,     0,
     272,     0,     0,     0,     0,     0,     0,   274,   157,    23,
      24,    30,    31,    32,    33,    25,    26,    34,    15,    36,
      37,   111,    39,    40,    41,    42,    43,     0,     0,     0,
       0,    23,    24,     0,    58,    32,    33,    25,    26,     0,
       0,     0,     0,     0,    23,    24,    30,    31,    32,    33,
      25,    26,    34,    35,    36,    37,    38,    39,    40,    41,
      42,    43,    76,     0,     0,     0,     0,     0,     0,   112,
       0,     0,     0,     0,     0,     0,    23,    24,    30,    31,
      32,    33,    25,    26,    34,    35,    36,    37,    38,    39,
      40,    41,    42,    43,     0,     0,     0,     0,     0,     0,
       0,    77,    23,    24,    30,    31,    32,    33,    25,    26,
      34,    35,    36,    37,    38,    39,    40,    41,    42,    43,
       0,     0,     0,     0,     0,     0,     0,    44,    23,    24,
      30,    31,    32,    33,    25,    26,    34,    35,    36,    37,
      38,    39,    40,    41,    42,    43,     0,     0,     0,     0,
       0,     0,     0,    93,    23,    24,    30,    31,    32,    33,
      25,    26,    34,     0,    36,    37,     0,    39,    40,    41,
      42,    43,     0,     0,     0,     0,     0,     0,     0,    58,
      23,    24,    30,    31,    32,    33,    25,    26,    34,     0,
      36,    37,     0,    39,    40,    41,    42,    43,     0,     0,
       0,     0,     0,     0,     0,    88,    23,    24,   160,   144,
      32,    33,    25,    26,    34,     0,    36,    37,     0,    39,
      40,    41,    42,    43,     0,     0,     0,     0,     0,     0,
       0,   161,    23,    24,   149,   163,    32,    33,    25,    26,
      34,     0,    36,    37,     0,    39,    40,    41,    42,    43,
       0,     0,     0,     0,     0,    23,    24,   164,   144,    32,
      33,    25,    26,    34,     0,    36,    37,     0,    39,    40,
      41,    42,    43,     0,     0,     0,     0,    23,    24,   149,
     145,    32,    33,    25,    26,    34,     0,    36,    37,     0,
      39,    40,    41,    42,    43,     0,     0,     0,     0,    23,
      24,     0,   150,    32,    33,    25,    26,    23,    24,     0,
       0,    32,    33,    25,    26,    23,    24,     0,     0,    32,
      33,    25,    26,     0,   166,     0,     0,     0,     0,    23,
      24,     0,   167,    32,    33,    25,    26,    23,    24,     0,
     184,    32,    33,    25,    26,     0,     0,     0,     0,     0,
       0,     0,     0,     0,   195,     0,     0,   190,   191,     0,
     192
};

static const short int yycheck[] =
{
       8,   106,    28,    18,    18,    28,     6,    15,     8,     8,
     222,     9,    20,    50,     3,    28,    15,     6,     0,     8,
      54,    20,   245,   246,    15,    16,    18,    53,    51,     7,
      48,    49,    69,    51,    68,   258,   259,   142,    16,     7,
      53,     6,    47,     8,     9,    53,    11,    52,    37,    57,
      65,    65,     4,    53,    53,    46,    79,    94,    57,    67,
     272,    84,   274,   198,    37,   200,   171,    59,    67,    61,
      62,    38,    64,    53,    28,    73,    28,    29,    43,    44,
      38,    53,    34,    35,    42,   100,   100,    47,    53,    97,
      50,   114,    52,    85,    86,    87,    47,     6,    97,     8,
       9,    52,    11,   131,   132,   120,   120,    38,    39,     5,
     118,    28,    29,    47,    68,   138,    50,    34,    35,   118,
     148,    75,    53,    30,    31,   153,    53,   135,   136,    46,
      84,   154,   155,    47,    43,    44,   135,   136,    52,    38,
      53,    28,    29,     5,    53,    32,    33,    34,    35,    28,
      29,   174,    47,    38,    41,    34,    35,    52,    10,   167,
     265,    28,   185,   117,    53,    41,    53,   175,   167,    21,
      22,    23,    24,    25,    26,    27,   175,   200,    17,    47,
     188,    38,    53,    53,    52,   139,    18,    19,    20,    28,
      29,    38,    53,    53,    46,    34,    35,    38,    30,    31,
      59,   155,    61,    62,    53,   159,    48,    46,   216,    53,
      38,    28,    29,    38,   168,    32,    33,    34,    35,    38,
      53,    48,    53,    53,   178,    28,   120,   181,   182,   213,
     254,   173,   240,     0,   257,   153,    53,   245,   246,    22,
     148,   240,   196,   197,   267,   199,     5,    65,   256,    65,
     258,   259,    -1,    -1,    -1,    -1,   210,   256,    -1,    -1,
     214,    -1,    28,    29,   218,   219,    32,    33,    34,    35,
      -1,    -1,    -1,   227,    -1,    -1,    -1,   231,   232,   233,
      -1,    -1,    -1,   237,    -1,    -1,    -1,    53,   242,    -1,
     244,    -1,    -1,   247,     3,    -1,    -1,     6,    -1,     8,
      -1,    28,    29,    -1,    -1,    32,    33,    34,    35,    -1,
     264,    -1,    -1,    -1,    -1,    -1,    -1,   271,    45,    28,
      29,    30,    31,    32,    33,    34,    35,    36,    37,    38,
      39,     5,    41,    42,    43,    44,    45,    -1,    -1,    -1,
      -1,    28,    29,    -1,    53,    32,    33,    34,    35,    -1,
      -1,    -1,    -1,    -1,    28,    29,    30,    31,    32,    33,
      34,    35,    36,    37,    38,    39,    40,    41,    42,    43,
      44,    45,    14,    -1,    -1,    -1,    -1,    -1,    -1,    53,
      -1,    -1,    -1,    -1,    -1,    -1,    28,    29,    30,    31,
      32,    33,    34,    35,    36,    37,    38,    39,    40,    41,
      42,    43,    44,    45,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    53,    28,    29,    30,    31,    32,    33,    34,    35,
      36,    37,    38,    39,    40,    41,    42,    43,    44,    45,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    53,    28,    29,
      30,    31,    32,    33,    34,    35,    36,    37,    38,    39,
      40,    41,    42,    43,    44,    45,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    53,    28,    29,    30,    31,    32,    33,
      34,    35,    36,    -1,    38,    39,    -1,    41,    42,    43,
      44,    45,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    53,
      28,    29,    30,    31,    32,    33,    34,    35,    36,    -1,
      38,    39,    -1,    41,    42,    43,    44,    45,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    53,    28,    29,    30,    31,
      32,    33,    34,    35,    36,    -1,    38,    39,    -1,    41,
      42,    43,    44,    45,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    53,    28,    29,    30,    31,    32,    33,    34,    35,
      36,    -1,    38,    39,    -1,    41,    42,    43,    44,    45,
      -1,    -1,    -1,    -1,    -1,    28,    29,    53,    31,    32,
      33,    34,    35,    36,    -1,    38,    39,    -1,    41,    42,
      43,    44,    45,    -1,    -1,    -1,    -1,    28,    29,    30,
      53,    32,    33,    34,    35,    36,    -1,    38,    39,    -1,
      41,    42,    43,    44,    45,    -1,    -1,    -1,    -1,    28,
      29,    -1,    53,    32,    33,    34,    35,    28,    29,    -1,
      -1,    32,    33,    34,    35,    28,    29,    -1,    -1,    32,
      33,    34,    35,    -1,    53,    -1,    -1,    -1,    -1,    28,
      29,    -1,    53,    32,    33,    34,    35,    28,    29,    -1,
      53,    32,    33,    34,    35,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    53,    -1,    -1,    48,    49,    -1,
      51
};

/* YYSTOS[STATE-NUM] -- The (internal number of the) accessing
   symbol of state STATE-NUM.  */
static const unsigned char yystos[] =
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
      47,    50,    52,    90,    71,    53,    47,    52,    48,    78,
     111,    47,    52,    47,    50,   111,   111,    76,   111,    18,
      19,    20,    77,   108,    73,    48,   111,    98,   111,   111,
      85,    85,   111,    47,    52,    53,    77,    98,    85,    85,
      47,    52,   111,   108,   111,    90,    90
};

#if ! defined (YYSIZE_T) && defined (__SIZE_TYPE__)
# define YYSIZE_T __SIZE_TYPE__
#endif
#if ! defined (YYSIZE_T) && defined (size_t)
# define YYSIZE_T size_t
#endif
#if ! defined (YYSIZE_T)
# if defined (__STDC__) || defined (__cplusplus)
#  include <stddef.h> /* INFRINGES ON USER NAME SPACE */
#  define YYSIZE_T size_t
# endif
#endif
#if ! defined (YYSIZE_T)
# define YYSIZE_T unsigned int
#endif

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
      YYPOPSTACK;						\
      goto yybackup;						\
    }								\
  else								\
    { 								\
      yyerror ("syntax error: cannot back up");\
      YYERROR;							\
    }								\
while (0)


#define YYTERROR	1
#define YYERRCODE	256


/* YYLLOC_DEFAULT -- Set CURRENT to span from RHS[1] to RHS[N].
   If N is 0, then set CURRENT to the empty location which ends
   the previous symbol: RHS[0] (always defined).  */

#define YYRHSLOC(Rhs, K) ((Rhs)[K])
#ifndef YYLLOC_DEFAULT
# define YYLLOC_DEFAULT(Current, Rhs, N)				\
    do									\
      if (N)								\
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
    while (0)
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
} while (0)

# define YY_SYMBOL_PRINT(Title, Type, Value, Location)		\
do {								\
  if (yydebug)							\
    {								\
      YYFPRINTF (stderr, "%s ", Title);				\
      yysymprint (stderr, 					\
                  Type, Value);	\
      YYFPRINTF (stderr, "\n");					\
    }								\
} while (0)

/*------------------------------------------------------------------.
| yy_stack_print -- Print the state stack from its BOTTOM up to its |
| TOP (included).                                                   |
`------------------------------------------------------------------*/

#if defined (__STDC__) || defined (__cplusplus)
static void
yy_stack_print (short int *bottom, short int *top)
#else
static void
yy_stack_print (bottom, top)
    short int *bottom;
    short int *top;
#endif
{
  YYFPRINTF (stderr, "Stack now");
  for (/* Nothing. */; bottom <= top; ++bottom)
    YYFPRINTF (stderr, " %d", *bottom);
  YYFPRINTF (stderr, "\n");
}

# define YY_STACK_PRINT(Bottom, Top)				\
do {								\
  if (yydebug)							\
    yy_stack_print ((Bottom), (Top));				\
} while (0)


/*------------------------------------------------.
| Report that the YYRULE is going to be reduced.  |
`------------------------------------------------*/

#if defined (__STDC__) || defined (__cplusplus)
static void
yy_reduce_print (int yyrule)
#else
static void
yy_reduce_print (yyrule)
    int yyrule;
#endif
{
  int yyi;
  unsigned int yylno = yyrline[yyrule];
  YYFPRINTF (stderr, "Reducing stack by rule %d (line %u), ",
             yyrule - 1, yylno);
  /* Print the symbols being reduced, and their result.  */
  for (yyi = yyprhs[yyrule]; 0 <= yyrhs[yyi]; yyi++)
    YYFPRINTF (stderr, "%s ", yytname [yyrhs[yyi]]);
  YYFPRINTF (stderr, "-> %s\n", yytname [yyr1[yyrule]]);
}

# define YY_REDUCE_PRINT(Rule)		\
do {					\
  if (yydebug)				\
    yy_reduce_print (Rule);		\
} while (0)

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
   SIZE_MAX < YYSTACK_BYTES (YYMAXDEPTH)
   evaluated with infinite-precision integer arithmetic.  */

#ifndef YYMAXDEPTH
# define YYMAXDEPTH 10000
#endif



#if YYERROR_VERBOSE

# ifndef yystrlen
#  if defined (__GLIBC__) && defined (_STRING_H)
#   define yystrlen strlen
#  else
/* Return the length of YYSTR.  */
static YYSIZE_T
#   if defined (__STDC__) || defined (__cplusplus)
yystrlen (const char *yystr)
#   else
yystrlen (yystr)
     const char *yystr;
#   endif
{
  register const char *yys = yystr;

  while (*yys++ != '\0')
    continue;

  return yys - yystr - 1;
}
#  endif
# endif

# ifndef yystpcpy
#  if defined (__GLIBC__) && defined (_STRING_H) && defined (_GNU_SOURCE)
#   define yystpcpy stpcpy
#  else
/* Copy YYSRC to YYDEST, returning the address of the terminating '\0' in
   YYDEST.  */
static char *
#   if defined (__STDC__) || defined (__cplusplus)
yystpcpy (char *yydest, const char *yysrc)
#   else
yystpcpy (yydest, yysrc)
     char *yydest;
     const char *yysrc;
#   endif
{
  register char *yyd = yydest;
  register const char *yys = yysrc;

  while ((*yyd++ = *yys++) != '\0')
    continue;

  return yyd - 1;
}
#  endif
# endif

#endif /* !YYERROR_VERBOSE */



#if YYDEBUG
/*--------------------------------.
| Print this symbol on YYOUTPUT.  |
`--------------------------------*/

#if defined (__STDC__) || defined (__cplusplus)
static void
yysymprint (FILE *yyoutput, int yytype, YYSTYPE *yyvaluep)
#else
static void
yysymprint (yyoutput, yytype, yyvaluep)
    FILE *yyoutput;
    int yytype;
    YYSTYPE *yyvaluep;
#endif
{
  /* Pacify ``unused variable'' warnings.  */
  (void) yyvaluep;

  if (yytype < YYNTOKENS)
    YYFPRINTF (yyoutput, "token %s (", yytname[yytype]);
  else
    YYFPRINTF (yyoutput, "nterm %s (", yytname[yytype]);


# ifdef YYPRINT
  if (yytype < YYNTOKENS)
    YYPRINT (yyoutput, yytoknum[yytype], *yyvaluep);
# endif
  switch (yytype)
    {
      default:
        break;
    }
  YYFPRINTF (yyoutput, ")");
}

#endif /* ! YYDEBUG */
/*-----------------------------------------------.
| Release the memory associated to this symbol.  |
`-----------------------------------------------*/

#if defined (__STDC__) || defined (__cplusplus)
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
  /* Pacify ``unused variable'' warnings.  */
  (void) yyvaluep;

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
# if defined (__STDC__) || defined (__cplusplus)
int yyparse (void *YYPARSE_PARAM);
# else
int yyparse ();
# endif
#else /* ! YYPARSE_PARAM */
#if defined (__STDC__) || defined (__cplusplus)
int yyparse (void);
#else
int yyparse ();
#endif
#endif /* ! YYPARSE_PARAM */






/*----------.
| yyparse.  |
`----------*/

#ifdef YYPARSE_PARAM
# if defined (__STDC__) || defined (__cplusplus)
int yyparse (void *YYPARSE_PARAM)
# else
int yyparse (YYPARSE_PARAM)
  void *YYPARSE_PARAM;
# endif
#else /* ! YYPARSE_PARAM */
#if defined (__STDC__) || defined (__cplusplus)
int
yyparse (void)
#else
int
yyparse ()

#endif
#endif
{
  /* The look-ahead symbol.  */
int yychar;

/* The semantic value of the look-ahead symbol.  */
YYSTYPE yylval;

/* Number of syntax errors so far.  */
int yynerrs;

  register int yystate;
  register int yyn;
  int yyresult;
  /* Number of tokens to shift before error messages enabled.  */
  int yyerrstatus;
  /* Look-ahead token as an internal (translated) token number.  */
  int yytoken = 0;

  /* Three stacks and their tools:
     `yyss': related to states,
     `yyvs': related to semantic values,
     `yyls': related to locations.

     Refer to the stacks thru separate pointers, to allow yyoverflow
     to reallocate them elsewhere.  */

  /* The state stack.  */
  short int yyssa[YYINITDEPTH];
  short int *yyss = yyssa;
  register short int *yyssp;

  /* The semantic value stack.  */
  YYSTYPE yyvsa[YYINITDEPTH];
  YYSTYPE *yyvs = yyvsa;
  register YYSTYPE *yyvsp;



#define YYPOPSTACK   (yyvsp--, yyssp--)

  YYSIZE_T yystacksize = YYINITDEPTH;

  /* The variables used to return semantic value and location from the
     action routines.  */
  YYSTYPE yyval;


  /* When reducing, the number of symbols on the RHS of the reduced
     rule.  */
  int yylen;

  YYDPRINTF ((stderr, "Starting parse\n"));

  yystate = 0;
  yyerrstatus = 0;
  yynerrs = 0;
  yychar = YYEMPTY;		/* Cause a token to be read.  */

  /* Initialize stack pointers.
     Waste one element of value and location stack
     so that they stay on the same level as the state stack.
     The wasted elements are never initialized.  */

  yyssp = yyss;
  yyvsp = yyvs;


  yyvsp[0] = yylval;

  goto yysetstate;

/*------------------------------------------------------------.
| yynewstate -- Push a new state, which is found in yystate.  |
`------------------------------------------------------------*/
 yynewstate:
  /* In all cases, when you get here, the value and location stacks
     have just been pushed. so pushing a state here evens the stacks.
     */
  yyssp++;

 yysetstate:
  *yyssp = yystate;

  if (yyss + yystacksize - 1 <= yyssp)
    {
      /* Get the current used size of the three stacks, in elements.  */
      YYSIZE_T yysize = yyssp - yyss + 1;

#ifdef yyoverflow
      {
	/* Give user a chance to reallocate the stack. Use copies of
	   these so that the &'s don't force the real ones into
	   memory.  */
	YYSTYPE *yyvs1 = yyvs;
	short int *yyss1 = yyss;


	/* Each stack pointer address is followed by the size of the
	   data in use in that stack, in bytes.  This used to be a
	   conditional around just the two extra args, but that might
	   be undefined if yyoverflow is a macro.  */
	yyoverflow ("parser stack overflow",
		    &yyss1, yysize * sizeof (*yyssp),
		    &yyvs1, yysize * sizeof (*yyvsp),

		    &yystacksize);

	yyss = yyss1;
	yyvs = yyvs1;
      }
#else /* no yyoverflow */
# ifndef YYSTACK_RELOCATE
      goto yyoverflowlab;
# else
      /* Extend the stack our own way.  */
      if (YYMAXDEPTH <= yystacksize)
	goto yyoverflowlab;
      yystacksize *= 2;
      if (YYMAXDEPTH < yystacksize)
	yystacksize = YYMAXDEPTH;

      {
	short int *yyss1 = yyss;
	union yyalloc *yyptr =
	  (union yyalloc *) YYSTACK_ALLOC (YYSTACK_BYTES (yystacksize));
	if (! yyptr)
	  goto yyoverflowlab;
	YYSTACK_RELOCATE (yyss);
	YYSTACK_RELOCATE (yyvs);

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

  goto yybackup;

/*-----------.
| yybackup.  |
`-----------*/
yybackup:

/* Do appropriate processing given the current state.  */
/* Read a look-ahead token if we need one and don't already have one.  */
/* yyresume: */

  /* First try to decide what to do without reference to look-ahead token.  */

  yyn = yypact[yystate];
  if (yyn == YYPACT_NINF)
    goto yydefault;

  /* Not known => get a look-ahead token if don't already have one.  */

  /* YYCHAR is either YYEMPTY or YYEOF or a valid look-ahead symbol.  */
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

  if (yyn == YYFINAL)
    YYACCEPT;

  /* Shift the look-ahead token.  */
  YY_SYMBOL_PRINT ("Shifting", yytoken, &yylval, &yylloc);

  /* Discard the token being shifted unless it is eof.  */
  if (yychar != YYEOF)
    yychar = YYEMPTY;

  *++yyvsp = yylval;


  /* Count tokens shifted since error; after three, turn off error
     status.  */
  if (yyerrstatus)
    yyerrstatus--;

  yystate = yyn;
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
#line 93 "xml_parse.y"
    {;}
    break;

  case 11:
#line 94 "xml_parse.y"
    {free_string((yyvsp[-1].str));;}
    break;

  case 12:
#line 95 "xml_parse.y"
    {free_string((yyvsp[-2].str)); free_string((yyvsp[0].str));;}
    break;

  case 13:
#line 96 "xml_parse.y"
    {free_string((yyvsp[0].str));;}
    break;

  case 14:
#line 98 "xml_parse.y"
    {;}
    break;

  case 17:
#line 102 "xml_parse.y"
    {;}
    break;

  case 18:
#line 103 "xml_parse.y"
    {if(!check_dtd_tok((yyvsp[-3].str),"SYSTEM")) printf("System '%s'\n",get_cstring((yyvsp[-1].str)));;}
    break;

  case 19:
#line 104 "xml_parse.y"
    {if(!check_dtd_tok((yyvsp[-5].str),"PUBLIC")) printf("Public '%s' '%s'\n",get_cstring((yyvsp[-3].str)),get_cstring((yyvsp[-1].str)));;}
    break;

  case 20:
#line 106 "xml_parse.y"
    {printf("DTD name %s\n",get_cstring((yyvsp[0].str)));;}
    break;

  case 21:
#line 108 "xml_parse.y"
    {;}
    break;

  case 22:
#line 109 "xml_parse.y"
    {;}
    break;

  case 23:
#line 111 "xml_parse.y"
    {;}
    break;

  case 24:
#line 112 "xml_parse.y"
    {;}
    break;

  case 25:
#line 114 "xml_parse.y"
    {;}
    break;

  case 26:
#line 115 "xml_parse.y"
    {;}
    break;

  case 27:
#line 117 "xml_parse.y"
    {;}
    break;

  case 28:
#line 119 "xml_parse.y"
    {element_state=1;;}
    break;

  case 29:
#line 119 "xml_parse.y"
    {element_state=0;;}
    break;

  case 32:
#line 122 "xml_parse.y"
    {;}
    break;

  case 33:
#line 123 "xml_parse.y"
    {;}
    break;

  case 34:
#line 125 "xml_parse.y"
    {attlist_state=1;;}
    break;

  case 35:
#line 125 "xml_parse.y"
    {attlist_state=2;;}
    break;

  case 36:
#line 125 "xml_parse.y"
    {attlist_state=0;;}
    break;

  case 37:
#line 126 "xml_parse.y"
    {attlist_state=1;;}
    break;

  case 38:
#line 126 "xml_parse.y"
    {attlist_state=2;;}
    break;

  case 39:
#line 126 "xml_parse.y"
    {attlist_state=0;;}
    break;

  case 40:
#line 128 "xml_parse.y"
    {;}
    break;

  case 41:
#line 129 "xml_parse.y"
    {;}
    break;

  case 42:
#line 130 "xml_parse.y"
    {;}
    break;

  case 43:
#line 131 "xml_parse.y"
    {;}
    break;

  case 44:
#line 132 "xml_parse.y"
    {;}
    break;

  case 45:
#line 133 "xml_parse.y"
    {;}
    break;

  case 46:
#line 134 "xml_parse.y"
    {;}
    break;

  case 47:
#line 135 "xml_parse.y"
    {;}
    break;

  case 48:
#line 136 "xml_parse.y"
    {;}
    break;

  case 49:
#line 138 "xml_parse.y"
    {;}
    break;

  case 50:
#line 139 "xml_parse.y"
    {;}
    break;

  case 51:
#line 141 "xml_parse.y"
    {;}
    break;

  case 52:
#line 142 "xml_parse.y"
    {;}
    break;

  case 53:
#line 144 "xml_parse.y"
    {;}
    break;

  case 54:
#line 145 "xml_parse.y"
    {;}
    break;

  case 55:
#line 146 "xml_parse.y"
    {;}
    break;

  case 56:
#line 147 "xml_parse.y"
    {;}
    break;

  case 58:
#line 150 "xml_parse.y"
    {;}
    break;

  case 59:
#line 152 "xml_parse.y"
    {;}
    break;

  case 60:
#line 153 "xml_parse.y"
    {;}
    break;

  case 61:
#line 154 "xml_parse.y"
    {;}
    break;

  case 62:
#line 155 "xml_parse.y"
    {;}
    break;

  case 64:
#line 158 "xml_parse.y"
    {;}
    break;

  case 65:
#line 159 "xml_parse.y"
    {;}
    break;

  case 67:
#line 162 "xml_parse.y"
    {;}
    break;

  case 68:
#line 163 "xml_parse.y"
    {;}
    break;

  case 69:
#line 164 "xml_parse.y"
    {;}
    break;

  case 70:
#line 166 "xml_parse.y"
    {;}
    break;

  case 71:
#line 167 "xml_parse.y"
    {;}
    break;

  case 72:
#line 168 "xml_parse.y"
    {;}
    break;

  case 73:
#line 170 "xml_parse.y"
    {;}
    break;

  case 74:
#line 171 "xml_parse.y"
    {;}
    break;

  case 75:
#line 173 "xml_parse.y"
    {;}
    break;

  case 76:
#line 174 "xml_parse.y"
    {;}
    break;

  case 77:
#line 176 "xml_parse.y"
    {;}
    break;

  case 78:
#line 177 "xml_parse.y"
    {;}
    break;

  case 79:
#line 179 "xml_parse.y"
    {pi_state=1;;}
    break;

  case 80:
#line 181 "xml_parse.y"
    {pi_state=0; check_declaration((yyvsp[-2].arg));;}
    break;

  case 81:
#line 183 "xml_parse.y"
    {pi_state=0; free_string((yyvsp[-2].str)); if(call->pi) call->pi((yyvsp[-3].str),(yyvsp[-1].str)); free_string((yyvsp[-3].str)); free_string((yyvsp[-1].str)); ;}
    break;

  case 82:
#line 184 "xml_parse.y"
    {pi_state=0; if(call->pi) call->pi((yyvsp[-2].str),0); free_string((yyvsp[-2].str)); ;}
    break;

  case 83:
#line 186 "xml_parse.y"
    {(yyval.str)=add_to_string(0,(yyvsp[0].c));;}
    break;

  case 85:
#line 188 "xml_parse.y"
    {(yyval.str)=add_to_string((yyvsp[-1].str),(yyvsp[0].c));;}
    break;

  case 86:
#line 189 "xml_parse.y"
    {(yyval.str)=add_strings((yyvsp[-1].str),(yyvsp[0].str));;}
    break;

  case 87:
#line 191 "xml_parse.y"
    {(yyval.str)=add_to_string(0,(yyvsp[0].c));;}
    break;

  case 88:
#line 192 "xml_parse.y"
    {(yyval.str)=add_to_string((yyvsp[-1].str),(yyvsp[0].c));;}
    break;

  case 115:
#line 227 "xml_parse.y"
    {(yyval.str)=add_to_string(0,(yyvsp[0].c));;}
    break;

  case 116:
#line 228 "xml_parse.y"
    {(yyval.str)=add_to_string((yyvsp[-1].str),(yyvsp[0].c));;}
    break;

  case 117:
#line 230 "xml_parse.y"
    {;}
    break;

  case 118:
#line 231 "xml_parse.y"
    {;}
    break;

  case 119:
#line 232 "xml_parse.y"
    {free_string((yyvsp[0].str));;}
    break;

  case 120:
#line 234 "xml_parse.y"
    {(yyval.str)=add_to_string(0,(yyvsp[0].c));;}
    break;

  case 122:
#line 236 "xml_parse.y"
    {(yyval.str)=add_to_string((yyvsp[-1].str),(yyvsp[0].c));;}
    break;

  case 123:
#line 237 "xml_parse.y"
    {(yyval.str)=add_strings((yyvsp[-1].str),(yyvsp[0].str));;}
    break;

  case 124:
#line 239 "xml_parse.y"
    {in_comment=1;;}
    break;

  case 125:
#line 239 "xml_parse.y"
    {in_comment=0; if(call->comment) call->comment((yyvsp[-2].str)); free_string((yyvsp[-2].str));;}
    break;

  case 129:
#line 245 "xml_parse.y"
    {(yyval.str)=(yyvsp[0].str);;}
    break;

  case 130:
#line 247 "xml_parse.y"
    {check_start((yyvsp[-2].str),(yyvsp[-1].arg)); ;}
    break;

  case 131:
#line 248 "xml_parse.y"
    {check_start((yyvsp[-2].str),0); ;}
    break;

  case 133:
#line 251 "xml_parse.y"
    {(yyvsp[0].arg)->next=(yyvsp[-1].arg); (yyval.arg)=(yyvsp[0].arg);;}
    break;

  case 134:
#line 253 "xml_parse.y"
    {free_string((yyvsp[-3].str)); (yyval.arg)=make_arg((yyvsp[-2].str),(yyvsp[0].str));;}
    break;

  case 135:
#line 255 "xml_parse.y"
    {(yyval.str)=(yyvsp[-1].str);;}
    break;

  case 136:
#line 256 "xml_parse.y"
    {(yyval.str)=(yyvsp[-1].str);;}
    break;

  case 137:
#line 258 "xml_parse.y"
    {(yyval.str)=add_to_string(0,(yyvsp[0].c));;}
    break;

  case 139:
#line 260 "xml_parse.y"
    {(yyval.str)=add_to_string((yyvsp[-1].str),(yyvsp[0].c));;}
    break;

  case 140:
#line 261 "xml_parse.y"
    {(yyval.str)=add_strings((yyvsp[-1].str),(yyvsp[0].str));;}
    break;

  case 141:
#line 263 "xml_parse.y"
    {(yyval.str)=add_to_string(0,(yyvsp[0].c));;}
    break;

  case 143:
#line 265 "xml_parse.y"
    {(yyval.str)=add_to_string((yyvsp[-1].str),(yyvsp[0].c));;}
    break;

  case 144:
#line 266 "xml_parse.y"
    {(yyval.str)=add_strings((yyvsp[-1].str),(yyvsp[0].str));;}
    break;

  case 146:
#line 269 "xml_parse.y"
    {free_string((yyvsp[0].str));;}
    break;

  case 147:
#line 271 "xml_parse.y"
    {check_end((yyvsp[-2].str),0);;}
    break;

  case 148:
#line 273 "xml_parse.y"
    {check_start((yyvsp[-3].str),(yyvsp[-2].arg)); check_end((yyvsp[-3].str),1);;}
    break;

  case 149:
#line 275 "xml_parse.y"
    {(yyval.str)=add_to_string(0,(yyvsp[0].c));;}
    break;

  case 151:
#line 277 "xml_parse.y"
    {(yyval.str)=add_to_string((yyvsp[-1].str),(yyvsp[0].c));;}
    break;

  case 152:
#line 278 "xml_parse.y"
    {(yyval.str)=add_strings((yyvsp[-1].str),(yyvsp[0].str));;}
    break;

  case 153:
#line 280 "xml_parse.y"
    { if(check_ws((yyvsp[0].str)) && call->content) call->content((yyvsp[0].str)); free_string((yyvsp[0].str)); ;}
    break;

  case 154:
#line 281 "xml_parse.y"
    { if(check_ws((yyvsp[0].str)) && call->content) call->content((yyvsp[0].str)); free_string((yyvsp[0].str)); ;}
    break;

  case 155:
#line 282 "xml_parse.y"
    { if(check_ws((yyvsp[0].str)) && call->content) call->content((yyvsp[0].str)); free_string((yyvsp[0].str)); ;}
    break;

  case 156:
#line 284 "xml_parse.y"
    {(yyval.str)=0;;}
    break;

  case 157:
#line 285 "xml_parse.y"
    {(yyval.str)=(yyvsp[0].str);;}
    break;

  case 158:
#line 286 "xml_parse.y"
    {(yyval.str)=0;;}
    break;

  case 159:
#line 287 "xml_parse.y"
    {(yyval.str)=(yyvsp[0].str);;}
    break;

  case 160:
#line 288 "xml_parse.y"
    {(yyval.str)=0;;}
    break;

  case 161:
#line 289 "xml_parse.y"
    {(yyval.str)=(yyvsp[0].str);;}
    break;


    }

/* Line 1037 of yacc.c.  */
#line 2067 "xml_parse.tab.c"

  yyvsp -= yylen;
  yyssp -= yylen;


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
#if YYERROR_VERBOSE
      yyn = yypact[yystate];

      if (YYPACT_NINF < yyn && yyn < YYLAST)
	{
	  YYSIZE_T yysize = 0;
	  int yytype = YYTRANSLATE (yychar);
	  const char* yyprefix;
	  char *yymsg;
	  int yyx;

	  /* Start YYX at -YYN if negative to avoid negative indexes in
	     YYCHECK.  */
	  int yyxbegin = yyn < 0 ? -yyn : 0;

	  /* Stay within bounds of both yycheck and yytname.  */
	  int yychecklim = YYLAST - yyn;
	  int yyxend = yychecklim < YYNTOKENS ? yychecklim : YYNTOKENS;
	  int yycount = 0;

	  yyprefix = ", expecting ";
	  for (yyx = yyxbegin; yyx < yyxend; ++yyx)
	    if (yycheck[yyx + yyn] == yyx && yyx != YYTERROR)
	      {
		yysize += yystrlen (yyprefix) + yystrlen (yytname [yyx]);
		yycount += 1;
		if (yycount == 5)
		  {
		    yysize = 0;
		    break;
		  }
	      }
	  yysize += (sizeof ("syntax error, unexpected ")
		     + yystrlen (yytname[yytype]));
	  yymsg = (char *) YYSTACK_ALLOC (yysize);
	  if (yymsg != 0)
	    {
	      char *yyp = yystpcpy (yymsg, "syntax error, unexpected ");
	      yyp = yystpcpy (yyp, yytname[yytype]);

	      if (yycount < 5)
		{
		  yyprefix = ", expecting ";
		  for (yyx = yyxbegin; yyx < yyxend; ++yyx)
		    if (yycheck[yyx + yyn] == yyx && yyx != YYTERROR)
		      {
			yyp = yystpcpy (yyp, yyprefix);
			yyp = yystpcpy (yyp, yytname[yyx]);
			yyprefix = " or ";
		      }
		}
	      yyerror (yymsg);
	      YYSTACK_FREE (yymsg);
	    }
	  else
	    yyerror ("syntax error; also virtual memory exhausted");
	}
      else
#endif /* YYERROR_VERBOSE */
	yyerror ("syntax error");
    }



  if (yyerrstatus == 3)
    {
      /* If just tried and failed to reuse look-ahead token after an
	 error, discard it.  */

      if (yychar <= YYEOF)
        {
          /* If at end of input, pop the error token,
	     then the rest of the stack, then return failure.  */
	  if (yychar == YYEOF)
	     for (;;)
	       {

		 YYPOPSTACK;
		 if (yyssp == yyss)
		   YYABORT;
		 yydestruct ("Error: popping",
                             yystos[*yyssp], yyvsp);
	       }
        }
      else
	{
	  yydestruct ("Error: discarding", yytoken, &yylval);
	  yychar = YYEMPTY;
	}
    }

  /* Else will try to reuse look-ahead token after shifting the error
     token.  */
  goto yyerrlab1;


/*---------------------------------------------------.
| yyerrorlab -- error raised explicitly by YYERROR.  |
`---------------------------------------------------*/
yyerrorlab:

#ifdef __GNUC__
  /* Pacify GCC when the user code never invokes YYERROR and the label
     yyerrorlab therefore never appears in user code.  */
  if (0)
     goto yyerrorlab;
#endif

yyvsp -= yylen;
  yyssp -= yylen;
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


      yydestruct ("Error: popping", yystos[yystate], yyvsp);
      YYPOPSTACK;
      yystate = *yyssp;
      YY_STACK_PRINT (yyss, yyssp);
    }

  if (yyn == YYFINAL)
    YYACCEPT;

  *++yyvsp = yylval;


  /* Shift the error token. */
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
  yydestruct ("Error: discarding lookahead",
              yytoken, &yylval);
  yychar = YYEMPTY;
  yyresult = 1;
  goto yyreturn;

#ifndef yyoverflow
/*----------------------------------------------.
| yyoverflowlab -- parser overflow comes here.  |
`----------------------------------------------*/
yyoverflowlab:
  yyerror ("parser stack overflow");
  yyresult = 2;
  /* Fall through.  */
#endif

yyreturn:
#ifndef yyoverflow
  if (yyss != yyssa)
    YYSTACK_FREE (yyss);
#endif
  return yyresult;
}


#line 291 "xml_parse.y"


static int check_ws(string *s)
{
  int i=0;
  char *p;

  if(s) {
    p=get_cstring(s);
    while(*p) {
      if(!isspace(*p)) break;
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
  int i,j,c1,err=0,buflen=16;
  
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
      attlist_state=0;
    }
  } else if(attlist_state==2 && c=='#') {
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


