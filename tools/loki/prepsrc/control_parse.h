/* A Bison parser, made by GNU Bison 1.875.  */

/* Skeleton parser for Yacc-like parsing with Bison,
   Copyright (C) 1984, 1989, 1990, 2000, 2001, 2002 Free Software Foundation, Inc.

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




#if ! defined (YYSTYPE) && ! defined (YYSTYPE_IS_DECLARED)
#line 48 "control_parse.y"
typedef union YYSTYPE {
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
} YYSTYPE;
/* Line 1248 of yacc.c.  */
#line 170 "y.tab.h"
# define yystype YYSTYPE /* obsolescent; will be withdrawn */
# define YYSTYPE_IS_DECLARED 1
# define YYSTYPE_IS_TRIVIAL 1
#endif

extern YYSTYPE yylval;



