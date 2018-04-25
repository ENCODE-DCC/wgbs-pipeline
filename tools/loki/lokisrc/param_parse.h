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




#if ! defined (YYSTYPE) && ! defined (YYSTYPE_IS_DECLARED)
#line 87 "param_parse.y"
typedef union YYSTYPE {
	char *string;
	int value;
	double rvalue;
	struct IBD_List *rlist;
	struct lk_variable *lk_var;
	struct num_array *num_array;
	struct clause_atom *clause;
	struct string_list *string_list;
} YYSTYPE;
/* Line 1248 of yacc.c.  */
#line 155 "y.tab.h"
# define yystype YYSTYPE /* obsolescent; will be withdrawn */
# define YYSTYPE_IS_DECLARED 1
# define YYSTYPE_IS_TRIVIAL 1
#endif

extern YYSTYPE yylval;



