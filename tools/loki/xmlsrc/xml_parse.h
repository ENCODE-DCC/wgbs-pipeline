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
     END = 258,
     XMLTOK = 259,
     DOUBLEDASH = 260,
     QCLOSE = 261,
     QOPEN = 262,
     BANGOPEN = 263,
     LETTER = 264,
     DIGIT = 265,
     SQ = 266,
     DQ = 267,
     DASH = 268,
     DOT = 269,
     USCORE = 270,
     COLON = 271,
     CHAR = 272,
     OPEN = 273,
     CLOSE = 274,
     SLASH = 275,
     AMPERSAND = 276,
     EQ = 277,
     SPACE = 278
   };
#endif
#define END 258
#define XMLTOK 259
#define DOUBLEDASH 260
#define QCLOSE 261
#define QOPEN 262
#define BANGOPEN 263
#define LETTER 264
#define DIGIT 265
#define SQ 266
#define DQ 267
#define DASH 268
#define DOT 269
#define USCORE 270
#define COLON 271
#define CHAR 272
#define OPEN 273
#define CLOSE 274
#define SLASH 275
#define AMPERSAND 276
#define EQ 277
#define SPACE 278




#if ! defined (YYSTYPE) && ! defined (YYSTYPE_IS_DECLARED)
#line 60 "xml_parse.y"
typedef union YYSTYPE {
	string *str;
	args *arg;
	char c;
	int i;
} YYSTYPE;
/* Line 1248 of yacc.c.  */
#line 89 "y.tab.h"
# define yystype YYSTYPE /* obsolescent; will be withdrawn */
# define YYSTYPE_IS_DECLARED 1
# define YYSTYPE_IS_TRIVIAL 1
#endif





