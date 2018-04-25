
/* A Bison parser, made by GNU Bison 2.4.1.  */

/* Skeleton interface for Bison's Yacc-like parsers in C
   
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

/* Line 1676 of yacc.c  */
#line 62 "xml_parse.y"

  string *str;
  args *arg;
  char c;
  int i;



/* Line 1676 of yacc.c  */
#line 167 "y.tab.h"
} YYSTYPE;
# define YYSTYPE_IS_TRIVIAL 1
# define yystype YYSTYPE /* obsolescent; will be withdrawn */
# define YYSTYPE_IS_DECLARED 1
#endif




