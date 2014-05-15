/* A Bison parser, made by GNU Bison 2.3.  */

/* Skeleton interface for Bison's Yacc-like parsers in C

   Copyright (C) 1984, 1989, 1990, 2000, 2001, 2002, 2003, 2004, 2005, 2006
   Free Software Foundation, Inc.

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
   Foundation, Inc., 51 Franklin Street, Fifth Floor,
   Boston, MA 02110-1301, USA.  */

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
     DECLSET = 258,
     DECLPAR = 259,
     DECLVAR = 260,
     DECLMIN = 261,
     DECLMAX = 262,
     DECLSUB = 263,
     DECLSOS = 264,
     DEFNUMB = 265,
     DEFSTRG = 266,
     DEFBOOL = 267,
     DEFSET = 268,
     PRINT = 269,
     CHECK = 270,
     BINARY = 271,
     INTEGER = 272,
     REAL = 273,
     IMPLICIT = 274,
     ASGN = 275,
     DO = 276,
     WITH = 277,
     IN = 278,
     TO = 279,
     UNTIL = 280,
     BY = 281,
     FORALL = 282,
     EXISTS = 283,
     PRIORITY = 284,
     STARTVAL = 285,
     DEFAULT = 286,
     CMP_LE = 287,
     CMP_GE = 288,
     CMP_EQ = 289,
     CMP_LT = 290,
     CMP_GT = 291,
     CMP_NE = 292,
     INFTY = 293,
     AND = 294,
     OR = 295,
     XOR = 296,
     NOT = 297,
     SUM = 298,
     MIN = 299,
     MAX = 300,
     ARGMIN = 301,
     ARGMAX = 302,
     PROD = 303,
     IF = 304,
     THEN = 305,
     ELSE = 306,
     END = 307,
     INTER = 308,
     UNION = 309,
     CROSS = 310,
     SYMDIFF = 311,
     WITHOUT = 312,
     PROJ = 313,
     MOD = 314,
     DIV = 315,
     POW = 316,
     FAC = 317,
     CARD = 318,
     ROUND = 319,
     FLOOR = 320,
     CEIL = 321,
     RANDOM = 322,
     ORD = 323,
     ABS = 324,
     SGN = 325,
     LOG = 326,
     LN = 327,
     EXP = 328,
     SQRT = 329,
     SIN = 330,
     COS = 331,
     TAN = 332,
     POWER = 333,
     SGNPOW = 334,
     READ = 335,
     AS = 336,
     SKIP = 337,
     USE = 338,
     COMMENT = 339,
     MATCH = 340,
     SUBSETS = 341,
     INDEXSET = 342,
     POWERSET = 343,
     VIF = 344,
     VABS = 345,
     TYPE1 = 346,
     TYPE2 = 347,
     LENGTH = 348,
     SUBSTR = 349,
     NUMBSYM = 350,
     STRGSYM = 351,
     VARSYM = 352,
     SETSYM = 353,
     NUMBDEF = 354,
     STRGDEF = 355,
     BOOLDEF = 356,
     SETDEF = 357,
     DEFNAME = 358,
     NAME = 359,
     STRG = 360,
     NUMB = 361,
     SCALE = 362,
     SEPARATE = 363,
     CHECKONLY = 364,
     INDICATOR = 365
   };
#endif
/* Tokens.  */
#define DECLSET 258
#define DECLPAR 259
#define DECLVAR 260
#define DECLMIN 261
#define DECLMAX 262
#define DECLSUB 263
#define DECLSOS 264
#define DEFNUMB 265
#define DEFSTRG 266
#define DEFBOOL 267
#define DEFSET 268
#define PRINT 269
#define CHECK 270
#define BINARY 271
#define INTEGER 272
#define REAL 273
#define IMPLICIT 274
#define ASGN 275
#define DO 276
#define WITH 277
#define IN 278
#define TO 279
#define UNTIL 280
#define BY 281
#define FORALL 282
#define EXISTS 283
#define PRIORITY 284
#define STARTVAL 285
#define DEFAULT 286
#define CMP_LE 287
#define CMP_GE 288
#define CMP_EQ 289
#define CMP_LT 290
#define CMP_GT 291
#define CMP_NE 292
#define INFTY 293
#define AND 294
#define OR 295
#define XOR 296
#define NOT 297
#define SUM 298
#define MIN 299
#define MAX 300
#define ARGMIN 301
#define ARGMAX 302
#define PROD 303
#define IF 304
#define THEN 305
#define ELSE 306
#define END 307
#define INTER 308
#define UNION 309
#define CROSS 310
#define SYMDIFF 311
#define WITHOUT 312
#define PROJ 313
#define MOD 314
#define DIV 315
#define POW 316
#define FAC 317
#define CARD 318
#define ROUND 319
#define FLOOR 320
#define CEIL 321
#define RANDOM 322
#define ORD 323
#define ABS 324
#define SGN 325
#define LOG 326
#define LN 327
#define EXP 328
#define SQRT 329
#define SIN 330
#define COS 331
#define TAN 332
#define POWER 333
#define SGNPOW 334
#define READ 335
#define AS 336
#define SKIP 337
#define USE 338
#define COMMENT 339
#define MATCH 340
#define SUBSETS 341
#define INDEXSET 342
#define POWERSET 343
#define VIF 344
#define VABS 345
#define TYPE1 346
#define TYPE2 347
#define LENGTH 348
#define SUBSTR 349
#define NUMBSYM 350
#define STRGSYM 351
#define VARSYM 352
#define SETSYM 353
#define NUMBDEF 354
#define STRGDEF 355
#define BOOLDEF 356
#define SETDEF 357
#define DEFNAME 358
#define NAME 359
#define STRG 360
#define NUMB 361
#define SCALE 362
#define SEPARATE 363
#define CHECKONLY 364
#define INDICATOR 365




#if ! defined YYSTYPE && ! defined YYSTYPE_IS_DECLARED
typedef union YYSTYPE
#line 79 "src/mmlparse2.y"
{
   unsigned int bits;
   Numb*        numb;
   const char*  strg;
   const char*  name;
   Symbol*      sym;
   Define*      def;
   CodeNode*    code;
}
/* Line 1489 of yacc.c.  */
#line 279 "src/mmlparse2.h"
	YYSTYPE;
# define yystype YYSTYPE /* obsolescent; will be withdrawn */
# define YYSTYPE_IS_DECLARED 1
# define YYSTYPE_IS_TRIVIAL 1
#endif



