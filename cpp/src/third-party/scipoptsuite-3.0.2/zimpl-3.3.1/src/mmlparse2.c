/* A Bison parser, made by GNU Bison 2.3.  */

/* Skeleton implementation for Bison's Yacc-like parsers in C

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
#define YYBISON_VERSION "2.3"

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




/* Copy the first part of user declarations.  */
#line 1 "src/mmlparse2.y"

#pragma ident "@(#) $Id: mmlparse2.y,v 1.12 2012/07/29 15:09:27 bzfkocht Exp $"
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*   File....: mmlparse2.y                                                   */
/*   Name....: MML Parser                                                    */
/*   Author..: Thorsten Koch                                                 */
/*   Copyright by Author, All rights reserved                                */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*
 * Copyright (C) 2001-2012 by Thorsten Koch <koch@zib.de>
 * 
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 * 
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 */
/*lint -e428 -e433 -e506 -e514 -e525 -e527 -e537 -e568 -e574 */
/*lint -e639 -e659 -e661 -e662 -e676 -e685 */
/*lint -e713 -e717 -e732 -e734 -e737 -e744 -e750 -e751 -e753 -e762 -e764 -e774 -e778 */
/*lint -e810 -e818 -e825 -e830 */
/*lint -esym(530,yylen) */
/*lint -esym(563,yyerrorlab) */   
/*lint -esym(746,__yy_memcpy) -esym(516,__yy_memcpy) */
/*lint -esym(718,yylex) -esym(746,yylex) */
/*lint -esym(644,yyval,yylval) -esym(645,yylval) -esym(550,yynerrs) */
/*lint -esym(553,__GNUC__)  -esym(578,yylen) */
/*lint -esym(768,bits) -esym(553,YYSTACK_USE_ALLOCA) */
/*lint -esym(593,yymsg) Custodial pointer possibly not freed */
/*lint -esym(426,mem_malloc) call violates semantics */
   
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "bool.h"
#include "mshell.h"
#include "ratlptypes.h"
#include "numb.h"
#include "elem.h"
#include "tuple.h"
#include "mme.h"
#include "set.h"
#include "symbol.h"
#include "entry.h"
#include "idxset.h"
#include "rdefpar.h"
#include "bound.h"
#include "define.h"
#include "mono.h"
#include "term.h"
#include "list.h"
#include "stmt.h"
#include "local.h"
#include "code.h"
#include "inst.h"   
        
#define YYERROR_VERBOSE 1

/*lint -sem(yyerror, 1p && nulterm(1), r_no) */ 
extern void yyerror(const char* s);
 


/* Enabling traces.  */
#ifndef YYDEBUG
# define YYDEBUG 1
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
/* Line 187 of yacc.c.  */
#line 401 "src/mmlparse2.c"
	YYSTYPE;
# define yystype YYSTYPE /* obsolescent; will be withdrawn */
# define YYSTYPE_IS_DECLARED 1
# define YYSTYPE_IS_TRIVIAL 1
#endif



/* Copy the second part of user declarations.  */


/* Line 216 of yacc.c.  */
#line 414 "src/mmlparse2.c"

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
YYID (int i)
#else
static int
YYID (i)
    int i;
#endif
{
  return i;
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
  yytype_int16 yyss;
  YYSTYPE yyvs;
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
# define YYSTACK_RELOCATE(Stack)					\
    do									\
      {									\
	YYSIZE_T yynewbytes;						\
	YYCOPY (&yyptr->Stack, Stack, yysize);				\
	Stack = &yyptr->Stack;						\
	yynewbytes = yystacksize * sizeof (*Stack) + YYSTACK_GAP_MAXIMUM; \
	yyptr += yynewbytes / sizeof (*yyptr);				\
      }									\
    while (YYID (0))

#endif

/* YYFINAL -- State number of the termination state.  */
#define YYFINAL  40
/* YYLAST -- Last index in YYTABLE.  */
#define YYLAST   3142

/* YYNTOKENS -- Number of terminals.  */
#define YYNTOKENS  123
/* YYNNTS -- Number of nonterminals.  */
#define YYNNTS  59
/* YYNRULES -- Number of rules.  */
#define YYNRULES  301
/* YYNRULES -- Number of states.  */
#define YYNSTATES  887

/* YYTRANSLATE(YYLEX) -- Bison symbol number corresponding to YYLEX.  */
#define YYUNDEFTOK  2
#define YYMAXUTOK   365

#define YYTRANSLATE(YYX)						\
  ((unsigned int) (YYX) <= YYMAXUTOK ? yytranslate[YYX] : YYUNDEFTOK)

/* YYTRANSLATE[YYLEX] -- Bison symbol number corresponding to YYLEX.  */
static const yytype_uint8 yytranslate[] =
{
       0,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
     118,   119,   113,   111,   117,   112,     2,   120,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,   114,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,   115,     2,   116,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,   121,     2,   122,     2,     2,     2,     2,
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
      55,    56,    57,    58,    59,    60,    61,    62,    63,    64,
      65,    66,    67,    68,    69,    70,    71,    72,    73,    74,
      75,    76,    77,    78,    79,    80,    81,    82,    83,    84,
      85,    86,    87,    88,    89,    90,    91,    92,    93,    94,
      95,    96,    97,    98,    99,   100,   101,   102,   103,   104,
     105,   106,   107,   108,   109,   110
};

#if YYDEBUG
/* YYPRHS[YYN] -- Index of the first RHS symbol of rule number YYN in
   YYRHS.  */
static const yytype_uint16 yyprhs[] =
{
       0,     0,     3,     5,     7,     9,    11,    13,    15,    17,
      19,    21,    23,    25,    31,    40,    49,    57,    59,    63,
      70,    79,    84,    87,    96,   105,   114,   123,   125,   129,
     139,   148,   154,   156,   158,   160,   161,   164,   174,   181,
     190,   196,   206,   213,   225,   234,   235,   237,   240,   241,
     244,   248,   258,   268,   269,   272,   275,   284,   293,   294,
     297,   298,   301,   303,   307,   309,   312,   315,   319,   323,
     328,   334,   340,   346,   348,   352,   357,   363,   371,   376,
     381,   386,   391,   398,   405,   412,   419,   432,   445,   458,
     471,   484,   497,   510,   523,   536,   549,   562,   575,   588,
     601,   614,   627,   636,   645,   654,   663,   667,   671,   675,
     679,   683,   687,   691,   695,   699,   703,   707,   711,   715,
     719,   723,   727,   731,   735,   739,   743,   747,   750,   754,
     755,   759,   761,   763,   765,   767,   769,   771,   773,   775,
     779,   783,   787,   791,   795,   799,   801,   805,   809,   813,
     817,   819,   822,   825,   827,   831,   836,   839,   844,   849,
     854,   859,   864,   869,   874,   879,   884,   889,   896,   903,
     911,   915,   921,   926,   931,   932,   934,   936,   940,   943,
     946,   949,   952,   955,   960,   962,   964,   970,   974,   976,
     980,   984,   988,   992,   996,  1000,  1002,  1007,  1009,  1013,
    1017,  1022,  1025,  1030,  1033,  1041,  1047,  1055,  1061,  1066,
    1074,  1079,  1087,  1091,  1095,  1099,  1103,  1109,  1115,  1122,
    1127,  1135,  1140,  1143,  1146,  1149,  1152,  1155,  1157,  1161,
    1163,  1167,  1171,  1175,  1179,  1183,  1187,  1191,  1195,  1199,
    1203,  1207,  1211,  1215,  1219,  1223,  1226,  1230,  1234,  1239,
    1244,  1252,  1255,  1259,  1260,  1264,  1266,  1270,  1272,  1276,
    1280,  1282,  1286,  1290,  1294,  1298,  1303,  1305,  1308,  1311,
    1313,  1317,  1322,  1327,  1332,  1337,  1342,  1344,  1346,  1348,
    1351,  1354,  1359,  1364,  1367,  1372,  1377,  1382,  1387,  1392,
    1397,  1402,  1407,  1412,  1417,  1421,  1426,  1435,  1442,  1450,
    1459,  1464
};

/* YYRHS -- A `-1'-separated list of the rules' RHS.  */
static const yytype_int16 yyrhs[] =
{
     124,     0,    -1,   125,    -1,   133,    -1,   136,    -1,   146,
      -1,   147,    -1,   159,    -1,   128,    -1,   129,    -1,   130,
      -1,   131,    -1,   162,    -1,     3,   104,    20,   166,   114,
      -1,     3,   104,   115,   164,   116,    20,   166,   114,    -1,
       3,   104,   115,   164,   116,    20,   126,   114,    -1,     3,
     104,   115,   116,    20,   126,   114,    -1,   127,    -1,   126,
     117,   127,    -1,    86,   118,   166,   117,   177,   119,    -1,
      86,   118,   166,   117,   177,   117,   177,   119,    -1,    88,
     118,   166,   119,    -1,   174,   166,    -1,    10,   103,   118,
     132,   119,    20,   177,   114,    -1,    11,   103,   118,   132,
     119,    20,   177,   114,    -1,    12,   103,   118,   132,   119,
      20,   173,   114,    -1,    13,   103,   118,   132,   119,    20,
     166,   114,    -1,   104,    -1,   132,   117,   104,    -1,     4,
     104,   115,   164,   116,    20,   142,   135,   114,    -1,     4,
     104,   115,   164,   116,    20,   177,   114,    -1,     4,   104,
      20,   134,   114,    -1,     4,    -1,   142,    -1,   177,    -1,
      -1,    31,   177,    -1,     5,   104,   115,   164,   116,   137,
     138,   139,   114,    -1,     5,   104,   137,   138,   139,   114,
      -1,     5,   104,   115,   164,   116,    19,    16,   114,    -1,
       5,   104,    19,    16,   114,    -1,     5,   104,   115,   164,
     116,    16,   140,   141,   114,    -1,     5,   104,    16,   140,
     141,   114,    -1,     5,   104,   115,   164,   116,    17,   138,
     139,   140,   141,   114,    -1,     5,   104,    17,   138,   139,
     140,   141,   114,    -1,    -1,    18,    -1,    19,    17,    -1,
      -1,    33,   177,    -1,    33,   112,    38,    -1,    33,    49,
     173,    50,   177,    51,   112,    38,    52,    -1,    33,    49,
     173,    50,   112,    38,    51,   177,    52,    -1,    -1,    32,
     177,    -1,    32,    38,    -1,    32,    49,   173,    50,   177,
      51,    38,    52,    -1,    32,    49,   173,    50,    38,    51,
     177,    52,    -1,    -1,    29,   177,    -1,    -1,    30,   177,
      -1,   143,    -1,   142,   117,   143,    -1,   170,    -1,   144,
     145,    -1,   174,   177,    -1,    22,   176,    22,    -1,   144,
     176,    22,    -1,   145,   144,   176,    22,    -1,     6,   104,
      21,   154,   114,    -1,     7,   104,    21,   154,   114,    -1,
       8,   104,    21,   148,   114,    -1,   149,    -1,   148,    39,
     149,    -1,    27,   164,    21,   148,    -1,    49,   173,    50,
     148,    52,    -1,    49,   173,    50,   148,    51,   148,    52,
      -1,   154,   153,   154,   151,    -1,   154,   153,   177,   151,
      -1,   177,   153,   154,   151,    -1,   177,   153,   177,   151,
      -1,   177,   153,   154,    32,   177,   151,    -1,   177,   153,
     177,    32,   177,   151,    -1,   177,   153,   154,    33,   177,
     151,    -1,   177,   153,   177,    33,   177,   151,    -1,    89,
     150,    50,   154,   153,   154,    51,   154,   153,   154,    52,
     151,    -1,    89,   150,    50,   177,   153,   154,    51,   154,
     153,   154,    52,   151,    -1,    89,   150,    50,   154,   153,
     177,    51,   154,   153,   154,    52,   151,    -1,    89,   150,
      50,   154,   153,   154,    51,   177,   153,   154,    52,   151,
      -1,    89,   150,    50,   154,   153,   154,    51,   154,   153,
     177,    52,   151,    -1,    89,   150,    50,   177,   153,   177,
      51,   154,   153,   154,    52,   151,    -1,    89,   150,    50,
     177,   153,   154,    51,   177,   153,   154,    52,   151,    -1,
      89,   150,    50,   177,   153,   154,    51,   154,   153,   177,
      52,   151,    -1,    89,   150,    50,   154,   153,   177,    51,
     177,   153,   154,    52,   151,    -1,    89,   150,    50,   154,
     153,   177,    51,   154,   153,   177,    52,   151,    -1,    89,
     150,    50,   154,   153,   154,    51,   177,   153,   177,    52,
     151,    -1,    89,   150,    50,   177,   153,   177,    51,   177,
     153,   154,    52,   151,    -1,    89,   150,    50,   177,   153,
     177,    51,   154,   153,   177,    52,   151,    -1,    89,   150,
      50,   177,   153,   154,    51,   177,   153,   177,    52,   151,
      -1,    89,   150,    50,   154,   153,   177,    51,   177,   153,
     177,    52,   151,    -1,    89,   150,    50,   177,   153,   177,
      51,   177,   153,   177,    52,   151,    -1,    89,   150,    50,
     154,   153,   154,    52,   151,    -1,    89,   150,    50,   177,
     153,   154,    52,   151,    -1,    89,   150,    50,   154,   153,
     177,    52,   151,    -1,    89,   150,    50,   177,   153,   177,
      52,   151,    -1,   154,    37,   154,    -1,   177,    37,   154,
      -1,   154,    37,   177,    -1,   154,    34,   154,    -1,   177,
      34,   154,    -1,   154,    34,   177,    -1,   154,    32,   154,
      -1,   177,    32,   154,    -1,   154,    32,   177,    -1,   154,
      33,   154,    -1,   177,    33,   154,    -1,   154,    33,   177,
      -1,   154,    35,   154,    -1,   177,    35,   154,    -1,   154,
      35,   177,    -1,   154,    36,   154,    -1,   177,    36,   154,
      -1,   154,    36,   177,    -1,   150,    39,   150,    -1,   150,
      40,   150,    -1,   150,    41,   150,    -1,    42,   150,    -1,
     118,   150,   119,    -1,    -1,   151,   117,   152,    -1,   107,
      -1,   108,    -1,   109,    -1,   110,    -1,    32,    -1,    33,
      -1,    34,    -1,   155,    -1,   154,   111,   155,    -1,   154,
     112,   155,    -1,   154,   111,   178,    -1,   154,   112,   178,
      -1,   177,   111,   155,    -1,   177,   112,   155,    -1,   156,
      -1,   155,   113,   179,    -1,   155,   120,   179,    -1,   178,
     113,   156,    -1,   155,   113,   156,    -1,   157,    -1,   111,
     156,    -1,   112,   156,    -1,   158,    -1,   158,    61,   179,
      -1,    43,   164,    21,   155,    -1,    97,   175,    -1,    90,
     118,   154,   119,    -1,    74,   118,   154,   119,    -1,    71,
     118,   154,   119,    -1,    73,   118,   154,   119,    -1,    72,
     118,   154,   119,    -1,    75,   118,   154,   119,    -1,    76,
     118,   154,   119,    -1,    77,   118,   154,   119,    -1,    69,
     118,   154,   119,    -1,    70,   118,   154,   119,    -1,    78,
     118,   154,   117,   177,   119,    -1,    79,   118,   154,   117,
     177,   119,    -1,    49,   173,    50,   154,    51,   154,    52,
      -1,   118,   154,   119,    -1,     9,   104,    21,   160,   114,
      -1,   161,   140,    21,   154,    -1,    27,   164,    21,   160,
      -1,    -1,    91,    -1,    92,    -1,    21,   163,   114,    -1,
      14,   176,    -1,    14,   174,    -1,    14,   166,    -1,    14,
      97,    -1,    15,   173,    -1,    27,   164,    21,   163,    -1,
     165,    -1,   166,    -1,   174,    23,   166,    22,   173,    -1,
     174,    23,   166,    -1,   167,    -1,   166,    54,   167,    -1,
     166,   111,   167,    -1,   166,    56,   167,    -1,   166,   112,
     167,    -1,   166,    57,   167,    -1,   166,    53,   167,    -1,
     168,    -1,    54,   164,    21,   168,    -1,   169,    -1,   168,
      55,   169,    -1,   168,   113,   169,    -1,    53,   164,    21,
     169,    -1,    98,   175,    -1,   102,   118,   176,   119,    -1,
     121,   122,    -1,   121,   177,    24,   177,    26,   177,   122,
      -1,   121,   177,    24,   177,   122,    -1,   121,   177,    25,
     177,    26,   177,   122,    -1,   121,   177,    25,   177,   122,
      -1,    46,   164,    21,   177,    -1,    46,   118,   177,   119,
     164,    21,   177,    -1,    47,   164,    21,   177,    -1,    47,
     118,   177,   119,   164,    21,   177,    -1,   118,   166,   119,
      -1,   121,   172,   122,    -1,   121,   176,   122,    -1,   121,
     164,   122,    -1,   121,   164,    21,   177,   122,    -1,   121,
     164,    21,   174,   122,    -1,    58,   118,   166,   117,   174,
     119,    -1,    87,   118,    98,   119,    -1,    49,   173,    50,
     166,    51,   166,    52,    -1,    80,   177,    81,   177,    -1,
     170,   171,    -1,    82,   177,    -1,    83,   177,    -1,    84,
     177,    -1,    85,   177,    -1,   174,    -1,   172,   117,   174,
      -1,   170,    -1,   177,    34,   177,    -1,   177,    37,   177,
      -1,   177,    36,   177,    -1,   177,    33,   177,    -1,   177,
      35,   177,    -1,   177,    32,   177,    -1,   166,    34,   166,
      -1,   166,    37,   166,    -1,   166,    36,   166,    -1,   166,
      33,   166,    -1,   166,    35,   166,    -1,   166,    32,   166,
      -1,   173,    39,   173,    -1,   173,    40,   173,    -1,   173,
      41,   173,    -1,    42,   173,    -1,   118,   173,   119,    -1,
     174,    23,   166,    -1,    28,   118,   164,   119,    -1,   101,
     118,   176,   119,    -1,    49,   173,    50,   173,    51,   173,
      52,    -1,    35,    36,    -1,    35,   176,    36,    -1,    -1,
     115,   176,   116,    -1,   177,    -1,   176,   117,   177,    -1,
     178,    -1,   177,   111,   178,    -1,   177,   112,   178,    -1,
     179,    -1,   178,   113,   179,    -1,   178,   120,   179,    -1,
     178,    59,   179,    -1,   178,    60,   179,    -1,    48,   164,
      21,   179,    -1,   180,    -1,   111,   180,    -1,   112,   180,
      -1,   181,    -1,   181,    61,   179,    -1,    43,   164,    21,
     178,    -1,    44,   165,    21,   179,    -1,    45,   165,    21,
     179,    -1,    44,   118,   164,   119,    -1,    45,   118,   164,
     119,    -1,   106,    -1,   105,    -1,   104,    -1,    95,   175,
      -1,    96,   175,    -1,    99,   118,   176,   119,    -1,   100,
     118,   176,   119,    -1,   181,    62,    -1,    63,   118,   166,
     119,    -1,    69,   118,   177,   119,    -1,    70,   118,   177,
     119,    -1,    64,   118,   177,   119,    -1,    65,   118,   177,
     119,    -1,    66,   118,   177,   119,    -1,    71,   118,   177,
     119,    -1,    72,   118,   177,   119,    -1,    73,   118,   177,
     119,    -1,    74,   118,   177,   119,    -1,   118,   177,   119,
      -1,    93,   118,   177,   119,    -1,    94,   118,   177,   117,
     177,   117,   177,   119,    -1,    67,   118,   177,   117,   177,
     119,    -1,    49,   173,    50,   177,    51,   177,    52,    -1,
      68,   118,   166,   117,   177,   117,   177,   119,    -1,    44,
     118,   176,   119,    -1,    45,   118,   176,   119,    -1
};

/* YYRLINE[YYN] -- source line where rule number YYN was defined.  */
static const yytype_uint16 yyrline[] =
{
       0,   148,   148,   149,   150,   151,   152,   153,   154,   155,
     156,   157,   158,   166,   173,   179,   185,   195,   196,   199,
     202,   205,   211,   220,   229,   238,   247,   256,   259,   269,
     272,   275,   282,   285,   286,   293,   294,   302,   309,   318,
     328,   339,   348,   358,   362,   372,   373,   374,   378,   381,
     382,   383,   388,   396,   397,   398,   399,   404,   412,   413,
     417,   418,   426,   427,   430,   431,   435,   439,   443,   446,
     458,   461,   471,   477,   480,   483,   488,   493,   501,   504,
     509,   514,   521,   525,   530,   534,   540,   543,   548,   553,
     558,   562,   569,   576,   582,   588,   594,   599,   607,   616,
     625,   633,   644,   647,   651,   656,   664,   665,   668,   671,
     672,   675,   678,   679,   682,   685,   686,   689,   692,   693,
     696,   699,   700,   703,   706,   707,   708,   709,   710,   714,
     715,   719,   720,   721,   722,   726,   727,   728,   732,   733,
     734,   735,   736,   739,   740,   748,   749,   750,   754,   755,
     759,   760,   761,   767,   768,   771,   777,   780,   781,   782,
     783,   784,   785,   786,   787,   788,   789,   790,   793,   796,
     799,   807,   813,   816,   822,   823,   824,   832,   836,   837,
     838,   839,   840,   841,   851,   852,   859,   862,   868,   869,
     870,   873,   874,   877,   878,   881,   882,   886,   887,   888,
     891,   895,   898,   903,   904,   907,   910,   913,   916,   919,
     922,   925,   928,   929,   930,   931,   932,   933,   934,   937,
     940,   946,   947,   951,   952,   953,   954,   958,   961,   964,
     968,   969,   970,   971,   972,   973,   974,   975,   976,   977,
     978,   979,   980,   981,   982,   983,   984,   985,   986,   987,
     992,   998,   999,  1003,  1006,  1012,  1015,  1021,  1022,  1023,
    1027,  1028,  1029,  1030,  1031,  1032,  1038,  1039,  1040,  1044,
    1045,  1046,  1049,  1052,  1055,  1058,  1064,  1065,  1066,  1069,
    1072,  1075,  1080,  1085,  1086,  1087,  1088,  1089,  1090,  1091,
    1092,  1093,  1094,  1095,  1097,  1098,  1099,  1102,  1105,  1108,
    1111,  1114
};
#endif

#if YYDEBUG || YYERROR_VERBOSE || YYTOKEN_TABLE
/* YYTNAME[SYMBOL-NUM] -- String name of the symbol SYMBOL-NUM.
   First, the terminals, then, starting at YYNTOKENS, nonterminals.  */
static const char *const yytname[] =
{
  "$end", "error", "$undefined", "DECLSET", "DECLPAR", "DECLVAR",
  "DECLMIN", "DECLMAX", "DECLSUB", "DECLSOS", "DEFNUMB", "DEFSTRG",
  "DEFBOOL", "DEFSET", "PRINT", "CHECK", "BINARY", "INTEGER", "REAL",
  "IMPLICIT", "ASGN", "DO", "WITH", "IN", "TO", "UNTIL", "BY", "FORALL",
  "EXISTS", "PRIORITY", "STARTVAL", "DEFAULT", "CMP_LE", "CMP_GE",
  "CMP_EQ", "CMP_LT", "CMP_GT", "CMP_NE", "INFTY", "AND", "OR", "XOR",
  "NOT", "SUM", "MIN", "MAX", "ARGMIN", "ARGMAX", "PROD", "IF", "THEN",
  "ELSE", "END", "INTER", "UNION", "CROSS", "SYMDIFF", "WITHOUT", "PROJ",
  "MOD", "DIV", "POW", "FAC", "CARD", "ROUND", "FLOOR", "CEIL", "RANDOM",
  "ORD", "ABS", "SGN", "LOG", "LN", "EXP", "SQRT", "SIN", "COS", "TAN",
  "POWER", "SGNPOW", "READ", "AS", "SKIP", "USE", "COMMENT", "MATCH",
  "SUBSETS", "INDEXSET", "POWERSET", "VIF", "VABS", "TYPE1", "TYPE2",
  "LENGTH", "SUBSTR", "NUMBSYM", "STRGSYM", "VARSYM", "SETSYM", "NUMBDEF",
  "STRGDEF", "BOOLDEF", "SETDEF", "DEFNAME", "NAME", "STRG", "NUMB",
  "SCALE", "SEPARATE", "CHECKONLY", "INDICATOR", "'+'", "'-'", "'*'",
  "';'", "'['", "']'", "','", "'('", "')'", "'/'", "'{'", "'}'", "$accept",
  "stmt", "decl_set", "set_entry_list", "set_entry", "def_numb",
  "def_strg", "def_bool", "def_set", "name_list", "decl_par",
  "par_singleton", "par_default", "decl_var", "var_type", "lower", "upper",
  "priority", "startval", "cexpr_entry_list", "cexpr_entry", "matrix_head",
  "matrix_body", "decl_obj", "decl_sub", "constraint_list", "constraint",
  "vbool", "con_attr_list", "con_attr", "con_type", "vexpr", "vproduct",
  "vfactor", "vexpo", "vval", "decl_sos", "soset", "sos_type", "exec_do",
  "command", "idxset", "pure_idxset", "sexpr", "sunion", "sproduct",
  "sval", "read", "read_par", "tuple_list", "lexpr", "tuple", "symidx",
  "cexpr_list", "cexpr", "cproduct", "cfactor", "cexpo", "cval", 0
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
     315,   316,   317,   318,   319,   320,   321,   322,   323,   324,
     325,   326,   327,   328,   329,   330,   331,   332,   333,   334,
     335,   336,   337,   338,   339,   340,   341,   342,   343,   344,
     345,   346,   347,   348,   349,   350,   351,   352,   353,   354,
     355,   356,   357,   358,   359,   360,   361,   362,   363,   364,
     365,    43,    45,    42,    59,    91,    93,    44,    40,    41,
      47,   123,   125
};
# endif

/* YYR1[YYN] -- Symbol number of symbol that rule YYN derives.  */
static const yytype_uint8 yyr1[] =
{
       0,   123,   124,   124,   124,   124,   124,   124,   124,   124,
     124,   124,   124,   125,   125,   125,   125,   126,   126,   126,
     126,   126,   127,   128,   129,   130,   131,   132,   132,   133,
     133,   133,   133,   134,   134,   135,   135,   136,   136,   136,
     136,   136,   136,   136,   136,   137,   137,   137,   138,   138,
     138,   138,   138,   139,   139,   139,   139,   139,   140,   140,
     141,   141,   142,   142,   142,   142,   143,   144,   145,   145,
     146,   146,   147,   148,   148,   148,   148,   148,   149,   149,
     149,   149,   149,   149,   149,   149,   149,   149,   149,   149,
     149,   149,   149,   149,   149,   149,   149,   149,   149,   149,
     149,   149,   149,   149,   149,   149,   150,   150,   150,   150,
     150,   150,   150,   150,   150,   150,   150,   150,   150,   150,
     150,   150,   150,   150,   150,   150,   150,   150,   150,   151,
     151,   152,   152,   152,   152,   153,   153,   153,   154,   154,
     154,   154,   154,   154,   154,   155,   155,   155,   155,   155,
     156,   156,   156,   157,   157,   157,   158,   158,   158,   158,
     158,   158,   158,   158,   158,   158,   158,   158,   158,   158,
     158,   159,   160,   160,   161,   161,   161,   162,   163,   163,
     163,   163,   163,   163,   164,   164,   165,   165,   166,   166,
     166,   166,   166,   166,   166,   167,   167,   168,   168,   168,
     168,   169,   169,   169,   169,   169,   169,   169,   169,   169,
     169,   169,   169,   169,   169,   169,   169,   169,   169,   169,
     169,   170,   170,   171,   171,   171,   171,   172,   172,   172,
     173,   173,   173,   173,   173,   173,   173,   173,   173,   173,
     173,   173,   173,   173,   173,   173,   173,   173,   173,   173,
     173,   174,   174,   175,   175,   176,   176,   177,   177,   177,
     178,   178,   178,   178,   178,   178,   179,   179,   179,   180,
     180,   180,   180,   180,   180,   180,   181,   181,   181,   181,
     181,   181,   181,   181,   181,   181,   181,   181,   181,   181,
     181,   181,   181,   181,   181,   181,   181,   181,   181,   181,
     181,   181
};

/* YYR2[YYN] -- Number of symbols composing right hand side of rule YYN.  */
static const yytype_uint8 yyr2[] =
{
       0,     2,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     5,     8,     8,     7,     1,     3,     6,
       8,     4,     2,     8,     8,     8,     8,     1,     3,     9,
       8,     5,     1,     1,     1,     0,     2,     9,     6,     8,
       5,     9,     6,    11,     8,     0,     1,     2,     0,     2,
       3,     9,     9,     0,     2,     2,     8,     8,     0,     2,
       0,     2,     1,     3,     1,     2,     2,     3,     3,     4,
       5,     5,     5,     1,     3,     4,     5,     7,     4,     4,
       4,     4,     6,     6,     6,     6,    12,    12,    12,    12,
      12,    12,    12,    12,    12,    12,    12,    12,    12,    12,
      12,    12,     8,     8,     8,     8,     3,     3,     3,     3,
       3,     3,     3,     3,     3,     3,     3,     3,     3,     3,
       3,     3,     3,     3,     3,     3,     3,     2,     3,     0,
       3,     1,     1,     1,     1,     1,     1,     1,     1,     3,
       3,     3,     3,     3,     3,     1,     3,     3,     3,     3,
       1,     2,     2,     1,     3,     4,     2,     4,     4,     4,
       4,     4,     4,     4,     4,     4,     4,     6,     6,     7,
       3,     5,     4,     4,     0,     1,     1,     3,     2,     2,
       2,     2,     2,     4,     1,     1,     5,     3,     1,     3,
       3,     3,     3,     3,     3,     1,     4,     1,     3,     3,
       4,     2,     4,     2,     7,     5,     7,     5,     4,     7,
       4,     7,     3,     3,     3,     3,     5,     5,     6,     4,
       7,     4,     2,     2,     2,     2,     2,     1,     3,     1,
       3,     3,     3,     3,     3,     3,     3,     3,     3,     3,
       3,     3,     3,     3,     3,     2,     3,     3,     4,     4,
       7,     2,     3,     0,     3,     1,     3,     1,     3,     3,
       1,     3,     3,     3,     3,     4,     1,     2,     2,     1,
       3,     4,     4,     4,     4,     4,     1,     1,     1,     2,
       2,     4,     4,     2,     4,     4,     4,     4,     4,     4,
       4,     4,     4,     4,     3,     4,     8,     6,     7,     8,
       4,     4
};

/* YYDEFACT[STATE-NAME] -- Default rule to reduce with in state
   STATE-NUM when YYTABLE doesn't specify something else to do.  Zero
   means the default is an error.  */
static const yytype_uint16 yydefact[] =
{
       0,     0,    32,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     2,     8,     9,    10,    11,     3,
       4,     5,     6,     7,    12,     0,     0,    45,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       1,     0,     0,     0,     0,    58,    48,    46,     0,     0,
      48,     0,     0,     0,   174,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,   253,   253,   181,   253,     0,
       0,     0,   278,   277,   276,     0,     0,     0,     0,   180,
     188,   195,   197,   179,   178,   255,   257,   260,   266,   269,
       0,     0,     0,     0,     0,     0,   182,     0,     0,     0,
       0,     0,   184,   185,     0,   177,     0,     0,     0,     0,
       0,     0,     0,     0,    33,    62,     0,    64,     0,    34,
       0,     0,    60,     0,    53,     0,    47,     0,    53,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,   253,     0,     0,     0,     0,   138,   145,
     150,   153,     0,   257,     0,     0,     0,     0,     0,    73,
       0,     0,     0,   175,   176,     0,    58,    27,     0,     0,
       0,     0,   251,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,   279,   280,   201,     0,     0,     0,
     267,   268,     0,     0,   203,     0,   229,     0,   227,     0,
     255,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,   283,     0,   245,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,    13,     0,     0,     0,     0,     0,
      31,     0,     0,    65,     0,     0,     0,     0,   222,    66,
       0,    59,     0,     0,     0,     0,    49,     0,    58,    40,
      45,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,   156,     0,     0,   151,
     152,     0,     0,     0,     0,    70,     0,     0,     0,     0,
       0,     0,    71,     0,     0,     0,     0,     0,     0,     0,
       0,    72,   135,   136,   137,     0,     0,     0,   171,     0,
       0,     0,     0,     0,     0,   252,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,   212,   294,     0,   215,     0,   213,   214,     0,     0,
     194,   189,   191,   193,   190,   192,   198,   199,   256,   258,
     259,   263,   264,   261,   262,   270,     0,     0,     0,   246,
     241,   239,   236,   240,   238,   237,   242,   243,   244,   247,
     235,   233,   230,   234,   232,   231,     0,   183,   187,     0,
       0,     0,    17,     0,     0,    67,     0,     0,    63,     0,
       0,   223,   224,   225,   226,     0,    61,    42,     0,    50,
      55,     0,    54,    60,    58,    48,     0,    48,    38,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,   170,   139,
     141,   140,   142,   149,   146,   147,   154,   143,   258,   144,
     259,   148,     0,     0,   127,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,    74,   129,   129,   129,   129,   174,
       0,    28,     0,     0,     0,     0,   271,   274,   300,   272,
     275,   301,   273,     0,   208,     0,   210,   265,     0,     0,
     200,   196,     0,   284,   287,   288,   289,     0,     0,   285,
     286,   290,   291,   292,   293,   219,   295,     0,   254,   281,
     282,   202,     0,     0,   228,     0,     0,   248,     0,     0,
       0,   249,     0,     0,     0,    16,     0,    22,     0,     0,
     221,    68,     0,    35,     0,     0,     0,     0,    60,    53,
       0,    53,   155,   271,     0,     0,   165,   166,   159,   161,
     160,   158,   162,   163,   164,     0,     0,   157,     0,     0,
       0,     0,     0,     0,     0,     0,    75,     0,     0,     0,
     128,   124,   125,   126,     0,     0,   112,   114,   115,   117,
     109,   111,   118,   120,   121,   123,   106,   108,   113,   116,
     110,   119,   122,   107,    78,    79,     0,     0,    80,     0,
       0,    81,   173,   172,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,   217,   216,     0,   205,
       0,   207,     0,   186,     0,     0,    18,    15,    14,    69,
       0,     0,    30,     0,     0,     0,    44,     0,    58,    39,
       0,     0,     0,     0,     0,     0,     0,    76,     0,     0,
       0,   129,   129,   129,   129,    23,    24,    25,    26,     0,
       0,     0,     0,   218,   297,     0,     0,     0,     0,     0,
       0,    21,    36,    29,     0,     0,     0,     0,    41,    60,
      37,     0,   167,   168,     0,     0,     0,     0,     0,     0,
     131,   132,   133,   134,   130,    82,    84,    83,    85,   209,
     211,   220,   298,     0,     0,   204,   206,   250,     0,     0,
       0,     0,     0,     0,   169,    77,     0,   129,     0,   129,
       0,   129,     0,   129,   299,   296,     0,    19,     0,     0,
       0,     0,    43,     0,     0,   102,     0,     0,   104,     0,
       0,   103,     0,     0,   105,     0,    52,    51,    57,    56,
       0,     0,     0,     0,     0,     0,     0,     0,    20,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,   129,   129,   129,   129,   129,
     129,   129,   129,   129,   129,   129,   129,   129,   129,   129,
     129,    86,    90,    89,    96,    88,    95,    94,   100,    87,
      93,    92,    99,    91,    98,    97,   101
};

/* YYDEFGOTO[NTERM-NUM].  */
static const yytype_int16 yydefgoto[] =
{
      -1,    13,    14,   451,   452,    15,    16,    17,    18,   188,
      19,   133,   711,    20,    50,   144,   308,   142,   303,   134,
     135,   136,   293,    21,    22,   178,   179,   347,   674,   774,
     355,   348,   168,   169,   170,   171,    23,   185,   186,    24,
      39,   121,   122,   123,   100,   101,   102,   137,   298,   237,
     116,   124,   224,   104,   172,   106,   107,   108,   109
};

/* YYPACT[STATE-NUM] -- Index in YYTABLE of the portion describing
   STATE-NUM.  */
#define YYPACT_NINF -513
static const yytype_int16 yypact[] =
{
    1259,   -55,   -20,     2,    31,    64,    73,    80,    -3,    16,
      36,    93,   422,   206,  -513,  -513,  -513,  -513,  -513,  -513,
    -513,  -513,  -513,  -513,  -513,    -4,    22,    21,   217,   232,
     236,   315,   221,   235,   238,   256,  1779,  1512,   500,   300,
    -513,   959,   358,  1331,   500,   388,   392,  -513,   580,   500,
     392,  2520,  2520,  1621,   171,   330,   330,   330,   330,  1426,
     500,     6,    11,   524,   578,   500,  1512,   500,   500,   345,
     351,   412,   434,   447,   467,   485,   499,   501,   505,   512,
     541,   545,   557,   573,   603,   495,   495,  -513,   495,   613,
     617,   627,  -513,  -513,  -513,  3012,  3012,  2325,  1699,  1149,
    -513,     9,  -513,  -513,   550,   657,   210,  -513,  -513,   719,
     630,  1512,  1512,   643,  1512,   787,   661,   655,   311,  1512,
     959,   674,  -513,  1149,   701,  -513,   837,   738,   623,  2701,
    1512,  2701,  2701,   704,   716,  -513,   817,  1213,  2701,   657,
     730,  2701,   822,  2766,   824,   768,  -513,   769,   824,   500,
    1512,   774,   860,   880,   883,   886,   889,   898,   901,   931,
     942,   944,   948,   495,  2578,  2578,  2520,   370,   -23,  -513,
    -513,   849,   691,   295,   396,   500,  1512,  2390,    34,  -513,
     190,   258,   500,  -513,  -513,   823,   388,  -513,   148,   321,
     419,   431,  -513,    17,   929,  1858,   955,  1858,  1048,  2325,
    1055,  2325,  1086,  1088,   732,  1103,  1105,   959,   959,  2701,
    2701,  2701,  2701,   959,  2701,  2701,  2701,  2701,  2701,  2701,
     854,  2701,  2701,  2701,  -513,  -513,  -513,  2701,  2701,  2701,
    -513,  -513,   307,    96,  -513,    -7,  1213,   109,   701,   110,
      20,   959,   959,   959,   959,   959,   959,  1069,  1069,  2701,
    2701,  2701,  2954,  2954,  2954,  2954,  2954,  -513,   500,  -513,
     795,  2701,   343,    19,   184,   959,   959,   959,   959,   959,
     959,  1512,  1512,  1512,   959,  2701,  2701,  2701,  2701,  2701,
    2701,   810,   422,   959,  -513,   249,  1081,    25,   927,   -30,
    -513,  1109,  2701,   817,  2701,  2701,  2701,  2701,  -513,   657,
    1127,   657,  2701,  1041,  1512,  2131,   657,   891,   388,  -513,
    1365,  1049,  1144,  1000,  2520,  2520,  2520,  2520,  2520,  2520,
    2520,  2520,  2520,  2520,  2520,  2520,  -513,  3024,  3024,  -513,
    -513,   188,   400,  2520,  2520,  -513,  2636,  2954,  2954,  2520,
    2520,  2636,  -513,  1145,  1033,  2390,  2390,  1080,   751,   832,
    2455,  -513,  -513,  -513,  -513,  2520,  2520,  1154,  -513,  1160,
    1066,  1169,  1171,  1172,  1174,  -513,  2701,  1076,   551,  2954,
    1090,   567,  2954,   640,  2701,   777,  2701,  2954,  2325,  1069,
     636,   830,   394,   789,   819,   834,   254,   926,   862,   881,
     899,   993,  1050,  1096,  1095,  1100,   361,   816,   648,   711,
     778,  -513,  -513,  1937,  -513,  1109,  -513,  -513,  2701,  2701,
    -513,  -513,  -513,  -513,  -513,  -513,  -513,  -513,   657,   210,
     210,  -513,  -513,  -513,  -513,  -513,  1106,  1512,   788,  -513,
    1149,  1149,  1149,  1149,  1149,  1149,  -513,  1187,  1187,  1149,
     657,   657,   657,   657,   657,   657,   959,  -513,    35,  1121,
    1122,   -51,  -513,   959,   474,  -513,  2701,  2701,  -513,    26,
    2701,   657,   657,   657,   657,  1331,   657,  -513,  1102,  -513,
    -513,  1512,   657,   822,   388,   392,  1112,   392,  -513,  2520,
    2520,  1111,  1117,  1132,  1134,  1137,  1164,  1167,  1199,  1201,
    1220,  1222,  1224,  1238,  1249,  1251,   469,   792,  1295,   500,
    1512,  1136,  1139,  1140,  1141,  1155,  1156,  2520,  -513,   -23,
     295,   -23,   295,  -513,  -513,  -513,  -513,   -23,   295,   -23,
     295,  -513,  1621,  1621,  -513,   115,   214,   240,  2390,  2390,
    2390,  2520,  2520,  2520,  2520,  2520,  2520,  2520,  2520,  2520,
    2520,  2520,  2520,  2520,  -513,   974,   691,   118,   272,   171,
    2520,  -513,  2701,  2701,  1512,   959,   210,  -513,  -513,  -513,
    -513,  -513,  -513,   500,   657,   500,   657,  -513,  1002,    -8,
    -513,     9,  1109,  -513,  -513,  -513,  -513,  2701,  2701,  -513,
    -513,  -513,  -513,  -513,  -513,  -513,  -513,  2701,  -513,  -513,
    -513,  -513,  1130,   -43,  -513,    -2,     1,  -513,   693,   598,
     742,  -513,  1512,   959,   959,  -513,  1109,  1149,   -40,  1011,
     657,  -513,    28,    -9,   606,  2831,  1110,  1113,   822,   824,
    1173,   824,   -23,   295,   143,   290,  -513,  -513,  -513,  -513,
    -513,  -513,  -513,  -513,  -513,  2701,  2701,  -513,  1234,  1118,
    2520,  2520,  2520,  2520,  2520,  2520,  1250,    47,   177,   283,
    -513,  -513,  1252,  1252,   190,   258,   974,   691,   974,   691,
     974,   691,   974,   691,   974,   691,   974,   691,   974,   974,
     974,   974,   974,   974,  1176,  1176,  2701,  2701,  1176,  2701,
    2701,  1176,  -513,   974,   940,   970,   131,  1120,  1288,  1293,
     959,  2701,  1196,  1297,   858,   877,  -513,  -513,  2701,  -513,
    2701,  -513,  1512,   661,   991,   414,  -513,  -513,  -513,  -513,
    2701,  1202,  -513,  2195,   341,  2002,  -513,  1205,   388,  -513,
    1207,  2520,  1309,  1321,  2520,  2520,  1621,  -513,  2520,  2520,
    1349,   657,   657,   657,   657,  -513,  -513,  -513,  -513,  2701,
    2701,  1126,   200,  -513,  -513,  2701,  2701,   445,   576,    15,
    2701,  -513,   657,  -513,  1279,  2896,  1286,   379,  -513,   822,
    -513,   316,  -513,  -513,   295,    13,   246,   270,   298,   391,
    -513,  -513,  -513,  -513,  -513,  1176,  1176,  1176,  1176,   657,
     657,  -513,  -513,  1327,  1333,  -513,  -513,  -513,   596,  2701,
    2259,  2701,  2067,  1226,  -513,  -513,  2520,  -513,  2520,  -513,
    2520,  -513,  2520,  -513,  -513,  -513,  2701,  -513,   346,  1290,
     354,  1292,  -513,   190,   258,  1176,   190,   258,  1176,   190,
     258,  1176,   190,   258,  1176,  1336,  -513,  -513,  -513,  -513,
    2520,  2520,  2520,  2520,  2520,  2520,  2520,  2520,  -513,   428,
     452,   463,   472,   477,   479,   482,   503,   517,   549,   560,
     594,   600,   621,   625,   629,  -513,  -513,  -513,  -513,  -513,
    -513,  -513,  -513,  -513,  -513,  -513,  -513,  -513,  -513,  -513,
    -513,  1176,  1176,  1176,  1176,  1176,  1176,  1176,  1176,  1176,
    1176,  1176,  1176,  1176,  1176,  1176,  1176
};

/* YYPGOTO[NTERM-NUM].  */
static const yytype_int16 yypgoto[] =
{
    -513,  -513,  -513,   894,   745,  -513,  -513,  -513,  -513,   577,
    -513,  -513,  -513,  -513,  1045,   -44,  -144,  -181,  -465,   902,
    1073,  -127,  -513,  -513,  -513,  -512,  1019,  -325,    58,  -513,
    -166,   556,  -321,  -139,  -513,  -513,  -513,   828,  -513,  -513,
    1089,   264,  1077,  1034,   144,   992,  -219,  1275,  -513,  -513,
      14,  1074,    -1,   -28,   -36,   492,  -220,   -93,  -513
};

/* YYTABLE[YYPACT[STATE-NUM]].  What to do in state STATE-NUM.  If
   positive, shift that token.  If negative, reduce the rule which
   number is the opposite.  If zero, do what YYDEFACT says.
   If YYTABLE_NINF, syntax error.  */
#define YYTABLE_NINF -1
static const yytype_uint16 yytable[] =
{
     105,   118,   230,   231,   311,   359,   148,   139,   617,   292,
     646,   647,   509,   511,   403,   356,    41,   181,   517,   519,
     524,   525,   710,   105,   698,   329,   330,   700,   416,   417,
     118,   193,   421,   422,   423,   424,   425,    45,    46,    47,
      48,    59,    43,   691,   408,   409,    59,   455,   611,    25,
     709,   457,   350,   365,   271,   272,   273,   602,   271,   272,
     273,   233,   240,   605,   247,   795,   606,   787,   250,   251,
     239,   230,   231,   350,   707,   118,   118,   606,   264,   697,
     204,   250,   251,   118,    26,   225,   350,   226,   241,   242,
     336,   243,   244,   105,   118,   289,   233,   337,   726,   727,
      32,   287,   299,   250,   251,   301,    27,   306,   291,   250,
     251,    42,   250,   251,   118,   404,   514,   515,   516,    33,
     699,   423,   248,   701,   195,   259,   260,   473,   263,   197,
     332,   250,   251,   281,   249,    28,    49,    44,   429,    34,
     118,   349,   249,   249,   288,   249,   245,   246,   351,   559,
     676,   677,   562,   717,   528,   529,   530,   567,   622,   105,
     570,   105,   326,   373,   313,   375,   460,   368,    29,   371,
     271,   272,   273,   383,   384,   385,   386,    30,   388,   389,
     390,   391,   392,   393,    31,   395,   396,   105,   329,   330,
     344,   105,   105,   105,   721,   397,    35,   513,   182,   398,
     399,   400,   521,   651,   652,   653,    40,   250,   251,   352,
     353,   354,   231,   418,   765,   402,   275,   276,   277,   278,
     279,   280,   352,   353,   354,   105,   405,   249,   721,   333,
     334,   406,   407,   428,   650,   118,   118,   118,    51,   440,
     441,   442,   443,   444,   445,   737,   532,   533,   534,   535,
     536,   537,   782,    52,   333,   334,   105,    53,   461,   462,
     463,   464,   183,   184,   459,   360,   466,   361,   118,   252,
     253,   472,   538,   539,   540,   541,   542,   543,   482,   484,
     486,   488,   490,   492,    59,   436,   437,   438,   333,   334,
     352,   353,   354,   618,   793,   250,   251,   796,   797,   333,
     334,   333,   334,   402,   679,   680,   128,   508,   140,   349,
     527,   250,   251,   147,   181,   352,   353,   354,   468,   546,
     548,   798,   799,   254,   194,   333,   334,   200,   202,   203,
     255,   205,   206,   508,   691,   449,    54,   450,   564,    55,
     566,   691,   569,   275,   276,   277,   278,   279,   280,   800,
     801,   339,   340,    56,   252,   253,    57,   333,   334,   402,
     241,   242,   235,   243,   244,   250,   251,   593,   794,   339,
     340,   577,   595,   596,    58,   265,   266,   267,   268,   269,
     270,   339,   340,   339,   340,   410,   411,   412,   413,   414,
     415,   600,   755,    59,   339,   340,   241,   242,   826,   243,
     244,   339,   340,   622,    63,    64,   828,   119,   341,   333,
     334,    67,    68,   312,   125,   255,    69,   141,   245,   246,
     569,   610,   250,   251,   105,   143,   401,   333,   334,   614,
     792,   619,   612,   621,   187,   118,    36,    37,   360,   343,
     362,   599,   802,   803,   625,    82,   357,   241,   242,    38,
     243,   244,   250,   251,   245,   246,    88,   250,   251,   367,
      91,   370,   401,   207,   118,   250,   251,   241,   242,   208,
     243,   244,   250,   251,   127,   718,   120,   720,   587,    98,
     855,   333,   334,   356,   335,   616,   181,   649,   728,   729,
     250,   251,   349,   349,   349,   655,   657,   659,   661,   663,
     665,   667,   339,   340,   856,   245,   246,   333,   334,    59,
     342,   339,   340,   573,   639,   857,   684,   685,   118,   402,
      63,    64,   426,   119,   858,   245,   246,    67,    68,   859,
     209,   860,    69,   751,   861,    59,   360,   759,   363,   333,
     334,   693,   694,   173,   173,   173,    63,    64,   360,   119,
     364,   695,   210,    67,    68,   862,   250,   251,    69,    59,
     449,    82,   450,   339,   340,   211,   118,   785,   686,   863,
      63,    64,    88,   119,   333,   334,    91,    67,    68,   714,
     333,   334,    69,   339,   340,   212,   635,    82,   333,   334,
     339,   340,   120,   333,   334,    98,   145,   146,    88,   722,
     723,   864,    91,   213,   675,   678,   681,   167,   174,   180,
     223,    82,   865,    59,   339,   340,   703,   214,   120,   215,
     231,    98,    88,   216,    63,    64,    91,   119,   333,   334,
     217,    67,    68,   189,   190,   191,    69,   271,   272,   273,
     731,   732,   199,   733,   734,    98,   866,   830,   831,   702,
     832,   833,   867,   834,   835,   742,   836,   837,   173,   218,
     339,   340,   747,   219,   748,    82,   118,   249,   249,   173,
     558,   333,   334,   868,   752,   220,    88,   869,   274,   757,
      91,   870,    63,    64,   249,   119,   561,   250,   251,    67,
     181,   221,   767,   769,    69,   282,   201,   231,   786,    98,
     271,   272,   273,   779,   780,   339,   340,   250,   251,   783,
     784,   333,   334,   806,   788,   807,   749,   250,   251,   742,
     712,   222,   331,    82,   283,   265,   266,   267,   268,   269,
     270,   227,   339,   340,    88,   228,   333,   334,    91,   286,
     339,   340,   419,   420,   690,   229,   241,   242,   258,   243,
     244,   250,   251,   808,   120,   810,   742,    98,   285,   563,
     814,   261,   817,   638,   820,   249,   823,   589,   250,   251,
     825,   271,   272,   273,   275,   276,   277,   278,   279,   280,
     256,   257,   378,   532,   533,   534,   535,   536,   537,   775,
     776,   777,   778,   691,   840,   842,   844,   846,   848,   850,
     852,   854,   339,   340,   245,   246,   173,   173,   173,   173,
     173,   173,   173,   173,   173,   173,   173,   173,   290,   265,
     266,   267,   268,   269,   270,   510,   512,   688,   249,   689,
     590,   518,   520,   291,   271,   272,   273,   173,   173,   129,
     241,   242,   173,   243,   244,   427,   300,   173,   173,   271,
     272,   273,   302,   250,   251,   815,   307,   818,   556,   821,
     446,   824,   333,   334,   538,   539,   540,   541,   542,   543,
     481,   483,   485,   487,   489,   491,   493,   494,   495,   496,
     497,   498,   309,   241,   242,   310,   243,   244,   250,   251,
     241,   242,   314,   243,   244,   249,   565,   591,   245,   246,
     250,   251,   526,   333,   334,   249,   180,   601,   574,   636,
     338,   545,   547,   871,   872,   873,   874,   875,   876,   877,
     878,   879,   880,   881,   882,   883,   884,   885,   886,   470,
     250,   251,   588,   249,    60,    61,    62,   358,   575,    65,
     471,   245,   246,   339,   340,   250,   251,   572,   245,   246,
     366,   284,   394,   576,    70,    71,    72,    73,    74,    75,
      76,    77,    78,    79,    80,    81,   271,   272,   273,   250,
     251,   623,   173,   250,   251,   745,   369,   456,   315,   241,
     242,   579,   243,   244,    83,    84,    85,    86,   250,   251,
      89,    90,   250,   251,   746,    92,    93,    94,   316,   173,
     580,   317,    95,    96,   318,    63,    64,   319,   119,   132,
     250,   251,    67,    68,   173,   173,   320,    69,   581,   321,
     173,   173,   173,   173,   173,   173,   173,   173,   173,   173,
     173,   173,   173,   173,   173,   173,   624,   245,   246,   271,
     272,   273,   173,   578,   241,   242,    82,   243,   244,   322,
     480,   250,   251,   690,   735,   241,   242,    88,   243,   244,
     323,    91,   324,   331,   241,   242,   325,   243,   244,   372,
      99,   115,   271,   272,   273,   126,   374,   120,   180,   648,
      98,   250,   251,   523,   736,   333,   334,   654,   656,   658,
     660,   662,   664,   666,   668,   669,   670,   671,   672,   673,
     115,   454,   245,   246,   250,   251,   683,   376,   750,   377,
     103,   117,   582,   245,   246,    63,    64,   138,   119,   528,
     529,   530,   245,   246,   379,   708,   380,    69,   620,   146,
     531,   232,   173,   173,   173,   173,   173,   173,   196,   198,
     117,   271,   272,   273,    59,   115,   115,   465,   262,   271,
     272,   273,   615,   115,   232,   467,    82,   271,   272,   273,
     715,   250,   251,   478,   115,   479,   522,    88,   725,   583,
     551,    91,   238,   241,   242,   549,   243,   244,   781,   241,
     242,   550,   243,   244,   115,   117,   117,   120,   117,   552,
      98,   553,   554,   117,   555,   557,   481,   483,   485,   487,
     489,   491,   241,   242,   117,   243,   244,   250,   251,   560,
     115,   250,   251,   173,   585,   584,   764,   173,   173,   586,
     173,   173,   333,   334,   117,   597,   271,   716,   339,   340,
     626,   245,   246,   232,   738,   232,   579,   245,   246,   603,
     604,   381,   382,   333,   334,   339,   340,   387,   333,   334,
     117,   627,   696,   580,   640,   724,   628,   641,   642,   643,
     245,   246,     1,     2,     3,     4,     5,     6,     7,     8,
       9,    10,    11,   644,   645,   339,   340,   761,   333,   334,
      12,   624,   180,   581,   766,   768,   629,   719,   173,   350,
     173,   528,   173,   730,   173,   294,   295,   296,   297,   430,
     431,   432,   433,   434,   435,   115,   115,   115,   439,   739,
     339,   340,   333,   334,   740,   743,   753,   448,   582,   758,
     630,   760,   173,   173,   173,   173,   173,   173,   173,   173,
     789,   339,   340,   333,   334,   339,   340,   791,   115,   583,
     812,   631,   827,   584,   829,   117,   117,   117,   608,   333,
     334,   706,   813,   129,   816,   477,   819,   632,   822,   453,
     333,   334,   333,   334,   458,   138,    59,   613,   633,   544,
     634,   447,   571,   236,    60,    61,    62,   682,   117,    65,
     130,   474,   475,    47,   476,     0,   839,   841,   843,   845,
     847,   849,   851,   853,    70,    71,    72,    73,    74,    75,
      76,    77,    78,    79,    80,    81,   333,   334,   250,   251,
       0,   131,   568,     0,   637,     0,   744,     0,     0,     0,
     250,   251,     0,     0,    83,    84,    85,    86,   762,     0,
      89,    90,   250,   251,     0,    92,    93,    94,   250,   251,
     763,     0,    95,    96,   250,   251,   804,   250,   251,   132,
       0,     0,   805,     0,     0,   838,   770,   771,   772,   773,
       0,   598,   192,     0,     0,     0,     0,     0,     0,    60,
      61,    62,     0,     0,    65,   130,     0,   592,     0,   594,
     568,     0,     0,     0,     0,     0,     0,   607,   609,    70,
      71,    72,    73,    74,    75,    76,    77,    78,    79,    80,
      81,   117,     0,     0,     0,   115,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,    83,
      84,    85,    86,     0,     0,    89,    90,     0,   453,     0,
      92,    93,    94,     0,   115,     0,     0,    95,    96,   138,
     110,     0,     0,     0,   132,   117,     0,    59,     0,     0,
       0,     0,     0,     0,   111,    60,    61,    62,    63,    64,
      65,   112,     0,     0,     0,    67,    68,     0,     0,     0,
      69,     0,     0,     0,   117,    70,    71,    72,    73,    74,
      75,    76,    77,    78,    79,    80,    81,     0,   115,   687,
       0,     0,     0,     0,     0,     0,     0,     0,     0,    82,
       0,     0,     0,     0,     0,    83,    84,    85,    86,     0,
      88,    89,    90,   113,    91,     0,    92,    93,    94,     0,
       0,     0,     0,    95,    96,     0,     0,     0,   117,     0,
     114,     0,     0,    98,     0,     0,   115,   704,   705,     0,
       0,     0,     0,     0,     0,     0,   692,     0,   175,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,   149,    61,    62,     0,     0,    65,
     176,     0,     0,     0,     0,     0,   117,     0,     0,     0,
     453,     0,     0,     0,    70,    71,    72,    73,    74,    75,
     151,   152,   153,   154,   155,   156,   157,   158,   159,   160,
     161,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     177,   162,     0,     0,    83,    84,    85,    86,   163,     0,
      89,    90,     0,     0,   741,    92,    93,    94,     0,     0,
       0,     0,   164,   165,    59,     0,   115,     0,     0,   166,
       0,     0,    60,    61,    62,    63,    64,    65,    66,     0,
       0,     0,    67,    68,     0,     0,     0,    69,     0,     0,
       0,     0,    70,    71,    72,    73,    74,    75,    76,    77,
      78,    79,    80,    81,     0,     0,   117,     0,     0,   131,
       0,     0,     0,     0,     0,     0,    82,     0,     0,     0,
       0,     0,    83,    84,    85,    86,     0,    88,    89,    90,
       0,    91,     0,    92,    93,    94,     0,     0,     0,     0,
      95,    96,     0,     0,    59,     0,     0,    97,     0,     0,
      98,   234,    60,    61,    62,    63,    64,    65,    66,     0,
       0,     0,    67,    68,     0,     0,     0,    69,     0,     0,
       0,     0,    70,    71,    72,    73,    74,    75,    76,    77,
      78,    79,    80,    81,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,    82,     0,     0,     0,
       0,     0,    83,    84,    85,    86,    87,    88,    89,    90,
       0,    91,     0,    92,    93,    94,     0,     0,     0,     0,
      95,    96,     0,    59,     0,     0,     0,    97,     0,     0,
      98,    60,    61,    62,    63,    64,    65,    66,     0,     0,
       0,    67,    68,     0,     0,     0,    69,     0,     0,     0,
       0,    70,    71,    72,    73,    74,    75,    76,    77,    78,
      79,    80,    81,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,    82,     0,     0,     0,     0,
       0,    83,    84,    85,    86,     0,    88,    89,    90,     0,
      91,     0,    92,    93,    94,     0,     0,     0,     0,    95,
      96,     0,    59,     0,     0,     0,    97,     0,     0,    98,
      60,    61,    62,     0,     0,    65,   130,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
      70,    71,    72,    73,    74,    75,    76,    77,    78,    79,
      80,    81,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
      83,    84,    85,    86,     0,     0,    89,    90,     0,     0,
     756,    92,    93,    94,     0,    60,    61,    62,    95,    96,
      65,   130,     0,     0,     0,   132,     0,     0,     0,     0,
       0,     0,     0,     0,     0,    70,    71,    72,    73,    74,
      75,    76,    77,    78,    79,    80,    81,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,    83,    84,    85,    86,     0,
       0,    89,    90,     0,     0,   811,    92,    93,    94,     0,
      60,    61,    62,    95,    96,    65,   130,     0,     0,     0,
     132,     0,     0,     0,     0,     0,     0,     0,     0,     0,
      70,    71,    72,    73,    74,    75,    76,    77,    78,    79,
      80,    81,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
      83,    84,    85,    86,     0,     0,    89,    90,     0,   469,
       0,    92,    93,    94,    60,    61,    62,     0,    95,    96,
     130,     0,     0,     0,     0,   132,     0,     0,     0,     0,
       0,     0,     0,     0,    70,    71,    72,    73,    74,    75,
      76,    77,    78,    79,    80,    81,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,    83,    84,    85,    86,     0,     0,
      89,    90,     0,   754,     0,    92,    93,    94,    60,    61,
      62,     0,     0,     0,   130,     0,     0,     0,     0,   132,
       0,     0,     0,     0,     0,     0,     0,     0,    70,    71,
      72,    73,    74,    75,    76,    77,    78,    79,    80,    81,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,    83,    84,
      85,    86,     0,     0,    89,    90,     0,   809,     0,    92,
      93,    94,    60,    61,    62,     0,     0,     0,   130,     0,
       0,     0,     0,   132,     0,     0,     0,     0,     0,     0,
       0,     0,    70,    71,    72,    73,    74,    75,    76,    77,
      78,    79,    80,    81,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,    83,    84,    85,    86,     0,     0,    89,    90,
       0,     0,     0,    92,    93,    94,     0,     0,    60,    61,
      62,    63,    64,    65,    66,     0,     0,   132,    67,    68,
       0,     0,     0,    69,     0,     0,     0,     0,    70,    71,
      72,    73,    74,    75,    76,    77,    78,    79,    80,    81,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,    82,     0,     0,     0,     0,     0,    83,    84,
      85,    86,     0,    88,    89,    90,     0,    91,     0,    92,
      93,    94,   345,   149,    61,    62,    95,    96,    65,   150,
       0,     0,     0,    97,     0,     0,    98,     0,     0,     0,
       0,     0,     0,    70,    71,    72,    73,    74,    75,   151,
     152,   153,   154,   155,   156,   157,   158,   159,   160,   161,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     162,     0,     0,    83,    84,    85,    86,   163,     0,    89,
      90,     0,     0,     0,    92,    93,    94,     0,   149,    61,
      62,   164,   165,    65,   150,     0,     0,     0,   346,     0,
       0,     0,     0,     0,     0,     0,     0,     0,    70,    71,
      72,    73,    74,    75,   151,   152,   153,   154,   155,   156,
     157,   158,   159,   160,   161,     0,     0,     0,     0,     0,
       0,     0,     0,     0,   177,   162,     0,     0,    83,    84,
      85,    86,   163,     0,    89,    90,     0,     0,     0,    92,
      93,    94,     0,   149,    61,    62,   164,   165,    65,   150,
       0,     0,     0,   166,     0,     0,     0,     0,     0,     0,
       0,     0,     0,    70,    71,    72,    73,    74,    75,   151,
     152,   153,   154,   155,   156,   157,   158,   159,   160,   161,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     162,     0,     0,    83,    84,    85,    86,   163,     0,    89,
      90,   149,    61,    62,    92,    93,    94,   150,     0,     0,
       0,   164,   165,     0,     0,     0,     0,     0,   166,     0,
       0,    70,    71,    72,    73,    74,    75,   151,   152,   153,
     154,   155,   156,   157,   158,   159,   160,   161,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,   162,     0,
       0,    83,    84,    85,    86,   163,     0,    89,    90,   149,
      61,    62,    92,    93,    94,   150,     0,     0,     0,   327,
     328,     0,     0,     0,     0,     0,   166,     0,     0,    70,
      71,    72,    73,    74,    75,   151,   152,   153,   154,   155,
     156,   157,   158,   159,   160,   161,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,   162,     0,     0,    83,
      84,    85,    86,   163,     0,    89,    90,     0,     0,     0,
      92,    93,    94,     0,    60,    61,    62,   164,   165,    65,
     130,     0,     0,     0,   166,     0,     0,     0,     0,     0,
       0,     0,     0,     0,    70,    71,    72,    73,    74,    75,
      76,    77,    78,    79,    80,    81,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,    83,    84,    85,    86,     0,     0,
      89,    90,     0,     0,     0,    92,    93,    94,     0,    60,
      61,    62,    95,    96,    65,   304,     0,     0,     0,   132,
       0,     0,     0,     0,     0,     0,     0,     0,     0,    70,
      71,    72,    73,    74,    75,    76,    77,    78,    79,    80,
      81,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,    83,
      84,    85,    86,     0,     0,    89,    90,     0,     0,     0,
      92,    93,    94,     0,    60,    61,    62,    95,   305,    65,
     130,     0,     0,     0,   132,     0,     0,     0,     0,     0,
       0,     0,     0,     0,    70,    71,    72,    73,    74,    75,
      76,    77,    78,    79,    80,    81,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,    83,    84,    85,    86,     0,     0,
      89,    90,     0,     0,     0,    92,    93,    94,     0,    60,
      61,    62,    95,   713,    65,   130,     0,     0,     0,   132,
       0,     0,     0,     0,     0,     0,     0,     0,     0,    70,
      71,    72,    73,    74,    75,    76,    77,    78,    79,    80,
      81,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,    83,
      84,    85,    86,     0,     0,    89,    90,    60,    61,    62,
      92,    93,    94,   130,     0,     0,     0,    95,   790,     0,
       0,     0,     0,     0,   132,     0,     0,    70,    71,    72,
      73,    74,    75,    76,    77,    78,    79,    80,    81,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,    83,    84,    85,
      86,     0,     0,    89,    90,    60,    61,    62,    92,    93,
      94,   130,     0,     0,     0,    95,    96,   499,     0,     0,
       0,     0,   132,   500,     0,    70,    71,    72,    73,    74,
      75,    76,    77,    78,    79,    80,    81,     0,     0,     0,
       0,     0,     0,   501,   502,   503,   504,   505,   506,   157,
     158,   159,   160,   161,     0,    83,    84,    85,    86,     0,
       0,    89,    90,     0,   162,     0,    92,    93,    94,     0,
       0,   163,     0,     0,     0,     0,     0,     0,     0,     0,
     132,     0,     0,     0,     0,   327,   328,     0,     0,     0,
       0,     0,   507
};

static const yytype_int16 yycheck[] =
{
      36,    37,    95,    96,   148,   186,    50,    43,   473,   136,
     522,   523,   333,   334,    21,   181,    20,    53,   339,   340,
     345,   346,    31,    59,    26,   164,   165,    26,   247,   248,
      66,    59,   252,   253,   254,   255,   256,    16,    17,    18,
      19,    35,    20,    51,    24,    25,    35,    22,    22,   104,
      22,    81,    39,    36,    39,    40,    41,    22,    39,    40,
      41,    97,    98,   114,    55,    52,   117,    52,   111,   112,
      98,   164,   165,    39,   114,   111,   112,   117,   114,   122,
      66,   111,   112,   119,   104,    86,    39,    88,    53,    54,
     113,    56,    57,   129,   130,   131,   132,   120,    51,    52,
     103,   129,   138,   111,   112,   141,   104,   143,   117,   111,
     112,   115,   111,   112,   150,   122,   336,   337,   338,   103,
     122,   341,   113,   122,   118,   111,   112,   308,   114,   118,
     166,   111,   112,   119,   117,   104,   115,   115,   119,   103,
     176,   177,   117,   117,   130,   117,   111,   112,   114,   369,
      32,    33,   372,   618,    39,    40,    41,   377,   479,   195,
     379,   197,   163,   199,   150,   201,   293,   195,   104,   197,
      39,    40,    41,   209,   210,   211,   212,   104,   214,   215,
     216,   217,   218,   219,   104,   221,   222,   223,   327,   328,
     176,   227,   228,   229,    51,   223,   103,   336,    27,   227,
     228,   229,   341,   528,   529,   530,     0,   111,   112,    32,
      33,    34,   305,   249,   726,   119,    32,    33,    34,    35,
      36,    37,    32,    33,    34,   261,   117,   117,    51,   111,
     112,   122,   122,   261,   119,   271,   272,   273,    21,   275,
     276,   277,   278,   279,   280,   114,    32,    33,    34,    35,
      36,    37,    52,    21,   111,   112,   292,    21,   294,   295,
     296,   297,    91,    92,   292,   117,   302,   119,   304,    59,
      60,   307,    32,    33,    34,    35,    36,    37,   314,   315,
     316,   317,   318,   319,    35,   271,   272,   273,   111,   112,
      32,    33,    34,   474,   759,   111,   112,    51,    52,   111,
     112,   111,   112,   119,    32,    33,    42,   119,    44,   345,
     346,   111,   112,    49,   350,    32,    33,    34,   304,   355,
     356,    51,    52,   113,    60,   111,   112,    63,    64,    65,
     120,    67,    68,   119,    51,    86,    21,    88,   374,   118,
     376,    51,   378,    32,    33,    34,    35,    36,    37,    51,
      52,   111,   112,   118,    59,    60,   118,   111,   112,   119,
      53,    54,    98,    56,    57,   111,   112,   403,    52,   111,
     112,   117,   408,   409,   118,    32,    33,    34,    35,    36,
      37,   111,   112,   111,   112,   241,   242,   243,   244,   245,
     246,   427,    51,    35,   111,   112,    53,    54,    52,    56,
      57,   111,   112,   724,    46,    47,    52,    49,   113,   111,
     112,    53,    54,   149,   114,   120,    58,    29,   111,   112,
     456,   457,   111,   112,   460,    33,   119,   111,   112,   465,
      51,   475,   460,   477,   104,   471,    14,    15,   117,   175,
     119,   427,    51,    52,   480,    87,   182,    53,    54,    27,
      56,    57,   111,   112,   111,   112,    98,   111,   112,   195,
     102,   197,   119,   118,   500,   111,   112,    53,    54,   118,
      56,    57,   111,   112,   116,   619,   118,   621,   117,   121,
      52,   111,   112,   649,   114,   471,   522,   523,   654,   655,
     111,   112,   528,   529,   530,   531,   532,   533,   534,   535,
     536,   537,   111,   112,    52,   111,   112,   111,   112,    35,
     114,   111,   112,   119,   500,    52,   552,   553,   554,   119,
      46,    47,   258,    49,    52,   111,   112,    53,    54,    52,
     118,    52,    58,   119,    52,    35,   117,   718,   119,   111,
     112,   577,   578,    51,    52,    53,    46,    47,   117,    49,
     119,   587,   118,    53,    54,    52,   111,   112,    58,    35,
      86,    87,    88,   111,   112,   118,   602,   122,   554,    52,
      46,    47,    98,    49,   111,   112,   102,    53,    54,   615,
     111,   112,    58,   111,   112,   118,   117,    87,   111,   112,
     111,   112,   118,   111,   112,   121,    16,    17,    98,   635,
     636,    52,   102,   118,   546,   547,   548,    51,    52,    53,
     115,    87,    52,    35,   111,   112,   602,   118,   118,   118,
     713,   121,    98,   118,    46,    47,   102,    49,   111,   112,
     118,    53,    54,    56,    57,    58,    58,    39,    40,    41,
     676,   677,   118,   679,   680,   121,    52,   813,   814,    51,
     816,   817,    52,   819,   820,   691,   822,   823,   166,   118,
     111,   112,   698,   118,   700,    87,   702,   117,   117,   177,
     119,   111,   112,    52,   710,   118,    98,    52,    23,   715,
     102,    52,    46,    47,   117,    49,   119,   111,   112,    53,
     726,   118,   728,   729,    58,    21,   118,   790,   122,   121,
      39,    40,    41,   739,   740,   111,   112,   111,   112,   745,
     746,   111,   112,   117,   750,   119,   702,   111,   112,   755,
     114,   118,   166,    87,    23,    32,    33,    34,    35,    36,
      37,   118,   111,   112,    98,   118,   111,   112,   102,   116,
     111,   112,   250,   251,    51,   118,    53,    54,   118,    56,
      57,   111,   112,   789,   118,   791,   792,   121,    20,   119,
     796,   118,   798,   499,   800,   117,   802,   119,   111,   112,
     806,    39,    40,    41,    32,    33,    34,    35,    36,    37,
      61,    62,    50,    32,    33,    34,    35,    36,    37,   731,
     732,   733,   734,    51,   830,   831,   832,   833,   834,   835,
     836,   837,   111,   112,   111,   112,   314,   315,   316,   317,
     318,   319,   320,   321,   322,   323,   324,   325,   114,    32,
      33,    34,    35,    36,    37,   333,   334,   563,   117,   565,
     119,   339,   340,   117,    39,    40,    41,   345,   346,    22,
      53,    54,   350,    56,    57,    50,   116,   355,   356,    39,
      40,    41,    30,   111,   112,   797,    32,   799,   366,   801,
      50,   803,   111,   112,    32,    33,    34,    35,    36,    37,
     314,   315,   316,   317,   318,   319,   320,   321,   322,   323,
     324,   325,   114,    53,    54,   116,    56,    57,   111,   112,
      53,    54,   118,    56,    57,   117,   119,   119,   111,   112,
     111,   112,   346,   111,   112,   117,   350,   119,   119,   117,
      61,   355,   356,   855,   856,   857,   858,   859,   860,   861,
     862,   863,   864,   865,   866,   867,   868,   869,   870,    38,
     111,   112,   116,   117,    43,    44,    45,   114,   119,    48,
      49,   111,   112,   111,   112,   111,   112,   117,   111,   112,
      21,   114,    98,   119,    63,    64,    65,    66,    67,    68,
      69,    70,    71,    72,    73,    74,    39,    40,    41,   111,
     112,   479,   480,   111,   112,   117,    21,    50,   118,    53,
      54,   119,    56,    57,    93,    94,    95,    96,   111,   112,
      99,   100,   111,   112,   117,   104,   105,   106,   118,   507,
     119,   118,   111,   112,   118,    46,    47,   118,    49,   118,
     111,   112,    53,    54,   522,   523,   118,    58,   119,   118,
     528,   529,   530,   531,   532,   533,   534,   535,   536,   537,
     538,   539,   540,   541,   542,   543,   480,   111,   112,    39,
      40,    41,   550,   117,    53,    54,    87,    56,    57,   118,
      50,   111,   112,    51,   114,    53,    54,    98,    56,    57,
     118,   102,   118,   507,    53,    54,   118,    56,    57,    21,
      36,    37,    39,    40,    41,    41,    21,   118,   522,   523,
     121,   111,   112,    50,   114,   111,   112,   531,   532,   533,
     534,   535,   536,   537,   538,   539,   540,   541,   542,   543,
      66,    20,   111,   112,   111,   112,   550,    21,   117,    21,
      36,    37,   119,   111,   112,    46,    47,    43,    49,    39,
      40,    41,   111,   112,    21,   114,    21,    58,    16,    17,
      50,    97,   640,   641,   642,   643,   644,   645,    61,    62,
      66,    39,    40,    41,    35,   111,   112,    20,   114,    39,
      40,    41,    50,   119,   120,   114,    87,    39,    40,    41,
      50,   111,   112,   114,   130,    21,    21,    98,    50,   119,
     104,   102,    98,    53,    54,    21,    56,    57,    52,    53,
      54,    21,    56,    57,   150,   111,   112,   118,   114,    20,
     121,    20,    20,   119,    20,   119,   640,   641,   642,   643,
     644,   645,    53,    54,   130,    56,    57,   111,   112,   119,
     176,   111,   112,   721,   119,   119,   724,   725,   726,   119,
     728,   729,   111,   112,   150,   119,    39,   114,   111,   112,
     119,   111,   112,   199,   114,   201,   119,   111,   112,   118,
     118,   207,   208,   111,   112,   111,   112,   213,   111,   112,
     176,   119,   122,   119,   118,    21,   119,   118,   118,   118,
     111,   112,     3,     4,     5,     6,     7,     8,     9,    10,
      11,    12,    13,   118,   118,   111,   112,   721,   111,   112,
      21,   725,   726,   119,   728,   729,   119,   114,   796,    39,
     798,    39,   800,   117,   802,    82,    83,    84,    85,   265,
     266,   267,   268,   269,   270,   271,   272,   273,   274,    21,
     111,   112,   111,   112,    21,   119,   114,   283,   119,   114,
     119,   114,   830,   831,   832,   833,   834,   835,   836,   837,
      51,   111,   112,   111,   112,   111,   112,    51,   304,   119,
     114,   119,    52,   119,    52,   271,   272,   273,   454,   111,
     112,   606,   796,    22,   798,   310,   800,   119,   802,   285,
     111,   112,   111,   112,   291,   291,    35,   465,   119,   350,
     119,   282,   380,    98,    43,    44,    45,   549,   304,    48,
      49,    16,    17,    18,    19,    -1,   830,   831,   832,   833,
     834,   835,   836,   837,    63,    64,    65,    66,    67,    68,
      69,    70,    71,    72,    73,    74,   111,   112,   111,   112,
      -1,    80,   378,    -1,   119,    -1,   119,    -1,    -1,    -1,
     111,   112,    -1,    -1,    93,    94,    95,    96,   119,    -1,
      99,   100,   111,   112,    -1,   104,   105,   106,   111,   112,
     119,    -1,   111,   112,   111,   112,   119,   111,   112,   118,
      -1,    -1,   119,    -1,    -1,   119,   107,   108,   109,   110,
      -1,   427,    36,    -1,    -1,    -1,    -1,    -1,    -1,    43,
      44,    45,    -1,    -1,    48,    49,    -1,   403,    -1,   405,
     446,    -1,    -1,    -1,    -1,    -1,    -1,   453,   454,    63,
      64,    65,    66,    67,    68,    69,    70,    71,    72,    73,
      74,   427,    -1,    -1,    -1,   471,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    93,
      94,    95,    96,    -1,    -1,    99,   100,    -1,   454,    -1,
     104,   105,   106,    -1,   500,    -1,    -1,   111,   112,   465,
      28,    -1,    -1,    -1,   118,   471,    -1,    35,    -1,    -1,
      -1,    -1,    -1,    -1,    42,    43,    44,    45,    46,    47,
      48,    49,    -1,    -1,    -1,    53,    54,    -1,    -1,    -1,
      58,    -1,    -1,    -1,   500,    63,    64,    65,    66,    67,
      68,    69,    70,    71,    72,    73,    74,    -1,   554,   555,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    87,
      -1,    -1,    -1,    -1,    -1,    93,    94,    95,    96,    -1,
      98,    99,   100,   101,   102,    -1,   104,   105,   106,    -1,
      -1,    -1,    -1,   111,   112,    -1,    -1,    -1,   554,    -1,
     118,    -1,    -1,   121,    -1,    -1,   602,   603,   604,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,   572,    -1,    27,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    43,    44,    45,    -1,    -1,    48,
      49,    -1,    -1,    -1,    -1,    -1,   602,    -1,    -1,    -1,
     606,    -1,    -1,    -1,    63,    64,    65,    66,    67,    68,
      69,    70,    71,    72,    73,    74,    75,    76,    77,    78,
      79,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      89,    90,    -1,    -1,    93,    94,    95,    96,    97,    -1,
      99,   100,    -1,    -1,   690,   104,   105,   106,    -1,    -1,
      -1,    -1,   111,   112,    35,    -1,   702,    -1,    -1,   118,
      -1,    -1,    43,    44,    45,    46,    47,    48,    49,    -1,
      -1,    -1,    53,    54,    -1,    -1,    -1,    58,    -1,    -1,
      -1,    -1,    63,    64,    65,    66,    67,    68,    69,    70,
      71,    72,    73,    74,    -1,    -1,   702,    -1,    -1,    80,
      -1,    -1,    -1,    -1,    -1,    -1,    87,    -1,    -1,    -1,
      -1,    -1,    93,    94,    95,    96,    -1,    98,    99,   100,
      -1,   102,    -1,   104,   105,   106,    -1,    -1,    -1,    -1,
     111,   112,    -1,    -1,    35,    -1,    -1,   118,    -1,    -1,
     121,   122,    43,    44,    45,    46,    47,    48,    49,    -1,
      -1,    -1,    53,    54,    -1,    -1,    -1,    58,    -1,    -1,
      -1,    -1,    63,    64,    65,    66,    67,    68,    69,    70,
      71,    72,    73,    74,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    87,    -1,    -1,    -1,
      -1,    -1,    93,    94,    95,    96,    97,    98,    99,   100,
      -1,   102,    -1,   104,   105,   106,    -1,    -1,    -1,    -1,
     111,   112,    -1,    35,    -1,    -1,    -1,   118,    -1,    -1,
     121,    43,    44,    45,    46,    47,    48,    49,    -1,    -1,
      -1,    53,    54,    -1,    -1,    -1,    58,    -1,    -1,    -1,
      -1,    63,    64,    65,    66,    67,    68,    69,    70,    71,
      72,    73,    74,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    87,    -1,    -1,    -1,    -1,
      -1,    93,    94,    95,    96,    -1,    98,    99,   100,    -1,
     102,    -1,   104,   105,   106,    -1,    -1,    -1,    -1,   111,
     112,    -1,    35,    -1,    -1,    -1,   118,    -1,    -1,   121,
      43,    44,    45,    -1,    -1,    48,    49,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      63,    64,    65,    66,    67,    68,    69,    70,    71,    72,
      73,    74,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      93,    94,    95,    96,    -1,    -1,    99,   100,    -1,    -1,
      38,   104,   105,   106,    -1,    43,    44,    45,   111,   112,
      48,    49,    -1,    -1,    -1,   118,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    63,    64,    65,    66,    67,
      68,    69,    70,    71,    72,    73,    74,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    93,    94,    95,    96,    -1,
      -1,    99,   100,    -1,    -1,    38,   104,   105,   106,    -1,
      43,    44,    45,   111,   112,    48,    49,    -1,    -1,    -1,
     118,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      63,    64,    65,    66,    67,    68,    69,    70,    71,    72,
      73,    74,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      93,    94,    95,    96,    -1,    -1,    99,   100,    -1,    38,
      -1,   104,   105,   106,    43,    44,    45,    -1,   111,   112,
      49,    -1,    -1,    -1,    -1,   118,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    63,    64,    65,    66,    67,    68,
      69,    70,    71,    72,    73,    74,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    93,    94,    95,    96,    -1,    -1,
      99,   100,    -1,    38,    -1,   104,   105,   106,    43,    44,
      45,    -1,    -1,    -1,    49,    -1,    -1,    -1,    -1,   118,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    63,    64,
      65,    66,    67,    68,    69,    70,    71,    72,    73,    74,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    93,    94,
      95,    96,    -1,    -1,    99,   100,    -1,    38,    -1,   104,
     105,   106,    43,    44,    45,    -1,    -1,    -1,    49,    -1,
      -1,    -1,    -1,   118,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    63,    64,    65,    66,    67,    68,    69,    70,
      71,    72,    73,    74,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    93,    94,    95,    96,    -1,    -1,    99,   100,
      -1,    -1,    -1,   104,   105,   106,    -1,    -1,    43,    44,
      45,    46,    47,    48,    49,    -1,    -1,   118,    53,    54,
      -1,    -1,    -1,    58,    -1,    -1,    -1,    -1,    63,    64,
      65,    66,    67,    68,    69,    70,    71,    72,    73,    74,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    87,    -1,    -1,    -1,    -1,    -1,    93,    94,
      95,    96,    -1,    98,    99,   100,    -1,   102,    -1,   104,
     105,   106,    42,    43,    44,    45,   111,   112,    48,    49,
      -1,    -1,    -1,   118,    -1,    -1,   121,    -1,    -1,    -1,
      -1,    -1,    -1,    63,    64,    65,    66,    67,    68,    69,
      70,    71,    72,    73,    74,    75,    76,    77,    78,    79,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      90,    -1,    -1,    93,    94,    95,    96,    97,    -1,    99,
     100,    -1,    -1,    -1,   104,   105,   106,    -1,    43,    44,
      45,   111,   112,    48,    49,    -1,    -1,    -1,   118,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    63,    64,
      65,    66,    67,    68,    69,    70,    71,    72,    73,    74,
      75,    76,    77,    78,    79,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    89,    90,    -1,    -1,    93,    94,
      95,    96,    97,    -1,    99,   100,    -1,    -1,    -1,   104,
     105,   106,    -1,    43,    44,    45,   111,   112,    48,    49,
      -1,    -1,    -1,   118,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    63,    64,    65,    66,    67,    68,    69,
      70,    71,    72,    73,    74,    75,    76,    77,    78,    79,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      90,    -1,    -1,    93,    94,    95,    96,    97,    -1,    99,
     100,    43,    44,    45,   104,   105,   106,    49,    -1,    -1,
      -1,   111,   112,    -1,    -1,    -1,    -1,    -1,   118,    -1,
      -1,    63,    64,    65,    66,    67,    68,    69,    70,    71,
      72,    73,    74,    75,    76,    77,    78,    79,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    90,    -1,
      -1,    93,    94,    95,    96,    97,    -1,    99,   100,    43,
      44,    45,   104,   105,   106,    49,    -1,    -1,    -1,   111,
     112,    -1,    -1,    -1,    -1,    -1,   118,    -1,    -1,    63,
      64,    65,    66,    67,    68,    69,    70,    71,    72,    73,
      74,    75,    76,    77,    78,    79,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    90,    -1,    -1,    93,
      94,    95,    96,    97,    -1,    99,   100,    -1,    -1,    -1,
     104,   105,   106,    -1,    43,    44,    45,   111,   112,    48,
      49,    -1,    -1,    -1,   118,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    63,    64,    65,    66,    67,    68,
      69,    70,    71,    72,    73,    74,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    93,    94,    95,    96,    -1,    -1,
      99,   100,    -1,    -1,    -1,   104,   105,   106,    -1,    43,
      44,    45,   111,   112,    48,    49,    -1,    -1,    -1,   118,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    63,
      64,    65,    66,    67,    68,    69,    70,    71,    72,    73,
      74,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    93,
      94,    95,    96,    -1,    -1,    99,   100,    -1,    -1,    -1,
     104,   105,   106,    -1,    43,    44,    45,   111,   112,    48,
      49,    -1,    -1,    -1,   118,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    63,    64,    65,    66,    67,    68,
      69,    70,    71,    72,    73,    74,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    93,    94,    95,    96,    -1,    -1,
      99,   100,    -1,    -1,    -1,   104,   105,   106,    -1,    43,
      44,    45,   111,   112,    48,    49,    -1,    -1,    -1,   118,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    63,
      64,    65,    66,    67,    68,    69,    70,    71,    72,    73,
      74,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    93,
      94,    95,    96,    -1,    -1,    99,   100,    43,    44,    45,
     104,   105,   106,    49,    -1,    -1,    -1,   111,   112,    -1,
      -1,    -1,    -1,    -1,   118,    -1,    -1,    63,    64,    65,
      66,    67,    68,    69,    70,    71,    72,    73,    74,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    93,    94,    95,
      96,    -1,    -1,    99,   100,    43,    44,    45,   104,   105,
     106,    49,    -1,    -1,    -1,   111,   112,    43,    -1,    -1,
      -1,    -1,   118,    49,    -1,    63,    64,    65,    66,    67,
      68,    69,    70,    71,    72,    73,    74,    -1,    -1,    -1,
      -1,    -1,    -1,    69,    70,    71,    72,    73,    74,    75,
      76,    77,    78,    79,    -1,    93,    94,    95,    96,    -1,
      -1,    99,   100,    -1,    90,    -1,   104,   105,   106,    -1,
      -1,    97,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
     118,    -1,    -1,    -1,    -1,   111,   112,    -1,    -1,    -1,
      -1,    -1,   118
};

/* YYSTOS[STATE-NUM] -- The (internal number of the) accessing
   symbol of state STATE-NUM.  */
static const yytype_uint8 yystos[] =
{
       0,     3,     4,     5,     6,     7,     8,     9,    10,    11,
      12,    13,    21,   124,   125,   128,   129,   130,   131,   133,
     136,   146,   147,   159,   162,   104,   104,   104,   104,   104,
     104,   104,   103,   103,   103,   103,    14,    15,    27,   163,
       0,    20,   115,    20,   115,    16,    17,    18,    19,   115,
     137,    21,    21,    21,    21,   118,   118,   118,   118,    35,
      43,    44,    45,    46,    47,    48,    49,    53,    54,    58,
      63,    64,    65,    66,    67,    68,    69,    70,    71,    72,
      73,    74,    87,    93,    94,    95,    96,    97,    98,    99,
     100,   102,   104,   105,   106,   111,   112,   118,   121,   166,
     167,   168,   169,   174,   176,   177,   178,   179,   180,   181,
      28,    42,    49,   101,   118,   166,   173,   174,   177,    49,
     118,   164,   165,   166,   174,   114,   166,   116,   164,    22,
      49,    80,   118,   134,   142,   143,   144,   170,   174,   177,
     164,    29,   140,    33,   138,    16,    17,   164,   138,    43,
      49,    69,    70,    71,    72,    73,    74,    75,    76,    77,
      78,    79,    90,    97,   111,   112,   118,   154,   155,   156,
     157,   158,   177,   178,   154,    27,    49,    89,   148,   149,
     154,   177,    27,    91,    92,   160,   161,   104,   132,   132,
     132,   132,    36,   176,   164,   118,   165,   118,   165,   118,
     164,   118,   164,   164,   173,   164,   164,   118,   118,   118,
     118,   118,   118,   118,   118,   118,   118,   118,   118,   118,
     118,   118,   118,   115,   175,   175,   175,   118,   118,   118,
     180,   180,   166,   177,   122,   164,   170,   172,   174,   176,
     177,    53,    54,    56,    57,   111,   112,    55,   113,   117,
     111,   112,    59,    60,   113,   120,    61,    62,   118,   173,
     173,   118,   166,   173,   177,    32,    33,    34,    35,    36,
      37,    39,    40,    41,    23,    32,    33,    34,    35,    36,
      37,   173,    21,    23,   114,    20,   116,   176,   173,   177,
     114,   117,   144,   145,    82,    83,    84,    85,   171,   177,
     116,   177,    30,   141,    49,   112,   177,    32,   139,   114,
     116,   139,   164,   173,   118,   118,   118,   118,   118,   118,
     118,   118,   118,   118,   118,   118,   175,   111,   112,   156,
     156,   154,   177,   111,   112,   114,   113,   120,    61,   111,
     112,   113,   114,   164,   173,    42,   118,   150,   154,   177,
      39,   114,    32,    33,    34,   153,   153,   164,   114,   140,
     117,   119,   119,   119,   119,    36,    21,   164,   176,    21,
     164,   176,    21,   177,    21,   177,    21,    21,    50,    21,
      21,   166,   166,   177,   177,   177,   177,   166,   177,   177,
     177,   177,   177,   177,    98,   177,   177,   176,   176,   176,
     176,   119,   119,    21,   122,   117,   122,   122,    24,    25,
     167,   167,   167,   167,   167,   167,   169,   169,   177,   178,
     178,   179,   179,   179,   179,   179,   164,    50,   176,   119,
     166,   166,   166,   166,   166,   166,   173,   173,   173,   166,
     177,   177,   177,   177,   177,   177,    50,   163,   166,    86,
      88,   126,   127,   174,    20,    22,    50,    81,   143,   176,
     144,   177,   177,   177,   177,    20,   177,   114,   173,    38,
      38,    49,   177,   140,    16,    17,    19,   137,   114,    21,
      50,   154,   177,   154,   177,   154,   177,   154,   177,   154,
     177,   154,   177,   154,   154,   154,   154,   154,   154,    43,
      49,    69,    70,    71,    72,    73,    74,   118,   119,   155,
     178,   155,   178,   156,   179,   179,   179,   155,   178,   155,
     178,   156,    21,    50,   150,   150,   154,   177,    39,    40,
      41,    50,    32,    33,    34,    35,    36,    37,    32,    33,
      34,    35,    36,    37,   149,   154,   177,   154,   177,    21,
      21,   104,    20,    20,    20,    20,   178,   119,   119,   179,
     119,   119,   179,   119,   177,   119,   177,   179,   166,   177,
     169,   168,   117,   119,   119,   119,   119,   117,   117,   119,
     119,   119,   119,   119,   119,   119,   119,   117,   116,   119,
     119,   119,   174,   177,   174,   177,   177,   119,   166,   173,
     177,   119,    22,   118,   118,   114,   117,   166,   126,   166,
     177,    22,   176,   142,   177,    50,   173,   141,   140,   138,
      16,   138,   155,   178,   154,   177,   119,   119,   119,   119,
     119,   119,   119,   119,   119,   117,   117,   119,   164,   173,
     118,   118,   118,   118,   118,   118,   148,   148,   154,   177,
     119,   150,   150,   150,   154,   177,   154,   177,   154,   177,
     154,   177,   154,   177,   154,   177,   154,   177,   154,   154,
     154,   154,   154,   154,   151,   151,    32,    33,   151,    32,
      33,   151,   160,   154,   177,   177,   173,   166,   164,   164,
      51,    51,   174,   177,   177,   177,   122,   122,    26,   122,
      26,   122,    51,   173,   166,   166,   127,   114,   114,    22,
      31,   135,   114,   112,   177,    50,   114,   141,   139,   114,
     139,    51,   177,   177,    21,    50,    51,    52,   153,   153,
     117,   177,   177,   177,   177,   114,   114,   114,   114,    21,
      21,   166,   177,   119,   119,   117,   117,   177,   177,   173,
     117,   119,   177,   114,    38,    51,    38,   177,   114,   140,
     114,   154,   119,   119,   178,   148,   154,   177,   154,   177,
     107,   108,   109,   110,   152,   151,   151,   151,   151,   177,
     177,    52,    52,   177,   177,   122,   122,    52,   177,    51,
     112,    51,    51,   141,    52,    52,    51,    52,    51,    52,
      51,    52,    51,    52,   119,   119,   117,   119,   177,    38,
     177,    38,   114,   154,   177,   151,   154,   177,   151,   154,
     177,   151,   154,   177,   151,   177,    52,    52,    52,    52,
     153,   153,   153,   153,   153,   153,   153,   153,   119,   154,
     177,   154,   177,   154,   177,   154,   177,   154,   177,   154,
     177,   154,   177,   154,   177,    52,    52,    52,    52,    52,
      52,    52,    52,    52,    52,    52,    52,    52,    52,    52,
      52,   151,   151,   151,   151,   151,   151,   151,   151,   151,
     151,   151,   151,   151,   151,   151,   151
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
yy_stack_print (yytype_int16 *bottom, yytype_int16 *top)
#else
static void
yy_stack_print (bottom, top)
    yytype_int16 *bottom;
    yytype_int16 *top;
#endif
{
  YYFPRINTF (stderr, "Stack now");
  for (; bottom <= top; ++bottom)
    YYFPRINTF (stderr, " %d", *bottom);
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
      fprintf (stderr, "   $%d = ", yyi + 1);
      yy_symbol_print (stderr, yyrhs[yyprhs[yyrule] + yyi],
		       &(yyvsp[(yyi + 1) - (yynrhs)])
		       		       );
      fprintf (stderr, "\n");
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






/*----------.
| yyparse.  |
`----------*/

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
  /* The look-ahead symbol.  */
int yychar;

/* The semantic value of the look-ahead symbol.  */
YYSTYPE yylval;

/* Number of syntax errors so far.  */
int yynerrs;

  int yystate;
  int yyn;
  int yyresult;
  /* Number of tokens to shift before error messages enabled.  */
  int yyerrstatus;
  /* Look-ahead token as an internal (translated) token number.  */
  int yytoken = 0;
#if YYERROR_VERBOSE
  /* Buffer for error messages, and its allocated size.  */
  char yymsgbuf[128];
  char *yymsg = yymsgbuf;
  YYSIZE_T yymsg_alloc = sizeof yymsgbuf;
#endif

  /* Three stacks and their tools:
     `yyss': related to states,
     `yyvs': related to semantic values,
     `yyls': related to locations.

     Refer to the stacks thru separate pointers, to allow yyoverflow
     to reallocate them elsewhere.  */

  /* The state stack.  */
  yytype_int16 yyssa[YYINITDEPTH];
  yytype_int16 *yyss = yyssa;
  yytype_int16 *yyssp;

  /* The semantic value stack.  */
  YYSTYPE yyvsa[YYINITDEPTH];
  YYSTYPE *yyvs = yyvsa;
  YYSTYPE *yyvsp;



#define YYPOPSTACK(N)   (yyvsp -= (N), yyssp -= (N))

  YYSIZE_T yystacksize = YYINITDEPTH;

  /* The variables used to return semantic value and location from the
     action routines.  */
  YYSTYPE yyval;


  /* The number of symbols on the RHS of the reduced rule.
     Keep to zero when no symbol should be popped.  */
  int yylen = 0;

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

  /* Do appropriate processing given the current state.  Read a
     look-ahead token if we need one and don't already have one.  */

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

  /* Count tokens shifted since error; after three, turn off error
     status.  */
  if (yyerrstatus)
    yyerrstatus--;

  /* Shift the look-ahead token.  */
  YY_SYMBOL_PRINT ("Shifting", yytoken, &yylval, &yylloc);

  /* Discard the shifted token unless it is eof.  */
  if (yychar != YYEOF)
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
#line 148 "src/mmlparse2.y"
    { code_set_root((yyvsp[(1) - (1)].code)); ;}
    break;

  case 3:
#line 149 "src/mmlparse2.y"
    { code_set_root((yyvsp[(1) - (1)].code)); ;}
    break;

  case 4:
#line 150 "src/mmlparse2.y"
    { code_set_root((yyvsp[(1) - (1)].code)); ;}
    break;

  case 5:
#line 151 "src/mmlparse2.y"
    { code_set_root((yyvsp[(1) - (1)].code)); ;}
    break;

  case 6:
#line 152 "src/mmlparse2.y"
    { code_set_root((yyvsp[(1) - (1)].code)); ;}
    break;

  case 7:
#line 153 "src/mmlparse2.y"
    { code_set_root((yyvsp[(1) - (1)].code)); ;}
    break;

  case 8:
#line 154 "src/mmlparse2.y"
    { code_set_root((yyvsp[(1) - (1)].code)); ;}
    break;

  case 9:
#line 155 "src/mmlparse2.y"
    { code_set_root((yyvsp[(1) - (1)].code)); ;}
    break;

  case 10:
#line 156 "src/mmlparse2.y"
    { code_set_root((yyvsp[(1) - (1)].code)); ;}
    break;

  case 11:
#line 157 "src/mmlparse2.y"
    { code_set_root((yyvsp[(1) - (1)].code)); ;}
    break;

  case 12:
#line 158 "src/mmlparse2.y"
    { code_set_root((yyvsp[(1) - (1)].code)); ;}
    break;

  case 13:
#line 166 "src/mmlparse2.y"
    {
         (yyval.code) = code_new_inst(i_newsym_set1, 3,
            code_new_name((yyvsp[(2) - (5)].name)),                                       /* Name */
            code_new_inst(i_idxset_pseudo_new, 1,               /* index set */
               code_new_inst(i_bool_true, 0)),              
            (yyvsp[(4) - (5)].code));                                              /* initial set */
      ;}
    break;

  case 14:
#line 173 "src/mmlparse2.y"
    {
         (yyval.code) = code_new_inst(i_newsym_set1, 3,
            code_new_name((yyvsp[(2) - (8)].name)),                                       /* Name */
            (yyvsp[(4) - (8)].code),                                                 /* index set */
            (yyvsp[(7) - (8)].code));                                                      /* set */
      ;}
    break;

  case 15:
#line 179 "src/mmlparse2.y"
    {
         (yyval.code) = code_new_inst(i_newsym_set2, 3,
            code_new_name((yyvsp[(2) - (8)].name)),                                       /* Name */
            (yyvsp[(4) - (8)].code),                                                 /* index set */
            (yyvsp[(7) - (8)].code));                                   /* initial set_entry_list */
      ;}
    break;

  case 16:
#line 185 "src/mmlparse2.y"
    {
         (yyval.code) = code_new_inst(i_newsym_set2, 3,
            code_new_name((yyvsp[(2) - (7)].name)),                                       /* Name */
            code_new_inst(i_idxset_pseudo_new, 1,               /* index set */
               code_new_inst(i_bool_true, 0)),              
            (yyvsp[(6) - (7)].code));                                   /* initial set_entry_list */
      ;}
    break;

  case 17:
#line 195 "src/mmlparse2.y"
    { (yyval.code) = code_new_inst(i_entry_list_new, 1, (yyvsp[(1) - (1)].code)); ;}
    break;

  case 18:
#line 196 "src/mmlparse2.y"
    {
         (yyval.code) = code_new_inst(i_entry_list_add, 2, (yyvsp[(1) - (3)].code), (yyvsp[(3) - (3)].code));
      ;}
    break;

  case 19:
#line 199 "src/mmlparse2.y"
    {
         (yyval.code) = code_new_inst(i_entry_list_subsets, 3, (yyvsp[(3) - (6)].code), (yyvsp[(5) - (6)].code), code_new_numb(numb_new_integer(-1)));
      ;}
    break;

  case 20:
#line 202 "src/mmlparse2.y"
    {
         (yyval.code) = code_new_inst(i_entry_list_subsets, 3, (yyvsp[(3) - (8)].code), (yyvsp[(5) - (8)].code), (yyvsp[(7) - (8)].code));
      ;}
    break;

  case 21:
#line 205 "src/mmlparse2.y"
    {
         (yyval.code) = code_new_inst(i_entry_list_powerset, 1, (yyvsp[(3) - (4)].code));
      ;}
    break;

  case 22:
#line 211 "src/mmlparse2.y"
    { (yyval.code) = code_new_inst(i_entry, 2, (yyvsp[(1) - (2)].code), (yyvsp[(2) - (2)].code)); ;}
    break;

  case 23:
#line 220 "src/mmlparse2.y"
    {
         (yyval.code) = code_new_inst(i_newdef, 3,
            code_new_define((yyvsp[(2) - (8)].def)),
            code_new_inst(i_tuple_new, 1, (yyvsp[(4) - (8)].code)),
            (yyvsp[(7) - (8)].code));
      ;}
    break;

  case 24:
#line 229 "src/mmlparse2.y"
    {
         (yyval.code) = code_new_inst(i_newdef, 3,
            code_new_define((yyvsp[(2) - (8)].def)),
            code_new_inst(i_tuple_new, 1, (yyvsp[(4) - (8)].code)),
            (yyvsp[(7) - (8)].code));
      ;}
    break;

  case 25:
#line 238 "src/mmlparse2.y"
    {
         (yyval.code) = code_new_inst(i_newdef, 3,
            code_new_define((yyvsp[(2) - (8)].def)),
            code_new_inst(i_tuple_new, 1, (yyvsp[(4) - (8)].code)),
            (yyvsp[(7) - (8)].code));
      ;}
    break;

  case 26:
#line 247 "src/mmlparse2.y"
    {
         (yyval.code) = code_new_inst(i_newdef, 3,
            code_new_define((yyvsp[(2) - (8)].def)),
            code_new_inst(i_tuple_new, 1, (yyvsp[(4) - (8)].code)),
            (yyvsp[(7) - (8)].code));
      ;}
    break;

  case 27:
#line 256 "src/mmlparse2.y"
    {
         (yyval.code) = code_new_inst(i_elem_list_new, 1, code_new_name((yyvsp[(1) - (1)].name)));
      ;}
    break;

  case 28:
#line 259 "src/mmlparse2.y"
    {
         (yyval.code) = code_new_inst(i_elem_list_add, 2, (yyvsp[(1) - (3)].code), code_new_name((yyvsp[(3) - (3)].name)));
      ;}
    break;

  case 29:
#line 269 "src/mmlparse2.y"
    {
         (yyval.code) = code_new_inst(i_newsym_para1, 4, code_new_name((yyvsp[(2) - (9)].name)), (yyvsp[(4) - (9)].code), (yyvsp[(7) - (9)].code), (yyvsp[(8) - (9)].code));
      ;}
    break;

  case 30:
#line 272 "src/mmlparse2.y"
    {
         (yyval.code) = code_new_inst(i_newsym_para2, 4, code_new_name((yyvsp[(2) - (8)].name)), (yyvsp[(4) - (8)].code), (yyvsp[(7) - (8)].code), code_new_inst(i_nop, 0));
      ;}
    break;

  case 31:
#line 275 "src/mmlparse2.y"
    {
         (yyval.code) = code_new_inst(i_newsym_para1, 4,
            code_new_name((yyvsp[(2) - (5)].name)),
            code_new_inst(i_idxset_pseudo_new, 1, code_new_inst(i_bool_true, 0)),
            (yyvsp[(4) - (5)].code),
            code_new_inst(i_nop, 0));
      ;}
    break;

  case 32:
#line 282 "src/mmlparse2.y"
    { (yyval.code) = code_new_inst(i_nop, 0); ;}
    break;

  case 33:
#line 285 "src/mmlparse2.y"
    { (yyval.code) = (yyvsp[(1) - (1)].code); ;}
    break;

  case 34:
#line 286 "src/mmlparse2.y"
    {
         (yyval.code) = code_new_inst(i_entry_list_new, 1,
            code_new_inst(i_entry, 2, code_new_inst(i_tuple_empty, 0), (yyvsp[(1) - (1)].code)));
      ;}
    break;

  case 35:
#line 293 "src/mmlparse2.y"
    { (yyval.code) = code_new_inst(i_nop, 0); ;}
    break;

  case 36:
#line 294 "src/mmlparse2.y"
    { (yyval.code) = code_new_inst(i_entry, 2, code_new_inst(i_tuple_empty, 0), (yyvsp[(2) - (2)].code)); ;}
    break;

  case 37:
#line 302 "src/mmlparse2.y"
    {
         (yyval.code) = code_new_inst(i_newsym_var, 7,
            code_new_name((yyvsp[(2) - (9)].name)),
            (yyvsp[(4) - (9)].code), (yyvsp[(6) - (9)].code), (yyvsp[(7) - (9)].code), (yyvsp[(8) - (9)].code),
            code_new_numb(numb_copy(numb_unknown())),
            code_new_numb(numb_copy(numb_unknown())));
      ;}
    break;

  case 38:
#line 309 "src/mmlparse2.y"
    {
         (yyval.code) = code_new_inst(i_newsym_var, 7,
            code_new_name((yyvsp[(2) - (6)].name)),
            code_new_inst(i_idxset_pseudo_new, 1,
               code_new_inst(i_bool_true, 0)),              
            (yyvsp[(3) - (6)].code), (yyvsp[(4) - (6)].code), (yyvsp[(5) - (6)].code),
            code_new_numb(numb_copy(numb_unknown())),
            code_new_numb(numb_copy(numb_unknown())));
      ;}
    break;

  case 39:
#line 318 "src/mmlparse2.y"
    {
         (yyval.code) = code_new_inst(i_newsym_var, 7,
            code_new_name((yyvsp[(2) - (8)].name)),
            (yyvsp[(4) - (8)].code),
            code_new_varclass(VAR_IMP),
            code_new_inst(i_bound_new, 1, code_new_numb(numb_new_integer(0))),
            code_new_inst(i_bound_new, 1, code_new_numb(numb_new_integer(1))),
            code_new_numb(numb_copy(numb_unknown())),
            code_new_numb(numb_copy(numb_unknown())));
      ;}
    break;

  case 40:
#line 328 "src/mmlparse2.y"
    {
         (yyval.code) = code_new_inst(i_newsym_var, 7,
            code_new_name((yyvsp[(2) - (5)].name)),
            code_new_inst(i_idxset_pseudo_new, 1,
               code_new_inst(i_bool_true, 0)),              
            code_new_varclass(VAR_IMP),
            code_new_inst(i_bound_new, 1, code_new_numb(numb_new_integer(0))),
            code_new_inst(i_bound_new, 1, code_new_numb(numb_new_integer(1))),
            code_new_numb(numb_copy(numb_unknown())),
            code_new_numb(numb_copy(numb_unknown())));
      ;}
    break;

  case 41:
#line 339 "src/mmlparse2.y"
    {
         (yyval.code) = code_new_inst(i_newsym_var, 7,
            code_new_name((yyvsp[(2) - (9)].name)),
            (yyvsp[(4) - (9)].code),
            code_new_varclass(VAR_INT),
            code_new_inst(i_bound_new, 1, code_new_numb(numb_new_integer(0))),
            code_new_inst(i_bound_new, 1, code_new_numb(numb_new_integer(1))),
            (yyvsp[(7) - (9)].code), (yyvsp[(8) - (9)].code));
      ;}
    break;

  case 42:
#line 348 "src/mmlparse2.y"
    {
         (yyval.code) = code_new_inst(i_newsym_var, 7,
            code_new_name((yyvsp[(2) - (6)].name)),
            code_new_inst(i_idxset_pseudo_new, 1,
               code_new_inst(i_bool_true, 0)),              
            code_new_varclass(VAR_INT),
            code_new_inst(i_bound_new, 1, code_new_numb(numb_new_integer(0))),
            code_new_inst(i_bound_new, 1, code_new_numb(numb_new_integer(1))),
            (yyvsp[(4) - (6)].code), (yyvsp[(5) - (6)].code));
      ;}
    break;

  case 43:
#line 358 "src/mmlparse2.y"
    {
         (yyval.code) = code_new_inst(i_newsym_var, 7,
            code_new_name((yyvsp[(2) - (11)].name)), (yyvsp[(4) - (11)].code), code_new_varclass(VAR_INT), (yyvsp[(7) - (11)].code), (yyvsp[(8) - (11)].code), (yyvsp[(9) - (11)].code), (yyvsp[(10) - (11)].code));
      ;}
    break;

  case 44:
#line 362 "src/mmlparse2.y"
    {
         (yyval.code) = code_new_inst(i_newsym_var, 7,
            code_new_name((yyvsp[(2) - (8)].name)),
            code_new_inst(i_idxset_pseudo_new, 1,
               code_new_inst(i_bool_true, 0)),              
            code_new_varclass(VAR_INT), (yyvsp[(4) - (8)].code), (yyvsp[(5) - (8)].code), (yyvsp[(6) - (8)].code), (yyvsp[(7) - (8)].code));
      ;}
    break;

  case 45:
#line 372 "src/mmlparse2.y"
    { (yyval.code) = code_new_varclass(VAR_CON); ;}
    break;

  case 46:
#line 373 "src/mmlparse2.y"
    { (yyval.code) = code_new_varclass(VAR_CON); ;}
    break;

  case 47:
#line 374 "src/mmlparse2.y"
    { (yyval.code) = code_new_varclass(VAR_IMP); ;}
    break;

  case 48:
#line 378 "src/mmlparse2.y"
    {
         (yyval.code) = code_new_inst(i_bound_new, 1, code_new_numb(numb_new_integer(0)));
      ;}
    break;

  case 49:
#line 381 "src/mmlparse2.y"
    { (yyval.code) = code_new_inst(i_bound_new, 1, (yyvsp[(2) - (2)].code)); ;}
    break;

  case 50:
#line 382 "src/mmlparse2.y"
    { (yyval.code) = code_new_bound(BOUND_MINUS_INFTY); ;}
    break;

  case 51:
#line 383 "src/mmlparse2.y"
    {
         (yyval.code) = code_new_inst(i_expr_if_else, 3, (yyvsp[(3) - (9)].code),
            code_new_inst(i_bound_new, 1, (yyvsp[(5) - (9)].code)),
            code_new_bound(BOUND_MINUS_INFTY));
      ;}
    break;

  case 52:
#line 388 "src/mmlparse2.y"
    {
         (yyval.code) = code_new_inst(i_expr_if_else, 3, (yyvsp[(3) - (9)].code),
            code_new_bound(BOUND_MINUS_INFTY),
            code_new_inst(i_bound_new, 1, (yyvsp[(8) - (9)].code)));
      ;}
    break;

  case 53:
#line 396 "src/mmlparse2.y"
    { (yyval.code) = code_new_bound(BOUND_INFTY); ;}
    break;

  case 54:
#line 397 "src/mmlparse2.y"
    { (yyval.code) = code_new_inst(i_bound_new, 1, (yyvsp[(2) - (2)].code)); ;}
    break;

  case 55:
#line 398 "src/mmlparse2.y"
    { (yyval.code) = code_new_bound(BOUND_INFTY); ;}
    break;

  case 56:
#line 399 "src/mmlparse2.y"
    {
         (yyval.code) = code_new_inst(i_expr_if_else, 3, (yyvsp[(3) - (8)].code),
            code_new_inst(i_bound_new, 1, (yyvsp[(5) - (8)].code)),
            code_new_bound(BOUND_INFTY));
      ;}
    break;

  case 57:
#line 404 "src/mmlparse2.y"
    {
         (yyval.code) = code_new_inst(i_expr_if_else, 3, (yyvsp[(3) - (8)].code),
            code_new_bound(BOUND_INFTY),
            code_new_inst(i_bound_new, 1, (yyvsp[(7) - (8)].code)));
      ;}
    break;

  case 58:
#line 412 "src/mmlparse2.y"
    { (yyval.code) = code_new_numb(numb_new_integer(0)); ;}
    break;

  case 59:
#line 413 "src/mmlparse2.y"
    { (yyval.code) = (yyvsp[(2) - (2)].code); ;}
    break;

  case 60:
#line 417 "src/mmlparse2.y"
    { (yyval.code) = code_new_numb(numb_copy(numb_unknown())); ;}
    break;

  case 61:
#line 418 "src/mmlparse2.y"
    { (yyval.code) = (yyvsp[(2) - (2)].code); ;}
    break;

  case 62:
#line 426 "src/mmlparse2.y"
    { (yyval.code) = code_new_inst(i_entry_list_new, 1, (yyvsp[(1) - (1)].code)); ;}
    break;

  case 63:
#line 427 "src/mmlparse2.y"
    {
         (yyval.code) = code_new_inst(i_entry_list_add, 2, (yyvsp[(1) - (3)].code), (yyvsp[(3) - (3)].code));
      ;}
    break;

  case 64:
#line 430 "src/mmlparse2.y"
    { (yyval.code) = code_new_inst(i_read, 1, (yyvsp[(1) - (1)].code)); ;}
    break;

  case 65:
#line 431 "src/mmlparse2.y"
    { (yyval.code) = code_new_inst(i_list_matrix, 2, (yyvsp[(1) - (2)].code), (yyvsp[(2) - (2)].code)); ;}
    break;

  case 66:
#line 435 "src/mmlparse2.y"
    { (yyval.code) = code_new_inst(i_entry, 2, (yyvsp[(1) - (2)].code), (yyvsp[(2) - (2)].code)); ;}
    break;

  case 67:
#line 439 "src/mmlparse2.y"
    { (yyval.code) = (yyvsp[(2) - (3)].code); ;}
    break;

  case 68:
#line 443 "src/mmlparse2.y"
    {
         (yyval.code) = code_new_inst(i_matrix_list_new, 2, (yyvsp[(1) - (3)].code), (yyvsp[(2) - (3)].code));
      ;}
    break;

  case 69:
#line 446 "src/mmlparse2.y"
    {
         (yyval.code) = code_new_inst(i_matrix_list_add, 3, (yyvsp[(1) - (4)].code), (yyvsp[(2) - (4)].code), (yyvsp[(3) - (4)].code));
      ;}
    break;

  case 70:
#line 458 "src/mmlparse2.y"
    {
         (yyval.code) = code_new_inst(i_object_min, 2, code_new_name((yyvsp[(2) - (5)].name)), (yyvsp[(4) - (5)].code));
      ;}
    break;

  case 71:
#line 461 "src/mmlparse2.y"
    {
         (yyval.code) = code_new_inst(i_object_max, 2, code_new_name((yyvsp[(2) - (5)].name)), (yyvsp[(4) - (5)].code));
      ;}
    break;

  case 72:
#line 471 "src/mmlparse2.y"
    {
        (yyval.code) = code_new_inst(i_subto, 2, code_new_name((yyvsp[(2) - (5)].name)), (yyvsp[(4) - (5)].code));
     ;}
    break;

  case 73:
#line 477 "src/mmlparse2.y"
    {
        (yyval.code) = code_new_inst(i_constraint_list, 2, (yyvsp[(1) - (1)].code), code_new_inst(i_nop, 0));
     ;}
    break;

  case 74:
#line 480 "src/mmlparse2.y"
    {
        (yyval.code) = code_new_inst(i_constraint_list, 2, (yyvsp[(1) - (3)].code), (yyvsp[(3) - (3)].code));
     ;}
    break;

  case 75:
#line 483 "src/mmlparse2.y"
    {
        (yyval.code) = code_new_inst(i_constraint_list, 2, 
           code_new_inst(i_forall, 2, (yyvsp[(2) - (4)].code), (yyvsp[(4) - (4)].code)),
           code_new_inst(i_nop, 0));
     ;}
    break;

  case 76:
#line 488 "src/mmlparse2.y"
    {
        (yyval.code) = code_new_inst(i_constraint_list, 2, 
           code_new_inst(i_expr_if_else, 3, (yyvsp[(2) - (5)].code), (yyvsp[(4) - (5)].code), code_new_inst(i_nop, 0)),
           code_new_inst(i_nop, 0));
      ;}
    break;

  case 77:
#line 493 "src/mmlparse2.y"
    {
        (yyval.code) = code_new_inst(i_constraint_list, 2, 
           code_new_inst(i_expr_if_else, 3, (yyvsp[(2) - (7)].code), (yyvsp[(4) - (7)].code), (yyvsp[(6) - (7)].code)),
           code_new_inst(i_nop, 0));
      ;}
    break;

  case 78:
#line 501 "src/mmlparse2.y"
    {
        (yyval.code) = code_new_inst(i_constraint, 4, (yyvsp[(1) - (4)].code), (yyvsp[(2) - (4)].code), (yyvsp[(3) - (4)].code), code_new_bits((yyvsp[(4) - (4)].bits)));
     ;}
    break;

  case 79:
#line 504 "src/mmlparse2.y"
    {
        (yyval.code) = code_new_inst(i_constraint, 4, (yyvsp[(1) - (4)].code), (yyvsp[(2) - (4)].code),
           code_new_inst(i_term_expr, 1, (yyvsp[(3) - (4)].code)),
           code_new_bits((yyvsp[(4) - (4)].bits)));
     ;}
    break;

  case 80:
#line 509 "src/mmlparse2.y"
    {
        (yyval.code) = code_new_inst(i_constraint, 4,
           code_new_inst(i_term_expr, 1, (yyvsp[(1) - (4)].code)),
           (yyvsp[(2) - (4)].code), (yyvsp[(3) - (4)].code), code_new_bits((yyvsp[(4) - (4)].bits)));
     ;}
    break;

  case 81:
#line 514 "src/mmlparse2.y"
    { 
        (yyval.code) = code_new_inst(i_constraint, 4,
           code_new_inst(i_term_expr, 1, (yyvsp[(1) - (4)].code)),
           (yyvsp[(2) - (4)].code),
           code_new_inst(i_term_expr, 1, (yyvsp[(3) - (4)].code)),
           code_new_bits((yyvsp[(4) - (4)].bits)));
     ;}
    break;

  case 82:
#line 521 "src/mmlparse2.y"
    {
        (yyval.code) = code_new_inst(i_rangeconst, 6, (yyvsp[(1) - (6)].code), (yyvsp[(3) - (6)].code), (yyvsp[(5) - (6)].code), (yyvsp[(2) - (6)].code),
           code_new_contype(CON_RHS), code_new_bits((yyvsp[(6) - (6)].bits))); 
     ;}
    break;

  case 83:
#line 525 "src/mmlparse2.y"
    {
        (yyval.code) = code_new_inst(i_rangeconst, 6, (yyvsp[(1) - (6)].code),
           code_new_inst(i_term_expr, 1, (yyvsp[(3) - (6)].code)), (yyvsp[(5) - (6)].code), (yyvsp[(2) - (6)].code),
           code_new_contype(CON_RHS), code_new_bits((yyvsp[(6) - (6)].bits))); 
     ;}
    break;

  case 84:
#line 530 "src/mmlparse2.y"
    {
        (yyval.code) = code_new_inst(i_rangeconst, 6, (yyvsp[(5) - (6)].code), (yyvsp[(3) - (6)].code), (yyvsp[(1) - (6)].code), (yyvsp[(2) - (6)].code),
           code_new_contype(CON_LHS), code_new_bits((yyvsp[(6) - (6)].bits))); 
     ;}
    break;

  case 85:
#line 534 "src/mmlparse2.y"
    {
        (yyval.code) = code_new_inst(i_rangeconst, 6, (yyvsp[(5) - (6)].code),
           code_new_inst(i_term_expr, 1, (yyvsp[(3) - (6)].code)),
           (yyvsp[(1) - (6)].code), (yyvsp[(2) - (6)].code),
           code_new_contype(CON_LHS), code_new_bits((yyvsp[(6) - (6)].bits))); 
     ;}
    break;

  case 86:
#line 540 "src/mmlparse2.y"
    {
         (yyval.code) = code_new_inst(i_vif_else, 8, (yyvsp[(2) - (12)].code), (yyvsp[(4) - (12)].code), (yyvsp[(5) - (12)].code), (yyvsp[(6) - (12)].code), (yyvsp[(8) - (12)].code), (yyvsp[(9) - (12)].code), (yyvsp[(10) - (12)].code), code_new_bits((yyvsp[(12) - (12)].bits)));
      ;}
    break;

  case 87:
#line 543 "src/mmlparse2.y"
    {
         (yyval.code) = code_new_inst(i_vif_else, 8, (yyvsp[(2) - (12)].code),
            code_new_inst(i_term_expr, 1, (yyvsp[(4) - (12)].code)),
            (yyvsp[(5) - (12)].code), (yyvsp[(6) - (12)].code), (yyvsp[(8) - (12)].code), (yyvsp[(9) - (12)].code), (yyvsp[(10) - (12)].code), code_new_bits((yyvsp[(12) - (12)].bits)));
      ;}
    break;

  case 88:
#line 548 "src/mmlparse2.y"
    { 
         (yyval.code) = code_new_inst(i_vif_else, 8, (yyvsp[(2) - (12)].code), (yyvsp[(4) - (12)].code), (yyvsp[(5) - (12)].code),
            code_new_inst(i_term_expr, 1, (yyvsp[(6) - (12)].code)),
            (yyvsp[(8) - (12)].code), (yyvsp[(9) - (12)].code), (yyvsp[(10) - (12)].code), code_new_bits((yyvsp[(12) - (12)].bits)));
      ;}
    break;

  case 89:
#line 553 "src/mmlparse2.y"
    { 
         (yyval.code) = code_new_inst(i_vif_else, 8, (yyvsp[(2) - (12)].code), (yyvsp[(4) - (12)].code), (yyvsp[(5) - (12)].code), (yyvsp[(6) - (12)].code),
            code_new_inst(i_term_expr, 1, (yyvsp[(8) - (12)].code)),
            (yyvsp[(9) - (12)].code), (yyvsp[(10) - (12)].code), code_new_bits((yyvsp[(12) - (12)].bits)));
      ;}
    break;

  case 90:
#line 558 "src/mmlparse2.y"
    { 
         (yyval.code) = code_new_inst(i_vif_else, 8, (yyvsp[(2) - (12)].code), (yyvsp[(4) - (12)].code), (yyvsp[(5) - (12)].code), (yyvsp[(6) - (12)].code), (yyvsp[(8) - (12)].code), (yyvsp[(9) - (12)].code),
            code_new_inst(i_term_expr, 1, (yyvsp[(10) - (12)].code)), code_new_bits((yyvsp[(12) - (12)].bits)));
      ;}
    break;

  case 91:
#line 562 "src/mmlparse2.y"
    { /* ??? This is an error */
         (yyval.code) = code_new_inst(i_vif_else, 8, (yyvsp[(2) - (12)].code),
            code_new_inst(i_term_expr, 1, (yyvsp[(4) - (12)].code)),
            (yyvsp[(5) - (12)].code),
            code_new_inst(i_term_expr, 1, (yyvsp[(6) - (12)].code)),
            (yyvsp[(8) - (12)].code), (yyvsp[(9) - (12)].code), (yyvsp[(10) - (12)].code), code_new_bits((yyvsp[(12) - (12)].bits)));
      ;}
    break;

  case 92:
#line 569 "src/mmlparse2.y"
    { 
         (yyval.code) = code_new_inst(i_vif_else, 8, (yyvsp[(2) - (12)].code),
            code_new_inst(i_term_expr, 1, (yyvsp[(4) - (12)].code)),
            (yyvsp[(5) - (12)].code), (yyvsp[(6) - (12)].code),
            code_new_inst(i_term_expr, 1, (yyvsp[(8) - (12)].code)),
            (yyvsp[(9) - (12)].code), (yyvsp[(10) - (12)].code), code_new_bits((yyvsp[(12) - (12)].bits)));
      ;}
    break;

  case 93:
#line 576 "src/mmlparse2.y"
    { 
         (yyval.code) = code_new_inst(i_vif_else, 8, (yyvsp[(2) - (12)].code),
            code_new_inst(i_term_expr, 1, (yyvsp[(4) - (12)].code)),
            (yyvsp[(5) - (12)].code), (yyvsp[(6) - (12)].code), (yyvsp[(8) - (12)].code), (yyvsp[(9) - (12)].code),
            code_new_inst(i_term_expr, 1, (yyvsp[(10) - (12)].code)), code_new_bits((yyvsp[(12) - (12)].bits)));
      ;}
    break;

  case 94:
#line 582 "src/mmlparse2.y"
    { 
         (yyval.code) = code_new_inst(i_vif_else, 8, (yyvsp[(2) - (12)].code), (yyvsp[(4) - (12)].code), (yyvsp[(5) - (12)].code),
            code_new_inst(i_term_expr, 1, (yyvsp[(6) - (12)].code)),
            code_new_inst(i_term_expr, 1, (yyvsp[(8) - (12)].code)),
            (yyvsp[(9) - (12)].code), (yyvsp[(10) - (12)].code), code_new_bits((yyvsp[(12) - (12)].bits)));
      ;}
    break;

  case 95:
#line 588 "src/mmlparse2.y"
    { 
         (yyval.code) = code_new_inst(i_vif_else, 8, (yyvsp[(2) - (12)].code), (yyvsp[(4) - (12)].code), (yyvsp[(5) - (12)].code),
            code_new_inst(i_term_expr, 1, (yyvsp[(6) - (12)].code)),
            (yyvsp[(8) - (12)].code), (yyvsp[(9) - (12)].code),
            code_new_inst(i_term_expr, 1, (yyvsp[(10) - (12)].code)), code_new_bits((yyvsp[(12) - (12)].bits)));
      ;}
    break;

  case 96:
#line 594 "src/mmlparse2.y"
    { /* ??? This is an error */
         (yyval.code) = code_new_inst(i_vif_else, 8, (yyvsp[(2) - (12)].code), (yyvsp[(4) - (12)].code), (yyvsp[(5) - (12)].code), (yyvsp[(6) - (12)].code),
            code_new_inst(i_term_expr, 1, (yyvsp[(8) - (12)].code)), (yyvsp[(9) - (12)].code),
            code_new_inst(i_term_expr, 1, (yyvsp[(10) - (12)].code)), code_new_bits((yyvsp[(12) - (12)].bits)));
      ;}
    break;

  case 97:
#line 599 "src/mmlparse2.y"
    { /* ??? This is an error */
         (yyval.code) = code_new_inst(i_vif_else, 8, (yyvsp[(2) - (12)].code),
            code_new_inst(i_term_expr, 1, (yyvsp[(4) - (12)].code)),
            (yyvsp[(5) - (12)].code),
            code_new_inst(i_term_expr, 1, (yyvsp[(6) - (12)].code)),
            code_new_inst(i_term_expr, 1, (yyvsp[(8) - (12)].code)),
            (yyvsp[(9) - (12)].code), (yyvsp[(10) - (12)].code), code_new_bits((yyvsp[(12) - (12)].bits)));
      ;}
    break;

  case 98:
#line 607 "src/mmlparse2.y"
    { /* ??? This is an error */
         (yyval.code) = code_new_inst(i_vif_else, 8, (yyvsp[(2) - (12)].code),
            code_new_inst(i_term_expr, 1, (yyvsp[(4) - (12)].code)),
            (yyvsp[(5) - (12)].code),
            code_new_inst(i_term_expr, 1, (yyvsp[(6) - (12)].code)),
            (yyvsp[(8) - (12)].code), (yyvsp[(9) - (12)].code),
            code_new_inst(i_term_expr, 1, (yyvsp[(10) - (12)].code)), 
            code_new_bits((yyvsp[(12) - (12)].bits)));
      ;}
    break;

  case 99:
#line 616 "src/mmlparse2.y"
    { /* ??? This is an error */
         (yyval.code) = code_new_inst(i_vif_else, 8, (yyvsp[(2) - (12)].code),
            code_new_inst(i_term_expr, 1, (yyvsp[(4) - (12)].code)),
            (yyvsp[(5) - (12)].code), (yyvsp[(6) - (12)].code),
            code_new_inst(i_term_expr, 1, (yyvsp[(8) - (12)].code)),
            (yyvsp[(9) - (12)].code),
            code_new_inst(i_term_expr, 1, (yyvsp[(10) - (12)].code)), 
            code_new_bits((yyvsp[(12) - (12)].bits)));
      ;}
    break;

  case 100:
#line 625 "src/mmlparse2.y"
    { /* ??? This is an error */
         (yyval.code) = code_new_inst(i_vif_else, 8, (yyvsp[(2) - (12)].code), (yyvsp[(4) - (12)].code), (yyvsp[(5) - (12)].code),
            code_new_inst(i_term_expr, 1, (yyvsp[(6) - (12)].code)),
            code_new_inst(i_term_expr, 1, (yyvsp[(8) - (12)].code)),
            (yyvsp[(9) - (12)].code),
            code_new_inst(i_term_expr, 1, (yyvsp[(10) - (12)].code)), 
            code_new_bits((yyvsp[(12) - (12)].bits)));
      ;}
    break;

  case 101:
#line 633 "src/mmlparse2.y"
    { /* ??? This is an error */
         (yyval.code) = code_new_inst(i_vif_else, 8, (yyvsp[(2) - (12)].code),
            code_new_inst(i_term_expr, 1, (yyvsp[(4) - (12)].code)),
            (yyvsp[(5) - (12)].code),
            code_new_inst(i_term_expr, 1, (yyvsp[(6) - (12)].code)),
            code_new_inst(i_term_expr, 1, (yyvsp[(8) - (12)].code)),
            (yyvsp[(9) - (12)].code),
            code_new_inst(i_term_expr, 1, (yyvsp[(10) - (12)].code)), 
            code_new_bits((yyvsp[(12) - (12)].bits)));
      ;}
    break;

  case 102:
#line 644 "src/mmlparse2.y"
    {
         (yyval.code) = code_new_inst(i_vif, 5, (yyvsp[(2) - (8)].code), (yyvsp[(4) - (8)].code), (yyvsp[(5) - (8)].code), (yyvsp[(6) - (8)].code), code_new_bits((yyvsp[(8) - (8)].bits)));
      ;}
    break;

  case 103:
#line 647 "src/mmlparse2.y"
    {
         (yyval.code) = code_new_inst(i_vif, 5, (yyvsp[(2) - (8)].code), 
            code_new_inst(i_term_expr, 1, (yyvsp[(4) - (8)].code)), (yyvsp[(5) - (8)].code), (yyvsp[(6) - (8)].code), code_new_bits((yyvsp[(8) - (8)].bits)));
      ;}
    break;

  case 104:
#line 651 "src/mmlparse2.y"
    {
         (yyval.code) = code_new_inst(i_vif, 5, (yyvsp[(2) - (8)].code), 
            (yyvsp[(4) - (8)].code), (yyvsp[(5) - (8)].code), code_new_inst(i_term_expr, 1, (yyvsp[(6) - (8)].code)), 
            code_new_bits((yyvsp[(8) - (8)].bits)));
      ;}
    break;

  case 105:
#line 656 "src/mmlparse2.y"
    { /* ??? This is an error */
         (yyval.code) = code_new_inst(i_vif, 5, (yyvsp[(2) - (8)].code),
            code_new_inst(i_term_expr, 1, (yyvsp[(4) - (8)].code)), (yyvsp[(5) - (8)].code), 
            code_new_inst(i_term_expr, 1, (yyvsp[(6) - (8)].code)), code_new_bits((yyvsp[(8) - (8)].bits)));
      ;}
    break;

  case 106:
#line 664 "src/mmlparse2.y"
    { (yyval.code) = code_new_inst(i_vbool_ne, 2, (yyvsp[(1) - (3)].code), (yyvsp[(3) - (3)].code)); ;}
    break;

  case 107:
#line 665 "src/mmlparse2.y"
    {
         (yyval.code) = code_new_inst(i_vbool_ne, 2, code_new_inst(i_term_expr, 1, (yyvsp[(1) - (3)].code)), (yyvsp[(3) - (3)].code));
      ;}
    break;

  case 108:
#line 668 "src/mmlparse2.y"
    {
         (yyval.code) = code_new_inst(i_vbool_ne, 2, (yyvsp[(1) - (3)].code), code_new_inst(i_term_expr, 1, (yyvsp[(3) - (3)].code)));
      ;}
    break;

  case 109:
#line 671 "src/mmlparse2.y"
    { (yyval.code) = code_new_inst(i_vbool_eq, 2, (yyvsp[(1) - (3)].code), (yyvsp[(3) - (3)].code)); ;}
    break;

  case 110:
#line 672 "src/mmlparse2.y"
    {
         (yyval.code) = code_new_inst(i_vbool_eq, 2, code_new_inst(i_term_expr, 1, (yyvsp[(1) - (3)].code)), (yyvsp[(3) - (3)].code));
      ;}
    break;

  case 111:
#line 675 "src/mmlparse2.y"
    {
         (yyval.code) = code_new_inst(i_vbool_eq, 2, (yyvsp[(1) - (3)].code), code_new_inst(i_term_expr, 1, (yyvsp[(3) - (3)].code)));
      ;}
    break;

  case 112:
#line 678 "src/mmlparse2.y"
    { (yyval.code) = code_new_inst(i_vbool_le, 2, (yyvsp[(1) - (3)].code), (yyvsp[(3) - (3)].code)); ;}
    break;

  case 113:
#line 679 "src/mmlparse2.y"
    {
         (yyval.code) = code_new_inst(i_vbool_le, 2, code_new_inst(i_term_expr, 1, (yyvsp[(1) - (3)].code)), (yyvsp[(3) - (3)].code));
      ;}
    break;

  case 114:
#line 682 "src/mmlparse2.y"
    {
         (yyval.code) = code_new_inst(i_vbool_le, 2, (yyvsp[(1) - (3)].code), code_new_inst(i_term_expr, 1, (yyvsp[(3) - (3)].code)));
      ;}
    break;

  case 115:
#line 685 "src/mmlparse2.y"
    { (yyval.code) = code_new_inst(i_vbool_ge, 2, (yyvsp[(1) - (3)].code), (yyvsp[(3) - (3)].code)); ;}
    break;

  case 116:
#line 686 "src/mmlparse2.y"
    {
         (yyval.code) = code_new_inst(i_vbool_ge, 2, code_new_inst(i_term_expr, 1, (yyvsp[(1) - (3)].code)), (yyvsp[(3) - (3)].code));
      ;}
    break;

  case 117:
#line 689 "src/mmlparse2.y"
    {
         (yyval.code) = code_new_inst(i_vbool_ge, 2, (yyvsp[(1) - (3)].code), code_new_inst(i_term_expr, 1, (yyvsp[(3) - (3)].code)));
      ;}
    break;

  case 118:
#line 692 "src/mmlparse2.y"
    { (yyval.code) = code_new_inst(i_vbool_lt, 2, (yyvsp[(1) - (3)].code), (yyvsp[(3) - (3)].code)); ;}
    break;

  case 119:
#line 693 "src/mmlparse2.y"
    {
         (yyval.code) = code_new_inst(i_vbool_lt, 2, code_new_inst(i_term_expr, 1, (yyvsp[(1) - (3)].code)), (yyvsp[(3) - (3)].code));
      ;}
    break;

  case 120:
#line 696 "src/mmlparse2.y"
    {
         (yyval.code) = code_new_inst(i_vbool_lt, 2, (yyvsp[(1) - (3)].code), code_new_inst(i_term_expr, 1, (yyvsp[(3) - (3)].code)));
      ;}
    break;

  case 121:
#line 699 "src/mmlparse2.y"
    { (yyval.code) = code_new_inst(i_vbool_gt, 2, (yyvsp[(1) - (3)].code), (yyvsp[(3) - (3)].code)); ;}
    break;

  case 122:
#line 700 "src/mmlparse2.y"
    {
         (yyval.code) = code_new_inst(i_vbool_gt, 2, code_new_inst(i_term_expr, 1, (yyvsp[(1) - (3)].code)), (yyvsp[(3) - (3)].code));
      ;}
    break;

  case 123:
#line 703 "src/mmlparse2.y"
    {
         (yyval.code) = code_new_inst(i_vbool_gt, 2, (yyvsp[(1) - (3)].code), code_new_inst(i_term_expr, 1, (yyvsp[(3) - (3)].code)));
      ;}
    break;

  case 124:
#line 706 "src/mmlparse2.y"
    { (yyval.code) = code_new_inst(i_vbool_and, 2, (yyvsp[(1) - (3)].code), (yyvsp[(3) - (3)].code)); ;}
    break;

  case 125:
#line 707 "src/mmlparse2.y"
    { (yyval.code) = code_new_inst(i_vbool_or,  2, (yyvsp[(1) - (3)].code), (yyvsp[(3) - (3)].code)); ;}
    break;

  case 126:
#line 708 "src/mmlparse2.y"
    { (yyval.code) = code_new_inst(i_vbool_xor, 2, (yyvsp[(1) - (3)].code), (yyvsp[(3) - (3)].code)); ;}
    break;

  case 127:
#line 709 "src/mmlparse2.y"
    { (yyval.code) = code_new_inst(i_vbool_not, 1, (yyvsp[(2) - (2)].code)); ;}
    break;

  case 128:
#line 710 "src/mmlparse2.y"
    { (yyval.code) = (yyvsp[(2) - (3)].code); ;}
    break;

  case 129:
#line 714 "src/mmlparse2.y"
    { (yyval.bits) = 0; ;}
    break;

  case 130:
#line 715 "src/mmlparse2.y"
    { (yyval.bits) = (yyvsp[(1) - (3)].bits) | (yyvsp[(3) - (3)].bits); ;}
    break;

  case 131:
#line 719 "src/mmlparse2.y"
    { (yyval.bits) = LP_FLAG_CON_SCALE; ;}
    break;

  case 132:
#line 720 "src/mmlparse2.y"
    { (yyval.bits) = LP_FLAG_CON_SEPAR; ;}
    break;

  case 133:
#line 721 "src/mmlparse2.y"
    { (yyval.bits) = LP_FLAG_CON_CHECK; ;}
    break;

  case 134:
#line 722 "src/mmlparse2.y"
    { (yyval.bits) = LP_FLAG_CON_INDIC; ;}
    break;

  case 135:
#line 726 "src/mmlparse2.y"
    { (yyval.code) = code_new_contype(CON_RHS); ;}
    break;

  case 136:
#line 727 "src/mmlparse2.y"
    { (yyval.code) = code_new_contype(CON_LHS); ;}
    break;

  case 137:
#line 728 "src/mmlparse2.y"
    { (yyval.code) = code_new_contype(CON_EQUAL); ;}
    break;

  case 138:
#line 732 "src/mmlparse2.y"
    { (yyval.code) = (yyvsp[(1) - (1)].code); ;}
    break;

  case 139:
#line 733 "src/mmlparse2.y"
    { (yyval.code) = code_new_inst(i_term_add, 2, (yyvsp[(1) - (3)].code), (yyvsp[(3) - (3)].code)); ;}
    break;

  case 140:
#line 734 "src/mmlparse2.y"
    { (yyval.code) = code_new_inst(i_term_sub, 2, (yyvsp[(1) - (3)].code), (yyvsp[(3) - (3)].code)); ;}
    break;

  case 141:
#line 735 "src/mmlparse2.y"
    { (yyval.code) = code_new_inst(i_term_const, 2, (yyvsp[(1) - (3)].code), (yyvsp[(3) - (3)].code)); ;}
    break;

  case 142:
#line 736 "src/mmlparse2.y"
    {
         (yyval.code) = code_new_inst(i_term_sub, 2, (yyvsp[(1) - (3)].code), code_new_inst(i_term_expr, 1, (yyvsp[(3) - (3)].code)));
      ;}
    break;

  case 143:
#line 739 "src/mmlparse2.y"
    { (yyval.code) = code_new_inst(i_term_const, 2, (yyvsp[(3) - (3)].code), (yyvsp[(1) - (3)].code)); ;}
    break;

  case 144:
#line 740 "src/mmlparse2.y"
    {
         (yyval.code) = code_new_inst(i_term_sub, 2,
            code_new_inst(i_term_expr, 1, (yyvsp[(1) - (3)].code)),
            (yyvsp[(3) - (3)].code));
      ;}
    break;

  case 145:
#line 748 "src/mmlparse2.y"
    { (yyval.code) = (yyvsp[(1) - (1)].code); ;}
    break;

  case 146:
#line 749 "src/mmlparse2.y"
    { (yyval.code) = code_new_inst(i_term_coeff, 2, (yyvsp[(1) - (3)].code), (yyvsp[(3) - (3)].code));  ;}
    break;

  case 147:
#line 750 "src/mmlparse2.y"
    {
         (yyval.code) = code_new_inst(i_term_coeff, 2, (yyvsp[(1) - (3)].code),
            code_new_inst(i_expr_div, 2, code_new_numb(numb_new_integer(1)), (yyvsp[(3) - (3)].code)));
      ;}
    break;

  case 148:
#line 754 "src/mmlparse2.y"
    { (yyval.code) = code_new_inst(i_term_coeff, 2, (yyvsp[(3) - (3)].code), (yyvsp[(1) - (3)].code)); ;}
    break;

  case 149:
#line 755 "src/mmlparse2.y"
    { (yyval.code) = code_new_inst(i_term_mul, 2, (yyvsp[(1) - (3)].code), (yyvsp[(3) - (3)].code)); ;}
    break;

  case 151:
#line 760 "src/mmlparse2.y"
    { (yyval.code) = (yyvsp[(2) - (2)].code); ;}
    break;

  case 152:
#line 761 "src/mmlparse2.y"
    { 
         (yyval.code) = code_new_inst(i_term_coeff, 2, (yyvsp[(2) - (2)].code), code_new_numb(numb_new_integer(-1)));
      ;}
    break;

  case 153:
#line 767 "src/mmlparse2.y"
    { (yyval.code) = (yyvsp[(1) - (1)].code); ;}
    break;

  case 154:
#line 768 "src/mmlparse2.y"
    { 
         (yyval.code) = code_new_inst(i_term_power, 2, (yyvsp[(1) - (3)].code), (yyvsp[(3) - (3)].code));
      ;}
    break;

  case 155:
#line 771 "src/mmlparse2.y"
    {
         (yyval.code) = code_new_inst(i_term_sum, 2, (yyvsp[(2) - (4)].code), (yyvsp[(4) - (4)].code));
      ;}
    break;

  case 156:
#line 777 "src/mmlparse2.y"
    {
         (yyval.code) = code_new_inst(i_symbol_deref, 2, code_new_symbol((yyvsp[(1) - (2)].sym)), (yyvsp[(2) - (2)].code));
      ;}
    break;

  case 157:
#line 780 "src/mmlparse2.y"
    { (yyval.code) = code_new_inst(i_vabs, 1, (yyvsp[(3) - (4)].code)); ;}
    break;

  case 158:
#line 781 "src/mmlparse2.y"
    { (yyval.code) = code_new_inst(i_vexpr_fun, 2, code_new_numb(numb_new_integer(-2)), (yyvsp[(3) - (4)].code)); ;}
    break;

  case 159:
#line 782 "src/mmlparse2.y"
    { (yyval.code) = code_new_inst(i_vexpr_fun, 2, code_new_numb(numb_new_integer(3)), (yyvsp[(3) - (4)].code)); ;}
    break;

  case 160:
#line 783 "src/mmlparse2.y"
    { (yyval.code) = code_new_inst(i_vexpr_fun, 2, code_new_numb(numb_new_integer(4)), (yyvsp[(3) - (4)].code)); ;}
    break;

  case 161:
#line 784 "src/mmlparse2.y"
    { (yyval.code) = code_new_inst(i_vexpr_fun, 2, code_new_numb(numb_new_integer(5)), (yyvsp[(3) - (4)].code)); ;}
    break;

  case 162:
#line 785 "src/mmlparse2.y"
    { (yyval.code) = code_new_inst(i_vexpr_fun, 2, code_new_numb(numb_new_integer(6)), (yyvsp[(3) - (4)].code)); ;}
    break;

  case 163:
#line 786 "src/mmlparse2.y"
    { (yyval.code) = code_new_inst(i_vexpr_fun, 2, code_new_numb(numb_new_integer(7)), (yyvsp[(3) - (4)].code)); ;}
    break;

  case 164:
#line 787 "src/mmlparse2.y"
    { (yyval.code) = code_new_inst(i_vexpr_fun, 2, code_new_numb(numb_new_integer(8)), (yyvsp[(3) - (4)].code)); ;}
    break;

  case 165:
#line 788 "src/mmlparse2.y"
    { (yyval.code) = code_new_inst(i_vexpr_fun, 2, code_new_numb(numb_new_integer(9)), (yyvsp[(3) - (4)].code)); ;}
    break;

  case 166:
#line 789 "src/mmlparse2.y"
    { (yyval.code) = code_new_inst(i_vexpr_fun, 2, code_new_numb(numb_new_integer(10)), (yyvsp[(3) - (4)].code)); ;}
    break;

  case 167:
#line 790 "src/mmlparse2.y"
    {
         (yyval.code) = code_new_inst(i_vexpr_fun, 3, code_new_numb(numb_new_integer(11)), (yyvsp[(3) - (6)].code), (yyvsp[(5) - (6)].code));
      ;}
    break;

  case 168:
#line 793 "src/mmlparse2.y"
    {
         (yyval.code) = code_new_inst(i_vexpr_fun, 3, code_new_numb(numb_new_integer(12)), (yyvsp[(3) - (6)].code), (yyvsp[(5) - (6)].code));
      ;}
    break;

  case 169:
#line 796 "src/mmlparse2.y"
    {
         (yyval.code) = code_new_inst(i_expr_if_else, 3, (yyvsp[(2) - (7)].code), (yyvsp[(4) - (7)].code), (yyvsp[(6) - (7)].code));
      ;}
    break;

  case 170:
#line 799 "src/mmlparse2.y"
    { (yyval.code) = (yyvsp[(2) - (3)].code); ;}
    break;

  case 171:
#line 807 "src/mmlparse2.y"
    {
        (yyval.code) = code_new_inst(i_sos, 2, code_new_name((yyvsp[(2) - (5)].name)), (yyvsp[(4) - (5)].code));
     ;}
    break;

  case 172:
#line 813 "src/mmlparse2.y"
    {
        (yyval.code) = code_new_inst(i_soset, 3, (yyvsp[(4) - (4)].code), (yyvsp[(1) - (4)].code), (yyvsp[(2) - (4)].code));
     ;}
    break;

  case 173:
#line 816 "src/mmlparse2.y"
    {
         (yyval.code) = code_new_inst(i_forall, 2, (yyvsp[(2) - (4)].code), (yyvsp[(4) - (4)].code));
      ;}
    break;

  case 174:
#line 822 "src/mmlparse2.y"
    { (yyval.code) = code_new_numb(numb_new_integer(1)); ;}
    break;

  case 175:
#line 823 "src/mmlparse2.y"
    { (yyval.code) = code_new_numb(numb_new_integer(1)); ;}
    break;

  case 176:
#line 824 "src/mmlparse2.y"
    { (yyval.code) = code_new_numb(numb_new_integer(2)); ;}
    break;

  case 177:
#line 832 "src/mmlparse2.y"
    { (yyval.code) = (yyvsp[(2) - (3)].code); ;}
    break;

  case 178:
#line 836 "src/mmlparse2.y"
    { (yyval.code) = code_new_inst(i_print, 1, (yyvsp[(2) - (2)].code)); ;}
    break;

  case 179:
#line 837 "src/mmlparse2.y"
    { (yyval.code) = code_new_inst(i_print, 1, (yyvsp[(2) - (2)].code)); ;}
    break;

  case 180:
#line 838 "src/mmlparse2.y"
    { (yyval.code) = code_new_inst(i_print, 1, (yyvsp[(2) - (2)].code)); ;}
    break;

  case 181:
#line 839 "src/mmlparse2.y"
    { (yyval.code) = code_new_inst(i_print, 1, code_new_symbol((yyvsp[(2) - (2)].sym))); ;}
    break;

  case 182:
#line 840 "src/mmlparse2.y"
    { (yyval.code) = code_new_inst(i_check, 1, (yyvsp[(2) - (2)].code)); ;}
    break;

  case 183:
#line 841 "src/mmlparse2.y"
    {
        (yyval.code) = code_new_inst(i_forall, 2, (yyvsp[(2) - (4)].code), (yyvsp[(4) - (4)].code));
     ;}
    break;

  case 184:
#line 851 "src/mmlparse2.y"
    { (yyval.code) = (yyvsp[(1) - (1)].code); ;}
    break;

  case 185:
#line 852 "src/mmlparse2.y"
    {
         (yyval.code) = code_new_inst(i_idxset_new, 3,
            code_new_inst(i_tuple_empty, 0), (yyvsp[(1) - (1)].code), code_new_inst(i_bool_true, 0));
      ;}
    break;

  case 186:
#line 859 "src/mmlparse2.y"
    {
         (yyval.code) = code_new_inst(i_idxset_new, 3, (yyvsp[(1) - (5)].code), (yyvsp[(3) - (5)].code), (yyvsp[(5) - (5)].code));
      ;}
    break;

  case 187:
#line 862 "src/mmlparse2.y"
    {
         (yyval.code) = code_new_inst(i_idxset_new, 3, (yyvsp[(1) - (3)].code), (yyvsp[(3) - (3)].code), code_new_inst(i_bool_true, 0));
      ;}
    break;

  case 189:
#line 869 "src/mmlparse2.y"
    { (yyval.code) = code_new_inst(i_set_union, 2, (yyvsp[(1) - (3)].code), (yyvsp[(3) - (3)].code)); ;}
    break;

  case 190:
#line 870 "src/mmlparse2.y"
    {
         (yyval.code) = code_new_inst(i_set_union, 2, (yyvsp[(1) - (3)].code), (yyvsp[(3) - (3)].code));
      ;}
    break;

  case 191:
#line 873 "src/mmlparse2.y"
    { (yyval.code) = code_new_inst(i_set_sdiff, 2, (yyvsp[(1) - (3)].code), (yyvsp[(3) - (3)].code)); ;}
    break;

  case 192:
#line 874 "src/mmlparse2.y"
    {
         (yyval.code) = code_new_inst(i_set_minus, 2, (yyvsp[(1) - (3)].code), (yyvsp[(3) - (3)].code));
      ;}
    break;

  case 193:
#line 877 "src/mmlparse2.y"
    { (yyval.code) = code_new_inst(i_set_minus, 2, (yyvsp[(1) - (3)].code), (yyvsp[(3) - (3)].code)); ;}
    break;

  case 194:
#line 878 "src/mmlparse2.y"
    { (yyval.code) = code_new_inst(i_set_inter, 2, (yyvsp[(1) - (3)].code), (yyvsp[(3) - (3)].code)); ;}
    break;

  case 196:
#line 882 "src/mmlparse2.y"
    { (yyval.code) = code_new_inst(i_set_union2, 2, (yyvsp[(2) - (4)].code), (yyvsp[(4) - (4)].code)); ;}
    break;

  case 198:
#line 887 "src/mmlparse2.y"
    { (yyval.code) = code_new_inst(i_set_cross, 2, (yyvsp[(1) - (3)].code), (yyvsp[(3) - (3)].code)); ;}
    break;

  case 199:
#line 888 "src/mmlparse2.y"
    {
         (yyval.code) = code_new_inst(i_set_cross, 2, (yyvsp[(1) - (3)].code), (yyvsp[(3) - (3)].code));
      ;}
    break;

  case 200:
#line 891 "src/mmlparse2.y"
    { (yyval.code) = code_new_inst(i_set_inter2, 2, (yyvsp[(2) - (4)].code), (yyvsp[(4) - (4)].code)); ;}
    break;

  case 201:
#line 895 "src/mmlparse2.y"
    {
         (yyval.code) = code_new_inst(i_symbol_deref, 2, code_new_symbol((yyvsp[(1) - (2)].sym)), (yyvsp[(2) - (2)].code));
      ;}
    break;

  case 202:
#line 898 "src/mmlparse2.y"
    {
         (yyval.code) = code_new_inst(i_define_deref, 2,
            code_new_define((yyvsp[(1) - (4)].def)),
            code_new_inst(i_tuple_new, 1, (yyvsp[(3) - (4)].code)));
      ;}
    break;

  case 203:
#line 903 "src/mmlparse2.y"
    { (yyval.code) = code_new_inst(i_set_empty, 1, code_new_size(0)); ;}
    break;

  case 204:
#line 904 "src/mmlparse2.y"
    {
         (yyval.code) = code_new_inst(i_set_range2, 3, (yyvsp[(2) - (7)].code), (yyvsp[(4) - (7)].code), (yyvsp[(6) - (7)].code));
      ;}
    break;

  case 205:
#line 907 "src/mmlparse2.y"
    {
         (yyval.code) = code_new_inst(i_set_range2, 3, (yyvsp[(2) - (5)].code), (yyvsp[(4) - (5)].code), code_new_numb(numb_new_integer(1)));
      ;}
    break;

  case 206:
#line 910 "src/mmlparse2.y"
    {
         (yyval.code) = code_new_inst(i_set_range, 3, (yyvsp[(2) - (7)].code), (yyvsp[(4) - (7)].code), (yyvsp[(6) - (7)].code));
      ;}
    break;

  case 207:
#line 913 "src/mmlparse2.y"
    {
         (yyval.code) = code_new_inst(i_set_range, 3, (yyvsp[(2) - (5)].code), (yyvsp[(4) - (5)].code), code_new_numb(numb_new_integer(1)));
      ;}
    break;

  case 208:
#line 916 "src/mmlparse2.y"
    {
         (yyval.code) = code_new_inst(i_set_argmin, 3, code_new_numb(numb_new_integer(1)), (yyvsp[(2) - (4)].code), (yyvsp[(4) - (4)].code));
      ;}
    break;

  case 209:
#line 919 "src/mmlparse2.y"
    {
         (yyval.code) = code_new_inst(i_set_argmin, 3, (yyvsp[(3) - (7)].code), (yyvsp[(5) - (7)].code), (yyvsp[(7) - (7)].code));
      ;}
    break;

  case 210:
#line 922 "src/mmlparse2.y"
    {
         (yyval.code) = code_new_inst(i_set_argmax, 3, code_new_numb(numb_new_integer(1)), (yyvsp[(2) - (4)].code), (yyvsp[(4) - (4)].code));
      ;}
    break;

  case 211:
#line 925 "src/mmlparse2.y"
    {
         (yyval.code) = code_new_inst(i_set_argmax, 3, (yyvsp[(3) - (7)].code), (yyvsp[(5) - (7)].code), (yyvsp[(7) - (7)].code));
      ;}
    break;

  case 212:
#line 928 "src/mmlparse2.y"
    { (yyval.code) = (yyvsp[(2) - (3)].code); ;}
    break;

  case 213:
#line 929 "src/mmlparse2.y"
    { (yyval.code) = code_new_inst(i_set_new_tuple, 1, (yyvsp[(2) - (3)].code)); ;}
    break;

  case 214:
#line 930 "src/mmlparse2.y"
    { (yyval.code) = code_new_inst(i_set_new_elem, 1, (yyvsp[(2) - (3)].code)); ;}
    break;

  case 215:
#line 931 "src/mmlparse2.y"
    { (yyval.code) = code_new_inst(i_set_idxset, 1, (yyvsp[(2) - (3)].code)); ;}
    break;

  case 216:
#line 932 "src/mmlparse2.y"
    { (yyval.code) = code_new_inst(i_set_expr, 2, (yyvsp[(2) - (5)].code), (yyvsp[(4) - (5)].code)); ;}
    break;

  case 217:
#line 933 "src/mmlparse2.y"
    { (yyval.code) = code_new_inst(i_set_expr, 2, (yyvsp[(2) - (5)].code), (yyvsp[(4) - (5)].code)); ;}
    break;

  case 218:
#line 934 "src/mmlparse2.y"
    {
         (yyval.code) = code_new_inst(i_set_proj, 2, (yyvsp[(3) - (6)].code), (yyvsp[(5) - (6)].code));
       ;}
    break;

  case 219:
#line 937 "src/mmlparse2.y"
    {
          (yyval.code) = code_new_inst(i_set_indexset, 1, code_new_symbol((yyvsp[(3) - (4)].sym)));
       ;}
    break;

  case 220:
#line 940 "src/mmlparse2.y"
    {
         (yyval.code) = code_new_inst(i_expr_if_else, 3, (yyvsp[(2) - (7)].code), (yyvsp[(4) - (7)].code), (yyvsp[(6) - (7)].code));
      ;}
    break;

  case 221:
#line 946 "src/mmlparse2.y"
    { (yyval.code) = code_new_inst(i_read_new, 2, (yyvsp[(2) - (4)].code), (yyvsp[(4) - (4)].code)); ;}
    break;

  case 222:
#line 947 "src/mmlparse2.y"
    { (yyval.code) = code_new_inst(i_read_param, 2, (yyvsp[(1) - (2)].code), (yyvsp[(2) - (2)].code)); ;}
    break;

  case 223:
#line 951 "src/mmlparse2.y"
    { (yyval.code) = code_new_inst(i_read_skip, 1, (yyvsp[(2) - (2)].code)); ;}
    break;

  case 224:
#line 952 "src/mmlparse2.y"
    { (yyval.code) = code_new_inst(i_read_use, 1, (yyvsp[(2) - (2)].code)); ;}
    break;

  case 225:
#line 953 "src/mmlparse2.y"
    { (yyval.code) = code_new_inst(i_read_comment, 1, (yyvsp[(2) - (2)].code)); ;}
    break;

  case 226:
#line 954 "src/mmlparse2.y"
    { (yyval.code) = code_new_inst(i_read_match, 1, (yyvsp[(2) - (2)].code)); ;}
    break;

  case 227:
#line 958 "src/mmlparse2.y"
    {
         (yyval.code) = code_new_inst(i_tuple_list_new, 1, (yyvsp[(1) - (1)].code));
      ;}
    break;

  case 228:
#line 961 "src/mmlparse2.y"
    {
         (yyval.code) = code_new_inst(i_tuple_list_add, 2, (yyvsp[(1) - (3)].code), (yyvsp[(3) - (3)].code));
      ;}
    break;

  case 229:
#line 964 "src/mmlparse2.y"
    { (yyval.code) = code_new_inst(i_read, 1, (yyvsp[(1) - (1)].code)); ;}
    break;

  case 230:
#line 968 "src/mmlparse2.y"
    { (yyval.code) = code_new_inst(i_bool_eq, 2, (yyvsp[(1) - (3)].code), (yyvsp[(3) - (3)].code)); ;}
    break;

  case 231:
#line 969 "src/mmlparse2.y"
    { (yyval.code) = code_new_inst(i_bool_ne, 2, (yyvsp[(1) - (3)].code), (yyvsp[(3) - (3)].code)); ;}
    break;

  case 232:
#line 970 "src/mmlparse2.y"
    { (yyval.code) = code_new_inst(i_bool_gt, 2, (yyvsp[(1) - (3)].code), (yyvsp[(3) - (3)].code)); ;}
    break;

  case 233:
#line 971 "src/mmlparse2.y"
    { (yyval.code) = code_new_inst(i_bool_ge, 2, (yyvsp[(1) - (3)].code), (yyvsp[(3) - (3)].code)); ;}
    break;

  case 234:
#line 972 "src/mmlparse2.y"
    { (yyval.code) = code_new_inst(i_bool_lt, 2, (yyvsp[(1) - (3)].code), (yyvsp[(3) - (3)].code)); ;}
    break;

  case 235:
#line 973 "src/mmlparse2.y"
    { (yyval.code) = code_new_inst(i_bool_le, 2, (yyvsp[(1) - (3)].code), (yyvsp[(3) - (3)].code)); ;}
    break;

  case 236:
#line 974 "src/mmlparse2.y"
    { (yyval.code) = code_new_inst(i_bool_seq, 2, (yyvsp[(1) - (3)].code), (yyvsp[(3) - (3)].code)); ;}
    break;

  case 237:
#line 975 "src/mmlparse2.y"
    { (yyval.code) = code_new_inst(i_bool_sneq, 2, (yyvsp[(1) - (3)].code), (yyvsp[(3) - (3)].code)); ;}
    break;

  case 238:
#line 976 "src/mmlparse2.y"
    { (yyval.code) = code_new_inst(i_bool_subs, 2, (yyvsp[(3) - (3)].code), (yyvsp[(1) - (3)].code)); ;}
    break;

  case 239:
#line 977 "src/mmlparse2.y"
    { (yyval.code) = code_new_inst(i_bool_sseq, 2, (yyvsp[(3) - (3)].code), (yyvsp[(1) - (3)].code)); ;}
    break;

  case 240:
#line 978 "src/mmlparse2.y"
    { (yyval.code) = code_new_inst(i_bool_subs, 2, (yyvsp[(1) - (3)].code), (yyvsp[(3) - (3)].code)); ;}
    break;

  case 241:
#line 979 "src/mmlparse2.y"
    { (yyval.code) = code_new_inst(i_bool_sseq, 2, (yyvsp[(1) - (3)].code), (yyvsp[(3) - (3)].code)); ;}
    break;

  case 242:
#line 980 "src/mmlparse2.y"
    { (yyval.code) = code_new_inst(i_bool_and, 2, (yyvsp[(1) - (3)].code), (yyvsp[(3) - (3)].code)); ;}
    break;

  case 243:
#line 981 "src/mmlparse2.y"
    { (yyval.code) = code_new_inst(i_bool_or,  2, (yyvsp[(1) - (3)].code), (yyvsp[(3) - (3)].code)); ;}
    break;

  case 244:
#line 982 "src/mmlparse2.y"
    { (yyval.code) = code_new_inst(i_bool_xor, 2, (yyvsp[(1) - (3)].code), (yyvsp[(3) - (3)].code)); ;}
    break;

  case 245:
#line 983 "src/mmlparse2.y"
    { (yyval.code) = code_new_inst(i_bool_not, 1, (yyvsp[(2) - (2)].code)); ;}
    break;

  case 246:
#line 984 "src/mmlparse2.y"
    { (yyval.code) = (yyvsp[(2) - (3)].code); ;}
    break;

  case 247:
#line 985 "src/mmlparse2.y"
    { (yyval.code) = code_new_inst(i_bool_is_elem, 2, (yyvsp[(1) - (3)].code), (yyvsp[(3) - (3)].code)); ;}
    break;

  case 248:
#line 986 "src/mmlparse2.y"
    { (yyval.code) = code_new_inst(i_bool_exists, 1, (yyvsp[(3) - (4)].code)); ;}
    break;

  case 249:
#line 987 "src/mmlparse2.y"
    {
         (yyval.code) = code_new_inst(i_define_deref, 2,
            code_new_define((yyvsp[(1) - (4)].def)),
            code_new_inst(i_tuple_new, 1, (yyvsp[(3) - (4)].code)));
      ;}
    break;

  case 250:
#line 992 "src/mmlparse2.y"
    {
        (yyval.code) = code_new_inst(i_expr_if_else, 3, (yyvsp[(2) - (7)].code), (yyvsp[(4) - (7)].code), (yyvsp[(6) - (7)].code));
     ;}
    break;

  case 251:
#line 998 "src/mmlparse2.y"
    { (yyval.code) = code_new_inst(i_tuple_empty, 0); ;}
    break;

  case 252:
#line 999 "src/mmlparse2.y"
    { (yyval.code) = code_new_inst(i_tuple_new, 1, (yyvsp[(2) - (3)].code));  ;}
    break;

  case 253:
#line 1003 "src/mmlparse2.y"
    {
         (yyval.code) = code_new_inst(i_tuple_empty, 0);
      ;}
    break;

  case 254:
#line 1006 "src/mmlparse2.y"
    {
         (yyval.code) = code_new_inst(i_tuple_new, 1, (yyvsp[(2) - (3)].code));
      ;}
    break;

  case 255:
#line 1012 "src/mmlparse2.y"
    {
         (yyval.code) = code_new_inst(i_elem_list_new, 1, (yyvsp[(1) - (1)].code));
      ;}
    break;

  case 256:
#line 1015 "src/mmlparse2.y"
    {
         (yyval.code) = code_new_inst(i_elem_list_add, 2, (yyvsp[(1) - (3)].code), (yyvsp[(3) - (3)].code));
      ;}
    break;

  case 257:
#line 1021 "src/mmlparse2.y"
    { (yyval.code) = (yyvsp[(1) - (1)].code); ;}
    break;

  case 258:
#line 1022 "src/mmlparse2.y"
    { (yyval.code) = code_new_inst(i_expr_add, 2, (yyvsp[(1) - (3)].code), (yyvsp[(3) - (3)].code)); ;}
    break;

  case 259:
#line 1023 "src/mmlparse2.y"
    { (yyval.code) = code_new_inst(i_expr_sub, 2, (yyvsp[(1) - (3)].code), (yyvsp[(3) - (3)].code)); ;}
    break;

  case 260:
#line 1027 "src/mmlparse2.y"
    { (yyval.code) = (yyvsp[(1) - (1)].code); ;}
    break;

  case 261:
#line 1028 "src/mmlparse2.y"
    { (yyval.code) = code_new_inst(i_expr_mul, 2, (yyvsp[(1) - (3)].code), (yyvsp[(3) - (3)].code)); ;}
    break;

  case 262:
#line 1029 "src/mmlparse2.y"
    { (yyval.code) = code_new_inst(i_expr_div, 2, (yyvsp[(1) - (3)].code), (yyvsp[(3) - (3)].code)); ;}
    break;

  case 263:
#line 1030 "src/mmlparse2.y"
    { (yyval.code) = code_new_inst(i_expr_mod, 2, (yyvsp[(1) - (3)].code), (yyvsp[(3) - (3)].code)); ;}
    break;

  case 264:
#line 1031 "src/mmlparse2.y"
    { (yyval.code) = code_new_inst(i_expr_intdiv, 2, (yyvsp[(1) - (3)].code), (yyvsp[(3) - (3)].code)); ;}
    break;

  case 265:
#line 1032 "src/mmlparse2.y"
    {
         (yyval.code) = code_new_inst(i_expr_prod, 2, (yyvsp[(2) - (4)].code), (yyvsp[(4) - (4)].code));
      ;}
    break;

  case 267:
#line 1039 "src/mmlparse2.y"
    { (yyval.code) = (yyvsp[(2) - (2)].code); ;}
    break;

  case 268:
#line 1040 "src/mmlparse2.y"
    { (yyval.code) = code_new_inst(i_expr_neg, 1, (yyvsp[(2) - (2)].code)); ;}
    break;

  case 270:
#line 1045 "src/mmlparse2.y"
    { (yyval.code) = code_new_inst(i_expr_pow, 2, (yyvsp[(1) - (3)].code), (yyvsp[(3) - (3)].code)); ;}
    break;

  case 271:
#line 1046 "src/mmlparse2.y"
    {
         (yyval.code) = code_new_inst(i_expr_sum, 2, (yyvsp[(2) - (4)].code), (yyvsp[(4) - (4)].code));
      ;}
    break;

  case 272:
#line 1049 "src/mmlparse2.y"
    {
         (yyval.code) = code_new_inst(i_expr_min, 2, (yyvsp[(2) - (4)].code), (yyvsp[(4) - (4)].code));
      ;}
    break;

  case 273:
#line 1052 "src/mmlparse2.y"
    {
         (yyval.code) = code_new_inst(i_expr_max, 2, (yyvsp[(2) - (4)].code), (yyvsp[(4) - (4)].code));
      ;}
    break;

  case 274:
#line 1055 "src/mmlparse2.y"
    {
         (yyval.code) = code_new_inst(i_expr_sglmin, 1, (yyvsp[(3) - (4)].code));
         ;}
    break;

  case 275:
#line 1058 "src/mmlparse2.y"
    {
         (yyval.code) = code_new_inst(i_expr_sglmax, 1, (yyvsp[(3) - (4)].code));
      ;}
    break;

  case 276:
#line 1064 "src/mmlparse2.y"
    { (yyval.code) = code_new_numb((yyvsp[(1) - (1)].numb)); ;}
    break;

  case 277:
#line 1065 "src/mmlparse2.y"
    { (yyval.code) = code_new_strg((yyvsp[(1) - (1)].strg));  ;}
    break;

  case 278:
#line 1066 "src/mmlparse2.y"
    {
         (yyval.code) = code_new_inst(i_local_deref, 1, code_new_name((yyvsp[(1) - (1)].name)));
      ;}
    break;

  case 279:
#line 1069 "src/mmlparse2.y"
    { 
         (yyval.code) = code_new_inst(i_symbol_deref, 2, code_new_symbol((yyvsp[(1) - (2)].sym)), (yyvsp[(2) - (2)].code));
      ;}
    break;

  case 280:
#line 1072 "src/mmlparse2.y"
    { 
         (yyval.code) = code_new_inst(i_symbol_deref, 2, code_new_symbol((yyvsp[(1) - (2)].sym)), (yyvsp[(2) - (2)].code));
      ;}
    break;

  case 281:
#line 1075 "src/mmlparse2.y"
    {
         (yyval.code) = code_new_inst(i_define_deref, 2,
            code_new_define((yyvsp[(1) - (4)].def)),
            code_new_inst(i_tuple_new, 1, (yyvsp[(3) - (4)].code)));
      ;}
    break;

  case 282:
#line 1080 "src/mmlparse2.y"
    {
         (yyval.code) = code_new_inst(i_define_deref, 2,
            code_new_define((yyvsp[(1) - (4)].def)),
            code_new_inst(i_tuple_new, 1, (yyvsp[(3) - (4)].code)));
      ;}
    break;

  case 283:
#line 1085 "src/mmlparse2.y"
    { (yyval.code) = code_new_inst(i_expr_fac, 1, (yyvsp[(1) - (2)].code)); ;}
    break;

  case 284:
#line 1086 "src/mmlparse2.y"
    { (yyval.code) = code_new_inst(i_expr_card, 1, (yyvsp[(3) - (4)].code)); ;}
    break;

  case 285:
#line 1087 "src/mmlparse2.y"
    { (yyval.code) = code_new_inst(i_expr_abs, 1, (yyvsp[(3) - (4)].code)); ;}
    break;

  case 286:
#line 1088 "src/mmlparse2.y"
    { (yyval.code) = code_new_inst(i_expr_sgn, 1, (yyvsp[(3) - (4)].code)); ;}
    break;

  case 287:
#line 1089 "src/mmlparse2.y"
    { (yyval.code) = code_new_inst(i_expr_round, 1, (yyvsp[(3) - (4)].code)); ;}
    break;

  case 288:
#line 1090 "src/mmlparse2.y"
    { (yyval.code) = code_new_inst(i_expr_floor, 1, (yyvsp[(3) - (4)].code)); ;}
    break;

  case 289:
#line 1091 "src/mmlparse2.y"
    { (yyval.code) = code_new_inst(i_expr_ceil, 1, (yyvsp[(3) - (4)].code)); ;}
    break;

  case 290:
#line 1092 "src/mmlparse2.y"
    { (yyval.code) = code_new_inst(i_expr_log, 1, (yyvsp[(3) - (4)].code)); ;}
    break;

  case 291:
#line 1093 "src/mmlparse2.y"
    { (yyval.code) = code_new_inst(i_expr_ln, 1, (yyvsp[(3) - (4)].code)); ;}
    break;

  case 292:
#line 1094 "src/mmlparse2.y"
    { (yyval.code) = code_new_inst(i_expr_exp, 1, (yyvsp[(3) - (4)].code)); ;}
    break;

  case 293:
#line 1095 "src/mmlparse2.y"
    { (yyval.code) = code_new_inst(i_expr_sqrt, 1, (yyvsp[(3) - (4)].code)); ;}
    break;

  case 294:
#line 1097 "src/mmlparse2.y"
    { (yyval.code) = (yyvsp[(2) - (3)].code); ;}
    break;

  case 295:
#line 1098 "src/mmlparse2.y"
    { (yyval.code) = code_new_inst(i_expr_length, 1, (yyvsp[(3) - (4)].code)); ;}
    break;

  case 296:
#line 1099 "src/mmlparse2.y"
    {
         (yyval.code) = code_new_inst(i_expr_substr, 3, (yyvsp[(3) - (8)].code), (yyvsp[(5) - (8)].code), (yyvsp[(7) - (8)].code));
      ;}
    break;

  case 297:
#line 1102 "src/mmlparse2.y"
    {
         (yyval.code) = code_new_inst(i_expr_rand, 2, (yyvsp[(3) - (6)].code), (yyvsp[(5) - (6)].code));
      ;}
    break;

  case 298:
#line 1105 "src/mmlparse2.y"
    {
         (yyval.code) = code_new_inst(i_expr_if_else, 3, (yyvsp[(2) - (7)].code), (yyvsp[(4) - (7)].code), (yyvsp[(6) - (7)].code));
      ;}
    break;

  case 299:
#line 1108 "src/mmlparse2.y"
    {
         (yyval.code) = code_new_inst(i_expr_ord, 3, (yyvsp[(3) - (8)].code), (yyvsp[(5) - (8)].code), (yyvsp[(7) - (8)].code));
      ;}
    break;

  case 300:
#line 1111 "src/mmlparse2.y"
    {
         (yyval.code) = code_new_inst(i_expr_min2, 1, (yyvsp[(3) - (4)].code));
      ;}
    break;

  case 301:
#line 1114 "src/mmlparse2.y"
    {
         (yyval.code) = code_new_inst(i_expr_max2, 1, (yyvsp[(3) - (4)].code));
      ;}
    break;


/* Line 1267 of yacc.c.  */
#line 4733 "src/mmlparse2.c"
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
      /* If just tried and failed to reuse look-ahead token after an
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

  /* Else will try to reuse look-ahead token after shifting the error
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

  if (yyn == YYFINAL)
    YYACCEPT;

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

#ifndef yyoverflow
/*-------------------------------------------------.
| yyexhaustedlab -- memory exhaustion comes here.  |
`-------------------------------------------------*/
yyexhaustedlab:
  yyerror (YY_("memory exhausted"));
  yyresult = 2;
  /* Fall through.  */
#endif

yyreturn:
  if (yychar != YYEOF && yychar != YYEMPTY)
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



