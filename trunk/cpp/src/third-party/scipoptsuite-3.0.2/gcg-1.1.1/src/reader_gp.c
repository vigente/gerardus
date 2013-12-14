/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* Copyright (C) 2010-2013 Operations Research, RWTH Aachen University       */
/*                         Zuse Institute Berlin (ZIB)                       */
/*                                                                           */
/* This program is free software; you can redistribute it and/or             */
/* modify it under the terms of the GNU Lesser General Public License        */
/* as published by the Free Software Foundation; either version 3            */
/* of the License, or (at your option) any later version.                    */
/*                                                                           */
/* This program is distributed in the hope that it will be useful,           */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of            */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             */
/* GNU Lesser General Public License for more details.                       */
/*                                                                           */
/* You should have received a copy of the GNU Lesser General Public License  */
/* along with this program; if not, write to the Free Software               */
/* Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA.*/
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   reader_gp.c
 * @brief  GP file reader writing gnuplot files
 * @author Martin Bergner
 * @todo change output file type based on parameter
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "reader_gp.h"
#include "scip_misc.h"
#include "struct_decomp.h"
#include "cons_decomp.h"

#define READER_NAME             "gpreader"
#define READER_DESC             "gnuplot file writer for matrix visualization"
#define READER_EXTENSION        "gp"

#define READERGP_GNUPLOT_BOXTEMPLATE(i, x1, y1, x2, y2) "set object %d rect from %.1f,%.1f to %.1f,%.1f fc rgb \"grey\"\n", (i), (x1), (y1), (x2), (y2)
#define READERGP_GNUPLOT_HEADER(outputname) "set terminal pdf\nset output \"%s.pdf\"\nunset xtics\nunset ytics\nunset border\nunset key\nset style fill solid 1.0 noborder\nset size ratio -1\n", (outputname)
#define READERGP_GNUPLOT_RANGES(xmax, ymax) "set xrange [0:%d]\nset yrange[%d:0]\n", (xmax), (ymax)
#define READERGP_GNUPLOT_PLOTCMD "plot \"-\" using 1:2:3 with circles fc rgb \"black\"\n"

/*
 * Local methods
 */

/** write file header with terminal etc. */
static
SCIP_RETCODE writeFileHeader(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file,               /**< File pointer to write to */
   const char*           outname             /**< the name of the gnuplot outputname */
   )
{

   SCIPinfoMessage(scip, file, READERGP_GNUPLOT_HEADER(outname));
   SCIPinfoMessage(scip, file, READERGP_GNUPLOT_RANGES(SCIPgetNVars(scip), SCIPgetNConss(scip)));
   return SCIP_OKAY;
}

/** write decomposition header such as rectangles for blocks etc. */
static
SCIP_RETCODE writeDecompositionHeader(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file,               /**< File pointer to write to */
   DEC_DECOMP*           decdecomp           /**< Decomposition pointer */
   )
{

   int i;
   int startx;
   int starty;
   int endx;
   int endy;
   assert(scip != NULL);
   assert(file != NULL);
   assert(decdecomp != NULL);
   if( decdecomp->type == DEC_DECTYPE_UNKNOWN || decdecomp->nblocks == 0 )
   {
      return SCIP_OKAY;
   }

   if( decdecomp->type == DEC_DECTYPE_ARROWHEAD || decdecomp->type == DEC_DECTYPE_BORDERED )
   {
      startx = 0;
      starty = 0;
      endx = 0;
      endy = 0;

      for( i = 0; i < decdecomp->nblocks; ++i )
      {
         endx += decdecomp->nsubscipvars[i];
         endy += decdecomp->nsubscipconss[i];
         SCIPinfoMessage(scip, file, READERGP_GNUPLOT_BOXTEMPLATE(i+1, startx+0.5, starty+0.5, endx+0.5, endy+0.5));
         startx = endx;
         starty = endy;
      }
      endx += decdecomp->nlinkingvars;
      endy += decdecomp->nlinkingconss;
      SCIPinfoMessage(scip, file, READERGP_GNUPLOT_BOXTEMPLATE(i+2, 0.5, starty+0.5, endx+0.5, endy+0.5));
      SCIPinfoMessage(scip, file, READERGP_GNUPLOT_BOXTEMPLATE(i+3, startx+0.5, +0.5, endx+0.5, endy+0.5));
      SCIPinfoMessage(scip, file, READERGP_GNUPLOT_BOXTEMPLATE(i+4, startx+0.5, starty+0.5, endx+0.5, endy+0.5));
   }
   return SCIP_OKAY;
}

/** write the plot commands */
static
SCIP_RETCODE writePlotCommands(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file                /**< File pointer to write to */
   )
{
   assert(scip != NULL);
   assert(file != NULL);

   SCIPinfoMessage(scip, file, READERGP_GNUPLOT_PLOTCMD);
   return SCIP_OKAY;
}

/** write the data optionally using the decomposition data */
static
SCIP_RETCODE writeData(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file,               /**< File pointer to write to */
   DEC_DECOMP*           decdecomp           /**< Decomposition pointer */
   )
{
   SCIP_VAR** vars;
   SCIP_CONS** conss;

   SCIP_HASHMAP* varindexmap;
   SCIP_HASHMAP* consindexmap;
   int nvars;
   int nconss;

   int i;
   int j;
   size_t varindex;
   size_t consindex;

   assert(scip != NULL);
   assert(file != NULL);

   conss = SCIPgetConss(scip);
   nconss = SCIPgetNConss(scip);

   varindexmap = NULL;
   consindexmap = NULL;

   SCIP_CALL( SCIPhashmapCreate(&varindexmap, SCIPblkmem(scip), SCIPgetNVars(scip)) );
   SCIP_CALL( SCIPhashmapCreate(&consindexmap, SCIPblkmem(scip), SCIPgetNConss(scip)) );

   if( decdecomp != NULL )
   {
      assert(decdecomp->type == DEC_DECTYPE_ARROWHEAD
               || decdecomp->type == DEC_DECTYPE_BORDERED
               || decdecomp->type == DEC_DECTYPE_DIAGONAL
               || decdecomp->type == DEC_DECTYPE_UNKNOWN
               || decdecomp->type == DEC_DECTYPE_STAIRCASE);

      /* if we don't have staicase, but something else, go through the blocks and create the indices */
      if( decdecomp->type == DEC_DECTYPE_ARROWHEAD || decdecomp->type == DEC_DECTYPE_BORDERED || decdecomp->type == DEC_DECTYPE_DIAGONAL )
      {
         SCIPdebugMessage("Block information:\n");
         varindex = 1;
         consindex = 1;
         for( i = 0; i < decdecomp->nblocks; ++i )
         {
            SCIPdebugPrintf("Block %d:\n", i+1);
            SCIPdebugPrintf("\tVars: %d", decdecomp->nsubscipvars[i]);
            SCIPdebugPrintf("\tConss: %d\n", decdecomp->nsubscipconss[i]);
            for( j = 0; j < decdecomp->nsubscipvars[i]; ++j )
            {
               assert(decdecomp->subscipvars[i][j] != NULL);
               SCIP_CALL( SCIPhashmapInsert(varindexmap, decdecomp->subscipvars[i][j], (void*)varindex) );
               varindex++;
            }
            for( j = 0; j < decdecomp->nsubscipconss[i]; ++j )
            {
               assert(decdecomp->subscipconss[i][j] != NULL);
               SCIP_CALL( SCIPhashmapInsert(consindexmap, decdecomp->subscipconss[i][j], (void*)consindex) );
               consindex++;
            }
         }

         SCIPdebugPrintf("Linking:\n");
         SCIPdebugPrintf("\tVars: %d", decdecomp->nlinkingvars);
         SCIPdebugPrintf("\tConss: %d\n\n", decdecomp->nlinkingconss);

         for( j = 0; j < decdecomp->nlinkingvars; ++j )
         {
            assert(decdecomp->linkingvars[j]);
            SCIP_CALL( SCIPhashmapInsert(varindexmap, decdecomp->linkingvars[j], (void*)varindex) );
            varindex++;
         }
         for( j = 0; j < decdecomp->nlinkingconss; ++j )
         {
            assert(decdecomp->linkingconss[j]);
            SCIP_CALL( SCIPhashmapInsert(consindexmap, decdecomp->linkingconss[j], (void*)consindex) );
            consindex++;
         }
      }
      else if( decdecomp->type == DEC_DECTYPE_STAIRCASE )
      {
         varindexmap = decdecomp->varindex;
         consindexmap = decdecomp->consindex;
      }
   }

   for( i = 0; i < nconss; i++ )
   {
      nvars = SCIPgetNVarsXXX(scip, conss[i]);
      vars = NULL;

      if( nvars > 0 )
      {
         SCIP_CALL( SCIPallocBufferArray( scip, &vars, nvars) );
         SCIP_CALL( SCIPgetVarsXXX(scip, conss[i], vars, nvars) );
      }

      for( j = 0; j < nvars; j++ )
      {
         assert(vars != NULL);

         /* if the problem has been created, output the whole model */
         if( SCIPgetStage(scip) == SCIP_STAGE_PROBLEM )
         {
            SCIPinfoMessage(scip, file, "%d %d 0.5\n", SCIPvarGetIndex(vars[j]), i);
            continue;
         }

         /* if there is no decomposition, output the presolved model! */
         if( decdecomp == NULL || decdecomp->type == DEC_DECTYPE_UNKNOWN )
         {
            SCIPinfoMessage(scip, file, "%d %d 0.5\n", SCIPvarGetIndex(vars[j]), i);
         }
         /* if there is a decomposition, output the indices derived from the decomposition above*/
         else
         {
            assert(varindexmap != NULL);
            assert(consindexmap != NULL);
            assert(SCIPhashmapGetImage(varindexmap, SCIPvarGetProbvar(vars[j])) != NULL);
            assert(SCIPhashmapGetImage(consindexmap, conss[i]) != NULL);

            SCIPinfoMessage(scip, file, "%d %d 0.5\n",
               SCIPhashmapGetImage(varindexmap, SCIPvarGetProbvar(vars[j])),
               SCIPhashmapGetImage(consindexmap, conss[i])
               );
         }
      }

      SCIPfreeBufferArrayNull(scip, &vars);
   }

   if( decdecomp != NULL && decdecomp->type != DEC_DECTYPE_STAIRCASE )
   {
      SCIPhashmapFree(&varindexmap);
      SCIPhashmapFree(&consindexmap);
   }

   return SCIP_OKAY;
}


/** write trailer of the file */
static
SCIP_RETCODE writeFileTrailer(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file                /**< File pointer to write to */
   )
{
   SCIPinfoMessage(scip, file, "e\n");
   return SCIP_OKAY;
}


/*
 * Callback methods of reader
 */

#define readerCopyGp NULL

/** destructor of reader to free user data (called when SCIP is exiting) */
static
SCIP_DECL_READERFREE(readerFreeGp)
{
   assert(strcmp(SCIPreaderGetName(reader), READER_NAME) == 0);
   return SCIP_OKAY;
}


/** problem reading method of reader */
#define readerReadGp NULL



/** problem writing method of reader */
static
SCIP_DECL_READERWRITE(readerWriteGp)
{
   /*lint --e{715}*/
   assert(scip != NULL);

   SCIP_CALL( SCIPwriteGp(scip, file, DECgetBestDecomp(scip), TRUE) );

   *result = SCIP_SUCCESS;
   return SCIP_OKAY;
}


/*
 * reader specific interface methods
 */

/** writes the decomposition to the specific file */
SCIP_RETCODE SCIPwriteGp(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file,               /**< File pointer to write to */
   DEC_DECOMP*           decdecomp,          /**< Decomposition pointer */
   SCIP_Bool             writeDecomposition  /**< whether to write decomposed problem */
   )
{
   char probname[SCIP_MAXSTRLEN];
   char outname[SCIP_MAXSTRLEN];
   char *name;

   assert(scip != NULL);
   assert(file != NULL);

   if( writeDecomposition && decdecomp == NULL )
   {
      SCIPwarningMessage(scip, "Cannot write decomposed problem if decomposition structure empty!");
      writeDecomposition = FALSE;
   }
   /* sanitize filename */
   (void) SCIPsnprintf(probname, SCIP_MAXSTRLEN, "%s", SCIPgetProbName(scip));
   SCIPsplitFilename(probname, NULL, &name, NULL, NULL);

   /* print header */
   if( decdecomp == NULL )
      (void) SCIPsnprintf(outname, SCIP_MAXSTRLEN, "%s", name);
   else
      (void) SCIPsnprintf(outname, SCIP_MAXSTRLEN, "%s_%d", name, decdecomp->nblocks);

   SCIP_CALL( writeFileHeader(scip, file, outname) );

   /* write decomp information such as rectangles */
   if( writeDecomposition )
      SCIP_CALL( writeDecompositionHeader(scip, file, decdecomp) );

   /* write the plot header*/
   SCIP_CALL( writePlotCommands(scip, file) );

   /* write data */
   SCIP_CALL( writeData(scip, file, decdecomp) );

   /* write file end */
   SCIP_CALL( writeFileTrailer(scip, file) );
   return SCIP_OKAY;
}


/** includes the gp file reader into SCIP */
SCIP_RETCODE SCIPincludeReaderGp(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   /* include gp reader */
   SCIP_CALL( SCIPincludeReader(scip, READER_NAME, READER_DESC, READER_EXTENSION,
         readerCopyGp, readerFreeGp, readerReadGp, readerWriteGp, NULL) );

   return SCIP_OKAY;
}
