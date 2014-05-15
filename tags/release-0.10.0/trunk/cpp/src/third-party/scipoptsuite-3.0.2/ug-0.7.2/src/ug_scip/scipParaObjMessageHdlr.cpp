/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*             This file is part of the program and software framework       */
/*                  UG --- Ubquity Generator Framework                       */
/*                                                                           */
/*    Copyright (C) 2010-2013 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  UG is distributed under the terms of the ZIB Academic Licence.           */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with UG; see the file COPYING. If not email to scip@zib.de.        */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file    scipParaObjMessageHdlr.cpp
 * @brief   SCIP message handler for ParaSCIP and FiberSCIP.
 * @author  Yuji Shinano
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <cassert>
#include <iostream>

#include "scipParaObjMessageHdlr.h"

using namespace ParaSCIP;


/** constructor */
ScipParaObjMessageHdlr::ScipParaObjMessageHdlr(
   UG::ParaComm           *inComm,
   FILE*              inFile,
   SCIP_Bool          inQuiet,
   SCIP_Bool          inBufferedoutput      /**< should the output be buffered up to the next newline? */
   )  : ObjMessagehdlr(inBufferedoutput)
{
   comm = inComm;     // for debugging
   logfile = inFile;
   quiet = inQuiet;
}

ScipParaObjMessageHdlr::~ScipParaObjMessageHdlr(
      )
{
}


void ScipParaObjMessageHdlr::logMessage(
   FILE*                 file,               /**< file stream to print message into */
   const char*           msg                 /**< message to print */
   )
{

   if( !quiet || (file != stdout && file != stderr) )
   {
      fputs(msg, file);
      fflush(file);
   }
   if( logfile != NULL && (file == stdout || file == stderr) )
   {
      fputs(msg, logfile);
      fflush(logfile);
   }
}

/** error message print method of message handler
*
*  This method is invoked, if SCIP wants to display an error message to the screen or a file
*/
void ScipParaObjMessageHdlr::scip_error(
   SCIP_MESSAGEHDLR*  messagehdlr,        /**< the message handler itself */
   FILE*              file,               /**< file stream to print into */
   const char*        msg                 /**< string to output into the file */
   )
{
   logMessage(file, msg);
}

/** warning message print method of message handler
 *
 *  This method is invoked, if SCIP wants to display a warning message to the screen or a file
 */
void ScipParaObjMessageHdlr::scip_warning(
   SCIP_MESSAGEHDLR*  messagehdlr,        /**< the message handler itself */
   FILE*              file,               /**< file stream to print into */
   const char*        msg                 /**< string to output into the file */
   )
{
   // logMessage(mymessagehdlrdata, file, msg);
   logMessage(file, msg);
}

/** dialog message print method of message handler
 *
 *  This method is invoked, if SCIP wants to display a dialog message to the screen or a file
 */
void ScipParaObjMessageHdlr::scip_dialog(
   SCIP_MESSAGEHDLR*  messagehdlr,        /**< the message handler itself */
   FILE*              file,               /**< file stream to print into */
   const char*        msg                 /**< string to output into the file */
   )
{
   logMessage(file, msg);
}

/** info message print method of message handler
 *
 *  This method is invoked, if SCIP wants to display an information message to the screen or a file
 */
void ScipParaObjMessageHdlr::scip_info(
   SCIP_MESSAGEHDLR*  messagehdlr,        /**< the message handler itself */
   FILE*              file,               /**< file stream to print into */
   const char*        msg                 /**< string to output into the file */
   )
{
   logMessage(file, msg);
}


