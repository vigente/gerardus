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

/**@file    scipParaObjMessageHdlr.h
 * @brief   SCIP message handler for ParaSCIP and FiberSCIP.
 * @author  Yuji Shinano
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_PARA_OBJ_MESSAGE_HDLR_H__
#define __SCIP_PARA_OBJ_MESSAGE_HDLR_H__

#include "objscip/objmessagehdlr.h"
#include "scip/scip.h"
#include "ug/paraComm.h"

namespace ParaSCIP
{

/** C++ wrapper object for file readers */
class ScipParaObjMessageHdlr: public scip::ObjMessagehdlr
{
   // MYMESSAGEHDLRDATA *mymessagehdlrdata;
   UG::ParaComm          *comm;
   FILE*                 logfile;            /**< log file where to copy messages into */
   SCIP_Bool             quiet;              /**< should screen messages be suppressed? */

   void logMessage(
         FILE*                 file,               /**< file stream to print message into */
         const char*           msg                 /**< message to print */
         );

public:

   /** default constructor */
   ScipParaObjMessageHdlr(
      UG::ParaComm       *comm,
      FILE*              file,
      SCIP_Bool          quiet,
      SCIP_Bool          bufferedoutput      /**< should the output be buffered up to the next newline? */
      );

   /** destructor */
   virtual ~ScipParaObjMessageHdlr();

   /** error message print method of message handler
    *
    *  This method is invoked, if SCIP wants to display an error message to the screen or a file
    */
   virtual void scip_error(
      SCIP_MESSAGEHDLR*  messagehdlr,        /**< the message handler itself */
      FILE*              file,               /**< file stream to print into */
      const char*        msg                 /**< string to output into the file */
      );

   /** warning message print method of message handler
    *
    *  This method is invoked, if SCIP wants to display a warning message to the screen or a file
    */
   virtual void scip_warning(
      SCIP_MESSAGEHDLR*  messagehdlr,        /**< the message handler itself */
      FILE*              file,               /**< file stream to print into */
      const char*        msg                 /**< string to output into the file */
      );

   /** dialog message print method of message handler
    *
    *  This method is invoked, if SCIP wants to display a dialog message to the screen or a file
    */
   virtual void scip_dialog(
      SCIP_MESSAGEHDLR*  messagehdlr,        /**< the message handler itself */
      FILE*              file,               /**< file stream to print into */
      const char*        msg                 /**< string to output into the file */
      );

   /** info message print method of message handler
    *
    *  This method is invoked, if SCIP wants to display an information message to the screen or a file
    */
   virtual void scip_info(
      SCIP_MESSAGEHDLR*  messagehdlr,        /**< the message handler itself */
      FILE*              file,               /**< file stream to print into */
      const char*        msg                 /**< string to output into the file */
      );
};

}

#endif  // __SCIP_PARA_OBJ_MESSAGE_HDLR_H__
