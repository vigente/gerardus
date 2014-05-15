/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the class library                   */
/*       SoPlex --- the Sequential object-oriented simPlex.                  */
/*                                                                           */
/*    Copyright (C) 1996-2013 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SoPlex is distributed under the terms of the ZIB Academic Licence.       */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SoPlex; see the file COPYING. If not email to soplex@zib.de.  */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */


/**@file  timer.h
 * @brief TraceMethod class.
 */

#ifndef _TRACEMETHOD_H_
#define _TRACEMETHOD_H_


#include <iostream>


#define FILE_NAME_COL  60

namespace soplex
{

/**@class TraceMethod 
   @brief Helper class to trace the calling of methods.

   Prints out a message with the current filename and line to indicate
   that some method has been called.
*/
class TraceMethod
{
private:

   //------------------------------------
   /**@name Data */
   //@{
   /// the current indentation
   static int s_indent;
   //@}

public:

   //------------------------------------
   /**@name Construction / destruction */
   //@{
   /// constructor
   TraceMethod(const char* s, const char* file, int line );
 
  /// destructor
   virtual ~TraceMethod()
   {
      s_indent--;
   }
   //@}
};


}


#endif // _TRACEMETHOD_H_
