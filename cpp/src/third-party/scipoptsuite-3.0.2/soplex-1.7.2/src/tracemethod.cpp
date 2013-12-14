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

#include "tracemethod.h"
#include "spxout.h"


namespace soplex
{
#if defined(TRACE_METHOD)

   int TraceMethod::s_indent = 0;

   /// constructor
   TraceMethod::TraceMethod(const char* s, const char* file, int line )
   {
      int i;
 
      spxout << "\t";
      
      for(i = 0; i < s_indent; i++)
         spxout << ".";      
      
      spxout << s;
      
      for(i = strlen(s) + s_indent; i < FILE_NAME_COL - 8; i++)
         spxout << "_";             
      spxout << "[" << file << ":" << line << "]" << std::endl; 
      s_indent++;
   }
#endif //TRACE_METHOD
}
