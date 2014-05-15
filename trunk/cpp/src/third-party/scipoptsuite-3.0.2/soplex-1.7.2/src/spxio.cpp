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


//#define DEBUGGING 1

#include <iostream>
#include <assert.h>
#include <string.h>

#include "spxdefines.h"
#include "spxlp.h"

#include "dvector.h"
#include "dataarray.h"
#include "lprow.h"
#include "lpcol.h"
#include "lprowset.h"
#include "lpcolset.h"
#include "nameset.h"
#include "spxout.h"
#include "spxfileio.h"

namespace soplex
{
/**@param is       input stream. 
 * @param rowNames contains after the call the names of the constraints 
 *                 (rows) in the same order as the rows in the LP.
 *                 Constraints without a name (only possible with LPF 
 *                 files) are automatically assigned a name.
 *                 Maybe 0 if the names are not needed.
 * @param colNames contains after the call the names of the variables
 *                 (columns) in the same order as the columns in the LP.
 *                 Maybe 0 if the names are not needed.
 * @param intVars  contains after the call the indices of those variables
 *                 that where marked as beeing integer in the file.
 *                 Maybe 0 if the information is not needed.
 * @todo  Make sure the Id's in the NameSet%s are the same as in the LP.
 */
bool SPxLP::read(std::istream&   is, 
                 NameSet*        rowNames,
                 NameSet*        colNames,
                 DIdxSet*        intVars)
{
   bool ok;
   char c;

   is.get(c);
   is.putback(c);

   /* MPS starts either with a comment mark '*' or with the keyword
    * 'NAME' at the first column.
    * LPF starts either with blanks, a comment mark '\' or with
    * the keyword "MAX" or "MIN" in upper or lower case.
    * There is no possible valid LPF file starting with a '*' or 'N'.
    */
   ok = ((c == '*') || (c == 'N'))
      ? readMPS(is, rowNames, colNames, intVars)
      : readLPF(is, rowNames, colNames, intVars);

   MSG_DEBUG( spxout << "DSPXIO01\n" << *this; );

   return ok;
}


bool SPxLP::readFile( 
   const char* filename, 
   NameSet*    rowNames,
   NameSet*    colNames, 
   DIdxSet*    intVars)
{
   METHOD( "SPxLP::readFile()" );

   spxifstream file(filename);

   if (!file)
      return false;

   return read(file, rowNames, colNames, intVars);
}


/** write loaded LP to \p filename
 */
void SPxLP::writeFile(const char* filename,
      const NameSet* rowNames,
      const NameSet* colNames, 
      const DIdxSet* p_intvars ) const
{
   std::ofstream tmp(filename);
   size_t len_f = strlen(filename);
   if (len_f > 4 && filename[len_f-1] == 's' && filename[len_f-2] == 'p' && filename[len_f-3] == 'm' && filename[len_f-4] == '.')
   {
      writeMPS(tmp, rowNames, colNames, p_intvars);
   }
   else
   {
      writeLPF(tmp, rowNames, colNames, p_intvars);
   }
}

} // namespace soplex

//-----------------------------------------------------------------------------
//Emacs Local Variables:
//Emacs mode:c++
//Emacs c-basic-offset:3
//Emacs tab-width:8
//Emacs indent-tabs-mode:nil
//Emacs End:
//-----------------------------------------------------------------------------
