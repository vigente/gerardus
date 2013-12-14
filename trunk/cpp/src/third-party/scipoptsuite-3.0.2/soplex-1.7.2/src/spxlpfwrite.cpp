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

/**@file  spxlpfwrite.cpp
 * @brief Write LP as LPF format file.
 */
#define NUM_ENTRIES_PER_LINE 5 

#include <assert.h>
#include <stdio.h>
#include <fstream>

#include "spxdefines.h"
#include "spxlp.h"

namespace soplex
{
   namespace {

      // get the name of a row or construct one
      const char* getRowName(
         const SPxLP&    p_lp,
         int             p_idx,
         const NameSet*  p_rnames, 
         char*           p_buf,
         int             p_num_written_rows)
      {
         assert(p_buf != 0);
         assert(p_idx >= 0);
         assert(p_idx <  p_lp.nRows());

         if (p_rnames != 0) 
         {
            DataKey key = p_lp.rId(p_idx);

            if (p_rnames->has(key))
               return (*p_rnames)[key];
         }
         sprintf(p_buf, "C%d", p_num_written_rows);
         
         return p_buf;
      }
         
      // get the name of a column or construct one
      const char* getColName(
         const SPxLP&    p_lp,
         int             p_idx,
         const NameSet*  p_cnames, 
         char*           p_buf)
      {
         assert(p_buf != 0);
         assert(p_idx >= 0);
         assert(p_idx <  p_lp.nCols());

         if (p_cnames != 0) 
         {
            DataKey key = p_lp.cId(p_idx);

            if (p_cnames->has(key))
               return (*p_cnames)[key];
         }
         sprintf(p_buf, "x%d", p_idx);
         
         return p_buf;
      }
         
      // write an SVector
      void writeSVector( const SPxLP&     p_lp,       ///< the LP
                         std::ostream&    p_output,   ///< output stream
                         const NameSet*   p_cnames,   ///< column names
                         const SVector&   p_svec )    ///< vector to write
      {
         METHOD("writeSVector");

         char name[16];
         int  num_coeffs = 0;

         for (int j = 0; j < p_lp.nCols(); ++j)
         {
            const Real coeff = p_svec[j];
            if (coeff == 0)
               continue;

            if (num_coeffs == 0)
               p_output << coeff << " "
                        << getColName(p_lp, j, p_cnames, name);
            else
            {
               // insert a line break every NUM_ENTRIES_PER_LINE columns
               if (num_coeffs % NUM_ENTRIES_PER_LINE == 0)
                  p_output << "\n\t";

               if (coeff < 0)
                  p_output << " - " << -coeff;
               else
                  p_output << " + " << coeff;

               p_output << " " << getColName(p_lp, j, p_cnames, name);
            }
            ++num_coeffs;
         }
      }

      // write the objective
      void writeObjective( const SPxLP&   p_lp,       ///< the LP
                           std::ostream&  p_output,   ///< output stream
                           const NameSet* p_cnames )  ///< column names
      {
         METHOD("writeObjective");

         const int sense = p_lp.spxSense();

         p_output << ((sense == SPxLP::MINIMIZE) ? "Minimize\n" : "Maximize\n");
         p_output << "  obj: ";

         const Vector&  obj  = p_lp.maxObj();
         DSVector svec(obj.dim());
         svec.operator=(obj);
         svec *= Real(sense);
         writeSVector(p_lp, p_output, p_cnames, svec);
         p_output << "\n";
      } 

      // write non-ranged rows
      void writeRow( const SPxLP&     p_lp,       ///< the LP
                     std::ostream&    p_output,   ///< output stream
                     const NameSet*   p_cnames,   ///< column names
                     const SVector&   p_svec,     ///< vector of the row
                     const Real&      p_lhs,      ///< lhs of the row
                     const Real&      p_rhs )     ///< rhs of the row
      {
         METHOD("writeRow");

         writeSVector(p_lp, p_output, p_cnames, p_svec);

         if ( p_lhs == p_rhs )
            p_output << " = " << p_rhs;
         else if ( p_lhs <= -infinity )
            p_output << " <= " << p_rhs;
         else {
            assert(p_rhs >= infinity);
            p_output << " >= " << p_lhs;
         }

         p_output << "\n";
      }
        
      // write all rows
      void writeRows( const SPxLP&   p_lp,       ///< the LP
                      std::ostream&  p_output,   ///< output stream
                      const NameSet* p_rnames,   ///< row names
                      const NameSet* p_cnames )  ///< column names
      {
         METHOD("writeRows");
         char      name[16];

         p_output << "Subject To\n";

         int num_written_rows = 0;  // num_written_rows > nRows with ranged rows
         for (int i = 0; i < p_lp.nRows(); ++i)
         {
            const Real lhs = p_lp.lhs(i);
            const Real rhs = p_lp.rhs(i);

            if (lhs > -infinity && rhs < infinity && lhs != rhs) {
               // ranged row -> write two non-ranged rows
               p_output << " " << getRowName(p_lp, i, p_rnames, 
                                             name, ++num_written_rows) << "_1 : ";
               writeRow (p_lp, p_output, p_cnames, p_lp.rowVector(i), lhs, infinity);

               p_output << " " << getRowName(p_lp, i, p_rnames, 
                                             name, ++num_written_rows) << "_2 : ";
               writeRow (p_lp, p_output, p_cnames, p_lp.rowVector(i), -infinity, rhs);
            }
            else {
               p_output << " " << getRowName(p_lp, i, p_rnames, 
                                             name, ++num_written_rows) << " : ";
               writeRow (p_lp, p_output, p_cnames, p_lp.rowVector(i), lhs, rhs);
            }
         }
      }

      // write the variable bounds
      // (the default bounds 0 <= x <= infinity are not written)
      void writeBounds( const SPxLP&   p_lp,       ///< the LP to write
                        std::ostream&  p_output,   ///< output stream
                        const NameSet* p_cnames )  ///< column names
      {
         METHOD("writeBounds");
         char      name[16];

         p_output << "Bounds\n";
         for (int j = 0; j < p_lp.nCols(); ++j)
         {
            const Real lower = p_lp.lower(j);
            const Real upper = p_lp.upper(j);

            if (lower == upper) {
               p_output << "  "   << getColName(p_lp, j, p_cnames, name) 
                        << " = "  << upper << '\n';
            }
            else if (lower > -infinity) {
               if (upper < infinity) {
                  // range bound
                  if (lower != 0)
                     p_output << "  "   << lower << " <= "
                              << getColName(p_lp, j, p_cnames, name) 
                              << " <= " << upper << '\n';
                  else
                     p_output << "  "   << getColName(p_lp, j, p_cnames, name) 
                              << " <= " << upper << '\n';
               }
               else if (lower != 0)
                  p_output << "  " << lower << " <= "
                           << getColName(p_lp, j, p_cnames, name)
                           << '\n';
            }
            else if (upper < infinity)
               p_output << "   -Inf <= "
                        << getColName(p_lp, j, p_cnames, name) 
                        << " <= " << upper << '\n';
            else
               p_output << "  "   << getColName(p_lp, j, p_cnames, name) 
                        << " free\n";
         }
      }
         
      // write the generals section 
      void writeGenerals( const SPxLP&   p_lp,         ///< the LP to write
                          std::ostream&  p_output,     ///< output stream
                          const NameSet* p_cnames,     ///< column names
                          const DIdxSet* p_intvars )   ///< integer variables
      {
         METHOD("writeGenerals");
         char  name[16];

         if (p_intvars == NULL || p_intvars->size() <= 0) 
            return;  // no integer variables
         
         p_output << "Generals\n";
         for (int j = 0; j < p_lp.nCols(); ++j)
            if (p_intvars->number(j) >= 0)
               p_output << "  " << getColName(p_lp, j, p_cnames, name) << "\n";
      }
   } // anonymous namespace
   
/// Write LP in "LPF File Format".
void SPxLP::writeLPF( std::ostream&  p_output,          ///< output stream
                      const NameSet* p_rnames,          ///< row names
                      const NameSet* p_cnames,          ///< column names
                      const DIdxSet* p_intvars )        ///< integer variables
   const
{
   METHOD("writeLPF");

   p_output << std::setprecision(15);
   writeObjective (*this, p_output, p_cnames);
   writeRows      (*this, p_output, p_rnames, p_cnames);
   writeBounds    (*this, p_output, p_cnames);
   writeGenerals  (*this, p_output, p_cnames, p_intvars);
   p_output << "End" << std::endl;
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

