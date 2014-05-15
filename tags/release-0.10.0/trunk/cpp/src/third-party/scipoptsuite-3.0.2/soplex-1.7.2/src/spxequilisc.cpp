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

/**@file  spxequilisc.cpp
 * @brief Equilibrium row/column scaling.
 */
#include <assert.h>

#include "spxequilisc.h"
#include "spxout.h"

namespace soplex
{
static const char* makename(bool doBoth)
{
   return doBoth ? "bi-Equilibrium" : "uni-Equilibrium";
}

SPxEquiliSC::SPxEquiliSC(bool doBoth)
   : SPxScaler(makename(doBoth), false, doBoth)
{}

SPxEquiliSC::SPxEquiliSC(const SPxEquiliSC& old)
   : SPxScaler(old)
{}

SPxEquiliSC& SPxEquiliSC::operator=(const SPxEquiliSC& rhs)
{
   if(this != &rhs)
   {
      SPxScaler::operator=(rhs);
   }

   return *this;
}

Real SPxEquiliSC::computeScale(Real /*mini*/, Real maxi) const
{
   METHOD( "SPxEquiliSC::computeScale()" );

   return maxi;
}

void SPxEquiliSC::scale(SPxLP& lp) 
{
   METHOD( "SPxEquiliSC::scale()" );

   MSG_INFO1( spxout << "IEQUSC01 Equilibrium scaling LP" << std::endl; )

   setup(lp);

   /* We want to do that direction first, with the lower ratio.
    * Reason:           
    *                               Rowratio
    *            0.04  0.02  0.01      4
    *            4000    20  1000    200
    * Colratio    1e5   1e3   1e5
    *
    * Row first =>                  Col next =>
    *               1   0.5  0.25         1   1   1 
    *               1   0.05 0.25         1  0.1  1
    *
    * Col first =>                  Row next =>
    *            1e-5  1e-3  1e-5        0.01  1  0.01
    *               1     1     1          1   1    1
    *
    */
   Real colratio = maxColRatio(lp);
   Real rowratio = maxRowRatio(lp);

   bool colFirst = colratio < rowratio;

   MSG_INFO2( spxout << "IEQUSC02 LP scaling statistics:" 
                        << " min= " << lp.minAbsNzo()
                        << " max= " << lp.maxAbsNzo()
                        << " col-ratio= " << colratio 
                        << " row-ratio= " << rowratio
                        << std::endl; )
   if (colFirst)
   {
      computeScalingVecs(lp.colSet(), m_rowscale, m_colscale);

      if (m_doBoth)
         computeScalingVecs(lp.rowSet(), m_colscale, m_rowscale);
   }
   else
   {
      computeScalingVecs(lp.rowSet(), m_colscale, m_rowscale);

      if (m_doBoth)
         computeScalingVecs(lp.colSet(), m_rowscale, m_colscale);
   }
   applyScaling(lp);

   MSG_INFO3( spxout << "IEQUSC03 \tRow scaling min= " << minAbsRowscale()
                        << " max= " << maxAbsRowscale()
                        << std::endl
                        << "\tCol scaling min= " << minAbsColscale()
                        << " max= " << maxAbsColscale()
                        << std::endl; )

   MSG_INFO2( spxout << "IEQUSC04 LP scaling statistics:" 
                        << " min= " << lp.minAbsNzo()
                        << " max= " << lp.maxAbsNzo()
                        << " col-ratio= " << maxColRatio(lp) 
                        << " row-ratio= " << maxRowRatio(lp) 
                        << std::endl; )
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





