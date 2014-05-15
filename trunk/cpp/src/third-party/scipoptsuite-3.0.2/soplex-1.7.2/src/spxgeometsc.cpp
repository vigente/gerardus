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

/**@file  spxgeometsc.cpp
 * @brief Geometric mean row/column scaling.
 */
#include <assert.h>

#include "spxgeometsc.h"
#include "spxout.h"

namespace soplex
{
/**@param maxIters   arbitrary small number, we choose 8
   @param minImpr    Bixby said Fourer said in MP 23, 274 ff. that 0.9 is a good value.
   @param goodEnough if the max/min ratio is allready less then 1000/1 we do not scale.
*/ 
SPxGeometSC::SPxGeometSC(int maxIters, Real minImpr, Real goodEnough)
   : SPxScaler("Geometric")
   , m_maxIterations(maxIters)
   , m_minImprovement(minImpr)
   , m_goodEnoughRatio(goodEnough)
{}

SPxGeometSC::SPxGeometSC(const SPxGeometSC& old)
   : SPxScaler(old)
   , m_maxIterations(old.m_maxIterations)
   , m_minImprovement(old.m_minImprovement)
   , m_goodEnoughRatio(old.m_goodEnoughRatio)
{}

SPxGeometSC& SPxGeometSC::operator=(const SPxGeometSC& rhs)
{
   if (this != &rhs)
   {
      SPxScaler::operator=(rhs);
   }

   return *this;
}

Real SPxGeometSC::computeScale(Real mini, Real maxi) const
{
   METHOD( "SPxGeometSC::computeScale()" );

   return sqrt(mini * maxi);
}

void SPxGeometSC::scale(SPxLP& lp) 
{
   METHOD( "SPxGeometSC::scale()" );

   MSG_INFO1( spxout << "IGEOSC01 Geometric scaling LP" << std::endl; )

   Real pstart = 0.0;
   Real p0     = 0.0;
   Real p1     = 0.0;

   setup(lp);

   /* We want to do that direction first, with the lower ratio.
    * See SPxEquiliSC::scale() for a reasoning.
    */
   Real colratio = maxColRatio(lp);
   Real rowratio = maxRowRatio(lp);

   bool colFirst = colratio < rowratio;

   MSG_INFO2( spxout << "IGEOSC02 LP scaling statistics:" 
                        << " min= " << lp.minAbsNzo()
                        << " max= " << lp.maxAbsNzo()
                        << " col-ratio= " << colratio
                        << " row-ratio= " << rowratio
                        << std::endl; )

   // We make at most m_maxIterations. 
   for(int count = 0; count < m_maxIterations; count++)
   {
      if (colFirst)
      {
         p0 = computeScalingVecs(lp.colSet(), m_rowscale, m_colscale);
         p1 = computeScalingVecs(lp.rowSet(), m_colscale, m_rowscale);
      }
      else
      {
         p0 = computeScalingVecs(lp.rowSet(), m_colscale, m_rowscale);
         p1 = computeScalingVecs(lp.colSet(), m_rowscale, m_colscale);
      }
      MSG_INFO3( spxout << "IGEOSC03 Geometric scaling round " << count
                           << " col-ratio= " << (colFirst ? p0 : p1)
                           << " row-ratio= " << (colFirst ? p1 : p0)
                           << std::endl; )

      // record start value, this is done with m_col/rowscale = 1.0, so it is the
      // value frome the "original" (as passed to the scaler) LP.
      if (count == 0)
      {
         pstart = p0;
         // are we allready good enough ?
         if (pstart < m_goodEnoughRatio)
            break;
      }
      else // do not test at the first iteration, then abort if no improvement.
         if (p1 > m_minImprovement * p0)
            break;
   }      
   
   // we scale only if either:
   // - we had at the beginng a ratio worse then 1000/1
   // - we have at least a 15% improvement.
   if (pstart < m_goodEnoughRatio || p1 > pstart * m_minImprovement)
   {
      // reset m_colscale/m_rowscale to 1.0
      setup(lp);

      MSG_INFO2( spxout << "IGEOSC04 No scaling done." << std::endl; )
   }
   else
   {
      applyScaling(lp);

      MSG_INFO3( spxout << "IGEOSC05 Row scaling min= " << minAbsRowscale()
                           << " max= " << maxAbsRowscale()
                           << std::endl
                           << "IGEOSC06 Col scaling min= " << minAbsColscale()
                           << " max= " << maxAbsColscale()
                           << std::endl; )

      MSG_INFO2( spxout << "IGEOSC07 LP scaling statistics:" 
                           << " min= " << lp.minAbsNzo()
                           << " max= " << lp.maxAbsNzo()
                           << " col-ratio= " << maxColRatio(lp) 
                           << " row-ratio= " << maxRowRatio(lp) 
                           << std::endl; )
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



