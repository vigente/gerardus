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

/**@file  spxscaler.cpp
 * @brief LP scaling base class.
 */
#include <iostream>
#include <assert.h>

#include "spxscaler.h"
#include "spxlp.h"

namespace soplex
{

std::ostream& operator<<(std::ostream& s, const SPxScaler& sc)
{
   s << sc.getName() << " scaler:" << std::endl;
   s << "colscale = [ ";
   for(int ci = 0; ci < sc.m_colscale.size(); ++ci )
      s << sc.m_colscale[ci] << " ";
   s << "]" << std::endl;
      
   s << "rowscale = [ ";
   for(int ri = 0; ri < sc.m_rowscale.size(); ++ri )
      s << sc.m_rowscale[ri] << " ";
   s << "]" << std::endl;

   return s;
}

SPxScaler::SPxScaler(
   const char* name, 
   bool        colFirst, 
   bool        doBoth) 
   : m_name(name)
   , m_colFirst(colFirst)
   , m_doBoth(doBoth)
{
   assert(SPxScaler::isConsistent());
}

SPxScaler::SPxScaler(const SPxScaler& old)
   : m_name(old.m_name)
   , m_colscale(old.m_colscale)
   , m_rowscale(old.m_rowscale)
   , m_colFirst(old.m_colFirst)
   , m_doBoth(old.m_doBoth)
{
   assert(SPxScaler::isConsistent());
}

SPxScaler::~SPxScaler()
{
   m_name = 0;
}   

SPxScaler& SPxScaler::operator=(const SPxScaler& rhs)
{
   if (this != &rhs)
   {
      m_name     = rhs.m_name;
      m_colscale = rhs.m_colscale;
      m_rowscale = rhs.m_rowscale;
      m_colFirst = rhs.m_colFirst;
      m_doBoth   = rhs.m_doBoth;

      assert(SPxScaler::isConsistent());
   }
   return *this;
}


const char* SPxScaler::getName() const
{
   METHOD( "SPxScaler::getName()" );

   return m_name;
}

void SPxScaler::setOrder(bool colFirst)
{
   METHOD( "SPxScaler::setOrder()" );

   m_colFirst = colFirst;
}

void SPxScaler::setBoth(bool both)
{
   METHOD( "SPxScaler::setBoth()" );

   m_doBoth = both;
}

void SPxScaler::setup(SPxLP& lp)
{
   METHOD( "SPxScaler::setup()" );

   assert(lp.isConsistent());

   m_colscale.reSize(lp.nCols());
   m_rowscale.reSize(lp.nRows());

   int i;

   for(i = 0; i < lp.nCols(); ++i )
      m_colscale[i] = 1.0;

   for(i = 0; i < lp.nRows(); ++i )
      m_rowscale[i] = 1.0;
}

/** This function is used by computeScaleVecs and has to be overridden.
 */
Real SPxScaler::computeScale(Real /*mini*/, Real /*maxi*/) const
{
   METHOD( "SPxScaler::computeScale" );

   return 1.0;
}

Real SPxScaler::computeScalingVecs(
   const SVSet*           vecset, 
   const DataArray<Real>& coScaleval, 
   DataArray<Real>&       scaleval) 
{
   METHOD( "SPxScaler::computeScalingVecs()" );

   Real pmax = 0.0;

   for(int i = 0; i < vecset->num(); ++i )
   {
      const SVector& vec = (*vecset)[i];
            
      Real maxi = 0.0;
      Real mini = infinity;
            
      for( int j = 0; j < vec.size(); ++j)
      {
         Real x = fabs(vec.value(j) * coScaleval[vec.index(j)]);
               
         if (!isZero(x))
         {
            if (x > maxi)
               maxi = x;
            if (x < mini)
               mini = x;
         }
      }
      // empty rows/cols are possible
      if (mini == infinity || maxi == 0.0)
      {
         mini = 1.0;
         maxi = 1.0;
      }
      assert(mini < infinity);
      assert(maxi > 0.0);

      scaleval[i] = 1.0 / computeScale(mini, maxi);
            
      Real p = maxi / mini;
            
      if (p > pmax)
         pmax = p;
   }
   return pmax;
}

void SPxScaler::applyScaling(SPxLP& lp)
{
   METHOD( "SPxScaler::applyScaling()" );

   int i;

   for(i = 0; i < lp.nRows(); ++i )
   {
      SVector& vec = lp.rowVector_w(i);

      for( int j = 0; j < vec.size(); ++j)
         vec.value(j) *= m_colscale[vec.index(j)] * m_rowscale[i];

      if (lp.rhs(i) < infinity)
         lp.rhs_w(i) *= m_rowscale[i];
      if (lp.lhs(i) > -infinity)
         lp.lhs_w(i) *= m_rowscale[i];
   }
   for(i = 0; i < lp.nCols(); ++i )
   {
      SVector& vec = lp.colVector_w(i);
      
      for( int j = 0; j < vec.size(); ++j)
         vec.value(j) *= m_rowscale[vec.index(j)] * m_colscale[i];
      
      lp.maxObj_w(i) *= m_colscale[i];
      
      if (lp.upper(i) < infinity)
         lp.upper_w(i) /= m_colscale[i];
      if (lp.lower(i) > -infinity)
         lp.lower_w(i) /= m_colscale[i];
   }
   assert(lp.isConsistent());
}

void SPxScaler::unscalePrimal(Vector& x) const
{
   METHOD( "SPxScaler::unscalePrimal()" );

   assert(x.dim() == m_colscale.size());

   for(int j = 0; j < x.dim(); ++j)
      x[j] *= m_colscale[j];
}

void SPxScaler::unscaleSlacks(Vector& s) const
{
   METHOD( "SPxScaler::unscaleSlacks()" );

   assert(s.dim() == m_rowscale.size());

   for(int i = 0; i < s.dim(); ++i)
      s[i] /= m_rowscale[i];
}

void SPxScaler::unscaleDual(Vector& pi) const
{
   METHOD( "SPxScaler::unscaleDual()" );

   assert(pi.dim() == m_rowscale.size());

   for(int i = 0; i < pi.dim(); ++i)
      pi[i] *= m_rowscale[i];
}

void SPxScaler::unscaleRedCost(Vector& r) const
{
   METHOD( "SPxScaler::unscaleRedCost()" );

   assert(r.dim() == m_colscale.size());

   for(int j = 0; j < r.dim(); ++j)
      r[j] /= m_colscale[j];
}

Real SPxScaler::minAbsColscale() const
{
   METHOD( "SPxScaler::minAbsColscale()" );

   Real mini = infinity;

   for(int i = 0; i < m_colscale.size(); ++i)
      if (fabs(m_colscale[i]) < mini)
         mini = fabs(m_colscale[i]);

   return mini;
}

Real SPxScaler::maxAbsColscale() const
{
   METHOD( "SPxScaler::maxAbsColscale()" );

   Real maxi = 0.0;

   for(int i = 0; i < m_colscale.size(); ++i)
      if (fabs(m_colscale[i]) > maxi)
         maxi = fabs(m_colscale[i]);

   return maxi;
}

Real SPxScaler::minAbsRowscale() const
{
   METHOD( "SPxScaler::minAbsRowscale()" );

   Real mini = infinity;

   for(int i = 0; i < m_rowscale.size(); ++i)
      if (fabs(m_rowscale[i]) < mini)
         mini = fabs(m_rowscale[i]);

   return mini;
}

Real SPxScaler::maxAbsRowscale() const
{
   METHOD( "SPxScaler::maxAbsRowscale()" );

   Real maxi = 0.0;

   for(int i = 0; i < m_rowscale.size(); ++i)
      if (fabs(m_rowscale[i]) > maxi)
         maxi = fabs(m_rowscale[i]);

   return maxi;
}

/** \f$\max_{j\in\mbox{ cols}}
 *   \left(\frac{\max_{i\in\mbox{ rows}}|a_ij|}
 *              {\min_{i\in\mbox{ rows}}|a_ij|}\right)\f$
 */
Real SPxScaler::maxColRatio(const SPxLP& lp) const
{
   METHOD( "SPxScaler::maxColRatio()" );

   Real pmax = 0.0;

   for(int i = 0; i < lp.nCols(); ++i )
   {
      const SVector& vec  = lp.colVector(i);
      Real           mini = infinity;
      Real           maxi = 0.0;

      for(int j = 0; j < vec.size(); ++j)
      {
         Real x = fabs(vec.value(j));

         if (x < mini)
            mini = x;
         if (x > maxi)
            maxi = x;
      }
      Real p = maxi / mini;

      if (p > pmax)
         pmax = p;
   }   
   return pmax;
}

/** \f$\max_{i\in\mbox{ rows}}
 *   \left(\frac{\max_{j\in\mbox{ cols}}|a_ij|}
 *              {\min_{j\in\mbox{ cols}}|a_ij|}\right)\f$
 */
Real SPxScaler::maxRowRatio(const SPxLP& lp) const
{
   METHOD( "SPxScaler::maxRowRatio()" );

   Real pmax = 0.0;

   for(int i = 0; i < lp.nRows(); ++i )
   {
      const SVector& vec  = lp.rowVector(i);
      Real           mini = infinity;
      Real           maxi = 0.0;

      for(int j = 0; j < vec.size(); ++j)
      {
         Real x = fabs(vec.value(j));

         if (x < mini)
            mini = x;
         if (x > maxi)
            maxi = x;
      }
      Real p = maxi / mini;

      if (p > pmax)
         pmax = p;
   }   
   return pmax;
}

bool SPxScaler::isConsistent() const
{
#ifdef ENABLE_CONSISTENCY_CHECKS
   METHOD( "SPxScaler::isConsistent()" );

   return m_colscale.isConsistent() && m_rowscale.isConsistent();
#else
   return true;
#endif
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


