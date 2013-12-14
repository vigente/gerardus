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

/**@file  spxproof.cpp
 * @brief provable bounds
 */

#ifdef WITH_EASYVAL

#include <assert.h>
#include <iostream>

#include "Easyval.hh"

#include "spxdefines.h"
#include "spxsolver.h"

//#undef DEBUG
//#define DEBUG(x) {x}

namespace soplex
{

Real SPxSolver::provedBound(Vector& dualsol, const Vector& objvec) const
{
   DVector lhsvec(dim());
   DVector rhsvec(dim());
   Real eps = opttol();
   Easyval ytb;
   Easyval scalprod;
   int sen = (spxSense() == MINIMIZE ? +1 : -1);

   assert(rep() == COLUMN);

   /* for minimization problems: calculate y^Tb + min{(c^T - y^TA)x | l <= x <= u}
    * for maximization problems: calculate y^Tb + max{(c^T - y^TA)x | l <= x <= u}
    */

   getLhs(lhsvec);
   getRhs(rhsvec);

   /* 1. calculate y^Tb */
   ytb = 0;
   for( int i = 0; i < dim(); ++i )
   {
      Easyval b;

      if( sen * dualsol[i] > eps )
         b = lhsvec[i];
      else if( sen * dualsol[i] < -eps )
         b = rhsvec[i];
      else
         dualsol[i] = 0.0;

      Easyval y(dualsol[i]);

      ytb += y*b;
      MSG_DEBUG( spxout << "DPROOF02 y[" << i << "] = " << dualsol[i]
                        << ", lhs[" << i << "] = " << lhsvec[i]
                        << ", rhs[" << i << "] = " << rhsvec[i]
                        << " -> ytb = " << ytb << std::endl; )
   }

   /* 2. calculate min/max{(c^T - y^TA)x} */
   scalprod = 0;
   for( int j = 0; j < coDim(); ++j )
   {
      LPCol col;
      Easyval x(lower(j), upper(j)); /* add/sub machine epsilon? */

      getCol(j, col);
      SVector colvec = col.colVector();  /* make this faster !!! */
      Easyval diff = objvec[j];
      MSG_DEBUG( spxout << "DPROOF03 obj[" << j << "] = " << objvec[j]
                        << " -> diff = " << diff << std::endl; )

      for( int i = 0; i < colvec.size(); ++i )
      {
         Easyval y(dualsol[colvec.index(i)]);
         Easyval a(colvec.value(i));
         diff -= y*a;
         MSG_DEBUG( spxout << "DPROOF04 y[" << colvec.index(i) << "] = " 
                           << dualsol[colvec.index(i)]
                           << ", a[" << colvec.index(i) << "] = "
                           << colvec.value(i)
                           << " -> diff = " << diff << std::endl; )
      }

      diff *= x;
      scalprod += diff;
      MSG_DEBUG( spxout << "DPROOF05 x[" << j << "] = " << x 
                        << " -> diff = " << diff << std::endl
                        << "   ------------> scalprod = " 
                        << scalprod << std::endl; )
   }

   /* add y^Tb */
   scalprod += ytb;
   MSG_INFO1( spxout << "IPROOF01 proved bound = " << scalprod << std::endl; );

   /* depending on the objective sense, choose min or max */
   if( spxSense() == MINIMIZE )
      return scalprod.getInf();
   else
      return scalprod.getSup();
}

Real SPxSolver::provedDualbound() const
{
   DVector dualsol(dim());
   DVector objvec(coDim());

   getDual(dualsol);
   getObj(objvec);

   return provedBound(dualsol, objvec);
}

bool SPxSolver::isProvenInfeasible() const
{
   DVector dualfarkas(dim());
   DVector zerovec(coDim());

   getDualfarkas(dualfarkas);
   zerovec.clear();

   return (provedBound(dualfarkas, zerovec) > 0.0);
}

}

#endif // WITH_EASYVAL
