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

/**@file  lprow.h
 * @brief (In)equality for LPs.
 */
#ifndef _LPROW_H_
#define _LPROW_H_

#include <assert.h>

#include "spxdefines.h"
#include "dsvector.h"

namespace soplex
{
/**@brief   (In)equality for LPs.
   @ingroup Algo

   Class LPRow provides constraints for linear programs in the form
   \f[
                       l \le a^Tx \le r,
   \f]
   where \em a is a DSVector. \em l is referred to as 
   %left hand side,
   \em r as %right hand side and \em a as \em row \em vector or 
   the constraint vector. \em l and \em r may also take values 
   \f$\pm\f$ #infinity. 
   This static member is predefined, but may be overridden to meet 
   the needs of the LP solver to be used.
 
   LPRows allow to specify regular inequalities of the form 
   \f[
                           a^Tx \sim \alpha,
   \f]
   where \f$\sim\f$ can take any value 
   of \f$\le, =, \ge\f$, by setting rhs and
   lhs to the same value or setting one of them to \f$\infty\f$.

   Since constraints in the regular form occur often, LPRows offers methods
   type() and value() for retreiving \f$\sim\f$ and \f$\alpha\f$ of 
   an LPRow in this form, respectively. Also, a constructor for 
   LPRows given in regular form is provided.
*/
class LPRow
{
private:

   //------------------------------------
   /**@name Data */
   //@{
   Real   left;      ///< left-hand side of the constraint
   Real   right;     ///< right-hand side of the constraint
   DSVector vec;     ///< the row vector
   //@}

public:

   //------------------------------------
   /**@name Types */
   //@{
   /// (In)Equality of an LP row.
   /** #LPRow%s may be of one of the above Types. This datatype may be
    *  used for constructing new #LPRow%s in the regular form.
    */
   enum Type
   {                          
      LESS_EQUAL,          ///< \f$a^Tx \le \alpha\f$.   
      EQUAL,               ///< \f$a^Tx = \alpha\f$.   
      GREATER_EQUAL,       ///< \f$a^Tx \ge \alpha\f$.    
      RANGE                ///< \f$\lambda \le a^Tx \le \rho\f$.
   };
   //@}

   //------------------------------------
   /**@name Construction / destruction */
   //@{
   /// Construct LPRow with a vector ready to hold \p defDim nonzeros
   explicit LPRow(int defDim = 0) 
      : left(0), right(infinity), vec(defDim)
   {
      assert(isConsistent());
   }

   /// copy constructor
   LPRow(const LPRow& row) 
      : left(row.left), right(row.right), vec(row.vec)
   {
      assert(isConsistent());
   }

   /// Construct LPRow with the given left-hand side, right-hand side
   /// and rowVector.
   LPRow(Real p_lhs, const SVector& p_rowVector, Real p_rhs)
      : left(p_lhs), right(p_rhs), vec(p_rowVector)
   {
      assert(isConsistent());
   }

   /// Construct LPRow from passed \p rowVector, \p type and \p value
   LPRow(const SVector& p_rowVector, Type p_type, Real p_value);

   /// destructor
   ~LPRow()
   {}
   //@}

   //------------------------------------
   /**@name Access / modification */
   //@{
   /// get type of row.
   Type type() const;

   /// set type of (in)equality
   void setType(Type p_type);

   /// Right hand side value of (in)equality.
   /** This method returns \f$\alpha\f$ for a LPRow in regular form.
    *  However, value() may only be called for LPRow%s with
    *  type() != \c RANGE.
    */
   Real value() const;

   /// get left-hand side value.
   Real lhs() const
   {
      return left;
   }

   /// access left-hand side value.
   void setLhs(Real p_left)
   {
      left = p_left;
   }

   /// get right-hand side value.
   Real rhs() const
   {
      return right;
   }

   /// access right-hand side value.
   void setRhs(Real p_right)
   {
      right = p_right;
   }

   /// get constraint row vector
   const SVector& rowVector() const
   {
      return vec;
   }

   /// access constraint row vector.
   void setRowVector(const DSVector& p_vec)
   {
      vec = p_vec;
   }
   //@}

   //------------------------------------
   /**@name Consistency check */
   //@{
   /// check consistency.
   bool isConsistent() const
   {
#ifdef ENABLE_CONSISTENCY_CHECKS
      return vec.isConsistent();
#else
      return true;
#endif
   }
   //@}
};

} // namespace soplex
#endif // _LPROW_H_

//-----------------------------------------------------------------------------
//Emacs Local Variables:
//Emacs mode:c++
//Emacs c-basic-offset:3
//Emacs tab-width:8
//Emacs indent-tabs-mode:nil
//Emacs End:
//-----------------------------------------------------------------------------
