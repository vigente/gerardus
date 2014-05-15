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

/**@file  lprowset.h
 * @brief Set of LP columns.
 */
#ifndef _LPROWSET_H_
#define _LPROWSET_H_


#include <assert.h>

#include "spxdefines.h"
#include "lprow.h"
#include "dvector.h"
#include "svset.h"
#include "datakey.h"

namespace soplex
{
/**@brief   Set of LP rows.
   @ingroup Algebra
   
   Class LPRowSet implements a set of \ref LPRow "LPRows". Unless for memory
   limitations, any number of LPRows may be #add%ed to an LPRowSet. Single
   or multiple LPRows may be added to an LPRowSet, where each method
   add() comes with two different signatures. One with and one without a
   parameter, used for returning the Keys assigned to the new LPRows
   by the set. See DataKey for a more detailed description of the
   concept of keys. For the concept of renumbering LPRows within an
   LPRowSet after removal of some LPRows see DataSet.
   
   @see        DataSet, DataKey
*/
class LPRowSet : protected SVSet
{
private:

   //-----------------------------------
   /**@name Data */
   //@{
   DVector left;  ///< vector of left hand sides (lower bounds) of LPRows.
   DVector right; ///< vector of right hand sides (upper bounds) of LPRows.
   //@}

protected:

   //-----------------------------------
   /**@name Helpers */
   //@{
   /// return the complete SVSet.
   const SVSet* rowSet() const 
   {
      return this;
   }
   //@}

public:

   //-----------------------------------
   /**@name Access / modification */
   //@{
   /// returns the number of LPRows in LPRowSet.
   int num() const
   {
      return SVSet::num();
   }
   /// returns the maximum number of LPRows that fit.
   int max() const
   {
      return SVSet::max();
   }
   /// returns the vector of lhs values.
   const Vector& lhs() const
   {
      return left;
   }
   /// returns the vector of lhs values.
   Vector& lhs_w()
   {
      return left;
   }
   /// returns the lhs of the \p i 'th LPRow.
   Real lhs(int i) const
   {
      return left[i];
   }
   /// returns the lhs of the \p i 'th LPRow.
   Real& lhs_w(int i)
   {
      return left[i];
   }
   /// returns the lhs of the LPRow with DataKey \p k in LPRowSet.
   Real lhs(const DataKey& k) const
   {
      return left[number(k)];
   }
   /// returns the lhs of the LPRow with DataKey \p k in LPRowSet.
   Real& lhs_w(const DataKey& k)
   {
      return left[number(k)];
   }
   /// returns the vector of rhs values.
   const Vector& rhs() const
   {
      return right;
   }
   /// returns the vector of rhs values (writeable).
   Vector& rhs_w()
   {
      return right;
   }
   /// returns the rhs of the \p i 'th LPRow.
   Real rhs(int i) const
   {
      return right[i];
   }
   /// returns the rhs of the \p i 'th LPRow (writeable).
   Real& rhs_w(int i)
   {
      return right[i];
   }
   /// returns the rhs of the LPRow with DataKey \p k in LPRowSet.
   Real rhs(const DataKey& k) const
   {
      return right[number(k)];
   }
   /// returns the rhs of the LPRow with DataKey \p k in LPRowSet (writeable).
   Real& rhs_w(const DataKey& k)
   {
      return right[number(k)];
   }
   /// returns a writable rowVector of the \p i 'th LPRow.
   SVector& rowVector_w(int i)
   {
      return operator[](i);
   }
   /// returns the rowVector of the \p i 'th LPRow.
   const SVector& rowVector(int i) const
   {
      return operator[](i);
   }
   /// returns a writable rowVector of the LPRow with DataKey \p k.
   SVector& rowVector_w(const DataKey& k)
   {
      return operator[](k);
   }
   /// returns the rowVector of the LPRow with DataKey \p k.
   const SVector& rowVector(const DataKey& k) const
   {
      return operator[](k);
   }

   /// returns the inequalitiy type of the \p i 'th LPRow.
   LPRow::Type type(int i) const
   {
      if (rhs(i) >= infinity)
         return LPRow::GREATER_EQUAL;
      if (lhs(i) <= -infinity)
         return LPRow::LESS_EQUAL;
      if (lhs(i) == rhs(i))
         return LPRow::EQUAL;
      return LPRow::RANGE;
   }

   /// returns the inequality type of the LPRow with DataKey \p k.
   LPRow::Type type(const DataKey& k) const
   {
      return type(number(k));
   }

   /// changes the inequality type of row \p i to \p type.
   void setType(int i, LPRow::Type type);

   /// returns the value of the \p i'th LPRow.
   Real value(int i) const
   {
      if (rhs(i) < infinity)
         return rhs(i);
      else
      {
         assert(lhs(i) > -infinity);
         return lhs(i);
      }
   }

   /// returns the value of the LPRow with DataKey \p k.
   /** The \em value of a row depends on its type: if the inequality is of
       type "greater or equal", the value is the lhs of the row. Otherwise,
       the value is the rhs.
   */
   Real value(const DataKey& k) const
   {
      return value(number(k));
   }

   /// returns the DataKey of the \p i 'th LPRow in LPRowSet.
   DataKey key(int i) const
   {
      return SVSet::key(i);
   }

   /// returns the number of the LPRow with DataKey \p k in LPRowSet.
   int number(const DataKey& k) const
   {
      return SVSet::number(k);
   }

   /// does DataKey \p k belong to LPRowSet ?
   bool has(const DataKey& k) const
   {
      return SVSet::has(k);
   }
   //@}


   //-----------------------------------
   /**@name Extension
      Extension methods come with two signatures, one of them providing a
      parameter to return the assigned DataKey(s). See DataSet for a more
      detailed description. All extension methods will automatically rearrange
      or allocate more memory if required.
   */
   //@{
   ///
   void add(const LPRow& row)
   {
      DataKey k;
      add(k, row);
   }
   /// adds \p row to LPRowSet.
   void add(DataKey& pkey, const LPRow& prow)
   {
      add(pkey, prow.lhs(), prow.rowVector(), prow.rhs());
   }

   /// adds LPRow consisting of left hand side \p lhs, row vector \p rowVector, and right hand side \p rhs to LPRowSet.
   void add(Real plhs, const SVector& prowVector, Real prhs)
   {
      DataKey k;
      add(k, plhs, prowVector, prhs);
   }
   /// adds LPRow consisting of left hand side \p lhs, row vector \p
   /// rowVector, and right hand side \p rhs to LPRowSet, with DataKey
   /// \p key.
   void add(DataKey& key, Real lhs, const SVector& rowVector, Real rhs);

   ///
   void add(const LPRowSet& set);
   /// adds all LPRows of \p set to LPRowSet.
   void add(DataKey key[], const LPRowSet& set);

   /// extends row \p n to fit \p newmax nonzeros.
   void xtend(int n, int newmax)
   {
      SVSet::xtend(rowVector_w(n), newmax);
   }

   /// extend row with DataKey \p key to fit \p newmax nonzeros.
   void xtend(const DataKey& pkey, int pnewmax)
   {
      SVSet::xtend(rowVector_w(pkey), pnewmax);
   }

   /// adds \p n nonzero (\p idx, \p val)-pairs to rowVector with DataKey \p k.
   void add2(const DataKey& k, int n, const int idx[], const Real val[])
   {
      SVSet::add2(rowVector_w(k), n, idx, val);
   }

   /// adds \p n nonzero (\p idx, \p val)-pairs to \p i 'th rowVector.
   void add2(int i, int n, const int idx[], const Real val[])
   {
      SVSet::add2(rowVector_w(i), n, idx, val);
   }

   /// creates new LPRow with specified parameters and returns a reference to its row vector.
   SVector& create(int pnonzeros = 0, Real plhs = 0, Real prhs = 1)
   {
      DataKey k;
      return create(k, pnonzeros, plhs, prhs);
   }
   /// creates new LPRow with specified parameters and returns a reference to its row vector.
   SVector& create(DataKey& nkey, int nonzeros = 0, Real lhs = 0, Real rhs = 1);
   //@}


   //-----------------------------------
   /**@name Shrinking
       See DataSet for a description of the renumbering of the remaining
       LPRows in a LPRowSet after the call of a removal method.
    */
   //@{
   /// removes \p i 'th LPRow.
   void remove(int i);
   /// removes LPRow with DataKey \p k.
   void remove(const DataKey& k)
   {
      remove(number(k));
   }


   /// removes multiple LPRows.
   void remove(int perm[]);

   /// removes \p n LPRows with row numbers given by \p nums.
   void remove(const int nums[], int n)
   {
      DataArray<int> perm(num());
      remove(nums, n, perm.get_ptr());
   }
   /// removes \p n LPRows with row numbers given by \p nums, 
   /// Stores permutation of row indices in \p perm.
   void remove(const int nums[], int n, int* perm);

   /// removes all LPRows.
   void clear();
   //@}


   //-----------------------------------
   /**@name Memory Management
       For a description of the memory management methods, see the
       documentation of SVSet, which has been used for implementating
       LPRowSet.
    */
   //@{
   /// reallocates memory to be able to store \p newmax LPRows.
   void reMax(int newmax = 0)
   {
      SVSet::reMax(newmax);
      left.reSize (max());
      right.reSize(max());
   }

   /// returns number of used nonzero entries.
   int memSize() const
   {
      return SVSet::memSize();
   }

   /// returns length of nonzero memory.
   int memMax() const
   {
      return SVSet::memMax();
   }

   /// reallocates memory to be able to store \p newmax nonzeros.
   void memRemax(int newmax)
   {
      SVSet::memRemax(newmax);
   }

   /// garbage collection in nonzero memory.
   void memPack()
   {
      SVSet::memPack();
   }
   //@}

   //-----------------------------------
   /**@name Consistency check */
   /// check consistency.
   bool isConsistent() const;
   //@}

   //------------------------------------
   /**@name Construction / Destruction */
   //@{
   /// default constructor.
   /** The user can specify the initial maximum number of rows \p max
       and the initial maximum number of nonzero entries \p memmax. If these
       parameters are omitted, a default size is used. However, one can add
       an arbitrary number of rows to the LPRowSet, which may result in
       automated memory realllocation.
   */
   explicit
   LPRowSet(int pmax = -1, int pmemmax = -1)
      : SVSet(pmax, pmemmax), left(0), right(0)
   { 
      assert(isConsistent());
   }

   /// assignment operator.
   LPRowSet& operator=(const LPRowSet& rs)
   {
      if (this != &rs)
      {
         SVSet::operator=(rs);
         left  = rs.left;
         right = rs.right;

         assert(isConsistent());
      }
      return *this;
   }

   /// copy constructor.
   LPRowSet(const LPRowSet& rs)
      : SVSet( rs )
      , left ( rs.left )
      , right( rs.right )
   {
      assert(isConsistent());
   }

   /// destructor
   ~LPRowSet()
   {}
   //@}

};
} // namespace soplex
#endif // _LPROWSET_H_

//-----------------------------------------------------------------------------
//Emacs Local Variables:
//Emacs mode:c++
//Emacs c-basic-offset:3
//Emacs tab-width:8
//Emacs indent-tabs-mode:nil
//Emacs End:
//-----------------------------------------------------------------------------
