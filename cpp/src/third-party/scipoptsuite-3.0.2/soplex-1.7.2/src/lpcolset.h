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

/**@file  lpcolset.h
 * @brief Set of LP columns.
 */
#ifndef _LPCOLSET_H_
#define _LPCOLSET_H_

#include <assert.h>

#include "spxdefines.h"
#include "lpcol.h"
#include "dvector.h"
#include "svset.h"
#include "datakey.h"

namespace soplex
{
/**@brief   Set of LP columns.
   @ingroup Algebra

   Class LPColSet implements a set of \ref LPCol "LPCols". Unless for
   memory limitations, any number of LPCols may be #add%ed to an
   LPColSet. Single or multiple LPCols may be #add%ed to an LPColSet,
   where each method add() comes with two different signatures. One with
   and one without a parameter, used for returning the \ref DataKey
   "DataKeys" assigned to the new LPCols by the set. See DataKey for a
   more detailed description of the concept of keys. For the concept of
   renumbering LPCols within an LPColSet after removal of some LPCols,
   see DataSet.

   @see        DataSet, DataKey
 */
class LPColSet : protected SVSet
{
private:

   //------------------------------------
   /**@name Data */
   //@{
   DVector low;     ///< vector of lower bounds.
   DVector up;      ///< vector of upper bounds.
   DVector object;  ///< vector of objective coefficients.
   //@}

protected:

   //------------------------------------
   /**@name Protected helpers */
   //@{
   /// return the complete SVSet.
   const SVSet* colSet() const 
   {
      return this;
   }
   //@}

public:

   //------------------------------------
   /**@name Inquiry */
   //@{
   /// returns the number of LPCols currently in LPColSet.
   int num() const
   {
      return SVSet::num();
   }

   /// returns maximum number of LPCols currently fitting into LPColSet.
   int max() const
   {
      return SVSet::max();
   }

   ///
   const Vector& maxObj() const
   {
      return object;
   }
   /// returns vector of objective values w.r.t. maximization.
   Vector& maxObj_w()
   {
      return object;
   }
   ///
   Real maxObj(int i) const
   {
      return object[i];
   }
   /// returns objective value (w.r.t. maximization) of \p i 'th LPCol in LPColSet.
   Real& maxObj_w(int i)
   {
      return object[i];
   }
   ///
   Real maxObj(const DataKey& k) const
   {
      return object[number(k)];
   }
   /// returns objective value (w.r.t. maximization) of LPCol with DataKey \p k in LPColSet.
   Real& maxObj_w(const DataKey& k)
   {
      return object[number(k)];
   }
   ///
   const Vector& lower() const
   {
      return low;
   }
   /// returns vector of lower bound values.
   Vector& lower_w()
   {
      return low;
   }
   ///
   Real lower(int i) const
   {
      return low[i];
   }
   /// returns lower bound of \p i 'th LPCol in LPColSet.
   Real& lower_w(int i)
   {
      return low[i];
   }
   ///
   Real lower(const DataKey& k) const
   {
      return low[number(k)];
   }
   /// returns lower bound of LPCol with DataKey \p k in LPColSet.
   Real& lower_w(const DataKey& k)
   {
      return low[number(k)];
   }
   ///
   const Vector& upper() const
   {
      return up;
   }
   /// returns vector of upper bound values.
   Vector& upper_w()
   {
      return up;
   }
   ///
   Real upper(int i) const
   {
      return up[i];
   }
   /// returns upper bound of \p i 'th LPCol in LPColSet.
   Real& upper_w(int i)
   {
      return up[i];
   }
   ///
   Real upper(const DataKey& k) const
   {
      return up[number(k)];
   }
   /// returns upper bound of LPCol with DataKey \p k in LPColSet.
   Real& upper_w(const DataKey& k)
   {
      return up[number(k)];
   }
   ///
   SVector& colVector_w(int i)
   {
      return operator[](i);
   }
   /// returns colVector of \p i 'th LPCol in LPColSet.
   const SVector& colVector(int i) const
   {
      return operator[](i);
   }
   /// returns writeable colVector of LPCol with DataKey \p k in LPColSet.
   SVector& colVector_w(const DataKey& k)
   {
      return operator[](k);
   }
   /// returns colVector of LPCol with DataKey \p k in LPColSet.
   const SVector& colVector(const DataKey& k) const
   {
      return operator[](k);
   }

   /// returns DataKey of \p i 'th LPCol in LPColSet.
   DataKey key(int i) const
   {
      return SVSet::key(i);
   }

   /// returns number of LPCol with DataKey \p k in LPColSet.
   int number(const DataKey& k) const
   {
      return SVSet::number(k);
   }

   /// does DataKey \p k belong to LPColSet ?
   bool has(const DataKey& k) const
   {
      return SVSet::has(k);
   }
   //@}


   //------------------------------------
   /**@name Extension
      All extension methods come with two signatures, one of which providing a
      parameter to return the assigned DataKey(s). See DataSet for a more
      detailed description. All extension methods are designed to
      automatically reallocate memory if required.
   */
   //@{
   ///
   void add(const LPCol& pcol)
   {
      DataKey k;
      add(k, pcol);
   }
   /// adds p pcol to LPColSet.
   void add(DataKey& pkey, const LPCol& pcol)
   {
      add(pkey, pcol.obj(), pcol.lower(),
           pcol.colVector(), pcol.upper());
   }

   ///
   void add(Real pobj, Real plower, const SVector& pcolVector, Real pupper)
   {
      DataKey k;
      add(k, pobj, plower, pcolVector, pupper);
   }
   /// adds LPCol consisting of objective value \p obj, lower bound \p lower, column vector \p colVector and upper bound \p upper to LPColSet.
   void add (DataKey& key,
             Real obj,
             Real lower,
             const SVector& colVector,
             Real upper);

   ///
   void add(const LPColSet& set);
   /// adds all LPCols of \p set to LPColSet.
   void add(DataKey key[], const LPColSet& set);

   /// extends column \p n to fit \p newmax nonzeros.
   void xtend(int n, int newmax)
   {
      SVSet::xtend(colVector_w(n), newmax);
   }

   /// extend column with DataKey \p key to fit \p newmax nonzeros.
   void xtend(const DataKey& pkey, int pnewmax)
   {
      SVSet::xtend(colVector_w(pkey), pnewmax);
   }

   ///
   void add2(const DataKey& k, int n, const int idx[], const Real val[])
   {
      SVSet::add2(colVector_w(k), n, idx, val);
   }
   /// adds \p n nonzero (\p idx, \p val)-pairs to \p i 'th colVector.
   void add2(int i, int n, const int idx[], const Real val[])
   {
      SVSet::add2(colVector_w(i), n, idx, val);
   }

   ///
   SVector& create(int pnonzeros = 0, Real pobj = 1, Real plw = 0, Real pupp = 1)
   {
      DataKey k;
      return create(k, pnonzeros, pobj, plw, pupp);
   }
   /// creates new LPCol with specified arguments and returns a reference to its column vector.
   SVector& create(DataKey& nkey, int nonzeros = 0, Real obj = 1, Real low = 0, Real up = 1);
   //@}


   //------------------------------------
   /**@name Shrinking
      See DataSet for a description of the renumbering of the remaining
      LPCols in a LPColSet after the call of a removal method.
   */
   //@{
   /// removes \p i 'th LPCol.
   void remove(int i);

   /// removes LPCol with DataKey \p k.
   void remove(const DataKey& k)
   {
      remove(number(k));
   }

   /// removes multiple elements.
   void remove(int perm[]);

   /// removes LPCols with numbers \p nums, where \p n is the length of the array \p nums
   void remove(const int nums[], int n)
   {
      DataArray < int > perm(num());
      remove(nums, n, perm.get_ptr());
   }
   /// removes LPCols with numbers \p nums, where \p n is the length of the array \p nums, and stores the index permutation in array \p perm.
   void remove(const int nums[], int n, int* perm);

   /// removes all LPCols from the set.
   void clear();
   //@}


   //------------------------------------
   /**@name Memory Management
      See SVSet for a description of the memory management methods.
   */
   //@{
   /// reallocates memory to be able to store \p newmax LPCols.
   void reMax(int newmax = 0)
   {
      SVSet::reMax (newmax);
      up.reSize (max());
      low.reSize (max());
      object.reSize(max());
   }

   /// returns used nonzero memory.
   int memSize() const
   {
      return SVSet::memSize();
   }

   /// returns length of nonzero memory.
   int memMax() const
   {
      return SVSet::memMax();
   }

   /// resets length of nonzero memory.
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

   //------------------------------------
   /**@name Miscellaneous */
   //@{
   ///
   bool isConsistent() const;
   //@}

   //------------------------------------
   /**@name Constructors / Destructors */
   //@{
   /// default constructor.
   /** The user can specify the initial maximum number of columns \p max
       and the initial maximum number of nonzero entries \p memmax. If these
       parameters are omitted, a default size is used. However, one can add
       an arbitrary number of columns to the LPColSet, which may result in
       automated memory realllocation.
   */
   explicit
   LPColSet(int pmax = -1, int pmemmax = -1)
      : SVSet(pmax, pmemmax), low(0), up(0), object(0)
   { 
      assert(isConsistent());
   }

   /// assignment operator.
   LPColSet& operator=(const LPColSet& rs)
   {
      if (this != &rs)
      {
         SVSet::operator=(rs);
         low = rs.low;
         up = rs.up;
         object = rs.object;

         assert(isConsistent());
      }
      return *this;
   }
   /// copy constructor.
   LPColSet(const LPColSet& rs)
      : SVSet ( rs )
      , low   ( rs.low )
      , up    ( rs.up )
      , object( rs.object )
   {
      assert(isConsistent());
   }

   /// destructor
   ~LPColSet()
   {}
   //@}
};


} // namespace soplex
#endif // _LPCOLSET_H_

//-----------------------------------------------------------------------------
//Emacs Local Variables:
//Emacs mode:c++
//Emacs c-basic-offset:3
//Emacs tab-width:8
//Emacs indent-tabs-mode:nil
//Emacs End:
//-----------------------------------------------------------------------------
