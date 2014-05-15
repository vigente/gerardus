/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the class library                   */
/*       SoPlex --- the Sequential object-oriented simPlex.                  */
/*                                                                           */
/*    Copyright (C) 1996      Roland Wunderling                              */
/*                  1996-2013 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SoPlex is distributed under the terms of the ZIB Academic Licence.       */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SoPlex; see the file COPYING. If not email to soplex@zib.de.  */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file  vector_exact.h
 * @brief Dense vector of MpqReal.
 */
#ifndef _VECTOR_EXACT_H_
#define _VECTOR_EXACT_H_

#include <assert.h>
#include <string.h>
#include <math.h>
#include <iostream>

#include "mpqreal.h"

namespace soplex
{
class Vector;
class SVector;
class SSVector;

/**@brief   Dense vector of MpqReal.
 * @ingroup Algebra
 *
 * Class Vector_exact is a copy of class Vector replacing the floating point type Real with the exact MpqReal.
 */
class Vector_exact
{
protected:

   //------------------------------------
   /**@name Data */
   //@{
   /// dimension of vector
   int dimen;
   /// values of a vector
   /** The memory block pointed to by val must at least have size
    *  dimen * sizeof(Real)
    */
   MpqReal* val;
   //@}

public:

   //------------------------------------
   /**@name Construction and assignment */
   //@{
   /// construction
   /** There is no default constructor since the storage for a
    *  Vector must be provided externally.
    *  Storage must be passed as a memory block val at construction. It
    *  must be large enough to fit at least dimen Real values.
    */
   Vector_exact(int p_dimen, MpqReal *p_val)
      : dimen(p_dimen)
      , val(p_val)
   {
      assert(dimen >= 0);
      assert(isConsistent());
   }
   /// Assignment operator.
   Vector_exact& operator=(const Vector_exact& vec);

   /// Assignment operator.
   Vector_exact& operator=(const Vector& vec);

   /// Assignment operator.
   /** Assigning a SVector to a Vector using operator=()
    *  will set all values to 0 except the nonzeros of \p vec.
    *  This is diffent in method assign().
    */
   Vector_exact& operator=(const SVector& vec);
   //@}

   //------------------------------------
   /**@name Access */
   //@{
   /// dimension of vector
   inline int dim() const
   {
      return dimen;
   }
   /// return \p n 'th value by reference
   MpqReal& operator[](int n)
   {
      assert(n >= 0 && n < dimen);
      return val[n];
   }

   /// return \p n 'th value
   MpqReal operator[](int n) const
   {
      assert(n >= 0 && n < dimen);
      return val[n];
   }

   /// equality operator
   friend bool operator==(const Vector_exact& vec1, const Vector_exact& vec2);
   //@}

   //------------------------------------
   /**@name Algebraic methods */
   //@{
   /// vector addition
   Vector_exact& operator+=(const Vector_exact& vec);
   /// vector addition
   Vector_exact& operator+=(const Vector& vec);
   /// vector addition
   Vector_exact& operator+=(const SVector& vec);

   /// vector difference
   Vector_exact& operator-=(const Vector_exact& vec);
   /// vector difference
   Vector_exact& operator-=(const Vector& vec);
   /// vector difference
   Vector_exact& operator-=(const SVector& vec);

   /// scaling
   Vector_exact& operator*=(MpqReal x);

   /// absolute biggest element (infinity norm).
   MpqReal maxAbs() const;
   /// absolute smallest element.
   MpqReal minAbs() const;

   /// addition of scaled vector
   Vector_exact& multAdd(MpqReal x, const Vector_exact& vec);
   ///  addition of scaled vector
   Vector_exact& multAdd(MpqReal x, const Vector& vec);
   /// addition of scaled vector
   Vector_exact& multAdd(MpqReal x, const SVector& vec);
   //@}

   //------------------------------------
   /**@name Utilities */
   //@{
   /// Conversion to C-style pointer.
   /** This function serves for using a Vector in an C-style
    *  function. It returns a pointer to the first value of the array.
    *
    *  @todo check whether this non-const c-style acces should indeed be public
    */
   MpqReal* get_ptr()
   {
      return val;
   }
   /// Conversion to C-style pointer.
   /** This function serves for using a Vector in an C-style
    *  function. It returns a pointer to the first value of the array.
    */
   const MpqReal* get_const_ptr() const
   {
      return val;
   }
   /// output operator.
   friend std::ostream& operator<<(std::ostream& s, const Vector_exact& vec);

   /// consistency check.
   bool isConsistent() const;

   /// set vector to 0.
   void clear()
   {
      if( dimen > 0 )
      {
         for( int i = 0; i < dimen; i++ )
            val[i] = 0;
      }
   }
   //@}

private:

   //------------------------------------
   /**@name Blocked */
   //@{
   /// we have no default constructor.
   Vector_exact();
   //@}
};
} // namespace soplex
#endif // _VECTOR_EXACT_H_

//-----------------------------------------------------------------------------
//Emacs Local Variables:
//Emacs mode:c++
//Emacs c-basic-offset:3
//Emacs tab-width:8
//Emacs indent-tabs-mode:nil
//Emacs End:
//-----------------------------------------------------------------------------
