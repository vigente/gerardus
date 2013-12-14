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

/**@file  vector.h
 * @brief Dense vector for linear algebra.
 */
#ifndef _VECTOR_H_
#define _VECTOR_H_

#include <assert.h>
#include <string.h>
#include <math.h>
#include <iostream>

#include "spxdefines.h"

namespace soplex
{
class SVector;
class SSVector;
class Vector_exact;

/**@brief   Dense vector
   @ingroup Algebra
 
   Class Vector provides dense linear algebra vectors. It does not
   provide memory management for the %array of values. Instead, the
   constructor requires a pointer to a memory block large enough to
   fit the desired dimension of Real values.
 
   After construction, the values of a Vector can be accessed with
   the subscript operator[]() .  Safety is provided by
   - checking of array bound when accessing elements with the
     subscript operator[]() (only when compiled without \c -DNDEBUG).
 
   A Vector is distinguished from a simple array of #Real%s by
   providing a set of mathematical operations. Since Vector does
   not provide any memory management features, no operations are
   available that would require allocation of temporary memory
   space.

   The following mathematical operations are provided by class Vector
   (Vector \p a, \p b; Real \p x): 
 
   <TABLE>
   <TR><TD>Operation</TD><TD>Description   </TD><TD></TD>&nbsp;</TR>
   <TR><TD>\c -=    </TD><TD>subtraction   </TD><TD>\c a \c -= \c b </TD></TR>
   <TR><TD>\c +=    </TD><TD>addition      </TD><TD>\c a \c += \c b </TD></TR>
   <TR><TD>\c *     </TD><TD>scalar product</TD>
       <TD>\c x = \c a \c * \c b </TD></TR>
   <TR><TD>\c *=    </TD><TD>scaling       </TD><TD>\c a \c *= \c x </TD></TR>
   <TR><TD>maxAbs() </TD><TD>infinity norm </TD>
       <TD>\c a.maxAbs() == \f$\|a\|_{\infty}\f$ </TD></TR>
   <TR><TD>minAbs() </TD><TD> </TD>
       <TD>\c a.minAbs() == \f$\min |a_i|\f$ </TD></TR>

   <TR><TD>length() </TD><TD>euclidian norm</TD>
       <TD>\c a.length() == \f$\sqrt{a^2}\f$ </TD></TR>
   <TR><TD>length2()</TD><TD>square norm   </TD>
       <TD>\c a.length2() == \f$a^2\f$ </TD></TR>
   <TR><TD>multAdd(\c x,\c b)</TD><TD>add scaled vector</TD>
       <TD> \c a +=  \c x * \c b </TD></TR>
   </TABLE>
 
   When using any of these operations, the vectors involved must be of
   the same dimension. Also an SVector \c b is allowed if it
   does not contain nonzeros with index greater than the dimension of
   \c a.
*/
class Vector
{
protected:

   //------------------------------------
   /**@name Data */
   //@{
   /// dimension of vector
   int dimen;
   /// values of a vector
   /** The memory block pointed to by val must at least have size
    *  dimen * sizeof(Real).  
    */
   Real* val;
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
   Vector(int p_dimen, Real *p_val)
      : dimen(p_dimen)
      , val(p_val)
   {
      assert(dimen >= 0);
      assert(isConsistent());
   }
   /// Assignment operator.
   Vector& operator=(const Vector& vec);

   /// Assignment operator.
   /** Assigning a SVector to a Vector using operator=()
    *  will set all values to 0 except the nonzeros of \p vec. 
    *  This is diffent in method assign().
    */
   Vector& operator=(const SVector& vec);

   /// Assignment operator.
   /** Assigning a SSVector to a Vector using operator=()
    *  will set all values to 0 except the nonzeros of \p vec. 
    *  This is diffent in method assign().
    */
   Vector& operator=(const SSVector& vec);

   /// Assignment operator.
   Vector& operator=(const Vector_exact& vec);

   /// Assign values of \p sv.
   /** Assigns all nonzeros of \p sv to the vector. 
    *  All other values remain unchanged.
   */
   Vector& assign(const SVector& sv);

   /// Assign values of \p sv.
   /** Assigns all nonzeros of \p sv to the vector. 
    *  All other values remain unchanged.
    */
   Vector& assign(const SSVector& sv);
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
   Real& operator[](int n)
   {
      assert(n >= 0 && n < dimen);
      return val[n];
   }

   /// return \p n 'th value
   Real operator[](int n) const
   {
      assert(n >= 0 && n < dimen);
      return val[n];
   }
   //@}

   //------------------------------------
   /**@name Algebraic methods */
   //@{   
   /// vector addition
   Vector& operator+=(const Vector& vec);
   /// vector addition
   Vector& operator+=(const SVector& vec);
   /// vector addition
   Vector& operator+=(const SSVector& vec);

   /// vector difference
   Vector& operator-=(const Vector& vec);
   /// vector difference
   Vector& operator-=(const SVector& vec);
   /// vector difference
   Vector& operator-=(const SSVector& vec);

   /// scaling
   Vector& operator*=(Real x);

   /// inner product.
   Real operator*(const SSVector& v) const;
   /// inner product.
   Real operator*(const SVector& v) const;
   /// inner product.
   Real operator*(const Vector& v) const
   {
      assert(v.dim() == dimen);
      Real x = 0;
      for(int i = 0; i < dimen; i++)
         x += val[i] * v.val[i];
      return x;
   }

   /// absolute biggest element (infinity norm). 
   Real maxAbs() const;
   /// absolute smallest element.
   Real minAbs() const;
   /// euclidian norm.
   Real length() const;
   /// squared norm.
   Real length2() const;

   /// addition of scaled vector
   Vector& multAdd(Real x, const SVector& vec);
   /// addition of scaled vector
   Vector& multAdd(Real x, const SSVector& svec);
   ///  addition of scaled vector
   Vector& multAdd(Real x, const Vector& vec)
   {
      assert(vec.dim() == dimen);

      for(int i = 0; i < dimen; i++)
         val[i] += x * vec.val[i];

      return *this;
   }
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
   Real* get_ptr()
   {
      return val;
   }
   /// Conversion to C-style pointer.
   /** This function serves for using a Vector in an C-style
    *  function. It returns a pointer to the first value of the array.
    */
   const Real* get_const_ptr() const
   {
      return val;
   }
   /// output operator.
   friend std::ostream& operator<<(std::ostream& s, const Vector& vec);

   /// consistency check.
   bool isConsistent() const;

   /// set vector to 0.
   void clear()
   {
      if (dimen > 0)
         memset(val, 0, dimen * sizeof(Real));
   }
   //@}

private:

   //------------------------------------
   /**@name Blocked */
   //@{
   /// we have no default constructor.
   Vector();
   //@}
};
} // namespace soplex
#endif // _VECTOR_H_

//-----------------------------------------------------------------------------
//Emacs Local Variables:
//Emacs mode:c++
//Emacs c-basic-offset:3
//Emacs tab-width:8
//Emacs indent-tabs-mode:nil
//Emacs End:
//-----------------------------------------------------------------------------
