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

/**@file  dvector.h
 * @brief Dynamic vectors.
 */
#ifndef _DVECTOR_H_
#define _DVECTOR_H_

#include <iostream>
#include <assert.h>

#include "spxdefines.h"
#include "vector.h"
#include "svector.h"

namespace soplex
{
class Vector_exact;

/**@brief   Dynamic vectors.
   @ingroup Algebra

   Class DVector is a derived class of Vector adding automatic
   memory management to such objects. This allows to implement maths
   operations operator+() and operator-(). Further, it is possible
   to reset the dimension of a DVector via method reDim(). However,
   this may render all references to values of a #reDim()%ed DVector
   invalid.
   
   For vectors that are often subject to reDim() it may be
   unconvenient to reallocate the required memory every
   time. Instead, an array of values of length memSize() is kept,
   where only the first dim() elements are used.  Initially,
   memSize() == dim(). However, if the dimension is increased,
   memSize() will be increased at least by a factor of 1.2 to be
   prepared for future (small) #reDim()%s. Finally, one can explicitly
   set memSize() with method reSize(), but not lower than
   dim().  
*/
class DVector : public Vector
{
   //-----------------------------------
   /**@name Data */
   //@{
   int   memsize;        ///< length of array of values \ref soplex::DVector::mem "mem"
   Real* mem;            ///< value array to be used
   //@}

public:

   //-----------------------------------------
   /**@name Arithmetic operations (friends) */
   //@{
   /// adding vectors.
   friend DVector operator+(const Vector& v, const Vector& w);
   /// adding vectors.
   friend DVector operator+(const SVector& v, const Vector& w);
   /// adding vectors.
   friend DVector operator+(const Vector& v, const SVector& w);

   /// subtracting vectors.
   friend DVector operator-(const Vector& v, const Vector& w);
   /// subtracting vectors.
   friend DVector operator-(const SVector& v, const Vector& w);
   /// subtracting vectors.
   friend DVector operator-(const Vector& v, const SVector& w);

   /// negation operator.
   friend DVector operator-(const Vector& vec);
   /// negation operator.
   friend DVector operator-(const SVector& vec);

   /// scaling vectors with a real number.
   friend DVector operator*(const Vector& v, Real x);
   /// scaling vectors with a real number.
   friend DVector operator*(Real x, const Vector& v);
   //@}

   //-----------------------------------
   /**@name I/O */
   //@{
   /// output operator.
   friend std::istream& operator>>(std::istream& s, DVector& vec);
   //@}

#if 0
   //--------------------------------------------------
   /**@name Arithmetic operations (member functions) */
   //@{
   /**@todo Do we really have to reimplement the following operators 
    *       (inheritance from Vector?). All the more as these methods are
    *       not virtual.
    */
   /// adds \p vec to %vector.
   DVector& operator+=(const Vector& vec)
   {
      Vector::operator+=(vec);
      return *this;
   }
   /// adds \p vec to %vector.
   DVector& operator+=(const SVector& vec)
   {
      Vector::operator+=(vec);
      return *this;
   }

   /// subtracts \p vec from %vector.
   DVector& operator-=(const Vector& vec)
   {
      Vector::operator-=(vec);
      return *this;
   }
   /// subtracts \p vec from %vector.
   DVector& operator-=(const SVector& vec)
   {
      Vector::operator-=(vec);
      return *this;
   }

   /// scales vector with factor \p x
   DVector& operator*=(Real x)
   {
      Vector::operator*=(x);
      return *this;
   }
   //@}
#endif

   //--------------------------------------------------
   /**@name Access */
   //@{
   /// resets  \ref soplex::DVector "DVector"'s dimension to \p newdim.
   void reDim(int newdim);

   /// resets  \ref soplex::DVector "DVector"'s memory size to \p newsize.
   void reSize(int newsize);

   /// resets  \ref soplex::DVector "DVector"'s memory size to \p newsize and dimension to \p newdim.
   void reSize(int newsize, int newdim);

   /// returns \ref soplex::DVector "DVector"'s memory size.
   int memSize() const
   {
      return memsize;
   }

   /// consistency check.
   bool isConsistent() const;
   //@}

   //--------------------------------------------------
   /**@name Construction / destruction */
   //@{
   /// default constructor. \p dim is the initial dimension.
   explicit
   DVector(int dim = 0);
   /// copy constructor.
   explicit DVector(const Vector& old);
   /// cast from exact vector.
   explicit DVector(const Vector_exact& old);
   /// copy constructor.
   DVector(const DVector& old);
   /// assignment operator.
   DVector& operator=(const Vector& vec)
   {
      if (vec.dim() != dim())
         reDim(vec.dim());
      Vector::operator=(vec);
      
      assert(isConsistent());

      return *this;
   }
   /// assignment operator.
   DVector& operator=(const DVector& vec)
   {
      if (this != &vec)
      {
         if (vec.dim() != dim())
            reDim(vec.dim());
         Vector::operator=(vec);

         assert(isConsistent());
      }
      return *this;
   }
   /// assingment operator.
   DVector& operator=(const SVector& vec)
   {
      if (vec.dim() != dim())
         reDim(vec.dim()); 
      Vector::operator=(vec);

      assert(isConsistent());

      return *this;
   }

   /// destructor.
   ~DVector();
};
} // namespace soplex
#endif // _DVECTOR_H_

//-----------------------------------------------------------------------------
//Emacs Local Variables:
//Emacs mode:c++
//Emacs c-basic-offset:3
//Emacs tab-width:8
//Emacs indent-tabs-mode:nil
//Emacs End:
//-----------------------------------------------------------------------------
