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

/**@file  dvector_exact.h
 * @brief Dynamic vectors of MpqReal.
 */
#ifndef _DVECTOR_EXACT_H_
#define _DVECTOR_EXACT_H_

#include <iostream>
#include <assert.h>

#include "spxdefines.h"
#include "mpqreal.h"
#include "vector_exact.h"
#include "vector.h"
#include "svector.h"

namespace soplex
{
/**@brief   Dynamic vectors of MpqReal.
 * @ingroup Algebra
 *
 * Class DVector_exact is a copy of class DVector replacing the floating point type Real with the exact MpqReal.
 */
class DVector_exact : public Vector_exact
{
   //-----------------------------------
   /**@name Data */
   //@{
   int   memsize;        ///< length of array of values \ref soplex::DVector::mem "mem"
   MpqReal* mem;         ///< value array to be used
   //@}

public:

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
   DVector_exact(int dim = 0);
   /// copy constructor.
   explicit DVector_exact(const Vector& old);
   /// copy constructor.
   explicit DVector_exact(const Vector_exact& old);
   /// copy constructor.
   DVector_exact(const DVector_exact& old);
   /// assignment operator.
   DVector_exact& operator=(const DVector_exact& vec)
   {
      if( this != &vec )
      {
         if ( vec.dim() != dim() )
            reDim(vec.dim());
         Vector_exact::operator=(vec);

         assert(isConsistent());
      }

      return *this;
   }
   /// assignment operator.
   DVector_exact& operator=(const Vector_exact& vec)
   {
      if( vec.dim() != dim() )
         reDim(vec.dim());
      Vector_exact::operator=(vec);

      assert(isConsistent());

      return *this;
   }
   /// assignment operator.
   DVector_exact& operator=(const Vector& vec)
   {
      if( vec.dim() != dim() )
         reDim(vec.dim());
      Vector_exact::operator=(vec);

      assert(isConsistent());

      return *this;
   }
   /// assignment operator.
   DVector_exact& operator=(const SVector& vec)
   {
      if( vec.dim() != dim() )
         reDim(vec.dim());
      Vector_exact::operator=(vec);

      assert(isConsistent());

      return *this;
   }

   /// destructor.
   ~DVector_exact();
};
} // namespace soplex
#endif // _DVECTOR_EXACT_H_

//-----------------------------------------------------------------------------
//Emacs Local Variables:
//Emacs mode:c++
//Emacs c-basic-offset:3
//Emacs tab-width:8
//Emacs indent-tabs-mode:nil
//Emacs End:
//-----------------------------------------------------------------------------
