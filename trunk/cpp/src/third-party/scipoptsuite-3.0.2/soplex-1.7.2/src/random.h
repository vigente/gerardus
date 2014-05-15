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

/**@file  random.h
 * @brief Random numbers.
 */
#ifndef _RANDOM_H_
#define _RANDOM_H_

#include <assert.h>

namespace soplex
{
#define  RSTEP    1103515245
#define  RDIVIDE  65536
#define  RADD     12345
#define  RMULT    32768

/**@brief   Random numbers.
   @ingroup Elementary

   Class Random provides random Real variables, i.e. a value variable that
   gives another value each time it is accessed. It may be used just like an
   ordinary Real by means of an overloaded cast operator Real()%.
*/
class Random
{
private:

   //--------------------------------------
   /**@name Data */
   //@{
   Real themin;           ///< minimum random number to be returned
   Real themax;           ///< maximum random number to be returned
   unsigned long next;    ///< random seed.
   //@}

   //--------------------------------------
   /**@name Helpers */
   //@{
   /// increases rand seed and returns a pseudo random Real value in [0,1).
   Real next_random ()
   {
      next = next * RSTEP + RADD;
      return last_random();
   }

   /// returns the last used random value in [0,1).
   Real last_random() const
   {
      Real i = int ((next / RDIVIDE) % RMULT);
      return ( i / RMULT );
   }
   //@}

public:

   //--------------------------------------
   /**@name Access */
   //@{
   /// returns lower bound of random numbers.
   Real min() const
   {
      return themin;
   }
   /// returns upper bound of random numbers.
   Real max() const
   {
      return themax;
   }

   /// returns next random number.
   /** When a Random variable is used where a Real value is
       expected, a new random number within the range specified in the
       constructor is retured.
    */
   operator Real()
   {
      return (themin + (themax - themin) * next_random());
   }

   /// returns last random number or seed for next one.
   Real last() const
   {
      return (themin + (themax - themin) * last_random());
   }
   //@}

   //--------------------------------------
   /**@name Modification */
   //@{
   /// resets lower bound for random numbers.
   void setMin(Real p_min)
   {
      themin = p_min;
   }

   /// resets upper bound for random numbers.
   void setMax(Real p_max)
   {
      themax = p_max;
   }

   /// resets seed for next random number.
   void setSeed(Real seed)
   {
      seed = (seed - themin) / (themax - themin);
      next = static_cast<unsigned int>(seed * RMULT * RDIVIDE);
   }
   //@}

   //--------------------------------------
   /**@name Debugging */
   //@{
   /// consistency check.
   bool isConsistent() const
   {
#ifdef ENABLE_CONSISTENCY_CHECKS
      return themin <= themax;
#else
      return true;
#endif
   }
   //@}

   //--------------------------------------
   /**@name Constructors / destructors */
   //@{
   /// default constructor.
   /** Constructs a new (pseudo) Random variable returning values between
       \p p_min and \p p_max and using \p p_seed as seed for the random
       variable's sequence.
   */
   explicit
   Random(Real p_min = 0, Real p_max = 1, Real p_seed = 0.5)
      : themin(p_min), themax(p_max)
   {
      if (p_seed < p_min || p_seed > p_max)
         p_seed = (p_min + p_max) / 2;
      setSeed(p_seed);

      assert(isConsistent());
   }
   /// destructor
   ~Random() 
   {}
   //@}
};

} // namespace soplex
#endif // _RANDOM_H_

//-----------------------------------------------------------------------------
//Emacs Local Variables:
//Emacs mode:c++
//Emacs c-basic-offset:3
//Emacs tab-width:8
//Emacs indent-tabs-mode:nil
//Emacs End:
//-----------------------------------------------------------------------------
