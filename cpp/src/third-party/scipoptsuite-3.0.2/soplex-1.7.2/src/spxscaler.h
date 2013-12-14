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

/**@file  spxscaler.h
 * @brief LP scaling base class.
 */
#ifndef _SPXSCALER_H_
#define _SPXSCALER_H_

#include <assert.h>

#include "spxdefines.h"
#include "dataarray.h"
#include "spxlp.h"

namespace soplex
{
/**@brief   LP scaler abstract base class.
   @ingroup Algo

   Instances of classes derived from SPxScaler may be loaded to SoPlex in
   order to scale LPs before solving them. SoPlex will load() itself to
   the SPxScaler and then call #scale(). Generally any SPxLP can be
   loaded to a SPxScaler for #scale()%ing it. The scaling can
   be undone by calling unscale().
*/
class SPxScaler
{
protected:

   //-------------------------------------
   /**@name Data */
   //@{
   const char*        m_name;      ///< Name of the scaler
   DataArray < Real > m_colscale;  ///< column scaling factors
   DataArray < Real > m_rowscale;  ///< row scaling factors
   bool               m_colFirst;  ///< do column scaling first 
   bool               m_doBoth;    ///< do columns and rows
   //@}

   //-------------------------------------
   /**@name Protected helpers */
   //@{
   /// setup scale array for the LP.
   virtual void setup(SPxLP& lp);
   /// computes scaling value for a minimum and maximum pair.
   virtual Real computeScale(Real mini, Real maxi) const;
   /// iterates through vecset and calls computeScale() for each vector.
   /**@return maximum ratio between absolute biggest and smallest element for any vector.
    */
   virtual Real computeScalingVecs( const SVSet* vecset, 
                                    const DataArray<Real>& coScaleval, 
                                    DataArray<Real>& scaleval );
   /// applies m_colscale and m_rowscale to the \p lp.
   virtual void applyScaling(SPxLP& lp);
   //@}

public:

   friend std::ostream& operator<<(std::ostream& s, const SPxScaler& sc);

   //-------------------------------------
   /**@name Construction / destruction */
   //@{
   /// constructor
   explicit SPxScaler(const char* name, bool colFirst = false, bool doBoth = true);
   /// copy constructor
   SPxScaler(const SPxScaler& );
   /// assignment operator
   SPxScaler& operator=(const SPxScaler& );
   /// destructor.
   virtual ~SPxScaler();
   /// clone function for polymorphism
   virtual SPxScaler* clone() const = 0;
   //@}

   //-------------------------------------
   /**@name Access / modification */
   //@{
   /// get name of scaler.
   virtual const char* getName() const;
   /// set scaling order.
   virtual void setOrder(bool colFirst); 
   /// set wether column and row scaling should be performed.
   virtual void setBoth(bool both); 
   //@}

   //-------------------------------------
   /**@name Scaling */
   //@{
   /// scale SPxLP. 
   virtual void scale(SPxLP& lp) = 0;
   /// unscale dense primal solution vector given in \p x. 
   virtual void unscalePrimal(Vector& x) const;
   /// unscale dense slack vector given in \p s.
   virtual void unscaleSlacks(Vector& s) const;
   /// unscale dense dual solution vector given in \p pi. 
   virtual void unscaleDual(Vector& pi) const;
   /// unscale dense reduced cost vector given in \p r.
   virtual void unscaleRedCost(Vector& r) const;

   /// absolute smallest column scaling factor
   virtual Real minAbsColscale() const;
   /// absolute biggest column scaling factor
   virtual Real maxAbsColscale() const;
   /// absolute smallest row scaling factor
   virtual Real minAbsRowscale() const;
   /// absolute biggest row scaling factor
   virtual Real maxAbsRowscale() const;
   /// maximum ratio between absolute biggest and smallest element in any column.
   virtual Real maxColRatio(const SPxLP& lp) const;
   /// maximum ratio between absolute biggest and smallest element in any row.
   virtual Real maxRowRatio(const SPxLP& lp) const;
   //@}

   //-------------------------------------
   /**@name Debugging */
   //@{
   /// consistency check
   virtual bool isConsistent() const;
   //@}
};
} // namespace soplex
#endif // _SPXSCALER_H_

//-----------------------------------------------------------------------------
//Emacs Local Variables:
//Emacs mode:c++
//Emacs c-basic-offset:3
//Emacs tab-width:8
//Emacs indent-tabs-mode:nil
//Emacs End:
//-----------------------------------------------------------------------------

