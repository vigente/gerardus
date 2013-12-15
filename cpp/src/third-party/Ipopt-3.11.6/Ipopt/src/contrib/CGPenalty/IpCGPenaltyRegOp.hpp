// Copyright (C) 2007 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// $Id: IpCGPenaltyRegOp.hpp 1861 2010-12-21 21:34:47Z andreasw $
//
// Authors:  Andreas Waechter        IBM       2007-06-01

#ifndef __IPCGPENALTYREGOP_HPP__
#define __IPCGPENALTYREGOP_HPP__

#include "IpSmartPtr.hpp"

namespace Ipopt
{
  class RegisteredOptions;

  void RegisterOptions_CGPenalty(const SmartPtr<RegisteredOptions>& roptions);

} // namespace Ipopt

#endif
