// Copyright (C) 2004, 2006 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// $Id: IpObserver.cpp 1861 2010-12-21 21:34:47Z andreasw $
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

#include "IpObserver.hpp"

namespace Ipopt
{
#ifdef IP_DEBUG_OBSERVER
  const Index Observer::dbg_verbosity = 0;
  const Index Subject::dbg_verbosity = 0;
#endif
} // namespace Ipopt
