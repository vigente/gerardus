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

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

#include "spxdefines.h"
#include "spxsolver.h"
#include "spxpricer.h"
#include "spxratiotester.h"
#include "spxstarter.h"
#include "slinsolver.h"
#include "slufactor.h"

namespace soplex
{
bool SPxSolver::writeState(
   const char*    filename,
   const NameSet* rowNames,
   const NameSet* colNames ) const
{
   METHOD( "SPxSolver::writeState()" );

   std::string ofname;
   std::ofstream ofs;

   // write parameter settings
   ofname = std::string(filename) + ".set";
   ofs.open(ofname.c_str());
   if (!ofs)
      return false;

   std::stringstream table, commandline;
   table
      << "Delta            = " << std::setw(8) << delta()
      << std::endl
      << "Epsilon Zero     = " << std::setw(8) << Param::epsilon()
      << std::endl
      << "Epsilon Factor   = " << std::setw(8) << Param::epsilonFactorization()
      << std::endl
      << "Epsilon Update   = " << std::setw(8) << Param::epsilonUpdate()
      << std::endl
      << "Verbosity        = " << std::setw(8) << Param::verbose()
      << std::endl << std::endl
      << "Algorithm        = " << (type() == SPxSolver::ENTER ? "Entering" : "Leaving")
      << std::endl
      << "Representation   = " << (rep() == SPxSolver::ROW ? "Row" : "Column")
      << std::endl
      << "Update           = " << slinSolver()->getName()
      << std::endl
      << "Pricer           = " << pricer()->getName()
#ifdef PARTIAL_PRICING
      << " (partial, size = " << MAX_PRICING_CANDIDATES << ")"
#endif
      << std::endl
      << "Starter          = " << ((starter() == 0) ? "no" : starter()->getName())
      << std::endl
      << "Ratiotest        = " << ratiotester()->getName()
      << std::endl << std::endl;

   commandline
      << "bin/soplex -g0 -s0"
      << " -f" << feastol()
      << " -o" << opttol()
      << (type() == SPxSolver::ENTER ? " -e" : "")
      << (rep()  == SPxSolver::ROW   ? " -r" : "")
      << (!strcmp(slinSolver()->getName(), "SLU-Eta") ? " -i" : "");
   if (!strcmp(pricer()->getName(), "Dantzig"))
      commandline << " -p0";
   else if (!strcmp(pricer()->getName(), "ParMult"))
      commandline << " -p1";
   else if (!strcmp(pricer()->getName(), "Devex"))
      commandline << " -p2";
   else if (!strcmp(pricer()->getName(), "Hybrid"))
      commandline << " -p3";
   else if (!strcmp(pricer()->getName(), "Steep"))
      commandline << " -p4";
   else if (!strcmp(pricer()->getName(), "Weight"))
      commandline << " -p5";
   else if (!strcmp(pricer()->getName(), "SteepEx"))
      commandline << " -p6";
   if (starter() != 0)
   {
      if (!strcmp(starter()->getName(), "Weight"))
         commandline << " -s1";
      else if (!strcmp(starter()->getName(), "Sum"))
         commandline << " -s2";
      else if (!strcmp(starter()->getName(), "Vector"))
         commandline << " -s3";
   }
   if (!strcmp(ratiotester()->getName(), "Default"))
      commandline << " -t0";
   else if (!strcmp(ratiotester()->getName(), "Harris"))
      commandline << " -t1";
   else if (!strcmp(ratiotester()->getName(), "Fast"))
      commandline << " -t2";
   else if (!strcmp(ratiotester()->getName(), "Bound Flipping"))
      commandline << " -t3";
   commandline  << " -br " << filename << ".mps " << filename << ".bas";
   ofs << "SoPlex Parameters:\n\n"  << table.str() << "Command line     > " << commandline.str();
   ofs.close();

   // write LP
   ofname = std::string(filename) + ".mps";
   ofs.open(ofname.c_str());
   if (!ofs)
      return false;

   writeMPS(ofs, rowNames, colNames, NULL);
   ofs.close();

   // write basis
   ofname = std::string(filename) + ".bas";
   return writeBasisFile(ofname.c_str(), rowNames, colNames);
}

} // namespace soplex


//-----------------------------------------------------------------------------
//Emacs Local Variables:
//Emacs mode:c++
//Emacs c-basic-offset:3
//Emacs tab-width:8
//Emacs indent-tabs-mode:nil
//Emacs End:
//-----------------------------------------------------------------------------
