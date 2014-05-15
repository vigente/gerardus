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

#include <assert.h>
#include <math.h>
#include <string.h>
#include <iostream>
#include <iomanip>
#include <fstream>

#include "spxdefines.h"
#include "soplex.h"
#include "spxsolver.h"

#include "timer.h"
#include "spxgithash.h"
#include "spxpricer.h"
#include "spxdantzigpr.h"
#include "spxparmultpr.h"
#include "spxdevexpr.h"
#include "spxhybridpr.h"
#include "spxsteeppr.h"
#include "spxsteepexpr.h"
#include "spxweightpr.h"
#include "spxratiotester.h"
#include "spxharrisrt.h"
#include "spxdefaultrt.h"
#include "spxfastrt.h"
#include "spxboundflippingrt.h"
#include "spxsimplifier.h"
#include "spxmainsm.h"
#include "spxscaler.h"
#include "spxequilisc.h"
#include "spxgeometsc.h"
#include "spxsumst.h"
#include "spxweightst.h"
#include "spxvectorst.h"
#include "slufactor.h"
#include "spxout.h"

using namespace soplex;


//------------------------------------------------------------------------
// for simplicity: store whether we are in check mode:
static bool checkMode = false;
//------------------------------------------------------------------------


//------------------------------------------------------------------------
//    class MySoPlex
//------------------------------------------------------------------------

/** LP solver class for the command line. */
class MySoPlex : public SoPlex
{
public:
   /// default constructor
   MySoPlex( SPxSolver::Type           p_type = SPxSolver::LEAVE,
             SPxSolver::Representation p_rep  = SPxSolver::COLUMN )
      : SoPlex(p_type, p_rep)
   {}
   //------------------------------------------------------------------------
   /// virtual destructor
   virtual ~MySoPlex()
   {}
   //------------------------------------------------------------------------
   void displayQuality() const
   {
      Real maxviol;
      Real sumviol;

      if ( checkMode )
      {
	 MSG_INFO1( spxout << "IEXAMP05 Violations (max/sum)" << std::endl; )

	 m_solver.qualConstraintViolation(maxviol, sumviol);

	 MSG_INFO1( spxout << "IEXAMP06 Constraints      :"
	    << std::setw(16) << maxviol << "  "
	    << std::setw(16) << sumviol << std::endl; )

	 qualConstraintViolation(maxviol, sumviol);

	 MSG_INFO1( spxout << "IEXAMP07       (unscaled) :"
	    << std::setw(16) << maxviol << "  "
	    << std::setw(16) << sumviol << std::endl; )

	 m_solver.qualBoundViolation(maxviol, sumviol);

	 MSG_INFO1( spxout << "IEXAMP08 Bounds           :"
	    << std::setw(16) << maxviol << "  "
	    << std::setw(16) << sumviol << std::endl; )

	 qualBoundViolation(maxviol, sumviol);

	 MSG_INFO1( spxout << "IEXAMP09       (unscaled) :"
	    << std::setw(16) << maxviol << "  "
	    << std::setw(16) << sumviol << std::endl; )

	 if (!m_vanished)
	 {
	    m_solver.qualSlackViolation(maxviol, sumviol);

	    MSG_INFO1( spxout << "IEXAMP10 Slacks           :"
	       << std::setw(16) << maxviol << "  "
	       << std::setw(16) << sumviol << std::endl; )

	    m_solver.qualRedCostViolation(maxviol, sumviol);

	    MSG_INFO1( spxout << "IEXAMP11 Reduced costs    :"
	       << std::setw(16) << maxviol << "  "
	       << std::setw(16) << sumviol << std::endl; )
#if 0
	    MSG_INFO1( spxout << "IEXAMP12 Proven dual bound:"
	       << std::setw(20)
	       << std::setprecision(20)
	       << m_solver.provedDualbound() << std::endl; )
#endif
	  }
      }
      else
      {
	 MSG_INFO1( spxout << "Violations (max/sum)" << std::endl; )

	 m_solver.qualConstraintViolation(maxviol, sumviol);

	 MSG_INFO1( spxout << "Constraints      :"
	    << std::setw(16) << maxviol << "  "
	    << std::setw(16) << sumviol << std::endl; )

	 qualConstraintViolation(maxviol, sumviol);

	 MSG_INFO1( spxout << "      (unscaled) :"
	    << std::setw(16) << maxviol << "  "
	    << std::setw(16) << sumviol << std::endl; )

	 m_solver.qualBoundViolation(maxviol, sumviol);

	 MSG_INFO1( spxout << "Bounds           :"
	    << std::setw(16) << maxviol << "  "
	    << std::setw(16) << sumviol << std::endl; )

	 qualBoundViolation(maxviol, sumviol);

	 MSG_INFO1( spxout << "      (unscaled) :"
	    << std::setw(16) << maxviol << "  "
	    << std::setw(16) << sumviol << std::endl; )

	 if (!m_vanished)
	 {
	    m_solver.qualSlackViolation(maxviol, sumviol);

	    MSG_INFO1( spxout << "Slacks           :"
	       << std::setw(16) << maxviol << "  "
	       << std::setw(16) << sumviol << std::endl; )

	    m_solver.qualRedCostViolation(maxviol, sumviol);

	    MSG_INFO1( spxout << "Reduced costs    :"
	       << std::setw(16) << maxviol << "  "
	       << std::setw(16) << sumviol << std::endl; )
#if 0
	    MSG_INFO1( spxout << "Proven dual bound:"
	       << std::setw(20)
	       << std::setprecision(20)
	       << m_solver.provedDualbound() << std::endl; )
#endif
	  }
      }
   }
   //------------------------------------------------------------------------
   void displayInfeasibility() const
   {
      assert(m_solver.status() == SPxSolver::INFEASIBLE);

#if 0
      if ( checkMode )
      {
	 if( m_solver.isProvenInfeasible() )
	    MSG_INFO1( spxout << "IEXAMP13 Infeasibility is proven." << std::endl; )
         else
	    MSG_INFO1( spxout << "IEXAMP13 Infeasibility could not be proven!" << std::endl; )
      }
      else
      {
         if ( m_solver.isProvenInfeasible() )
	 {
	    MSG_INFO1( spxout << "Infeasibility is proven." << std::endl; )
	 }
	 else
	 {
	    MSG_INFO1( spxout << "Infeasibility could not be proven!" << std::endl; )
	 }
      }
#endif
   }
};


//------------------------------------------------------------------------
//    Helpers
//------------------------------------------------------------------------

static
void print_version_info()
{
   const char* banner1 =
   "************************************************************************\n"
   "*                                                                      *\n"
   "*       SoPlex --- the Sequential object-oriented simPlex.             *\n"
   ;

   const char* banner2 = 
   "*                                                                      *\n"
   "*    Copyright (C) 1996-2013 Konrad-Zuse-Zentrum                       *\n"
   "*                            fuer Informationstechnik Berlin           *\n"
   "*                                                                      *\n"
   "*  SoPlex is distributed under the terms of the ZIB Academic Licence.  *\n"
   "*  You should have received a copy of the ZIB Academic License         *\n"
   "*  along with SoPlex; If not email to soplex@zib.de.                   *\n"
   "*                                                                      *\n"
   "************************************************************************\n"
   ;

   if( !checkMode )
      std::cout << banner1;

#if (SOPLEX_SUBVERSION > 0)
   if( !checkMode )
      std::cout <<    "*                  Version ";
   else
      std::cout << "SoPlex version ";
   std::cout << SOPLEX_VERSION/100 << "."
             << (SOPLEX_VERSION % 100)/10 << "."
             << SOPLEX_VERSION % 10 << "."
             << SOPLEX_SUBVERSION
             << " - Githash "
             << std::setw(13) << std::setiosflags(std::ios::left) << getGitHash();
   if( !checkMode )
      std::cout << "             *\n" << banner2 << std::endl;
   else
      std::cout << "\n";
#else
   if( !checkMode )
      std::cout <<    "*                  Release ";
   else
      std::cout << "SoPlex release ";
   std::cout << SOPLEX_VERSION/100 << "."
             << (SOPLEX_VERSION % 100)/10 << "."
             << SOPLEX_VERSION % 10
             << " - Githash "
             << std::setw(13) << std::setiosflags(std::ios::left) << getGitHash();
   if( !checkMode )
      std::cout << "               *\n" << banner2 << std::endl;
   else
      std::cout << "\n";
#endif

   /// The following code block is tests and shows compilation parameters.
   std::cout << "[NDEBUG:"
#ifdef NDEBUG
             << "YES"
#else
             << "NO"
#endif
             << "]";

   std::cout << "[WITH_WARNINGS:"
#ifdef WITH_WARNINGS
             << "YES"
#else
             << "NO"
#endif
             << "]";

   std::cout << "[ENABLE_ADDITIONAL_CHECKS:"
#ifdef ENABLE_ADDITIONAL_CHECKS
             << "YES"
#else
             << "NO"
#endif
             << "]";

   std::cout << "[ENABLE_CONSISTENCY_CHECKS:"
#ifdef ENABLE_CONSISTENCY_CHECKS
             << "YES"
#else
             << "NO"
#endif
             << "]";

   std::cout << "[SOPLEX_WITH_GMP:"
#ifdef SOPLEX_WITH_GMP
             << "YES"
#else
             << "NO"
#endif
             << "]" << std::endl;

   std::cout << std::endl;
}

#if 0
static
void print_short_version_info()
{
   const char* banner1 =
   "************************************************************************\n"
   "* SoPlex --- the Sequential object-oriented simPlex. ";
   const char* banner2 =
   "* Copyright (C)  1996-2013 Zuse Institute Berlin                       *\n"
   "************************************************************************\n";

   std::cout << banner1;
#if (SOPLEX_SUBVERSION > 0)
   std::cout <<    "Version "
             << SOPLEX_VERSION/100 << "."
             << (SOPLEX_VERSION % 100)/10 << "."
             << SOPLEX_VERSION % 10 << "."
             << SOPLEX_SUBVERSION
             << "   *\n";
#else
   std::cout <<    "Release "
             << SOPLEX_VERSION/100 << "."
             << (SOPLEX_VERSION % 100)/10 << "."
             << SOPLEX_VERSION % 10
             << "     *\n";
#endif
   std::cout << banner2 << std::endl;
}
#endif

//------------------------------------------------------------------------
static
void print_usage_and_exit( const char* const argv[] )
{
   const char* usage =
      "[options] LPfile [Basfile]\n\n"
      "          LPfile can be either in MPS or LPF format\n\n"
      "options:  (*) indicates default\n"
      "          (!) indicates experimental features which may give wrong results\n"
      " -e        select entering algorithm (default is leaving)\n"
      " -r        select row wise representation (default is column)\n"
      " -i        select Eta-update (default is Forest-Tomlin)\n"
      " -x        output solution vector\n"
      " -y        output dual multipliers\n"
      " -q        display solution quality\n"
      " -br       read file with starting basis from Basfile\n"
      " -bw       write file with optimal basis to Basfile\n"
      " -l        set time limit in seconds\n"
      " -L        set iteration limit\n"
      " -f        set primal feasibility tolerance\n"
      " -o        set optimality, i.e., dual feasibility tolerance\n"
      " -d        set primal and dual feasibility tolerance to same value\n"
      " -R        set threshold for tolerances below which iterative refinement is applied\n"
      " -zz       set general zero tolerance\n"
      " -zf       set factorization zero tolerance\n"
      " -zu       set update zero tolerance\n"
      " -v        set verbosity Level: from 0 (ERROR) to 5 (INFO3), default 3 (INFO1)\n"
      " -V        show program version\n"
      " -C        check mode (for check scripts)\n"
      " -h        show this help\n\n"
      "Simplifier:  Scaler:           Starter:    Pricer:        Ratiotester:\n"
      " -s0 none     -g0 none          -c0 none*   -p0 Textbook   -t0 Textbook\n"
      " -s1 Main*    -g1 uni-Equi      -c1 Weight  -p1 ParMult    -t1 Harris\n"
      "              -g2 bi-Equi*      -c2 Sum     -p2 Devex      -t2 Fast*\n"
      "              -g3 bi-Equi+Geo1  -c3 Vector  -p3 Hybrid!    -t3 Bound Flipping\n"
      "              -g4 bi-Equi+Geo8              -p4 Steep*\n"
      "                                            -p5 Weight\n"
      "                                            -p6 SteepExactSetup\n"
      ;

   std::cerr << "usage: " << argv[0] << " " << usage << std::endl;
   exit(0);
}

//------------------------------------------------------------------------
static
void check_parameter(const char param, const char* const argv[])
{
   if (param == '\0')
      print_usage_and_exit( argv );
}

//------------------------------------------------------------------------
static
void print_algorithm_parameters(
   MySoPlex&                       work,
   const SPxSolver::Representation representation,
   const SLUFactor::UpdateType     update
   )
{
   if ( checkMode )
   {
      MSG_INFO1( spxout
	 << "IEXAMP12 Feastol        = "
	 << std::setw(16) << work.feastol() << std::endl
	 << "IEXAMP52 Opttol         = "
	 << std::setw(16) << work.opttol() << std::endl
	 << "IEXAMP53 Irthreshold    = "
	 << std::setw(16) << work.irthreshold() << std::endl
	 << "IEXAMP13 Epsilon Zero   = "
	 << std::setw(16) << Param::epsilon() << std::endl
	 << "IEXAMP37 Epsilon Factor = "
	 << std::setw(16) << Param::epsilonFactorization() << std::endl
	 << "IEXAMP38 Epsilon Update = "
	 << std::setw(16) << Param::epsilonUpdate() << std::endl
	 << "IEXAMP14 "
	 << (work.type() == SPxSolver::ENTER ? "Entering" : "Leaving")
	 << " algorithm" << std::endl
	 << "IEXAMP15 "
	 << (representation == SPxSolver::ROW ? "Row" : "Column")
	 << " representation" << std::endl
	 << "IEXAMP16 "
	 << (update == SLUFactor::ETA ? "Eta" : "Forest-Tomlin")
	 << " update" << std::endl; )
   }
   else
   {
      MSG_INFO1( spxout
	 << "SoPlex parameters: " << std::endl
	 << "Feastol        = "
	 << std::setw(16) << work.feastol() << std::endl
	 << "Opttol         = "
	 << std::setw(16) << work.opttol() << std::endl
	 << "Irthreshold    = "
	 << std::setw(16) << work.irthreshold() << std::endl
	 << "Epsilon Zero   = "
	 << std::setw(16) << Param::epsilon() << std::endl
	 << "Epsilon Factor = "
	 << std::setw(16) << Param::epsilonFactorization() << std::endl
	 << "Epsilon Update = "
	 << std::setw(16) << Param::epsilonUpdate() << std::endl
	 << std::endl
	 << "algorithm      = " << (work.type() == SPxSolver::ENTER ? "Entering" : "Leaving")
	 << std::endl
	 << "representation = " << (representation == SPxSolver::ROW ? "Row" : "Column")
	 << std::endl
	 << "update         = " << (update == SLUFactor::ETA ? "Eta" : "Forest-Tomlin")
	 << std::endl; )
   }
}

//------------------------------------------------------------------------
static
SPxPricer* get_pricer(const int pricing)
{
   SPxPricer* pricer = 0;
   switch(pricing)
   {
   case 6 :
      pricer = new SPxSteepExPR;
      break;
   case 5 :
      pricer = new SPxWeightPR;
      break;
   case 4 :
      pricer = new SPxSteepPR;
      break;
   case 3 :
      pricer = new SPxHybridPR;
      break;
   case 2 :
      pricer = new SPxDevexPR;
      break;
   case 1 :
      pricer = new SPxParMultPR;
      break;
   case 0 :
      /*FALLTHROUGH*/
   default :
      pricer = new SPxDantzigPR;
      break;
   }

   assert(pricer != 0);
   if ( checkMode )
#ifdef PARTIAL_PRICING
      MSG_INFO1( spxout << "IEXAMP17 " << pricer->getName() << " pricing"
                        << " (partial, size = " << MAX_PRICING_CANDIDATES << ")"
                        << std::endl; )
#else
      MSG_INFO1( spxout << "IEXAMP17 " << pricer->getName() << " pricing"
                        << std::endl; )
#endif
   else
#ifdef PARTIAL_PRICING
      MSG_INFO1( spxout << "pricing        = " << pricer->getName()
                        << " (partial, size = " << MAX_PRICING_CANDIDATES << ")"
                        << std::endl; )
#else
      MSG_INFO1( spxout << "pricing        = " << pricer->getName()
                        << std::endl; )
#endif
   return pricer;
}

//------------------------------------------------------------------------
static
SPxRatioTester* get_ratio_tester(const int ratiotest)
{
   SPxRatioTester* ratiotester = 0;
   switch(ratiotest)
   {
   case 3 :
      ratiotester = new SPxBoundFlippingRT;
      break;
   case 2 :
      ratiotester = new SPxFastRT;
      break;
   case 1 :
      ratiotester = new SPxHarrisRT;
      break;
   case 0 :
      /*FALLTHROUGH*/
   default:
      ratiotester = new SPxDefaultRT;
      break;
   }

   assert(ratiotester != 0);
   if ( checkMode )
      MSG_INFO1( spxout << "IEXAMP18 " << ratiotester->getName() << " ratiotest" << std::endl; )
   else
      MSG_INFO1( spxout << "ratiotest      = " << ratiotester->getName() << std::endl; )
   return ratiotester;
}

//------------------------------------------------------------------------
static
void get_scalers(
   SPxScaler*& prescaler,
   SPxScaler*& postscaler,
   const int   scaling
   )
{
   switch(scaling)
   {
   case 4:
      prescaler  = new SPxEquiliSC(true);
      postscaler = new SPxGeometSC(8);
      break;
   case 3:
      prescaler  = new SPxEquiliSC(true);
      postscaler = new SPxGeometSC(1);
      break;
   case 2 :
      prescaler  = new SPxEquiliSC(true);
      postscaler = 0;
      break;
   case 1 :
      prescaler  = new SPxEquiliSC(false);
      postscaler = 0;
      break;
   case 0 :
      /*FALLTHROUGH*/
   default :
      prescaler  = 0;
      postscaler = 0;
      break;
   }

   if ( checkMode )
   {
      MSG_INFO1( spxout << "IEXAMP19 "
	 << ((prescaler != 0) ? prescaler->getName() : "no")
	 << " / "
	 << ((postscaler != 0) ? postscaler->getName() : "no")
	 << " scaling" << std::endl; )
   }
   else
   {
      MSG_INFO1( spxout << "scaling        = "
	 << ((prescaler != 0) ? prescaler->getName() : "no")
	 << " / "
	 << ((postscaler != 0) ? postscaler->getName() : "no")
	 << std::endl; )
   }
}

//------------------------------------------------------------------------
static
SPxSimplifier* get_simplifier(const int simplifying)
{
   SPxSimplifier* simplifier = 0;
   switch(simplifying)
   {
   case 1 :
      simplifier = new SPxMainSM;
      break;
   case 0  :
      /*FALLTHROUGH*/
   default :
      assert(simplifier == 0);
      break;
   }

   if ( checkMode )
      MSG_INFO1( spxout << "IEXAMP20 " << ((simplifier == 0) ? "no" : simplifier->getName()) << " simplifier" << std::endl; )
   else
      MSG_INFO1( spxout << "simplifier     = " << ((simplifier == 0) ? "no" : simplifier->getName()) << std::endl; )
   return simplifier;
}

//------------------------------------------------------------------------
static
SPxStarter* get_starter(const int starting)
{
   SPxStarter* starter = 0;
   switch(starting)
   {
   case 3 :
      starter = new SPxVectorST;
      break;
   case 2 :
      starter = new SPxSumST;
      break;
   case 1 :
      starter = new SPxWeightST;
      break;
   case 0 :
      /*FALLTHROUGH*/
   default :
      break;
   }

   if ( checkMode )
      MSG_INFO1( spxout << "IEXAMP21 " << ((starter == 0) ? "no" : starter->getName()) << " starter" << std::endl; )
   else
      MSG_INFO1( spxout << "starter        = " << ((starter == 0) ? "no" : starter->getName()) << std::endl; )

   return starter;
}

//------------------------------------------------------------------------
#ifdef SEND_ALL_OUTPUT_TO_FILES
static
void redirect_output(
   std::ostream&  myerrstream,
   std::ostream&  myinfostream
   )
{
   myerrstream .setf( std::ios::scientific | std::ios::showpoint );
   myinfostream.setf( std::ios::scientific | std::ios::showpoint );
   spxout.setStream( SPxOut::ERROR,    myerrstream );
   spxout.setStream( SPxOut::WARNING,  myerrstream );
   spxout.setStream( SPxOut::INFO1,    myinfostream );
   spxout.setStream( SPxOut::INFO2,    myinfostream );
   spxout.setStream( SPxOut::INFO3,    myinfostream );
   spxout.setStream( SPxOut::DEBUG,    myinfostream );
}
#endif
//------------------------------------------------------------------------
static
void read_input_file(
   MySoPlex&      work,
   const char*    filename,
   NameSet&       rownames,
   NameSet&       colnames)
{
   if ( checkMode )
      MSG_INFO1( spxout << "IEXAMP22 loading LP file " << filename << std::endl; )
   else
      MSG_INFO1( spxout << "\nLoading LP file " << filename << std::endl; )

   Timer timer;
   timer.start();

   if ( ! work.readFile(filename, &rownames, &colnames, 0) )
   {
      if ( checkMode )
	 MSG_INFO1( spxout << "EEXAMP23 error while reading file \"" << filename << "\"" << std::endl; )
      else
	 MSG_INFO1( spxout << "error while reading file \""  << filename << "\"" << std::endl; )
      exit(1);
   }
   assert(work.isConsistent());

   timer.stop();

   if ( checkMode )
   {
      MSG_INFO1( spxout << "IEXAMP24 LP has "
	 << work.nRows() << " rows "
	 << work.nCols() << " columns "
	 << work.nNzos() << " nonzeros"
	 << std::endl; )

      MSG_INFO1( spxout << "IEXAMP41 LP reading time: " << timer.userTime() << std::endl; )
   }
   else
   {
      MSG_INFO1( spxout << "LP has "
	 << work.nRows() << " rows "
	 << work.nCols() << " columns "
	 << work.nNzos() << " nonzeros"
	 << std::endl; )

      MSG_INFO1(
	 std::streamsize prec = spxout.precision();
	 spxout << "LP reading time: " << std::fixed << std::setprecision(2) << timer.userTime();
	 spxout << std::scientific << std::setprecision(int(prec)) << std::endl; )
   }
}

//------------------------------------------------------------------------
static
void read_basis_file(
   MySoPlex&      work    ,
   const char*    filename,
   const NameSet* rownames,
   const NameSet* colnames)
{
   MSG_INFO1( spxout << "Reading basis from file (disables simplifier)" << std::endl; )
   if (!work.readBasisFile(filename, rownames, colnames))
   {
      if ( checkMode )
         MSG_INFO1( spxout << "EEXAMP25 error while reading file \"" << filename << "\"" << std::endl; )
      else
         MSG_INFO1( spxout << "Error while reading file \"" << filename << "\"" << std::endl; )
      exit(1);
   }
}

//------------------------------------------------------------------------
static
void solve_LP(MySoPlex& work)
{
   Timer timer;
   timer.start();

   if ( checkMode )
      MSG_INFO1( spxout << "IEXAMP26 solving LP" << std::endl; )
   else
      MSG_INFO1( spxout << "\nSolving LP ..." << std::endl; )

   work.solve();
   timer.stop();

   MSG_INFO1( spxout << "\nSoPlex statistics:\n" << work.statistics(); )
}

//------------------------------------------------------------------------
static
void print_solution_and_status(
   MySoPlex&            work,
   const NameSet&       rownames,
   const NameSet&       colnames,
   const int            precision,
   const bool           print_quality,
   const bool           print_solution,
   const bool           print_dual,
   const bool           write_basis,
   const char*          basisname
   )
{
   // get the solution status
   SPxSolver::Status stat = work.status();

   if ( ! checkMode )
      MSG_INFO1( spxout << std::endl; )
   switch (stat)
   {
   case SPxSolver::OPTIMAL:
      if ( checkMode )
	 MSG_INFO1( spxout << "IEXAMP29 solution value is: " << std::setprecision( precision ) << work.objValue() << std::endl; )
      else
	 MSG_INFO1( spxout << "Solution value is: " << std::setprecision( precision ) << work.objValue() << std::endl; )

      if ( print_quality )
         work.displayQuality();

      if ( print_solution )
      {
         DVector objx(work.nCols());

         if( work.getPrimal(objx) != SPxSolver::ERROR )
         {
            MSG_INFO1( spxout << std::endl << "Primal solution (name, id, value):" << std::endl; )
            for( int i = 0; i < work.nCols(); ++i )
            {
               if ( isNotZero( objx[i], 0.001 * work.feastol() ) )
                  MSG_INFO1( spxout << colnames[ work.cId(i) ] << "\t"
                                    << i << "\t"
                                    << std::setw(17)
                                    << std::setprecision( precision )
                                    << objx[i] << std::endl; )
            }
            MSG_INFO1( spxout << "All other variables are zero (within " << std::setprecision(1) << 0.001*work.feastol() << ")." << std::endl; )
         }
      }
      if ( print_dual )
      {
         DVector objy(work.nRows());
         bool allzero = true;

         if( work.getDual(objy) != SPxSolver::ERROR )
         {
            MSG_INFO1( spxout << std::endl << "Dual multipliers (name, id, value):" << std::endl; )
            for( int i = 0; i < work.nRows(); ++i )
            {
               if ( isNotZero( objy[i] , 0.001 * work.opttol() ) )
               {
                  MSG_INFO1( spxout << rownames[ work.rId(i) ] << "\t"
                                    << i << "\t"
                                    << std::setw(17)
                                    << std::setprecision( precision )
                                    << objy[i] << std::endl; )
                  allzero = false;
               }
            }

            MSG_INFO1( spxout << "All " << (allzero ? "" : "other ") << "dual values are zero (within "
                              << std::setprecision(1) << 0.001*work.opttol() << ")." << std::endl; )

            if( !allzero )
            {
               if( work.spxSense() == SPxLP::MINIMIZE )
               {
                  MSG_INFO1( spxout << "Minimizing: a positive/negative value corresponds to left-hand (>=) resp. right-hand (<=) side."
                                    << std::endl; )
               }
               else
               {
                  MSG_INFO1( spxout << "Maximizing: a positive/negative value corresponds to right-hand (<=) resp. left-hand (>=) side."
                                    << std::endl; )
               }
            }
         }
      }
      if ( write_basis )
      {
         MSG_INFO1( spxout << "Writing basis of original problem to file " << basisname << std::endl; )
         if ( ! work.writeBasisFile( basisname, &rownames, &colnames ) )
         {
            if ( checkMode )
               MSG_INFO1( spxout << "EEXAMP30 error while writing file \"" << basisname << "\"" << std::endl; )
            else
               MSG_INFO1( spxout << "Error while writing file \"" << basisname << "\"" << std::endl; )
         }
      }
      break;
   case SPxSolver::UNBOUNDED:
      if ( checkMode )
	 MSG_INFO1( spxout << "IEXAMP31 LP is unbounded" << std::endl; )
      else
	 MSG_INFO1( spxout << "LP is unbounded" << std::endl; )

      if ( print_solution )
      {
         DVector objx(work.nCols());
         if( work.getPrimal(objx) != SPxSolver::ERROR )
         {
            MSG_INFO1( spxout << std::endl << "Primal solution (name, id, value):" << std::endl; )
            for( int i = 0; i < work.nCols(); ++i )
            {
               if ( isNotZero( objx[i], 0.001 * work.feastol() ) )
                  MSG_INFO1( spxout << colnames[ work.cId(i) ] << "\t"
                                    << i << "\t"
                                    << std::setw(17)
                                    << std::setprecision( precision )
                                    << objx[i] << std::endl; )
            }
            MSG_INFO1( spxout << "All other variables are zero (within " << std::setprecision(1) << 0.001*work.feastol() << ")." << std::endl; )
         }

         DVector objcoef(work.nCols());
         DVector ray(work.nCols());
         if( work.getPrimalray(ray) != SPxSolver::ERROR )
         {
            Real rayobjval = 0.0;

            work.getObj(objcoef);

            MSG_INFO1( spxout << std::endl << "Primal ray (name, id, value):" << std::endl; )
            for( int i = 0; i < work.nCols(); ++i )
            {
               if ( isNotZero( ray[i], 0.001 * work.feastol() ) )
               {
                  rayobjval += ray[i] * objcoef[i];

                  MSG_INFO1( spxout << colnames[ work.cId(i) ] << "\t"
                                    << i << "\t"
                                    << std::setw(17)
                                    << std::setprecision( precision )
                                    << ray[i] << std::endl; )
               }
            }
            MSG_INFO1( spxout << "All other variables have zero value (within " << std::setprecision(1) << 0.001*work.feastol() << ")." << std::endl; )
            MSG_INFO1( spxout << "Objective change per unit along primal ray is " << rayobjval << "." << std::endl; )
         }
      }
      break;
   case SPxSolver::INFEASIBLE:
      if ( checkMode )
	 MSG_INFO1( spxout << "IEXAMP32 LP is infeasible" << std::endl; )
      else
	 MSG_INFO1( spxout << "LP is infeasible" << std::endl; )
      if ( print_solution )
      {
         DVector farkasx(work.nRows());

         if( work.getDualfarkas(farkasx) != SPxSolver::ERROR )
         {
            DVector proofvec(work.nCols());
            double lhs;
            double rhs;

            lhs = 0.0;
            rhs = 0.0;
            proofvec.clear();
            for( int i = 0; i < work.nRows(); ++i )
            {
               if ( isNotZero( farkasx[i], 0.001 * work.opttol() ) )
               {
                  MSG_INFO1( spxout << rownames[ work.rId(i) ] << "\t"
                                    << i << "\t"
                                    << std::setw(16)
                                    << std::setprecision( precision )
                                    << farkasx[i] << "\t"; )
                  LPRow row;
                  work.getRow(i, row);
                  if( row.lhs() > -soplex::infinity )
                  {
                     MSG_INFO1( spxout << row.lhs() << " <= "; );
                  }
                  for( int j = 0; j < row.rowVector().size(); ++j )
                  {
                     if( row.rowVector().value(j) > 0 )
                     {
                        MSG_INFO1( spxout << "+"; )
                     }
                     MSG_INFO1( spxout
                        << row.rowVector().value(j) << " "
                        << colnames[ work.cId(row.rowVector().index(j)) ]
                        << " "; );
                  }
                  if( row.rhs() < soplex::infinity )
                  {
                     MSG_INFO1( spxout << "<= " << row.rhs(); );
                  }
                  MSG_INFO1( spxout << std::endl; )
                  if( farkasx[i] > 0.0 )
                  {
                     lhs += farkasx[i] * row.lhs();
                     rhs += farkasx[i] * row.rhs();
                  }
                  else
                  {
                     lhs += farkasx[i] * row.rhs();
                     rhs += farkasx[i] * row.lhs();
                  }
                  SVector vec(row.rowVector());
                  vec *= farkasx[i];
                  proofvec += vec;
               }
            }

            MSG_INFO1( spxout << "All other row multipliers are zero (within " << std::setprecision(1) << 0.001*work.opttol() << ")." << std::endl; )
            MSG_INFO1( spxout << "Farkas infeasibility proof: \t"; )
            MSG_INFO1( spxout << lhs << " <= "; )

            bool nonzerofound = false;
            for( int i = 0; i < work.nCols(); ++i )
            {
               if ( isNotZero( proofvec[i], 0.001 * work.opttol() ) )
               {
                  if( proofvec[i] > 0 )
                  {
                     MSG_INFO1( spxout << "+"; )
                  }
                  MSG_INFO1( spxout << proofvec[i] << " " << colnames[ work.cId(i) ] << " "; )
                  nonzerofound = true;
               }
            }
            if( !nonzerofound )
            {
               MSG_INFO1( spxout << "0 "; );
            }
            MSG_INFO1( spxout << "<= " << rhs << std::endl; );
         }
      }
      if ( print_quality )
         work.displayInfeasibility();
      if ( write_basis )  // write basis even if we are infeasible
         if ( ! work.writeBasisFile( basisname, &rownames, &colnames ) )
         {
	    if ( checkMode )
	       MSG_INFO1( spxout << "EEXAMP30 error while writing file \"" << basisname << "\"" << std::endl; )
	    else
	       MSG_INFO1( spxout << "Error while writing file \"" << basisname << "\"" << std::endl; )
         }
      break;
   case SPxSolver::ABORT_CYCLING:
      if ( checkMode )
	 MSG_INFO1( spxout << "EEXAMP40 aborted due to cycling" << std::endl; )
      else
	 MSG_INFO1( spxout << "Aborted due to cycling" << std::endl; )
      break;
   case SPxSolver::ABORT_TIME:
      if ( checkMode )
	 MSG_INFO1( spxout << "IEXAMP33 aborted due to time limit" << std::endl; )
      else
	 MSG_INFO1( spxout << "Aborted due to time limit" << std::endl; )
      break;
   case SPxSolver::ABORT_ITER:
      if ( checkMode )
	 MSG_INFO1( spxout << "IEXAMP34 aborted due to iteration limit" << std::endl; )
      else
	 MSG_INFO1( spxout << "Aborted due to iteration limit" << std::endl; )
      break;
   case SPxSolver::ABORT_VALUE:
      if ( checkMode )
	 MSG_INFO1( spxout << "IEXAMP35 aborted due to objective value limit" << std::endl; )
      else
	 MSG_INFO1( spxout << "Aborted due to objective value limit" << std::endl; )
      break;
   case SPxSolver::SINGULAR:
      if ( checkMode )
	 MSG_INFO1( spxout << "EEXAMP39 basis is singular" << std::endl; )
      else
	 MSG_INFO1( spxout << "Basis is singular" << std::endl; )
      break;
   default:
      if ( checkMode )
	 MSG_INFO1( spxout << "EEXAMP36 An error occurred during " << "the solution process" << std::endl; )
      else
	 MSG_INFO1( spxout << "An error occurred during " << "the solution process" << std::endl; )
      break;
   }
   MSG_INFO1( spxout << std::endl; )
}

//------------------------------------------------------------------------
static
void clean_up(
   SPxScaler*&       prescaler,
   SPxScaler*&       postscaler,
   SPxSimplifier*&   simplifier,
   SPxStarter*&      starter,
   SPxPricer*&       pricer,
   SPxRatioTester*&  ratiotester,
   char*&            basisname
   )
{
   if ( prescaler != 0 )
   {
      delete prescaler;
      prescaler = 0;
   }
   if ( postscaler != 0 )
   {
      delete postscaler;
      postscaler = 0;
   }
   if ( simplifier != 0 )
   {
      delete simplifier;
      simplifier = 0;
   }
   if ( starter != 0 )
   {
      delete starter;
      starter = 0;
   }

   assert( pricer != 0 );
   delete pricer;
   pricer = 0;

   assert( ratiotester != 0 );
   delete ratiotester;
   ratiotester = 0;

   if ( basisname != 0 )
      delete [] basisname;
   basisname = 0;
}

//------------------------------------------------------------------------
//    main program
//------------------------------------------------------------------------

int main(int argc, char* argv[])
{
   const char*               filename;
   char*                     basisname      = 0;
   SPxSolver::Type           type           = SPxSolver::LEAVE;
   SPxSolver::Representation representation = SPxSolver::COLUMN;
   SLUFactor::UpdateType     update         = SLUFactor::FOREST_TOMLIN;
   SPxSimplifier*            simplifier     = 0;
   SPxStarter*               starter        = 0;
   SPxPricer*                pricer         = 0;
   SPxRatioTester*           ratiotester    = 0;
   SPxScaler*                prescaler      = 0;
   SPxScaler*                postscaler     = 0;

   try {
      NameSet                   rownames;
      NameSet                   colnames;
      int                       starting       = 0;
      int                       pricing        = 4;
      int                       ratiotest      = 2;
      int                       scaling        = 2;
      int                       simplifying    = 1;
      int                       iterlimit      = -1;
      Real                      timelimit      = -1.0;
      Real                      delta          = DEFAULT_BND_VIOL;
      Real                      feastol        = DEFAULT_BND_VIOL;
      Real                      opttol         = DEFAULT_BND_VIOL;
      Real                      irthreshold    = DEFAULT_BND_VIOL * 1e-6;
      Real                      epsilon        = DEFAULT_EPS_ZERO;
      Real                      epsilon_factor = DEFAULT_EPS_FACTOR;
      Real                      epsilon_update = DEFAULT_EPS_UPDATE;
      int                       verbose        = SPxOut::INFO1;
      bool                      print_solution = false;
      bool                      print_dual     = false;
      bool                      print_quality  = false;
      bool                      read_basis     = false;
      bool                      write_basis    = false;
      int                       precision;
      int                       optidx;

      for(optidx = 1; optidx < argc; optidx++)
      {
         if (*argv[optidx] != '-')
            break;

         switch(argv[optidx][1])
         {
         case 'b' :
            check_parameter(argv[optidx][2], argv); // use -b{r,w}, not -b
            if (argv[optidx][2] == 'r')
               read_basis = true;
            if (argv[optidx][2] == 'w')
               write_basis = true;
            break;
         case 'c' :
            check_parameter(argv[optidx][2], argv); // use -c[0-3], not -c
            starting = atoi(&argv[optidx][2]);
            break;
         case 'd' :
            check_parameter(argv[optidx][2], argv); // use -dx, not -d
            delta = atof(&argv[optidx][2]);
            break;
         case 'f' :
            check_parameter(argv[optidx][2], argv); // use -fx, not -f
            feastol = atof(&argv[optidx][2]);
            break;
         case 'o' :
            check_parameter(argv[optidx][2], argv); // use -ox, not -o
            opttol = atof(&argv[optidx][2]);
            break;
         case 'R' :
            check_parameter(argv[optidx][2], argv); // use -Rx, not -R
            irthreshold = atof(&argv[optidx][2]);
            break;
         case 'e':
            type = SPxSolver::ENTER;
            break;
         case 'g' :
            check_parameter(argv[optidx][2], argv); // use -g[0-5], not -g
            scaling = atoi(&argv[optidx][2]);
            break;
         case 'i' :
            update = SLUFactor::ETA;
            break;
         case 'l' :
            if (argv[optidx][2] == '\0' )  // use -lx, not -l
               print_usage_and_exit( argv );
            timelimit = atof(&argv[optidx][2]);
            break;
         case 'L' :
            if (argv[optidx][2] == '\0' )  // use -Lx, not -L
               print_usage_and_exit( argv );
            iterlimit = atoi(&argv[optidx][2]);
            break;
         case 'p' :
            check_parameter(argv[optidx][2], argv); // use -p[0-5], not -p
            pricing = atoi(&argv[optidx][2]);
            break;
         case 'q' :
            print_quality = true;
            break;
         case 'r' :
            representation = SPxSolver::ROW;
            break;
         case 's' :
            check_parameter(argv[optidx][2], argv); // use -s[0-4], not -s
            simplifying = atoi(&argv[optidx][2]);
            break;
         case 't' :
            check_parameter(argv[optidx][2], argv); // use -r[0-2], not -r
            ratiotest = atoi(&argv[optidx][2]);
            break;
         case 'v' :
            check_parameter(argv[optidx][2], argv); // use -v[0-5], not -v
            if (argv[optidx][2] >= '0' && argv[optidx][2] <= '9')
               verbose = argv[optidx][2] - '0';
            break;
         case 'V' :
            print_version_info();
            exit(0);
         case 'x' :
            print_solution = true;
            break;
         case 'y' :
            print_dual = true;
            break;
         case 'z' :
            check_parameter(argv[optidx][2], argv); // must not be empty
            check_parameter(argv[optidx][3], argv); // must not be empty
            switch(argv[optidx][2])
            {
            case 'z' :
               epsilon = atof(&argv[optidx][3]);
               break;
            case 'f' :
               epsilon_factor = atof(&argv[optidx][3]);
               break;
            case 'u' :
               epsilon_update = atof(&argv[optidx][3]);
               break;
            default :
               print_usage_and_exit( argv );
            }
            break;
         case 'C' :
            checkMode = true;
            break;
         case 'h' :
         case '?' :
            print_version_info();
            //lint -fallthrough
         default :
            print_usage_and_exit( argv );
         }
      }

      // print version
      print_version_info();

      // enough arguments?
      if ((argc - optidx) < 1 + (read_basis ? 1 : 0) + (write_basis ? 1 : 0))
         print_usage_and_exit( argv );
      filename  = argv[optidx];

      ++optidx;

      // switch off simplifier when using a starting basis
      if ( read_basis )
         simplifying = 0;

      if ( read_basis || write_basis )
         basisname = strcpy( new char[strlen(argv[optidx]) + 1], argv[optidx] );

      // set some algorithm parameters
      Param::setEpsilon             ( epsilon );
      Param::setEpsilonFactorization( epsilon_factor );
      Param::setEpsilonUpdate       ( epsilon_update );
      Param::setVerbose             ( verbose );

      // Set the output precision.
      precision = int(-log10(std::min(feastol, opttol))) + 1;

      std::cout.setf( std::ios::scientific | std::ios::showpoint );
      std::cerr.setf( std::ios::scientific | std::ios::showpoint );

#ifdef SEND_ALL_OUTPUT_TO_FILES
      // Example of redirecting output to different files.
      // Default is cerr for errors and warnings, cout for everything else.
      std::ofstream  myerrstream ( "errwarn.txt" );
      std::ofstream  myinfostream( "infos.txt" );
      redirect_output(myerrstream, myinfostream);
#endif

      // create an instance of MySoPlex
      MySoPlex work( type, representation );
      work.setUtype             ( update );
      work.setFeastol           ( std::min(feastol, delta) );
      work.setOpttol            ( std::min(opttol, delta) );
      work.setIrthreshold       ( irthreshold );
      work.setTerminationTime   ( timelimit );
      work.setTerminationIter   ( iterlimit );
      print_algorithm_parameters( work, representation, update );
      assert( work.isConsistent() );

      // set pricer, starter, simplifier, and ratio tester
      work.setPricer    ( pricer      = get_pricer      (pricing) );
      work.setStarter   ( starter     = get_starter     (starting) );
      work.setSimplifier( simplifier  = get_simplifier  (simplifying) );
      work.setTester    ( ratiotester = get_ratio_tester(ratiotest) );
      assert(work.isConsistent());

      // set pre- and postscaler
      get_scalers(prescaler, postscaler, scaling);
      work.setPreScaler (prescaler);
      work.setPostScaler(postscaler);
      assert(work.isConsistent());

      // read the LP from an input file (.lp or .mps)
      read_input_file(work, filename, rownames, colnames);

      // read a basis file if specified
      if (read_basis)
         read_basis_file(work, basisname, &rownames, &colnames);

      // solve the LP
      solve_LP(work);

      // print solution, status, infeasibility system,...
      print_solution_and_status(work, rownames, colnames, precision, print_quality,
                                print_solution, print_dual, write_basis, basisname);

      // clean up
      clean_up(prescaler, postscaler, simplifier, starter, pricer, ratiotester, basisname);

      return 0;
   }
   catch(SPxException& x) {
      std::cout << "exception caught : " << x.what() << std::endl;
      delete [] basisname;
      if (simplifier)
         delete simplifier;
      delete starter;
      delete pricer;
      delete ratiotester;
      delete prescaler;
      delete postscaler;
   }
}
