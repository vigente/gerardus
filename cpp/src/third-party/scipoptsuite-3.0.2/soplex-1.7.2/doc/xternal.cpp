//-----------------------------------------------------------------------------
/**@mainpage Overview
   @version  1.7.2
   @author   Roland Wunderling
   @author   Tobias Achterberg
   @author   Timo Berthold
   @author   Andreas Bley
   @author   Ambros Gleixner
   @author   Wei Huang
   @author   Benjamin Hiller
   @author   Thorsten Koch
   @author   Matthias Miltenberger
   @author   Sebastian Orlowski
   @author   Marc Pfetsch
   @author   Eva Ramlow
   @author   Andreas Tuchscherer

   @section Main The Sequential object-oriented simplex class library.

   This software has been implemented as a part of Roland Wunderling's 
   Ph.D. thesis "Paralleler und Objektorientierter Simplex-Algorithmus"
   which can be found at http://www.zib.de/PaperWeb/abstracts/TR-96-09
   (in German).

   SoPlex is part of the SCIP Optmization Suite.  A tutorial article for
   getting started with the SCIP Optimization Suite is available as ZIB
   technical report 12-27 <a href="http://scip.zib.de/doc/ZR-12-27.pdf">here</a>.

   SoPlex is implemented in C++. The code should be compliant with the 
   current ANSI standard. RTTI and STL (other then iostream) are not used. 
   Everything is in one namespace \em soplex.

   - \ref RUN      "Where does it run"

   - \ref INST     "Installation"

   - \ref FLAGS    "How to use the flags in the example program"

   - \ref FAQ      "Frequently asked questions"

   - \ref IR       "Iterative Refinement"

   - \ref PROG     "Programming with SoPlex"
*/
//-----------------------------------------------------------------------------
/**@namespace soplex
   @brief     Everything should be within this namespace.

   We have put the whole class library in the namespace soplex.
   If anything here is defined outside, this is a mistake and 
   should be reported. 
*/
//-----------------------------------------------------------------------------
/**@defgroup Elementary Elementary Classes
   @brief    General purpose classes.
   
   Elementary classes are provided for general purpose use in
   projects way beyond the scope of numerical software or linear
   programming.
*/
//-----------------------------------------------------------------------------
/**@defgroup Algebra Linear Algebra Classes
   @brief Basic data types for linear algebra computations.
   
   Linear algebra classes provide basic data types for (sparse)
   linear algebra computations. However, their functionality is
   restricted to simple operations such as addition and scaling.
   For complex tasks, such as solving linear systems of equations,
   algorithmic classes are provided instead.
*/
//-----------------------------------------------------------------------------
/**@defgroup Algo Algorithmic Classes
   @brief Implementation of numerical algorithms.   
   
   Algorithmic classes serve for implementing a variety of
   algorithms for solving numerical (sub-)problems.
*/
//-----------------------------------------------------------------------------
/**@page DataObjects Data Objects 

    \em Data \em objects refer to C++ objects that do not allocate any
    resources, particularly that do not allocate any memory.  This
    makes them behave just like ordinary C structures, in that both,
    the copy constructor and assignment operator are equivalent to a
    memcopy of the memory taken by the object. Examples for data
    objects are all builtin types such as \c int or \c double or
    \e simple classes such as \c complex.
 
    We distinguish \em data \em objects from general C++ objects that
    may include some allocation of resources. (Note that for general
    C++ objects that do allocate resources, this must be respected by
    providing appropriate copy constructor and assignment operators.)
    An example for a general C++ class is class DataArray.
 
    The distinction between data and general C++ objects becomes
    relevant when using such objects in container classes such as
    DataArray or Array.  
*/
//-----------------------------------------------------------------------------
/**@page RUN On which Platforms is SoPlex running

   The current release of SoPlex was tested to compile with the following
   compilers / architectures:
   <TABLE>
   <TR><TD>Architecture</TD><TD>CPU              </TD><TD>OS                </TD><TD>Compiler             </TD></TR>
   <TR><TD>x86         </TD><TD>Pentium 4 3.2 GHz</TD><TD>Suse Linux 10.0   </TD><TD>GCC 4.0.2            </TD></TR>
   <TR><TD>x86         </TD><TD>Pentium 4 3.8 GHz</TD><TD>Suse Linux 10.0   </TD><TD>Intel 9.0.026        </TD></TR>
   <TR><TD>x86         </TD><TD>Pentium 4 3.2 GHz</TD><TD>Windows XP SP2    </TD><TD>MS Visual Studio 2003</TD></TR>
   <TR><TD>x86_64      </TD><TD>AMD Opteron 875  </TD><TD>SLES 9            </TD><TD>GCC 4.0.2            </TD></TR>
   <TR><TD>PowerPC     </TD><TD>G5 2.3           </TD><TD>MacOS-X 10.4.4    </TD><TD>GCC 4.0.0            </TD></TR>
   <TR><TD>Power4      </TD><TD>Power4 1.7 GHz   </TD><TD>AIX 5.1           </TD><TD>VisualAge 6          </TD></TR>
   <TR><TD>Alpha       </TD><TD>21264a 750 MHz   </TD><TD>OSF1 V5.1         </TD><TD>Compaq 6.5-39        </TD></TR>
   <TR><TD>Sparc       </TD><TD>Ultra 60         </TD><TD>Solaris 9         </TD><TD>Sun Studio 11        </TD></TR>
   </TABLE>

   We expect SoPlex to be compilable using any sufficiently recent compiler.
   Remember, your mileage may vary.
 */
//-----------------------------------------------------------------------------
/**@page INST Installation
 
 \section Prerequisites Prerequisites 
 You need the following programs to compile SoPlex:

  - C++ compiler (e.g. http://www.gnu.org/software/gcc)
  - gnu make http://www.gnu.org/software/make
  - gnu awk http://www.gnu.org/software/gawk (if you want the \c check target)
  - doxygen http://www.doxygen.org (if you want to generate the documentation)

 After receiving SoPlex you have to uncompress it with \c gzip 
 (http://www.gnu.org/software/gzip) and 
 unpack it with \c tar (http://www.gnu.org/software/tar)
 into a directory. Then change to that directory.

 \section Tested Linux/x86, Linux/AXP, Darwin, Tru64, Solaris, IRIX, HP-UX
 If you are working under one those OSes you should try the following:
 
 \c gmake \c COMP=xxx \c OPT=yyy

 with \c xxx being one of \c gnu, \c intel, \c sun, \c hp, \c sgi, \c compaq
 or \c ibm
 and \c yyy one of \c dbg (if you want the debug version) or \c opt 
 (if you want the optimized version).

 This should generate a binary in the \c bin subdirectory. 

 \section Others If the previous section was not for you
 First the Makefile tries to find out the OS and the 
 architecture by using the \c uname command. 
 If your OS or architecture is missing, please update the 
 Makefile and submit the change to me.

 Then a submakefile from \c make/make.OS.ARCH.COMP.OPT is included in
 the main Makefile. You should adapt the compiler flags as needed.
 Be especially careful with the \c AR setting since some C++ compilers do
 not like using the standard \c \c ar program together with code that uses
 templates.

 If this all does not work, change to the \c src directory and type

 \c CC \c *.cpp 

 This should do the trick. Adding \c -DNDEBUG gives you a non debugging 
 version. Add flags as needed. 

 \section Testing Testing the Binary
 After you compiled the binary you should download the Netlib LP files at
 http://www.zib.de/Optimization/Software/Soplex/netlib.tar.gz and
 unpack them in the \c check directory. Then you can try

 \c gmake \c COMP=xxx \c OPT=dbg \c quick

 \c gmake \c COMP=xxx \c OPT=opt \c quick

 \c gmake \c COMP=xxx \c OPT=opt \c check

 Use the \c check target together with \c OPT=dbg only if you have really 
 a lot of time. \c quick should run in a few minutes and \c check will 
 need between less than one hour and a day, depending on your machine.

 \c quick should report no fails at all. \c check should report no fails in the
 \c LC and \c EC columns and in the \c LR and \c ER columns only with 
 the instances greenbea, greenbeb, pilot-ja or pilot87. 
 One or two fails is normal, above four is probably a problem with the 
 compiler. Look how big the error is. Try again with less optimization. 
 Our results are available at 
 http://www.zib.de/Optimization/Software/Soplex/results
 to compare with.

 \section Documentation  Generating the documentation
 If you have \c doxygen (and \c dot) installed, you just can say 

 \c gmake \c doc

 After that the documentation should be in doc/html.

 \section Installation Installation
 The binary is in the \c bin directory, the library in \c lib and all
 headers are in \c src. Feel free to install them at a suitable place.

 \section Naming Naming of the OPT Variable

  - dbg (debugging) -DNDEBUG is \b not set. Optimization is mostly off
    and debugging info is generated.

  - std (standard) -DNDEBUG is set. All optimizations may be switched on
    that do \b not \b alter the floating point behaviour or may 
    otherwise by result in wrong code.

  - opt (optimized) -DNDEBUG is set. All optimizations may be switched on,
    as long as the code seems to run correctly. The code should run
    on the relevant architectures. Best is something like -fast that
    uses the right optimizations for the architecture that is used.

  - opt-XXX (optimized for XXX) like opt-p4. This includes optimization
    for a specific processor. 

  - prf (profile) like opt, but generate profile data.

  - XXX-ld (long-double) same as XXX but with long doubles.
 */
//-----------------------------------------------------------------------------
/**@page FAQ Frequently Asked Questions
 * \htmlinclude faq.inc
 * \htmlinclude faqcss.inc
 */           
//-----------------------------------------------------------------------------
/**@page FLAGS How to use the flags in the example program

   Here are some tips on which flags to use with the example program:

   If you have more constraints (rows) than variables (cols) it is 
   a good idea to try the \c -r flag to choose a row-wise representation
   of the basis.

   Setting the different epsilons to a smaller value like 1e-18 or 1e-20
   (using \c -zz or \c -zu) might improve the quality of the
   solution, but may also slow down the program. Setting the epsilons to
   bigger values may speed up the algorithm, but values greater than 1e-12 
   are definitely a bad idea.

   Setting \c -d to smaller values like 1e-7 or 1e-8 will improve the
   quality of the solution, but it will take longer. Values smaller
   then 1e-9 are not recommended. The \c -d value should be
   substantial bigger then the \c -z values.

   If the default settings are too slow, using \c -e eventually together 
   with \c -p1 might improve the running time.  
*/
//-----------------------------------------------------------------------------
/**@page PROG Programming with SoPlex 
   
   The SoPlex class library comprises classes that may be categorized into
   three different types:

   - Elementary classes are provided for general purpose use in
     projects way beyond the scope of numerical software or linear
     programming.
   - Linear algebra classes provide basic data types for (sparse)
     linear algebra computations. However, their functionality is
     restricted to simple operations such as addition and scaling.
     For complex tasks, such as solving linear systems of equations,
     algorithmic classes are provided instead.
   - Algorithmic classes serve for implementing maybe a variety of
     algorithms for solving numerical (sub-)problems.

   The following sections are dedicated to users who want to
   provide own pricers, ratio test, start basis generation codes or
   LP simplifiers to use with SoPlex or who want to derive own
   implementations (e.g. parallel versions) using SoPlex.

   @section Representation Virtualizing the Representation
   The primal Simplex on the columnwise representation is
   structurally equivalent to the dual Simplex on the rowwise
   representation and vice versa (see below). Hence, it is
   desirable to treat both cases in a very similar manner. This
   is supported by the programmer's interface of SoPlex which
   provides access methods for all internal data in two ways: one
   is relative to the "physical" representation of the LP in
   rows and columns, while the other is relative to the chosen
   basis representation. If e.g. a soplex::SPxPricer is
   written using the second type of methods only (which will
   generally be the case), the same code can be used for running
   SoPlex's simplex algorithm for both representations. 
   We will now give two examples for this
   abstraction from the chosen representation.

   Methods \c vector() will return a column or a row vector,
   corresponding to the chosen basis representation. 
   The other "vectors" will be referred to as \em covectors:
     
   <TABLE>
   <TR><TD>&nbsp;  </TD><TD>ROW      </TD><TD>COLUMN   </TD></TR>
   <TR><TD>vector  </TD><TD>rowVector</TD><TD>colVector</TD></TR>
   <TR><TD>coVector</TD><TD>colVector</TD><TD>rowVector</TD></TR>
   </TABLE>
    
   Whether the soplex::SPxBasis::Desc::Status of a variable indicates that the
   corresponding vector is in the basis matrix or not also depends on the
   chosen representation. Hence, methods \c isBasic() are provided to get the
   correct answer for both representations.  
   
   @section Simplex Vectors and Bounds
   The Simplex algorithms keeps three vectors which are associated to each basis.
   Two of them are required for the pricing, while the third one is needed for
   detecting feasibility of the basis. For all three vectors, bounds are
   defined. The Simplex algorithm changes the basis until all three vectors
   satisfy their bounds, which means that the optimal solution has been found.
    
   With each update of the basis, also the three vectors need to be
   updated. This is best supported by the use of \c UpdateVectors.
    
   @subsection Variables
   The Simplex algorithm works with two types of variables, primals and
   duals.  The primal variables are associated with each column of
   an LP, whereas the dual variables are associated with each row.
   However, for each row a slack variable must be added to the set of
   primals (to represent inequalities), and a reduced cost variable must be
   added for each column (to represent upper or lower bounds). Note, that
   mathematically, one dual variable for each bound (upper and lower) should
   be added. However, this variable would always yield the same value and
   can, hence, be implemented as one.
    
   To summarize, we have a primal variable for each LP column and row
   (i.e., its slack) as well as a dual variable for each LP row and column
   (i.e., its bounds). However, not all these values need to be stored and
   computed, since the structure of the Simplex algorithms allow to
   keep them implicitly.
      
   If the SPxBasis's Status of a row or column is one of \c P_ON_LOWER,
   \c P_ON_UPPER, \c P_FIXED or \c P_FREE, the value of the corresponding
   primal variable is the lower, upper or both bound(s) or 0, respectively.
   The corresponding dual variable needs to be computed. Equivalently, for
   a Status of \c D_FREE, \c D_ON_UPPER, \c D_ON_LOWER, \c D_ON_BOTH or
   \c D_UNDEFINED, the corresponding dual variable is 0, whereas the primal 
   one needs to be computed.

   The following vectors are declared for holding the values to be computed:
   \c primRhs, \c primVec (with dimension \c nCols()) for the primal
   variables, and \c dualRhs, \c dualVec (with dimension \c nRows()) for the 
   dual variables. The additional variable \c addvec (with dimension \c coDim())
   depends on the representation.

   @subsection Bounds 
   Primal and dual variables are bounded (including \f$\pm\infty\f$ as
   bounds).  If all primal variables are within their bounds, the
   Simplex basis is said to be primal feasible. Analogously, if all
   dual variables are within their bounds, its is called dual
   feasible.  If a basis is both, primal and dual feasible, the
   optimal solution has been found.

   In the dual Simplex, the basis is maintained dual feasible, while
   primal feasibility is improved via basis updates. However, for
   numerical reasons dual feasibility must be relaxed from time to time.
   Equivalently, primal feasibility will be relaxed to
   retain numerical stability in the primal Simplex algorithm.

   Relaxation of (dual or primal) feasibility is achieved by
   relaxing the bounds of primal or dual variables. However, for each
   type of Simplex only the corresponding bounds need to be
   relaxed. Hence, we define only one vector of upper and lower bound
   for each row and column and initialize it with primal or dual
   bound, depending on the Simplex type (see \c theURbound,
   \c theLRbound, \c theUCbound, \c theLCbound). 
*/
//-----------------------------------------------------------------------------
/**@page IR Iterative Refinement

   Since version 1.7, SoPlex provides the new feature \em iterative \em refinement that
   allows for computing extended-precision solutions beyond the limits of
   standard floating-point arithmetic.  It may be particularly helpful for
   numerically troublesome LPs and applications that require solutions within
   tight feasibility tolerances.

   Iterative refinement is still under development and deactivated by default.
   In the following, we explain in detail how to install and use iterative
   refinement in SoPlex 1.7.  Be aware that the interface may change in future
   releases.

   @section IR-Installation Installation

   For iterative refinement, SoPlex must be built with support for
   GMP, the GNU Multiple Precision library (http://www.gmplib.org/).  Since
   SoPlex uses the C++ interface of GMP, both libgmp and libgmpxx must be
   available on your system.

   When using SoPlex's Makefile build system, it suffices to make with option
   GMP=true, see the INSTALL file for details.  If you use a different build
   system than the provided Makefile, you need to define the preprocessor flag
   SOPLEX_WITH_GMP and link with libgmp and libgmpxx.

   @section IR-Usage Usage

   SoPlex 1.7 comes with a new parameter called "iterative refinement
   threshold".  At the command line, this threshold can be changed via option
   -R; the callable library provides the methods \ref soplex::SoPlex::setIrthreshold() and
   \ref soplex::SPxSolver::setIrthreshold().  By default, this threshold is 1e-12.

   If GMP support is not available, the primal and dual feasibility tolerance
   cannot be set below the iterative refinement threshold.  With GMP support
   enabled, as soon as either the primal or dual feasibility tolerance (command
   line options -f and -o, respectively) is set to a value smaller than this
   threshold, iterative refinement is performed automatically.

   For a detailed explanation of the iterative refinement algorithm see the ZIB
   technical report 12-19 available online at

   http://opus4.kobv.de/opus4-zib/frontdoor/index/index/docId/1545

   Vaguely speaking, at each round of the iterative refinement process, SoPlex
   computes the primal and dual violation of the current solution in exact
   arithmetic and uses this to set up an LP with modified bounds, sides, and
   objective function.  The solution to this LP is then used to correct the
   solution and reduce primal and dual violation.  This is repeated until the
   requested feasibility and optimality tolerance is reached.

   @section IR-Limitations Current limitations

   Currently, iterative refinement is only implemented for LPs with
   an optimal solution; unbounded or infeasible LPs cannot be handled, yet.
   Furthermore, SoPlex supports exact arithmetic only internally.  Input and
   output cannot yet be performed in standard floating-point precision.  This
   leads to the following two limitations.

   First, SoPlex only accepts LPs in double-precision, hence when parsing an LP
   or MPS file, roundoff errors may be introduced and hence the LP in SoPlex may
   be slightly perturbed.  Be aware that the result returned by SoPlex hence
   might apply for a slightly modified problem.  In extreme cases, this can even
   lead to the internally stored LP becoming slightly infeasible.  While a
   standard flaoting-point LP solver would not detect these minimal
   infeasibilities, the iterative refinement procedure will do so if high
   feasibility tolerances are used.  Because of this, we have currently decided
   that in case infeasibility is detected at higher rounds of iterative
   refinement, SoPlex will return the solution of the previous round.

   Second, SoPlex can currently return only a floating-point rounding of the
   high-precision solution computed internally.  Be aware that even if iterative
   refinement terminates and claims that the requested tolerances have been
   reached, this may not hold for the solution returned by the interface, e.g.,
   the options -x and -y on the command line.  Note that the basis information
   returned by SoPlex (command line option -bw) is exact in any case and tools
   like PerPlex or QSopt_ex can be used to recompute the corresponding
   primal-dual solution in exact arithmetic.

   @section IR-Feedback Feedback

   These limitations will be overcome in future release.  We appreciate any
   feedback, comments, and questions on this new feature.  They can be directed
   to soplex@zib.de and will help improve the future development of SoPlex.
 */
//-----------------------------------------------------------------------------


