<div id="faq">    
   <p class="q">How can I contribute?</p>

   <p class="a">Please report
    <a href="bugs.html">bugs and problems</a>.<br/>

     In fact, a couple of improvements of SoPlex are
     possible. If you are interested in helping please contact
     us.
   </p>


   <p class="q">
     So, how do I start working with SoPlex?
   </p>

   <p class="a">
     To just read LPs from a file and solve them without modifying the LPs, 
     compile SoPlex and use it on the command line. Calling the binary
     without parameters gives you some usage information. To use SoPlex as a
     callable library, see src/example.cpp or src/soplexmain.cpp (from which 
     the SoPlex binary is compiled) on how to construct an LP, solve it, and 
     extract the solution. For further details, you have to consult the doxygen
     documentation of the code (or even the code itself).<br>
     Another good starting point might also be
     the SoPlex interface of <a href="http://scip.zib.de">SCIP</a>.
   </p>

   <p class="q">
     Can I use SoPlex with Branch-and-cut or MIP solvers?
   </p>

   <p class="a">
     Yes, but the corresponding software needs an interface.

     Currently, the following frameworks have an
     interface to SoPlex.</p>
     <ul>
     <li><a  href="http://scip.zib.de">SCIP</a> - Solving
     Constraint Integer Programs.</li>
     <li><a href="http://www.informatik.uni-koeln.de/abacus/">ABACUS</a>
     has an older interface to 1.2.1 only.</li>
     <li>There is also an interface to <a
     href="https://projects.coin-or.org/Osi">COIN
     OSI</a>.</li>
     </ul>

   <p class="q">
     Why is &lt;iostream&gt; used but &lt;assert.h&gt; and not &lt;cassert&gt;?
   </p>

   <p class="a">
     The reason is twofold. From the theoretical point we were not
     able to exactly find out in TC++PL in what namespace cassert 
     should load its declarations. Surely in std. But since these are
     normally functions with C linkage, this won't work.
     Then some of them, like assert, are macros which have no namespace.
     The practical point is that the compiler vendors seem to be 
     unsure as well. Most of them put everything in both namespaces std 
     and global. So there is no advantage in using &lt;cassert&gt;. Compaq 
     even left them off because it seemed unclear to them.
     So our reasoning was: If it is a C++ header we use the correct form
     without the .h and in std. If it is a C header, according to the
     standard the .h header has to be there and uses the global namespace.
     That seems acceptable to us, especially for C functions.
   </p>
      
   <p class="q">
     Why is malloc/free sometimes used and not new/delete?
   </p>

   <p class="a">
     Because there is no realloc with new/delete. Because malloc
     is faster. And we only use it for builtin types or so called 
     "Data Objects" .
     If you do not like this decision, it is quite easy to change 
     spxalloc.h such as to use new/delete.
   </p>

   <p class="q">
     Can SoPlex solve Integer Programs (IPs)?
   </p>

   <p class="a">
     No. You need an IP-Solver for this. Most IP-Solver use LP-Solvers
     as a subroutine and do some kind of Branch-and-Bound. For
     instance, you can use  SCIP (Solving Constraint Integer Programs) 
     together with SoPlex to solve IPs.  SCIP can be obtained at 
     <a href="http://scip.zib.de/">http://scip.zib.de/</a> and is 
     distributed under the ZIB academic license, like SoPlex.
   </p>

   <p class="q">
     Is there a Windows version?
   </p>

   <p class="a">
     The code is tested to compile under Microsoft Visual Studio 2008 and 2010.
     It should also be possible to create a Windows version by compiling it
     with the GCC under <a href="http://www.cygwin.com/">Cygwin</a> or minGW
   </p>

   <p class="q">
     I want a primal and a dual simplex, where are they?
   </p>

   <p class="a">
    That is nearly easy. You can set SoPlex to use the  
    ENTERing and LEAVEing algorithm and
    COLUMN and ROW representation.
   </p>

   <center>
    <table border="1" rules="all">
    <tr><td>&nbsp;</td><td>ENTER </td><td>LEAVE </td></tr>
    <tr><td>COLUMN</td><td>Primal</td><td>Dual  </td></tr>
    <tr><td>ROW   </td><td>Dual  </td><td>Primal</td></tr>
    </table>
   </center>

   <p class="a">
     COLUMN oriented is the "usual" representation.
     Then Entering is the Primal and Leaving is the Dual algorithm.
     In ROW oriented representation, we have in principle the
     explicit dual and then the algorithms are reversed.

     The only problem is that SoPlex is a composite simplex algorithm.
     That means it switches between entering and leaving algorithm
     as it needs. So all you can choose is which algorithm is used first,
     but then an arbitrary number of switches may occur. (Even so, often
     no switch at all happens.)
   </p>

   <p class="q">
     I got a segment violation or a signal 11.
   </p>

   <p class="a">
     If all of the test instances from Netlib work, but your LP gives this
     problem, mail your LP in as an gzip'ed MPS of LP file and we will
     check. If you have this problem also with the test instances, check your
     stack space:  ulimit  -s will report the current size in kilobytes.
     Try a higher value. If this doesn't help, maybe your compiler is broken.
     Try compiling without optimization. 
   </p>
       
   <p class="q">
     I got messages or exceptions stating &quot;This shall not be.&quot; or similar.    
   </p>

   <p class="a">
     Such messages are an indication of numerical trouble, which may happen
     from time to time especially on LPs with weird numbers. If the problem
     occurs with the binary, try using another scaler (call soplex without
     parameters to see how to specify this). If the problem occurs while using
     the callable library of SoPlex, you can try to add a scaler and a
     preprocessor yourself (see the SoPlex class) unless you solve the same LP
     repeatedly with small changes (as in a Branch-and-cut process, for
     instance). If you send us an LP or MPS file with the offending linear
     program together with a description of your settings, we may have a
     closer look at it to see if we can do something about it.
   </p>

   <p class="q">
     Can SoPlex handle basis files?
   </p>
   <p class="a">
     Yes, the binary version allows to read and write basis files using the
     options <code>-br</code> and <code>-bw</code>, respectively. The library 
     version allows to read or write a basis file from within the code as well. 
     This can be useful for debugging. If you cannot reconstruct an error occurring 
     in the library code with the binary version, try to save a basis and try again. 
     Note that the LP needs to be in MPS format in order to ensure compatability of the 
     basis and the LP (the LP format cannot store ranged rows).
  </p>
       
   <p class="q">
     SoPlex means &quot;Sequential Simplex&quot;. Is there a parallel version
     available? 
   </p>

   <p class="a">
     No. Some research was done in this direction. You can find most of 
     the results in
     <a href="http://www.zib.de/PaperWeb/abstracts/TR-96-09">
     http://www.zib.de/PaperWeb/abstracts/TR-96-09</a> and 
     <a href="http://www.zib.de/PaperWeb/abstracts/SC-95-45">
     http://www.zib.de/PaperWeb/abstracts/SC-95-45</a>.
   </p>

   <p class="q">
     Is there a wrapper class/library to use SoPlex instead of CPLEX ? 
   </p>

   <p class="a">
     No. 
   </p>

   <p class="q">
     How can I make LP generation easier?
   </p>

   <p class="a">
     You can use <a href="http://zimpl.zib.de/">ZIMPL</a>, available at 
     under the ZIB academic license. It takes a (human readable) file
     describing the linear program together with a data file as input and 
     generates LPs or MIPs in LP- or MPS-format.
   </p>

   <p class="q">
     What kind of preprocessing is done?
   </p>

   <p class="a">
     In the default setting SoPlex presolves the given LP in order to
     simplify the problem as far as possible by removing rows,
     columns, and bounds. Moreover, infeasibility or unboundedness may
     be detected. The simplification cannot be undone, but given a
     primal/dual solution for the simplified LP, the simplifier can
     reconstruct the primal/dual solution and basis of the
     unsimplified LP.<br/><br/>
     Handled are:</p>
     <ul>
     <li> empty rows / columns                                  </li>
     <li> unconstrained rows                                    </li>
     <li> row singletons                                        </li>
     <li> forcing rows                                          </li>
     <li> zero objective column singletons                      </li>
     <li> (implied) free column singletons                      </li>
     <li> doubleton equations combined with a column singleton  </li>
     <li> (implicitly) fixed columns                            </li>
     <li> redundant lhs / rhs                                   </li>
     <li> redundant variable bounds                             </li>
     <li> variables that are free in one direction              </li>
     <li> (weakly) dominated columns                            </li>
     <li> duplicate rows / columns                              </li>
     </ul>

   <p class="a">
     The preprocessing can be switched off by using the option  <code>-s0</code>.
   </p>
</div>
