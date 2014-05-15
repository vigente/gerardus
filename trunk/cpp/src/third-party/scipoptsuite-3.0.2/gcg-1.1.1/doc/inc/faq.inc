<div id="faq">

    <h3><a name="Section0">Sections</a></h3>
    <ul>
     <li><a href="#Section1">General Questions about GCG</a></li>
     <li><a href="#Section2">Decomposition file related questions</a></li>
     <li><a href="#Section3">Miscellaneous</a></li>
    </ul>

    <h4><a name="Section1">General Questions about GCG</a></h4>
    <p class="q">
      <a name="Q1">What is GCG?</a>
    </p>

    <p class="a">
      GCG is a generic solver for Dantzig-Wolfe reformulation
      using <b>branch-cut-and-price</b>. It contains all code to serve as
      a <b>standalone solver</b> for structured mixed integer programs (MIPs)
      and the accompanying structure information.

      This FAQ contains common questions as well as starting information and known bugs.
   </p>

    <p class="q">
      <a name="Q2">What MIPs can GCG solve?</a>
    </p>
    <p class="a">
      GCG can solve all problems for which a structure for Dantzig-Wolfe
      reformulation is known. The structure is either automatically detected by
      GCG (see  <a href="#Q3">next question</a>) or supplied as a decomposition
      file (see <a href="group__FILEREADERS.html">filereaders</a>).
   </p>

    <p class="q">
      <a name="Q3">What are the detected structures?</a>
    </p>
    <p class="a">
      GCG has built in detection code for block diagonal problems and for
      classical set partitioning or setcovering master problems. In order to
      benefit from the automatic detection, your problem must contain only set
      partitioning or set covering constraints in the master and non set
      covering or set partitioning constraints in each of the pricing
      subproblems. <br />
      After removing all set partitioning or covering constraints from the
      problem, the pricing problems need to be disjoint, e.g. they share no
      common variables. In these cases, GCG can detect the structure
      automatically. <br />
      This works for bin packing, coloring, capacitated p-median, generalized
      assignment problems as well as for cutting stock, in which the set
      partitioning is relaxed to set partitioning with integer right hand side.
   </p>

    <h4><a name="Section2">Decomposition file related questions</a></h4>
    <p class="q">
      <a name="Q4">Which filereaders can read decompositions</a>
    </p>
    <p class="a">
      Currently GCG reads three different decomposition structure information:
      <ul>
	<li>dec (<a href="reader__dec_8h.html">documentation</a>)</li>
	<li>blk (<a href="reader__blk_8h.html">documentation</a>)</li>
	<li>ref (<a href="reader__ref_8h.html">documentation</a>)</li>
      </ul>
    </p>

    <p class="q">
      <a name="Q5">Why is presolving important for the decomposition</a>
    </p>
    <p class="a">
      As GCG uses several presolving methods from SCIP, the transformed problem
      (see also SCIP FAQ) may change significantly from the original
      problem. Variables or constraints may be deleted or added which renders
      the current decomposition invalid or unsuitable in some cases. GCG does
      some basic sanity checks, however, it doesn't handle all problems and may
      crash if the decomposition is read at the wrong time (e.g. a
      decomposition found after presolving is read before the problem is
      presolved).<br />
      Use only non presolved decompositions and disable presolving if you are
      in doubt!
    </p>

    <p class="q">
      <a name="Q6">After reading the decomposition file, GCG tells me that the
	behaviour is undefined</a>
    </p>
    <p class="a">
      If GCG can not detect whether your decomposition is for the presolved or
      the original problem, it cannot guard you against errors. See the
      <a href="#Q5">previous question</a> for more information.
    </p>

    <h4><a name="Section3">Miscellaneous</a></h4>
    <p class="q">
      <a name="Q7">Why is CTRL-C unsafe to use?</a>
    </p>
    <p class="a">
       In its current version, SCIP will finish the current node when pressing
    CTRL-C. This is a problem if the master problem is not completely finished
    solving. The current node will be marked as finished and two branching
    children are created. Depending on when you press CTRL-C, this may and
    will lead to different solving processes. If an optimal solution is found,
    they will be the same in all runs, however, it might take a substantially
    different amount of time to do so.
    </p>

    <p class="q">
      <a name="Q8">How does aggregation work in GCG?</a>
    </p>
    <p class="a">
       Currently, aggregation in GCG is very simple. If the same variables
      appear in the same order with the same coefficients in the same
      constraints for all pricing problems to be aggregated, GCG will detect
      that. However, it currently does not do more then that e.g., trying to
      permute the pricing problems in order to find out if they can be
      aggregated.<br />

      Furthermore, aggregation is only supported if the master problem
      consists of set partitioning or (in case of a minimization problem) set
      covering constraints. We are currently implementing a generic branching
      scheme as well as more sophisticated detection methods in order to
      improve the aggregation in the future.
    </p>
</div>

