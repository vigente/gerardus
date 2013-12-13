%% Pre-Distribution-Release Test File
% Variety of tests to check opti OK for distribution
clc
clear all

%Run Install Test
opti_Install_Test;

%Check MEX Interfaces all working
vers = checkSolver('ver');

%Run Benchmarks
optiBench('LP');
optiBench('BILP');
optiBench('MILP');
optiBench('QP');
optiBench('SDP');
optiBench('NLS');
optiBench('NLP');

%Run Demos
Basic_demo;
Overload_demo;
LP_demo;
MILP_demo;
QP_demo;
GNLP_demo;
NLP_demo;
Differentiation_demo;
FileIO_demo;
AMPL_demo;
Benchmark_demo;

%Run MATLAB Wrapper Soak Tests
test_probs_mw;
test_probs_mw;
test_probs_mw;

%Run MPS Test
test_mps;

%Run SOS Test
test_sos;

%QC Rest
test_new_qc;

%Run Row Bounds Test
test_rowcon

%Run NumDiff Test
test_numdiff;

%Run Global Test Set
test_global;

%Run iterfun Test
test_iterfunc;

%Run hessian Test
test_hessian;

%Run AMPL Soak Tests
test_ampl;
test_ampl;
test_ampl;

%Run Solver Tests
test_clp;
test_filterSD;
test_nlopt;
test_ooqp;
test_scip_formal;
test_lipsol;

%Run linear constraint tests
test_ipopt_linearcon;
test_bonmin_linearcon;

%Run sdp tests
test_csdp;
test_dsdp;
test_sedumi;
test_sdpio;

%Run other tests
test_objc;
test_sym_diff;

%Plot tests
test_1dplots;
test_3dplots;

%Multisolve Tests
test_multisolve;
auto2fit_problems;

%Dynamic Optimization Tests
test_dyn_nls;
dnls_examples;

%Doc Test
opti_doc_examples;

%Run General Test Set
test_probs;

clc
fprintf('\nAll Appears OK!\n\n');
