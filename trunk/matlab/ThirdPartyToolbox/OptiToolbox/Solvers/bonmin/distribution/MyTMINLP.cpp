// (C) Copyright Carnegie Mellon University 2006
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Authors :
// P. Bonami, Carnegie Mellon University
//
// Date :  03/17/2006
#include "MyTMINLP.hpp"
#include "mex.h"

bool 
MyTMINLP::get_variables_types(Index n, VariableType* var_types)
{
  var_types[0] = BINARY;
  var_types[1] = CONTINUOUS;
  var_types[2] = CONTINUOUS;
  var_types[3] = INTEGER;
  return true;
}


bool 
MyTMINLP::get_variables_linearity(Index n, Ipopt::TNLP::LinearityType* var_types)
{
  var_types[0] = Ipopt::TNLP::LINEAR;
  var_types[1] = Ipopt::TNLP::NON_LINEAR;
  var_types[2] = Ipopt::TNLP::NON_LINEAR;
  var_types[3] = Ipopt::TNLP::LINEAR;
  return true;
}


bool 
MyTMINLP::get_constraints_linearity(Index m, Ipopt::TNLP::LinearityType* const_types)
{
  assert (m==3);
  const_types[0] = Ipopt::TNLP::NON_LINEAR;
  const_types[1] = Ipopt::TNLP::LINEAR;
  const_types[2] = Ipopt::TNLP::LINEAR;
  return true;
}
bool 
MyTMINLP::get_nlp_info(Index& n, Index&m, Index& nnz_jac_g,
                       Index& nnz_h_lag, TNLP::IndexStyleEnum& index_style)
{
  n = 4;//number of variable
  m = 3;//number of constraints
  nnz_jac_g = 7;//number of non zeroes in Jacobian
  nnz_h_lag = 2;//number of non zeroes in Hessian of Lagrangean
  index_style = TNLP::FORTRAN_STYLE;
  return true;
}

bool 
MyTMINLP::get_bounds_info(Index n, Number* x_l, Number* x_u,
                            Index m, Number* g_l, Number* g_u)
{
  assert(n==4);
  assert(m==3);
  x_l[0] = 0.;
  x_u[0] = 1.;
  
  x_l[1] = 0.;
  x_u[1] = DBL_MAX;
  
  x_l[2] =0.;
  x_u[2] = DBL_MAX;
  
  x_l[3] = 0;
  x_u[3] = 5;
  
  g_l[0] = -DBL_MAX;
  g_u[0] = 1./4.;

  g_l[1] = -DBL_MAX;
  g_u[1] = 0;
  
  g_l[2] = -DBL_MAX;
  g_u[2] = 2;
  return true;
}

bool 
MyTMINLP::get_starting_point(Index n, bool init_x, Number* x,
                             bool init_z, Number* z_L, Number* z_U,
                             Index m, bool init_lambda,
                             Number* lambda)
{
  assert(n==4);
  assert(m==3);
  
  assert(init_x);
  assert(!init_lambda);
  x[0] = 0;
  x[1] = 0;
  x[2] = 0;
  x[3] = 0;
  return true;
}

bool 
MyTMINLP::eval_f(Index n, const Number* x, bool new_x, Number& obj_value)
{
  assert(n==4);
  obj_value = - x[0] - x[1] - x[2];
  return true;
}

bool
MyTMINLP::eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f)
{
  assert(n==4);
  grad_f[0] = -1.;
  grad_f[1] = -1.;  
  grad_f[2] = -1.;
  grad_f[3] = 0.;
  return true;
}

bool
MyTMINLP::eval_g(Index n, const Number* x, bool new_x, Index m, Number* g)
{
  assert(n==4);
  assert(m==3);
  
  g[0] = (x[1] - 1./2.)*(x[1] - 1./2.) + (x[2] - 1./2.)*(x[2] - 1./2.);
  g[1] = x[0] - x[1];
  g[2] = x[0] + x[2] + x[3];
  
  return true;
}

bool
MyTMINLP::eval_jac_g(Index n, const Number* x, bool new_x,
                     Index m, Index nnz_jac, Index* iRow, Index *jCol,
                     Number* values)
{
  assert(n==4);
  assert(nnz_jac == 7);
  if(values == NULL) {
    iRow[0] = 2;
    jCol[0] = 1;

    iRow[1] = 3;
    jCol[1] = 1;
    
    iRow[2] = 1;
    jCol[2] = 2;
    
    iRow[3] = 2;
    jCol[3] = 2;
    
    iRow[4] = 1;
    jCol[4] = 3;
        
    iRow[5] = 3;
    jCol[5] = 3;
    
    iRow[6] = 3;
    jCol[6] = 4;
    return true;
  }
  else {
    values[0] = 1.;
    values[1] = 1;

    values[2] = 2*x[1] - 1;
    values[3] = -1.;

    values[4] = 2*x[2] - 1;
    values[5] = 1.;
    
    values[6] = 1.;
    
    return true;
  }
}

bool
MyTMINLP::eval_h(Index n, const Number* x, bool new_x,
                 Number obj_factor, Index m, const Number* lambda,
                 bool new_lambda, Index nele_hess, Index* iRow,
                 Index* jCol, Number* values)
{
  assert (n==4);
  assert (m==3);
  assert(nele_hess==2);
  if(values==NULL)
  {
    iRow[0] = 2;
    jCol[0] = 2;
    
    iRow[1] = 3;
    jCol[1] = 3;
  }
  else {
    values[0] = 2*lambda[0];
    values[1] = 2*lambda[0];
  }
  return true;
}

void
MyTMINLP::finalize_solution(TMINLP::SolverReturn status,
                            Index n, const Number* x, Number obj_value)
{
  //std::cout<<"Problem status: "<<status<<std::endl;
    mexPrintf("Problem status: %d\n",status);
  //std::cout<<"Objective value: "<<obj_value<<std::endl;
    mexPrintf("Problem value: %f\n",obj_value);
  if(printSol_ && x != NULL){
    //std::cout<<"Solution:"<<std::endl;
    mexPrintf("Solution: \n");
    for(int i = 0 ; i < n ; i++)
        mexPrintf("x[%d] = %f\n",i,x[i]);
        
      //std::cout<<"x["<<i<<"] = "<<x[i];
      //if(i < n-1) std::cout<<", ";}
    //std::cout<<std::endl;
  }
}
