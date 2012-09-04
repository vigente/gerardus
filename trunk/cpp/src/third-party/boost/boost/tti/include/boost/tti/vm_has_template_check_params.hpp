
//  (C) Copyright Edward Diener 2011
//  Use, modification and distribution are subject to the Boost Software License,
//  Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
//  http://www.boost.org/LICENSE_1_0.txt).

#if !defined(TTI_VM_HAS_TEMPLATE_CHECK_PARAMS_HPP)
#define TTI_VM_HAS_TEMPLATE_CHECK_PARAMS_HPP

#include <boost/preprocessor/config/config.hpp>

#if BOOST_PP_VARIADICS

#include <boost/preprocessor/cat.hpp>
#include <boost/tti/gen/vm_has_template_check_params_gen.hpp>
#include <boost/tti/detail/dvm_template_params.hpp>

/*

  The succeeding comments in this file are in doxygen format.

*/

/** \file
*/

/// Expands to a metafunction which tests whether an inner class template with a particular name and signature exists.
/**

    trait = the name of the metafunction within the tti namespace.
    
    name  = the name of the inner class template.
    
    ...   = variadic macro data which has the class template parameters.

    generates a metafunction called "trait" where 'trait' is the macro parameter.
    
              template<class TTI_T>
              struct trait
                {
                static const value = unspecified;
                typedef mpl::bool_<true-or-false> type;
                };

              The metafunction types and return:
    
                TTI_T = the enclosing type in which to look for our 'name'.
                
                returns = 'value' is true if the 'name' class template, with the signature
                          as defined by the '...' variadic macro data, exists within the enclosing TTI_T type,
                          otherwise 'value' is false.
    
*/
#define BOOST_TTI_VM_TRAIT_HAS_TEMPLATE_CHECK_PARAMS(trait,name,...) \
  TTI_VM_DETAIL_TRAIT_HAS_TEMPLATE_CHECK_PARAMS(BOOST_PP_CAT(trait,_detail),name,__VA_ARGS__) \
  template<class TTI_T> \
  struct trait \
    { \
    typedef typename BOOST_PP_CAT(trait,_detail)<TTI_T>::type type; \
    \
    BOOST_STATIC_CONSTANT(bool,value=type::value); \
    }; \
/**/

/// Expands to a metafunction which tests whether an inner class template with a particular name and signature exists.
/**

    name  = the name of the inner class template.
    
    ...   = variadic macro data which has the class template parameters.

    generates a metafunction called "has_template_check_params_name" where 'name' is the macro parameter.
    
              template<class TTI_T>
              struct has_template_check_params_name
                {
                static const value = unspecified;
                typedef mpl::bool_<true-or-false> type;
                };

              The metafunction types and return:
    
                TTI_T = the enclosing type in which to look for our 'name'.
                
                returns = 'value' is true if the 'name' class template, with the signature
                          as defined by the '...' variadic macro data, exists within the enclosing TTI_T type,
                          otherwise 'value' is false.
    
*/
#define BOOST_TTI_VM_HAS_TEMPLATE_CHECK_PARAMS(name,...) \
  BOOST_TTI_VM_TRAIT_HAS_TEMPLATE_CHECK_PARAMS \
  ( \
  BOOST_TTI_VM_HAS_TEMPLATE_CHECK_PARAMS_GEN(name), \
  name, \
  __VA_ARGS__ \
  ) \
/**/

#endif // BOOST_PP_VARIADICS
#endif // TTI_VM_HAS_TEMPLATE_CHECK_PARAMS_HPP
