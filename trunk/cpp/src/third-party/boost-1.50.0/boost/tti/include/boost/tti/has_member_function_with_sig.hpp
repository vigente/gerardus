
//  (C) Copyright Edward Diener 2011
//  Use, modification and distribution are subject to the Boost Software License,
//  Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
//  http://www.boost.org/LICENSE_1_0.txt).

#if !defined(TTI_HAS_MEMBER_FUNCTION_WITH_SIG_HPP)
#define TTI_HAS_MEMBER_FUNCTION_WITH_SIG_HPP

#include <boost/config.hpp>
#include <boost/mpl/apply.hpp>
#include <boost/mpl/identity.hpp>
#include <boost/preprocessor/cat.hpp>
#include <boost/tti/gen/has_member_function_with_sig_gen.hpp>
#include <boost/tti/detail/dcomp_mem_fun.hpp>

/*

  The succeeding comments in this file are in doxygen format.

*/

/** \file
*/

/// Expands to a metafunction which tests whether a member function with a particular name and composite type exists.
/**

    trait = the name of the metafunction within the tti namespace.
    
    name  = the name of the inner member.

    generates a metafunction called "trait" where 'trait' is the macro parameter.
    
              template<class TTI_T>
              struct trait
                {
                static const value = unspecified;
                typedef mpl::bool_<true-or-false> type;
                };

              The metafunction types and return:
    
                TTI_T = the member function type, 
                        in the form of a member function pointer - 'return_type (Class::*)(parameter_types...)', 
                        in which to look for our 'name'.
                    
                returns = 'value' is true if the 'name' exists, 
                          with the appropriate member function type as defined by TTI_T,
                          otherwise 'value' is false.
                          
*/
#define BOOST_TTI_TRAIT_HAS_MEMBER_FUNCTION_WITH_SIG(trait,name) \
  TTI_DETAIL_TRAIT_HAS_COMP_MEMBER_FUNCTION(trait,name) \
  template<class TTI_T> \
  struct trait : \
    BOOST_PP_CAT(trait,_detail)<TTI_T> \
    { \
    }; \
/**/

/// Expands to a metafunction which tests whether a member function with a particular name and composite type exists.
/**

    name  = the name of the inner member.

    generates a metafunction called "has_member_function_with_sig_name" where 'name' is the macro parameter.
    
              template<class TTI_T>
              struct has_member_function_with_sig_name
                {
                static const value = unspecified;
                typedef mpl::bool_<true-or-false> type;
                };

              The metafunction types and return:
    
                TTI_T = the member function type, 
                        in the form of a member function pointer - 'return_type (Class::*)(parameter_types...)', 
                        in which to look for our 'name'.
                    
                returns = 'value' is true if the 'name' exists, 
                          with the appropriate member function type as defined by TTI_T,
                          otherwise 'value' is false.
                          
*/
#define BOOST_TTI_HAS_MEMBER_FUNCTION_WITH_SIG(name) \
  BOOST_TTI_TRAIT_HAS_MEMBER_FUNCTION_WITH_SIG \
  ( \
  BOOST_TTI_HAS_MEMBER_FUNCTION_WITH_SIG_GEN(name), \
  name \
  ) \
/**/

#endif // TTI_HAS_MEMBER_FUNCTION_WITH_SIG_HPP
