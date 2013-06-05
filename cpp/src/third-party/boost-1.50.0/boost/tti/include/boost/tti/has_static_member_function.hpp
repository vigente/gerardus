
//  (C) Copyright Edward Diener 2011
//  Use, modification and distribution are subject to the Boost Software License,
//  Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
//  http://www.boost.org/LICENSE_1_0.txt).

#if !defined(TTI_HAS_STATIC_MEMBER_FUNCTION_HPP)
#define TTI_HAS_STATIC_MEMBER_FUNCTION_HPP

#include <boost/config.hpp>
#include <boost/function_types/property_tags.hpp>
#include <boost/mpl/vector.hpp>
#include <boost/preprocessor/cat.hpp>
#include <boost/tti/gen/has_static_member_function_gen.hpp>
#include <boost/tti/gen/namespace_gen.hpp>
#include <boost/tti/detail/dstatic_mem_fun.hpp>
#include <boost/tti/detail/dtfunction.hpp>

/*

  The succeeding comments in this file are in doxygen format.

*/

/** \file
*/

/// Expands to a metafunction which tests whether a static member function with a particular name and signature exists.
/**

    trait = the name of the metafunction within the tti namespace.
    
    name  = the name of the inner member.

    generates a metafunction called "trait" where 'trait' is the macro parameter.
    
              template<class TTI_T,class TTI_R,class TTI_FS,class TTI_TAG>
              struct trait
                {
                static const value = unspecified;
                typedef mpl::bool_<true-or-false> type;
                };

              The metafunction types and return:
    
                TTI_T   = the enclosing type in which to look for our 'name'.
                
                TTI_R   = the return type of the static member function.
                
                TTI_FS  = an optional parameter which are the parameters of the static member function as a boost::mpl forward sequence.
                
                TTI_TAG = an optional parameter which is a boost::function_types tag to apply to the static member function.
                
                returns = 'value' is true if the 'name' exists, 
                          with the appropriate static member function type,
                          as defined by TTI_R, TTI_FS, and TTI_TAG,
                          within the enclosing TTI_T type, 
                          otherwise 'value' is false.
                          
*/
#define BOOST_TTI_TRAIT_HAS_STATIC_MEMBER_FUNCTION(trait,name) \
  TTI_DETAIL_TRAIT_HAS_STATIC_MEMBER_FUNCTION(trait,name) \
  template<class TTI_T,class TTI_R,class TTI_FS = boost::mpl::vector<>,class TTI_TAG = boost::function_types::null_tag> \
  struct trait : \
    BOOST_PP_CAT(trait,_detail)<TTI_T,typename BOOST_TTI_NAMESPACE::detail::tfunction_seq<TTI_R,TTI_FS,TTI_TAG>::type> \
    { \
    }; \
/**/

/// Expands to a metafunction which tests whether a static member function with a particular name and signature exists.
/**

    name  = the name of the inner member.

    generates a metafunction called "has_static_member_function_name" where 'name' is the macro parameter.
    
              template<class TTI_T,class TTI_R,class TTI_FS,class TTI_TAG>
              struct trait
                {
                static const value = unspecified;
                typedef mpl::bool_<true-or-false> type;
                };

              The metafunction types and return:
    
                TTI_T   = the enclosing type in which to look for our 'name'.
                
                TTI_R   = the return type of the static member function.
                
                TTI_FS  = an optional parameter which are the parameters of the static member function as a boost::mpl forward sequence.
                
                TTI_TAG = an optional parameter which is a boost::function_types tag to apply to the static member function.
                
                returns = 'value' is true if the 'name' exists, 
                          with the appropriate static member function type,
                          as defined by TTI_R, TTI_FS, and TTI_TAG,
                          within the enclosing TTI_T type, 
                          otherwise 'value' is false.
                          
*/
#define BOOST_TTI_HAS_STATIC_MEMBER_FUNCTION(name) \
  BOOST_TTI_TRAIT_HAS_STATIC_MEMBER_FUNCTION \
  ( \
  BOOST_TTI_HAS_STATIC_MEMBER_FUNCTION_GEN(name), \
  name \
  ) \
/**/

#endif // TTI_HAS_STATIC_MEMBER_FUNCTION_HPP
