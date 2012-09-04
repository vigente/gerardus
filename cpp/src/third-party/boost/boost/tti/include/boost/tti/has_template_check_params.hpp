
//  (C) Copyright Edward Diener 2011
//  Use, modification and distribution are subject to the Boost Software License,
//  Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
//  http://www.boost.org/LICENSE_1_0.txt).

#if !defined(TTI_HAS_TEMPLATE_CHECK_PARAMS_HPP)
#define TTI_HAS_TEMPLATE_CHECK_PARAMS_HPP

#include <boost/config.hpp>
#include <boost/preprocessor/cat.hpp>
#include <boost/tti/gen/has_template_check_params_gen.hpp>
#include <boost/tti/detail/dtemplate_params.hpp>

/*

  The succeeding comments in this file are in doxygen format.

*/

/** \file
*/

/// Expands to a metafunction which tests whether an inner class template with a particular name and signature exists.
/**

    trait = the name of the metafunction within the tti namespace.
    
    name  = the name of the inner class template.
    
    tpSeq = a Boost PP sequence which has the class template parameters.
            Each part of the template parameters separated by a comma ( , )
            is put in a separate sequence element.

    generates a metafunction called "trait" where 'trait' is the macro parameter.
    
              template<class TTI_T>
              struct trait
                {
                static const value = unspecified;
                typedef mpl::bool_<true-or-false> type;
                };

              The metafunction types and return:
    
                TTI_T = the enclosing type in which to look for our 'name'.
                
                returns = 'value' is true if the 'name' class template with the signature
                          as defined by the 'tpSeq' exists within the enclosing TTI_T type,
                          otherwise 'value' is false.
    
*/
#define BOOST_TTI_TRAIT_HAS_TEMPLATE_CHECK_PARAMS(trait,name,tpSeq) \
  TTI_DETAIL_TRAIT_HAS_TEMPLATE_CHECK_PARAMS(BOOST_PP_CAT(trait,_detail),name,tpSeq) \
  template<class TTI_T> \
  struct trait : \
    BOOST_PP_CAT(trait,_detail)<TTI_T> \
    { \
    }; \
/**/

/// Expands to a metafunction which tests whether an inner class template with a particular name and signature exists.
/**

    name  = the name of the inner class template.
    
    tpSeq = a Boost PP sequence which has the class template parameters.
            Each part of the template parameters separated by a comma ( , )
            is put in a separate sequence element.

    generates a metafunction called "has_template_check_params_name" where 'name' is the macro parameter.
    
              template<class TTI_T>
              struct has_template_check_params_name
                {
                static const value = unspecified;
                typedef mpl::bool_<true-or-false> type;
                };

              The metafunction types and return:
    
                TTI_T = the enclosing type in which to look for our 'name'.
                
                returns = 'value' is true if the 'name' class template with the signature
                          as defined by the 'tpSeq' exists within the enclosing TTI_T type,
                          otherwise 'value' is false.
    
*/
#define BOOST_TTI_HAS_TEMPLATE_CHECK_PARAMS(name,tpSeq) \
  BOOST_TTI_TRAIT_HAS_TEMPLATE_CHECK_PARAMS \
  ( \
  BOOST_TTI_HAS_TEMPLATE_CHECK_PARAMS_GEN(name), \
  name, \
  tpSeq \
  ) \
/**/

#endif // TTI_HAS_TEMPLATE_CHECK_PARAMS_HPP
