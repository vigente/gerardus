
//  (C) Copyright Edward Diener 2011
//  Use, modification and distribution are subject to the Boost Software License,
//  Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
//  http://www.boost.org/LICENSE_1_0.txt).

#if !defined(TTI_HAS_TEMPLATE_HPP)
#define TTI_HAS_TEMPLATE_HPP

#include <boost/config.hpp>
#include <boost/mpl/has_xxx.hpp>
#include <boost/preprocessor/cat.hpp>
#include <boost/tti/gen/has_template_gen.hpp>

/*

  The succeeding comments in this file are in doxygen format.

*/

/** \file
*/

/// Expands to a metafunction which tests whether an inner class template with a particular name exists.
/**

    trait = the name of the metafunction within the tti namespace.
    
    name  = the name of the inner template.

    generates a metafunction called "trait" where 'trait' is the macro parameter.
    
              template<class TTI_T>
              struct trait
                {
                static const value = unspecified;
                typedef mpl::bool_<true-or-false> type;
                };

              The metafunction types and return:
    
                TTI_T = the enclosing type in which to look for our 'name'.
                
                returns = 'value' is true if the 'name' template exists within the enclosing type,
                          otherwise 'value' is false.
                          
                The template must have all 'class' ( or 'typename' ) parameters types.
    
*/
#define BOOST_TTI_TRAIT_HAS_TEMPLATE(trait,name) \
  BOOST_MPL_HAS_XXX_TEMPLATE_NAMED_DEF(BOOST_PP_CAT(trait,_detail_mpl), name, false) \
  template<class TTI_T> \
  struct trait : \
    BOOST_PP_CAT(trait,_detail_mpl)<TTI_T> \
    { \
    }; \
/**/
  
/// Expands to a metafunction which tests whether an inner class template with a particular name exists.
/**

    name  = the name of the inner template.

    generates a metafunction called "has_template_name" where 'name' is the macro parameter.
    
              template<class TTI_T>
              struct has_template_name
                {
                static const value = unspecified;
                typedef mpl::bool_<true-or-false> type;
                };

              The metafunction types and return:
    
                TTI_T = the enclosing type in which to look for our 'name'.
                
                returns = 'value' is true if the 'name' template exists within the enclosing type,
                          otherwise 'value' is false.
                          
                The template must have all 'class' ( or 'typename' ) parameters types.
    
*/
#define BOOST_TTI_HAS_TEMPLATE(name) \
  BOOST_TTI_TRAIT_HAS_TEMPLATE \
  ( \
  BOOST_TTI_HAS_TEMPLATE_GEN(name), \
  name \
  ) \
/**/

#endif // TTI_HAS_TEMPLATE_HPP
