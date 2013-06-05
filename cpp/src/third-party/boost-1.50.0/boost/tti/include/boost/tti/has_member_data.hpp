
//  (C) Copyright Edward Diener 2011
//  Use, modification and distribution are subject to the Boost Software License,
//  Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
//  http://www.boost.org/LICENSE_1_0.txt).

#if !defined(TTI_HAS_MEMBER_DATA_HPP)
#define TTI_HAS_MEMBER_DATA_HPP

#include <boost/config.hpp>
#include <boost/preprocessor/cat.hpp>
#include <boost/type_traits/remove_const.hpp>
#include <boost/tti/gen/has_member_data_gen.hpp>
#include <boost/tti/gen/namespace_gen.hpp>
#include <boost/tti/detail/dmem_data.hpp>

/*

  The succeeding comments in this file are in doxygen format.

*/

/** \file
*/

/// Expands to a metafunction which tests whether a member data with a particular name and type exists.
/**

    trait = the name of the metafunction within the tti namespace.
    
    name  = the name of the inner member to introspect.

    generates a metafunction called "trait" where 'trait' is the macro parameter.
    
              template<class TTI_T,class TTI_Type>
              struct trait
                {
                static const value = unspecified;
                typedef mpl::bool_<true-or-false> type;
                };

              The metafunction types and return:
    
                TTI_T    = the enclosing type in which to look for our 'name'.
                
                TTI_Type = the type of the member data.
                
                returns  = 'value' is true if the 'name' exists, with the TTI_Type data type,
                           within the enclosing TTI_T type, otherwise 'value' is false.
                          
*/
#define BOOST_TTI_TRAIT_HAS_MEMBER_DATA(trait,name) \
  TTI_DETAIL_TRAIT_HAS_MEMBER_DATA(trait,name) \
  template<class TTI_T,class TTI_Type> \
  struct trait : \
    BOOST_PP_CAT(trait,_detail)<typename BOOST_TTI_NAMESPACE::detail::ptmd<TTI_T,TTI_Type>::type,typename boost::remove_const<TTI_T>::type> \
    { \
    }; \
/**/

/// Expands to a metafunction which tests whether a member data with a particular name and type exists.
/**

    name  = the name of the inner member.

    generates a metafunction called "has_member_data_name" where 'name' is the macro parameter.
    
              template<class TTI_T,class TTI_Type>
              struct has_member_data_name
                {
                static const value = unspecified;
                typedef mpl::bool_<true-or-false> type;
                };

              The metafunction types and return:
    
                TTI_T    = the enclosing type in which to look for our 'name'.
                
                TTI_Type = the type of the member data.
                
                returns  = 'value' is true if the 'name' exists, with the TTI_Type data type,
                           within the enclosing TTI_T type, otherwise 'value' is false.
                          
*/
#define BOOST_TTI_HAS_MEMBER_DATA(name) \
  BOOST_TTI_TRAIT_HAS_MEMBER_DATA \
  ( \
  BOOST_TTI_HAS_MEMBER_DATA_GEN(name), \
  name \
  ) \
/**/

#endif // TTI_HAS_MEMBER_DATA_HPP
