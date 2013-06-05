
//  (C) Copyright Edward Diener 2011
//  Use, modification and distribution are subject to the Boost Software License,
//  Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
//  http://www.boost.org/LICENSE_1_0.txt).

#if !defined(TTI_DETAIL_TEMPLATE_PARAMS_HPP)
#define TTI_DETAIL_TEMPLATE_PARAMS_HPP

#include <boost/config.hpp>
#include <boost/mpl/bool.hpp>
#include <boost/mpl/has_xxx.hpp>
#include <boost/preprocessor/arithmetic/add.hpp>
#include <boost/preprocessor/arithmetic/sub.hpp>
#include <boost/preprocessor/array/elem.hpp>
#include <boost/preprocessor/array/size.hpp>
#include <boost/preprocessor/punctuation/comma_if.hpp>
#include <boost/preprocessor/repetition/repeat.hpp>
#include <boost/preprocessor/repetition/enum.hpp>
#include <boost/preprocessor/seq/enum.hpp>
#include <boost/preprocessor/seq/size.hpp>

#if !defined(BOOST_MPL_CFG_NO_HAS_XXX_TEMPLATE)

#define TTI_DETAIL_TEMPLATE_PARAMETERS(z,n,args) \
BOOST_PP_ARRAY_ELEM(BOOST_PP_ADD(4,n),args) \
/**/

#define TTI_DETAIL_HAS_MEMBER_IMPLEMENTATION(args,introspect_macro) \
   template \
     < \
     typename T, \
     typename fallback_ \
       = boost::mpl::bool_< BOOST_PP_ARRAY_ELEM(3, args) > \
     > \
   class BOOST_PP_ARRAY_ELEM(0, args) \
     { \
     introspect_macro(args) \
     public: \
       static const bool value \
         = BOOST_MPL_HAS_MEMBER_INTROSPECTION_NAME(args)< T >::value; \
       typedef typename BOOST_MPL_HAS_MEMBER_INTROSPECTION_NAME(args) \
         < \
         T \
         >::type type; \
     }; \
/**/

#if !BOOST_WORKAROUND(BOOST_MSVC, <= 1400)

#define TTI_DETAIL_HAS_MEMBER_MULTI_SUBSTITUTE(z,n,args) \
  template \
    < \
    template \
      < \
      BOOST_PP_ENUM_ ## z \
        ( \
        BOOST_PP_SUB \
          ( \
          BOOST_PP_ARRAY_SIZE(args), \
          4 \
          ), \
        TTI_DETAIL_TEMPLATE_PARAMETERS, \
        args \
        ) \
      > \
    class V \
    > \
  struct BOOST_MPL_HAS_MEMBER_INTROSPECTION_SUBSTITUTE_NAME(args, n) \
    { \
    }; \
/**/

#define TTI_DETAIL_HAS_MEMBER_SUBSTITUTE(args) \
  BOOST_PP_REPEAT \
    ( \
    BOOST_PP_ARRAY_ELEM(2, args), \
    TTI_DETAIL_HAS_MEMBER_MULTI_SUBSTITUTE, \
    args \
    ) \
/**/

#define TTI_DETAIL_HAS_MEMBER_INTROSPECT(args) \
  template< typename U > \
  struct BOOST_MPL_HAS_MEMBER_INTROSPECTION_NAME(args) \
    { \
    TTI_DETAIL_HAS_MEMBER_SUBSTITUTE(args) \
    BOOST_MPL_HAS_MEMBER_REJECT(args, BOOST_PP_NIL) \
    BOOST_MPL_HAS_MEMBER_ACCEPT(args, BOOST_PP_NIL) \
    BOOST_STATIC_CONSTANT \
      ( \
      bool, value = BOOST_MPL_HAS_MEMBER_TEST(args) \
      ); \
    typedef boost::mpl::bool_< value > type; \
    }; \
/**/

#define TTI_DETAIL_HAS_MEMBER_WITH_FUNCTION_SFINAE(args) \
  TTI_DETAIL_HAS_MEMBER_IMPLEMENTATION \
    ( \
    args, \
    TTI_DETAIL_HAS_MEMBER_INTROSPECT \
    ) \
/**/

#else // !!BOOST_WORKAROUND(BOOST_MSVC, <= 1400)

#define TTI_DETAIL_HAS_MEMBER_MULTI_SUBSTITUTE_WITH_TEMPLATE_SFINAE(z,n,args) \
  template \
    < \
    template \
      < \
      BOOST_PP_ENUM_ ## z \
        ( \
        BOOST_PP_SUB \
          ( \
          BOOST_PP_ARRAY_SIZE(args), \
          4 \
          ), \
        TTI_DETAIL_TEMPLATE_PARAMETERS, \
        args \
        ) \
      > \
    class U \
    > \
  struct BOOST_MPL_HAS_MEMBER_INTROSPECTION_SUBSTITUTE_NAME_WITH_TEMPLATE_SFINAE \
    ( \
    args, \
    n \
    ) \
    { \
    typedef \
      BOOST_MPL_HAS_MEMBER_INTROSPECTION_SUBSTITUTE_TAG_NAME(args) \
      type; \
    }; \
/**/

#define TTI_DETAIL_HAS_MEMBER_SUBSTITUTE_WITH_TEMPLATE_SFINAE(args) \
  typedef void \
      BOOST_MPL_HAS_MEMBER_INTROSPECTION_SUBSTITUTE_TAG_NAME(args); \
  BOOST_PP_REPEAT \
    ( \
    BOOST_PP_ARRAY_ELEM(2, args), \
    TTI_DETAIL_HAS_MEMBER_MULTI_SUBSTITUTE_WITH_TEMPLATE_SFINAE, \
    args \
    ) \
/**/

#define TTI_DETAIL_HAS_MEMBER_INTROSPECT_WITH_TEMPLATE_SFINAE(args) \
  BOOST_MPL_HAS_MEMBER_REJECT_WITH_TEMPLATE_SFINAE(args,BOOST_PP_NIL) \
  BOOST_MPL_HAS_MEMBER_ACCEPT_WITH_TEMPLATE_SFINAE(args,BOOST_PP_NIL) \
  template< typename U > \
  struct BOOST_MPL_HAS_MEMBER_INTROSPECTION_NAME(args) \
      : BOOST_MPL_HAS_MEMBER_INTROSPECTION_TEST_NAME(args)< U > { \
  }; \
/**/

#define TTI_DETAIL_HAS_MEMBER_WITH_TEMPLATE_SFINAE(args) \
  TTI_DETAIL_HAS_MEMBER_SUBSTITUTE_WITH_TEMPLATE_SFINAE \
    ( \
    args \
    ) \
  TTI_DETAIL_HAS_MEMBER_IMPLEMENTATION \
    ( \
    args, \
    TTI_DETAIL_HAS_MEMBER_INTROSPECT_WITH_TEMPLATE_SFINAE \
    ) \
/**/

#endif // !BOOST_WORKAROUND(BOOST_MSVC, <= 1400)
#else // !!defined(BOOST_MPL_CFG_NO_HAS_XXX_TEMPLATE)

#define TTI_DETAIL_SAME(trait,name) \
  BOOST_MPL_HAS_XXX_TEMPLATE_NAMED_DEF \
    ( \
    trait, \
    name, \
    false \
    ) \
/**/

#endif // !defined(BOOST_MPL_CFG_NO_HAS_XXX_TEMPLATE)

#if !defined(BOOST_MPL_CFG_NO_HAS_XXX_TEMPLATE)
#if !BOOST_WORKAROUND(BOOST_MSVC, <= 1400)

#define TTI_DETAIL_TRAIT_HAS_TEMPLATE_CHECK_PARAMS(trait,name,tpSeq) \
  TTI_DETAIL_HAS_MEMBER_WITH_FUNCTION_SFINAE \
    (  \
      ( BOOST_PP_ADD(BOOST_PP_SEQ_SIZE(tpSeq),4), ( trait, name, 1, false, BOOST_PP_SEQ_ENUM(tpSeq) ) )  \
    )  \
/**/

#else // !!BOOST_WORKAROUND(BOOST_MSVC, <= 1400)

#define TTI_DETAIL_TRAIT_HAS_TEMPLATE_CHECK_PARAMS(trait,name,tpSeq) \
  TTI_DETAIL_HAS_MEMBER_WITH_TEMPLATE_SFINAE \
    ( \
      ( BOOST_PP_ADD(BOOST_PP_SEQ_SIZE(tpSeq),4), ( trait, name, 1, false, BOOST_PP_SEQ_ENUM(tpSeq) ) )  \
    ) \
/**/

#endif // !BOOST_WORKAROUND(BOOST_MSVC, <= 1400)
#else // !!defined(BOOST_MPL_CFG_NO_HAS_XXX_TEMPLATE)

#define TTI_DETAIL_TRAIT_HAS_TEMPLATE_CHECK_PARAMS(trait,name,tpSeq) \
  TTI_DETAIL_SAME(trait,name) \
/**/

#endif // !defined(BOOST_MPL_CFG_NO_HAS_XXX_TEMPLATE)

#endif // TTI_DETAIL_TEMPLATE_PARAMS_HPP
