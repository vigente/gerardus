
//  (C) Copyright Edward Diener 2011
//  Use, modification and distribution are subject to the Boost Software License,
//  Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
//  http://www.boost.org/LICENSE_1_0.txt).

#if !defined(TTI_VM_DETAIL_TEMPLATE_PARAMS_HPP)
#define TTI_VM_DETAIL_TEMPLATE_PARAMS_HPP

#include <boost/preprocessor/config/config.hpp>

#if BOOST_PP_VARIADICS

#include <boost/mpl/has_xxx.hpp>
#include <boost/mpl/identity.hpp>
#include <boost/preprocessor/arithmetic/add.hpp>
#include <boost/preprocessor/variadic/size.hpp>
#include "dtemplate_params.hpp"

#if !defined(BOOST_MPL_CFG_NO_HAS_XXX_TEMPLATE)
#if !BOOST_WORKAROUND(BOOST_MSVC, <= 1400)

#define TTI_VM_DETAIL_TRAIT_HAS_TEMPLATE_CHECK_PARAMS(trait,name,...) \
  TTI_DETAIL_HAS_MEMBER_WITH_FUNCTION_SFINAE \
    (  \
      ( BOOST_PP_ADD(BOOST_PP_VARIADIC_SIZE(__VA_ARGS__),4), ( trait, name, 1, false, __VA_ARGS__ ) )  \
    )  \
/**/

#else // !!BOOST_WORKAROUND(BOOST_MSVC, <= 1400)

#define TTI_VM_DETAIL_TRAIT_HAS_TEMPLATE_CHECK_PARAMS(trait,name,...) \
  TTI_DETAIL_HAS_MEMBER_WITH_TEMPLATE_SFINAE \
    ( \
      ( BOOST_PP_ADD(BOOST_PP_VARIADIC_SIZE(__VA_ARGS__),4), ( trait, name, 1, false, __VA_ARGS__ ) )  \
    ) \
/**/

#endif // !BOOST_WORKAROUND(BOOST_MSVC, <= 1400)
#else // !!defined(BOOST_MPL_CFG_NO_HAS_XXX_TEMPLATE)

#define TTI_VM_DETAIL_TRAIT_HAS_TEMPLATE_CHECK_PARAMS(trait,name,...) \
  TTI_DETAIL_SAME(trait,name) \
/**/

#endif // !defined(BOOST_MPL_CFG_NO_HAS_XXX_TEMPLATE)
#endif // BOOST_PP_VARIADICS
#endif // TTI_VM_DETAIL_TEMPLATE_PARAMS_HPP
