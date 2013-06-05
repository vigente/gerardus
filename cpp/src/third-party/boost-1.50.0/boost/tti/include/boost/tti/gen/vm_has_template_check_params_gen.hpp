
//  (C) Copyright Edward Diener 2011
//  Use, modification and distribution are subject to the Boost Software License,
//  Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
//  http://www.boost.org/LICENSE_1_0.txt).

#if !defined(TTI_VM_TEMPLATE_PARAMS_GEN_HPP)
#define TTI_VM_TEMPLATE_PARAMS_GEN_HPP

#include <boost/preprocessor/config/config.hpp>

#if BOOST_PP_VARIADICS

#include <boost/preprocessor/cat.hpp>

/*

  The succeeding comments in this file are in doxygen format.

*/

/** \file
*/

/// Generates the macro metafunction name for BOOST_TTI_VM_HAS_TEMPLATE_CHECK_PARAMS.
/**
    name  = the name of the class template.

    returns = the generated macro metafunction name.
*/
#define BOOST_TTI_VM_HAS_TEMPLATE_CHECK_PARAMS_GEN(name) \
  BOOST_PP_CAT(has_template_check_params_,name) \
/**/

#endif // BOOST_PP_VARIADICS
#endif // TTI_VM_TEMPLATE_PARAMS_GEN_HPP
