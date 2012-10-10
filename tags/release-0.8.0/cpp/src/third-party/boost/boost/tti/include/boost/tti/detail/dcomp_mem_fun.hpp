
//  (C) Copyright Edward Diener 2011
//  Use, modification and distribution are subject to the Boost Software License,
//  Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
//  http://www.boost.org/LICENSE_1_0.txt).

#if !defined(TTI_DETAIL_COMP_MEM_FUN_HPP)
#define TTI_DETAIL_COMP_MEM_FUN_HPP

#include <boost/config.hpp>
#include <boost/mpl/at.hpp>
#include <boost/mpl/bool.hpp>
#include <boost/mpl/identity.hpp>
#include <boost/mpl/int.hpp>
#include <boost/mpl/quote.hpp>
#include <boost/function_types/parameter_types.hpp>
#include <boost/preprocessor/cat.hpp>
#include <boost/type_traits/remove_const.hpp>
#include <boost/type_traits/detail/yes_no_type.hpp>

#if defined(BOOST_NO_NULLPTR)

#define TTI_DETAIL_TRAIT_HAS_COMP_MEMBER_FUNCTION(trait,name) \
  template<class T> \
  struct BOOST_PP_CAT(trait,_detail) \
    { \
    template<class F> \
    struct class_type \
      { \
      typedef typename \
      boost::remove_const \
        < \
        typename \
        boost::mpl::at \
          < \
          typename \
          boost::function_types::parameter_types \
            < \
            F, \
            boost::mpl::quote1 \
              < \
              boost::mpl::identity \
              > \
            > \
            ::type, \
            boost::mpl::int_<0> \
          >::type \
        >::type \
      type; \
      }; \
    \
    template<T> \
    struct helper; \
    \
    template<class U> \
    static ::boost::type_traits::yes_type check(helper<&U::name> *); \
    \
    template<class U> \
    static ::boost::type_traits::no_type check(...); \
    \
    BOOST_STATIC_CONSTANT(bool,value=sizeof(check<typename class_type<T>::type>(0))==sizeof(::boost::type_traits::yes_type)); \
    \
    typedef boost::mpl::bool_<value> type; \
    }; \
/**/

#else // !defined(BOOST_NO_NULLPTR)

#define TTI_DETAIL_TRAIT_HAS_COMP_MEMBER_FUNCTION(trait,name) \
  template<class T> \
  struct BOOST_PP_CAT(trait,_detail) \
    { \
    template<class F> \
    struct class_type \
      { \
      typedef typename \
      boost::remove_const \
        < \
        typename \
        boost::mpl::at \
          < \
          typename \
          boost::function_types::parameter_types \
            < \
            F, \
            boost::mpl::quote1 \
              < \
              boost::mpl::identity \
              > \
            > \
            ::type, \
            boost::mpl::int_<0> \
          >::type \
        >::type \
      type; \
      }; \
    \
    template<T> \
    struct helper; \
    \
    template<class U> \
    static ::boost::type_traits::yes_type check(helper<&U::name> *); \
    \
    template<class U> \
    static ::boost::type_traits::no_type check(...); \
    \
    BOOST_STATIC_CONSTANT(bool,value=sizeof(check<typename class_type<T>::type>(nullptr))==sizeof(::boost::type_traits::yes_type)); \
    \
    typedef boost::mpl::bool_<value> type; \
    }; \
/**/

#endif // defined(BOOST_NO_NULLPTR)

#endif // TTI_DETAIL_COMP_MEM_FUN_HPP
