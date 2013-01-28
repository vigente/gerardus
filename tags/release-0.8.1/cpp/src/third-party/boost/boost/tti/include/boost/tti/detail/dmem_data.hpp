
//  (C) Copyright Edward Diener 2011
//  Use, modification and distribution are subject to the Boost Software License,
//  Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
//  http://www.boost.org/LICENSE_1_0.txt).

#if !defined(TTI_DETAIL_MEM_DATA_HPP)
#define TTI_DETAIL_MEM_DATA_HPP

#include <boost/config.hpp>
#include <boost/mpl/bool.hpp>
#include <boost/preprocessor/cat.hpp>
#include <boost/type_traits/detail/yes_no_type.hpp>
#if !defined(BOOST_MSVC)
#include <boost/tti/detail/dmem_fun.hpp>
#endif

#if defined(BOOST_MSVC)

#define TTI_DETAIL_TRAIT_HAS_MEMBER_DATA(trait,name) \
  template<class T,class C> \
  struct BOOST_PP_CAT(trait,_detail) \
    { \
    template<class> \
    struct return_of; \
    \
    template<class R,class C> \
    struct return_of<R C::*> \
      { \
      typedef R type; \
      }; \
    \
    template<bool,typename U> \
    struct menable_if; \
    \
    template<typename U> \
    struct menable_if<true,U> \
      { \
      typedef U type; \
      }; \
    \
    template<typename U,typename V> \
    static ::boost::type_traits::yes_type check2(V U::*); \
    \
    template<typename U,typename V> \
    static ::boost::type_traits::no_type check2(U); \
    \
    template<typename U,typename V> \
    static typename \
      menable_if \
        < \
        sizeof(check2<U,V>(&U::name))==sizeof(::boost::type_traits::yes_type), \
        ::boost::type_traits::yes_type \
        > \
      ::type \
    has_matching_member(int); \
    \
    template<typename U,typename V> \
    static ::boost::type_traits::no_type has_matching_member(...); \
    \
    template<class U,class V> \
    struct ttc_md \
      { \
      typedef boost::mpl::bool_<sizeof(has_matching_member<V,typename return_of<U>::type>(0))==sizeof(::boost::type_traits::yes_type)> type; \
      }; \
    \
    typedef typename ttc_md<T,C>::type type; \
    \
    BOOST_STATIC_CONSTANT(bool,value=type::value); \
    \
    }; \
/**/

#else // !defined(BOOST_MSVC)

#define TTI_DETAIL_TRAIT_HAS_MEMBER_DATA(trait,name) \
  TTI_DETAIL_TRAIT_HAS_MEMBER_FUNCTION(trait,name) \
/**/

#endif // defined(BOOST_MSVC)

namespace boost
  {
  namespace tti
    {
    namespace detail
      {
      template
        <
        class T,
        class R
        >
      struct ptmd
        {
        typedef R T::* type;
        };
      }
    }
  }
  
#endif // TTI_DETAIL_MEM_DATA_HPP
