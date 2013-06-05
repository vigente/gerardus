
//  (C) Copyright Edward Diener 2011
//  Use, modification and distribution are subject to the Boost Software License,
//  Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
//  http://www.boost.org/LICENSE_1_0.txt).

#if !defined(TTI_DETAIL_TCLASS_HPP)
#define TTI_DETAIL_TCLASS_HPP

#include <boost/mpl/eval_if.hpp>
#include <boost/mpl/identity.hpp>
#include <boost/type_traits/is_class.hpp>

namespace boost
  {
  namespace tti
    {
    namespace detail
      {
      template <class T>
      struct tclass :
        boost::mpl::eval_if
          <
          boost::is_class<T>,
          T,
          boost::mpl::identity<T>
          >
        {
        };
      }
    }
  }
  
#endif // TTI_DETAIL_TCLASS_HPP
