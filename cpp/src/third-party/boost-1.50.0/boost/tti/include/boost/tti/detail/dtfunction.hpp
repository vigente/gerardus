
//  (C) Copyright Edward Diener 2011
//  Use, modification and distribution are subject to the Boost Software License,
//  Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
//  http://www.boost.org/LICENSE_1_0.txt).

#if !defined(TTI_DETAIL_TFUNCTION_HPP)
#define TTI_DETAIL_TFUNCTION_HPP

#include <boost/config.hpp>
#include <boost/mpl/push_front.hpp>
#include <boost/function_types/function_type.hpp>

namespace boost
  {
  namespace tti
    {
    namespace detail
      {
      template
        <
        class R,
        class FS,
        class TAG
        >
      struct tfunction_seq
        {
        typedef typename boost::mpl::push_front<FS,R>::type ftseq;
        typedef typename boost::function_types::function_type<ftseq,TAG>::type type;
        };
      }
    }
  }
  
#endif // TTI_DETAIL_TFUNCTION_HPP
