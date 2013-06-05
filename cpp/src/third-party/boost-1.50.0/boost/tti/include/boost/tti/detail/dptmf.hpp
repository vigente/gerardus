
//  (C) Copyright Edward Diener 2011
//  Use, modification and distribution are subject to the Boost Software License,
//  Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
//  http://www.boost.org/LICENSE_1_0.txt).

#if !defined(TTI_DETAIL_PTMF_HPP)
#define TTI_DETAIL_PTMF_HPP

#include <boost/config.hpp>
#include <boost/mpl/push_front.hpp>
#include <boost/function_types/member_function_pointer.hpp>

namespace boost
  {
  namespace tti
    {
    namespace detail
      {
      template
        <
        class T,
        class R,
        class FS,
        class TAG
        >
      struct ptmf_seq
        {
        typedef typename boost::mpl::push_front<FS,T>::type tfs1;
        typedef typename boost::mpl::push_front<tfs1,R>::type tfs2;
        typedef typename boost::function_types::member_function_pointer<tfs2,TAG>::type type;
        };
      }
    }
  }
  
#endif // TTI_DETAIL_PTMF_HPP
