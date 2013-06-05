//  (C) Copyright Gennadiy Rozental 2011.
//  Distributed under the Boost Software License, Version 1.0.
//  (See accompanying file LICENSE_1_0.txt or copy at
//  http://www.boost.org/LICENSE_1_0.txt)

//  See http://www.boost.org/libs/test for the library home page.
//
//  File        : $RCSfile$
//
//  Version     : $Revision: 74248 $
//
//  Description : implementation of test tools
// ***************************************************************************

#ifndef BOOST_TEST_TEST_TOOLS_IMPL_HPP_012705GER
#define BOOST_TEST_TEST_TOOLS_IMPL_HPP_012705GER

// Boost.Test
#include <boost/test/tools/predicate_result.hpp>
#ifndef BOOST_TEST_PROD
#include <boost/test/unit_test_log.hpp>
#endif

#include <boost/test/tools/floating_point_comparison.hpp>
#include <boost/test/tools/assertion.hpp>

#include <boost/test/detail/config.hpp>
#include <boost/test/detail/global_typedef.hpp>
#include <boost/test/detail/workaround.hpp>

#include <boost/test/utils/wrap_stringstream.hpp>
#include <boost/test/utils/basic_cstring/io.hpp>
#include <boost/test/utils/lazy_ostream.hpp>

// Boost
#include <boost/limits.hpp>

#include <boost/type_traits/is_array.hpp>
#include <boost/type_traits/is_function.hpp>
#include <boost/type_traits/is_abstract.hpp>

#include <boost/mpl/or.hpp>

// STL
#include <cstddef>          // for std::size_t
#include <iosfwd>
#include <ios>              // for std::boolalpha
#include <climits>          // for CHAR_BIT

#ifdef BOOST_MSVC
# pragma warning(disable: 4127) // conditional expression is constant
#endif

#include <boost/test/detail/suppress_warnings.hpp>

//____________________________________________________________________________//

// ************************************************************************** //
// **************                    TOOL BOX                  ************** //
// ************************************************************************** //

namespace boost {
namespace test_tools {

typedef unit_test::const_string      const_string;

namespace { bool dummy_cond = false; }

// ************************************************************************** //
// **************                print_log_value               ************** //
// ************************************************************************** //

template<typename T>
struct print_log_value {
    void    operator()( std::ostream& ostr, T const& t )
    {
        // avoid warning: 'boost::test_tools::<unnamed>::dummy_cond' defined but not used
        if( ::boost::test_tools::dummy_cond ) {}

        typedef typename mpl::or_<is_array<T>,is_function<T>,is_abstract<T> >::type cant_use_nl;

        set_precision( ostr, cant_use_nl() );

        ostr << t; // by default print the value
    }

    void set_precision( std::ostream& ostr, mpl::false_ )
    {
        if( std::numeric_limits<T>::is_specialized && std::numeric_limits<T>::radix == 2 )
        {
            ostr.precision( 2 + std::numeric_limits<T>::digits * 301/1000 );
        }
        else if ( std::numeric_limits<T>::is_specialized && std::numeric_limits<T>::radix == 10 )
        {
#ifdef BOOST_NO_CXX11_NUMERIC_LIMITS
// (was BOOST_NO_NUMERIC_LIMITS_LOWEST but now deprecated).
// No support for std::numeric_limits<double>::max_digits10,
// so guess that a couple of guard digits more than digits10 will display any difference.
           ostr.precision( 2 + std::numeric_limits<T>::digits10 );
#else
  // std::numeric_limits<double>::max_digits10; IS supported.
  // Any noisy or guard digits needed to display any difference are included in max_digits10.
           ostr.precision( std::numeric_limits<T>::max_digits10 );
#endif
       }
       // else if T is not specialized for std::numeric_limits<>,
       // then will just get the default precision of 6 digits.
    }

    void set_precision( std::ostream&, mpl::true_ ) {}
};

//____________________________________________________________________________//

#if BOOST_WORKAROUND(__BORLANDC__, BOOST_TESTED_AT(0x564))
template<typename T, std::size_t N >
struct print_log_value< T[N] > {
    void    operator()( std::ostream& ostr, T const* t )
    {
        ostr << t;
    }
};
#endif

//____________________________________________________________________________//

template<>
struct BOOST_TEST_DECL print_log_value<bool> {
    void    operator()( std::ostream& ostr, bool t )
    {
         ostr << std::boolalpha << t;
    }
};

//____________________________________________________________________________//

template<>
struct BOOST_TEST_DECL print_log_value<char> {
    void    operator()( std::ostream& ostr, char t );
};

//____________________________________________________________________________//

template<>
struct BOOST_TEST_DECL print_log_value<unsigned char> {
    void    operator()( std::ostream& ostr, unsigned char t );
};

//____________________________________________________________________________//

template<>
struct BOOST_TEST_DECL print_log_value<char const*> {
    void    operator()( std::ostream& ostr, char const* t );
};

//____________________________________________________________________________//

template<>
struct BOOST_TEST_DECL print_log_value<wchar_t const*> {
    void    operator()( std::ostream& ostr, wchar_t const* t );
};

//____________________________________________________________________________//

namespace tt_detail {

// ************************************************************************** //
// **************              tools classification            ************** //
// ************************************************************************** //

enum check_type {
    CHECK_PRED,
    CHECK_MSG,
    CHECK_EQUAL,
    CHECK_NE,
    CHECK_LT,
    CHECK_LE,
    CHECK_GT,
    CHECK_GE,
    CHECK_CLOSE,
    CHECK_CLOSE_FRACTION,
    CHECK_SMALL,
    CHECK_BITWISE_EQUAL,
    CHECK_PRED_WITH_ARGS,
    CHECK_EQUAL_COLL,
    CHECK_BUILT_ASSERTION
};

enum tool_level {
    WARN, CHECK, REQUIRE, PASS
};

// ************************************************************************** //
// **************                 print_helper                 ************** //
// ************************************************************************** //
// Adds level of indirection to the output operation, allowing us to customize
// it for types that do not support operator << directly or for any other reason

template<typename T>
struct print_helper_t {
    explicit    print_helper_t( T const& t ) : m_t( t ) {}

    T const&    m_t;
};

//____________________________________________________________________________//

#if BOOST_WORKAROUND(__BORLANDC__, BOOST_TESTED_AT(0x564))
// Borland suffers premature pointer decay passing arrays by reference
template<typename T, std::size_t N >
struct print_helper_t< T[N] > {
    explicit    print_helper_t( T const * t ) : m_t( t ) {}

    T const *   m_t;
};
#endif

//____________________________________________________________________________//

template<typename T>
inline print_helper_t<T> print_helper( T const& t )
{
    return print_helper_t<T>( t );
}

//____________________________________________________________________________//

template<typename T>
inline std::ostream&
operator<<( std::ostream& ostr, print_helper_t<T> const& ph )
{
    print_log_value<T>()( ostr, ph.m_t );

    return ostr;
}

//____________________________________________________________________________//

// ************************************************************************** //
// **************            TOOL BOX Implementation           ************** //
// ************************************************************************** //

BOOST_TEST_DECL
bool check_impl( predicate_result const& pr, ::boost::unit_test::lazy_ostream const& assertion_descr,
                 const_string file_name, std::size_t line_num,
                 tool_level tl, check_type ct,
                 std::size_t num_args, ... );

//____________________________________________________________________________//

// This function adds level of indirection, but it makes sure we evaluate predicate
// arguments only once

#ifndef BOOST_TEST_PROD
#define TEMPL_PARAMS( z, m, dummy ) , typename BOOST_JOIN( Arg, m )

#define FUNC_PARAMS( z, m, dummy )                                                  \
 , BOOST_JOIN( Arg, m ) const& BOOST_JOIN( arg, m )                                 \
 , char const* BOOST_JOIN( BOOST_JOIN( arg, m ), _descr )                           \
/**/

#define PRED_PARAMS( z, m, dummy ) BOOST_PP_COMMA_IF( m ) BOOST_JOIN( arg, m )

#define ARG_INFO( z, m, dummy )                                                     \
 , BOOST_JOIN( BOOST_JOIN( arg, m ), _descr )                                       \
 , &static_cast<const unit_test::lazy_ostream&>(unit_test::lazy_ostream::instance() \
        << ::boost::test_tools::tt_detail::print_helper( BOOST_JOIN( arg, m ) ))    \
/**/

#define IMPL_FRWD( z, n, dummy )                                                    \
template<typename Pred                                                              \
         BOOST_PP_REPEAT_ ## z( BOOST_PP_ADD( n, 1 ), TEMPL_PARAMS, _ )>            \
inline bool                                                                         \
check_frwd( Pred P, unit_test::lazy_ostream const& assertion_descr,                 \
            const_string file_name, std::size_t line_num,                           \
            tool_level tl, check_type ct                                            \
            BOOST_PP_REPEAT_ ## z( BOOST_PP_ADD( n, 1 ), FUNC_PARAMS, _ )           \
)                                                                                   \
{                                                                                   \
    return                                                                          \
    check_impl( P( BOOST_PP_REPEAT_ ## z( BOOST_PP_ADD( n, 1 ), PRED_PARAMS, _ ) ), \
                assertion_descr, file_name, line_num, tl, ct,                       \
                BOOST_PP_ADD( n, 1 )                                                \
                BOOST_PP_REPEAT_ ## z( BOOST_PP_ADD( n, 1 ), ARG_INFO, _ )          \
    );                                                                              \
}                                                                                   \
/**/

#ifndef BOOST_TEST_MAX_PREDICATE_ARITY
#define BOOST_TEST_MAX_PREDICATE_ARITY 5
#endif

BOOST_PP_REPEAT( BOOST_TEST_MAX_PREDICATE_ARITY, IMPL_FRWD, _ )

#undef TEMPL_PARAMS
#undef FUNC_PARAMS
#undef PRED_INFO
#undef ARG_INFO
#undef IMPL_FRWD

#endif

//____________________________________________________________________________//

template <class Left, class Right>
predicate_result equal_impl( Left const& left, Right const& right )
{
    return left == right;
}

//____________________________________________________________________________//

predicate_result        BOOST_TEST_DECL equal_impl( char const* left, char const* right );
inline predicate_result equal_impl( char* left, char const* right ) { return equal_impl( static_cast<char const*>(left), static_cast<char const*>(right) ); }
inline predicate_result equal_impl( char const* left, char* right ) { return equal_impl( static_cast<char const*>(left), static_cast<char const*>(right) ); }
inline predicate_result equal_impl( char* left, char* right )       { return equal_impl( static_cast<char const*>(left), static_cast<char const*>(right) ); }

#if !defined( BOOST_NO_CWCHAR )
predicate_result        BOOST_TEST_DECL equal_impl( wchar_t const* left, wchar_t const* right );
inline predicate_result equal_impl( wchar_t* left, wchar_t const* right ) { return equal_impl( static_cast<wchar_t const*>(left), static_cast<wchar_t const*>(right) ); }
inline predicate_result equal_impl( wchar_t const* left, wchar_t* right ) { return equal_impl( static_cast<wchar_t const*>(left), static_cast<wchar_t const*>(right) ); }
inline predicate_result equal_impl( wchar_t* left, wchar_t* right )       { return equal_impl( static_cast<wchar_t const*>(left), static_cast<wchar_t const*>(right) ); }
#endif

//____________________________________________________________________________//

struct equal_impl_frwd {
    template <typename Left, typename Right>
    inline predicate_result
    call_impl( Left const& left, Right const& right, mpl::false_ ) const
    {
        return equal_impl( left, right );
    }

    template <typename Left, typename Right>
    inline predicate_result
    call_impl( Left const& left, Right const& right, mpl::true_ ) const
    {
        return (*this)( right, &left[0] );
    }

    template <typename Left, typename Right>
    inline predicate_result
    operator()( Left const& left, Right const& right ) const
    {
        typedef typename is_array<Left>::type left_is_array;
        return call_impl( left, right, left_is_array() );
    }
};

//____________________________________________________________________________//

struct ne_impl {
    template <class Left, class Right>
    predicate_result operator()( Left const& left, Right const& right )
    {
        return !equal_impl_frwd()( left, right );
    }
};

//____________________________________________________________________________//

struct lt_impl {
    template <class Left, class Right>
    predicate_result operator()( Left const& left, Right const& right )
    {
        return left < right;
    }
};

//____________________________________________________________________________//

struct le_impl {
    template <class Left, class Right>
    predicate_result operator()( Left const& left, Right const& right )
    {
        return left <= right;
    }
};

//____________________________________________________________________________//

struct gt_impl {
    template <class Left, class Right>
    predicate_result operator()( Left const& left, Right const& right )
    {
        return left > right;
    }
};

//____________________________________________________________________________//

struct ge_impl {
    template <class Left, class Right>
    predicate_result operator()( Left const& left, Right const& right )
    {
        return left >= right;
    }
};

//____________________________________________________________________________//

struct equal_coll_impl {
    template <typename Left, typename Right>
    predicate_result operator()( Left left_begin, Left left_end, Right right_begin, Right right_end )
    {
        predicate_result    pr( true );
        std::size_t         pos = 0;

        for( ; left_begin != left_end && right_begin != right_end; ++left_begin, ++right_begin, ++pos ) {
            if( *left_begin != *right_begin ) {
                pr = false;
                pr.message() << "\nMismatch in a position " << pos << ": "  << *left_begin << " != " << *right_begin;
            }
        }

        if( left_begin != left_end ) {
            std::size_t r_size = pos;
            while( left_begin != left_end ) {
                ++pos;
                ++left_begin;
            }

            pr = false;
            pr.message() << "\nCollections size mismatch: " << pos << " != " << r_size;
        }

        if( right_begin != right_end ) {
            std::size_t l_size = pos;
            while( right_begin != right_end ) {
                ++pos;
                ++right_begin;
            }

            pr = false;
            pr.message() << "\nCollections size mismatch: " << l_size << " != " << pos;
        }

        return pr;
    }
};

//____________________________________________________________________________//

struct bitwise_equal_impl {
    template <class Left, class Right>
    predicate_result    operator()( Left const& left, Right const& right )
    {
        predicate_result    pr( true );

        std::size_t left_bit_size  = sizeof(Left)*CHAR_BIT;
        std::size_t right_bit_size = sizeof(Right)*CHAR_BIT;

        static Left const leftOne( 1 );
        static Right const rightOne( 1 );

        std::size_t total_bits = left_bit_size < right_bit_size ? left_bit_size : right_bit_size;

        for( std::size_t counter = 0; counter < total_bits; ++counter ) {
            if( ( left & ( leftOne << counter ) ) != ( right & ( rightOne << counter ) ) ) {
                pr = false;
                pr.message() << "\nMismatch in a position " << counter;
            }
        }

        if( left_bit_size != right_bit_size ) {
            pr = false;
            pr.message() << "\nOperands bit sizes mismatch: " << left_bit_size << " != " << right_bit_size;
        }

        return pr;
    }
};

//____________________________________________________________________________//

bool BOOST_TEST_DECL is_defined_impl( const_string symbol_name, const_string symbol_value );

//____________________________________________________________________________//

// ************************************************************************** //
// **************                 context_frame                ************** //
// ************************************************************************** //

struct BOOST_TEST_DECL context_frame {
    explicit    context_frame( ::boost::unit_test::lazy_ostream const& context_descr );
    ~context_frame();

    operator    bool();

private:
    // Data members
    int         m_frame_id;
};

} // namespace tt_detail
} // namespace test_tools
} // namespace boost

#include <boost/test/detail/enable_warnings.hpp>

#endif // BOOST_TEST_TEST_TOOLS_IMPL_HPP_012705GER
