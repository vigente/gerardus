//  (C) Copyright Gennadiy Rozental 2011.
//  Distributed under the Boost Software License, Version 1.0.
//  (See accompanying file LICENSE_1_0.txt or copy at 
//  http://www.boost.org/LICENSE_1_0.txt)

//  See http://www.boost.org/libs/test for the library home page.
//
//  File        : $RCSfile$
//
//  Version     : $Revision: 74663 $
//
//  Description : defines framework for automated assertion construction
// ***************************************************************************

#ifndef BOOST_TEST_TOOLS_ASSERTION_HPP_100911GER
#define BOOST_TEST_TOOLS_ASSERTION_HPP_100911GER

// Boost.Test

// Boost
#include <boost/mpl/assert.hpp>
#include <boost/utility/declval.hpp>

#include <boost/test/detail/suppress_warnings.hpp>

//____________________________________________________________________________//

namespace boost {
namespace test_tools {
namespace assertion {

// ************************************************************************** //
// **************             assertion::expression            ************** //
// ************************************************************************** //

class expression {
public:
    // expression interface
    virtual predicate_result    evaluate() const = 0;
};

// ************************************************************************** //
// **************             assertion::operators             ************** //
// ************************************************************************** //

namespace op {

enum id {
    // precedence 4: ->*, .*
    MEMP,
    // precedence 5: *, /, %
    MUL, DIV, MOD,
    // precedence 6: +, -
    ADD, SUB,
    // precedence 7: << , >>
    LSH, RSH,
    // precedence 8: <, <=, > and >=
    LT, LE, GT, GE,
    // precedence 9: == and !=
    EQ, NE,
    // precedence 10: bitwise AND
    BAND,
    // precedence 11: bitwise XOR
    XOR,
    // precedence 12: bitwise OR
    BOR,
    // precedence 13: logical AND
    //  disabled
    // precedence 14: logical OR
    //  disabled
    // precedence 15: ternary conditional
    //  disabled

    // precedence 16: = and OP= operators
    SET, IADD, ISUB, IMUL, IDIV, IMOD,
    ILSH, IRSH,
    IAND, IXOR, IOR
    // precedence 17: throw operator
    // not supported
    // precedence 18: comma
    // not supported
};

#ifndef BOOST_NO_DECLTYPE

#define BOOST_TEST_FOR_EACH_OP( action )    \
    action(->*, MEMP )                      \
                                            \
    action( * , MUL  )                      \
    action( / , DIV  )                      \
    action( % , MOD  )                      \
                                            \
    action( + , ADD  )                      \
    action( - , SUB  )                      \
                                            \
    action( <<, LSH  )                      \
    action( >>, RSH  )                      \
                                            \
    action( < , LT   )                      \
    action( <=, LE   )                      \
    action( > , GT   )                      \
    action( >=, GE   )                      \
                                            \
    action( ==, EQ   )                      \
    action( !=, NE   )                      \
                                            \
    action( & , BAND )                      \
    action( ^ , XOR  )                      \
    action( | , BOR  )                      \
/**/

#else

#define BOOST_TEST_FOR_EACH_OP( action )    \
    action( < , LT   )                      \
    action( <=, LE   )                      \
    action( > , GT   )                      \
    action( >=, GE   )                      \
                                            \
    action( ==, EQ   )                      \
    action( !=, NE   )                      \
/**/

#endif

#define BOOST_TEST_FOR_EACH_MUT_OP( action )\
    action( = , SET  )                      \
    action( +=, IADD )                      \
    action( -=, ISUB )                      \
    action( *=, IMUL )                      \
    action( /=, IDIV )                      \
    action( %=, IMOD )                      \
    action(<<=, ILSH )                      \
    action(>>=, IRSH )                      \
    action( &=, IAND )                      \
    action( ^=, IXOR )                      \
    action( |=, IOR  )                      \
/**/

// ************************************************************************** //
// **************          assertion::operator traits          ************** //
// ************************************************************************** //

template<id ID>
struct traits {
    static char const* reverse( char const* direct ) { return direct; }
};

template<> struct traits<EQ> { static char const* reverse( char const* ) { return "!="; } };
template<> struct traits<NE> { static char const* reverse( char const* ) { return "=="; } };

template<> struct traits<LT> { static char const* reverse( char const* ) { return ">="; } };
template<> struct traits<LE> { static char const* reverse( char const* ) { return ">"; } };
template<> struct traits<GT> { static char const* reverse( char const* ) { return "<="; } };
template<> struct traits<GE> { static char const* reverse( char const* ) { return "<"; } };

} // anmespace op

// ************************************************************************** //
// **************         assertion::expression_result         ************** //
// ************************************************************************** //

namespace detail {

template<typename PrevExprType,typename Rhs,op::id OP>
class expression_result;

#define DEFINE_OP_EXPRESSION_RESULT( OPER, ID )             \
template<typename PrevExprType,typename Rhs>                \
class expression_result<PrevExprType,Rhs,op::ID> {          \
    typedef typename PrevExprType::result_type Lhs;         \
public:                                                     \
                                                            \
    typedef DEDUCE_RESULT_TYPE( OPER ) type;                \
                                                            \
    static type                                             \
    produce( Lhs const& lhs, Rhs const& rhs)                \
    {                                                       \
        return lhs OPER rhs;                                \
    }                                                       \
                                                            \
    static void                                             \
    report( std::ostream& ostr,                             \
            PrevExprType const& lhs, Rhs const& rhs)        \
    {                                                       \
        lhs.report( ostr );                                 \
        ostr << op::traits<op::ID>::reverse( #OPER )        \
             << rhs;                                        \
    }                                                       \
};                                                          \
/**/

////////////////////////////////////////////////////////////////

#ifndef BOOST_NO_DECLTYPE
#define DEDUCE_RESULT_TYPE( OPER ) decltype(boost::declval<Lhs>() OPER boost::declval<Rhs>() ) 
#else
#define DEDUCE_RESULT_TYPE( OPER ) predicate_result 
#endif

////////////////////////////////////////////////////////////////

BOOST_TEST_FOR_EACH_OP( DEFINE_OP_EXPRESSION_RESULT )

#undef DEDUCE_RESULT_TYPE
#undef DEFINE_OP_EXPRESSION_RESULT

} // namespace detail

////////////////////////////////////////////////////////////////

// ************************************************************************** //
// **************           assertion::expression_op           ************** //
// ************************************************************************** //
// Defines expression operators

template<typename Lhs,typename Rhs,op::id OP> class binary_expr;

template<typename ExprType>
class expression_op {
public:

#define ADD_OP_SUPPORT( OPER, ID )                          \
    template<typename T>                                    \
    binary_expr<ExprType,T,op::ID>                          \
    operator OPER( T const& rhs ) const                     \
    {                                                       \
        return binary_expr<ExprType,T,op::ID>(              \
            *static_cast<ExprType const*>(this),            \
            rhs );                                          \
    }                                                       \
/**/

    BOOST_TEST_FOR_EACH_OP( ADD_OP_SUPPORT )
#undef ADD_OP_SUPPORT

    // Disabled operators
    template<typename T>
    ExprType&
    operator ||( T const& rhs )
    {
        BOOST_MPL_ASSERT_MSG(false, CANT_USE_LOGICAL_OPERATOR_OR_WITHIN_THIS_TESTING_TOOL, () );

        return *static_cast<ExprType*>(this);
    }

    template<typename T>
    ExprType&
    operator &&( T const& rhs )
    {
        BOOST_MPL_ASSERT_MSG(false, CANT_USE_LOGICAL_OPERATOR_AND_WITHIN_THIS_TESTING_TOOL, () );

        return *static_cast<ExprType*>(this);
    }

    operator bool()
    {
        BOOST_MPL_ASSERT_MSG(false, CANT_USE_TERNARY_OPERATOR_WITHIN_THIS_TESTING_TOOL, () );

        return false;
    }
};

// ************************************************************************** //
// **************            assertion::value_expr             ************** //
// ************************************************************************** //
// simple value expression

template<typename T>
class value_expr : public expression, public expression_op<value_expr<T> > {
public:
    // Public types
    typedef T                   result_type;

    // Constructor
#ifndef BOOST_NO_RVALUE_REFERENCES
    explicit                    value_expr( T&& val )
    : m_value( std::forward<T>(val) )
#else
    explicit                    value_expr( T const& val )
    : m_value( val )
#endif
    {}

    // Specific expresson interface
    T const&                    value() const
    {
        return m_value;
    }
    void                        report( std::ostream& ostr ) const
    {
        ostr << m_value;
    }

    // Mutating operators
#define ADD_OP_SUPPORT( OPER, ID )                          \
    template<typename U>                                    \
    value_expr<T>&                                          \
    operator OPER( U const& rhs )                           \
    {                                                       \
        m_value OPER rhs;                                   \
                                                            \
        return *this;                                       \
    }                                                       \
/**/

    BOOST_TEST_FOR_EACH_MUT_OP( ADD_OP_SUPPORT )
#undef ADD_OP_SUPPORT

private:
    template<typename U>
    static void format_message( wrap_stringstream& ostr, U const& v )    { ostr << "(bool)" << v << " is false"; }
    static void format_message( wrap_stringstream& ostr, bool v )        {}

    // expression interface
    virtual predicate_result    evaluate() const
    {
        predicate_result res( value() );
        if( !res )
            format_message( res.message(), value() );

        return res;
    }

    // Data members
    T                           m_value;
};

// ************************************************************************** //
// **************         assertion::binary_expr_expr          ************** //
// ************************************************************************** //
// binary expression

template<typename Lhs,typename Rhs,op::id OP>
class binary_expr : public expression, public expression_op<binary_expr<Lhs,Rhs,OP> > {
    typedef detail::expression_result<Lhs,Rhs,OP>   result;
public:
    // Public types
    typedef typename result::type result_type;

    // Constructor
    binary_expr( Lhs const& lhs, Rhs const& rhs )
    : m_lhs( lhs )
    , m_rhs( rhs )
    {}

    // Specifica expression interface
    result_type                 value() const
    {
        return result::produce( m_lhs.value(), m_rhs );
    }
    void                        report( std::ostream& ostr ) const
    {
        return result::report( ostr, m_lhs, m_rhs );
    }

private:
    // expression interface
    virtual predicate_result    evaluate() const
    {
        predicate_result res( value() );
        if( !res )
            report( res.message().stream() );

        return res;
    }

    // Data members
    Lhs                         m_lhs;
    Rhs                         m_rhs;
};

// ************************************************************************** //
// **************               assertion::seed                ************** //
// ************************************************************************** //
// seed added ot the input expression to form an assertion expression

class seed {
public:
    // ->* is highest precedence left to right operator
    template<typename T>
    value_expr<T>
#ifndef BOOST_NO_RVALUE_REFERENCES
    operator->*( T&& v ) const
    {
        return value_expr<T>( std::forward<T>( v ) );
    }
#else
    operator->*( T const& v )  const
    {
        return value_expr<T>( v );
    }
#endif

};

#undef BOOST_TEST_FOR_EACH_OP

} // namespace assertion
} // namespace test_tools
} // namespace boost

#include <boost/test/detail/enable_warnings.hpp>

#endif // BOOST_TEST_TOOLS_ASSERTION_HPP_100911GER

