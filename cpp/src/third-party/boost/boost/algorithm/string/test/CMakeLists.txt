# Copyright Dave Abrahams 2012. Distributed under the Boost
# Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

boost_add_multiple_tests(
  RUN
    trim_test.cpp
    conv_test.cpp
    predicate_test.cpp
    find_test.cpp
    split_test.cpp
    join_test.cpp
    replace_test.cpp
    find_format_test.cpp
)

boost_add_test(
  regex_test
  LINK
    regex_test.cpp 
  LINK_BOOST_LIBRARIES
    regex
)
