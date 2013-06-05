# Copyright Dave Abrahams 2012. Distributed under the Boost
# Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

boost_add_multiple_tests(
  LINK
    conv_example.cpp
    predicate_example.cpp
    find_example.cpp
    replace_example.cpp
    rle_example.cpp
    trim_example.cpp
    split_example.cpp
)

boost_add_test(
  regex_example
  LINK
    regex_example.cpp 
  LINK_BOOST_LIBRARIES
    regex
)
