# Distributed under the Boost Software License, Version 1.0.
# See http://www.boost.org/LICENSE_1_0.txt

include(Ryppl)
ryppl_find_and_use_package(Boost COMPONENTS test NO_MODULE)
include(${Boost_USE_FILE})
include(RypplTestSuite)

add_definitions(-DBOOST_ALL_DYN_LINK)
add_definitions(-DBOOST_DATE_TIME_POSIX_TIME_STD_CONFIG)
add_definitions(-DBOOST_ALL_NO_LIB)
#local DATE_TIME_DYNAMIC_PROPERTIES = <define>BOOST_ALL_DYN_LINK <runtime-link>shared #<define>BOOST_DATE_TIME_POSIX_TIME_STD_CONFIG <define>BOOST_ALL_NO_LIB ;

# date_time_core
ryppl_test_suite(
  RUN
    testint_adapter.cpp
		testtime_resolution_traits.cpp
		testwrapping_int.cpp
		testconstrained_value.cpp
		testgregorian_calendar.cpp
		testgeneric_period.cpp    
  )

# date_time_gregorian
ryppl_test_suite(
  RUN
		gregorian/testdate.cpp # $(DATE_TIME_PROPERTIES)
		gregorian/testdate_duration.cpp
		gregorian/testgreg_durations.cpp
		gregorian/testperiod.cpp
		gregorian/testdate_iterator.cpp
		gregorian/testformatters.cpp
###  streaming
		gregorian/testdate_facet_new.cpp
		gregorian/testdate_input_facet.cpp
###
		gregorian/testgenerators.cpp
		gregorian/testgreg_cal.cpp
		gregorian/testgreg_day.cpp
		gregorian/testgreg_month.cpp
		gregorian/testgreg_year.cpp
	LINK_LIBRARIES
    boost_date_time-static
  )

add_definitions(-DBOOST_DATE_TIME_POSIX_TIME_STD_CONFIG)
add_definitions(-DBOOST_DATE_TIME_STATIC_LINK)
add_definitions(-DBOOST_ALL_NO_LIB)
#local DATE_TIME_PROPERTIES = <define>BOOST_DATE_TIME_POSIX_TIME_STD_CONFIG <define>BOOST_DATE_TIME_STATIC_LINK
  # FIXME 
  #std::locale-support 
#  <define>BOOST_ALL_NO_LIB ;

# date_time_gregorian_dll
ryppl_test_suite(
  RUN
		gregorian/testdate.cpp # $(DATE_TIME_DYNAMIC_PROPERTIES)
		gregorian/testdate_duration.cpp
		gregorian/testperiod.cpp
		gregorian/testdate_iterator.cpp
		gregorian/testformatters.cpp
###  streaming 
		gregorian/testdate_facet_new.cpp
		gregorian/testdate_input_facet.cpp
###
		gregorian/testgenerators.cpp
		gregorian/testgreg_cal.cpp
		gregorian/testgreg_day.cpp
		gregorian/testgreg_month.cpp
		gregorian/testgreg_year.cpp
		gregorian/testgreg_durations.cpp
	LINK_LIBRARIES
    boost_date_time
  )

