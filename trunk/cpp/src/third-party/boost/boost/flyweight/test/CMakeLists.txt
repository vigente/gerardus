# Distributed under the Boost Software License, Version 1.0.
# See http://www.boost.org/LICENSE_1_0.txt

include(Ryppl)
include(RypplTestSuite)

ryppl_find_and_use_package(BoostContainer)
ryppl_find_and_use_package(BoostFlyweight)
ryppl_find_and_use_package(BoostInterprocess)
ryppl_find_and_use_package(BoostSerialization)
ryppl_find_package(Threads)

# assoc_cont_factory
ryppl_test_suite(
  RUN
    test_assoc_cont_factory.cpp
  ADDITIONAL_SOURCES
    test_assoc_cont_fact_main.cpp
  LINK_LIBRARIES
    ${CMAKE_THREAD_LIBS_INIT}
  )

# basic
ryppl_test_suite(
  RUN
    test_basic.cpp
  ADDITIONAL_SOURCES
    test_basic_main.cpp
  LINK_LIBRARIES
    ${CMAKE_THREAD_LIBS_INIT}
  )

# custom_factory
ryppl_test_suite(
  RUN
    test_custom_factory.cpp
  ADDITIONAL_SOURCES
    test_custom_factory_main.cpp
  LINK_LIBRARIES
    ${CMAKE_THREAD_LIBS_INIT}
  )

# init
ryppl_test_suite(
  RUN
    test_init.cpp
  ADDITIONAL_SOURCES
    test_init_main.cpp
  LINK_LIBRARIES
    ${CMAKE_THREAD_LIBS_INIT}
  )

# intermod_holder <threading>multi   
ryppl_test_suite(
  RUN
    test_intermod_holder.cpp
  ADDITIONAL_SOURCES
    test_intermod_holder_main.cpp
  LINK_LIBRARIES
    ${CMAKE_THREAD_LIBS_INIT}
    intermod_holder_dll
  )

# multictor
ryppl_test_suite(
  RUN
    test_multictor.cpp
  ADDITIONAL_SOURCES
    test_multictor_main.cpp
  LINK_LIBRARIES
    ${CMAKE_THREAD_LIBS_INIT}
  )

# no_locking
ryppl_test_suite(
  RUN
    test_no_locking.cpp
  ADDITIONAL_SOURCES
    test_no_locking_main.cpp
  LINK_LIBRARIES
    ${CMAKE_THREAD_LIBS_INIT}
  )

# no_tracking
ryppl_test_suite(
  RUN
    test_no_tracking.cpp
  ADDITIONAL_SOURCES
    test_no_tracking_main.cpp
  LINK_LIBRARIES
    ${CMAKE_THREAD_LIBS_INIT}
  )

# set_factory
ryppl_test_suite(
  RUN
    test_set_factory.cpp
  ADDITIONAL_SOURCES
    test_set_factory_main.cpp
  LINK_LIBRARIES
    ${CMAKE_THREAD_LIBS_INIT}
  )

add_library(intermod_holder_dll SHARED EXCLUDE_FROM_ALL
  intermod_holder_dll.cpp
  )
set_property(TARGET intermod_holder_dll APPEND PROPERTY
  COMPILE_DEFINITIONS "BOOST_FLYWEIGHT_TEST_INTERMOD_HOLDER_DLL_SOURCE=1"
  )
if(UNIX AND NOT APPLE)
  target_link_libraries(intermod_holder_dll rt)
endif()
