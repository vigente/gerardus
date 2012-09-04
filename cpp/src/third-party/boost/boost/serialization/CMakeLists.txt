# Distributed under the Boost Software License, Version 1.0.
# See http://www.boost.org/LICENSE_1_0.txt

cmake_minimum_required(VERSION 2.8.8)
project(BoostSerialization CXX)

include(Ryppl)

ryppl_find_and_use_package(BoostUtility)
ryppl_find_and_use_package(BoostSpiritClassic)
include_directories(include)

if(MSVC)
  add_definitions(
    /Gy
    /wd4244 # converting 'std::streamsize' in 'size_t'
    /D_CRT_SECURE_NO_WARNINGS
    /D_SCL_SECURE_NO_WARNINGS
    )
endif(MSVC)

# if(GCC)
#   add_definitions(-fvisibility=hidden)
# endif(GCC)

add_library(boost_serialization
  src/basic_archive
  src/basic_iarchive
  src/basic_iserializer
  src/basic_oarchive
  src/basic_oserializer
  src/basic_pointer_iserializer
  src/basic_pointer_oserializer
  src/basic_serializer_map
  src/basic_text_iprimitive
  src/basic_text_oprimitive
  src/basic_xml_archive
  src/binary_iarchive
  src/binary_oarchive
  src/extended_type_info
  src/extended_type_info_typeid
  src/extended_type_info_no_rtti
  src/polymorphic_iarchive
  src/polymorphic_oarchive
  src/stl_port
  src/text_iarchive
  src/text_oarchive
  src/void_cast
  src/archive_exception
  src/xml_grammar
  src/xml_iarchive
  src/xml_oarchive
  src/xml_archive_exception
  src/shared_ptr_helper
  )

add_library(boost_wserialization
  src/basic_text_wiprimitive
  src/basic_text_woprimitive
  src/text_wiarchive
  src/text_woarchive
  src/utf8_codecvt_facet
  src/xml_wgrammar
  src/xml_wiarchive
  src/xml_woarchive
  src/codecvt_null
  )

target_link_libraries(boost_wserialization
  boost_serialization
  )

ryppl_export(
  TARGETS
    boost_serialization
    boost_wserialization
  DEPENDS
    BoostArray
    BoostConfig
    BoostCore
    BoostIO
    BoostInteger
    BoostIterator
    BoostMPL
    BoostOptional
    BoostPreprocessor
    BoostSmartPtr
    BoostSpiritClassic
    BoostStaticAssert
    BoostTypeTraits
    BoostUtility
    BoostVariant
  INCLUDE_DIRECTORIES
    include
  )
