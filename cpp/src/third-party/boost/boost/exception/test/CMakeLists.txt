# Distributed under the Boost Software License, Version 1.0.
# See http://www.boost.org/LICENSE_1_0.txt

include(Ryppl)
ryppl_find_and_use_package(Boost COMPONENTS test utility NO_MODULE)
include(${Boost_USE_FILE})
include(RypplTestSuite)

# todo
#project
#    : requirements
#        <link>static
#        <exception-handling>on
#    ;

#to_string
ryppl_test_suite(
  RUN
    is_output_streamable_test.cpp
		has_to_string_test.cpp
		to_string_test.cpp
		to_string_stub_test.cpp
	COMPILE_FAIL
		to_string_fail.cpp
  )

#exception
ryppl_test_suite(
  RUN
    1-throw_exception_test.cpp
		2-throw_exception_no_exceptions_test.cpp
		3-throw_exception_no_integration_test.cpp
		4-throw_exception_no_both_test.cpp
		cloning_test.cpp
#		copy_exception_test.cpp /boost//thread : : : <threading>multi ;
		unknown_exception_test.cpp
		exception_test.cpp
#		enable_error_info_test.cpp helper1.cpp ;
#		throw_exception_test.cpp helper2.cpp ;
		errno_test.cpp
		error_info_test.cpp
		diagnostic_information_test.cpp
		refcount_ptr_test.cpp
		current_exception_cast_test.cpp
#		no_exceptions_test.cpp : : : <exception-handling>off ;
		errinfos_test.cpp
#		exception_ptr_test.cpp/<define>BOOST_ENABLE_NON_INTRUSIVE_EXCEPTION_PTR /boost/exception /boost//thread : : : <threading>multi : non_intrusive_exception_ptr_test ;
# 	exception_ptr_test.cpp /boost//thread : : : <threading>multi ;
	COMPILE_FAIL
		exception_fail.cpp
		throw_exception_fail.cpp
		error_info_const_fail.cpp
  )

#headers
ryppl_test_suite(
  COMPILE
		exception_ptr_hpp_test.cpp
		diagnostic_information_hpp_test.cpp
		error_info_hpp_test.cpp
		get_error_info_hpp_test.cpp
		info_hpp_test.cpp
		info_tuple_hpp_test.cpp
		to_string_hpp_test.cpp
		to_string_stub_hpp_test.cpp
		all_hpp_test.cpp
		current_exception_cast_hpp_test.cpp
		errinfo_api_function_hpp_test.cpp
		errinfo_at_line_hpp_test.cpp
		errinfo_errno_hpp_test.cpp
		errinfo_file_handle_hpp_test.cpp
		errinfo_file_name_hpp_test.cpp
		errinfo_file_open_mode_hpp_test.cpp
		errinfo_nested_exception_hpp_test.cpp
		errinfo_type_info_name_hpp_test.cpp
)
