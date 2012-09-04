# Distributed under the Boost Software License, Version 1.0.
# See http://www.boost.org/LICENSE_1_0.txt

include(Ryppl)
include(RypplTestSuite)

ryppl_find_and_use_package(BoostCore)
ryppl_find_and_use_package(BoostMPL)

include_directories(../example) # needed for Boost.PP file iteration with some compilers

ryppl_test_suite(

  COMPILE

    # Classification
    classification/is_function.cpp
    classification/is_function_pointer.cpp
    classification/is_function_reference.cpp
    classification/is_member_function_pointer.cpp
    classification/is_member_object_pointer.cpp
    classification/is_callable_builtin.cpp
    classification/is_nonmember_callable_builtin.cpp
    classification/is_member_pointer.cpp
    classification/is_cv_mem_func_ptr.cpp
    classification/is_variadic.cpp
    classification/is_cv_pointer.cpp
    #classification/is_cv_function.cpp

    # Decomposition
    decomposition/components.cpp
    decomposition/result_type.cpp
    decomposition/function_arity.cpp
    decomposition/parameter_types.cpp
    decomposition/components_seq.cpp
    decomposition/class_type_transform.cpp

    # Synthesis
    synthesis/function_type.cpp
    synthesis/function_pointer.cpp
    synthesis/function_reference.cpp
    synthesis/member_function_pointer.cpp
    synthesis/member_object_pointer.cpp
    synthesis/transformation.cpp
    synthesis/mem_func_ptr_cv1.cpp
    synthesis/mem_func_ptr_cv2.cpp
    synthesis/mem_func_ptr_cv_ptr_to_this.cpp
    synthesis/variadic_function_synthesis.cpp
    #synthesis/cv_function_synthesis.cpp

    # Reconfiguration
    reconfiguration/preprocessing_mode.cpp
    reconfiguration/partial_arity_preprocessing.cpp
    reconfiguration/cc_preprocessing.cpp

    # Custom calling conventions
    custom_ccs/nonmember_ccs.cpp
    custom_ccs/nonmember_ccs_exact.cpp
    #custom_ccs/member_ccs.cpp
    #custom_ccs/member_ccs_exact.cpp

    # Code from the examples
    ../example/interpreter_example.cpp
    ../example/result_of_example.cpp
    ../example/interface_example.cpp
    ../example/fast_mem_fn_example.cpp
    ../example/macro_type_args_example.cpp

  COMPILE_FAIL

    # Decomposition
    decomposition/result_type_fail.cpp
    decomposition/parameter_types_fail.cpp
    decomposition/function_arity_fail.cpp
  NO_SINGLE_TARGET
  )
