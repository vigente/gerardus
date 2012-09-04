# Distributed under the Boost Software License, Version 1.0.
# See http://www.boost.org/LICENSE_1_0.txt

include(Ryppl)
include(RypplTestSuite)

ryppl_find_and_use_package(BoostMPL)
ryppl_find_and_use_package(BoostAny)

ryppl_test_suite(
  RUN
    algorithm/all.cpp
    algorithm/any.cpp
    algorithm/clear.cpp
    algorithm/copy.cpp
    algorithm/count.cpp
    algorithm/count_if.cpp
    algorithm/erase.cpp
    algorithm/erase_key.cpp
    algorithm/filter.cpp
    algorithm/filter_if.cpp
    algorithm/find.cpp
    algorithm/find_if.cpp
    algorithm/fold.cpp
    algorithm/fold2.cpp
    algorithm/for_each.cpp
    algorithm/insert.cpp
    algorithm/insert_range.cpp
    algorithm/iter_fold.cpp
    algorithm/none.cpp
    algorithm/pop_back.cpp
    algorithm/pop_front.cpp
    algorithm/push_back.cpp
    algorithm/push_front.cpp
    algorithm/remove.cpp
    algorithm/remove_if.cpp
    algorithm/replace.cpp
    algorithm/replace_if.cpp
    algorithm/reverse_fold.cpp
    algorithm/reverse_iter_fold.cpp
    algorithm/reverse.cpp
    algorithm/transform.cpp
    algorithm/join.cpp
    algorithm/zip.cpp
    algorithm/zip2.cpp
    algorithm/zip_ignore.cpp

    sequence/as_list.cpp
    sequence/as_map.cpp
    sequence/as_set.cpp
    sequence/as_vector.cpp
    sequence/boost_tuple.cpp
    sequence/cons.cpp
    sequence/filter_view.cpp
    sequence/io.cpp
    sequence/iterator_range.cpp
    sequence/joint_view.cpp
    sequence/list_comparison.cpp
    sequence/list_construction.cpp
    sequence/list_copy.cpp
    sequence/list_iterator.cpp
    sequence/list_make.cpp
    sequence/list_misc.cpp
    sequence/list_mutate.cpp
    sequence/list_tie.cpp
    sequence/list_value_at.cpp
    sequence/deque_comparison.cpp
    sequence/deque_construction.cpp
    sequence/deque_copy.cpp
    sequence/deque_iterator.cpp
    sequence/deque_make.cpp
    sequence/deque_misc.cpp
    sequence/deque_mutate.cpp
    sequence/deque_tie.cpp
    sequence/deque_value_at.cpp
    sequence/front_extended_deque.cpp
    sequence/back_extended_deque.cpp
    sequence/make_list.cpp
    sequence/make_vector.cpp
    sequence/map.cpp
    sequence/map_tie.cpp
    sequence/nview.cpp
    sequence/reverse_view.cpp
    sequence/set.cpp
    sequence/single_view.cpp
    sequence/std_pair.cpp
    sequence/boost_array.cpp
    sequence/array.cpp
    sequence/tuple_comparison.cpp
    sequence/tuple_construction.cpp
    sequence/tuple_copy.cpp
    sequence/tuple_element.cpp
    sequence/tuple_make.cpp
    sequence/tuple_misc.cpp
    sequence/tuple_mutate.cpp
    sequence/tuple_tie.cpp
    sequence/tr1_tuple_auto_conv.cpp
    sequence/transform_view.cpp
    sequence/vector_comparison.cpp
    sequence/vector_construction.cpp
    sequence/vector_copy.cpp
    sequence/vector_iterator.cpp
    sequence/vector_make.cpp
    sequence/vector_misc.cpp
    sequence/vector_mutate.cpp
    sequence/vector_n.cpp
    sequence/vector_tie.cpp
    sequence/vector_value_at.cpp
    sequence/zip_view.cpp
    sequence/zip_view2.cpp
    sequence/zip_view_ignore.cpp
    sequence/repetitive_view.cpp
    sequence/deduce_sequence.cpp
    sequence/adapt_adt_named.cpp
    sequence/adapt_adt.cpp
    sequence/adapt_assoc_adt_named.cpp
    sequence/adapt_assoc_adt.cpp
    sequence/adapt_assoc_struct_named.cpp
    sequence/adapt_assoc_struct.cpp
    sequence/adapt_assoc_tpl_adt.cpp
    sequence/adapt_assoc_tpl_struct.cpp
    sequence/adapt_struct_named.cpp
    sequence/adapt_struct.cpp
    sequence/adapt_tpl_adt.cpp
    sequence/adapt_tpl_struct.cpp
    sequence/adt_attribute_proxy.cpp
    sequence/define_struct.cpp
    sequence/define_assoc_struct.cpp
    sequence/define_tpl_struct.cpp
    sequence/define_assoc_tpl_struct.cpp

    functional/fused.cpp
    functional/fused_function_object.cpp
    functional/fused_procedure.cpp
    functional/unfused.cpp
    functional/unfused_typed.cpp
    functional/make_fused.cpp
    functional/make_fused_function_object.cpp
    functional/make_fused_procedure.cpp
    functional/make_unfused.cpp
    functional/invoke.cpp
    functional/invoke_function_object.cpp
    functional/invoke_procedure.cpp
    sequence/swap.cpp

  COMPILE_FAIL
    #xxx.cpp
    
  NO_SINGLE_TARGET
  )

# Text for extension features, must be explicitly specified on the command line to be run
# TODO what does "explicit" mean in this context? How can we make them pass?

#boost_add_multiple_tests(RUN
#  algorithm/ext_/for_each_s.cpp
#  algorithm/ext_/find_if_s.cpp
#  sequence/ext_/iterator_range_s.cpp
#  )
