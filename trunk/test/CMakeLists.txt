# This example requires static library since it uses internal functions
# and they are not exposed by shared library.

if(NOT BUILD_SHARED_LIBS)
    add_executable(test_p_norm_sparsity_vector test_p_norm_sparsity_vector.cpp)
    target_link_libraries(test_p_norm_sparsity_vector TxSSA)
endif()
