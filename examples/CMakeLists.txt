
# HACK: So that C++ run-time is linked in if not building TxSSA as shared.
if(NOT BUILD_SHARED_LIBS)
    set_source_files_properties(ssa_real_hermitian.c PROPERTIES LANGUAGE CXX)
endif()

add_executable(ssa_real_hermitian ssa_real_hermitian.c)
add_executable(ssa_file_io_example ssa_file_io_example.cpp)

if(BUILD_SHARED_LIBS)
  target_link_libraries(ssa_real_hermitian TxSSA)
  target_link_libraries(ssa_file_io_example TxSSA)
  if(WIN32)
    add_definitions(-DTXSSA_DLL_IMPORTS)
  endif()
else()
  target_link_libraries(ssa_real_hermitian TxSSA ${BLAS_LAPACK_LIB_PATHS})
  target_link_libraries(ssa_file_io_example TxSSA ${BLAS_LAPACK_LIB_PATHS})
endif()

install(FILES ssa_file_io_example.cpp ssa_real_hermitian.c DESTINATION share/examples/)
install(TARGETS ssa_file_io_example ssa_real_hermitian RUNTIME DESTINATION bin)

