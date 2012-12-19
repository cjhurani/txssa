cmake_minimum_required(VERSION 2.6)

project(TxSSA)

if(WIN32)
    set(GREP_PROGRAM "findstr")
else()
    set(GREP_PROGRAM "grep")
endif()

execute_process(COMMAND ${GREP_PROGRAM} " ssa_intf_version" txssa.h
                WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/include
                OUTPUT_VARIABLE intf_version 
                OUTPUT_STRIP_TRAILING_WHITESPACE)

string(REGEX REPLACE ";" "" intf_version ${intf_version})
string(REGEX REPLACE ".* " "" intf_version ${intf_version})

execute_process(COMMAND ${GREP_PROGRAM} " ssa_impl_version" txssa.cpp
                WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/src/sparse_spectral_approximation 
                OUTPUT_VARIABLE impl_version 
                OUTPUT_STRIP_TRAILING_WHITESPACE)

string(REGEX REPLACE ";" "" impl_version ${impl_version})
string(REGEX REPLACE ".* " "" impl_version ${impl_version})

set(TXSSA_VERSION ${intf_version}.${impl_version})

if(NOT BLAS_LAPACK_LIB_PATHS)
    # Find BLAS and LAPACK using CMake

    message(STATUS "TxSSA: Will try to find BLAS and LAPACK libraries.")
    if(NOT WIN32)
        enable_language(Fortran)
    endif()

    include(${CMAKE_ROOT}/Modules/FindBLAS.cmake)

    if(NOT BLAS_FOUND)
        message(FATAL_ERROR "TxSSA: Could not find BLAS.")
    endif()

    include(${CMAKE_ROOT}/Modules/FindLAPACK.cmake)

    if(NOT LAPACK_FOUND)
        message(FATAL_ERROR "TxSSA: Could not find LAPACK.")
    endif()

    message(STATUS "TxSSA: Found BLAS: ${BLAS_LIBRARIES}")
    message(STATUS "TxSSA: Found LAPACK: ${LAPACK_LIBRARIES}")
    set(BLAS_LAPACK_LIB_PATHS "${BLAS_LIBRARIES}" "${LAPACK_LIBRARIES}")
endif()

set(CMAKE_REQUIRED_LIBRARIES ${CMAKE_REQUIRED_LIBRARIES} ${BLAS_LAPACK_LIB_PATHS})
list(REMOVE_DUPLICATES CMAKE_REQUIRED_LIBRARIES)

include(${CMAKE_ROOT}/Modules/CheckFunctionExists.cmake)

# Determine symbol convention
message(STATUS "TxSSA: Trying to find FORTRAN_SYMBOL_SCHEME using libraries ${CMAKE_REQUIRED_LIBRARIES}")

set(FORTRAN_SYMBOL_SCHEME 0)

check_function_exists(dgeqp3_ blas_lapack_name_found)
if(blas_lapack_name_found)
    set(FORTRAN_SYMBOL_SCHEME 1)
else()
    check_function_exists(dgeqp3 blas_lapack_name_found)
    if(blas_lapack_name_found)
        set(FORTRAN_SYMBOL_SCHEME 2)
    else()
        check_function_exists(DGEQP3_ blas_lapack_name_found)
        if(blas_lapack_name_found)
            set(FORTRAN_SYMBOL_SCHEME 3)
        else()
            check_function_exists(DGEQP3 blas_lapack_name_found)
            if(blas_lapack_name_found)
                set(FORTRAN_SYMBOL_SCHEME 4)
            else()
                message(FATAL_ERROR "Could not link with BLAS, LAPACK, and any extra libraries given when looking for symbol corresponding to 'dgeqp3'.")
            endif()
        endif()
    endif()
endif()

add_definitions(-DFORTRAN_SYMBOL_SCHEME=${FORTRAN_SYMBOL_SCHEME})

if(WIN32)
    add_definitions(-D_SCL_SECURE_NO_WARNINGS)
endif()

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src)

add_subdirectory(src)
add_subdirectory(test)
add_subdirectory(examples)

install(DIRECTORY include matlab DESTINATION . PATTERN ".svn" EXCLUDE)

# CPack

if(WIN32)
    set(CPACK_GENERATOR ZIP)
else()
    set(CPACK_GENERATOR TGZ)
endif()

set(CPACK_PACKAGE_NAME TxSSA)
set(CPACK_PACKAGE_FILE_NAME "${CPACK_PACKAGE_NAME}_${TXSSA_VERSION}_${CMAKE_SYSTEM_NAME}" )

include(CPack)

# End CPack

