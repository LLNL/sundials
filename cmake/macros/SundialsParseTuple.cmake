# ------------------------------------------------------------------------------
# Programmer(s): Shahbaj Sohal @ LLNL
# ------------------------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2002-2023, Lawrence Livermore National Security
# and Southern Methodist University.
# All rights reserved.
#
# See the top-level LICENSE and NOTICE files for details.
#
# SPDX-License-Identifier: BSD-3-Clause
# SUNDIALS Copyright End
# ------------------------------------------------------------------------------
#
# sundials_parse_tuple(<example_tuple>)
#
# CMake macro to parse the lists of example tuples. Keyword input arguments can be
# added after <examples_tuple> to put each tuple into its own list (see oneValueArgs below).

macro(sundials_parse_tuple example_tuple)
    set(options )
    set(oneValueArgs FILENAME_VAR PARAMS_VAR NODES_VAR NPROCS_VAR TYPE_VAR FLT_PRECISION_VAR INT_PRECISION_VAR)
    set(multiValueArgs)
    
    cmake_parse_arguments(sundials_parse_tuple
        "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})
    if(sundials_parse_tuple_FILENAME_VAR)
        list(GET example_tuple 0 ${sundials_parse_tuple_FILENAME_VAR})
    endif()
    if(sundials_parse_tuple_PARAMS_VAR)
        list(GET example_tuple 1 ${sundials_parse_tuple_PARAMS_VAR})
    endif()
    if(sundials_parse_tuple_NODES_VAR)
        list(GET example_tuple 2 ${sundials_parse_tuple_NODES_VAR})
    endif()
    if(sundials_parse_tuple_NPROCS_VAR)
        list(GET example_tuple 3 ${sundials_parse_tuple_NPROCS_VAR})
    endif()
    if(sundials_parse_tuple_TYPE_VAR)
        if(HAS_MPI)
            list(GET example_tuple 4 ${sundials_parse_tuple_TYPE_VAR})
        else()
            list(GET example_tuple 2 ${sundials_parse_tuple_TYPE_VAR})
        endif()
    endif()
    if(sundials_parse_tuple_FLT_PRECISION_VAR)
        if(HAS_OPENMP)
            list(GET example_tuple 3 ${sundials_parse_tuple_FLT_PRECISION_VAR})
        else()
            list(GET example_tuple 5 ${sundials_parse_tuple_FLT_PRECISION_VAR})
        endif()
    endif()
    if(sundials_parse_tuple_INT_PRECISION_VAR)
        list(GET example_tuple 4 ${sundials_parse_tuple_INT_PRECISION_VAR})
    endif()
endmacro()