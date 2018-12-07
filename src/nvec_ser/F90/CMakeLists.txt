# ---------------------------------------------------------------
# Programmer:  Cody J. Balos @ LLNL
# ---------------------------------------------------------------
# LLNS Copyright Start
# Copyright (c) 2017, Lawrence Livermore National Security
# This work was performed under the auspices of the U.S. Department
# of Energy by Lawrence Livermore National Laboratory in part under
# Contract W-7405-Eng-48 and in part under Contract DE-AC52-07NA27344.
# Produced at the Lawrence Livermore National Laboratory.
# All rights reserved.
# For details, see the LICENSE file.
# LLNS Copyright End
# ---------------------------------------------------------------
# CMakeLists.txt file for the F2003 serial NVECTOR object library

set(nvecserial_SOURCES fnvector_serial.f90)

if(BUILD_STATIC_LIBS)
  sundials_add_f2003_interface_library(sundials_fnvecserial_mod_static STATIC
    ${nvecserial_SOURCES}
  )
  set_target_properties(sundials_fnvecserial_mod_static
    PROPERTIES OUTPUT_NAME sundials_fnvecserial_mod CLEAN_DIRECT_OUTPUT 1)
  install(TARGETS sundials_fnvecserial_mod_static DESTINATION ${CMAKE_INSTALL_LIBDIR})
endif(BUILD_STATIC_LIBS)

if(BUILD_SHARED_LIBS)
  sundials_add_f2003_interface_library(sundials_fnvecserial_mod_shared SHARED
    ${nvecserial_SOURCES}
  )
  set_target_properties(sundials_fnvecserial_mod_shared
    PROPERTIES OUTPUT_NAME sundials_fnvecserial_mod CLEAN_DIRECT_OUTPUT 1)
  set_target_properties(sundials_fnvecserial_mod_shared
    PROPERTIES VERSION ${nveclib_VERSION} SOVERSION ${nveclib_SOVERSION})
  install(TARGETS sundials_fnvecserial_mod_shared DESTINATION ${CMAKE_INSTALL_LIBDIR})
endif(BUILD_SHARED_LIBS)

message(STATUS "Added NVECTOR_SERIAL F2003 interface")
