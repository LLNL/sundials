# Install script for directory: /home/maggul/STS/sundials/examples/arkode/C_serial

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/home/maggul/STS/sundials_install")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "RelWithDebInfo")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "1")
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

# Set default install directory permissions.
if(NOT DEFINED CMAKE_OBJDUMP)
  set(CMAKE_OBJDUMP "/usr/bin/objdump")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/maggul/STS/sundials_install/examples/arkode/C_serial/ark_analytic.c;/home/maggul/STS/sundials_install/examples/arkode/C_serial/ark_analytic.out")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  file(INSTALL DESTINATION "/home/maggul/STS/sundials_install/examples/arkode/C_serial" TYPE FILE FILES
    "/home/maggul/STS/sundials/examples/arkode/C_serial/ark_analytic.c"
    "/home/maggul/STS/sundials/examples/arkode/C_serial/ark_analytic.out"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/maggul/STS/sundials_install/examples/arkode/C_serial/lsrk_analytic.c;/home/maggul/STS/sundials_install/examples/arkode/C_serial/lsrk_analytic.out")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  file(INSTALL DESTINATION "/home/maggul/STS/sundials_install/examples/arkode/C_serial" TYPE FILE FILES
    "/home/maggul/STS/sundials/examples/arkode/C_serial/lsrk_analytic.c"
    "/home/maggul/STS/sundials/examples/arkode/C_serial/lsrk_analytic.out"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/maggul/STS/sundials_install/examples/arkode/C_serial/lsrk_analytic_VarJac.c;/home/maggul/STS/sundials_install/examples/arkode/C_serial/lsrk_analytic_VarJac.out")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  file(INSTALL DESTINATION "/home/maggul/STS/sundials_install/examples/arkode/C_serial" TYPE FILE FILES
    "/home/maggul/STS/sundials/examples/arkode/C_serial/lsrk_analytic_VarJac.c"
    "/home/maggul/STS/sundials/examples/arkode/C_serial/lsrk_analytic_VarJac.out"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/maggul/STS/sundials_install/examples/arkode/C_serial/erk_analytic.c;/home/maggul/STS/sundials_install/examples/arkode/C_serial/erk_analytic.out")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  file(INSTALL DESTINATION "/home/maggul/STS/sundials_install/examples/arkode/C_serial" TYPE FILE FILES
    "/home/maggul/STS/sundials/examples/arkode/C_serial/erk_analytic.c"
    "/home/maggul/STS/sundials/examples/arkode/C_serial/erk_analytic.out"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/maggul/STS/sundials_install/examples/arkode/C_serial/ark_analytic_mels.c;/home/maggul/STS/sundials_install/examples/arkode/C_serial/ark_analytic_mels.out")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  file(INSTALL DESTINATION "/home/maggul/STS/sundials_install/examples/arkode/C_serial" TYPE FILE FILES
    "/home/maggul/STS/sundials/examples/arkode/C_serial/ark_analytic_mels.c"
    "/home/maggul/STS/sundials/examples/arkode/C_serial/ark_analytic_mels.out"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/maggul/STS/sundials_install/examples/arkode/C_serial/ark_analytic_nonlin.c;/home/maggul/STS/sundials_install/examples/arkode/C_serial/ark_analytic_nonlin.out")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  file(INSTALL DESTINATION "/home/maggul/STS/sundials_install/examples/arkode/C_serial" TYPE FILE FILES
    "/home/maggul/STS/sundials/examples/arkode/C_serial/ark_analytic_nonlin.c"
    "/home/maggul/STS/sundials/examples/arkode/C_serial/ark_analytic_nonlin.out"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/maggul/STS/sundials_install/examples/arkode/C_serial/ark_brusselator_1D_mri.c;/home/maggul/STS/sundials_install/examples/arkode/C_serial/ark_brusselator_1D_mri.out")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  file(INSTALL DESTINATION "/home/maggul/STS/sundials_install/examples/arkode/C_serial" TYPE FILE FILES
    "/home/maggul/STS/sundials/examples/arkode/C_serial/ark_brusselator_1D_mri.c"
    "/home/maggul/STS/sundials/examples/arkode/C_serial/ark_brusselator_1D_mri.out"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/maggul/STS/sundials_install/examples/arkode/C_serial/ark_brusselator_fp.c;/home/maggul/STS/sundials_install/examples/arkode/C_serial/ark_brusselator_fp.out")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  file(INSTALL DESTINATION "/home/maggul/STS/sundials_install/examples/arkode/C_serial" TYPE FILE FILES
    "/home/maggul/STS/sundials/examples/arkode/C_serial/ark_brusselator_fp.c"
    "/home/maggul/STS/sundials/examples/arkode/C_serial/ark_brusselator_fp.out"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/maggul/STS/sundials_install/examples/arkode/C_serial/ark_brusselator_mri.c;/home/maggul/STS/sundials_install/examples/arkode/C_serial/ark_brusselator_mri.out")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  file(INSTALL DESTINATION "/home/maggul/STS/sundials_install/examples/arkode/C_serial" TYPE FILE FILES
    "/home/maggul/STS/sundials/examples/arkode/C_serial/ark_brusselator_mri.c"
    "/home/maggul/STS/sundials/examples/arkode/C_serial/ark_brusselator_mri.out"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/maggul/STS/sundials_install/examples/arkode/C_serial/ark_brusselator.c;/home/maggul/STS/sundials_install/examples/arkode/C_serial/ark_brusselator.out")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  file(INSTALL DESTINATION "/home/maggul/STS/sundials_install/examples/arkode/C_serial" TYPE FILE FILES
    "/home/maggul/STS/sundials/examples/arkode/C_serial/ark_brusselator.c"
    "/home/maggul/STS/sundials/examples/arkode/C_serial/ark_brusselator.out"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/maggul/STS/sundials_install/examples/arkode/C_serial/ark_brusselator1D_imexmri.c;/home/maggul/STS/sundials_install/examples/arkode/C_serial/ark_brusselator1D_imexmri_0_0.001.out")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  file(INSTALL DESTINATION "/home/maggul/STS/sundials_install/examples/arkode/C_serial" TYPE FILE FILES
    "/home/maggul/STS/sundials/examples/arkode/C_serial/ark_brusselator1D_imexmri.c"
    "/home/maggul/STS/sundials/examples/arkode/C_serial/ark_brusselator1D_imexmri_0_0.001.out"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/maggul/STS/sundials_install/examples/arkode/C_serial/ark_brusselator1D_imexmri.c;/home/maggul/STS/sundials_install/examples/arkode/C_serial/ark_brusselator1D_imexmri_2_0.001.out")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  file(INSTALL DESTINATION "/home/maggul/STS/sundials_install/examples/arkode/C_serial" TYPE FILE FILES
    "/home/maggul/STS/sundials/examples/arkode/C_serial/ark_brusselator1D_imexmri.c"
    "/home/maggul/STS/sundials/examples/arkode/C_serial/ark_brusselator1D_imexmri_2_0.001.out"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/maggul/STS/sundials_install/examples/arkode/C_serial/ark_brusselator1D_imexmri.c;/home/maggul/STS/sundials_install/examples/arkode/C_serial/ark_brusselator1D_imexmri_3_0.001.out")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  file(INSTALL DESTINATION "/home/maggul/STS/sundials_install/examples/arkode/C_serial" TYPE FILE FILES
    "/home/maggul/STS/sundials/examples/arkode/C_serial/ark_brusselator1D_imexmri.c"
    "/home/maggul/STS/sundials/examples/arkode/C_serial/ark_brusselator1D_imexmri_3_0.001.out"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/maggul/STS/sundials_install/examples/arkode/C_serial/ark_brusselator1D_imexmri.c;/home/maggul/STS/sundials_install/examples/arkode/C_serial/ark_brusselator1D_imexmri_4_0.001.out")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  file(INSTALL DESTINATION "/home/maggul/STS/sundials_install/examples/arkode/C_serial" TYPE FILE FILES
    "/home/maggul/STS/sundials/examples/arkode/C_serial/ark_brusselator1D_imexmri.c"
    "/home/maggul/STS/sundials/examples/arkode/C_serial/ark_brusselator1D_imexmri_4_0.001.out"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/maggul/STS/sundials_install/examples/arkode/C_serial/ark_brusselator1D_imexmri.c;/home/maggul/STS/sundials_install/examples/arkode/C_serial/ark_brusselator1D_imexmri_5_0.001.out")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  file(INSTALL DESTINATION "/home/maggul/STS/sundials_install/examples/arkode/C_serial" TYPE FILE FILES
    "/home/maggul/STS/sundials/examples/arkode/C_serial/ark_brusselator1D_imexmri.c"
    "/home/maggul/STS/sundials/examples/arkode/C_serial/ark_brusselator1D_imexmri_5_0.001.out"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/maggul/STS/sundials_install/examples/arkode/C_serial/ark_brusselator1D_imexmri.c;/home/maggul/STS/sundials_install/examples/arkode/C_serial/ark_brusselator1D_imexmri_6_0.001.out")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  file(INSTALL DESTINATION "/home/maggul/STS/sundials_install/examples/arkode/C_serial" TYPE FILE FILES
    "/home/maggul/STS/sundials/examples/arkode/C_serial/ark_brusselator1D_imexmri.c"
    "/home/maggul/STS/sundials/examples/arkode/C_serial/ark_brusselator1D_imexmri_6_0.001.out"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/maggul/STS/sundials_install/examples/arkode/C_serial/ark_brusselator1D_imexmri.c;/home/maggul/STS/sundials_install/examples/arkode/C_serial/ark_brusselator1D_imexmri_7_0.001.out")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  file(INSTALL DESTINATION "/home/maggul/STS/sundials_install/examples/arkode/C_serial" TYPE FILE FILES
    "/home/maggul/STS/sundials/examples/arkode/C_serial/ark_brusselator1D_imexmri.c"
    "/home/maggul/STS/sundials/examples/arkode/C_serial/ark_brusselator1D_imexmri_7_0.001.out"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/maggul/STS/sundials_install/examples/arkode/C_serial/ark_brusselator1D.c;/home/maggul/STS/sundials_install/examples/arkode/C_serial/ark_brusselator1D.out")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  file(INSTALL DESTINATION "/home/maggul/STS/sundials_install/examples/arkode/C_serial" TYPE FILE FILES
    "/home/maggul/STS/sundials/examples/arkode/C_serial/ark_brusselator1D.c"
    "/home/maggul/STS/sundials/examples/arkode/C_serial/ark_brusselator1D.out"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/maggul/STS/sundials_install/examples/arkode/C_serial/ark_conserved_exp_entropy_ark.c;/home/maggul/STS/sundials_install/examples/arkode/C_serial/ark_conserved_exp_entropy_ark_1_0.out")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  file(INSTALL DESTINATION "/home/maggul/STS/sundials_install/examples/arkode/C_serial" TYPE FILE FILES
    "/home/maggul/STS/sundials/examples/arkode/C_serial/ark_conserved_exp_entropy_ark.c"
    "/home/maggul/STS/sundials/examples/arkode/C_serial/ark_conserved_exp_entropy_ark_1_0.out"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/maggul/STS/sundials_install/examples/arkode/C_serial/ark_conserved_exp_entropy_ark.c;/home/maggul/STS/sundials_install/examples/arkode/C_serial/ark_conserved_exp_entropy_ark_1_1.out")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  file(INSTALL DESTINATION "/home/maggul/STS/sundials_install/examples/arkode/C_serial" TYPE FILE FILES
    "/home/maggul/STS/sundials/examples/arkode/C_serial/ark_conserved_exp_entropy_ark.c"
    "/home/maggul/STS/sundials/examples/arkode/C_serial/ark_conserved_exp_entropy_ark_1_1.out"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/maggul/STS/sundials_install/examples/arkode/C_serial/ark_conserved_exp_entropy_erk.c;/home/maggul/STS/sundials_install/examples/arkode/C_serial/ark_conserved_exp_entropy_erk_1.out")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  file(INSTALL DESTINATION "/home/maggul/STS/sundials_install/examples/arkode/C_serial" TYPE FILE FILES
    "/home/maggul/STS/sundials/examples/arkode/C_serial/ark_conserved_exp_entropy_erk.c"
    "/home/maggul/STS/sundials/examples/arkode/C_serial/ark_conserved_exp_entropy_erk_1.out"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/maggul/STS/sundials_install/examples/arkode/C_serial/ark_damped_harmonic_symplectic.c;/home/maggul/STS/sundials_install/examples/arkode/C_serial/ark_damped_harmonic_symplectic.out")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  file(INSTALL DESTINATION "/home/maggul/STS/sundials_install/examples/arkode/C_serial" TYPE FILE FILES
    "/home/maggul/STS/sundials/examples/arkode/C_serial/ark_damped_harmonic_symplectic.c"
    "/home/maggul/STS/sundials/examples/arkode/C_serial/ark_damped_harmonic_symplectic.out"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/maggul/STS/sundials_install/examples/arkode/C_serial/ark_dissipated_exp_entropy.c;/home/maggul/STS/sundials_install/examples/arkode/C_serial/ark_dissipated_exp_entropy_1_0.out")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  file(INSTALL DESTINATION "/home/maggul/STS/sundials_install/examples/arkode/C_serial" TYPE FILE FILES
    "/home/maggul/STS/sundials/examples/arkode/C_serial/ark_dissipated_exp_entropy.c"
    "/home/maggul/STS/sundials/examples/arkode/C_serial/ark_dissipated_exp_entropy_1_0.out"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/maggul/STS/sundials_install/examples/arkode/C_serial/ark_dissipated_exp_entropy.c;/home/maggul/STS/sundials_install/examples/arkode/C_serial/ark_dissipated_exp_entropy_1_1.out")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  file(INSTALL DESTINATION "/home/maggul/STS/sundials_install/examples/arkode/C_serial" TYPE FILE FILES
    "/home/maggul/STS/sundials/examples/arkode/C_serial/ark_dissipated_exp_entropy.c"
    "/home/maggul/STS/sundials/examples/arkode/C_serial/ark_dissipated_exp_entropy_1_1.out"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/maggul/STS/sundials_install/examples/arkode/C_serial/ark_harmonic_symplectic.c;/home/maggul/STS/sundials_install/examples/arkode/C_serial/ark_harmonic_symplectic.out")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  file(INSTALL DESTINATION "/home/maggul/STS/sundials_install/examples/arkode/C_serial" TYPE FILE FILES
    "/home/maggul/STS/sundials/examples/arkode/C_serial/ark_harmonic_symplectic.c"
    "/home/maggul/STS/sundials/examples/arkode/C_serial/ark_harmonic_symplectic.out"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/maggul/STS/sundials_install/examples/arkode/C_serial/ark_heat1D_adapt.c;/home/maggul/STS/sundials_install/examples/arkode/C_serial/ark_heat1D_adapt.out")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  file(INSTALL DESTINATION "/home/maggul/STS/sundials_install/examples/arkode/C_serial" TYPE FILE FILES
    "/home/maggul/STS/sundials/examples/arkode/C_serial/ark_heat1D_adapt.c"
    "/home/maggul/STS/sundials/examples/arkode/C_serial/ark_heat1D_adapt.out"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/maggul/STS/sundials_install/examples/arkode/C_serial/ark_heat1D.c;/home/maggul/STS/sundials_install/examples/arkode/C_serial/ark_heat1D.out")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  file(INSTALL DESTINATION "/home/maggul/STS/sundials_install/examples/arkode/C_serial" TYPE FILE FILES
    "/home/maggul/STS/sundials/examples/arkode/C_serial/ark_heat1D.c"
    "/home/maggul/STS/sundials/examples/arkode/C_serial/ark_heat1D.out"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/maggul/STS/sundials_install/examples/arkode/C_serial/ark_kepler.c;/home/maggul/STS/sundials_install/examples/arkode/C_serial/ark_kepler_--stepper_ERK_--step-mode_adapt.out")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  file(INSTALL DESTINATION "/home/maggul/STS/sundials_install/examples/arkode/C_serial" TYPE FILE FILES
    "/home/maggul/STS/sundials/examples/arkode/C_serial/ark_kepler.c"
    "/home/maggul/STS/sundials/examples/arkode/C_serial/ark_kepler_--stepper_ERK_--step-mode_adapt.out"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/maggul/STS/sundials_install/examples/arkode/C_serial/ark_kepler.c;/home/maggul/STS/sundials_install/examples/arkode/C_serial/ark_kepler_--stepper_ERK_--step-mode_fixed_--count-orbits.out")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  file(INSTALL DESTINATION "/home/maggul/STS/sundials_install/examples/arkode/C_serial" TYPE FILE FILES
    "/home/maggul/STS/sundials/examples/arkode/C_serial/ark_kepler.c"
    "/home/maggul/STS/sundials/examples/arkode/C_serial/ark_kepler_--stepper_ERK_--step-mode_fixed_--count-orbits.out"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/maggul/STS/sundials_install/examples/arkode/C_serial/ark_kepler.c;/home/maggul/STS/sundials_install/examples/arkode/C_serial/ark_kepler_--stepper_SPRK_--step-mode_fixed_--count-orbits_--use-compensated-sums.out")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  file(INSTALL DESTINATION "/home/maggul/STS/sundials_install/examples/arkode/C_serial" TYPE FILE FILES
    "/home/maggul/STS/sundials/examples/arkode/C_serial/ark_kepler.c"
    "/home/maggul/STS/sundials/examples/arkode/C_serial/ark_kepler_--stepper_SPRK_--step-mode_fixed_--count-orbits_--use-compensated-sums.out"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/maggul/STS/sundials_install/examples/arkode/C_serial/ark_kepler.c;/home/maggul/STS/sundials_install/examples/arkode/C_serial/ark_kepler_--stepper_SPRK_--step-mode_fixed_--method_ARKODE_SPRK_EULER_1_1_--tf_50_--check-order_--nout_1.out")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  file(INSTALL DESTINATION "/home/maggul/STS/sundials_install/examples/arkode/C_serial" TYPE FILE FILES
    "/home/maggul/STS/sundials/examples/arkode/C_serial/ark_kepler.c"
    "/home/maggul/STS/sundials/examples/arkode/C_serial/ark_kepler_--stepper_SPRK_--step-mode_fixed_--method_ARKODE_SPRK_EULER_1_1_--tf_50_--check-order_--nout_1.out"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/maggul/STS/sundials_install/examples/arkode/C_serial/ark_kepler.c;/home/maggul/STS/sundials_install/examples/arkode/C_serial/ark_kepler_--stepper_SPRK_--step-mode_fixed_--method_ARKODE_SPRK_LEAPFROG_2_2_--tf_50_--check-order_--nout_1.out")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  file(INSTALL DESTINATION "/home/maggul/STS/sundials_install/examples/arkode/C_serial" TYPE FILE FILES
    "/home/maggul/STS/sundials/examples/arkode/C_serial/ark_kepler.c"
    "/home/maggul/STS/sundials/examples/arkode/C_serial/ark_kepler_--stepper_SPRK_--step-mode_fixed_--method_ARKODE_SPRK_LEAPFROG_2_2_--tf_50_--check-order_--nout_1.out"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/maggul/STS/sundials_install/examples/arkode/C_serial/ark_kepler.c;/home/maggul/STS/sundials_install/examples/arkode/C_serial/ark_kepler_--stepper_SPRK_--step-mode_fixed_--method_ARKODE_SPRK_MCLACHLAN_2_2_--tf_50_--check-order_--nout_1.out")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  file(INSTALL DESTINATION "/home/maggul/STS/sundials_install/examples/arkode/C_serial" TYPE FILE FILES
    "/home/maggul/STS/sundials/examples/arkode/C_serial/ark_kepler.c"
    "/home/maggul/STS/sundials/examples/arkode/C_serial/ark_kepler_--stepper_SPRK_--step-mode_fixed_--method_ARKODE_SPRK_MCLACHLAN_2_2_--tf_50_--check-order_--nout_1.out"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/maggul/STS/sundials_install/examples/arkode/C_serial/ark_kepler.c;/home/maggul/STS/sundials_install/examples/arkode/C_serial/ark_kepler_--stepper_SPRK_--step-mode_fixed_--method_ARKODE_SPRK_MCLACHLAN_3_3_--tf_50_--check-order_--nout_1.out")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  file(INSTALL DESTINATION "/home/maggul/STS/sundials_install/examples/arkode/C_serial" TYPE FILE FILES
    "/home/maggul/STS/sundials/examples/arkode/C_serial/ark_kepler.c"
    "/home/maggul/STS/sundials/examples/arkode/C_serial/ark_kepler_--stepper_SPRK_--step-mode_fixed_--method_ARKODE_SPRK_MCLACHLAN_3_3_--tf_50_--check-order_--nout_1.out"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/maggul/STS/sundials_install/examples/arkode/C_serial/ark_kepler.c;/home/maggul/STS/sundials_install/examples/arkode/C_serial/ark_kepler_--stepper_SPRK_--step-mode_fixed_--method_ARKODE_SPRK_MCLACHLAN_4_4_--tf_50_--check-order_--nout_1.out")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  file(INSTALL DESTINATION "/home/maggul/STS/sundials_install/examples/arkode/C_serial" TYPE FILE FILES
    "/home/maggul/STS/sundials/examples/arkode/C_serial/ark_kepler.c"
    "/home/maggul/STS/sundials/examples/arkode/C_serial/ark_kepler_--stepper_SPRK_--step-mode_fixed_--method_ARKODE_SPRK_MCLACHLAN_4_4_--tf_50_--check-order_--nout_1.out"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/maggul/STS/sundials_install/examples/arkode/C_serial/ark_kepler.c;/home/maggul/STS/sundials_install/examples/arkode/C_serial/ark_kepler_--stepper_SPRK_--step-mode_fixed_--method_ARKODE_SPRK_MCLACHLAN_5_6_--tf_50_--check-order_--nout_1.out")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  file(INSTALL DESTINATION "/home/maggul/STS/sundials_install/examples/arkode/C_serial" TYPE FILE FILES
    "/home/maggul/STS/sundials/examples/arkode/C_serial/ark_kepler.c"
    "/home/maggul/STS/sundials/examples/arkode/C_serial/ark_kepler_--stepper_SPRK_--step-mode_fixed_--method_ARKODE_SPRK_MCLACHLAN_5_6_--tf_50_--check-order_--nout_1.out"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/maggul/STS/sundials_install/examples/arkode/C_serial/ark_kepler.c;/home/maggul/STS/sundials_install/examples/arkode/C_serial/ark_kepler_--stepper_SPRK_--step-mode_fixed_--method_ARKODE_SPRK_PSEUDO_LEAPFROG_2_2_--tf_50_--check-order_--nout_1.out")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  file(INSTALL DESTINATION "/home/maggul/STS/sundials_install/examples/arkode/C_serial" TYPE FILE FILES
    "/home/maggul/STS/sundials/examples/arkode/C_serial/ark_kepler.c"
    "/home/maggul/STS/sundials/examples/arkode/C_serial/ark_kepler_--stepper_SPRK_--step-mode_fixed_--method_ARKODE_SPRK_PSEUDO_LEAPFROG_2_2_--tf_50_--check-order_--nout_1.out"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/maggul/STS/sundials_install/examples/arkode/C_serial/ark_kepler.c;/home/maggul/STS/sundials_install/examples/arkode/C_serial/ark_kepler_--stepper_SPRK_--step-mode_fixed_--method_ARKODE_SPRK_RUTH_3_3_--tf_50_--check-order_--nout_1.out")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  file(INSTALL DESTINATION "/home/maggul/STS/sundials_install/examples/arkode/C_serial" TYPE FILE FILES
    "/home/maggul/STS/sundials/examples/arkode/C_serial/ark_kepler.c"
    "/home/maggul/STS/sundials/examples/arkode/C_serial/ark_kepler_--stepper_SPRK_--step-mode_fixed_--method_ARKODE_SPRK_RUTH_3_3_--tf_50_--check-order_--nout_1.out"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/maggul/STS/sundials_install/examples/arkode/C_serial/ark_kepler.c;/home/maggul/STS/sundials_install/examples/arkode/C_serial/ark_kepler_--stepper_SPRK_--step-mode_fixed_--method_ARKODE_SPRK_YOSHIDA_6_8_--tf_50_--check-order_--nout_1.out")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  file(INSTALL DESTINATION "/home/maggul/STS/sundials_install/examples/arkode/C_serial" TYPE FILE FILES
    "/home/maggul/STS/sundials/examples/arkode/C_serial/ark_kepler.c"
    "/home/maggul/STS/sundials/examples/arkode/C_serial/ark_kepler_--stepper_SPRK_--step-mode_fixed_--method_ARKODE_SPRK_YOSHIDA_6_8_--tf_50_--check-order_--nout_1.out"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/maggul/STS/sundials_install/examples/arkode/C_serial/ark_kepler.c;/home/maggul/STS/sundials_install/examples/arkode/C_serial/ark_kepler.out")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  file(INSTALL DESTINATION "/home/maggul/STS/sundials_install/examples/arkode/C_serial" TYPE FILE FILES
    "/home/maggul/STS/sundials/examples/arkode/C_serial/ark_kepler.c"
    "/home/maggul/STS/sundials/examples/arkode/C_serial/ark_kepler.out"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/maggul/STS/sundials_install/examples/arkode/C_serial/ark_kpr_mri.c;/home/maggul/STS/sundials_install/examples/arkode/C_serial/ark_kpr_mri.out")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  file(INSTALL DESTINATION "/home/maggul/STS/sundials_install/examples/arkode/C_serial" TYPE FILE FILES
    "/home/maggul/STS/sundials/examples/arkode/C_serial/ark_kpr_mri.c"
    "/home/maggul/STS/sundials/examples/arkode/C_serial/ark_kpr_mri.out"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/maggul/STS/sundials_install/examples/arkode/C_serial/ark_kpr_mri.c;/home/maggul/STS/sundials_install/examples/arkode/C_serial/ark_kpr_mri_0_0.002.out")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  file(INSTALL DESTINATION "/home/maggul/STS/sundials_install/examples/arkode/C_serial" TYPE FILE FILES
    "/home/maggul/STS/sundials/examples/arkode/C_serial/ark_kpr_mri.c"
    "/home/maggul/STS/sundials/examples/arkode/C_serial/ark_kpr_mri_0_0.002.out"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/maggul/STS/sundials_install/examples/arkode/C_serial/ark_kpr_mri.c;/home/maggul/STS/sundials_install/examples/arkode/C_serial/ark_kpr_mri_1_0.002.out")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  file(INSTALL DESTINATION "/home/maggul/STS/sundials_install/examples/arkode/C_serial" TYPE FILE FILES
    "/home/maggul/STS/sundials/examples/arkode/C_serial/ark_kpr_mri.c"
    "/home/maggul/STS/sundials/examples/arkode/C_serial/ark_kpr_mri_1_0.002.out"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/maggul/STS/sundials_install/examples/arkode/C_serial/ark_kpr_mri.c;/home/maggul/STS/sundials_install/examples/arkode/C_serial/ark_kpr_mri_2_0.005.out")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  file(INSTALL DESTINATION "/home/maggul/STS/sundials_install/examples/arkode/C_serial" TYPE FILE FILES
    "/home/maggul/STS/sundials/examples/arkode/C_serial/ark_kpr_mri.c"
    "/home/maggul/STS/sundials/examples/arkode/C_serial/ark_kpr_mri_2_0.005.out"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/maggul/STS/sundials_install/examples/arkode/C_serial/ark_kpr_mri.c;/home/maggul/STS/sundials_install/examples/arkode/C_serial/ark_kpr_mri_3_0.01.out")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  file(INSTALL DESTINATION "/home/maggul/STS/sundials_install/examples/arkode/C_serial" TYPE FILE FILES
    "/home/maggul/STS/sundials/examples/arkode/C_serial/ark_kpr_mri.c"
    "/home/maggul/STS/sundials/examples/arkode/C_serial/ark_kpr_mri_3_0.01.out"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/maggul/STS/sundials_install/examples/arkode/C_serial/ark_kpr_mri.c;/home/maggul/STS/sundials_install/examples/arkode/C_serial/ark_kpr_mri_4_0.002.out")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  file(INSTALL DESTINATION "/home/maggul/STS/sundials_install/examples/arkode/C_serial" TYPE FILE FILES
    "/home/maggul/STS/sundials/examples/arkode/C_serial/ark_kpr_mri.c"
    "/home/maggul/STS/sundials/examples/arkode/C_serial/ark_kpr_mri_4_0.002.out"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/maggul/STS/sundials_install/examples/arkode/C_serial/ark_kpr_mri.c;/home/maggul/STS/sundials_install/examples/arkode/C_serial/ark_kpr_mri_5_0.002.out")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  file(INSTALL DESTINATION "/home/maggul/STS/sundials_install/examples/arkode/C_serial" TYPE FILE FILES
    "/home/maggul/STS/sundials/examples/arkode/C_serial/ark_kpr_mri.c"
    "/home/maggul/STS/sundials/examples/arkode/C_serial/ark_kpr_mri_5_0.002.out"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/maggul/STS/sundials_install/examples/arkode/C_serial/ark_kpr_mri.c;/home/maggul/STS/sundials_install/examples/arkode/C_serial/ark_kpr_mri_6_0.005.out")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  file(INSTALL DESTINATION "/home/maggul/STS/sundials_install/examples/arkode/C_serial" TYPE FILE FILES
    "/home/maggul/STS/sundials/examples/arkode/C_serial/ark_kpr_mri.c"
    "/home/maggul/STS/sundials/examples/arkode/C_serial/ark_kpr_mri_6_0.005.out"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/maggul/STS/sundials_install/examples/arkode/C_serial/ark_kpr_mri.c;/home/maggul/STS/sundials_install/examples/arkode/C_serial/ark_kpr_mri_7_0.001_-100_100_0.5_1.out")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  file(INSTALL DESTINATION "/home/maggul/STS/sundials_install/examples/arkode/C_serial" TYPE FILE FILES
    "/home/maggul/STS/sundials/examples/arkode/C_serial/ark_kpr_mri.c"
    "/home/maggul/STS/sundials/examples/arkode/C_serial/ark_kpr_mri_7_0.001_-100_100_0.5_1.out"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/maggul/STS/sundials_install/examples/arkode/C_serial/ark_kpr_mri.c;/home/maggul/STS/sundials_install/examples/arkode/C_serial/ark_kpr_mri_7_0.001.out")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  file(INSTALL DESTINATION "/home/maggul/STS/sundials_install/examples/arkode/C_serial" TYPE FILE FILES
    "/home/maggul/STS/sundials/examples/arkode/C_serial/ark_kpr_mri.c"
    "/home/maggul/STS/sundials/examples/arkode/C_serial/ark_kpr_mri_7_0.001.out"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/maggul/STS/sundials_install/examples/arkode/C_serial/ark_kpr_mri.c;/home/maggul/STS/sundials_install/examples/arkode/C_serial/ark_kpr_mri_8_0.001_-100_100_0.5_1.out")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  file(INSTALL DESTINATION "/home/maggul/STS/sundials_install/examples/arkode/C_serial" TYPE FILE FILES
    "/home/maggul/STS/sundials/examples/arkode/C_serial/ark_kpr_mri.c"
    "/home/maggul/STS/sundials/examples/arkode/C_serial/ark_kpr_mri_8_0.001_-100_100_0.5_1.out"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/maggul/STS/sundials_install/examples/arkode/C_serial/ark_kpr_mri.c;/home/maggul/STS/sundials_install/examples/arkode/C_serial/ark_kpr_mri_8_0.001.out")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  file(INSTALL DESTINATION "/home/maggul/STS/sundials_install/examples/arkode/C_serial" TYPE FILE FILES
    "/home/maggul/STS/sundials/examples/arkode/C_serial/ark_kpr_mri.c"
    "/home/maggul/STS/sundials/examples/arkode/C_serial/ark_kpr_mri_8_0.001.out"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/maggul/STS/sundials_install/examples/arkode/C_serial/ark_kpr_mri.c;/home/maggul/STS/sundials_install/examples/arkode/C_serial/ark_kpr_mri_9_0.001_-100_100_0.5_1.out")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  file(INSTALL DESTINATION "/home/maggul/STS/sundials_install/examples/arkode/C_serial" TYPE FILE FILES
    "/home/maggul/STS/sundials/examples/arkode/C_serial/ark_kpr_mri.c"
    "/home/maggul/STS/sundials/examples/arkode/C_serial/ark_kpr_mri_9_0.001_-100_100_0.5_1.out"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/maggul/STS/sundials_install/examples/arkode/C_serial/ark_kpr_mri.c;/home/maggul/STS/sundials_install/examples/arkode/C_serial/ark_kpr_mri_9_0.001.out")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  file(INSTALL DESTINATION "/home/maggul/STS/sundials_install/examples/arkode/C_serial" TYPE FILE FILES
    "/home/maggul/STS/sundials/examples/arkode/C_serial/ark_kpr_mri.c"
    "/home/maggul/STS/sundials/examples/arkode/C_serial/ark_kpr_mri_9_0.001.out"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/maggul/STS/sundials_install/examples/arkode/C_serial/ark_KrylovDemo_prec.c;/home/maggul/STS/sundials_install/examples/arkode/C_serial/ark_KrylovDemo_prec.out")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  file(INSTALL DESTINATION "/home/maggul/STS/sundials_install/examples/arkode/C_serial" TYPE FILE FILES
    "/home/maggul/STS/sundials/examples/arkode/C_serial/ark_KrylovDemo_prec.c"
    "/home/maggul/STS/sundials/examples/arkode/C_serial/ark_KrylovDemo_prec.out"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/maggul/STS/sundials_install/examples/arkode/C_serial/ark_KrylovDemo_prec.c;/home/maggul/STS/sundials_install/examples/arkode/C_serial/ark_KrylovDemo_prec_1.out")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  file(INSTALL DESTINATION "/home/maggul/STS/sundials_install/examples/arkode/C_serial" TYPE FILE FILES
    "/home/maggul/STS/sundials/examples/arkode/C_serial/ark_KrylovDemo_prec.c"
    "/home/maggul/STS/sundials/examples/arkode/C_serial/ark_KrylovDemo_prec_1.out"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/maggul/STS/sundials_install/examples/arkode/C_serial/ark_KrylovDemo_prec.c;/home/maggul/STS/sundials_install/examples/arkode/C_serial/ark_KrylovDemo_prec_2.out")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  file(INSTALL DESTINATION "/home/maggul/STS/sundials_install/examples/arkode/C_serial" TYPE FILE FILES
    "/home/maggul/STS/sundials/examples/arkode/C_serial/ark_KrylovDemo_prec.c"
    "/home/maggul/STS/sundials/examples/arkode/C_serial/ark_KrylovDemo_prec_2.out"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/maggul/STS/sundials_install/examples/arkode/C_serial/ark_onewaycouple_mri.c;/home/maggul/STS/sundials_install/examples/arkode/C_serial/ark_onewaycouple_mri.out")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  file(INSTALL DESTINATION "/home/maggul/STS/sundials_install/examples/arkode/C_serial" TYPE FILE FILES
    "/home/maggul/STS/sundials/examples/arkode/C_serial/ark_onewaycouple_mri.c"
    "/home/maggul/STS/sundials/examples/arkode/C_serial/ark_onewaycouple_mri.out"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/maggul/STS/sundials_install/examples/arkode/C_serial/ark_reaction_diffusion_mri.c;/home/maggul/STS/sundials_install/examples/arkode/C_serial/ark_reaction_diffusion_mri.out")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  file(INSTALL DESTINATION "/home/maggul/STS/sundials_install/examples/arkode/C_serial" TYPE FILE FILES
    "/home/maggul/STS/sundials/examples/arkode/C_serial/ark_reaction_diffusion_mri.c"
    "/home/maggul/STS/sundials/examples/arkode/C_serial/ark_reaction_diffusion_mri.out"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/maggul/STS/sundials_install/examples/arkode/C_serial/ark_robertson_constraints.c;/home/maggul/STS/sundials_install/examples/arkode/C_serial/ark_robertson_constraints.out")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  file(INSTALL DESTINATION "/home/maggul/STS/sundials_install/examples/arkode/C_serial" TYPE FILE FILES
    "/home/maggul/STS/sundials/examples/arkode/C_serial/ark_robertson_constraints.c"
    "/home/maggul/STS/sundials/examples/arkode/C_serial/ark_robertson_constraints.out"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/maggul/STS/sundials_install/examples/arkode/C_serial/ark_robertson_root.c;/home/maggul/STS/sundials_install/examples/arkode/C_serial/ark_robertson_root.out")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  file(INSTALL DESTINATION "/home/maggul/STS/sundials_install/examples/arkode/C_serial" TYPE FILE FILES
    "/home/maggul/STS/sundials/examples/arkode/C_serial/ark_robertson_root.c"
    "/home/maggul/STS/sundials/examples/arkode/C_serial/ark_robertson_root.out"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/maggul/STS/sundials_install/examples/arkode/C_serial/ark_robertson.c;/home/maggul/STS/sundials_install/examples/arkode/C_serial/ark_robertson.out")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  file(INSTALL DESTINATION "/home/maggul/STS/sundials_install/examples/arkode/C_serial" TYPE FILE FILES
    "/home/maggul/STS/sundials/examples/arkode/C_serial/ark_robertson.c"
    "/home/maggul/STS/sundials/examples/arkode/C_serial/ark_robertson.out"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/maggul/STS/sundials_install/examples/arkode/C_serial/ark_twowaycouple_mri.c;/home/maggul/STS/sundials_install/examples/arkode/C_serial/ark_twowaycouple_mri.out")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  file(INSTALL DESTINATION "/home/maggul/STS/sundials_install/examples/arkode/C_serial" TYPE FILE FILES
    "/home/maggul/STS/sundials/examples/arkode/C_serial/ark_twowaycouple_mri.c"
    "/home/maggul/STS/sundials/examples/arkode/C_serial/ark_twowaycouple_mri.out"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/maggul/STS/sundials_install/examples/arkode/C_serial/README")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  file(INSTALL DESTINATION "/home/maggul/STS/sundials_install/examples/arkode/C_serial" TYPE FILE FILES "/home/maggul/STS/sundials/examples/arkode/C_serial/README")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/maggul/STS/sundials_install/examples/arkode/C_serial/ark_analytic_nonlin_stats.csv")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  file(INSTALL DESTINATION "/home/maggul/STS/sundials_install/examples/arkode/C_serial" TYPE FILE FILES "/home/maggul/STS/sundials/examples/arkode/C_serial/ark_analytic_nonlin_stats.csv")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/maggul/STS/sundials_install/examples/arkode/C_serial/ark_damped_harmonic_symplectic.h")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  file(INSTALL DESTINATION "/home/maggul/STS/sundials_install/examples/arkode/C_serial" TYPE FILE FILES "/home/maggul/STS/sundials/examples/arkode/C_serial/ark_damped_harmonic_symplectic.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/maggul/STS/sundials_install/examples/arkode/C_serial/ark_harmonic_symplectic.h")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  file(INSTALL DESTINATION "/home/maggul/STS/sundials_install/examples/arkode/C_serial" TYPE FILE FILES "/home/maggul/STS/sundials/examples/arkode/C_serial/ark_harmonic_symplectic.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/maggul/STS/sundials_install/examples/arkode/C_serial/ark_kepler_plot.py")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  file(INSTALL DESTINATION "/home/maggul/STS/sundials_install/examples/arkode/C_serial" TYPE FILE FILES "/home/maggul/STS/sundials/examples/arkode/C_serial/ark_kepler_plot.py")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/maggul/STS/sundials_install/examples/arkode/C_serial/ark_kepler.h")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  file(INSTALL DESTINATION "/home/maggul/STS/sundials_install/examples/arkode/C_serial" TYPE FILE FILES "/home/maggul/STS/sundials/examples/arkode/C_serial/ark_kepler.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/maggul/STS/sundials_install/examples/arkode/C_serial/ark_reaction_diffusion_mri_fast_stats.csv")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  file(INSTALL DESTINATION "/home/maggul/STS/sundials_install/examples/arkode/C_serial" TYPE FILE FILES "/home/maggul/STS/sundials/examples/arkode/C_serial/ark_reaction_diffusion_mri_fast_stats.csv")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/maggul/STS/sundials_install/examples/arkode/C_serial/ark_reaction_diffusion_mri_slow_stats.csv")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  file(INSTALL DESTINATION "/home/maggul/STS/sundials_install/examples/arkode/C_serial" TYPE FILE FILES "/home/maggul/STS/sundials/examples/arkode/C_serial/ark_reaction_diffusion_mri_slow_stats.csv")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/maggul/STS/sundials_install/examples/arkode/C_serial/ark_robertson_stats.csv")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  file(INSTALL DESTINATION "/home/maggul/STS/sundials_install/examples/arkode/C_serial" TYPE FILE FILES "/home/maggul/STS/sundials/examples/arkode/C_serial/ark_robertson_stats.csv")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/maggul/STS/sundials_install/examples/arkode/C_serial/plot_brusselator1D_FEM.py")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  file(INSTALL DESTINATION "/home/maggul/STS/sundials_install/examples/arkode/C_serial" TYPE FILE FILES "/home/maggul/STS/sundials/examples/arkode/C_serial/plot_brusselator1D_FEM.py")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/maggul/STS/sundials_install/examples/arkode/C_serial/plot_brusselator1D.py")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  file(INSTALL DESTINATION "/home/maggul/STS/sundials_install/examples/arkode/C_serial" TYPE FILE FILES "/home/maggul/STS/sundials/examples/arkode/C_serial/plot_brusselator1D.py")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/maggul/STS/sundials_install/examples/arkode/C_serial/plot_heat1D_adapt.py")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  file(INSTALL DESTINATION "/home/maggul/STS/sundials_install/examples/arkode/C_serial" TYPE FILE FILES "/home/maggul/STS/sundials/examples/arkode/C_serial/plot_heat1D_adapt.py")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/maggul/STS/sundials_install/examples/arkode/C_serial/plot_heat1D.py")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  file(INSTALL DESTINATION "/home/maggul/STS/sundials_install/examples/arkode/C_serial" TYPE FILE FILES "/home/maggul/STS/sundials/examples/arkode/C_serial/plot_heat1D.py")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/maggul/STS/sundials_install/examples/arkode/C_serial/plot_sol_log.py")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  file(INSTALL DESTINATION "/home/maggul/STS/sundials_install/examples/arkode/C_serial" TYPE FILE FILES "/home/maggul/STS/sundials/examples/arkode/C_serial/plot_sol_log.py")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/maggul/STS/sundials_install/examples/arkode/C_serial/plot_sol.py")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  file(INSTALL DESTINATION "/home/maggul/STS/sundials_install/examples/arkode/C_serial" TYPE FILE FILES "/home/maggul/STS/sundials/examples/arkode/C_serial/plot_sol.py")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/maggul/STS/sundials_install/examples/arkode/C_serial/CMakeLists.txt")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  file(INSTALL DESTINATION "/home/maggul/STS/sundials_install/examples/arkode/C_serial" TYPE FILE FILES "/home/maggul/STS/sundials/build/examples/arkode/C_serial/CMakeLists.txt")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/maggul/STS/sundials_install/examples/arkode/C_serial/Makefile")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  file(INSTALL DESTINATION "/home/maggul/STS/sundials_install/examples/arkode/C_serial" TYPE FILE RENAME "Makefile" FILES "/home/maggul/STS/sundials/build/examples/arkode/C_serial/Makefile_ex")
endif()

