# ---------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2002-2023, Lawrence Livermore National Security
# and Southern Methodist University.
# All rights reserved.
#
# See the top-level LICENSE and NOTICE files for details.
#
# SPDX-License-Identifier: BSD-3-Clause
# SUNDIALS Copyright End
# ---------------------------------------------------------------

add_subdirectory(sys)

if(ENABLE_CUDA)
  add_subdirectory(cuda)
endif()

if(ENABLE_HIP)
  add_subdirectory(hip)
endif()

if(ENABLE_SYCL)
  add_subdirectory(sycl)
endif()

