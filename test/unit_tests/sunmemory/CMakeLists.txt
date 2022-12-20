# ---------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2002-2022, Lawrence Livermore National Security
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

if(ENABLE_HIP AND (TARGET sundials_sunmemhip))
  add_subdirectory(hip)
endif()

if(ENABLE_SYCL AND (TARGET sundials_sunmemsycl))
  add_subdirectory(sycl)
endif()

