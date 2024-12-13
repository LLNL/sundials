/* -----------------------------------------------------------------------------
 * Programmer(s): David J. Gardner @ LLNL
 * -----------------------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2024, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------------------
 * Utility functions for checking function returns
 * ---------------------------------------------------------------------------*/

#include <iostream>
#include <string>

// Check function for non-zero return flag
static int check_flag(int flag, const std::string funcname)
{
  if (!flag) { return 0; }
  if (flag < 0) { std::cerr << "ERROR: "; }
  else { std::cerr << "WARNING: "; }
  std::cerr << funcname << " returned " << flag << std::endl;
  return (flag < 0) ? 1 : 0;
}

// Check if a function returned a NULL pointer
static int check_ptr(void* ptr, const std::string funcname)
{
  if (ptr) { return 0; }
  std::cerr << "ERROR: " << funcname << " returned NULL" << std::endl;
  return 1;
}
