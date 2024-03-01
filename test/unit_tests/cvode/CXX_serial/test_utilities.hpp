/* -----------------------------------------------------------------------------
 * Programmer(s): David J. Gardner @ LLNL
 * -----------------------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2023, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------------------
 * Utility functions
 * ---------------------------------------------------------------------------*/

#include <iostream>

#include "sundials/sundials_config.h"

// Check function return flag
int check_flag(int flag, const std::string funcname)
{
  if (!flag) return 0;
  if (flag < 0) std::cerr << "ERROR: ";
  std::cerr << funcname << " returned " << flag << std::endl;
  return 1;
}

// Check if a function returned a NULL pointer
int check_ptr(void *ptr, const std::string funcname)
{
  if (ptr) return 0;
  std::cerr << "ERROR: " << funcname << " returned NULL" << std::endl;
  return 1;
}

inline void find_arg(std::vector<std::string> &args, const std::string key,
                     sunrealtype &dest)
{
  auto it = find(args.begin(), args.end(), key);
  if (it != args.end())
  {
#if defined(SUNDIALS_SINGLE_PRECISION)
    dest = stof(*(it + 1));
#elif defined(SUNDIALS_DOUBLE_PRECISION)
    dest = stod(*(it + 1));
#elif defined(SUNDIALS_EXTENDED_PRECISION)
    dest = stold(*(it + 1));
#endif
    args.erase(it, it + 2);
  }
}

inline void find_arg(std::vector<std::string> &args, const std::string key,
                     long int &dest)
{
  auto it = find(args.begin(), args.end(), key);
  if (it != args.end())
  {
    dest = stoll(*(it + 1));
    args.erase(it, it + 2);
  }
}

inline void find_arg(std::vector<std::string> &args, const std::string key,
                     int &dest)
{
  auto it = find(args.begin(), args.end(), key);
  if (it != args.end())
  {
    dest = stoi(*(it + 1));
    args.erase(it, it + 2);
  }
}

inline void find_arg(std::vector<std::string> &args, const std::string key,
                     bool &dest, bool store = true)
{
  auto it = find(args.begin(), args.end(), key);
  if (it != args.end())
  {
    dest = store;
    args.erase(it);
  }
}

/*---- end of file ----*/
