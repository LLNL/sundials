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
 * This header file should *NOT* be included in user codes and exists *ONLY* to
 * reduce duplicate utility functions across example programs.
 * ---------------------------------------------------------------------------*/

#include <algorithm>
#include <iostream>
#include <string>
#include <vector>

// Check for an unrecoverable (negative) return value from a SUNDIALS function
int check_flag(const int flag, const std::string funcname)
{
  if (flag < 0)
  {
    std::cerr << "ERROR: " << funcname << " returned " << flag << std::endl;
    return 1;
  }
  return 0;
}

// Check if a function returned a NULL pointer
int check_ptr(const void* ptr, const std::string funcname)
{
  if (ptr) { return 0; }
  std::cerr << "ERROR: " << funcname << " returned NULL" << std::endl;
  return 1;
}

// Functions for parsing vectors of command line inputs
inline void find_arg(std::vector<std::string>& args, const std::string key,
                     float& dest)
{
  auto it = std::find(args.begin(), args.end(), key);
  if (it != args.end())
  {
    dest = stof(*(it + 1));
    args.erase(it, it + 2);
  }
}

inline void find_arg(std::vector<std::string>& args, const std::string key,
                     double& dest)
{
  auto it = std::find(args.begin(), args.end(), key);
  if (it != args.end())
  {
    dest = stod(*(it + 1));
    args.erase(it, it + 2);
  }
}

inline void find_arg(std::vector<std::string>& args, const std::string key,
                     long double& dest)
{
  auto it = std::find(args.begin(), args.end(), key);
  if (it != args.end())
  {
    dest = stold(*(it + 1));
    args.erase(it, it + 2);
  }
}

inline void find_arg(std::vector<std::string>& args, const std::string key,
                     long long& dest)
{
  auto it = std::find(args.begin(), args.end(), key);
  if (it != args.end())
  {
    dest = stoll(*(it + 1));
    args.erase(it, it + 2);
  }
}

inline void find_arg(std::vector<std::string>& args, const std::string key,
                     long int& dest)
{
  auto it = std::find(args.begin(), args.end(), key);
  if (it != args.end())
  {
    dest = stol(*(it + 1));
    args.erase(it, it + 2);
  }
}

inline void find_arg(std::vector<std::string>& args, const std::string key,
                     int& dest)
{
  auto it = std::find(args.begin(), args.end(), key);
  if (it != args.end())
  {
    dest = stoi(*(it + 1));
    args.erase(it, it + 2);
  }
}

inline void find_arg(std::vector<std::string>& args, const std::string key,
                     bool& dest, bool store = true)
{
  auto it = std::find(args.begin(), args.end(), key);
  if (it != args.end())
  {
    dest = store;
    args.erase(it);
  }
}
