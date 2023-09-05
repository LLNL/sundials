/* -----------------------------------------------------------------------------
 * Programmer(s): David J. Gardner @ LLNL
 * -----------------------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2022, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * ---------------------------------------------------------------------------*/

#include <iostream>
#include <ostream>
#include <string>
#include <chrono>
#include <thread>
#include <cstdio>

#include "sundials/sundials_profiler.h"
#include "sundials/sundials_types.h"

int sleep(SUNProfiler prof, int sec)
{
  int flag = SUNProfiler_Begin(prof, "sleep");
  if (flag) return flag;
  std::this_thread::sleep_for(std::chrono::seconds(sec));
  flag = SUNProfiler_End(prof, "sleep");
  if (flag) return flag;
  return 0;
}

int print_timings(SUNProfiler prof)
{
  // Output timing in default (table) format
  int flag = SUNProfiler_Print(prof, stdout);
  if (flag)
  {
    std::cerr << "SUNProfiler_Print returned " << flag << "\n";
    return 1;
  }

  flag = SUNProfiler_Print(prof, stdout);
  if (flag)
  {
    std::cerr << "SUNProfiler_Print returned " << flag << "\n";
    return 1;
  }

  return 0;
}

int main()
{
  // -----
  // Setup
  // -----

  std::cout << "Testing SUNProfiler\n";

  SUNProfiler prof = nullptr;
  int flag = SUNProfiler_Create(nullptr, "SUNProfiler Test", &prof);
  if (flag)
  {
    std::cerr << "SUNProfiler_Create returned " << flag << "\n";
    return 1;
  }

  // ------
  // Test 1
  // ------

  std::cout << "\nTest 1: sleep 1s, print timings \n";

  flag = sleep(prof, 1);
  if (flag)
  {
    std::cerr << "sleep returned " << flag << "\n";
    return 1;
  }

  flag = print_timings(prof);
  if (flag)
  {
    std::cerr << "print_timings returned " << flag << "\n";
    return 1;
  }

  // ------
  // Test 2
  // ------

  std::cout << "\nTest 2: reset, sleep 2s, print timings\n";

  // Reset timing
  flag = SUNProfiler_Reset(prof);
  if (flag)
  {
    std::cerr << "SUNProfiler_Reset returned " << flag << "\n";
    return 1;
  }

  flag = sleep(prof, 2);
  if (flag)
  {
    std::cerr << "sleep returned " << flag << "\n";
    return 1;
  }

  flag = print_timings(prof);
  if (flag)
  {
    std::cerr << "print_timings returned " << flag << "\n";
    return 1;
  }

  // ------
  // Test 3
  // ------

  std::cout << "\nTest 3: multiple outputs to a file for plot test\n";

  std::FILE* fout = std::fopen("profiling_test_output.txt", "w");
  if (!fout)
  {
    std::cerr << "fopen returned a null pointer\n";
    return 1;
  }

  for (int i = 1; i < 4; i++)
  {
    // Reset timing
    flag = SUNProfiler_Reset(prof);
    if (flag)
    {
      std::cerr << "SUNProfiler_Reset returned " << flag << "\n";
      return 1;
    }

    flag = sleep(prof, i);
    if (flag)
    {
      std::cerr << "sleep returned " << flag << "\n";
      return 1;
    }

    flag = SUNProfiler_Print(prof, fout);
    if (flag)
    {
      std::cerr << "SUNProfiler_Print returned " << flag << "\n";
      return 1;
    }
  }

  std::fclose(fout);

  // --------
  // Clean up
  // --------

  flag = SUNProfiler_Free(&prof);
  if (flag)
  {
    std::cerr << "SUNProfiler_Free returned " << flag << "\n";
    return 1;
  }

  std::cout << "\nTest complete\n";

  return 0;
}
