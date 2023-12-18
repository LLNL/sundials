/* -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2023, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------*/

#ifndef SUN_DUMPSTDERR_HPP_
#define SUN_DUMPSTDERR_HPP_

#include <fstream>
#include <string>
#include <sundials/sundials_core.hpp>

std::string dumpstderr(SUNContext sunctx, const std::string& errfile)
{
  SUNLogger logger = NULL;
  SUNContext_GetLogger(sunctx, &logger);
  SUNLogger_Flush(logger, SUN_LOGLEVEL_ERROR);
  std::ifstream file(errfile);
  std::string line;
  std::string file_contents;
  while (std::getline(file, line))
  {
    file_contents += line;
    file_contents.push_back('\n');
  }
  return file_contents;
}

#endif /* SUN_DUMPSTDERR_HPP_ */
