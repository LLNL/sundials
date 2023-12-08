/*---------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ SMU
 *---------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2023, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 *----------------------------------------------------------------
 * Routine to test all built-in Butcher tables to check their
 * order of accuracy.
 *-----------------------------------------------------------------*/

// Header files
#include <arkode/arkode.h>
#include <arkode/arkode_butcher_dirk.h>
#include <arkode/arkode_butcher_erk.h>
#include <iostream>
#include <ostream>
#include <string>
#include <sundials/sundials_types.h>
#include <vector>

// Main Program
int main()
{
  // set vectors of individual tables to test
  std::vector<std::string> Tables_ERK =
    {"ARKODE_HEUN_EULER_2_1_2",       "ARKODE_ARK2_ERK_3_1_2",
     "ARKODE_BOGACKI_SHAMPINE_4_2_3", "ARKODE_ARK324L2SA_ERK_4_2_3",
     "ARKODE_ZONNEVELD_5_3_4",        "ARKODE_ARK436L2SA_ERK_6_3_4",
     "ARKODE_SAYFY_ABURUB_6_3_4",     "ARKODE_CASH_KARP_6_4_5",
     "ARKODE_FEHLBERG_6_4_5",         "ARKODE_DORMAND_PRINCE_7_4_5",
     "ARKODE_ARK548L2SA_ERK_8_4_5",   "ARKODE_VERNER_8_5_6",
     "ARKODE_FEHLBERG_13_7_8",        "ARKODE_ARK437L2SA_ERK_7_3_4",
     "ARKODE_ARK548L2SAb_ERK_8_4_5",  "ARKODE_SOFRONIOU_SPALETTA_5_3_4",
     "ARKODE_SHU_OSHER_3_2_3",        "ARKODE_VERNER_9_5_6",
     "ARKODE_VERNER_10_6_7",          "ARKODE_VERNER_13_7_8",
     "ARKODE_VERNER_16_8_9"};
  std::vector<std::string> Tables_DIRK          = {"ARKODE_SDIRK_2_1_2",
                                                   "ARKODE_ARK2_DIRK_3_1_2",
                                                   "ARKODE_BILLINGTON_3_3_2",
                                                   "ARKODE_TRBDF2_3_3_2",
                                                   "ARKODE_KVAERNO_4_2_3",
                                                   "ARKODE_ARK324L2SA_DIRK_4_2_3",
                                                   "ARKODE_CASH_5_2_4",
                                                   "ARKODE_CASH_5_3_4",
                                                   "ARKODE_SDIRK_5_3_4",
                                                   "ARKODE_KVAERNO_5_3_4",
                                                   "ARKODE_ARK436L2SA_DIRK_6_3_4",
                                                   "ARKODE_KVAERNO_7_4_5",
                                                   "ARKODE_ARK548L2SA_DIRK_8_4_5",
                                                   "ARKODE_ARK437L2SA_DIRK_7_3_4",
                                                   "ARKODE_ARK548L2SAb_DIRK_8_4_5",
                                                   "ARKODE_ESDIRK324L2SA_4_2_3",
                                                   "ARKODE_ESDIRK325L2SA_5_2_3",
                                                   "ARKODE_ESDIRK32I5L2SA_5_2_3",
                                                   "ARKODE_ESDIRK436L2SA_6_3_4",
                                                   "ARKODE_ESDIRK43I6L2SA_6_3_4",
                                                   "ARKODE_QESDIRK436L2SA_6_3_4",
                                                   "ARKODE_ESDIRK437L2SA_7_3_4",
                                                   "ARKODE_ESDIRK547L2SA_7_4_5",
                                                   "ARKODE_ESDIRK547L2SA2_7_4_5"};
  std::vector<ARKODE_ERKTableID> Tables_ARK_ERK = {ARKODE_ARK2_ERK_3_1_2,
                                                   ARKODE_ARK324L2SA_ERK_4_2_3,
                                                   ARKODE_ARK436L2SA_ERK_6_3_4,
                                                   ARKODE_ARK437L2SA_ERK_7_3_4,
                                                   ARKODE_ARK548L2SA_ERK_8_4_5,
                                                   ARKODE_ARK548L2SAb_ERK_8_4_5};
  std::vector<ARKODE_DIRKTableID> Tables_ARK_DIRK =
    {ARKODE_ARK2_DIRK_3_1_2,       ARKODE_ARK324L2SA_DIRK_4_2_3,
     ARKODE_ARK436L2SA_DIRK_6_3_4, ARKODE_ARK437L2SA_DIRK_7_3_4,
     ARKODE_ARK548L2SA_DIRK_8_4_5, ARKODE_ARK548L2SAb_DIRK_8_4_5};
  std::vector<std::string> STables_ARK = {"ARKODE_ARK2_3_1_2",
                                          "ARKODE_ARK324L2SA_4_2_3",
                                          "ARKODE_ARK436L2SA_6_3_4",
                                          "ARKODE_ARK437L2SA_7_3_4",
                                          "ARKODE_ARK548L2SA_8_4_5",
                                          "ARKODE_ARK548L2SAb_8_4_5"};
  int numfails                         = 0;

  // loop over individual ERK tables
  std::cout << "\nTesting individual ERK methods:\n\n";
  for (std::string table : Tables_ERK)
  {
    std::cout << "Testing method " << table << ":";

    // load Butcher table
    ARKodeButcherTable B = ARKodeButcherTable_LoadERKByName(table.c_str());
    if (B == NULL)
    {
      std::cout << "  error retrieving table, aborting\n";
      return 1;
    }

    // test method, upon questionable results, re-test with statistics printing
    int q, p, retval;
    retval = ARKodeButcherTable_CheckOrder(B, &q, &p, NULL);
    if (retval == 0)
    {
      std::cout << "  table matches predicted method/embedding orders of " << q
                << "/" << p << "\n";
    }
    else if (retval > 0)
    {
      std::cout << "  WARNING:\n";
      retval = ARKodeButcherTable_CheckOrder(B, &q, &p, stdout);
    }
    else
    {
      std::cout << "  ERROR:\n";
      retval = ARKodeButcherTable_CheckOrder(B, &q, &p, stdout);
      numfails++;
    }

    // clean up after this test
    ARKodeButcherTable_Free(B);
  }

  // loop over individual DIRK tables
  std::cout << "\nTesting individual DIRK methods:\n\n";
  for (std::string table : Tables_DIRK)
  {
    std::cout << "Testing method " << table << ":";

    // load Butcher table
    ARKodeButcherTable B = ARKodeButcherTable_LoadDIRKByName(table.c_str());
    if (B == NULL)
    {
      std::cout << "  error retrieving table, aborting\n";
      return 1;
    }

    // test method, upon questionable results, re-test with statistics printing
    int q, p, retval;
    retval = ARKodeButcherTable_CheckOrder(B, &q, &p, NULL);
    if (retval == 0)
    {
      std::cout << "  table matches predicted method/embedding orders of " << q
                << "/" << p << "\n";
    }
    else if (retval > 0)
    {
      std::cout << "  WARNING:\n";
      retval = ARKodeButcherTable_CheckOrder(B, &q, &p, stdout);
    }
    else
    {
      std::cout << "  ERROR:\n";
      retval = ARKodeButcherTable_CheckOrder(B, &q, &p, stdout);
      numfails++;
    }

    // clean up after this test
    ARKodeButcherTable_Free(B);
  }

  // loop over ARK pairs
  std::cout << "\nTesting ARK pairs:\n\n";
  for (size_t i = 0; i < Tables_ARK_ERK.size(); i++)
  {
    std::cout << "Testing method " << STables_ARK[i] << ":";

    // load Butcher tables
    ARKodeButcherTable Be = ARKodeButcherTable_LoadERK(Tables_ARK_ERK[i]);
    if (Be == NULL)
    {
      std::cout << "  error retrieving explicit table, aborting\n";
      return 1;
    }
    ARKodeButcherTable Bi = ARKodeButcherTable_LoadDIRK(Tables_ARK_DIRK[i]);
    if (Bi == NULL)
    {
      std::cout << "  error retrieving implicit table, aborting";
      return 1;
    }

    // test method, upon questionable results, re-test with statistics printing
    int q, p, retval;
    retval = ARKodeButcherTable_CheckARKOrder(Be, Bi, &q, &p, NULL);
    if (retval == 0)
    {
      std::cout << "  Method/embedding match predicted orders of " << q << "/"
                << p << "\n";
    }
    else if (retval > 0)
    {
      std::cout << "  WARNING:\n";
      retval = ARKodeButcherTable_CheckARKOrder(Be, Bi, &q, &p, stdout);
    }
    else
    {
      std::cout << "  ERROR:\n";
      retval = ARKodeButcherTable_CheckARKOrder(Be, Bi, &q, &p, stdout);
      numfails++;
    }

    // clean up after this test
    ARKodeButcherTable_Free(Be);
    ARKodeButcherTable_Free(Bi);
  }

  // determine overall success/failure and return
  if (numfails == 0)
  {
    std::cout << "All Butcher tables passed\n";
    return 0;
  }
  else
  {
    std::cout << numfails << " Butcher tables failed\n";
    return 1;
  }
}

/*---- end of file ----*/
