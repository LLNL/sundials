/* -----------------------------------------------------------------
 * Programmer(s): Cody J. Balos @ LLNL
 * -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2025-2025, Lawrence Livermore National Security,
 * University of Maryland Baltimore County, and the SUNDIALS contributors.
 * Copyright (c) 2013-2025, Lawrence Livermore National Security
 * and Southern Methodist University.
 * Copyright (c) 2002-2013, Lawrence Livermore National Security.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------
 * This file defines the pysundials Python module and includes all
 * of the submodule pieces.
 * -----------------------------------------------------------------*/

#include <nanobind/nanobind.h>

namespace nb = nanobind;

namespace pysundials {

//
// Forward declarations of all of the binding functions
//

void bind_core(nb::module_& m);
void bind_arkode(nb::module_& m);
void bind_cvodes(nb::module_& m);
void bind_idas(nb::module_& m);
void bind_kinsol(nb::module_& m);

void bind_nvector_serial(nb::module_& m);

void bind_sunlinsol_spgmr(nb::module_& m);
void bind_sunlinsol_dense(nb::module_& m);
void bind_sunlinsol_band(nb::module_& m);
void bind_sunlinsol_spbcgs(nb::module_& m);
void bind_sunlinsol_spfgmr(nb::module_& m);
void bind_sunlinsol_sptfqmr(nb::module_& m);
void bind_sunlinsol_pcg(nb::module_& m);

void bind_sunmatrix_band(nb::module_& m);
void bind_sunmatrix_dense(nb::module_& m);
void bind_sunmatrix_sparse(nb::module_& m);

void bind_sunnonlinsol_newton(nb::module_& m);
void bind_sunnonlinsol_fixedpoint(nb::module_& m);

void bind_sunadaptcontroller_mrihtol(nb::module_& m);
void bind_sunadaptcontroller_soderlind(nb::module_& m);
void bind_sunadaptcontroller_imexgus(nb::module_& m);

// SUNDomEigEst bindings
void bind_sundomeigest_arnoldi(nb::module_& m);
void bind_sundomeigest_power(nb::module_& m);

// SUNAdjointCheckpointScheme bindings
void bind_sunadjointcheckpointscheme_fixed(nb::module_& m);

} // namespace pysundials

//
// Define main module, pysundials, and all of its submodules
//

NB_MODULE(pysundials, m)
{
  nb::module_ core_m = m.def_submodule("core", "A submodule of 'pysundials'");
  pysundials::bind_core(core_m);

  nb::module_ arkode_m = m.def_submodule("arkode",
                                         "A submodule of 'pysundials'");
  pysundials::bind_arkode(arkode_m);

  nb::module_ cvodes_m = m.def_submodule("cvodes",
                                         "A submodule of 'pysundials'");
  pysundials::bind_cvodes(cvodes_m);

  nb::module_ idas_m = m.def_submodule("idas", "A submodule of 'pysundials'");
  pysundials::bind_idas(idas_m);

  nb::module_ kinsol_m = m.def_submodule("kinsol",
                                         "A submodule of 'pysundials'");
  pysundials::bind_kinsol(kinsol_m);

  // Bind all implementation modules directly to core_m
  pysundials::bind_nvector_serial(core_m);
  pysundials::bind_sunlinsol_spgmr(core_m);
  pysundials::bind_sunlinsol_dense(core_m);
  pysundials::bind_sunlinsol_band(core_m);
  pysundials::bind_sunlinsol_spbcgs(core_m);
  pysundials::bind_sunlinsol_spfgmr(core_m);
  pysundials::bind_sunlinsol_sptfqmr(core_m);
  pysundials::bind_sunlinsol_pcg(core_m);

  pysundials::bind_sunmatrix_dense(core_m);
  pysundials::bind_sunmatrix_band(core_m);
  pysundials::bind_sunmatrix_sparse(core_m);

  pysundials::bind_sunnonlinsol_newton(core_m);
  pysundials::bind_sunnonlinsol_fixedpoint(core_m);

  pysundials::bind_sunadaptcontroller_mrihtol(core_m);
  pysundials::bind_sunadaptcontroller_soderlind(core_m);
  pysundials::bind_sunadaptcontroller_imexgus(core_m);

  // TODO(CJB): enable arnoldi once LAPACK is enabled for pysundials
  // pysundials::bind_sundomeigest_arnoldi(core_m);
  pysundials::bind_sundomeigest_power(core_m);

  pysundials::bind_sunadjointcheckpointscheme_fixed(core_m);
}