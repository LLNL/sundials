/* -----------------------------------------------------------------
 * Programmer(s): Cody J. Balos @ LLNL
 * -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2023, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------
 * SUNDIALS context class. A context object holds data that all
 * SUNDIALS objects in a simulation share.
 * ----------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sundials/sundials_context.h>
#include <sundials/sundials_logger.h>
#include <sundials/sundials_profiler.h>

#include "sundials_context_impl.h"
#include "sundials_debug.h"

#ifdef SUNDIALS_ADIAK_ENABLED
#include <adiak.h>
void sunAdiakCollectMetadata();
#endif

int SUNContext_Create(SUNComm comm, SUNContext* sunctx)
{
  SUNProfiler profiler = NULL;
  SUNLogger logger     = NULL;

#if defined(SUNDIALS_BUILD_WITH_PROFILING) && !defined(SUNDIALS_CALIPER_ENABLED)
  if (SUNProfiler_Create(comm, "SUNContext Default", &profiler)) return (-1);
#endif

#ifdef SUNDIALS_ADIAK_ENABLED 
  adiak_init(&comm);
  sunAdiakCollectMetadata();
#endif

#if SUNDIALS_LOGGING_LEVEL > 0 
#if SUNDIALS_MPI_ENABLED
  if (SUNLogger_CreateFromEnv(comm, &logger))
#else
  if (SUNLogger_CreateFromEnv(SUN_COMM_NULL, &logger))
#endif
  {
    return (-1);
  }
#else
#if SUNDIALS_MPI_ENABLED
  if (SUNLogger_Create(comm, 0, &logger))
#else
  if (SUNLogger_Create(SUN_COMM_NULL, 0, &logger))
#endif
  {
    return (-1);
  }
  SUNLogger_SetErrorFilename(logger, "");
  SUNLogger_SetWarningFilename(logger, "");
  SUNLogger_SetInfoFilename(logger, "");
  SUNLogger_SetDebugFilename(logger, "");
#endif

  *sunctx = NULL;
  *sunctx = (SUNContext)malloc(sizeof(struct _SUNContext));

  if (*sunctx == NULL)
  {
#if defined(SUNDIALS_BUILD_WITH_PROFILING) && !defined(SUNDIALS_CALIPER_ENABLED)
    SUNProfiler_Free(&profiler);
#endif
    SUNLogger_Destroy(&logger);
    return (-1);
  }

  (*sunctx)->logger       = logger;
  (*sunctx)->own_logger   = logger != NULL;
  (*sunctx)->profiler     = profiler;
  (*sunctx)->own_profiler = profiler != NULL;

  return (0);
}

int SUNContext_GetProfiler(SUNContext sunctx, SUNProfiler* profiler)
{
  if (sunctx == NULL)
  {
    return (-1);
  }

#ifdef SUNDIALS_BUILD_WITH_PROFILING
  /* get profiler */
  *profiler = sunctx->profiler;
#else
  *profiler = NULL;
#endif

  return (0);
}

int SUNContext_SetProfiler(SUNContext sunctx, SUNProfiler profiler)
{
  if (sunctx == NULL)
  {
    return (-1);
  }

#ifdef SUNDIALS_BUILD_WITH_PROFILING
  /* free any existing profiler */
  if (sunctx->profiler && sunctx->own_profiler)
  {
    if (SUNProfiler_Free(&(sunctx->profiler))) return (-1);
    sunctx->profiler = NULL;
  }

  /* set profiler */
  sunctx->profiler     = profiler;
  sunctx->own_profiler = SUNFALSE;
#endif

  return (0);
}

int SUNContext_GetLogger(SUNContext sunctx, SUNLogger* logger)
{
  if (sunctx == NULL)
  {
    return (-1);
  }

  /* get logger */
  *logger = sunctx->logger;

  return (0);
}

int SUNContext_SetLogger(SUNContext sunctx, SUNLogger logger)
{
  if (sunctx == NULL)
  {
    return (-1);
  }

  /* free any existing logger */
  if (sunctx->logger && sunctx->own_logger)
  {
    if (SUNLogger_Destroy(&(sunctx->logger)))
    {
      return (-1);
    }
    sunctx->logger = NULL;
  }

  /* set logger */
  sunctx->logger     = logger;
  sunctx->own_logger = SUNFALSE;

  return (0);
}

int SUNContext_Free(SUNContext* sunctx)
{
#if defined(SUNDIALS_BUILD_WITH_PROFILING) && !defined(SUNDIALS_CALIPER_ENABLED)
  FILE* fp;
  char* sunprofiler_print_env;
#endif

  if (!sunctx)
  {
    return (0);
  }
  if (!(*sunctx))
  {
    return (0);
  }

#if defined(SUNDIALS_BUILD_WITH_PROFILING) && !defined(SUNDIALS_CALIPER_ENABLED)
  /* Find out where we are printing to */
  sunprofiler_print_env = getenv("SUNPROFILER_PRINT");
  fp                    = NULL;
  if (sunprofiler_print_env)
  {
    if (!strcmp(sunprofiler_print_env, "0"))
      fp = NULL;
    else if (!strcmp(sunprofiler_print_env, "1") ||
             !strcmp(sunprofiler_print_env, "TRUE") ||
             !strcmp(sunprofiler_print_env, "stdout"))
      fp = stdout;
    else
      fp = fopen(sunprofiler_print_env, "a");
  }

  /* Enforce that the profiler is freed before finalizing,
     if it is not owned by the sunctx. */
  if ((*sunctx)->profiler)
  {
    if (fp) SUNProfiler_Print((*sunctx)->profiler, fp);
    if (fp) fclose(fp);
    if ((*sunctx)->own_profiler) SUNProfiler_Free(&(*sunctx)->profiler);
  }
#endif

#ifdef SUNDIALS_ADIAK_ENABLED
  adiak_fini();
#endif

  if ((*sunctx)->logger && (*sunctx)->own_logger)
  {
    SUNLogger_Destroy(&(*sunctx)->logger);
  }

  free(*sunctx);
  *sunctx = NULL;

  return (0);
}

#ifdef SUNDIALS_ADIAK_ENABLED
void sunAdiakCollectMetadata() {
  adiak_launchdate();
  adiak_executable();
  adiak_cmdline();
  adiak_clustername();

  adiak_job_size();
  adiak_num_hosts();

  adiak_namevalue("c_compiler", 2, NULL, "%s", SUN_C_COMPILER);
  adiak_namevalue("c_compiler_version", 2, NULL, "%s", SUN_C_COMPILER_VERSION);
  adiak_namevalue("c_compiler_flags", 2, NULL, "%s", SUN_C_COMPILER_FLAGS);

  adiak_namevalue("cxx_compiler", 2, NULL, "%s", SUN_CXX_COMPILER);
  adiak_namevalue("cxx_compiler_version", 2, NULL, "%s", SUN_CXX_COMPILER_VERSION);
  adiak_namevalue("cxx_compiler_flags", 2, NULL, "%s", SUN_CXX_COMPILER_FLAGS);

  adiak_namevalue("fortran_compiler", 2, NULL, "%s", SUN_FORTRAN_COMPILER);
  adiak_namevalue("fortran_compiler_version", 2, NULL, "%s", SUN_FORTRAN_COMPILER_VERSION);
  adiak_namevalue("fortran_compiler_flags", 2, NULL, "%s", SUN_FORTRAN_COMPILER_FLAGS);

  adiak_namevalue("sundials_version", 2, NULL, "%s", SUNDIALS_VERSION);
  adiak_namevalue("sundials_git_version", 2, NULL, "%s", SUNDIALS_GIT_VERSION);
  adiak_namevalue("build_type", 2, NULL, "%s", SUN_BUILD_TYPE);
  adiak_namevalue("third_party_libraries", 2, NULL, "%s", SUN_TPL_LIST);
#ifdef SUN_JOB_ID
  adiak_namevalue("job_id", 2, NULL, "%s", SUN_JOB_ID);
#endif
  adiak_namevalue("job_start_time", 2, NULL, "%s", SUN_JOB_START_TIME);

#ifdef SUNDIALS_SPACK_VERSION
  adiak_namevalue("spack_version", 2, NULL, "%s", SUNDIALS_SPACK_VERSION);
#endif

#ifdef SUNDIALS_GINKGO_ENABLED
  adiak_namevalue("ginkgo_version", 2, NULL, "%s", SUN_GINKGO_VERSION);
#endif

#ifdef SUNDIALS_HYPRE_ENABLED
  adiak_namevalue("hypre_version", 2, NULL, "%s", SUN_HYPRE_VERSION);
#endif

#ifdef SUNDIALS_KLU_ENABLED
  adiak_namevalue("klu_version", 2, NULL, "%s", SUN_KLU_VERSION);
#endif

#ifdef SUNDIALS_KOKKOS_ENABLED
  adiak_namevalue("kokkos_version", 2, NULL, "%s", SUN_KOKKOS_VERSION);
#endif

#ifdef SUNDIALS_KOKKOS_KERNELS_ENABLED
  adiak_namevalue("kokkos_kernels_version", 2, NULL, "%s", SUN_KOKKOS_KERNELS_VERSION);
#endif

#ifdef SUNDIALS_BLAS_LAPACK_ENABLED
  adiak_namevalue("lapack_version", 2, NULL, "%s", SUN_LAPACK_VERSION);
#endif

#ifdef SUNDIALS_MAGMA_ENABLED
  adiak_namevalue("magma_version", 2, NULL, "%s", SUN_MAGMA_VERSION);
#endif

#if SUNDIALS_MPI_ENABLED
  adiak_namevalue("mpi_c_compiler", 2, NULL, "%s", SUN_MPI_C_COMPILER);
  adiak_namevalue("mpi_c_version", 2, NULL, "%s", SUN_MPI_C_VERSION);

  adiak_namevalue("mpi_cxx_compiler", 2, NULL, "%s", SUN_MPI_CXX_COMPILER);
  adiak_namevalue("mpi_cxx_version", 2, NULL, "%s", SUN_MPI_CXX_VERSION);

  adiak_namevalue("mpi_fortran_compiler", 2, NULL, "%s", SUN_MPI_FORTRAN_COMPILER);
  adiak_namevalue("mpi_fortran_version", 2, NULL, "%s", SUN_MPI_FORTRAN_VERSION);
#endif

#ifdef SUNDIALS_ONEMKL_ENABLED
  adiak_namevalue("onemkl_version", 2, NULL, "%s", SUN_ONEMKL_VERSION);
#endif

#ifdef SUNDIALS_OPENMP_ENABLED
  adiak_namevalue("openmp_version", 2, NULL, "%s", SUN_OPENMP_VERSION);
#endif

#ifdef SUNDIALS_PETSC_ENABLED
  adiak_namevalue("petsc_version", 2, NULL, "%s", SUN_PETSC_VERSION);
#endif

#ifdef SUNDIALS_PTHREADS_ENABLED
  adiak_namevalue("pthreads_version", 2, NULL, "%s", SUN_PTHREADS_VERSION);
#endif

#ifdef SUNDIALS_RAJA_ENABLED
  adiak_namevalue("raja_version", 2, NULL, "%s", SUN_RAJA_VERSION);
#endif

#ifdef SUNDIALS_SUPERLUDIST_ENABLED
  adiak_namevalue("superludist_version", 2, NULL, "%s", SUN_SUPERLUDIST_VERSION);
#endif

#ifdef SUNDIALS_SUPERLUMT_ENABLED
  adiak_namevalue("superlumt_version", 2, NULL, "%s", SUN_SUPERLUMT_VERSION);
#endif

#ifdef SUNDIALS_TRILLINOS_ENABLED
  adiak_namevalue("trillinos_version", 2, NULL, "%s", SUN_TRILLINOS_VERSION);
#endif

#ifdef SUNDIALS_XBRAID_ENABLED
  adiak_namevalue("xbraid_version", 2, NULL, "%s", SUN_XBRAID_VERSION);
#endif

#ifdef SUNDIALS_CUDA_ENABLED
  adiak_namevalue("cuda_version", 2, NULL, "%s", SUN_CUDA_VERSION);
  adiak_namevalue("cuda_compiler", 2, NULL, "%s", SUN_CUDA_COMPILER);
  adiak_namevalue("cuda_architectures", 2, NULL, "%s", SUN_CUDA_ARCHITECTURES);
#endif

#ifdef SUNDIALS_HIP_ENABLED
  adiak_namevalue("hip_version", 2, NULL, "%s", SUN_HIP_VERSION);
  adiak_namevalue("amdgpu_targets", 2, NULL, "%s", SUN_AMDGPU_TARGETS);
#endif

}
#endif
