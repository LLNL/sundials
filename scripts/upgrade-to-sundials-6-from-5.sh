#!/bin/bash
# -----------------------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2002-2021, Lawrence Livermore National Security
# and Southern Methodist University.
# All rights reserved.
#
# See the top-level LICENSE and NOTICE files for details.
#
# SPDX-License-Identifier: BSD-3-Clause
# SUNDIALS Copyright End
# -----------------------------------------------------------------------------
# This script is a utility to upgrade codes from SUNDIALS v5.x.x to v6.0.0.
# It does a few things:
#   1. Updates all SUNDIALS constructor invocations to pass 'SUNCTX_PLACEHOLDER'
#      as the SUNContext argument. Users will need to then manually replace
#      'SUNCTX_PLACEHOLDER' with an actual SUNContext object.
#
# Usage:
#   ./upgrade-to-sundials-6-from-5.sh <files to update>
# -----------------------------------------------------------------------------

updateConstructors() {
  # All of the constructors updated to take a SUNContext.
  # Obtained with
  #   clang-format include/**/*.h | grep 'SUNContext sunctx)' | sed 's/SUNDIALS_EXPORT//' | awk -F ' ' '{print $2}' > constructors.txt
  # where the .clang-format file sets the ColumnLimit to 1000.
  declare -a constructors=(
    "ARKStepCreate"
    "CVodeCreate"
    "CVodeCreate"
    "ERKStepCreate"
    "IDACreate"
    "IDACreate"
    "KINCreate"
    "MRIStepCreate"
    "N_VMake_Cuda"
    "N_VMake_Hip"
    "N_VMake_MPIManyVector"
    "N_VMake_MPIPlusX"
    "N_VMake_OpenMP"
    "N_VMake_OpenMPDEV"
    "N_VMake_Parallel"
    "N_VMake_ParHyp"
    "N_VMake_Petsc"
    "N_VMake_Pthreads"
    "N_VMake_Raja"
    "N_VMake_Serial"
    "N_VMake_Sycl"
    "N_VMakeManaged_Cuda"
    "N_VMakeManaged_Hip"
    "N_VMakeManaged_Raja"
    "N_VMakeManaged_Sycl"
    "N_VNew_Cuda"
    "N_VNew_Hip"
    "N_VNew_ManyVector"
    "N_VNew_MPIManyVector"
    "N_VNew_OpenMP"
    "N_VNew_OpenMPDEV"
    "N_VNew_Parallel"
    "N_VNew_Pthreads"
    "N_VNew_Raja"
    "N_VNew_Serial"
    "N_VNew_Sycl"
    "N_VNewEmpty"
    "N_VNewEmpty_Cuda"
    "N_VNewEmpty_Hip"
    "N_VNewEmpty_OpenMP"
    "N_VNewEmpty_OpenMPDEV"
    "N_VNewEmpty_Parallel"
    "N_VNewEmpty_ParHyp"
    "N_VNewEmpty_Petsc"
    "N_VNewEmpty_Pthreads"
    "N_VNewEmpty_Raja"
    "N_VNewEmpty_SensWrapper"
    "N_VNewEmpty_Serial"
    "N_VNewEmpty_Sycl"
    "N_VNewEmpty_Trilinos"
    "N_VNewManaged_Cuda"
    "N_VNewManaged_Hip"
    "N_VNewManaged_Raja"
    "N_VNewManaged_Sycl"
    "N_VNewWithMemHelp_Cuda"
    "N_VNewWithMemHelp_Hip"
    "N_VNewWithMemHelp_Raja"
    "N_VNewWithMemHelp_Sycl"
    "N_VSetContext"
    "N_VSetContext_ManyVector"
    "N_VSetContext_MPIManyVector"
    "SUNBandMatrix"
    "SUNBandMatrixStorage"
    "SUNDenseMatrix"
    "SUNLinSol_Band"
    "SUNLinSol_cuSolverSp_batchQR"
    "SUNLinSol_Dense"
    "SUNLinSol_KLU"
    "SUNLinSol_LapackBand"
    "SUNLinSol_LapackDense"
    "SUNLinSol_MagmaDense"
    "SUNLinSol_OneMklDense"
    "SUNLinSol_PCG"
    "SUNLinSol_SPBCGS"
    "SUNLinSol_SPFGMR"
    "SUNLinSol_SPGMR"
    "SUNLinSol_SPTFQMR"
    "SUNLinSol_SuperLUDIST"
    "SUNLinSol_SuperLUMT"
    "SUNLinSolNewEmpty"
    "SUNLinSolSetContext"
    "SUNMatNewEmpty"
    "SUNMatrix_cuSparse_MakeCSR"
    "SUNMatrix_cuSparse_NewBlockCSR"
    "SUNMatrix_cuSparse_NewCSR"
    "SUNMatrix_MagmaDense"
    "SUNMatrix_MagmaDenseBlock"
    "SUNMatrix_OneMklDense"
    "SUNMatrix_OneMklDenseBlock"
    "SUNMatrix_SLUNRloc"
    "SUNMatSetContext"
    "SUNMemoryHelper_Cuda"
    "SUNMemoryHelper_Hip"
    "SUNMemoryHelper_NewEmpty"
    "SUNMemoryHelper_SetContext"
    "SUNMemoryHelper_Sycl"
    "SUNMemoryHelper_Sys"
    "SUNNonlinSol_FixedPoint"
    "SUNNonlinSol_FixedPointSens"
    "SUNNonlinSol_Newton"
    "SUNNonlinSol_NewtonSens"
    "SUNNonlinSol_PetscSNES"
    "SUNNonlinSolNewEmpty"
    "SUNNonlinSolSetContext"
    "SUNSparseMatrix"
    )

  echo "replacing in ${@}"
  for constructor in "${constructors[@]}"
  do
    perl -pi -e 's/'"${constructor}"'\((.*?)\)/'"${constructor}"'($1, SUNCTX_PLACEHOLDER)/g;' ${@}
  done
}

while true; do
    echo ""
    read -p "> This script will edit files in-place. Are you sure you want to proceed? [Y/N or ls] " yn
    case $yn in
        [Yy]* ) break;;
        [Nn]* ) exit;;
        ls ) echo "> The following files will be edited:" && printf "  %s\n" "${@}" && exit;;
        * ) echo "> Please answer 'ls' to see the list of files that will be edited, 'no' to cancel, or 'yes' to proceed.";;
    esac
done
updateConstructors "${@}"
echo "> SUNDIALS constructors have been updated with a new 'SUNCTX_PLACEHOLDER' argument. 'SUNCTX_PLACEHOLDER' MUST BE REPLACED with a 'SUNContext' object from a call to 'SUNContext_Create'."