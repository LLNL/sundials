#!/usr/bin/env bash

# Build and test first
source "$BUILD_ROOT/.gitlab/build_and_test.sh" --no-clean

benchmark_dir="${build_dir}/benchmarks"
ar3d_dir="${benchmark_dir}/advection_reaction_3D"
d2d_dir="${benchmark_dir}/diffusion_2D"

nresg=$((BENCHMARK_NNODES * 4)) # Lassen has 4 GPUs per node
nresc=$((BENCHMARK_NNODES * 40)) # Lassen has 40 cores per node

echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
echo "~~~~~ Benchmarking SUNDIALS"
echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"

if [[ ! -d ${benchmark_dir} ]]
then
    echo "ERROR: Benchmark directory not found : ${benchmark_dir}" && exit 1
fi

cd "${benchmark_dir}"

# We use synchronous kernel launches for benchmarking
# because we are interested in individual region timings.
export CUDA_LAUNCH_BLOCKING=1

if [[ -d ${ar3d_dir} ]]
then
    date

    if [[ -f ${ar3d_dir}/advection_reaction_3D_mpicuda ]]
    then
        export CALI_CONFIG="spot(output=${BUILD_ROOT}/ar3d_arkimex_tlnewton.cali),runtime-report(calc.inclusive)"
        jsrun --smpiargs="-gpu" -n${nresg} -a1 -c1 -g1 "${ar3d_dir}/advection_reaction_3D_mpicuda" --method ARK-IMEX --nls tl-newton --tf 10.0 --dont-save
        export CALI_CONFIG="spot(output=${BUILD_ROOT}/ar3d_arkdirk_newton.cali),runtime-report(calc.inclusive)"
        jsrun --smpiargs="-gpu" -n${nresg} -a1 -c1 -g1 "${ar3d_dir}/advection_reaction_3D_mpicuda" --method ARK-DIRK --nls newton --tf 10.0 --dont-save
        export CALI_CONFIG="spot(output=${BUILD_ROOT}/ar3d_cvbdf_newton.cali),runtime-report(calc.inclusive)"
        jsrun --smpiargs="-gpu" -n${nresg} -a1 -c1 -g1 "${ar3d_dir}/advection_reaction_3D_mpicuda" --method CV-BDF --nls newton --tf 10.0 --dont-save
        export CALI_CONFIG="spot(output=${BUILD_ROOT}/ar3d_ida_newton.cali),runtime-report(calc.inclusive)"
        jsrun --smpiargs="-gpu" -n${nresg} -a1 -c1 -g1 "${ar3d_dir}/advection_reaction_3D_mpicuda" --method IDA --nls newton --tf 10.0 --dont-save
    elif [[ -f ${ar3d_dir}/advection_reaction_3D ]]
    then
        export CALI_CONFIG="spot(output=${BUILD_ROOT}/ar3d_mpionly_arkimex_tlnewton.cali),runtime-report(calc.inclusive)"
        jsrun -n${nresc} -a1 -c1 "${ar3d_dir}/advection_reaction_3D" --method ARK-IMEX --nls tl-newton --tf 2.0 --dont-save
        export CALI_CONFIG="spot(output=${BUILD_ROOT}/ar3d_mpionly_arkdirk_newton.cali),runtime-report(calc.inclusive)"
        jsrun -n${nresc} -a1 -c1 "${ar3d_dir}/advection_reaction_3D" --method ARK-DIRK --nls newton --tf 2.0 --dont-save
        export CALI_CONFIG="spot(output=${BUILD_ROOT}/ar3d_mpionly_cvbdf_newton.cali),runtime-report(calc.inclusive)"
        jsrun -n${nresc} -a1 -c1 "${ar3d_dir}/advection_reaction_3D" --method CV-BDF --nls newton --tf 2.0 --dont-save
        export CALI_CONFIG="spot(output=${BUILD_ROOT}/ar3d_mpionly_ida_newton.cali),runtime-report(calc.inclusive)"
        jsrun -n${nresc} -a1 -c1 "${ar3d_dir}/advection_reaction_3D" --method IDA --nls newton --tf 2.0 --dont-save
    fi

    date
fi


if [[ -d ${d2d_dir} ]]
then
    date

    if [[ -d ${d2d_dir}/mpi_gpu ]]
    then
        export CALI_CONFIG="spot(output=${BUILD_ROOT}/d2d_arkode.cali),runtime-report(calc.inclusive)"
        jsrun --smpiargs="-gpu" -n${nresg} -a1 -c1 -g1 "${d2d_dir}/mpi_gpu/arkode_diffusion_2D_mpicuda"
        export CALI_CONFIG="spot(output=${BUILD_ROOT}/d2d_cvode.cali),runtime-report(calc.inclusive)"
        jsrun --smpiargs="-gpu" -n${nresg} -a1 -c1 -g1 "${d2d_dir}/mpi_gpu/cvode_diffusion_2D_mpicuda"
        export CALI_CONFIG="spot(output=${BUILD_ROOT}/d2d_ida.cali),runtime-report(calc.inclusive)"
        jsrun --smpiargs="-gpu" -n${nresg} -a1 -c1 -g1 "${d2d_dir}/mpi_gpu/ida_diffusion_2D_mpicuda"
    elif [[ -d ${d2d_dir}/mpi_serial ]]
    then
        export CALI_CONFIG="spot(output=${BUILD_ROOT}/d2d_arkode.cali),runtime-report(calc.inclusive)"
        jsrun -n${nresc} -a1 -c1 "${d2d_dir}/mpi_serial/arkode_diffusion_2D"
        export CALI_CONFIG="spot(output=${BUILD_ROOT}/d2d_cvode.cali),runtime-report(calc.inclusive)"
        jsrun -n${nresc} -a1 -c1 "${d2d_dir}/mpi_serial/arkode_diffusion_2D"
        export CALI_CONFIG="spot(output=${BUILD_ROOT}/d2d_ida.cali),runtime-report(calc.inclusive)"
        jsrun -n${nresc} -a1 -c1 "${d2d_dir}/mpi_serial/ida_diffusion_2D"
    fi

    date
fi

echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
echo "~~~~~ CLEAN UP"
echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
make clean

