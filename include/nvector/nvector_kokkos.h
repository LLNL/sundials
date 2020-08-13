#ifndef _NVECTOR_KOKKOS_H
#define _NVECTOR_KOKKOS_H

#include <stdio.h>
#include <sundials/sundials_nvector.h>
#include <sundials/sundials_config.h>

//TODO: determine best place to include this - do I need sundials_config/what is
#if defined(KOKKOS_ENABLE_CUDA_UVM)
  #define MemSpace Kokkos::CudaUVMSpace
#elif defined(KOKKOS_ENABLE_CUDA)
  #define MemSpace Kokkos::CudaSpace
#elif defined(KOKKOS_ENABLE_HIP)
  #define MemSpace Kokkos::Experimental::HIPSpace
#elif defined(KOKKOS_ENABLE_OPENMPTARGET)
  #define MemSpace Kokkos::OpenMPTargetSpace
#else
  #define MemSpace Kokkos::HostSpace
#endif

using ExecSpace = MemSpace::execution_space;
using range_policy = Kokkos::RangePolicy<ExecSpace>;

typedef Kokkos::View<realtype*, MemSpace> DeviceArrayView;
typedef Kokkos::View<realtype*, Kokkos::HostSpace> HostArrayView;

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

struct _N_VectorContent_Kokkos {
  sunindextype    length;
  HostArrayView   device_data;
  DeviceArrayView host_data;
};

typedef struct _N_VectorContent_Kokkos *N_VectorContent_Kokkos;

SUNDIALS_EXPORT N_Vector N_VNew_Kokkos(sunindextype length);

SUNDIALS_EXPORT N_Vector N_VNewEmpty_Kokkos();

SUNDIALS_EXPORT N_Vector N_VMake_Kokkos(sunindextype length, realtype *h_vdata, realtype *d_vdata);

SUNDIALS_EXPORT sunindextype N_VGetLength_Kokkos(N_Vector v);

SUNDIALS_EXPORT realtype *N_VGetHostArrayPointer_Kokkos(N_Vector v);

SUNDIALS_EXPORT realtype *N_VGetDeviceArrayPointer_Kokkos(N_Vector v);

SUNDIALS_EXPORT void N_VCopyToDevice_Kokkos(N_Vector v);

SUNDIALS_EXPORT void N_VCopyFromDevice_Kokkos(N_Vector v);

SUNDIALS_EXPORT void N_VPrint_Kokkos(N_Vector v);

SUNDIALS_EXPORT void N_VPrintFile_Kokkos(N_Vector v, FILE *outfile);

SUNDIALS_EXPORT N_Vector_ID N_VGetVectorID_Kokkos(N_Vector v);
SUNDIALS_EXPORT N_Vector N_VCloneEmpty_Kokkos(N_Vector w);
SUNDIALS_EXPORT N_Vector N_VClone_Kokkos(N_Vector w);
SUNDIALS_EXPORT void N_VDestroy_Kokkos(N_Vector v);
SUNDIALS_EXPORT void N_VSpace_Kokkos(N_Vector v, sunindextype *lrw, sunindextype *liw);
SUNDIALS_EXPORT realtype *N_VGetArrayPointer_Kokkos(N_Vector v);
SUNDIALS_EXPORT void N_VSetArrayPointer_Kokkos(realtype *v_data, N_Vector v);

/* standard vector operations */
SUNDIALS_EXPORT void N_VLinearSum_Kokkos(realtype a, N_Vector x, realtype b, N_Vector y, N_Vector z);
SUNDIALS_EXPORT void N_VConst_Kokkos(realtype c, N_Vector z);
SUNDIALS_EXPORT void N_VProd_Kokkos(N_Vector x, N_Vector y, N_Vector z);
SUNDIALS_EXPORT void N_VDiv_Kokkos(N_Vector x, N_Vector y, N_Vector z);
SUNDIALS_EXPORT void N_VScale_Kokkos(realtype c, N_Vector x, N_Vector z);
SUNDIALS_EXPORT void N_VAbs_Kokkos(N_Vector x, N_Vector z);
SUNDIALS_EXPORT void N_VInv_Kokkos(N_Vector x, N_Vector z);
SUNDIALS_EXPORT void N_VAddConst_Kokkos(N_Vector x, realtype b, N_Vector z);
SUNDIALS_EXPORT realtype N_VDotProd_Kokkos(N_Vector x, N_Vector y);
SUNDIALS_EXPORT realtype N_VMaxNorm_Kokkos(N_Vector x);
SUNDIALS_EXPORT realtype N_VWrmsNorm_Kokkos(N_Vector x, N_Vector w);
SUNDIALS_EXPORT realtype N_VWrmsNormMask_Kokkos(N_Vector x, N_Vector w, N_Vector id);
SUNDIALS_EXPORT realtype N_VMin_Kokkos(N_Vector x);
SUNDIALS_EXPORT realtype N_VWL2Norm_Kokkos(N_Vector x, N_Vector w);
SUNDIALS_EXPORT realtype N_VL1Norm_Kokkos(N_Vector x);
SUNDIALS_EXPORT void N_VCompare_Kokkos(realtype c, N_Vector x, N_Vector z);
SUNDIALS_EXPORT booleantype N_VInvTest_Kokkos(N_Vector x, N_Vector z);
SUNDIALS_EXPORT booleantype N_VConstrMask_Kokkos(N_Vector c, N_Vector x, N_Vector m);
SUNDIALS_EXPORT realtype N_VMinQuotient_Kokkos(N_Vector num, N_Vector denom);

/* fused vector operations */
SUNDIALS_EXPORT int N_VLinearCombination_Kokkos(int nvec, realtype* c, N_Vector* X,
                                              N_Vector z);
SUNDIALS_EXPORT int N_VScaleAddMulti_Kokkos(int nvec, realtype* c, N_Vector x,
                                          N_Vector* Y, N_Vector* Z);

/* vector array operations */
SUNDIALS_EXPORT int N_VLinearSumVectorArray_Kokkos(int nvec,
                                                 realtype a, N_Vector* X,
                                                 realtype b, N_Vector* Y,
                                                 N_Vector* Z);
SUNDIALS_EXPORT int N_VScaleVectorArray_Kokkos(int nvec, realtype* c, N_Vector* X,
                                             N_Vector* Z);
SUNDIALS_EXPORT int N_VConstVectorArray_Kokkos(int nvec, realtype c, N_Vector* Z);
SUNDIALS_EXPORT int N_VScaleAddMultiVectorArray_Kokkos(int nvec, int nsum,
                                                     realtype* a,
                                                     N_Vector* X, N_Vector** Y,
                                                     N_Vector** Z);
SUNDIALS_EXPORT int N_VLinearCombinationVectorArray_Kokkos(int nvec, int nsum,
                                                         realtype* c,
                                                         N_Vector** X,
                                                         N_Vector* Z);

/* OPTIONAL local reduction kernels (no parallel communication) */
SUNDIALS_EXPORT realtype N_VWSqrSumLocal_Kokkos(N_Vector x, N_Vector w);
SUNDIALS_EXPORT realtype N_VWSqrSumMaskLocal_Kokkos(N_Vector x, N_Vector w, N_Vector id);

/*
 * -----------------------------------------------------------------
 * Enable / disable fused vector operations
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT int N_VEnableFusedOps_Kokkos(N_Vector v, booleantype tf);

SUNDIALS_EXPORT int N_VEnableLinearCombination_Kokkos(N_Vector v, booleantype tf);
SUNDIALS_EXPORT int N_VEnableScaleAddMulti_Kokkos(N_Vector v, booleantype tf);

SUNDIALS_EXPORT int N_VEnableLinearSumVectorArray_Kokkos(N_Vector v, booleantype tf);
SUNDIALS_EXPORT int N_VEnableScaleVectorArray_Kokkos(N_Vector v, booleantype tf);
SUNDIALS_EXPORT int N_VEnableConstVectorArray_Kokkos(N_Vector v, booleantype tf);
SUNDIALS_EXPORT int N_VEnableScaleAddMultiVectorArray_Kokkos(N_Vector v, booleantype tf);
SUNDIALS_EXPORT int N_VEnableLinearCombinationVectorArray_Kokkos(N_Vector v, booleantype tf);

#ifdef __cplusplus
}
#endif

#endif
