#include <stdio.h>
#include <stdlib.h>

#include <nvector/nvector_kokkos.h>
#include <Kokkos_Core.hpp>

// TODO might not need debug include
#include "sundials_debug.h"

#define ZERO   RCONST(0.0)
#define HALF   RCONST(0.5)
#define ONE    RCONST(1.0)
#define ONEPT5 RCONST(1.5)

extern "C" {

//TODO fix indenting

//TODO find best location for this
using ExecSpace = MemSpace::execution_space;
using range_policy = Kokkos::RangePolicy<ExecSpace>;

// Static constants
static constexpr sunindextype zeroIdx = 0;

// Helpful macros
#define NVEC_KOKKOS_CONTENT(x) ((N_VectorContent_Kokkos)(x->content))
#define NVEC_KOKKOS_MEMSIZE(x) (NVEC_KOKKOS_CONTENT(x)->length * sizeof(realtype))
#define NVEC_KOKKOS_PRIVATE(x) ((N_PrivateVectorContent_Kokkos)(NVEC_KOKKOS_CONTENT(x)->priv))

// TODO Find best way to implement managed memory
// Will most likely only include cuda by including CUDAUVM memspace
struct _N_PrivateVectorContent_Kokkos
{
  booleantype use_managed_mem; /* indicates if the data pointers and buffer pointers are managed memory */
};

typedef struct _N_PrivateVectorContent_Kokkos*N_PrivateVectorContent_Kokkos;

/* ----------------------------------------------------------------
 * Returns vector type ID. Used to identify vector implementation
 * from abstract N_Vector interface.
 */
N_Vector_ID N_VGetVectorID_Kokkos(N_Vector v)
{
  return SUNDIALS_NVEC_KOKKOS;
}

N_Vector N_VNewEmpty_Kokkos()
{
  N_Vector v;

  /* Create an empty vector object */
  v = NULL;
  v = N_VNewEmpty();
  if (v == NULL) return(NULL);

  /* Attach operations */

  /* constructors, destructors, and utility operations */
  v->ops->nvgetvectorid     = N_VGetVectorID_Kokkos;
  v->ops->nvclone           = N_VClone_Kokkos;
  v->ops->nvcloneempty      = N_VCloneEmpty_Kokkos;
  v->ops->nvdestroy         = N_VDestroy_Kokkos;
  v->ops->nvspace           = N_VSpace_Kokkos;
  v->ops->nvgetlength       = N_VGetLength_Kokkos;

  /* standard vector operations */
  v->ops->nvlinearsum    = N_VLinearSum_Kokkos;
  v->ops->nvconst        = N_VConst_Kokkos;
  v->ops->nvprod         = N_VProd_Kokkos;
  v->ops->nvdiv          = N_VDiv_Kokkos;
  v->ops->nvscale        = N_VScale_Kokkos;
  v->ops->nvabs          = N_VAbs_Kokkos;
  v->ops->nvinv          = N_VInv_Kokkos;
  v->ops->nvaddconst     = N_VAddConst_Kokkos;
  v->ops->nvdotprod      = N_VDotProd_Kokkos;
  v->ops->nvmaxnorm      = N_VMaxNorm_Kokkos;
  v->ops->nvmin          = N_VMin_Kokkos;
  v->ops->nvl1norm       = N_VL1Norm_Kokkos;
  v->ops->nvinvtest      = N_VInvTest_Kokkos;
  v->ops->nvconstrmask   = N_VConstrMask_Kokkos;
  v->ops->nvminquotient  = N_VMinQuotient_Kokkos;
  v->ops->nvwrmsnormmask = N_VWrmsNormMask_Kokkos;
  v->ops->nvwrmsnorm     = N_VWrmsNorm_Kokkos;
  v->ops->nvwl2norm      = N_VWL2Norm_Kokkos;
  v->ops->nvcompare      = N_VCompare_Kokkos;

  /* fused and vector array operations are disabled (NULL) by default */

  /* local reduction operations */
  v->ops->nvwsqrsumlocal     = N_VWSqrSumLocal_Kokkos;
  v->ops->nvwsqrsummasklocal = N_VWSqrSumMaskLocal_Kokkos;
  v->ops->nvdotprodlocal     = N_VDotProd_Kokkos;
  v->ops->nvmaxnormlocal     = N_VMaxNorm_Kokkos;
  v->ops->nvminlocal         = N_VMin_Kokkos;
  v->ops->nvl1normlocal      = N_VL1Norm_Kokkos;
  v->ops->nvinvtestlocal     = N_VInvTest_Kokkos;
  v->ops->nvconstrmasklocal  = N_VConstrMask_Kokkos;
  v->ops->nvminquotientlocal = N_VMinQuotient_Kokkos;

  v->content = (N_VectorContent_Kokkos) malloc(sizeof(_N_VectorContent_Kokkos));
  if (v->content == NULL)
  {
    N_VDestroy(v);
    return NULL;
  }

  N_Vector N_VNew_Kokkos(sunindextype length)
  {
    N_Vector v;

    v = NULL;
    v = N_VNewEmpty_Kokkos();
    if (v == NULL) return(NULL);

    NVEC_KOKKOS_CONTENT(v)->length          = length;
    NVEC_KOKKOS_CONTENT(v)->own_data        = SUNTRUE;
    // Can you set the object to NULL?
    NVEC_KOKKOS_CONTENT(v)->host_data       = NULL;
    NVEC_KOKKOS_CONTENT(v)->device_data     = NULL;
    NVEC_KOKKOS_PRIVATE(v)->use_managed_mem = SUNFALSE;

    if (AllocateData(v))
    {
      SUNDIALS_DEBUG_PRINT("ERROR in N_VNew_Kokkos: AllocateData returned nonzero\n");
      N_VDestroy(v);
      return NULL;
    }

    return(v);
  }

  //TODO NewManaged will be inserter here based on managed memory implementation design

  N_Vector N_VMake_Kokkos(sunindextype length, realtype *h_vdata, realtype *d_vdata)
  {
    N_Vector v;

    if (h_vdata == NULL || d_vdata == NULL) return(NULL);

    v = NULL;
    v = N_VNewEmpty_Kokkos();
    if (v == NULL) return(NULL);

    int arrSize = sizeof(h_vdata)/sizeof(realtype);

    NVEC_KOKKOS_CONTENT(v)->length          = length;
    NVEC_KOKKOS_CONTENT(v)->own_data        = SUNFALSE;
    NVEC_KOKKOS_PRIVATE(v)->use_managed_mem = SUNFALSE;
    NVEC_KOKKOS_CONTENT(v)->host_data(h_vdata, arrSize);
    NVEC_KOKKOS_CONTENT(v)->device_data(d_vdata, arrSize);

    return(v);
  }

  //TODO MakeManaged will be inserted here

  /* -----------------------------------------------------------------
   * Function to return the global length of the vector.
   */
  sunindextype N_VGetLength_Kokkos(N_Vector v)
  {
    return NVEC_KOKKOS_CONTENT(v)->length;
  }

  /* ----------------------------------------------------------------------------
   * Return pointer to the raw host data
   */

  // Todo, this currently returns actual data pointer
  // Consider adding functions to return views?
  // Same with get device array
  realtype *N_VGetHostArrayPointer_Kokkos(N_Vector x)
  {
    return NVEC_KOKKOS_CONTENT(x)->host_data.data();
  }

  /* ----------------------------------------------------------------------------
   * Return pointer to the raw device data
   */

  realtype *N_VGetDeviceArrayPointer_Kokkos(N_Vector x)
  {
    return NVEC_KOKKOS_CONTENT(x)->device_data.data();;
  }

  /* ----------------------------------------------------------------------------
   * Return a flag indicating if the memory for the vector data is managed
   */
  booleantype N_VIsManagedMemory_Kokkos(N_Vector x)
  {
    return NVEC_KOKKOS_PRIVATE(x)->use_managed_mem;
  }

  /* ----------------------------------------------------------------------------
   * Copy vector data to the device
   */

  void N_VCopyToDevice_Kokkos(N_Vector x)
  {
    if(!NVEC_KOKKOS_PRIVATE(x)->use_managed_mem)
    {
      Kokkos::deep_copy(NVEC_KOKKOS_CONTENT(x)->device_data, NVEC_KOKKOS_CONTENT(x)->host_data);
      // TODO does this error check or will I need to error check some how
    }
  }

  /* ----------------------------------------------------------------------------
   * Copy vector data from the device to the host
   */

  void N_VCopyFromDevice_Kokkos(N_Vector x)
  {
    if(!NVEC_KOKKOS_PRIVATE(x)->use_managed_mem)
    {
      Kokkos::deep_copy(NVEC_KOKKOS_CONTENT(x)->host_data, NVEC_KOKKOS_CONTENT(x)->device_data);
      // TODO does this error check or will I need to do this
    }
  }

  /* ----------------------------------------------------------------------------
   * Function to print the a serial vector to stdout
   */

  void N_VPrint_Kokkos(N_Vector X)
  {
    N_VPrintFile_Kokkos(X, stdout);
  }

  /* ----------------------------------------------------------------------------
   * Function to print the a serial vector to outfile
   */

  void N_VPrintFile_Kokkos(N_Vector X, FILE *outfile)
  {
    const realtype *xd = NVEC_KOKKOS_CONTENT(X)->host_data.data();
    const sunindextype N = NVEC_KOKKOS_CONTENT(X)->length;
    sunindextype i;

    for (i = 0; i < N; ++i) {
  #if defined(SUNDIALS_EXTENDED_PRECISION)
      fprintf(outfile, "%35.32Lg\n", xd[i]);
  #elif defined(SUNDIALS_DOUBLE_PRECISION)
      fprintf(outfile, "%19.16g\n", xd[i]);
  #else
      fprintf(outfile, "%11.8g\n", xd[i]);
  #endif
    }
    fprintf(outfile, "\n");

    return;
  }

  /*
   * -----------------------------------------------------------------
   * implementation of vector operations
   * -----------------------------------------------------------------
   */

  N_Vector N_VCloneEmpty_Kokkos(N_Vector w)
  {
    N_Vector v;

    if (w == NULL) return(NULL);

    /* Create vector */
    v = NULL;
    v = N_VNewEmpty_Kokkos();
    if (v == NULL) return(NULL);

    /* Attach operations */
    if (N_VCopyOps(w, v)) { N_VDestroy(v); return(NULL); }

    /* Set content */
    NVEC_KOKKOS_CONTENT(v)->length          = NVEC_KOKKOS_CONTENT(w)->length;
    NVEC_KOKKOS_CONTENT(v)->own_data        = SUNFALSE;
    //TODO determine best way to set views to nothing
    NVEC_KOKKOS_CONTENT(v)->host_data       = NULL;
    NVEC_KOKKOS_CONTENT(v)->device_data     = NULL;
    NVEC_KOKKOS_PRIVATE(v)->use_managed_mem = NVEC_KOKKOS_PRIVATE(w)->use_managed_mem;


    return(v);
  }

  N_Vector N_VClone_Kokkos(N_Vector w)
  {
    N_Vector v;
    v = NULL;
    v = N_VCloneEmpty_Kokkos(w);
    if (v == NULL) return(NULL);

    NVEC_KOKKOS_CONTENT(v)->own_data = SUNTRUE;
    if (AllocateData(v))
    {
      SUNDIALS_DEBUG_PRINT("ERROR in N_VClone_Kokkos: AllocateData returned nonzero\n");
      N_VDestroy(v);
      return NULL;
    }

  return(v);

  }


  void N_VDestroy_Kokkos(N_Vector v)
  {
    if (v == NULL) return;

    N_VectorContent_Kokkosvc = NVEC_KOKKOS_CONTENT(v);
    if (vc == NULL)
    {
      free(v);
      v = NULL;
      return;
    }

    N_PrivateVectorContent_Kokkosvcp = NVEC_KOKKOS_PRIVATE(v);

    if (vc->own_data)
    {
      //TODO Figure out best way to free the views
      vc->host_data= NULL;
      vc->device_data= NULL;
    }

    /* free private content struct */
    if (vcp) free(vcp);
    vc->priv = NULL;

    /* free content struct */
    free(vc);
    v->content = NULL;

    /* free ops */
    if (v->ops != NULL)
    {
      free(v->ops);
      v->ops = NULL;
    }

    /* free vector */
    free(v);
    v = NULL;

    return;
  }


  void N_VSpace_Kokkos(N_Vector X, sunindextype *lrw, sunindextype *liw)
  {
    *lrw = NVEC_KOKKOS_CONTENT(X)->length;
    *liw = 2;
  }

  void N_VConst_Kokkos(realtype c, N_Vector Z)
  {
    const sunindextype N = NVEC_KOKKOS_CONTENT(Z)->length;
    auto zdata = NVEC_KOKKOS_CONTENT(Z)->device_data;

    Kokkos::parallel_for("N_VConst", range_policy(zeroIdx, N), KOKKOS_LAMBDA(sunindextype i){
       zdata(i) = c;
    });

  }

  void N_VLinearSum_Kokkos(realtype a, N_Vector X, realtype b, N_Vector Y, N_Vector Z)
  {
    const auto xdata = NVEC_KOKKOS_CONTENT(X)->device_data;
    const auto ydata = NVEC_KOKKOS_CONTENT(Y)->device_data;
    const sunindextype N = NVEC_KOKKOS_CONTENT(X)->length;
    auto zdata = NVEC_KOKKOS_CONTENT(Z)->device_data;


    Kokkos::parallel_for("N_VLinearSum", range_policy(zeroIdx, N), KOKKOS_LAMBDA(sunindextype i){
        zdata(i) = a*xdata(i) + b*ydata(i);
      }
    );
  }

  void N_VProd_Kokkos(N_Vector X, N_Vector Y, N_Vector Z)
  {
    const auto xdata = NVEC_KOKKOS_CONTENT(X)->device_data;
    const auto ydata = NVEC_KOKKOS_CONTENT(Y)->device_data;
    const sunindextype N = NVEC_KOKKOS_CONTENT(X)->length;
    auto zdata = NVEC_KOKKOS_CONTENT(Z)->device_data;

    Kokkos::parallel_for("N_VLinearSum", range_policy(zeroIdx, N), KOKKOS_LAMBDA(sunindextype i){
        zdata(i) = xdata(i) * ydata(i);
      }
    );
  }


  void N_VDiv_Kokkos(N_Vector X, N_Vector Y, N_Vector Z)
  {
    const auto xdata = NVEC_KOKKOS_CONTENT(X)->device_data;
    const auto ydata = NVEC_KOKKOS_CONTENT(Y)->device_data;
    const sunindextype N = NVEC_KOKKOS_CONTENT(X)->length;
    auto zdata = NVEC_KOKKOS_CONTENT(Z)->device_data;

    Kokkos::parallel_for("N_VDiv", range_policy(zeroIdx, N), KOKKOS_LAMBDA(sunindextype i){
        zdata(i) = xdata(i) / ydata(i);
      }
    );
  }

  void N_VScale_Kokkos(realtype c, N_Vector X, N_Vector Z)
  {
    const auto xdata = NVEC_KOKKOS_CONTENT(X)->device_data;
    const sunindextype N = NVEC_KOKKOS_CONTENT(X)->length;
    auto zdata = NVEC_KOKKOS_CONTENT(Z)->device_data;

    Kokkos::parallel_for("N_VScale", range_policy(zeroIdx, N), KOKKOS_LAMBDA(sunindextype i){
        zdata(i) = c * xdata(i);
      }
    );
  }

  void N_VAbs_Kokkos(N_Vector X, N_Vector Z)
  {
    const auto xdata = NVEC_KOKKOS_CONTENT(X)->device_data;
    const sunindextype N = NVEC_KOKKOS_CONTENT(X)->length;
    auto zdata = NVEC_KOKKOS_CONTENT(Z)->device_data;

    Kokkos::parallel_for("N_VAbs", range_policy(zeroIdx, N), KOKKOS_LAMBDA(sunindextype i){
        zdata(i) = abs(xdata(i));
      }
    );
  }

  void N_VInv_Kokkos(N_Vector X, N_Vector Z)
  {
    const auto xdata = NVEC_KOKKOS_CONTENT(X)->device_data;
    const sunindextype N = NVEC_KOKKOS_CONTENT(X)->length;
    auto zdata = NVEC_KOKKOS_CONTENT(Z)->device_data;

    Kokkos::parallel_for("N_VInv", range_policy(zeroIdx, N), KOKKOS_LAMBDA(sunindextype i){
        zdata(i) = ONE / xdata(i);
      }
    );
  }

  void N_VAddConst_Kokkos(N_Vector X, realtype b, N_Vector Z)
  {
    const auto xdata = NVEC_KOKKOS_CONTENT(X)->device_data;
    const sunindextype N = NVEC_KOKKOS_CONTENT(X)->length;
    auto zdata = NVEC_KOKKOS_CONTENT(Z)->device_data;

    Kokkos::parallel_for("N_VAddConst", range_policy(zeroIdx, N), KOKKOS_LAMBDA(sunindextype i){
        zdata(i) = xdata(i) + b;
      }
    );
  }

  realtype N_VDotProd_Kokkos(N_Vector X, N_Vector Y)
  {
    const auto xdata = NVEC_KOKKOS_CONTENT(X)->device_data;
    const auto ydata = NVEC_KOKKOS_CONTENT(Y)->device_data;
    const sunindextype N = NVEC_KOKKOS_CONTENT(X)->length;

    realtpe gpu_result = 0.0;

    Kokkos::parallel_reduce( "N_VDotProd", range_policy(zeroIdx, N), KOKKOS_LAMBDA(sunindextype i, realtype &update) {
      update += xdata(i) * ydata(i);
    }, gpu_result );

    return (static_cast<realtype>(gpu_result));
  }

  realtype N_VMaxNorm_Kokkos(N_Vector X)
  {
    const auto xdata = NVEC_KOKKOS_CONTENT(X)->device_data;
    const sunindextype N = NVEC_KOKKOS_CONTENT(X)->length;

    realtype gpu_resul= 0.0;

    Kokkos::parallel_reduce( "N_VMaxNorm", range_policy(zeroIdx, N), KOKKOS_LAMBDA(sunindextype i, realtype &update) {
      if(abs(xdata(i)) > update) update = abs(xdata(i));
    }, Kokkos::Max<realtype>(gpu_result));

    return (static_cast<realtype>(gpu_result));
  }

  realtype N_VWSqrSumLocal_Kokkos(N_Vector X, N_Vector W)
  {
    const auto xdata = NVEC_KOKKOS_CONTENT(X)->device_data;
    const auto wdata = NVEC_KOKKOS_CONTENT(W)->device_data;
    const sunindextype N = NVEC_KOKKOS_CONTENT(X)->length;

    realtpe gpu_result = 0.0;

    Kokkos::parallel_reduce( "N_VWSqrSumLocal", range_policy(zeroIdx, N), KOKKOS_LAMBDA(sunindextype i, realtype &update) {
      update += (xdata(i) * wdata(i) * xdata(i) * wdata(i));
    }, gpu_result );

    return (static_cast<realtype>(gpu_result));
  }

  realtype N_VWrmsNorm_Kokkos(N_Vector X, N_Vector W)
  {
    const realtype sum = N_VWSqrSumLocal_Kokkos(X, W);
    const sunindextype N = NVEC_KOKKOS_CONTENT(X)->length;
    return std::sqrt(sum/N);
  }

  realtype N_VWSqrSumMaskLocal_Kokkos(N_Vector X, N_Vector W, N_Vector ID)
  {
    const auto xdata = NVEC_KOKKOS_CONTENT(X)->device_data;
    const auto wdata = NVEC_KOKKOS_CONTENT(W)->device_data;
    const auto iddata = NVEC_KOKKOS_CONTENT(ID)->device_data;
    const sunindextype N = NVEC_KOKKOS_CONTENT(X)->length;

    realtpe gpu_result = 0.0;

    Kokkos::parallel_reduce( "N_VWSqrSumMaskLocal", range_policy(zeroIdx, N), KOKKOS_LAMBDA(sunindextype i, realtype &update) {
      if (iddata(i) > ZERO)
        update += (xdata(i) * wdata(i) * xdata(i) * wdata(i));
    }, gpu_result );

    return (static_cast<realtype>(gpu_result));
  }

  realtype N_VWrmsNormMask_Kokkos(N_Vector X, N_Vector W, N_Vector ID)
  {
    const realtype sum = N_VWSqrSumMaskLocal_Kokkos(X, W, ID);
    const sunindextype N = NVEC_KOKKOS_CONTENT(X)->length;
    return std::sqrt(sum/N);
  }

  realtype N_VMin_Kokkos(N_Vector X)
  {
    const auto xdata = NVEC_KOKKOS_CONTENT(X)->device_data;
    const sunindextype N = NVEC_KOKKOS_CONTENT(X)->length;

    realtype gpu_resul = std::numeric_limits<realtype>::max();

    Kokkos::parallel_reduce( "N_VMin", range_policy(zeroIdx, N), KOKKOS_LAMBDA(sunindextype i, realtype &update) {
      if(xdata(i) < update) update = xdata(i);
    }, Kokkos::Min<realtype>(gpu_result));

    return (static_cast<realtype>(gpu_result));

  }

  realtype N_VWL2Norm_Kokkos(N_Vector X, N_Vector W)
  {
    return std::sqrt(N_VWSqrSumLocal_Kokkos(X, W));
  }

  realtype N_VL1Norm_Kokkos(N_Vector X)
  {
    const auto xdata = NVEC_KOKKOS_CONTENT(X)->device_data;
    const sunindextype N = NVEC_KOKKOS_CONTENT(X)->length;

    realtpe gpu_result = 0.0;

    Kokkos::parallel_reduce( "N_VL1Norm", range_policy(zeroIdx, N), KOKKOS_LAMBDA(sunindextype i, realtype &update) {
        update += += (abs(xdata[i]));
    }, gpu_result );

    return (static_cast<realtype>(gpu_result));
  }

  void N_VCompare_Kokkos(realtype c, N_Vector X, N_Vector Z)
  {
    const auto xdata = NVEC_KOKKOS_CONTENT(X)->device_data;
    const sunindextype N = NVEC_KOKKOS_CONTENT(X)->length;
    auto zdata = NVEC_KOKKOS_CONTENT(Z)->device_data;

    Kokkos::parallel_for("N_VCompare", range_policy(zeroIdx, N), KOKKOS_LAMBDA(sunindextype i){
      zdata(i) = abs(xdata(i)) >= c ? ONE : ZERO;
    });

  }

  booleantype N_VInvTest_Kokkos(N_Vector x, N_Vector z)
  {
    const auto xdata = NVEC_KOKKOS_CONTENT(x)->device_data;
    const sunindextype N = NVEC_KOKKOS_CONTENT(x)->length;
    realtype *zdata = NVEC_KOKKOS_CONTENT(z)->device_data;

    realtpe gpu_result = 0.0;

    Kokkos::parallel_reduce( "N_VInvTest", range_policy(zeroIdx, N), KOKKOS_LAMBDA(sunindextype i, realtype &update) {
      if (xdata(i) == ZERO) {
        update += ONE;
      } else {
        zdata(i) = ONE/xdata(i);
      }
    }, gpu_result);

    realtype minimum = static_cast<realtype>(gpu_result);
    return (minimum < HALF);
  }

  booleantype N_VConstrMask_Kokkos(N_Vector c, N_Vector x, N_Vector m)
  {
    const auto cdata = NVEC_KOKKOS_CONTENT(c)->device_data;
    const auto xdata = NVEC_KOKKOS_CONTENT(x)->device_data;
    const sunindextype N = NVEC_KOKKOS_CONTENT(x)->length;
    auto mdata = NVEC_KOKKOS_CONTENT(m)->device_data;

    realtpe gpu_result = 0.0;

    Kokkos::parallel_reduce( "N_VConstrMask", range_policy(zeroIdx, N), KOKKOS_LAMBDA(sunindextype i, realtype &update) {
      KOKKOS_LAMBDA(sunindextype i) {
        bool test = (abs(cdata(i)) > ONEPT5 && cdata(i)*xdata(i) <= ZERO) ||
                    (abs(cdata(i)) > HALF   && cdata(i)*xdata(i) <  ZERO);
        mdata(i) = test ? ONE : ZERO;
        update += mdata(i);
      }
    }, gpu_result);

    realtype sum = static_cast<realtype>(gpu_result);
    return(sum < HALF);
  }


  realtype N_VMinQuotient_Kokkos(N_Vector num, N_Vector denom)
  {
    const auto ndata = NVEC_KOKKOS_CONTENT(num)->device_data;
    const auto ddata = NVEC_KOKKOS_CONTENT(denom)->device_data;
    const sunindextype N = NVEC_KOKKOS_CONTENT(num)->length;

    Kokkos::parallel_reduce( "N_VMinQuotient", range_policy(zeroIdx, N), KOKKOS_LAMBDA(sunindextype i, realtype &update) {
      if (ddata(i) != ZERO){
        if((ndata(i)/ddata(i)) < update) update = ndata(i)/ddata(i);
      }
    }, Kokkos::Min<realtype>(gpu_result));

    KOKKOS::ReduceMin< KOKKOS_REDUCE_TYPE, realtype> gpu_result(std::numeric_limits<realtype>::max());
    KOKKOS::forall< KOKKOS_NODE_TYPE >(KOKKOS::RangeSegment(zeroIdx, N),
      KOKKOS_LAMBDA(sunindextype i) {
        if (ddata[i] != ZERO)
          gpu_result.min(ndata[i]/ddata[i]);
      }
    );
    return (static_cast<realtype>(gpu_result));
  }
