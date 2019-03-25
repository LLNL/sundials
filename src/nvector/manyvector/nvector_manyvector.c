/* -----------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ SMU
 * -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2019, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------
 * This is the implementation file for the ManyVector implementation
 * of the NVECTOR package.
 * -----------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <nvector/nvector_manyvector.h>
#include <sundials/sundials_math.h>

#define ZERO RCONST(0.0)
#define ONE  RCONST(1.0)

/* -----------------------------------------------------------------
   ManyVector content accessor macros
   -----------------------------------------------------------------*/
#define MANYVECTOR_CONTENT(v)     ( (N_VectorContent_ManyVector)(v->content) )
#define MANYVECTOR_COMM(v)        ( MANYVECTOR_CONTENT(v)->comm )
#define MANYVECTOR_NUM_SUBVECS(v) ( MANYVECTOR_CONTENT(v)->num_subvectors )
#define MANYVECTOR_GLOBLENGTH(v)  ( MANYVECTOR_CONTENT(v)->global_length )
#define MANYVECTOR_SUBVECS(v)     ( MANYVECTOR_CONTENT(v)->subvec_array )
#define MANYVECTOR_SUBVEC(v,i)    ( MANYVECTOR_SUBVECS(v)[i] )
#define MANYVECTOR_OWN_DATA(v)    ( MANYVECTOR_CONTENT(v)->own_data )

/* -----------------------------------------------------------------
   Prototypes of utility routines
   -----------------------------------------------------------------*/
static N_Vector ManyVectorClone(N_Vector w, booleantype cloneempty);
static int SubvectorMPIRank(N_Vector w);


/* -----------------------------------------------------------------
   ManyVector API routines
   -----------------------------------------------------------------*/

/* This function creates a ManyVector from a set of existing
   N_Vector objects, along with a user-created MPI (inter/intra)communicator
   that can be used to couple all subvectors together. */
N_Vector N_VMake_ManyVector(void *comm, sunindextype num_subvectors,
                            N_Vector *vec_array)
{
  N_Vector v;
  N_Vector_Ops ops;
  N_VectorContent_ManyVector content;
  sunindextype i, local_length;
  realtype globlen;
  SUNMPI_Comm *vcomm;
  int rank, retval;

  /* Check that input N_Vectors are non-NULL */
  if (vec_array == NULL)  return(NULL);
  for (i=0; i<num_subvectors; i++)
    if (vec_array[i] == NULL)  return(NULL);

  /* Create vector */
  v = NULL;
  v = (N_Vector) malloc(sizeof *v);
  if (v == NULL) return(NULL);

  /* Create vector operation structure */
  ops = NULL;
  ops = (N_Vector_Ops) malloc(sizeof(struct _generic_N_Vector_Ops));
  if (ops == NULL) { free(v); return(NULL); }

  /* attach standard vector operations */
  ops->nvgetvectorid     = N_VGetVectorID_ManyVector;
  ops->nvcloneempty      = N_VCloneEmpty_ManyVector;
  ops->nvclone           = N_VClone_ManyVector;
  ops->nvdestroy         = N_VDestroy_ManyVector;
  ops->nvspace           = N_VSpace_ManyVector;
  ops->nvgetarraypointer = NULL;
  ops->nvsetarraypointer = NULL;
  ops->nvgetcommunicator = N_VGetCommunicator_ManyVector;
  ops->nvgetlength       = N_VGetLength_ManyVector;
  ops->nvlinearsum       = N_VLinearSum_ManyVector;
  ops->nvconst           = N_VConst_ManyVector;
  ops->nvprod            = N_VProd_ManyVector;
  ops->nvdiv             = N_VDiv_ManyVector;
  ops->nvscale           = N_VScale_ManyVector;
  ops->nvabs             = N_VAbs_ManyVector;
  ops->nvinv             = N_VInv_ManyVector;
  ops->nvaddconst        = N_VAddConst_ManyVector;
  ops->nvdotprod         = N_VDotProd_ManyVector;
  ops->nvmaxnorm         = N_VMaxNorm_ManyVector;
  ops->nvwrmsnorm        = N_VWrmsNorm_ManyVector;
  ops->nvwrmsnormmask    = N_VWrmsNormMask_ManyVector;
  ops->nvmin             = N_VMin_ManyVector;
  ops->nvwl2norm         = N_VWL2Norm_ManyVector;
  ops->nvl1norm          = N_VL1Norm_ManyVector;
  ops->nvcompare         = N_VCompare_ManyVector;
  ops->nvinvtest         = N_VInvTest_ManyVector;
  ops->nvconstrmask      = N_VConstrMask_ManyVector;
  ops->nvminquotient     = N_VMinQuotient_ManyVector;

  /* fused vector operations */
  ops->nvlinearcombination = N_VLinearCombination_ManyVector;
  ops->nvscaleaddmulti     = N_VScaleAddMulti_ManyVector;
  ops->nvdotprodmulti      = N_VDotProdMulti_ManyVector;

  /* vector array operations (optional, NULL means disabled by default) */
  ops->nvlinearsumvectorarray         = NULL;
  ops->nvscalevectorarray             = NULL;
  ops->nvconstvectorarray             = NULL;
  ops->nvscaleaddmultivectorarray     = NULL;
  ops->nvlinearcombinationvectorarray = NULL;
  ops->nvwrmsnormvectorarray          = N_VWrmsNormVectorArray_ManyVector;
  ops->nvwrmsnormmaskvectorarray      = N_VWrmsNormMaskVectorArray_ManyVector;

  /* local reduction kernels */
  ops->nvdotprodlocal     = N_VDotProdLocal_ManyVector;
  ops->nvmaxnormlocal     = N_VMaxNormLocal_ManyVector;
  ops->nvminlocal         = N_VMinLocal_ManyVector;
  ops->nvl1normlocal      = N_VL1NormLocal_ManyVector;
  ops->nvinvtestlocal     = N_VInvTestLocal_ManyVector;
  ops->nvconstrmasklocal  = N_VConstrMaskLocal_ManyVector;
  ops->nvminquotientlocal = N_VMinQuotientLocal_ManyVector;
  ops->nvwsqrsumlocal     = N_VWSqrSumLocal_ManyVector;
  ops->nvwsqrsummasklocal = N_VWSqrSumMaskLocal_ManyVector;

  /* Create content */
  content = NULL;
  content = (N_VectorContent_ManyVector) malloc(sizeof(struct _N_VectorContent_ManyVector));
  if (content == NULL) { free(ops); free(v); return(NULL); }


  /* Attach content components */

  /*    duplicate input communicator (if non-NULL) */
  content->comm = SUNMPI_COMM_NULL;
  if (comm != NULL) {
    vcomm = (SUNMPI_Comm *) comm;
    if (*vcomm != SUNMPI_COMM_NULL) {
      retval = SUNMPI_Comm_dup(*vcomm, &(content->comm));
      if (retval != SUNMPI_SUCCESS) { free(content); free(ops); free(v); return(NULL); }
    }
  }

  /*    allocate and set subvector array */
  content->num_subvectors = num_subvectors;
  content->own_data       = SUNFALSE;
  content->subvec_array   = NULL;

  content->subvec_array = (N_Vector *) malloc( content->num_subvectors * sizeof(N_Vector) );
  if (content->subvec_array == NULL) {
    free(ops); free(content); free(v); return(NULL); }

  for (i=0; i<content->num_subvectors; i++)
    content->subvec_array[i] = vec_array[i];

  /* Determine overall ManyVector length: sum contributions from all
     subvectors where this rank is the root, then perform reduction */
  local_length = 0;
  for (i=0; i<content->num_subvectors; i++) {

    /* determine rank of this task in subvector communicator
       (serial vectors default to rank=0) */
    rank = SubvectorMPIRank(vec_array[i]);
    if (rank < 0) {
      free(content); free(ops); free(v);
      free(content->subvec_array); return(NULL); }

    /* accumulate contribution from root tasks */
    if (rank == 0)  local_length += N_VGetLength(vec_array[i]);
  }
  retval = SUNMPI_Allreduce_scalar(local_length, &globlen, SUNMPI_SUM, content->comm);
  if (retval != SUNMPI_SUCCESS) {
    free(content); free(ops); free(v);
    free(content->subvec_array); return(NULL);
  }
  content->global_length = (sunindextype) globlen;

  /* Attach content and ops to new vector, and return */
  v->content = content;
  v->ops     = ops;

  return(v);
}


/* This function creates a ManyVector from a set of existing
   N_Vector objects, under the requirement that all MPI-aware
   sub-vectors use the same MPI communicator (this is verified
   internally).  If no sub-vector is MPI-aware, then this may be
   used to describe data partitioning within a single node, and
   a NULL communicator will be created. */
N_Vector N_VNew_ManyVector(sunindextype num_subvectors,
                           N_Vector *vec_array)
{
  sunindextype i;
  booleantype nocommfound;
  void* tmpcomm;
  SUNMPI_Comm comm, *vcomm;
  int retval, comparison;

  /* Check that all subvectors have identical MPI communicators (if present) */
  nocommfound = SUNTRUE;
  comm = SUNMPI_COMM_NULL;
  for (i=0; i<num_subvectors; i++) {

    /* access MPI communicator for subvector i (vcomm);
       if none is present then continue to next subvector */
    tmpcomm = N_VGetCommunicator(vec_array[i]);
    if (tmpcomm == NULL)  continue;
    vcomm = (SUNMPI_Comm *) tmpcomm;

    /* if this is the first communicator, create a copy */
    if (nocommfound) {

      /* set comm to duplicate this first subvector communicator */
      retval = SUNMPI_Comm_dup(*vcomm, &comm);
      if (retval != SUNMPI_SUCCESS)  return(NULL);
      nocommfound = SUNFALSE;

    /* otherwise, verify that vcomm matches stored comm */
    } else {

      retval = SUNMPI_Comm_compare(*vcomm, comm, &comparison);
      if ((comparison != SUNMPI_IDENT) && (comparison != SUNMPI_CONGRUENT))
        return(NULL);

    }
  }

  /* Create vector using "Make" routine and shared communicator (if non-NULL) */
  if (comm != SUNMPI_COMM_NULL)
    return(N_VMake_ManyVector(&comm, num_subvectors, vec_array));
  else
    return(N_VMake_ManyVector(NULL, num_subvectors, vec_array));
}


/* This function returns the vec_num sub-N_Vector from the N_Vector
   array.  If vec_num is outside of applicable bounds, NULL is returned. */
N_Vector N_VGetSubvector_ManyVector(N_Vector v, sunindextype vec_num)
{
  if ( (vec_num < 0) || (vec_num > MANYVECTOR_NUM_SUBVECS(v)) )
    return(NULL);
  return(MANYVECTOR_SUBVEC(v,vec_num));
}


/* This function returns the overall number of sub-vectors.
   It returns a locally stored integer, and is therefore a local call. */
sunindextype N_VGetNumSubvectors_ManyVector(N_Vector v)
{
  return(MANYVECTOR_NUM_SUBVECS(v));
}


/* -----------------------------------------------------------------
   ManyVector implementations of generic NVector routines
   -----------------------------------------------------------------*/

/* Returns vector type ID. Used to identify vector implementation
   from abstract N_Vector interface. */
N_Vector_ID N_VGetVectorID_ManyVector(N_Vector v)
{
  return(SUNDIALS_NVEC_MANYVECTOR);
}


/* Clones a ManyVector, calling CloneEmpty on subvectors. */
N_Vector N_VCloneEmpty_ManyVector(N_Vector w)
{
  return(ManyVectorClone(w, SUNTRUE));
}


/* Clones a ManyVector, calling Clone on subvectors. */
N_Vector N_VClone_ManyVector(N_Vector w)
{
  return(ManyVectorClone(w, SUNFALSE));
}


/* Destroys a ManyVector */
void N_VDestroy_ManyVector(N_Vector v)
{
  sunindextype i;

  /* if subvectors are owned by v, then Destroy those */
  if ((MANYVECTOR_OWN_DATA(v) == SUNTRUE) && (MANYVECTOR_SUBVECS(v) != NULL)) {
    for (i=0; i<MANYVECTOR_NUM_SUBVECS(v); i++)
      N_VDestroy(MANYVECTOR_SUBVEC(v,i));
  }

  /* free subvector array */
  free(MANYVECTOR_SUBVECS(v)); MANYVECTOR_SUBVECS(v) = NULL;

  /* free N_Vector structures */
  free(v->content); v->content = NULL;
  free(v->ops); v->ops = NULL;
  free(v); v = NULL;

  return;
}


/* Returns the space requirements for the ManyVector, by accumulating this
   information from all subvectors. */
void N_VSpace_ManyVector(N_Vector v, sunindextype *lrw, sunindextype *liw)
{
  sunindextype i, lrw1, liw1;
  *lrw = 0;
  *liw = 0;
  for (i=0; i<MANYVECTOR_NUM_SUBVECS(v); i++) {
    /* update space requirements for this subvector (if 'nvspace' is implemented) */
    if ((MANYVECTOR_SUBVEC(v,i))->ops->nvspace != NULL) {
      N_VSpace(MANYVECTOR_SUBVEC(v,i), &lrw1, &liw1);
      *lrw += lrw1;
      *liw += liw1;
    }
  }
  return;
}


/* This function retrieves the MPI Communicator from a ManyVector object. */
void *N_VGetCommunicator_ManyVector(N_Vector v)
{
  if (MANYVECTOR_COMM(v) == SUNMPI_COMM_NULL)
    return NULL;
  else
    return((void *) &MANYVECTOR_COMM(v));
}


/* This function retrieves the global length of a ManyVector object. */
sunindextype N_VGetLength_ManyVector(N_Vector v)
{
  return(MANYVECTOR_GLOBLENGTH(v));
}


/* Performs the linear sum z = a*x + b*y by calling N_VLinearSum on all subvectors;
   this routine does not check that x, y and z are ManyVectors, if they have the
   same number of subvectors, or if these subvectors are compatible. */
void N_VLinearSum_ManyVector(realtype a, N_Vector x, realtype b, N_Vector y, N_Vector z)
{
  sunindextype i;
  for (i=0; i<MANYVECTOR_NUM_SUBVECS(x); i++)
    N_VLinearSum(a, MANYVECTOR_SUBVEC(x,i), b, MANYVECTOR_SUBVEC(y,i),
                 MANYVECTOR_SUBVEC(z,i));
  return;
}


/* Performs the operation z = c by calling N_VConst on all subvectors. */
void N_VConst_ManyVector(realtype c, N_Vector z)
{
  sunindextype i;
  for (i=0; i<MANYVECTOR_NUM_SUBVECS(z); i++)
    N_VConst(c, MANYVECTOR_SUBVEC(z,i));
  return;
}


/* Performs the operation z_j = x_j*y_j by calling N_VProd on all subvectors;
   this routine does not check that x, y and z are ManyVectors, if they have the
   same number of subvectors, or if these subvectors are compatible. */
void N_VProd_ManyVector(N_Vector x, N_Vector y, N_Vector z)
{
  sunindextype i;
  for (i=0; i<MANYVECTOR_NUM_SUBVECS(x); i++)
    N_VProd(MANYVECTOR_SUBVEC(x,i), MANYVECTOR_SUBVEC(y,i),
            MANYVECTOR_SUBVEC(z,i));
  return;
}


/* Performs the operation z_j = x_j/y_j by calling N_VDiv on all subvectors;
   this routine does not check that x, y and z are ManyVectors, if they have the
   same number of subvectors, or if these subvectors are compatible. */
void N_VDiv_ManyVector(N_Vector x, N_Vector y, N_Vector z)
{
  sunindextype i;
  for (i=0; i<MANYVECTOR_NUM_SUBVECS(x); i++)
    N_VDiv(MANYVECTOR_SUBVEC(x,i), MANYVECTOR_SUBVEC(y,i),
           MANYVECTOR_SUBVEC(z,i));
  return;
}


/* Performs the operation z_j = c*x_j by calling N_VScale on all subvectors;
   this routine does not check that x and z are ManyVectors, if they have the
   same number of subvectors, or if these subvectors are compatible. */
void N_VScale_ManyVector(realtype c, N_Vector x, N_Vector z)
{
  sunindextype i;
  for (i=0; i<MANYVECTOR_NUM_SUBVECS(x); i++)
    N_VScale(c, MANYVECTOR_SUBVEC(x,i), MANYVECTOR_SUBVEC(z,i));
  return;
}


/* Performs the operation z_j = |x_j| by calling N_VAbs on all subvectors;
   this routine does not check that x and z are ManyVectors, if they have the
   same number of subvectors, or if these subvectors are compatible. */
void N_VAbs_ManyVector(N_Vector x, N_Vector z)
{
  sunindextype i;
  for (i=0; i<MANYVECTOR_NUM_SUBVECS(x); i++)
    N_VAbs(MANYVECTOR_SUBVEC(x,i), MANYVECTOR_SUBVEC(z,i));
  return;
}


/* Performs the operation z_j = 1/x_j by calling N_VInv on all subvectors;
   this routine does not check that x and z are ManyVectors, if they have the
   same number of subvectors, or if these subvectors are compatible. */
void N_VInv_ManyVector(N_Vector x, N_Vector z)
{
  sunindextype i;
  for (i=0; i<MANYVECTOR_NUM_SUBVECS(x); i++)
    N_VInv(MANYVECTOR_SUBVEC(x,i), MANYVECTOR_SUBVEC(z,i));
  return;
}


/* Performs the operation z_j = x_j + b by calling N_VAddConst on all subvectors;
   this routine does not check that x and z are ManyVectors, if they have the
   same number of subvectors, or if these subvectors are compatible. */
void N_VAddConst_ManyVector(N_Vector x, realtype b, N_Vector z)
{
  sunindextype i;
  for (i=0; i<MANYVECTOR_NUM_SUBVECS(x); i++)
    N_VAddConst(MANYVECTOR_SUBVEC(x,i), b, MANYVECTOR_SUBVEC(z,i));
  return;
}


/* Performs the MPI task-local dot product of two ManyVectors by calling
   N_VDotProdLocal on all subvectors; this routine does not check that x and
   y are ManyVectors, if they have the same number of subvectors, or if these
   subvectors are compatible.

   If any subvector does not implement the N_VDotProdLocal routine (NULL
   function pointer), then this routine will call N_VDotProd, but only
   accumulate the sum if this is the root task for that subvector's
   communicator (note: serial vectors are always root task). */
realtype N_VDotProdLocal_ManyVector(N_Vector x, N_Vector y)
{
  sunindextype i;
  realtype sum, contrib;
  int rank;

  /* initialize output*/
  sum = ZERO;

  for (i=0; i<MANYVECTOR_NUM_SUBVECS(x); i++) {

    /* check for nvdotprodlocal in subvector */
    if (MANYVECTOR_SUBVEC(x,i)->ops->nvdotprodlocal) {

      sum += N_VDotProdLocal(MANYVECTOR_SUBVEC(x,i), MANYVECTOR_SUBVEC(y,i));

    /* otherwise, call nvdotprod and root tasks accumulate to overall sum */
    } else {

      contrib = N_VDotProd(MANYVECTOR_SUBVEC(x,i), MANYVECTOR_SUBVEC(y,i));

      /* get this task's rank in subvector communicator (note: serial
         subvectors will result in rank==0) */
      rank = SubvectorMPIRank(MANYVECTOR_SUBVEC(x,i));
      if (rank < 0)   return(ZERO);
      if (rank == 0)  sum += contrib;

    }
  }

  return(sum);
}


/* Performs the dot product of two ManyVectors by calling N_VDotProdLocal and
   combining the results.  This routine does not check that x and y are
   ManyVectors, if they have the same number of subvectors, or if these
   subvectors are compatible. */
realtype N_VDotProd_ManyVector(N_Vector x, N_Vector y)
{
  realtype gsum;
  SUNMPI_Allreduce_scalar(N_VDotProdLocal_ManyVector(x,y),
                          &gsum, SUNMPI_SUM, MANYVECTOR_COMM(x));
  return(gsum);
}


/* Performs the MPI task-local maximum norm of a ManyVector by calling
   N_VMaxNormLocal on all subvectors.

   If any subvector does not implement the N_VMaxNormLocal routine (NULL
   function pointer), then this routine will call N_VMaxNorm instead. */
realtype N_VMaxNormLocal_ManyVector(N_Vector x)
{
  sunindextype i;
  realtype max, lmax;

  /* initialize output*/
  max = ZERO;

  for (i=0; i<MANYVECTOR_NUM_SUBVECS(x); i++) {

    /* check for nvmaxnormlocal in subvector */
    if (MANYVECTOR_SUBVEC(x,i)->ops->nvmaxnormlocal) {

      lmax = N_VMaxNormLocal(MANYVECTOR_SUBVEC(x,i));
      max = (max > lmax) ? max : lmax;

    /* otherwise, call nvmaxnorm and accumulate to overall max */
    } else {

      lmax = N_VMaxNorm(MANYVECTOR_SUBVEC(x,i));
      max = (max > lmax) ? max : lmax;

    }
  }

  return(max);
}


/* Performs the maximum norm of a ManyVector by calling N_VMaxNormLocal and
   combining the results. */
realtype N_VMaxNorm_ManyVector(N_Vector x)
{
  realtype gmax;
  SUNMPI_Allreduce_scalar(N_VMaxNormLocal_ManyVector(x),
                          &gmax, SUNMPI_MAX, MANYVECTOR_COMM(x));
  return(gmax);
}


/* Performs the MPI task-local weighted squared sum of a ManyVector by calling
   N_VWSqrSumLocal on all subvectors; this routine does not check that x and
   w are ManyVectors, if they have the same number of subvectors, or if these
   subvectors are compatible.

   If any subvector does not implement the N_VWSqrSumLocal routine (NULL
   function pointer), then this routine will call N_VWrmsNorm and N_VGetLength
   to unravel the squared sum of the subvector components.  It will then only
   accumulate this to the overall sum if this is the root task for that
   subvector's communicator (note: serial vectors are always root task). */
realtype N_VWSqrSumLocal_ManyVector(N_Vector x, N_Vector w)
{
  sunindextype i, N;
  realtype sum, contrib;
  int rank;

  /* initialize output*/
  sum = ZERO;

  for (i=0; i<MANYVECTOR_NUM_SUBVECS(x); i++) {

    /* check for nvwsqrsumlocal in subvector */
    if (MANYVECTOR_SUBVEC(x,i)->ops->nvwsqrsumlocal) {

      sum += N_VWSqrSumLocal(MANYVECTOR_SUBVEC(x,i), MANYVECTOR_SUBVEC(w,i));

    /* otherwise, call nvwrmsnorm, and accumulate to overall sum on root task */
    } else {

      contrib = N_VWrmsNorm(MANYVECTOR_SUBVEC(x,i), MANYVECTOR_SUBVEC(w,i));

      /* get this task's rank in subvector communicator (note: serial
         subvectors will result in rank==0) */
      rank = SubvectorMPIRank(MANYVECTOR_SUBVEC(x,i));
      if (rank < 0)  return(ZERO);
      if (rank == 0) {
        N = N_VGetLength(MANYVECTOR_SUBVEC(x,i));
        sum += (contrib*contrib*N);
      }
    }
  }

  return(sum);
}


/* Performs the WRMS norm of a ManyVector by calling N_VWSqrSumLocal and
   combining the results; this routine does not check that x and
   w are ManyVectors, if they have the same number of subvectors, or if these
   subvectors are compatible. */
realtype N_VWrmsNorm_ManyVector(N_Vector x, N_Vector w)
{
  realtype gsum;
  SUNMPI_Allreduce_scalar(N_VWSqrSumLocal_ManyVector(x, w),
                          &gsum, SUNMPI_SUM, MANYVECTOR_COMM(x));
  return(SUNRsqrt(gsum/(MANYVECTOR_GLOBLENGTH(x))));
}


/* Performs the MPI task-local masked weighted squared sum of a ManyVector by
   calling N_VWSqrSumMaskLocal on all subvectors; this routine does not check
   that x, w and id are ManyVectors, if they have the same number of
   subvectors, or if these subvectors are compatible.

   If any subvector does not implement the N_VWSqrSumLocal routine (NULL
   function pointer), then this routine will call N_VWrmsNormMask and
   N_VGetLength to unravel the masked squared sum of the subvector components.
   It will then only accumulate this to the overall sum if this is the root
   task for that subvector's communicator (note: serial vectors are always
   root task). */
realtype N_VWSqrSumMaskLocal_ManyVector(N_Vector x, N_Vector w, N_Vector id)
{
  sunindextype i, N;
  realtype sum, contrib;
  int rank;

  /* initialize output*/
  sum = ZERO;

  for (i=0; i<MANYVECTOR_NUM_SUBVECS(x); i++) {

    /* check for nvwsqrsummasklocal in subvector */
    if (MANYVECTOR_SUBVEC(x,i)->ops->nvwsqrsummasklocal) {

      sum += N_VWSqrSumMaskLocal(MANYVECTOR_SUBVEC(x,i), MANYVECTOR_SUBVEC(w,i),
                                 MANYVECTOR_SUBVEC(id,i));

    /* otherwise, call nvwrmsnormmask, and accumulate to overall sum on root task */
    } else {

      contrib = N_VWrmsNormMask(MANYVECTOR_SUBVEC(x,i), MANYVECTOR_SUBVEC(w,i),
                                MANYVECTOR_SUBVEC(id,i));

      /* get this task's rank in subvector communicator (note: serial
         subvectors will result in rank==0) */
      rank = SubvectorMPIRank(MANYVECTOR_SUBVEC(x,i));
      if (rank < 0)  return(ZERO);
      if (rank == 0) {
        N = N_VGetLength(MANYVECTOR_SUBVEC(x,i));
        sum += (contrib*contrib*N);
      }
    }
  }

  return(sum);
}


/* Performs the masked WRMS norm of a ManyVector by calling N_VWSqrSumMaskLocal
   and combining the results; this routine does not check that x, w and id are
   ManyVectors, if they have the same number of subvectors, or if these subvectors
   are compatible. */
realtype N_VWrmsNormMask_ManyVector(N_Vector x, N_Vector w, N_Vector id)
{
  realtype gsum;
  SUNMPI_Allreduce_scalar(N_VWSqrSumMaskLocal_ManyVector(x, w, id),
                          &gsum, SUNMPI_SUM, MANYVECTOR_COMM(x));
  return(SUNRsqrt(gsum/(MANYVECTOR_GLOBLENGTH(x))));
}


/* Computes the MPI task-local minimum entry of a ManyVector by calling
   N_VMinLocal on all subvectors.

   If any subvector does not implement the N_VMinLocal routine (NULL
   function pointer), then this routine will call N_VMin instead. */
realtype N_VMinLocal_ManyVector(N_Vector x)
{
  sunindextype i;
  realtype min, lmin;

  /* initialize output*/
  min = BIG_REAL;

  for (i=0; i<MANYVECTOR_NUM_SUBVECS(x); i++) {

    /* check for nvminlocal in subvector */
    if (MANYVECTOR_SUBVEC(x,i)->ops->nvminlocal) {

      lmin = N_VMinLocal(MANYVECTOR_SUBVEC(x,i));
      min = (min < lmin) ? min : lmin;

    /* otherwise, call nvmin and accumulate to overall min */
    } else {

      lmin = N_VMin(MANYVECTOR_SUBVEC(x,i));
      min = (min < lmin) ? min : lmin;

    }
  }

  return(min);
}


/* Computes the minimum entry of a ManyVector by calling N_VMinLocal and
   combining the results. */
realtype N_VMin_ManyVector(N_Vector x)
{
  realtype gmin;
  SUNMPI_Allreduce_scalar(N_VMinLocal_ManyVector(x),
                          &gmin, SUNMPI_MIN, MANYVECTOR_COMM(x));
  return(gmin);
}


/* Performs the WL2 norm of a ManyVector by calling N_VSqrSumLocal on all subvectors
   and combining the results.  This routine does not check that x and w are
   ManyVectors, if they have the same number of subvectors, or if these subvectors
   are compatible. */
realtype N_VWL2Norm_ManyVector(N_Vector x, N_Vector w)
{
  realtype gsum;
  SUNMPI_Allreduce_scalar(N_VWSqrSumLocal_ManyVector(x, w),
                          &gsum, SUNMPI_SUM, MANYVECTOR_COMM(x));
  return(SUNRsqrt(gsum));
}


/* Performs the MPI task-local L1 norm of a ManyVector by calling N_VL1NormLocal on
   all subvectors.  If any subvector does not implement the N_VL1NormLocal routine
   (NULL function pointer), then this routine will call N_VL1Norm, but only
   accumulate the sum if this is the root task for that subvector's
   communicator (note: serial vectors are always root task). */
realtype N_VL1NormLocal_ManyVector(N_Vector x)
{
  sunindextype i;
  realtype sum, contrib;
  int rank;

  /* initialize output*/
  sum = ZERO;

  for (i=0; i<MANYVECTOR_NUM_SUBVECS(x); i++) {

    /* check for nvl1normlocal in subvector */
    if (MANYVECTOR_SUBVEC(x,i)->ops->nvl1normlocal) {

      sum += N_VL1NormLocal(MANYVECTOR_SUBVEC(x,i));

    /* otherwise, call nvl1norm and root tasks accumulate to overall sum */
    } else {

      contrib = N_VL1Norm(MANYVECTOR_SUBVEC(x,i));

      /* get this task's rank in subvector communicator (note: serial
         subvectors will result in rank==0) */
      rank = SubvectorMPIRank(MANYVECTOR_SUBVEC(x,i));
      if (rank < 0)   return(ZERO);
      if (rank == 0)  sum += contrib;

    }
  }

  return(sum);
}


/* Performs the L1 norm of a ManyVector by calling N_VL1NormLocal and
   combining the results. */
realtype N_VL1Norm_ManyVector(N_Vector x)
{
  realtype gsum;
  SUNMPI_Allreduce_scalar(N_VL1NormLocal_ManyVector(x),
                          &gsum, SUNMPI_SUM, MANYVECTOR_COMM(x));
  return(gsum);
}


/* Performs N_VCompare on all subvectors; this routine does not check that x and z are
   ManyVectors, if they have the same number of subvectors, or if these subvectors are
   compatible. */
void N_VCompare_ManyVector(realtype c, N_Vector x, N_Vector z)
{
  sunindextype i;
  for (i=0; i<MANYVECTOR_NUM_SUBVECS(x); i++)
    N_VCompare(c, MANYVECTOR_SUBVEC(x,i), MANYVECTOR_SUBVEC(z,i));
  return;
}


/* Performs the MPI task-local InvTest for a ManyVector by calling N_VInvTestLocal
   on all subvectors and combining the results appropriately.  This routine does
   not check that x and z are ManyVectors, if they have the same number of
   subvectors, or if these subvectors are compatible.

   If any subvector does not implement the N_VInvTestLocal routine (NULL
   function pointer), then this routine will call N_VInvTest instead. */
booleantype N_VInvTestLocal_ManyVector(N_Vector x, N_Vector z)
{
  sunindextype i;
  booleantype val, subval;

  /* initialize output*/
  val = SUNTRUE;

  for (i=0; i<MANYVECTOR_NUM_SUBVECS(x); i++) {

    /* check for nvinvtestlocal in subvector */
    if (MANYVECTOR_SUBVEC(x,i)->ops->nvinvtestlocal) {

      subval = N_VInvTestLocal(MANYVECTOR_SUBVEC(x,i), MANYVECTOR_SUBVEC(z,i));
      val = (val && subval);

    /* otherwise, call nvinvtest and accumulate to overall val */
    } else {

      subval = N_VInvTest(MANYVECTOR_SUBVEC(x,i), MANYVECTOR_SUBVEC(z,i));
      val = (val && subval);

    }
  }

  return(val);
}


/* Performs the InvTest for a ManyVector by calling N_VInvTestLocal and
   combining the results. This routine does not check that x and z
   are ManyVectors, if they have the same number of subvectors, or if these
   subvectors are compatible. */
booleantype N_VInvTest_ManyVector(N_Vector x, N_Vector z)
{
  realtype val, gval;
  val = (N_VInvTestLocal_ManyVector(x, z)) ? ONE : ZERO;
  SUNMPI_Allreduce_scalar(val, &gval, SUNMPI_MIN, MANYVECTOR_COMM(x));
  return (gval != ZERO);
}


/* Performs the MPI task-local ConstrMask for a ManyVector by calling N_VConstrMaskLocal
   on all subvectors and combining the results appropriately.  This routine does not
   check that c, x and m are ManyVectors, if they have the same number of subvectors,
   or if these subvectors are compatible.

   If any subvector does not implement the N_VConstrMaskLocal routine (NULL
   function pointer), then this routine will call N_VConstrMask instead. */
booleantype N_VConstrMaskLocal_ManyVector(N_Vector c, N_Vector x, N_Vector m)
{
  sunindextype i;
  booleantype val, subval;

  /* initialize output*/
  val = SUNTRUE;

  for (i=0; i<MANYVECTOR_NUM_SUBVECS(x); i++) {

    /* check for nvconstrmasklocal in subvector */
    if (MANYVECTOR_SUBVEC(x,i)->ops->nvconstrmasklocal) {

      subval = N_VConstrMaskLocal(MANYVECTOR_SUBVEC(c,i), MANYVECTOR_SUBVEC(x,i),
                                  MANYVECTOR_SUBVEC(m,i));
      val = (val && subval);

    /* otherwise, call nvconstrmask and accumulate to overall val */
    } else {

      subval = N_VConstrMask(MANYVECTOR_SUBVEC(c,i), MANYVECTOR_SUBVEC(x,i),
                             MANYVECTOR_SUBVEC(m,i));
      val = (val && subval);

    }
  }

  return(val);
}


/* Performs the ConstrMask for a ManyVector by calling N_VConstrMaskLocal and
   combining the results.  This routine does not check that c, x and m
   are ManyVectors, if they have the same number of subvectors, or if these
   subvectors are compatible. */
booleantype N_VConstrMask_ManyVector(N_Vector c, N_Vector x, N_Vector m)
{
  realtype val, gval;
  val = (N_VConstrMaskLocal_ManyVector(c, x, m)) ? ONE : ZERO;
  SUNMPI_Allreduce_scalar(val, &gval, SUNMPI_MIN, MANYVECTOR_COMM(x));
  return (gval != ZERO);
}


/* Performs the MPI task-local MinQuotient for a ManyVector by calling N_VMinQuotientLocal
   on all subvectors and combining the results appropriately.  This routine does not check
   that num and denom are ManyVectors, if they have the same number of subvectors, or if
   these subvectors are compatible.

   If any subvector does not implement the N_VMinQuotientLocal routine (NULL
   function pointer), then this routine will call N_VMinQuotient instead. */
realtype N_VMinQuotientLocal_ManyVector(N_Vector num, N_Vector denom)
{
  sunindextype i;
  realtype min, lmin;

  /* initialize output*/
  min = BIG_REAL;

  for (i=0; i<MANYVECTOR_NUM_SUBVECS(num); i++) {

    /* check for nvminquotientlocal in subvector */
    if (MANYVECTOR_SUBVEC(num,i)->ops->nvminquotientlocal) {

      lmin = N_VMinQuotientLocal(MANYVECTOR_SUBVEC(num,i),
                                 MANYVECTOR_SUBVEC(denom,i));
      min = (min < lmin) ? min : lmin;

    /* otherwise, call nvmin and accumulate to overall min */
    } else {

      lmin = N_VMinQuotient(MANYVECTOR_SUBVEC(num,i),
                            MANYVECTOR_SUBVEC(denom,i));
      min = (min < lmin) ? min : lmin;

    }
  }

  return(min);
}


/* Performs the MinQuotient for a ManyVector by calling N_VMinQuotientLocal
   and combining the results.  This routine does not check that num and
   denom are ManyVectors, if they have the same number of subvectors, or if
   these subvectors are compatible. */
realtype N_VMinQuotient_ManyVector(N_Vector num, N_Vector denom)
{
  realtype gmin;
  SUNMPI_Allreduce_scalar(N_VMinQuotientLocal_ManyVector(num, denom),
                          &gmin, SUNMPI_MIN, MANYVECTOR_COMM(num));
  return(gmin);
}




/* -----------------------------------------------------------------
   Fused vector operations
   ----------------------------------------------------------------- */

/* Performs the linear combination z = sum_j c[j]*X[j] by calling
   N_VLinearCombination on all subvectors; this routine does not check that z
   or the components of X are ManyVectors, if they have the same number of
   subvectors, or if these subvectors are compatible.

   NOTE: implementation of this routine is more challenging, due to the
   array-of-arrays of N_Vectors that comprise X.  This routine will be
   passed an array of ManyVectors, so to call the subvector-specific routines
   we must unravel the subvectors while retaining an array of outer vectors. */
int N_VLinearCombination_ManyVector(int nvec, realtype* c, N_Vector* X, N_Vector z)
{
  sunindextype i, j;
  int retval;
  N_Vector *Xsub;

  /* create array of nvec N_Vector pointers for reuse within loop */
  Xsub = NULL;
  Xsub = (N_Vector *) malloc( nvec * sizeof(N_Vector) );
  if (Xsub == NULL)  return(1);

  /* perform operation by calling N_VLinearCombination for each subvector */
  for (i=0; i<MANYVECTOR_NUM_SUBVECS(z); i++) {

    /* for each subvector, create the array of subvectors of X */
    for (j=0; j<nvec; j++)  Xsub[j] = MANYVECTOR_SUBVEC(X[j],i);

    /* now call N_VLinearCombination for this array of subvectors */
    retval = N_VLinearCombination(nvec, c, Xsub, MANYVECTOR_SUBVEC(z,i));

    /* fail gracefully */
    if (retval) {
      free(Xsub);
      return(retval);
    }

  }

  /* clean up and return */
  free(Xsub);
  return(0);
}


/* Performs the ScaleAddMulti operation by calling N_VScaleAddMulti on all
   subvectors; this routine does not check that x, or the components of X and Z are
   ManyVectors, if they have the same number of subvectors, or if these subvectors
   are compatible.

   NOTE: this routine is more challenging, due to the array-of-arrays of
   N_Vectors that comprise Y and Z.  This routine will be passed an array of
   ManyVectors, so to call the subvector-specific routines we must unravel
   the subvectors while retaining an array of outer vectors. */
int N_VScaleAddMulti_ManyVector(int nvec, realtype* a, N_Vector x, N_Vector* Y, N_Vector* Z)
{
  sunindextype i, j;
  int retval;
  N_Vector *Ysub, *Zsub;

  /* create arrays of nvec N_Vector pointers for reuse within loop */
  Ysub = Zsub = NULL;
  Ysub = (N_Vector *) malloc( nvec * sizeof(N_Vector) );
  Zsub = (N_Vector *) malloc( nvec * sizeof(N_Vector) );
  if ( (Ysub == NULL) || (Zsub == NULL) )  return(1);

  /* perform operation by calling N_VScaleAddMulti for each subvector */
  for (i=0; i<MANYVECTOR_NUM_SUBVECS(x); i++) {

    /* for each subvector, create the array of subvectors of Y and Z */
    for (j=0; j<nvec; j++)  {
      Ysub[j] = MANYVECTOR_SUBVEC(Y[j],i);
      Zsub[j] = MANYVECTOR_SUBVEC(Z[j],i);
    }

    /* now call N_VScaleAddMulti for this array of subvectors */
    retval = N_VScaleAddMulti(nvec, a, MANYVECTOR_SUBVEC(x,i), Ysub, Zsub);

    /* fail gracefully */
    if (retval) {
      free(Ysub);
      free(Zsub);
      return(retval);
    }

  }

  /* clean up and return */
  free(Ysub);
  free(Zsub);
  return(0);
}


/* Performs the DotProdMulti operation by calling N_VDotProdLocal and combining results.
   This routine does not check that x, or the components of Y, are ManyVectors, if
   they have the same number of subvectors, or if these subvectors are compatible.

   NOTE: if all subvectors implement the N_VDotProdLocal routine, then this routine
   will require only a single array-valued reduction operation (in contrast to calling
   N_VDotProdMulti on all subvectors, where we would require num_subvectors separate
   reductions). */
int N_VDotProdMulti_ManyVector(int nvec, N_Vector x, N_Vector* Y, realtype* dotprods)
{
  sunindextype i;

  /* call N_VDotProdLocal for each <x,Y[i]> pair */
  for (i=0; i<nvec; i++)  dotprods[i] = N_VDotProdLocal_ManyVector(x,Y[i]);

  /* accumulate totals and return */
  return(SUNMPI_Allreduce(dotprods, nvec, SUNMPI_SUM, MANYVECTOR_COMM(x)));
}


/* -----------------------------------------------------------------
   Vector array operations
   ----------------------------------------------------------------- */

/* Performs the LinearSumVectorArray operation by calling N_VLinearSumVectorArray
   on all subvectors; this routine does not check that the components of X, Y or Z are
   ManyVectors, if they have the same number of subvectors, or if these subvectors
   are compatible.

   NOTE: this routine is more challenging, due to the array-of-arrays of
   N_Vectors that comprise X, Y and Z.  This routine will be passed arrays of
   ManyVectors, so to call the subvector-specific routines we must unravel
   the subvectors while retaining arrays of outer vectors. */
int N_VLinearSumVectorArray_ManyVector(int nvec, realtype a,
                                       N_Vector *X, realtype b,
                                       N_Vector *Y, N_Vector *Z)
{
  sunindextype i, j;
  int retval;
  N_Vector *Xsub, *Ysub, *Zsub;

  /* immediately return if nvec <= 0 */
  if (nvec <= 0)  return(0);

  /* create arrays of nvec N_Vector pointers for reuse within loop */
  Xsub = Ysub = Zsub = NULL;
  Xsub = (N_Vector *) malloc( nvec * sizeof(N_Vector) );
  Ysub = (N_Vector *) malloc( nvec * sizeof(N_Vector) );
  Zsub = (N_Vector *) malloc( nvec * sizeof(N_Vector) );
  if ( (Xsub == NULL) || (Ysub == NULL) || (Zsub == NULL) )  return(1);

  /* perform operation by calling N_VLinearSumVectorArray for each subvector */
  for (i=0; i<MANYVECTOR_NUM_SUBVECS(X[0]); i++) {

    /* for each subvector, create the array of subvectors of X, Y and Z */
    for (j=0; j<nvec; j++)  {
      Xsub[j] = MANYVECTOR_SUBVEC(X[j],i);
      Ysub[j] = MANYVECTOR_SUBVEC(Y[j],i);
      Zsub[j] = MANYVECTOR_SUBVEC(Z[j],i);
    }

    /* now call N_VLinearSumVectorArray for this array of subvectors */
    retval = N_VLinearSumVectorArray(nvec, a, Xsub, b, Ysub, Zsub);

    /* fail gracefully */
    if (retval) {
      free(Xsub);
      free(Ysub);
      free(Zsub);
      return(retval);
    }

  }

  /* clean up and return */
  free(Xsub);
  free(Ysub);
  free(Zsub);
  return(0);
}


/* Performs the ScaleVectorArray operation by calling N_VScaleVectorArray
   on all subvectors; this routine does not check that the components of X or Z are
   ManyVectors, if they have the same number of subvectors, or if these subvectors
   are compatible.

   NOTE: this routine is more challenging, due to the array-of-arrays of
   N_Vectors that comprise X and Z.  This routine will be passed arrays of
   ManyVectors, so to call the subvector-specific routines we must unravel
   the subvectors while retaining arrays of outer vectors. */
int N_VScaleVectorArray_ManyVector(int nvec, realtype* c, N_Vector* X, N_Vector* Z)
{
  sunindextype i, j;
  int retval;
  N_Vector *Xsub, *Zsub;

  /* immediately return if nvec <= 0 */
  if (nvec <= 0)  return(0);

  /* create arrays of nvec N_Vector pointers for reuse within loop */
  Xsub = Zsub = NULL;
  Xsub = (N_Vector *) malloc( nvec * sizeof(N_Vector) );
  Zsub = (N_Vector *) malloc( nvec * sizeof(N_Vector) );
  if ( (Xsub == NULL) || (Zsub == NULL) )  return(1);

  /* perform operation by calling N_VScaleVectorArray for each subvector */
  for (i=0; i<MANYVECTOR_NUM_SUBVECS(X[0]); i++) {

    /* for each subvector, create the array of subvectors of X, Y and Z */
    for (j=0; j<nvec; j++)  {
      Xsub[j] = MANYVECTOR_SUBVEC(X[j],i);
      Zsub[j] = MANYVECTOR_SUBVEC(Z[j],i);
    }

    /* now call N_VScaleVectorArray for this array of subvectors */
    retval = N_VScaleVectorArray(nvec, c, Xsub, Zsub);

    /* fail gracefully */
    if (retval) {
      free(Xsub);
      free(Zsub);
      return(retval);
    }

  }

  /* clean up and return */
  free(Xsub);
  free(Zsub);
  return(0);
}


/* Performs the ConstVectorArray operation by calling N_VConstVectorArray
   on all subvectors.

   NOTE: this routine is more challenging, due to the array-of-arrays of
   N_Vectors that comprise Z.  This routine will be passed an array of
   ManyVectors, so to call the subvector-specific routines we must unravel
   the subvectors while retaining an array of outer vectors. */
int N_VConstVectorArray_ManyVector(int nvec, realtype c, N_Vector* Z)
{
  sunindextype i, j;
  int retval;
  N_Vector *Zsub;

  /* immediately return if nvec <= 0 */
  if (nvec <= 0)  return(0);

  /* create array of N_Vector pointers for reuse within loop */
  Zsub = NULL;
  Zsub = (N_Vector *) malloc( nvec * sizeof(N_Vector) );
  if (Zsub == NULL)  return(1);

  /* perform operation by calling N_VConstVectorArray for each subvector */
  for (i=0; i<MANYVECTOR_NUM_SUBVECS(Z[0]); i++) {

    /* for each subvector, create the array of subvectors of X, Y and Z */
    for (j=0; j<nvec; j++)
      Zsub[j] = MANYVECTOR_SUBVEC(Z[j],i);

    /* now call N_VConstVectorArray for this array of subvectors */
    retval = N_VConstVectorArray(nvec, c, Zsub);

    /* fail gracefully */
    if (retval) {
      free(Zsub);
      return(retval);
    }

  }

  /* clean up and return */
  free(Zsub);
  return(0);
}


/* Performs the WrmsNormVectorArray operation by calling N_VWSqrSumLocal and combining
   results.  This routine does not check that the components of X or W are ManyVectors, if
   they have the same number of subvectors, or if these subvectors are compatible.

   NOTE: if all subvectors implement the N_VWSqrSumLocal routine, then this routine
   will require only a single array-valued reduction operation (in contrast to calling
   N_VWrmsNormVectorArray on all subvectors, where we would require num_subvectors
   separate reductions). */
int N_VWrmsNormVectorArray_ManyVector(int nvec, N_Vector* X, N_Vector* W, realtype* nrm)
{
  sunindextype i;
  int retval;

  /* immediately return if nvec <= 0 */
  if (nvec <= 0)  return(0);

  /* call N_VWSqrSumLocal for each (X[i],W[i]) pair */
  for (i=0; i<nvec; i++)  nrm[i] = N_VWSqrSumLocal_ManyVector(X[i], W[i]);

  /* accumulate totals */
  retval = SUNMPI_Allreduce(nrm, nvec, SUNMPI_SUM, MANYVECTOR_COMM(X[0]));

  /* finish off WRMS norms and return */
  for (i=0; i<nvec; i++)
    nrm[i] = SUNRsqrt(nrm[i]/(MANYVECTOR_GLOBLENGTH(X[i])));

  return (retval == SUNMPI_SUCCESS) ? 0 : -1;
}


/* Performs the WrmsNormMaskVectorArray operation by calling N_VWSqrSumMaskLocal and
   combining results.  This routine does not check that id or the components of X and
   W are ManyVectors, if they have the same number of subvectors, or if these
   subvectors are compatible.

   NOTE: if all subvectors implement the N_VWSqrSumMaskLocal routine, then this
   routine will require only a single array-valued reduction operation (in contrast
   to calling N_VWrmsNormMaskVectorArray on all subvectors, where we would require
   num_subvectors separate reductions). */
int N_VWrmsNormMaskVectorArray_ManyVector(int nvec, N_Vector* X, N_Vector* W,
                                          N_Vector id, realtype* nrm)
{
  sunindextype i;
  int retval;

  /* immediately return if nvec <= 0 */
  if (nvec <= 0)  return(0);

  /* call N_VWSqrSumMaskLocal for each (X[i],W[i]) pair */
  for (i=0; i<nvec; i++)  nrm[i] = N_VWSqrSumMaskLocal_ManyVector(X[i], W[i], id);

  /* accumulate totals */
  retval = SUNMPI_Allreduce(nrm, nvec, SUNMPI_SUM, MANYVECTOR_COMM(X[0]));

  /* finish off WRMS norms and return */
  for (i=0; i<nvec; i++)
    nrm[i] = SUNRsqrt(nrm[i]/(MANYVECTOR_GLOBLENGTH(X[i])));

  return (retval == SUNMPI_SUCCESS) ? 0 : -1;
}


/* -----------------------------------------------------------------
   Enable / Disable fused and vector array operations
   ----------------------------------------------------------------- */

int N_VEnableFusedOps_ManyVector(N_Vector v, booleantype tf)
{
  /* check that vector is non-NULL */
  if (v == NULL) return(-1);

  /* check that ops structure is non-NULL */
  if (v->ops == NULL) return(-1);

  if (tf) {
    /* enable all fused vector operations */
    v->ops->nvlinearcombination = N_VLinearCombination_ManyVector;
    v->ops->nvscaleaddmulti     = N_VScaleAddMulti_ManyVector;
    v->ops->nvdotprodmulti      = N_VDotProdMulti_ManyVector;
    /* enable all vector array operations */
    v->ops->nvlinearsumvectorarray         = N_VLinearSumVectorArray_ManyVector;
    v->ops->nvscalevectorarray             = N_VScaleVectorArray_ManyVector;
    v->ops->nvconstvectorarray             = N_VConstVectorArray_ManyVector;
    v->ops->nvwrmsnormvectorarray          = N_VWrmsNormVectorArray_ManyVector;
    v->ops->nvwrmsnormmaskvectorarray      = N_VWrmsNormMaskVectorArray_ManyVector;
    v->ops->nvscaleaddmultivectorarray     = NULL;
    v->ops->nvlinearcombinationvectorarray = NULL;
  } else {
    /* disable all fused vector operations */
    v->ops->nvlinearcombination = NULL;
    v->ops->nvscaleaddmulti     = NULL;
    v->ops->nvdotprodmulti      = NULL;
    /* disable all vector array operations */
    v->ops->nvlinearsumvectorarray         = NULL;
    v->ops->nvscalevectorarray             = NULL;
    v->ops->nvconstvectorarray             = NULL;
    v->ops->nvwrmsnormvectorarray          = NULL;
    v->ops->nvwrmsnormmaskvectorarray      = NULL;
    v->ops->nvscaleaddmultivectorarray     = NULL;
    v->ops->nvlinearcombinationvectorarray = NULL;
  }

  /* return success */
  return(0);
}


int N_VEnableLinearCombination_ManyVector(N_Vector v, booleantype tf)
{
  /* check that vector is non-NULL */
  if (v == NULL) return(-1);

  /* check that ops structure is non-NULL */
  if (v->ops == NULL) return(-1);

  /* enable/disable operation */
  if (tf)
    v->ops->nvlinearcombination = N_VLinearCombination_ManyVector;
  else
    v->ops->nvlinearcombination = NULL;

  /* return success */
  return(0);
}

int N_VEnableScaleAddMulti_ManyVector(N_Vector v, booleantype tf)
{
  /* check that vector is non-NULL */
  if (v == NULL) return(-1);

  /* check that ops structure is non-NULL */
  if (v->ops == NULL) return(-1);

  /* enable/disable operation */
  if (tf)
    v->ops->nvscaleaddmulti = N_VScaleAddMulti_ManyVector;
  else
    v->ops->nvscaleaddmulti = NULL;

  /* return success */
  return(0);
}

int N_VEnableDotProdMulti_ManyVector(N_Vector v, booleantype tf)
{
  /* check that vector is non-NULL */
  if (v == NULL) return(-1);

  /* check that ops structure is non-NULL */
  if (v->ops == NULL) return(-1);

  /* enable/disable operation */
  if (tf)
    v->ops->nvdotprodmulti = N_VDotProdMulti_ManyVector;
  else
    v->ops->nvdotprodmulti = NULL;

  /* return success */
  return(0);
}

int N_VEnableLinearSumVectorArray_ManyVector(N_Vector v, booleantype tf)
{
  /* check that vector is non-NULL */
  if (v == NULL) return(-1);

  /* check that ops structure is non-NULL */
  if (v->ops == NULL) return(-1);

  /* enable/disable operation */
  if (tf)
    v->ops->nvlinearsumvectorarray = N_VLinearSumVectorArray_ManyVector;
  else
    v->ops->nvlinearsumvectorarray = NULL;

  /* return success */
  return(0);
}

int N_VEnableScaleVectorArray_ManyVector(N_Vector v, booleantype tf)
{
  /* check that vector is non-NULL */
  if (v == NULL) return(-1);

  /* check that ops structure is non-NULL */
  if (v->ops == NULL) return(-1);

  /* enable/disable operation */
  if (tf)
    v->ops->nvscalevectorarray = N_VScaleVectorArray_ManyVector;
  else
    v->ops->nvscalevectorarray = NULL;

  /* return success */
  return(0);
}

int N_VEnableConstVectorArray_ManyVector(N_Vector v, booleantype tf)
{
  /* check that vector is non-NULL */
  if (v == NULL) return(-1);

  /* check that ops structure is non-NULL */
  if (v->ops == NULL) return(-1);

  /* enable/disable operation */
  if (tf)
    v->ops->nvconstvectorarray = N_VConstVectorArray_ManyVector;
  else
    v->ops->nvconstvectorarray = NULL;

  /* return success */
  return(0);
}

int N_VEnableWrmsNormVectorArray_ManyVector(N_Vector v, booleantype tf)
{
  /* check that vector is non-NULL */
  if (v == NULL) return(-1);

  /* check that ops structure is non-NULL */
  if (v->ops == NULL) return(-1);

  /* enable/disable operation */
  if (tf)
    v->ops->nvwrmsnormvectorarray = N_VWrmsNormVectorArray_ManyVector;
  else
    v->ops->nvwrmsnormvectorarray = NULL;

  /* return success */
  return(0);
}

int N_VEnableWrmsNormMaskVectorArray_ManyVector(N_Vector v, booleantype tf)
{
  /* check that vector is non-NULL */
  if (v == NULL) return(-1);

  /* check that ops structure is non-NULL */
  if (v->ops == NULL) return(-1);

  /* enable/disable operation */
  if (tf)
    v->ops->nvwrmsnormmaskvectorarray = N_VWrmsNormMaskVectorArray_ManyVector;
  else
    v->ops->nvwrmsnormmaskvectorarray = NULL;

  /* return success */
  return(0);
}


/* -----------------------------------------------------------------
   Implementation of utility routines
   -----------------------------------------------------------------*/

/* This function performs a generic clone operation on an input N_Vector.
   Based on the 'cloneempty' flag it will either call "nvclone" or
   "nvcloneempty" when creating subvectors in the cloned vector. */
static N_Vector ManyVectorClone(N_Vector w, booleantype cloneempty)
{
  N_Vector v;
  N_Vector_Ops ops;
  N_VectorContent_ManyVector content;
  booleantype fail_gracefully;
  sunindextype i;

  if (w == NULL) return(NULL);

  /* Create vector */
  v = NULL;
  v = (N_Vector) malloc(sizeof *v);
  if (v == NULL) return(NULL);

  /* Create vector operation structure */
  ops = NULL;
  ops = (N_Vector_Ops) malloc(sizeof(struct _generic_N_Vector_Ops));
  if (ops == NULL) { free(v); return(NULL); }

  /* standard vector operations */
  ops->nvgetvectorid     = w->ops->nvgetvectorid;
  ops->nvcloneempty      = w->ops->nvcloneempty;
  ops->nvclone           = w->ops->nvclone;
  ops->nvdestroy         = w->ops->nvdestroy;
  ops->nvspace           = w->ops->nvspace;
  ops->nvgetarraypointer = w->ops->nvgetarraypointer;
  ops->nvsetarraypointer = w->ops->nvsetarraypointer;
  ops->nvgetcommunicator = w->ops->nvgetcommunicator;
  ops->nvgetlength       = w->ops->nvgetlength;
  ops->nvlinearsum       = w->ops->nvlinearsum;
  ops->nvconst           = w->ops->nvconst;
  ops->nvprod            = w->ops->nvprod;
  ops->nvdiv             = w->ops->nvdiv;
  ops->nvscale           = w->ops->nvscale;
  ops->nvabs             = w->ops->nvabs;
  ops->nvinv             = w->ops->nvinv;
  ops->nvaddconst        = w->ops->nvaddconst;
  ops->nvdotprod         = w->ops->nvdotprod;
  ops->nvmaxnorm         = w->ops->nvmaxnorm;
  ops->nvwrmsnormmask    = w->ops->nvwrmsnormmask;
  ops->nvwrmsnorm        = w->ops->nvwrmsnorm;
  ops->nvmin             = w->ops->nvmin;
  ops->nvwl2norm         = w->ops->nvwl2norm;
  ops->nvl1norm          = w->ops->nvl1norm;
  ops->nvcompare         = w->ops->nvcompare;
  ops->nvinvtest         = w->ops->nvinvtest;
  ops->nvconstrmask      = w->ops->nvconstrmask;
  ops->nvminquotient     = w->ops->nvminquotient;

  /* fused vector operations */
  ops->nvlinearcombination = w->ops->nvlinearcombination;
  ops->nvscaleaddmulti     = w->ops->nvscaleaddmulti;
  ops->nvdotprodmulti      = w->ops->nvdotprodmulti;

  /* vector array operations */
  ops->nvlinearsumvectorarray         = w->ops->nvlinearsumvectorarray;
  ops->nvscalevectorarray             = w->ops->nvscalevectorarray;
  ops->nvconstvectorarray             = w->ops->nvconstvectorarray;
  ops->nvwrmsnormvectorarray          = w->ops->nvwrmsnormvectorarray;
  ops->nvwrmsnormmaskvectorarray      = w->ops->nvwrmsnormmaskvectorarray;
  ops->nvscaleaddmultivectorarray     = w->ops->nvscaleaddmultivectorarray;
  ops->nvlinearcombinationvectorarray = w->ops->nvlinearcombinationvectorarray;

  /* local reduction kernels */
  ops->nvdotprodlocal     = w->ops->nvdotprodlocal;
  ops->nvmaxnormlocal     = w->ops->nvmaxnormlocal;
  ops->nvminlocal         = w->ops->nvminlocal;
  ops->nvl1normlocal      = w->ops->nvl1normlocal;
  ops->nvinvtestlocal     = w->ops->nvinvtestlocal;
  ops->nvconstrmasklocal  = w->ops->nvconstrmasklocal;
  ops->nvminquotientlocal = w->ops->nvminquotientlocal;
  ops->nvwsqrsumlocal     = w->ops->nvwsqrsumlocal;
  ops->nvwsqrsummasklocal = w->ops->nvwsqrsummasklocal;

  /* Create content */
  content = NULL;
  content = (N_VectorContent_ManyVector) malloc(sizeof(struct _N_VectorContent_ManyVector));
  if (content == NULL) { free(ops); free(v); return(NULL); }


  /* Attach content components */

  /*   Set scalar components */
  content->comm = MANYVECTOR_COMM(w);
  content->num_subvectors = MANYVECTOR_NUM_SUBVECS(w);
  content->global_length = MANYVECTOR_GLOBLENGTH(w);
  content->own_data = SUNTRUE;

  /*   Set vector components */
  content->subvec_array = NULL;
  content->subvec_array = (N_Vector *) malloc( content->num_subvectors * sizeof(N_Vector) );
  if (content->subvec_array == NULL)
    { free(ops); free(content); free(v); return(NULL); }
  fail_gracefully = SUNFALSE;
  for (i=0; i<content->num_subvectors; i++) {
    if (cloneempty) {
      content->subvec_array[i] = N_VCloneEmpty(MANYVECTOR_SUBVEC(w,i));
    } else {
      content->subvec_array[i] = N_VClone(MANYVECTOR_SUBVEC(w,i));
    }
    if (content->subvec_array[i] == NULL) {
      fail_gracefully = SUNTRUE;
      break;
    }
  }

  /*   fail gracefully if any Clone failed */
  if (fail_gracefully) {
    for (i=0; i<content->num_subvectors; i++) {
      if (content->subvec_array[i] != NULL)
        N_VDestroy(content->subvec_array[i]);
    }
    free(content->subvec_array);
    free(ops);
    free(content);
    free(v);
    return(NULL);
  }

  /* Attach content and ops to new vector, and return */
  v->content = content;
  v->ops     = ops;

  return(v);
}


/* This function returns the rank of this task in the MPI communicator
   associated with the input N_Vector.  If the input N_Vector is MPI-unaware, it
   returns 0.  If an error occurs in the call to SUNMPI_Comm_Rank, it returns -1. */
static int SubvectorMPIRank(N_Vector x)
{
  void* tmpcomm;
  SUNMPI_Comm *comm;
  int rank, retval;
  tmpcomm = N_VGetCommunicator(x);
  if (tmpcomm == NULL) return(0);
  comm = (SUNMPI_Comm *) tmpcomm;
  retval = SUNMPI_Comm_rank(*comm, &rank);
  if (retval != SUNMPI_SUCCESS)  return(-1);
  return(rank);
}
