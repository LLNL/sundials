#ifndef _CHOMBOSUNDIALS_H
#define _CHOMBOSUNDIALS_H

#include <sundials/sundials_types.h>
#include <sundials/sundials_nvector.h>

#include "LevelData.H"
#include "FArrayBox.H"
#include "DisjointBoxLayout.H"
#include "CHArray.H"
#include "REAL.H"

#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

// define MPI data types
// real number
#if defined(SUNDIALS_SINGLE_PRECISION)
#define CH_REAL_MPI_TYPE MPI_FLOAT
#elif defined(SUNDIALS_DOUBLE_PRECISION)
#define CH_REAL_MPI_TYPE MPI_DOUBLE
#elif defined(SUNDIALS_EXTENDED_PRECISION)
#define CH_REAL_MPI_TYPE MPI_LONG_DOUBLE
#endif
// integer
#define CH_INTEGER_MPI_TYPE MPI_LONG

// Forward declare pointer-to type
typedef struct _N_VectorContent_Ch *N_VectorContent_Ch;

// content
struct _N_VectorContent_Ch {
  MPI_Comm comm;
  LevelData<FArrayBox>* Data;
  bool own_data;
};

// public API
N_Vector N_VAssociate(LevelData<FArrayBox>* Data);
N_Vector N_VNew(const DisjointBoxLayout& dp, const int nComp, const IntVect& ghost);
void N_VDataExchange(N_Vector v);
int N_VEquate(N_Vector x, N_Vector z);

#define CHSUN_SUCCESS    0
#define CHSUN_SIZE_ERROR 1

// accessor macros
#define NV_OPS_CH(v)     ( v->ops )
#define NV_CONTENT_CH(v) ( static_cast<N_VectorContent_Ch>(v->content) )
#define NV_COMM_CH(v)    ( NV_CONTENT_CH(v)->comm )
#define NV_OWNDATA_CH(v)   ( NV_CONTENT_CH(v)->own_data )
#define NV_DATA_CH(v)    ( (*NV_CONTENT_CH(v)->Data) )
#define NV_DBL_CH(v)     ( NV_DATA_CH(v).disjointBoxLayout() )
#define NV_DBLITER_CH(v) ( NV_DBL_CH(v).dataIterator() )
#define NV_GHOST_CH(v)   ( NV_DATA_CH(v).ghostVect() )
#define NV_NCOMP_CH(v)   ( NV_DATA_CH(v).nComp() )

#ifdef __cplusplus
} // extern "C"
#endif

#endif
