/*
 * -----------------------------------------------------------------
 * $Revision: 1.1 $
 * $Date: 2005-06-14 19:00:57 $
 * ----------------------------------------------------------------- 
 * Programmer(s): Daniel R. Reynolds and Radu Serban @LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2002, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see sundials/shared/LICENSE.
 *-----------------------------------------------------------------
 * This is the header file for an implementation of an NVECTOR
 * package for which the local data consists of 3 space dimensions 
 * + nspc species i.e. u(Xlo:Xhi, Ylo:Yhi, Zlo:Zhi, 1:nspc),
 * and the indices Xlo:Xhi, etc. may include ghost data for the
 * local domain.
 *
 * Part I contains declarations specific to the spcparallel
 * implementation of the supplied NVECTOR module.
 *
 * Part II defines accessor macros that allow the user to efficiently
 * use the type N_Vector without making explicit references to the
 * underlying data structure.
 *
 * Part III contains the prototype for the constructor
 * N_VNew_SpcParallel as well as implementation-specific prototypes
 * for various useful vector operations.
 *
 * Notes:
 *
 *   - The definition of the generic N_Vector structure can be
 *     found in the header file shared/include/nvector.h.
 *
 *   - The definition of the type realtype can be found in the
 *     header file shared/include/sundialstypes.h, and it may be
 *     changed (at the configuration stage) according to the user's
 *     needs. The sundialstypes.h file also contains the definition
 *     for the type booleantype.
 *
 *   - N_Vector arguments to arithmetic vector operations need not
 *     be distinct. For example, the following call:
 *
 *        N_VLinearSum_SpcParallel(a,x,b,y,y);
 *
 *     (which stores the result of the operation a*x+b*y in y)
 *     is legal.
 * -----------------------------------------------------------------
 */

#ifndef _NVECTOR_SPCPARALLEL_H
#define _NVECTOR_SPCPARALLEL_H

#ifdef __cplusplus     /* wrapper to enable C++ usage */
extern "C" {
#endif
  
#include "mpi.h"
#include "nvector.h"
#include "sundialstypes.h"

/*
 * -----------------------------------------------------------------
 * PART I: SPCPARALLEL implementation of N_Vector               
 * -----------------------------------------------------------------
 */
  
  /* define MPI data types */

#if defined(SUNDIALS_SINGLE_PRECISION)
#define SPVEC_REAL_MPI_TYPE MPI_FLOAT

#elif defined(SUNDIALS_DOUBLE_PRECISION)
#define SPVEC_REAL_MPI_TYPE MPI_DOUBLE

#elif defined(SUNDIALS_EXTENDED_PRECISION)
#define SPVEC_REAL_MPI_TYPE MPI_LONG_DOUBLE

#endif

#define SPVEC_INTEGER_MPI_TYPE MPI_LONG

/* The SPCPARALLEL implementation of the N_Vector 'content' structure 
 * contains the global length of the vector, local valid lengths and 
 * ghost dimensions in x, y and z-directions of the vector, a pointer 
 * to an array of real components, and the MPI communicator 
 */

struct _N_VectorContent_SpcParallel {
  int nspc;                /* number of species                               */
  long int Nx;           /* local x-mesh vector length                      */
  long int Ny;           /* local y-mesh vector length                      */
  long int Nz;           /* local z-mesh vector length                      */
  long int NGx;          /* x-width of ghost boundary                       */
  long int NGy;          /* y-width of ghost boundary                       */
  long int NGz;          /* z-width of ghost boundary                       */
  long int n1;           /* local computational vector length for 1 species */
  long int n1g;          /* local vector vector for 1 species               */
  long int n;            /* local computational vector length               */
  long int ng;           /* local vector length                             */
  long int Ng;           /* global vector length                            */
  realtype *data;        /* local data array                                */
  booleantype own_data;  /* flag for ownership of data                      */
  MPI_Comm comm;         /* pointer to MPI communicator                     */
};
  
typedef struct _N_VectorContent_SpcParallel *N_VectorContent_SpcParallel;
  
/*
 * -----------------------------------------------------------------
 * PART II: Macros                                              
 * --------------------------------------------------------------
 * Notes
 *
 * When looping over the components of an N_Vector v, it is
 * more efficient to first obtain the component array via
 * v_data = NV_DATA_MHD(v) and then access v_data[i] within
 * the loop than it is to use NV_Ith_MHD(v,i) within the
 * loop.
 *
 * -----------------------------------------------------------------
 */

  /* Vector content */

#define SPV_CONTENT(v) ( (N_VectorContent_SpcParallel)(v->content) )

  /* Number of species */

#define SPV_NSPECIES(v) (SPV_CONTENT(v)->nspc)

  /* Local grid dimensions */

#define SPV_XLENGTH(v) (SPV_CONTENT(v)->Nx)
#define SPV_YLENGTH(v) (SPV_CONTENT(v)->Ny)
#define SPV_ZLENGTH(v) (SPV_CONTENT(v)->Nz)

#define SPV_XGHOST(v) (SPV_CONTENT(v)->NGx)
#define SPV_YGHOST(v) (SPV_CONTENT(v)->NGy)
#define SPV_ZGHOST(v) (SPV_CONTENT(v)->NGz)

#define SPV_XGLENGTH(v) (SPV_XLENGTH(v)+2*SPV_XGHOST(v))
#define SPV_YGLENGTH(v) (SPV_YLENGTH(v)+2*SPV_YGHOST(v))
#define SPV_ZGLENGTH(v) (SPV_ZLENGTH(v)+2*SPV_ZGHOST(v))

  /* Local vector lengths */

#define SPV_LOCLENGTH1(v) (SPV_CONTENT(v)->n1g)
#define SPV_LOCLENGTH(v)  (SPV_CONTENT(v)->ng)

  /* Global length */

#define SPV_GLOBLENGTH(v) (SPV_CONTENT(v)->Ng)

  /* Data access */

#define SPV_DATA(v) (SPV_CONTENT(v)->data)

#define SPV_SDATA(v,s) (SPV_DATA(v) + ((s)*SPV_LOCLENGTH1(v)))

#define SPV_Ith(v,s,i) (SPV_SDATA(v,s)[i])

#define SPV_IJKth(v,s,i,j,k) (SPV_Ith(v,s,(i)+SPV_XGLENGTH(v)*((j)+SPV_YGLENGTH(v)*(k))))

  /* MPI communicator */

#define SPV_COMM(v) (SPV_CONTENT(v)->comm)

  /* Data ownership flag */

#define SPV_OWN_DATA(v) (SPV_CONTENT(v)->own_data)

/*
 * -----------------------------------------------------------------
 * PART III: functions exported by SPCPARALLEL
 * 
 * CONSTRUCTORS:
 *    N_VNew_SpcParallel
 *    N_VNewEmpty_SpcParallel
 *    N_VMake_SpcParallel
 *    N_VNewVectorArray_SpcParallel        NYI
 *    N_VNewVectorArrayEmpty_SpcParallel   NYI
 * DESTRUCTORS:
 *    N_VDestroy_SpcParallel              
 *    N_VDestroyVectorArray_SpcParallel    NYI
 * OTHER:
 *    N_VPrint_SpcParallel
 *    N_VPrintFile_SpcParallel
 * -----------------------------------------------------------------
 */
  
/*
 * -----------------------------------------------------------------
 * Function : N_VNew_SpcParallel
 * -----------------------------------------------------------------
 * This function creates and allocates memory for an 
 * SPCPARALLEL vector.
 * -----------------------------------------------------------------
 */

N_Vector N_VNew_SpcParallel(MPI_Comm comm, int nspc, 
                            long int Nx,  long int Ny,  long int Nz, 
                            long int NGx, long int NGy, long int NGz);
/*
 * -----------------------------------------------------------------
 * Function : N_VNewEmpty_SpcParallel
 * -----------------------------------------------------------------
 * This function creates a new SPCPARALLEL N_Vector with an empty
 * (NULL) data array.
 * -----------------------------------------------------------------
 */

N_Vector N_VNewEmpty_SpcParallel(MPI_Comm comm,  int nspc,
                                 long int Nx,  long int Ny,  long int Nz, 
                                 long int NGx, long int NGy, long int NGz);

/*
 * -----------------------------------------------------------------
 * Function : N_VAttach_SpcParallel
 * -----------------------------------------------------------------
 * This function creates an SPCPARALLEL vector and attaches to it
 * the user-supplied data (which must be a contiguous array)
 * -----------------------------------------------------------------
 */

N_Vector N_VAttach_SpcParallel(MPI_Comm comm, int nspc, 
                               long int Nx,  long int Ny,  long int Nz, 
                               long int NGx, long int NGy, long int NGz,
                               realtype *data);

/*
 * -----------------------------------------------------------------
 * Function : N_VLoad_SpcParallel
 * -----------------------------------------------------------------
 * This function creates an SPCPARALLEL vector and copiesinto it
 * the user-supplied data (which must be an array of nspc vectors)
 * -----------------------------------------------------------------
 */

N_Vector N_VLoad_SpcParallel(MPI_Comm comm, int nspc, 
                             long int Nx,  long int Ny,  long int Nz, 
                             long int NGx, long int NGy, long int NGz,
                             realtype **data);

/*
 * -----------------------------------------------------------------
 * Function : N_VNewVectorArray_SpcParallel
 * -----------------------------------------------------------------
 * This function creates an array of 'count' parallel vectors. This
 * array of N_Vectors can be freed using N_VDestroyVectorArray
 * (defined by the generic NVECTOR module).
 * -----------------------------------------------------------------
 */

N_Vector *N_VNewVectorArray_SpcParallel(int count, 
                                        MPI_Comm comm, int nspc, 
                                        long int Nx,  long int Ny,  long int Nz, 
                                        long int NGx, long int NGy, long int NGz);

/*
 * -----------------------------------------------------------------
 * Function : N_VNewVectorArrayEmpty_SpcParallel
 * -----------------------------------------------------------------
 * This function creates an array of 'count' parallel vectors each 
 * with an empty (NULL) data array.
 * -----------------------------------------------------------------
 */

N_Vector *N_VNewVectorArrayEmpty_SpcParallel(int count, 
                                             MPI_Comm comm, int nspc, 
                                             long int Nx,  long int Ny,  long int Nz, 
                                             long int NGx, long int NGy, long int NGz);

/*
 * -----------------------------------------------------------------
 * Function : N_VDestroyVectorArray_SpcParallel
 * -----------------------------------------------------------------
 * This function frees an array of N_Vector created with 
 * N_VNewVectorArray_SpcParallel.
 * -----------------------------------------------------------------
 */

void N_VDestroyVectorArray_SpcParallel(N_Vector *vs, int count);

/*
 * -----------------------------------------------------------------
 * Function : N_VPrint_SpcParallel
 * -----------------------------------------------------------------
 * This function prints the content of an SPCPARALLEL vector 
 * to stdout.
 * -----------------------------------------------------------------
 */

void N_VPrint_SpcParallel(N_Vector v);

/*
 * -----------------------------------------------------------------
 * Function : N_VPrintFile_SpcParallel
 * -----------------------------------------------------------------
 * This function prints the content of an SPCPARALLEL vector 
 * to a file.
 * -----------------------------------------------------------------
 */

void N_VPrintFile_SpcParallel(char *fname, N_Vector v);

/*
 * -----------------------------------------------------------------
 * SPCPARALLEL implementations of the vector operations
 * -----------------------------------------------------------------
 */

N_Vector    N_VCloneEmpty_SpcParallel(N_Vector);
N_Vector    N_VClone_SpcParallel(N_Vector);
void        N_VDestroy_SpcParallel(N_Vector);
void        N_VSpace_SpcParallel(N_Vector, long int *, long int *);
realtype   *N_VGetArrayPointer_SpcParallel(N_Vector);
void        N_VSetArrayPointer_SpcParallel(realtype *, N_Vector);
void        N_VLinearSum_SpcParallel(realtype, N_Vector, realtype, N_Vector, N_Vector);
void        N_VConst_SpcParallel(realtype, N_Vector);
void        N_VProd_SpcParallel(N_Vector, N_Vector, N_Vector);
void        N_VDiv_SpcParallel(N_Vector, N_Vector, N_Vector);
void        N_VScale_SpcParallel(realtype, N_Vector, N_Vector);
void        N_VAbs_SpcParallel(N_Vector, N_Vector);
void        N_VInv_SpcParallel(N_Vector, N_Vector);
void        N_VAddConst_SpcParallel(N_Vector, realtype, N_Vector);
realtype    N_VDotProd_SpcParallel(N_Vector, N_Vector);
realtype    N_VMaxNorm_SpcParallel(N_Vector);
realtype    N_VWrmsNorm_SpcParallel(N_Vector, N_Vector);
realtype    N_VWrmsNormMask_SpcParallel(N_Vector, N_Vector, N_Vector);
realtype    N_VMin_SpcParallel(N_Vector);
realtype    N_VWL2Norm_SpcParallel(N_Vector, N_Vector);
realtype    N_VL1Norm_SpcParallel(N_Vector);
void        N_VCompare_SpcParallel(realtype, N_Vector, N_Vector);
booleantype N_VInvTest_SpcParallel(N_Vector, N_Vector);
booleantype N_VConstrMask_SpcParallel(N_Vector, N_Vector, N_Vector);   
realtype    N_VMinQuotient_SpcParallel(N_Vector, N_Vector);



#ifdef __cplusplus
}
#endif

#endif

