/*
 * -----------------------------------------------------------------
 * $Revision: 1.12 $
 * $Date: 2004-07-22 22:21:59 $
 * ----------------------------------------------------------------- 
 * Programmers: Scott D. Cohen, Alan C. Hindmarsh, and 
 *              Radu Serban, LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2002, The Regents of the University of California
 * Produced at the Lawrence Livermore National Laboratory
 * All rights reserved
 * For details, see sundials/shared/LICENSE
 * -----------------------------------------------------------------
 * This is the header file for a parallel MPI (Message Passing   
 * Interface)implementation of NVECTOR package.                  
 *                                                                 
 * Part I of this file contains declarations which are specific    
 * to the particular parallel implementation.
 *                                                                 
 * Part II of this file defines accessor macros that allow the     
 * user to use efficiently the type N_Vector without making        
 * explicit references to its underlying representation.           
 *                                                                 
 * Part III of this file contains the prototype for the constructor 
 * N_VNew_Parallel, as well as prototypes for the vector kernels which 
 * operate on the parallel N_Vector. These prototypes are unique to this 
 * particular implementation of the vector package.           
 *                                                                 
 * NOTES:                                                          
 *                                                                 
 * The definition of the generic N_Vector structure is in the header 
 * file nvector.h.                               
 *                                                                 
 * The definition of the type realtype is in the header file       
 * sundialstypes.h and it may be changed (at the configuration     
 * stage) according to the user's needs. The sundialstypes.h file  
 * also contains the definition for the type booleantype.          
 *                                                                 
 * N_Vector arguments to arithmetic kernels need not be            
 * distinct. Thus, for example, the call                           
 *         N_VLinearSum_Parallel(a,x,b,y,y);   y <- ax+by            
 * is legal.                                                       
 * -----------------------------------------------------------------
 */

#ifndef included_nvector_parallel_h
#define included_nvector_parallel_h

#include "nvector.h"
#include "sundialstypes.h"
#include "mpi.h"

#ifdef __cplusplus     /* wrapper to enable C++ usage */
extern "C" {
#endif

/*
 * -----------------------------------------------------------------
 * PART I:                                                    
 * Parallel implementation of N_Vector               
 * -----------------------------------------------------------------
 */

/* 
 * Set types real and integer for MPI calls. 
 */

#if defined(SUNDIALS_SINGLE_PRECISION)
#define PVEC_REAL_MPI_TYPE MPI_FLOAT
#else
#define PVEC_REAL_MPI_TYPE MPI_DOUBLE
#endif

#define PVEC_INTEGER_MPI_TYPE MPI_LONG

/* 
 * The parallel implementation of the N_Vector 'content' 
 * structure contains the global and local lengths of the vector,
 * a pointer to an array of real components, and the MPI communicator.
 */

struct _N_VectorContent_Parallel {
  long int local_length;         /* local vector length         */
  long int global_length;        /* global vector length        */
  realtype   *data;              /* local data array            */
  MPI_Comm comm;                 /* pointer to MPI communicator */
};

typedef struct _N_VectorContent_Parallel *N_VectorContent_Parallel;

/*
 * -----------------------------------------------------------------
 *                                                            
 * PART II: Macros         
 * NV_SPEC_P, NV_DATA_P, NV_LOCLENGTH_P, NV_GLOBLENGTH_P, NV_Ith_P      
 * -----------------------------------------------------------------
 * In the descriptions below, the following user declarations 
 * are assumed:                                    
 *                                                              
 * N_Vector     v;                                         
 * long int     v_len, s_len, i;                                
 *                                                              
 * (1) NV_CONTENT_P                               
 *                                                              
 *     This routines gives access to the contents of the parallel
 *     vector N_Vector.
 *                                                              
 *     The assignment v_cont = NV_CONTENT_P(v) sets v_cont to be 
 *     a pointer to the parallel N_Vector content structure.
 *
 * (2) NV_DATA_P, NV_LOCLENGTH_P, NV_GLOBLENGTH_P               
 *                                                              
 *     These routines give individual access to the parts of    
 *     the content of a parallel N_Vector.                      
 *                                                              
 *     The assignment v_data=NV_DATA_P(v) sets v_data to be     
 *     a pointer to the first component of the local data for   
 *     the vector v. The assignment NV_DATA_P(v)=v_data sets    
 *     the component array of v to be v_data by storing the     
 *     pointer v_data.
 *                                                              
 *     The assignment v_llen=NV_LOCLENGTH_P(v) sets v_llen to   
 *     be the length of the local part of the vector v.         
 *     The call NV_LOCLENGTH_P(v)=llen_v sets the local length  
 *     of v to be llen_v.                                       
 *                                                              
 *     The assignment v_glen=NV_GLOBLENGTH_P(v) sets v_glen to  
 *     be the global length of the vector v.                    
 *     The call NV_GLOBLENGTH_P(v)=glen_v sets the global       
 *     length of v to be glen_v.                                
 *                                                              
 * (3) NV_Ith_P                                                 
 *                                                              
 *     In the following description, the components of the      
 *     local part of an N_Vector are numbered 0..n-1, where n   
 *     is the local length of (the local part of) v.            
 *                                                              
 *     The assignment r=NV_Ith_P(v,i) sets r to be the value    
 *     of the ith component of the local part of the vector v.  
 *     The assignment NV_Ith_P(v,i)=r sets the value of the     
 *     ith local component of v to be r.                        
 *                                                              
 * When looping over the components of an N_Vector v, it is     
 * more efficient to first obtain the component array via       
 * v_data=NV_DATA_P(v) and then access v_data[i] within the     
 * loop than it is to use NV_Ith_P(v,i) within the loop.        
 *                                                              
 * -----------------------------------------------------------------
 */                                                              

#define NV_CONTENT_P(v)    ( (N_VectorContent_Parallel)(v->content) )

#define NV_LOCLENGTH_P(v)  ( NV_CONTENT_P(v)->local_length )

#define NV_GLOBLENGTH_P(v) ( NV_CONTENT_P(v)->global_length )

#define NV_DATA_P(v)       ( NV_CONTENT_P(v)->data )

#define NV_COMM_P(v)       ( NV_CONTENT_P(v)->comm )

#define NV_Ith_P(v,i)      ( NV_DATA_P(v)[i] )
/*
 * -----------------------------------------------------------------
 * PART III: 
 * Functions exported by nvector_parallel
 * -----------------------------------------------------------------
 */

/*
 * N_VNew_Parallel
 * 
 * This function creates and allocates memory for a parallel vector.
 */

N_Vector N_VNew_Parallel(MPI_Comm comm, 
                         long int local_length,
                         long int global_length);

/*
 * N_VNewEmpty_Parallel
 *
 * This function creates a new parallel N_Vector with an empty (NULL)
 * data array.
 */

N_Vector N_VNewEmpty_Parallel(MPI_Comm comm, 
                              long int local_length,
                              long int global_length);

/*
 * N_VMake_Parallel
 * 
 * This function creates and allocates memory for a parallel vector
 * with user-provided data array.
 */

N_Vector N_VMake_Parallel(MPI_Comm comm, 
                          long int local_length,
                          long int global_length,
                          realtype *v_data);

/*
 * N_VNewVectorArray_Parallel
 *
 * This function creates an array of 'count' parallel vectors.
 * This array of N_Vectors can be freed with N_VDestroyVectorArray
 * (defined by the generic nvector module)
 */

N_Vector *N_VNewVectorArray_Parallel(int count, 
                                     MPI_Comm comm, 
                                     long int local_length,
                                     long int global_length);

/*
 * N_VDispose_Parallel
 *
 * This function frees an N_Vector created with N_VMake_Parallel.
 * Note that deallocation of the 'data' array is the user's
 * responsibility. In other words, N_VDispose_Parallel is identitical
 * to N_VDestroyEmpty_Parallel.
 */

void N_VDispose_Parallel(N_Vector v);

/*
 * N_VPrint_Parallel
 * 
 * This function prints the content of a parallel vector to stdout.
 */

void N_VPrint_Parallel(N_Vector v);

/*
 * Parallel implementations of the vector operations
 */

N_Vector N_VClone_Parallel(N_Vector w);
void N_VDestroy_Parallel(N_Vector v);
N_Vector N_VCloneEmpty_Parallel(N_Vector w);
void N_VDestroyEmpty_Parallel(N_Vector v);
void N_VSpace_Parallel(N_Vector v, long int *lrw, long int *liw);
realtype *N_VGetArrayPointer_Parallel(N_Vector v);
void N_VSetArrayPointer_Parallel(realtype *v_data, N_Vector v);
void N_VLinearSum_Parallel(realtype a, N_Vector x, realtype b, N_Vector y, N_Vector z);
void N_VConst_Parallel(realtype c, N_Vector z);
void N_VProd_Parallel(N_Vector x, N_Vector y, N_Vector z);
void N_VDiv_Parallel(N_Vector x, N_Vector y, N_Vector z);
void N_VScale_Parallel(realtype c, N_Vector x, N_Vector z);
void N_VAbs_Parallel(N_Vector x, N_Vector z);
void N_VInv_Parallel(N_Vector x, N_Vector z);
void N_VAddConst_Parallel(N_Vector x, realtype b, N_Vector z);
realtype N_VDotProd_Parallel(N_Vector x, N_Vector y);
realtype N_VMaxNorm_Parallel(N_Vector x);
realtype N_VWrmsNorm_Parallel(N_Vector x, N_Vector w);
realtype N_VWrmsNormMask_Parallel(N_Vector x, N_Vector w, N_Vector id);
realtype N_VMin_Parallel(N_Vector x);
realtype N_VWL2Norm_Parallel(N_Vector x, N_Vector w);
realtype N_VL1Norm_Parallel(N_Vector x);
void N_VCompare_Parallel(realtype c, N_Vector x, N_Vector z);
booleantype N_VInvTest_Parallel(N_Vector x, N_Vector z);
booleantype N_VConstrMask_Parallel(N_Vector c, N_Vector x, N_Vector m);
realtype N_VMinQuotient_Parallel(N_Vector num, N_Vector denom);


#ifdef __cplusplus
}
#endif

#endif
