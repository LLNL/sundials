/*
 * -----------------------------------------------------------------
 * $Revision: 1.10 $
 * $Date: 2004-07-22 21:10:21 $
 * ----------------------------------------------------------------- 
 * Programmers: Scott D. Cohen, Alan C. Hindmarsh, and 
 *              Radu Serban, LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2002, The Regents of the University of California
 * Produced at the Lawrence Livermore National Laboratory
 * All rights reserved
 * For details, see sundials/shared/LICENSE
 * -----------------------------------------------------------------
 * This is the header file for a serial implementation of the    
 * NVECTOR package.                                              
 *                                                               
 * Part I of this file contains declarations which are specific  
 * to the particular serial implementation.                      
 *                                                               
 * Part II of this file defines accessor macros that allow the   
 * user to use efficiently the type N_Vector without making      
 * explicit references to its underlying representation.         
 *                                                               
 * Part III of this file contains the prototype for the constructor 
 * N_VNew_Serial, as well as prototypes for the vector kernels which 
 * operate on the serial N_Vector. These prototypes are unique to this 
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
 *         N_VLinearSum_Serial(a,x,b,y,y);   y <- ax+by          
 * is legal.                                                     
 * -----------------------------------------------------------------
 */

#ifndef included_nvector_serial_h
#define included_nvector_serial_h

#include "nvector.h"
#include "sundialstypes.h"

#ifdef __cplusplus     /* wrapper to enable C++ usage */
extern "C" {
#endif

/*
 * -----------------------------------------------------------------
 * PART I:                                                    
 * Serial implementation of N_Vector               
 * -----------------------------------------------------------------
 */

/*
 * The serial implementation of the N_Vector 'content' structure 
 * contains the length of the vector and a pointer to an array of 
 * realtype components 
 */

struct _N_VectorContent_Serial {
  long int length;
  realtype   *data;
};

typedef struct _N_VectorContent_Serial *N_VectorContent_Serial;

/*
 * -----------------------------------------------------------------
 *                                                            
 * PART II: Macros                                            
 * NV_CONTENT_S, NV_DATA_S, NV_LENGTH_S, NV_Ith_S
 * -----------------------------------------------------------------
 * In the descriptions below, the following user declarations 
 * are assumed:                                               
 *                                                            
 * N_Vector  v;
 * long int  i;
 *                                                            
 * (1) NV_CONTENT_S                             
 *                                                            
 *     This routines gives access to the contents of the serial
 *     vector N_Vector.
 *                                                            
 *     The assignment v_cont = NV_CONTENT_S(v) sets           
 *     v_cont to be a pointer to the serial N_Vector content  
 *     structure.                                             
 *                                                            
 * (2) NV_DATA_S, NV_LENGTH_S                                 
 *                                                            
 *     These routines give individual access to the parts of  
 *     the content of a serial N_Vector.                      
 *                                                            
 *     The assignment v_data=NV_DATA_S(v) sets v_data to be   
 *     a pointer to the first component of v. The assignment  
 *     NV_DATA_S(v)=v_data sets the component array of v to   
 *     be v_data by storing the pointer v_data.
 *                                                            
 *     The assignment v_len=NV_LENGTH_S(v) sets v_len to be   
 *     the length of v. The call NV_LENGTH_S(v)=len_v sets    
 *     the length of v to be len_v.                           
 *                                                            
 * (3) NV_Ith_S                                               
 *                                                            
 *     In the following description, the components of an     
 *     N_Vector are numbered 0..N-1, where N is the length of v.
 *                                                            
 *     The assignment r=NV_Ith_S(v,i) sets r to be the value of
 *     the ith component of v. The assignment NV_Ith_S(v,i)=r 
 *     sets the value of the ith component of v to be r.      
 *                                                            
 * When looping over the components of an N_Vector v, it is   
 * more efficient to first obtain the component array via     
 * v_data=NV_DATA_S(v) and then access v_data[i] within the   
 * loop than it is to use NV_Ith_S(v,i) within the loop.      
 *                                                            
 * -----------------------------------------------------------------
 */                                                            

#define NV_CONTENT_S(v) ( (N_VectorContent_Serial)(v->content) )

#define NV_LENGTH_S(v) ( NV_CONTENT_S(v)->length )

#define NV_DATA_S(v) ( NV_CONTENT_S(v)->data )

#define NV_Ith_S(v,i) ( NV_DATA_S(v)[i] )

/*
 * -----------------------------------------------------------------
 * PART III: 
 * Functions exported by nvector_serial
 * -----------------------------------------------------------------
 */

/*
 * N_VNew_Serial
 * 
 * This function creates and allocates memory for a serial vector.
 * Its only argument is the vector length.
 */

N_Vector N_VNew_Serial(long int vec_length);

/*
 * N_VNewEmpty_Serial
 *
 * This function creates a new serial N_Vector with an empty (NULL)
 * data array.
 */

N_Vector N_VNewEmpty_Serial(long int vec_length);

/*
 * N_VMake_Serial
 *
 * This function creates and allocates memory for a serial vector
 * with user-provided data array.
 */

N_Vector N_VMake_Serial(long int vec_length, realtype *v_data);

/*
 * N_VCloneEmpty_Serial
 *
 * This function creates a new serial vector of the same type as an 
 * existing vector but with an empty (NULL) data array.
 */

N_Vector N_VCloneEmpty_Serial(N_Vector w);

/*
 * N_VNewVectorArray_Serial
 *
 * This function creates an array of 'count' serial vectors.
 * This array of N_Vectors can be freed with N_VDestroyVectorArray
 * (defined by the generic nvector module)
 */

N_Vector *N_VNewVectorArray_Serial(int count, long int vec_length);

/*
 * N_VDestroyEmpty_Serial
 *
 * This function frees an N_Vector created with N_VNewEmpty_Serial
 * or N_VCloneEmpty_Serial.
 */

void N_VDestroyEmpty_Serial(N_Vector v);

/*
 * N_VDispose_Serial
 *
 * This function frees an N_Vector created with N_VMake_Serial.
 * Note that deallocation of the 'data' array is the user's
 * responsibility. In other words, N_VDispose_Serial is identitical
 * to N_VDestroyEmpty_Serial.
 */

void N_VDispose_Serial(N_Vector v);

/*
 * N_VPrint_Serial
 * 
 * This function prints the content of a serial vector to stdout.
 */

void N_VPrint_Serial(N_Vector v);

/*
 * Serial implementations of the vector operations
 */

N_Vector N_VClone_Serial(N_Vector w);
void N_VDestroy_Serial(N_Vector v);
void N_VSpace_Serial(N_Vector v, long int *lrw, long int *liw);
realtype *N_VGetArrayPointer_Serial(N_Vector v);
void N_VSetArrayPointer_Serial(realtype *v_data, N_Vector v);
void N_VLinearSum_Serial(realtype a, N_Vector x, realtype b, N_Vector y, N_Vector z);
void N_VConst_Serial(realtype c, N_Vector z);
void N_VProd_Serial(N_Vector x, N_Vector y, N_Vector z);
void N_VDiv_Serial(N_Vector x, N_Vector y, N_Vector z);
void N_VScale_Serial(realtype c, N_Vector x, N_Vector z);
void N_VAbs_Serial(N_Vector x, N_Vector z);
void N_VInv_Serial(N_Vector x, N_Vector z);
void N_VAddConst_Serial(N_Vector x, realtype b, N_Vector z);
realtype N_VDotProd_Serial(N_Vector x, N_Vector y);
realtype N_VMaxNorm_Serial(N_Vector x);
realtype N_VWrmsNorm_Serial(N_Vector x, N_Vector w);
realtype N_VWrmsNormMask_Serial(N_Vector x, N_Vector w, N_Vector id);
realtype N_VMin_Serial(N_Vector x);
realtype N_VWL2Norm_Serial(N_Vector x, N_Vector w);
realtype N_VL1Norm_Serial(N_Vector x);
void N_VCompare_Serial(realtype c, N_Vector x, N_Vector z);
booleantype N_VInvTest_Serial(N_Vector x, N_Vector z);
booleantype N_VConstrMask_Serial(N_Vector c, N_Vector x, N_Vector m);
realtype N_VMinQuotient_Serial(N_Vector num, N_Vector denom);

#ifdef __cplusplus
}
#endif

#endif
