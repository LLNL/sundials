/*******************************************************************
 *                                                                 *
 * File          : nvector_parallel.h                              *
 * Programmers   : Scott D. Cohen, Alan C. Hindmarsh,              *
 *               : Radu Serban, and Allan G. Taylor, LLNL          *
 * Version of    : 07 February 2004                                *
 *-----------------------------------------------------------------*
 * Copyright (c) 2002, The Regents of the University of California *
 * Produced at the Lawrence Livermore National Laboratory          *
 * All rights reserved                                             *
 * For details, see sundials/shared/LICENSE                        *
 *-----------------------------------------------------------------*
 * This is the header file for a parallel MPI (Message Passing     *
 * Interface)implementation of NVECTOR package.                    *
 *                                                                 *
 * Part I of this file contains declarations which are specific    *
 * to the particular vector specification in which this version    *
 * of the NVECTOR module is to be used. This includes the          *
 * typedef for the 'content' fields of the structures NV_Spec and  *
 * N_Vector (NV_SpecContent_Serial and N_VectorContent_Parallel,   *
 * respectively).                                                  *
 *                                                                 *
 * Part II of this file defines accessor macros that allow the     *
 * user to use efficiently the type N_Vector without making        *
 * explicit references to its underlying representation.           *
 *                                                                 *
 * Part III of this file contains the prototype for the            *
 * initialization routine specific to this implementation          *
 * (NV_SpecInit_Parallel) as well as prototypes for the vector     *
 * kernels which operate on the parallel N_Vector. These           *
 * prototypes are unique to this particular implementation of      *
 * the vector package.                                             *
 *                                                                 *
 * NOTES:                                                          *
 *                                                                 *
 * The definitions of the generic NV_Spec and N_Vector structures  *
 * are in the header file nvector.h.                               *
 *                                                                 *
 * The definition of the type realtype is in the header file       *
 * sundialstypes.h and it may be changed (at the configuration     *
 * stage) according to the user's needs. The sundialstypes.h file  *
 * also contains the definition for the type booleantype.          *
 *                                                                 *
 * N_Vector arguments to arithmetic kernels need not be            *
 * distinct. Thus, for example, the call                           *
 *         N_VLinearSum_Serial(a,x,b,y,y);   y <- ax+by            *
 * is legal.                                                       *
 *                                                                 * 
 *******************************************************************/


#ifndef included_nvector_parallel_h
#define included_nvector_parallel_h

#include "nvector.h"
#include "sundialstypes.h"
#include "mpi.h"

#ifdef __cplusplus     /* wrapper to enable C++ usage */
extern "C" {
#endif

/****************************************************************
 * PART I:                                                      *
 * Parallel MPI implementaion of NV_Spec and N_Vector           *
 ****************************************************************/

/* Environment: MPI                          */
/* Set types real and integer for MPI calls. */

#if defined(SUNDIALS_SINGLE_PRECISION)
#define PVEC_REAL_MPI_TYPE MPI_FLOAT
#else
#define PVEC_REAL_MPI_TYPE MPI_DOUBLE
#endif

#define PVEC_INTEGER_MPI_TYPE MPI_LONG


/* The parallel implementation of the machine environment has 
   ID tag 'parallel' */
#define ID_TAG_P "parallel"

/* The parallel implementation of the vector specification 'content'
   structure contains the global and local lengths of vectors, a
   pointer to MPI communicator, and a flag showing if the user
   called MPI_Init */

struct _NV_SpecContent_Parallel {
  MPI_Comm comm;                 /* pointer to MPI communicator */
  long int local_vec_length;  /* local length of vectors */ 
  long int global_vec_length; /* global length of vectors */ 
  int init_by_user;              /* flag showing if user called MPI_Init */
};

typedef struct _NV_SpecContent_Parallel *NV_SpecContent_Parallel;
 
/* The parallel implementation of the N_Vector 'content' 
   structure contains the global and local lengths of the vector 
   and a pointer to an array of real components */

struct _N_VectorContent_Parallel {
  long int local_length;  /* local vector length  */
  long int global_length; /* global vector length */
  realtype   *data;          /* local data array     */
};

typedef struct _N_VectorContent_Parallel *N_VectorContent_Parallel;


/****************************************************************
 *                                                              *
 * PART II: Macros                                              *
 *    NV_MAKE_P, NV_DISPOSE_P, NVS_MAKE_P, NVS_DISPOSE_P        *
 *    NS_CONTENT_P, NV_CONTENT_P                                *
 *    NV_DATA_P, NV_LOCLENGTH_P, NV_GLOBLENGTH_P, NV_Ith_P      *
 *--------------------------------------------------------------*
 * In the descriptions below, the following user                *
 * declarations are assumed:                                    *
 *                                                              *
 * NV_Spec      nvspec;                                         *
 * N_Vector     v, *vs;                                         *
 * realtype     *v_data, **vs_data, r;                          *
 * long int  v_len, s_len, i;                                *
 *                                                              *
 * (1) NV_MAKE_P, NV_DISPOSE_P                                  *
 *                                                              *
 *     These companion routines are used to create and          *
 *     destroy an N_Vector with a component array v_data        *
 *     allocated by the user.                                   *
 *                                                              *
 *     The call NV_MAKE_P(v, v_data, nvspec) makes v an         *
 *     N_Vector with component array v_data.  The local and     *
 *     global vector lengths are taken from nvspec.             *
 *     NV_MAKE_P stores the pointer v_data so that              *
 *     changes made by the user to the elements of v_data are   *
 *     simultaneously reflected in v. There is no copying of    *
 *     elements.                                                *
 *                                                              *
 *     The call NV_DISPOSE_P(v) frees all memory associated     *
 *     with v except for its component array. This memory was   *
 *     allocated by the user and, therefore, should be          *
 *     deallocated by the user.                                 *
 *                                                              *
 * (2) NVS_MAKE_P, NVS_DISPOSE_P                                *
 *                                                              *
 *     These companion routines are used to create and          *
 *     destroy an array of N_Vectors with component vs_data     *
 *     allocated by the user.                                   *
 *                                                              * 
 *     The call NVS_MAKE_P(vs, vs_data, s_len, nvspec) makes    *
 *     vs an array of s_len N_Vectors, each with component      *
 *     array vs_data[i] and local and global vector lengths     *
 *     taken from nvspec.                                       *
 *     NVS_MAKE_P stores the pointers vs_data[i] so that        *
 *     changes made by the user to the elements of sdata are    *
 *     simultaneously reflected in vs. There is no copying of   *
 *     elements.                                                *
 *                                                              *
 *     The call NVS_DISPOSE_P(vs) frees all memory associated   *
 *     with vs except for its components' component array.      *
 *     This memory was allocated by the user and, therefore,    *
 *     should be deallocated by the user.                       *
 *                                                              *
 * (3) NS_CONTENT_P, NV_CONTENT_P                               *
 *                                                              *
 *     These routines give access to the contents of the        *
 *     parallel vector specification and N_Vector, respectively.* 
 *                                                              *
 *     The assignment ns_cont = NS_CONTENT_P(nvspec) sets       *
 *     ns_cont to be a pointer to the parallel vector           *
 *     specification content structure.                         * 
 *                                                              *
 *     The assignment v_cont = NV_CONTENT_P(v) sets             *
 *     v_cont to be a pointer to the parallel N_Vector content  *
 *     structure.                                               *
 *                                                              *
 * (4) NV_DATA_P, NV_LOCLENGTH_P, NV_GLOBLENGTH_P               *
 *                                                              *
 *     These routines give individual access to the parts of    *
 *     the content of a parallel N_Vector.                      *
 *                                                              *
 *     The assignment v_data=NV_DATA_P(v) sets v_data to be     *
 *     a pointer to the first component of the local data for   *
 *     the vector v. The assignment NV_DATA_P(v)=v_data sets    *
 *     the component array of v to be v_data by storing the     *
 *     pointer v_data.                                          *  
 *                                                              *
 *     The assignment v_llen=NV_LOCLENGTH_P(v) sets v_llen to   *
 *     be the length of the local part of the vector v.         *
 *     The call NV_LOCLENGTH_P(v)=llen_v sets the local length  *
 *     of v to be llen_v.                                       *
 *                                                              *
 *     The assignment v_glen=NV_GLOBLENGTH_P(v) sets v_glen to  *
 *     be the global length of the vector v.                    *
 *     The call NV_GLOBLENGTH_P(v)=glen_v sets the global       *
 *     length of v to be glen_v.                                *
 *                                                              *
 * (4) NV_Ith_P                                                 *
 *                                                              *
 *     In the following description, the components of the      *
 *     local part of an N_Vector are numbered 0..n-1, where n   *
 *     is the local length of (the local part of) v.            *
 *                                                              *
 *     The assignment r=NV_Ith_P(v,i) sets r to be the value    *
 *     of the ith component of the local part of the vector v.  *
 *     The assignment NV_Ith_P(v,i)=r sets the value of the     *
 *     ith local component of v to be r.                        *
 *                                                              *
 * Notes..                                                      *
 *                                                              *
 * Users who use the macros (1) and/or (2) must                 *
 * #include<stdlib.h> since these macros expand to calls to     *
 * malloc and free.                                             *
 *                                                              *
 * When looping over the components of an N_Vector v, it is     *
 * more efficient to first obtain the component array via       *
 * v_data=NV_DATA_P(v) and then access v_data[i] within the     *
 * loop than it is to use NV_Ith_P(v,i) within the loop.        *
 *                                                              *
 * NV_MAKE_P and NV_DISPOSE_P are similar to N_VNew_Parallel    *
 * and N_VFree_Parallel, while NVS_MAKE_P and NVS_DISPOSE_P     *
 * are similar to  N_VNew_S_Parallel and N_VFree_S_Parallel.    *
 * The difference is one of responsibility for component        *
 * memory allocation and deallocation. N_VNew_Parallel          *
 * allocates memory for the N_Vector components and             *
 * N_VFree_Parallel frees the component memory allocated by     *
 * N_VNew_Parallel. For NV_MAKE_P and NV_DISPOSE_P, the         *
 * component memory is allocated and freed by the user of this  *
 * package. Similar remarks hold for NVS_MAKE_P,                *
 * NVS_DISPOSE_P and N_VNew_S_Parallel, N_VFree_S_Parallel.     *
 *                                                              *
 ****************************************************************/ 

#define NV_MAKE_P(v, v_data, nvspec) \
        v = (N_Vector) malloc(sizeof(*v)); \
        v->content = (N_VectorContent_Parallel) malloc(sizeof(struct _N_VectorContent_Parallel)); \
        v->content->data = v_data; \
        v->content->local_length  = nvspec->content->local_vec_length; \
        v->content->global_length = nvspec->content->global_vec_length; \
        v->nvspec = nvspec

#define NV_DISPOSE_P(v) \
        free((N_VectorContent_Parallel)(v->content)); \
        free(v)

#define NVS_MAKE_P(vs, vs_data, s_len, nvspec) \
        vs = (N_Vector_S) malloc(s_len*sizeof(N_Vector *)); \
        for ((int)is=0; is<s_len; is++) { \
           NV_MAKE_P(vs[is], vs_data[is], nvspec); \
        }
#define NVS_DISPOSE_P(vs, s_len) \
        for ((int)is=0; is<s_len; is++) NV_DISPOSE_P(vs[i]); \
        free(vs);

#define NS_CONTENT_P(s) ( (NV_SpecContent_Parallel)(s->content) )

#define NV_CONTENT_P(v) ( (N_VectorContent_Parallel)(v->content) )

#define NV_LOCLENGTH_P(v) ( NV_CONTENT_P(v)->local_length )

#define NV_GLOBLENGTH_P(v) ( NV_CONTENT_P(v)->global_length )

#define NV_DATA_P(v) ( NV_CONTENT_P(v)->data )

#define NV_Ith_P(v,i) ( NV_DATA_P(v)[i] )


/****************************************************************
 * PART III:                                                    *
 * Functions exported by nvector_serial                         *
 ****************************************************************/

/*--------------------------------------------------------------*
 * Routine : NV_SpecInit_Parallel                               *
 *--------------------------------------------------------------*
 * This function sets the content field of the vector           *
 * specification for the parallel MPI implementation to a       *
 * structure of type _NV_SpecContent_Parallel and attaches the  *
 * vector operations defined for this implementation.           *
 *                                                              *
 * If successful, NV_SpecInit_Parallel returns a pointer of type*
 * NV_Spec. This pointer should in turn be passed in any user   *
 * calls to N_VNew, or uses of the macros NV_MAKE_P and         *
 * NVS_MAKE_P.                                                  *
 * If a memory allocation failure occurs, or if the global      *
 * length differs from the sum of the local lengths,            *
 * NV_SpecInit_Parallel returns NULL.  In the latter case, an   *
 * error message is printed to stdout.                          *
 *                                                              *
 *--------------------------------------------------------------*
 *                                                              *
 * comm              is a pointer to the MPI communicator,      *
 *                   of type MPI_Comm.  Must be non-NULL.       *
 *                                                              *
 * local_vec_length  is the length of the piece of the vectors  *
 *                   residing on this processor.                *
 *                   If the active processor set is a proper    *
 *                   subset of the full processor set assigned  *
 *                   to the job, the value of local_vec_length  *
 *                   should be 0 on the inactive processors.    *
 *                   (Otherwise, the two global length values   *
 *                   input and computed, may differ.)           *
 *                                                              *
 * global_vec_length is the global length of the vectors.       *
 *                   This must equal the sum of all local       *
 *                   lengths over the active processor set.     *
 *                   If not, a message is printed.              *
 *                                                              *
 * argc              is the command line arguments count from   *
 *                   the main program (or, a dummy if MPI_INIT  *
 *                   has already been called).                  *
 *                                                              *
 * argv              is the command line argument character     *
 *                   array from the main program (or, a dummy   *
 *                   if MPI_INIT has already been called)       *
 *                                                              *
 *--------------------------------------------------------------*/

NV_Spec NV_SpecInit_Parallel(MPI_Comm comm, 
                             long int local_vec_length,
                             long int global_vec_length, 
                             int *argc, char ***argv);

/*----------------------------------------------------------------*
 * Function NV_SpecFree_Parallel                                  *
 *----------------------------------------------------------------*
 * Function to free the block of machine-dependent environment    *
 * information created by NV_SpecInit_Parallel.                   *
 * Its only argument is the pointer nvspec returned by            *
 * NV_SpecInit_Parallel.                                          *
 * NOTE: if MPI is initialized by other than NV_SpecInit_Parallel *
 * it is necessary to call MPI_Finalize in addition to (after)    *
 * calling NV_SpecFree_Parallel.                                  *
 *----------------------------------------------------------------*/

void NV_SpecFree_Parallel(NV_Spec nvspec);

/*--------------------------------------------------------------*
 * Parallel implementations of the vector operations            *
 *                                                              *
 * For a complete description of each of the following routines *
 * see the header file nvector.h                                *
 *--------------------------------------------------------------*/

N_Vector N_VNew_Parallel(NV_Spec nvspec);
void N_VSpace_Parallel(NV_Spec nvspec, long int *lrw, long int *liw);
void N_VFree_Parallel(N_Vector v);
N_Vector N_VMake_Parallel(realtype *v_data, NV_Spec nvspec);
void N_VDispose_Parallel(N_Vector v);
realtype *N_VGetData_Parallel(N_Vector v);
void N_VSetData_Parallel(realtype *v_data, N_Vector v);
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
booleantype N_VConstrProdPos_Parallel(N_Vector c, N_Vector x);
booleantype N_VConstrMask_Parallel(N_Vector c, N_Vector x, N_Vector m);   
realtype N_VMinQuotient_Parallel(N_Vector num, N_Vector denom);
void N_VPrint_Parallel(N_Vector x);

#ifdef __cplusplus
}
#endif

#endif
