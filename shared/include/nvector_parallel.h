/*******************************************************************
 *                                                                 *
 * File          : nvector_parallel.h                              *
 * Programmers   : Scott D. Cohen, Alan C. Hindmarsh,              *
 *               : Radu Serban, and Allan G. Taylor, LLNL          *
 * Version of    : 26 June 2002                                    *
 *-----------------------------------------------------------------*
 * Copyright (c) 2002, The Regents of the University of California *
 * Produced at the Lawrence Livermore National Laboratory          *
 * All rights reserved                                             *
 * For details, see sundials/shared/LICENSE                        *
 *-----------------------------------------------------------------*
 * This is the header file for a parallel MPI implementation of    *
 * NVECTOR package.                                                *
 *                                                                 *
 * Part I of this file contains declarations which are specific    *
 * to the particular machine environment in which this version     *
 * of the NVECTOR module is to be used. This includes the          *
 * typedef for the 'content' fields of the structures M_Env and    *
 * N_Vector (M_EnvParallelContent and N_VectorParallelContent,     *
 * respectively).                                                  *
 *                                                                 *
 * Part II of this file defines accessor macros that allow the     *
 * user to use efficiently the type N_Vector without making        *
 * explicit references to its underlying representation.           *
 *                                                                 *
 * Part III of this file contains the prototype for the            *
 * initialization routine specific to this implementation          *
 * (M_EnvInit_Parallel) as well as prototypes for the vector       *
 * kernels which operate on the parallel N_Vector. These           *
 * prototypes are unique to this particular implementation of      *
 * the vector package.                                             *
 *                                                                 *
 * NOTES:                                                          *
 *                                                                 *
 * The definitions of the generic M_Env and N_Vector structures    *
 * are in the header file nvector.h.                               *
 *                                                                 *
 * The definitions of the types realtype and integertype are in    *
 * the header file sundialstypes.h and these may be changed        *
 * according to the user's needs. The sundialstypes.h file also    *
 * contains the definition for the type booleantype.               *
 *                                                                 *
 * N_Vector arguments to arithmetic kernels need not be            *
 * distinct. Thus, for example, the call                           *
 *         N_VLinearSum_Serial(a,x,b,y,y);   y <- ax+by            *
 * is legal.                                                       *
 *                                                                 * 
 * This version of nvector is for the MPI (Message Passing         *
 * Interface) machine environment. In the documentation given      *
 * below, N is the local length of all N_Vector parameters and     *
 * x[i] denotes the ith component of the local part of the         *
 * distributed N_Vector x,  where 0 <= i <= N-1.                   *
 *                                                                 *
 *******************************************************************/


#ifdef __cplusplus     /* wrapper to enable C++ usage */
extern "C" {
#endif

#ifndef included_nvector_parallel_h
#define included_nvector_parallel_h

#include "nvector.h"  /* Generic M_Env and N_Vector type definitions */
#include "sundialstypes.h"
#include "mpi.h"

/****************************************************************
 * PART I:                                                      *
 * Parallel MPI implementaion of M_Env and N_Vector             *
 ****************************************************************/

/* Environment: MPI                                  */
/* Set types realtype and integertype for MPI calls. */

#if (SUNDIALS_DOUBLE == 1)
#define PVEC_REAL_MPI_TYPE MPI_DOUBLE
#else
#define PVEC_REAL_MPI_TYPE MPI_FLOAT
#endif

#if (SUNDIALS_INT == 1)
#define PVEC_INTEGER_MPI_TYPE MPI_INT
#else
#define PVEC_INTEGER_MPI_TYPE MPI_LONG
#endif

/* The parallel implementation of the machine environment has 
   ID tag 'parallel' */
#define ID_TAG_P "parallel"

/* The parallel implementation of the machine environment 'content'
   structure contains the global and local lengths of vectors, a
   pointer to MPI communicator, and a flag showing if the user
   called MPI_Init */

struct _M_EnvParallelContent {
  MPI_Comm comm;                 /* pointer to MPI communicator */
  integertype local_vec_length;  /* local length of vectors */ 
  integertype global_vec_length; /* global length of vectors */ 
  int init_by_user;              /* flag showing if user called MPI_Init */
};

typedef struct _M_EnvParallelContent *M_EnvParallelContent;
 
/* The parallel implementation of the N_Vector 'content' 
   structure contains the global and local lengths of the vector 
   and a pointer to an array of real components */

struct _N_VectorParallelContent {
  integertype local_length;  /* local vector length  */
  integertype global_length; /* global vector length */
  realtype   *data;          /* local data array     */
};

typedef struct _N_VectorParallelContent *N_VectorParallelContent;


/****************************************************************
 *                                                              *
 * PART II: Macros                                              *
 *    NV_MAKE_P, NV_DISPOSE_P, NVS_MAKE_P, NVS_DISPOSE_P        *
 *    ME_CONTENT_P, NV_CONTENT_P                                *
 *    NV_DATA_P, NV_LOCLENGTH_P, NV_GLOBLENGTH_P, NV_Ith_P      *
 *--------------------------------------------------------------*
 * In the descriptions below, the following user                *
 * declarations are assumed:                                    *
 *                                                              *
 * M_Env        machenv;                                        *
 * N_Vector     v, *vs;                                         *
 * realtype     *v_data, **vs_data, r;                          *
 * integertype  v_len, s_len, i;                                *
 *                                                              *
 * (1) NV_MAKE_P, NV_DISPOSE_P                                  *
 *                                                              *
 *     These companion routines are used to create and          *
 *     destroy an N_Vector with a component array v_data        *
 *     allocated by the user.                                   *
 *                                                              *
 *     The call NV_MAKE_P(v, v_data, machenv) makes v an        *
 *     N_Vector with component array v_data.  The local and     *
 *     global vector lengths are taken from machenv.            *
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
 *     The call NVS_MAKE_P(vs, vs_data, s_len, machEnv) makes   *
 *     vs an array of s_len N_Vectors, each with component      *
 *     array vs_data[i] and local and global vector lengths     *
 *     taken from machenv.                                      *
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
 * (3) ME_CONTENT_P, NV_CONTENT_P                               *
 *                                                              *
 *     These routines give access to the contents of the        *
 *     parallel machine environment and N_Vector, respectively. * 
 *                                                              *
 *     The assignment m_cont = ME_CONTENT_P(machenv) sets       *
 *     m_cont to be a pointer to the parallel machine           *
 *     environment content structure.                           * 
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

#define NV_MAKE_P(v, v_data, machenv) \
        v = (N_Vector) malloc(sizeof(*v)); \
        v->content = (N_VectorParallelContent) malloc(sizeof(struct _N_VectorParallelContent)); \
        v->content->data = v_data; \
        v->content->local_length  = machenv->content->local_vec_length; \
        v->content->global_length = machenv->content->global_vec_length; \
        v->menv = machenv

#define NV_DISPOSE_P(v) \
        free((N_VectorParallelContent)(v->content)); \
        free(v)

#define NVS_MAKE_P(vs, vs_data, s_len, machenv) \
        vs = (N_Vector_S) malloc(s_len*sizeof(N_Vector *)); \
        for ((int)is=0; is<s_len; is++) { \
           NV_MAKE_P(vs[is], vs_data[is], machenv); \
        }
#define NVS_DISPOSE_P(vs, s_len) \
        for ((int)is=0; is<s_len; is++) NV_DISPOSE_P(vs[i]); \
        free(vs);

#define ME_CONTENT_P(m) ( (M_EnvParallelContent)(m->content) )

#define NV_CONTENT_P(v) ( (N_VectorParallelContent)(v->content) )

#define NV_LOCLENGTH_P(v) ( NV_CONTENT_P(v)->local_length )

#define NV_GLOBLENGTH_P(v) ( NV_CONTENT_P(v)->global_length )

#define NV_DATA_P(v) ( NV_CONTENT_P(v)->data )

#define NV_Ith_P(v,i) ( NV_DATA_P(v)[i] )


/****************************************************************
 * PART III:                                                    *
 * Functions exported by nvector_serial                         *
 ****************************************************************/

/*--------------------------------------------------------------*
 * Routine : M_EnvInit_Parallel                                 *
 *--------------------------------------------------------------*
 * This function sets the content field of the machine          *
 * environment for the parallel MPI implementation to a         *
 * structure of type _MEnvParallelContent and attaches the      *
 * vector operations defined for this implementation.           *
 *                                                              *
 * If successful, M_EnvInit_Parallel returns a pointer of type  *
 * M_Env. This pointer should in turn be passed in any user     *
 * calls to N_VNew, or uses of the macros NV_MAKE_P and         *
 * NVS_MAKE_P.                                                  *
 * If a memory allocation failure occurs, or if the global      *
 * length differs from the sum of the local lengths,            *
 * M_EnvInit_Parallel returns NULL.  In the latter case, an     *
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

M_Env M_EnvInit_Parallel(MPI_Comm comm, 
                         integertype local_vec_length,
                         integertype global_vec_length, 
                         int *argc, char ***argv);

/*--------------------------------------------------------------*
 * Function M_EnvFree_Parallel                                  *
 *--------------------------------------------------------------*
 * Function to free the block of machine-dependent environment  *
 * information created by N_VecInit_Parallel.                   *
 * Its only argument is the pointer machenv returned by         *
 * M_EnvInit_Parallel.                                          *
 * NOTE: if MPI is initialized by other than M_EnvInit_Parallel *
 * it is necessary to call MPI_Finalize in addition to (after)  *
 * calling M_EnvFree_Parallel.                                  *
 *--------------------------------------------------------------*/

void M_EnvFree_Parallel(M_Env machenv);

/*--------------------------------------------------------------*
 * Parallel implementations of the vector operations            *
 *                                                              *
 * For a complete description of each of the following routines *
 * see the header file nvector.h                                *
 *--------------------------------------------------------------*/

N_Vector N_VNew_Parallel(integertype n, M_Env machEnv);
N_Vector_S N_VNew_S_Parallel(integertype ns, integertype n, M_Env machEnv);
void N_VFree_Parallel(N_Vector v);
void N_VFree_S_Parallel(integertype ns, N_Vector_S vs);
N_Vector N_VMake_Parallel(integertype n, realtype *v_data, M_Env machEnv);
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
realtype N_VMin_Parallel(N_Vector x);
realtype N_VWL2Norm_Parallel(N_Vector x, N_Vector w);
realtype N_VL1Norm_Parallel(N_Vector x);
void N_VOneMask_Parallel(N_Vector x);
void N_VCompare_Parallel(realtype c, N_Vector x, N_Vector z);
booleantype N_VInvTest_Parallel(N_Vector x, N_Vector z);
booleantype N_VConstrProdPos_Parallel(N_Vector c, N_Vector x);
booleantype N_VConstrMask_Parallel(N_Vector c, N_Vector x, N_Vector m);   
realtype N_VMinQuotient_Parallel(N_Vector num, N_Vector denom);
void N_VPrint_Parallel(N_Vector x);


#endif
#ifdef __cplusplus
}
#endif
