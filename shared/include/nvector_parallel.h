/****************************************************************
 *                                                              *
 * File          : nvector.h                                    *
 * Programmers   : Scott D. Cohen, Alan C. Hindmarsh,           *
 *               : Radu Serban, and Allan G. Taylor, LLNL       *
 * Version of    : 20 December 2001                             *
 *--------------------------------------------------------------*
 *                                                              *
 * This is the header file for the MPI NVECTOR module.          *
 * It exports the type N_Vector.                                *
 *                                                              *
 * Part I of this file contains declarations which are specific *
 * to the particular machine environment in which this version  *
 * of the NVECTOR module is to be used. This includes the       *
 * typedef for the type machEnvType (machine environment data   *
 * block), type N_Vector, as well as accessor macros            *
 * that allow the user to use efficiently the type N_Vector     *
 * without making explicit references to its underlying         *
 * representation. The underlying type of N_Vector will always  *
 * be some pointer type.                                        *
 *                                                              *
 * Part II of this file contains the prototypes for the vector  *
 * kernels which operate on the type N_Vector. These prototypes *
 * are fixed for all implementations of the NVECTOR module. The *
 * definitions of the types real and integer are in the header  *
 * file llnltyps.h and these may be changed according to the    *
 * user's needs. The llnltyps.h file also contains the          *
 * definition for the type boole (short for boolean) that is    *
 * the return type for the routine N_VInvTest.                  *
 *                                                              *
 * Important Note: N_Vector arguments to arithmetic kernels     *
 * need not be distinct. Thus, for example, the call            *
 *         N_VLinearSum(a,x,b,y,y);    y <- ax+by               *
 * is legal.                                                    *
 *                                                              * 
 * This version of nvector.h is for the MPI (Message Passing    *
 * Interface) machine environment. In the documentation given   *
 * below, N is the local length of all N_Vector parameters and  *
 * x[i] denotes the ith component of the local part of the      *
 * distributed N_Vector x,  where 0 <= i <= N-1.                *
 *                                                              *
 ****************************************************************/


#ifdef __cplusplus     /* wrapper to enable C++ usage */
extern "C" {
#endif

#ifndef nvector_h
#define nvector_h


#include "llnltyps.h"
#include "mpi.h"


/* Part I: Machine Environment-Dependent Declarations */

/* Environment: MPI      */


/* Set types real and integer for MPI calls. */

#if (LLNL_DOUBLE == 1)
#define PVEC_REAL_MPI_TYPE MPI_DOUBLE
#else
#define PVEC_REAL_MPI_TYPE MPI_FLOAT
#endif

#if (LLNL_INT == 1)
#define PVEC_INTEGER_MPI_TYPE MPI_INT
#else
#define PVEC_INTEGER_MPI_TYPE MPI_LONG
#endif


/***************************************************************
 *                                                             *
 * Function PVecInitMPI                                        *
 *-------------------------------------------------------------*
 * Function to set block of machine-dependent environment      *
 * information.                                                *
 *                                                             *
 * comm              is a pointer to the MPI communicator,     *
 *                   of type MPI_Comm.  Must be non-NULL.      *
 *                                                             *
 * local_vec_length  is the length of the piece of the vectors *
 *                   residing on this processor.               *
 *                   If the active processor set is a proper   *
 *                   subset of the full processor set assigned *
 *                   to the job, the value of local_vec_length *
 *                   should be 0 on the inactive processors.   *
 *                   (Otherwise, the two global length values  *
 *                   input and computed, may differ.)          *
 *                                                             *
 * global_vec_length is the global length of the vectors.      *
 *                   This must equal the sum of all local      *
 *                   lengths over the active processor set.    *
 *                   If not, a message is printed.             *
 *                                                             *
 * argc              is the command line arguments count from  *
 *                   the main program (or, a dummy if MPI_INIT *
 *                   has already been called).                 *
 *                                                             *
 * argv              is the command line argument character    *
 *                   array from the main program (or, a dummy  *
 *                   if MPI_INIT has already been called)      *
 *                                                             *
 * If successful, PVecInitMPI returns a pointer to a block     *
 * of machine-environment information, of type *machEnvType.   *
 * This pointer should in turn be passed in any user calls     *
 * to N_VNew, or uses of the macro N_VMAKE.                    *
 * If a memory allocation failure occurs, or if the global     *
 * length differs from the sum of the local lengths,           *
 * PVecInitMPI returns NULL.  In the latter case, an           *
 * error message is printed to stdout.                         *
 *                                                             *
 ***************************************************************/

void *PVecInitMPI(MPI_Comm comm, integer local_vec_length,
		  integer global_vec_length, int *argc, char ***argv);


/***************************************************************
 *                                                             *
 * Function PVecFreeMPI                                        *
 *-------------------------------------------------------------*
 * Function to free the block of machine-dependent environment *
 * information created by PVecInitMPI.                         *
 * Its only argument is the pointer machEnv returned by        *
 * PVecInitMPI.                                                *
 * NOTE: if MPI is initialized by other than PVecInitMPI, it   *
 * is necessary to call MPI_Finalize in addition to (after)    *
 * calling PVecFreeMPI.                                        *
 ***************************************************************/

void PVecFreeMPI(void *machEnv);


/***************************************************************
 *                                                             *
 * Type: machEnvType                                           *
 *-------------------------------------------------------------*
 * The type machEnvType is a type for the block of             *
 * machine-dependent information required for parallel         *
 * implementations.  In this implementation, blocks of this    *
 * type are created by a user call to PVecInitMPI.             *
 * A structure of this type is a member of the type N_Vector.  *
 ***************************************************************/

typedef struct {
  MPI_Comm comm;             /* pointer to MPI communicator */
  integer local_vec_length;  /* local length of vectors */ 
  integer global_vec_length; /* global length of vectors */ 
  int init_by_user;          /* flag showing if user called MPI_Init */
} *machEnvType;

 
/***************************************************************
 *                                                             *
 * Type: N_Vector                                              *
 *-------------------------------------------------------------*
 * The type N_Vector is an abstract vector type. The fields of *
 * its concrete representation should not be accessed          *
 * directly, but rather through the macros given below.        *
 *                                                             *
 * A user may assume that the N components of an N_Vector      *
 * are stored contiguously. A pointer to the first component   *
 * can be obtained via the macro N_VDATA.                      *
 *                                                             *
 * Machine environment pointer machEnv is for parallel version.*
 *                                                             *
 ***************************************************************/

typedef struct {
  integer length;        /* local vector length */
  integer global_length; /* global vector length */
  real   *data;          /* local data array */
  machEnvType machEnv;   /* machine environment pointer */
} *N_Vector;
 
 
/***************************************************************
 *                                                             *
 * Macros: N_VMAKE, N_VDISPOSE, N_VMAKE_S, N_VDISPOSE_S,       *
 *         N_VDATA, N_VLOCLENGTH, N_VGLOBLENGTH, N_VIth        *
 *-------------------------------------------------------------*
 * In the descriptions below, the following user               *
 * declarations are assumed:                                   *
 *                                                             *
 * N_Vector v; real *v_data, r; integer v_len, i;              *
 *                                                             *
 * (1) N_VMAKE, N_VDISPOSE                                     *
 *                                                             *
 *     These companion routines are used to create and         *
 *     destroy an N_Vector with a component array v_data       *
 *     allocated by the user.                                  *
 *                                                             *
 *     The call N_VMAKE(v, v_data, machEnv) makes v an         *
 *     N_Vector with component array v_data.  The local and    *
 *     global vector lengths are taken from (*machEnv).        *
 *     N_VMAKE stores the pointer v_data so that               *
 *     changes made by the user to the elements of v_data are  *
 *     simultaneously reflected in v. There is no copying of   *
 *     elements.                                               *
 *                                                             *
 *     The call N_VDISPOSE(v) frees all memory associated      *
 *     with v except for its component array. This memory was  *
 *     allocated by the user and, therefore, should be         *
 *     deallocated by the user.                                *
 *                                                             *
 * (2) N_VMAKE_S, N_VDISPOSE_S                                 *
 *                                                             *
 *     These companion routines are used to create and         *
 *     destroy an array of N_Vectors with component sdata      *
 *     allocated by the user.                                  *
 *                                                             *
 *     The call N_VMAKE_S(vs, vs_data, s_len, machEnv) makes   *
 *     vs an array of s_len N_Vectors, each with component     *
 *     array vs_data[i] and local and global vector lengths    *
 *     taken from (*machEnv).                                  *
 *     N_VMAKE_S stores the pointers vs_data[i] so that        *
 *     changes made by the user to the elements of sdata are   *
 *     simultaneously reflected in vs. There is no copying of  *
 *     elements.                                               *
 *                                                             *
 *     The call N_VDISPOSE_S(vs) frees all memory associated   *
 *     with vs except for its components' component array.     *
 *     This memory was allocated by the user and, therefore,   *
 *     should be deallocated by the user.                      *
 *                                                             *
 * (3) N_VDATA, N_VLOCLENGTH, N_VGLOBLENGTH                    *
 *                                                             *
 *     These routines give individual access to the parts of   *
 *     an N_Vector.                                            *
 *                                                             *
 *     The assignment v_data=N_VDATA(v) sets v_data to be      *
 *     a pointer to the first component of the local data for  *
 *     the vector v. The assignment N_VDATA(v)=v_data sets the *
 *     component array of v to be v_data by storing the        *
 *     pointer v_data.                                         *  
 *                                                             *
 *     The assignment v_llen=N_VLOCLENGTH(v) sets v_llen to    *
 *     be the length of the local part of the vector v.        *
 *     The call N_VLOCLENGTH(v)=llen_v sets the local length   *
 *     of v to be llen_v.                                      *
 *                                                             *
 *     The assignment v_glen=N_VGLOBLENGTH(v) sets v_glen to   *
 *     be the global length of the vector v.                   *
 *     The call N_VGLOBLENGTH(v)=glen_v sets the global length *
 *     of v to be glen_v.                                      *
 *                                                             *
 * (4) N_VIth                                                  *
 *                                                             *
 *     In the following description, the components of the     *
 *     local part of an N_Vector are numbered 0..n-1, where n  *
 *     is the local length of (the local part of) v.           *
 *                                                             *
 *     The assignment r=N_VIth(v,i) sets r to be the value of  *
 *     the ith component of the local part of the vector v.    *
 *     The assignment N_VIth(v,i)=r sets the value of the      *
 *     ith local component of v to be r.                       *
 *                                                             *
 * Notes..                                                     *
 *                                                             *
 * Users who use the macros (1) and/or (2) must                *
 * #include<stdlib.h> since these macros expand to calls to    *
 * malloc and free.                                            *
 *                                                             *
 * When looping over the components of an N_Vector v, it is    *
 * more efficient to first obtain the component array via      *
 * v_data=N_VDATA(v) and then access v_data[i] within the      *
 * loop than it is to use N_VIth(v,i) within the loop.         *
 *                                                             *
 * N_VMAKE and N_VDISPOSE are similar to N_VNew and N_VFree.   *
 * The difference is one of responsibility for component       *
 * memory allocation and deallocation. N_VNew allocates memory *
 * for the N_Vector components and N_VFree frees the component *
 * memory allocated by N_VNew. For N_VMAKE and N_VDISPOSE, the *
 * component memory is allocated and freed by the user of      *
 * this module.                                                *
 *                                                             *
 * Similar remarks hold for N_VMAKE_S/N_VDISPOSE_S and         *
 * N_VNew_S/N_VFree_S.                                         *
 *                                                             *
 ***************************************************************/ 

#define N_VMAKE(v, v_data, machenv) \
        v = (N_Vector) malloc(sizeof(*v)); \
        v->data   = v_data; \
        v->length = machenv->local_vec_length; \
        v->global_length = machenv->global_vec_length; \
        v->machEnv = machenv

#define N_VDISPOSE(v) free(v)

#define N_VMAKE_S(vs, vs_data, s_len, machenv) \
        vs = (N_Vector *)malloc(s_len*sizeof(N_Vector *)); \
        for ((int)i=0; i<s_len; i++) { \
           vs[i] = (N_Vector) malloc(sizeof(N_Vector)); \
           vs[i]->data   = vs_data[i]; \
           vs[i]->length = machenv->local_vec_length; \
           vs[i]->global_length = machenv->global_vec_length; \
           vs[i]->machEnv = machenv; \
        }

#define N_VDISPOSE_S(vs, s_len) \
        for ((int)i=0; i<s_len; i++) N_VDISPOSE(vs[i]); \
        free(vs);

#define N_VDATA(v) (v->data)

#define N_VLOCLENGTH(v) (v->length)

#define N_VGLOBLENGTH(v) (v->global_length)

#define N_VIth(v,i) ((v->data)[i])


/* Part II: N_Vector Kernel Prototypes (Machine Environment-Independent) */

 
/***************************************************************
 *                                                             *
 * Memory Allocation and Deallocation: N_VNew, N_VFree         *
 *                                                             *
 ***************************************************************/


/***************************************************************
 *                                                             *
 * Function : N_VNew                                           *
 * Usage    : x = N_VNew(N, machEnv);                          *
 *-------------------------------------------------------------*
 *                                                             *
 * Returns a new N_Vector of length N. The parameter machEnv   *
 * is a pointer to machine environment-specific information.   *
 * It is ignored in the sequential machine environment and the *
 * user in this environment should simply pass NULL for this   *
 * argument. If there is not enough memory for a new N_Vector, *
 * then N_VNew returns NULL.                                   *
 *                                                             *
 ***************************************************************/

N_Vector N_VNew(integer n, machEnvType machEnv);


/***************************************************************
 *                                                             *
 * Function : N_VFree                                          *
 * Usage    : N_VFree(x);                                      *
 *-------------------------------------------------------------*
 *                                                             *
 * Frees the N_Vector x. It is illegal to use x after the call *
 * N_VFree(x).                                                 *
 *                                                             *
 ***************************************************************/

void N_VFree(N_Vector x);

/***************************************************************
 *                                                             *
 * Function : N_VNew_S                                         *
 * Usage    : x = N_VNew_S(Ns, N, machEnv);                    *
 *-------------------------------------------------------------*
 *                                                             *
 * Returns an array of Ns new N_Vectors of length N.           *
 * The parameter machEnv is a pointer to machine               *
 * environment-specific information. It is ignored in the      *
 * sequential machine environment and the user in this         *
 * environment should simply pass NULL for this argument.      * 
 * If there is not enough memory for a new array of N_Vectors  *
 * or for one of the components, then N_VNew_S returns NULL.   *
 *                                                             *
 ***************************************************************/

N_Vector *N_VNew_S(integer ns, integer n, machEnvType machEnv);

/***************************************************************
 *                                                             *
 * Function : N_VFree_S                                        *
 * Usage    : N_VFree_S(Ns, vs);                               *
 *-------------------------------------------------------------*
 *                                                             *
 * Frees the array of Ns N_Vectors vs.                         *
 * It is illegal to use vs after the call N_VFree_S(Ns,vs).    *
 *                                                             *
 ***************************************************************/

void N_VFree_S(integer ns, N_Vector *vs); 
 
/***************************************************************
 *                                                             *
 * N_Vector Arithmetic: N_VLinearSum, N_VConst, N_VProd,       *
 *                      N_VDiv, N_VScale, N_VAbs, N_VInv,      *
 *                      N_VAddConst                            *
 *                                                             *
 ***************************************************************/


/***************************************************************
 *                                                             *
 * Function  : N_VLinearSum                                    *
 * Operation : z = a x + b y                                   *
 *                                                             *
 ***************************************************************/

void N_VLinearSum(real a, N_Vector x, real b, N_Vector y, N_Vector z);


/***************************************************************
 *                                                             *
 * Function  : N_VConst                                        *
 * Operation : z[i] = c for i=0, 1, ..., N-1                   *
 *                                                             *
 ***************************************************************/

void N_VConst(real c, N_Vector z);


/***************************************************************
 *                                                             *
 * Function  : N_VProd                                         *
 * Operation : z[i] = x[i] * y[i] for i=0, 1, ..., N-1         *
 *                                                             *
 ***************************************************************/

void N_VProd(N_Vector x, N_Vector y, N_Vector z);


/***************************************************************
 *                                                             *
 * Function  : N_VDiv                                          *
 * Operation : z[i] = x[i] / y[i] for i=0, 1, ..., N-1         *
 *                                                             *
 ***************************************************************/

void N_VDiv(N_Vector x, N_Vector y, N_Vector z);


/***************************************************************
 *                                                             *
 * Function  : N_VScale                                        *
 * Operation : z = c x                                         *
 *                                                             *
 ***************************************************************/

void N_VScale(real c, N_Vector x, N_Vector z);


/***************************************************************
 *                                                             *
 * Function  : N_VAbs                                          *
 * Operation : z[i] = |x[i]|,   for i=0, 1, ..., N-1           *
 *                                                             *
 ***************************************************************/

void N_VAbs(N_Vector x, N_Vector z);


/***************************************************************
 *                                                             *
 * Function  : N_VInv                                          *
 * Operation : z[i] = 1.0 / x[i] for i = 0, 1, ..., N-1        *
 *-------------------------------------------------------------*
 *                                                             *
 * This routine does not check for division by 0. It should be *
 * called only with an N_Vector x which is guaranteed to have  *
 * all non-zero components.                                    *
 *                                                             *
 ***************************************************************/

void N_VInv(N_Vector x, N_Vector z);


/***************************************************************
 *                                                             *
 * Function  : N_VAddConst                                     *
 * Operation : z[i] = x[i] + b   for i = 0, 1, ..., N-1        *
 *                                                             *
 ***************************************************************/

void N_VAddConst(N_Vector x, real b, N_Vector z);
 
 
/***************************************************************
 *                                                             *
 * N_Vector Measures: N_VDotProd, N_VMaxNorm, VWrmsNorm,       *
 *                    N_VMin,  N_VWL2Norm, N_VL1Norm           *
 *                                                             *
 ***************************************************************/


/***************************************************************
 *                                                             *
 * Function : N_VDotProd                                       *
 * Usage    : dotprod = N_VDotProd(x, y);                      *
 *-------------------------------------------------------------*
 *                                                             *
 * Returns the value of the ordinary dot product of x and y:   *
 *                                                             *
 * -> sum (i=0 to N-1) {x[i] * y[i]}                           *
 *                                                             *
 * Returns 0.0 if N <= 0.                                      *
 *                                                             *
 ***************************************************************/

real N_VDotProd(N_Vector x, N_Vector y);


/***************************************************************
 *                                                             *
 * Function : N_VMaxNorm                                       *
 * Usage    : maxnorm = N_VMaxNorm(x);                         *
 *-------------------------------------------------------------*
 *                                                             *
 * Returns the maximum norm of x:                              *
 *                                                             *
 * -> max (i=0 to N-1) |x[i]|                                  *
 *                                                             *
 * Returns 0.0 if N <= 0.                                      *
 *                                                             *
 ***************************************************************/

real N_VMaxNorm(N_Vector x);


/***************************************************************
 *                                                             *
 * Function : N_VWrmsNorm                                      *
 * Usage    : wrmsnorm = N_VWrmsNorm(x, w);                    *
 *-------------------------------------------------------------*
 *                                                             *
 * Returns the weighted root mean square norm of x with        *
 * weight vector w:                                            *
 *                                                             *
 * -> sqrt [(sum (i=0 to N-1) {(x[i] * w[i])^2}) / N]          *
 *                                                             *
 * Returns 0.0 if N <= 0.                                      *
 *                                                             *
 ***************************************************************/

real N_VWrmsNorm(N_Vector x, N_Vector w);


/***************************************************************
 *                                                             *
 * Function : N_VMin                                           *
 * Usage    : min = N_VMin(x);                                 *
 *-------------------------------------------------------------*
 *                                                             *
 * Returns min x[i] if N > 0 and returns 0.0 if N <= 0.        *
 *          i                                                  *
 *                                                             *
 ***************************************************************/

real N_VMin(N_Vector x);

/***************************************************************
 *                                                             *
 * Function : N_VWL2Norm                                       *
 * Usage    : wl2norm = N_VWL2Norm(x, w);                      *
 *-------------------------------------------------------------*
 *                                                             *
 * Returns the weighted Euclidean L2 norm of x with            *
 * weight vector w:                                            *
 *                                                             *
 * -> sqrt [(sum (i=0 to N-1) {(x[i] * w[i])^2}) ]             *
 *                                                             *
 * Returns 0.0 if N <= 0.                                      *
 *                                                             *
 ***************************************************************/

real N_VWL2Norm(N_Vector x, N_Vector w);

 

/***************************************************************
 *                                                             *
 * Function : N_VL1Norm                                        *
 * Usage    : l1norm = N_VL1Norm(x);                           *
 *-------------------------------------------------------------*
 *                                                             *
 * Returns sum of ABS(x[i]) if N > 0 and returns 0.0 if N <= 0.*
 *          i                                                  *
 *                                                             *
 *     i.e., calculates and returns the L1 norm of x           *
 *                                                             *
 ***************************************************************/

real N_VL1Norm(N_Vector x);
 

/***************************************************************
 *                                                             *
 * Miscellaneous : N_VOneMask, N_VCompare, N_VInvTest,         *
 *        N_VConstrProdPos, N_VConstrMask, and N_VMinQuotient  *
 *                                                             *
 ***************************************************************/

/***************************************************************
 *                                                             *
 * Function  : N_VOneMask                                      *
 * Operation : x[i] = 1.0 if |x[i]| != 0.  i = 0, 1, ..., N-1  *
 *                    0.0 otherwise                            *
 *                                                             *
 ***************************************************************/

void N_VOneMask(N_Vector x);


/***************************************************************
 *                                                             *
 * Function  : N_VCompare                                      *
 * Operation : z[i] = 1.0 if |x[i]| >= c   i = 0, 1, ..., N-1  *
 *                    0.0 otherwise                            *
 *                                                             *
 ***************************************************************/


void N_VCompare(real c, N_Vector x, N_Vector z);


/***************************************************************
 *                                                             *
 * Function  : N_VInvTest                                      *
 * Operation : z[i] = 1.0 / x[i] with a test for x[i]==0.0     *
 *             before inverting x[i].                          *
 *-------------------------------------------------------------*
 *                                                             *
 * This routine returns TRUE if all components of x are        *
 * non-zero (successful inversion) and returns FALSE           *
 * otherwise.                                                  *
 *                                                             *
 ***************************************************************/

boole N_VInvTest(N_Vector x, N_Vector z);
 
 
/***************************************************************
 *                                                             *
 * Function : N_VConstrProdPos                                 *
 * Usage    : booltest = N_VConstrProdPos(c,x);                *
 *-------------------------------------------------------------*
 *                                                             *
 * Returns a boolean equal to                                  *
 *   FALSE if some c[i] != 0.0 and x[i]*c[i] <= 0.0,  or       *
 *   TRUE otherwise.                                           *
 *                                                             *
 * This routine is used for constraint checking.               *
 *                                                             *
 ***************************************************************/

boole N_VConstrProdPos(N_Vector c, N_Vector x);


/***************************************************************
 *                                                             *
 * Function  : N_VConstrMask                                   *
 * Operation : m[i] = 1.0 if constraint test fails for x[i]    *
 *             m[i] = 0.0 if constraint test passes for x[i]   *
 * where the constraint tests are as follows:                  *
 *    If c[i] = 2.0,  then x[i] must be > 0.0.                 *
 *    If c[i] = 1.0,  then x[i] must be >= 0.0.                *
 *    If c[i] = -1.0, then x[i] must be <= 0.0.                *
 *    If c[i] = -2.0, then x[i] must be < 0.0.                 *
 *-------------------------------------------------------------*
 *                                                             *
 * This routine returns a boole FALSE if any element failed    *
 * the constraint test, TRUE if all passed.  It also sets a    *
 * mask vector m, with elements equal to 1.0 where the         *
 * corresponding constraint test failed, and equal to 0.0      *
 * where the constraint test passed.                           *
 * This routine is specialized in that it is used only for     *
 * constraint checking.                                        *
 *                                                             *
 ***************************************************************/

boole N_VConstrMask(N_Vector c, N_Vector x, N_Vector m);   


/***************************************************************
 *                                                             *
 * Function  : N_VMinQuotient                                  *
 * Operation : minq  = min ( num[i]/denom[i]) over all i such  *
 *             that   denom[i] != 0.                           *
 *-------------------------------------------------------------*
 *                                                             *
 * This routine returns the minimum of the quotients obtained  *
 * by term-wise dividing num[i] by denom[i]. A zero element    *
 * in denom will be skipped. If no such quotients are found,   *
 * then the large value 1.e99 is returned.                     *
 *                                                             *
 ***************************************************************/

real N_VMinQuotient(N_Vector num, N_Vector denom);

 
/***************************************************************
 *                                                             *
 * Debugging Tools : N_VPrint                                  *
 *                                                             *
 ***************************************************************/

/***************************************************************
 *                                                             *
 * Function : N_VPrint                                         *
 * Usage    : N_VPrint(x);                                     *
 *-------------------------------------------------------------*
 *                                                             *
 * Prints the N_Vector x to stdout. Each component of x is     *
 * printed on a separate line using the %g specification. This *
 * routine is provided as an aid in debugging code which uses  *
 * this vector package.                                        *
 *                                                             *
 ***************************************************************/

void N_VPrint(N_Vector x);
 

#endif
#ifdef __cplusplus
}
#endif
