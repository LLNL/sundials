..
   Programmer(s): Daniel R. Reynolds @ SMU
   ----------------------------------------------------------------
   SUNDIALS Copyright Start
   Copyright (c) 2002-2024, Lawrence Livermore National Security
   and Southern Methodist University.
   All rights reserved.

   See the top-level LICENSE and NOTICE files for details.

   SPDX-License-Identifier: BSD-3-Clause
   SUNDIALS Copyright End
   ----------------------------------------------------------------

.. _SUNMatrix.Ops:

Description of the SUNMATRIX operations
=======================================

For each of the ``SUNMatrix`` operations, we give the name, usage
of the function, and a description of its mathematical operations
below.


.. c:function:: SUNMatrix_ID SUNMatGetID(SUNMatrix A)

   Returns the type identifier for the matrix *A*.  It is used to determine the
   matrix implementation type (e.g. dense, banded, sparse,...) from the abstract
   ``SUNMatrix`` interface.  This is used to assess compatibility with
   SUNDIALS-provided linear solver implementations.  Returned values
   are given in :numref:`SUNMatrix.Description.matrixIDs`

   Usage:

   .. code-block:: c

      id = SUNMatGetID(A);


.. c:function:: SUNMatrix SUNMatClone(SUNMatrix A)

   Creates a new ``SUNMatrix`` of the same type as an existing
   matrix *A* and sets the *ops* field.  It does not copy the matrix values,
   but rather allocates storage for the new matrix.

   Usage:

   .. code-block:: c

      B = SUNMatClone(A);


.. c:function:: void SUNMatDestroy(SUNMatrix A)

   Destroys the ``SUNMatrix`` *A* and frees memory allocated for its
   internal data.

   Usage:

   .. code-block:: c

      SUNMatDestroy(A);


.. c:function:: SUNErrCode SUNMatSpace(SUNMatrix A, long int *lrw, long int *liw)

   Returns the storage requirements for the matrix *A*.  *lrw*
   contains the number of sunrealtype words and *liw* contains the number
   of integer words.  The return value denotes success/failure of the
   operation.

   This function is advisory only, for use in determining a user's total
   space requirements; it could be a dummy function in a user-supplied
   ``SUNMatrix`` module if that information is not of interest.

   Usage:

   .. code-block:: c

      retval = SUNMatSpace(A, &lrw, &liw);


.. c:function:: SUNErrCode SUNMatZero(SUNMatrix A)

   Zeros all entries of the ``SUNMatrix`` *A*.  The return value
   denotes the success/failure of the operation:

   .. math::
      A_{i,j} = 0, \quad i=1,\ldots,m, \; j=1,\ldots,n.

   Usage:

   .. code-block:: c

      retval = SUNMatZero(A);


.. c:function:: SUNErrCode SUNMatCopy(SUNMatrix A, SUNMatrix B)

   Performs the operation *B \gets A* for all entries of the matrices *A*
   and *B*.  The return value denotes the success/failure of
   the operation:

   .. math::
      B_{i,j} = A_{i,j}, \quad i=1,\ldots,m, \; j=1,\ldots,n.

   Usage:

   .. code-block:: c

      retval = SUNMatCopy(A,B);


.. c:function:: SUNErrCode SUNMatScaleAdd(sunrealtype c, SUNMatrix A, SUNMatrix B)

   Performs the operation *A \gets cA + B*.  The return value
   denotes the success/failure of the operation:

   .. math::
      A_{i,j} = cA_{i,j} + B_{i,j}, \quad i=1,\ldots,m, \; j=1,\ldots,n.

   Usage:

   .. code-block:: c

      retval = SUNMatScaleAdd(c, A, B);


.. c:function:: SUNErrCode SUNMatScaleAddI(sunrealtype c, SUNMatrix A)

   Performs the operation *A \gets cA + I*.  The return value
   denotes the success/failure of the operation:

   .. math::
      A_{i,j} = cA_{i,j} + \delta_{i,j}, \quad i,j=1,\ldots,n.

   Usage:

   .. code-block:: c

      retval = SUNMatScaleAddI(c, A);


.. c:function:: SUNErrCode SUNMatMatvecSetup(SUNMatrix A)

   Performs any setup necessary to perform a matrix-vector product.
   The return value denotes the success/failure of the
   operation. It is useful for SUNMatrix implementations which need to
   prepare the matrix itself, or communication structures before performing
   the matrix-vector product.

   Usage:

   .. code-block:: c

      retval = SUNMatMatvecSetup(A);

.. c:function:: SUNErrCode SUNMatMatvec(SUNMatrix A, N_Vector x, N_Vector y)

   Performs the matrix-vector product *y \gets Ax*.  It should
   only be called with vectors *x* and *y* that are compatible with
   the matrix *A* -- both in storage type and dimensions.  The return
   value denotes the success/failure of the operation:

   .. math::
      y_i = \sum_{j=1}^n A_{i,j} x_j, \quad i=1,\ldots,m.

   Usage:

   .. code-block:: c

      retval = SUNMatMatvec(A, x, y);
