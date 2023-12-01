.. ----------------------------------------------------------------
   SUNDIALS Copyright Start
   Copyright (c) 2002-2023, Lawrence Livermore National Security
   and Southern Methodist University.
   All rights reserved.

   See the top-level LICENSE and NOTICE files for details.

   SPDX-License-Identifier: BSD-3-Clause
   SUNDIALS Copyright End
   ----------------------------------------------------------------

.. _SUNDIALS.Errors:

Error Checking
==============

.. versionadded:: 7.0.0

Errors that occur within SUNDIALS are dealt with through error codes and error handler callback
functions. The error code type :c:type:`SUNErrCode` is just a typedef to an ``int``:

.. c:type:: typedef int SUNErrCode 

Errors are always negative, while success is defined with the ``SUN_SUCCESS`` code and has the value ``0``.
To see all of the possible error codes, refer to the ``sundials/sundials_errors.h`` header file 
(`here <https://github.com/LLNL/sundials/blob/main/include/sundials/sundials_errors.h>_`).

Functions in the SUNDIALS core API (i.e., ``SUN`` or ``N_V`` functions only) either return a
:c:type:`SUNErrCode`, or (if they don't return a :c:type:`SUNErrCode`) they internally record an
error code (if an error occurs) within the :c:type:`SUNContext` for the execution stream. 
This "last error" is accessible via the :c:func:`SUNContext_GetLastError` or
 :c:func:`SUNContext_PeekLastError` functions.

Thus, in user code, SUNDIALS core functions can be checked for errors in one of two ways:

.. code-block:: C

  SUNContext sunctx;
  SUNErrCode sunerr; 
  N_Vector v;
  int length;
  sunrealtype dotprod;
  
  // Every code that uses SUNDIALS must create a SUNContext.
  sunctx = SUNContext_Create(...);

  // Create a SUNDIALS serial vector.
  // Some functions do not return an error code.
  // We have to check for errors in these functions using SUNContext_GetLastError.
  length = 2;
  v = N_VNew_Serial(length, sunctx);
  sunerr = SUNContext_GetLastError(sunctx);
  if (sunerr) {  /* an error occured, do something */ }

  // If the function returns a SUNErrCode, we can check it directly
  sunerr = N_VLinearCombination(...);
  if (sunerr) {  /* an error occured, do something */ }

  // Another function that does not return a SUNErrCode.
  dotprod = N_VDotProd(...);
  SUNContext_GetLastError(sunctx);
  if (sunerr) { 
   /* an error occured, do something */ 
  } else {
    print("dotprod = %.2f\n", dotprod);
  }


.. note::

  It is recommended in most cases that users check for an error after calling SUNDIALS functions.
  However, users concerned with getting the most performance might choose to exclude or limit these checks.


.. warning::

  If a function returns a :c:type:`SUNErrCode` then the return value is the only place the error is available. 
  I.e., these functions do not store their error code as the "last error" so it is invalid to use
  :c:func:`SUNContext_GetLastError` to check these functions for errors.


.. _SUNDIALS.Errors.Handlers:

Error Handler Functions
-----------------------

Errors that occur internally to SUNDIALS result in an error handler function being called. 
These error handler functions have the type

.. c:type:: typedef int (*SUNErrHandlerFn)(int line, const char* func, const char* file,
                                           const char* msg, SUNErrCode err_code,
                                           void* err_user_data, SUNContext sunctx);

SUNDIALS provides a few different error handlers that can be used, or a custom one defined by the
user can be provided (useful for linking SUNDIALS errors to your application's error handling).
The default error handler is :c:func:`SUNLogErrHandlerFn` which logs an error to a specified
file or ``stderr`` if no file is specified.

The error handlers provided in SUNDIALS are:

.. c:function:: void SUNLogErrHandlerFn(int line, const char* func, const char* file, \
                                        const char* msg, SUNErrCode err_code, \
                                        void* err_user_data, SUNContext sunctx)

  Logs the error that occurred using the :c:type:`SUNLogger` for ``sunctx``.
  This is the default error handler.

  :param line: the line number at which the error occured
  :param func: the function in which the error occured
  :param file: the file in which the error occured 
  :param msg: the message to log, if this is ``NULL`` then the default error message for the error code will be used
  :param err_code: the error code for the error that occured
  :param err_user_data: the user pointer provided to :c:func:`SUNContext_PushErrHandler`
  :param sunctx: pointer to a valid :c:type:`SUNContext` object

  :return: ``void``

.. c:function:: void SUNAbortErrHandlerFn(int line, const char* func, const char* file, \
                                          const char* msg, SUNErrCode err_code, \
                                          void* err_user_data, SUNContext sunctx)

  Logs the error and aborts the program if an error occured.

  :param line: the line number at which the error occured
  :param func: the function in which the error occured
  :param file: the file in which the error occured 
  :param msg: the message to log, if this is ``NULL`` then the default error message for the error code will be used
  :param err_code: the error code for the error that occured
  :param err_user_data: the user pointer provided to :c:func:`SUNContext_PushErrHandler`
  :param sunctx: pointer to a valid :c:type:`SUNContext` object

  :return: ``void``


.. c:function:: void SUNAssertErrHandlerFn(int line, const char* func, const char* file, \
                                           const char* stmt, SUNErrCode err_code, \
                                           void* err_user_data, SUNContext sunctx)

  Logs the error and aborts the program if an error occured, but with a messgae that displays
  the assertion that failed.

  :param line: the line number at which the error occured
  :param func: the function in which the error occured
  :param file: the file in which the error occured 
  :param stmt: the statement that failed the assertion
  :param err_code: the error code for the error that occured
  :param err_user_data: the user pointer provided to :c:func:`SUNContext_PushErrHandler`
  :param sunctx: pointer to a valid :c:type:`SUNContext` object

  :return: ``void``


.. c:function:: void SUNMPIAbortErrHandlerFn(int line, const char* func, const char* file, \
                                             const char* msg, SUNErrCode err_code, \
                                             void* err_user_data, SUNContext sunctx)

  Logs the error and calls ``MPI_Abort`` if an error occured.

  :param line: the line number at which the error occured
  :param func: the function in which the error occured
  :param file: the file in which the error occured 
  :param msg: the message to log, if this is ``NULL`` then the default error message for the error code will be used
  :param err_code: the error code for the error that occured
  :param err_user_data: the user pointer provided to :c:func:`SUNContext_PushErrHandler`
  :param sunctx: pointer to a valid :c:type:`SUNContext` object

  :return: ``void``


.. c:function:: void SUNMPIAssertErrHandlerFn(int line, const char* func, const char* file, \
                                              const char* stmt, SUNErrCode err_code, \
                                              void* err_user_data, SUNContext sunctx)

  Logs the error and calls ``MPI_Abort`` if an error occured, but with a messgae that displays
  the assertion that failed.

  :param line: the line number at which the error occured
  :param func: the function in which the error occured
  :param file: the file in which the error occured 
  :param stmt: the statement that failed the assertion
  :param err_code: the error code for the error that occured
  :param err_user_data: the user pointer provided to :c:func:`SUNContext_PushErrHandler`
  :param sunctx: pointer to a valid :c:type:`SUNContext` object

  :return: ``void``