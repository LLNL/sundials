.. ----------------------------------------------------------------
   SUNDIALS Copyright Start
   Copyright (c) 2002-2022, Lawrence Livermore National Security
   and Southern Methodist University.
   All rights reserved.

   See the top-level LICENSE and NOTICE files for details.

   SPDX-License-Identifier: BSD-3-Clause
   SUNDIALS Copyright End
   ----------------------------------------------------------------

.. _SUNDIALS.Logging:

SUNDIALS Status Logging
=======================

.. versionadded:: 6.2.0

SUNDIALS includes a built-in logging functionality which can be used to direct
error messages, warning messages, informational output, and debugging output to
specified files. This capability requires enabling both build-time and run-time
options to ensure the best possible performance is achieved.

.. _SUNDIALS.Logging.Enabling:

Enabling Logging
----------------

To enable logging, the CMake option :cmakeop:`SUNDIALS_LOGGING_LEVEL` must be
set to a value greater than ``0`` when configuring SUNDIALS. This option
specifies the maximum desired output level. See the documentation entry for
:cmakeop:`SUNDIALS_LOGGING_LEVEL` for the numeric values correspond to errors,
warnings, info output, and debug output where errors < warnings < info
output < debug output < extra debug output. If it is desired that the logger is
MPI-aware, then the option :cmakeop:`SUNDIALS_LOGGING_ENABLE_MPI` is set to
``TRUE``. More details in regards to configuring SUNDIALS with CMake can be
found in :numref:`Installation`.

When SUNDIALS is built with logging enabled, then the default logger (stored in
the :c:type:`SUNContext` object) may be configured through environment variables
without any changes to user code. The available environment variables are:

.. code-block::

   SUNLOGGER_ERROR_FILENAME
   SUNLOGGER_WARNING_FILENAME
   SUNLOGGER_INFO_FILENAME
   SUNLOGGER_DEBUG_FILENAME

These environment variables may be set to a filename string. There are two
special filenames: ``stdout`` and ``stderr``. These two filenames will
result in output going to the standard output file and standard error file.
The different variables may all be set to the same file, or to distinct files,
or some combination there of. To disable output for one of the streams, then
do not set the environment variable, or set it to an empty string.

.. warning::

   A non-default logger should be created prior to any other SUNDIALS calls
   in order to capture all log events.

.. note::

   If :cmakeop:`SUNDIALS_LOGGING_LEVEL` was set to ``1`` (corresponding to
   error-level output) at build-time, then setting the environment variable
   ``SUNLOGGER_INFO_FILENAME`` will do nothing.

.. note::

   Extra debugging output is turned on by setting :cmakeop:`SUNDIALS_LOGGING_LEVEL` to 5.
   This extra output includes vector-values (so long as the :c:type:`N_Vector` used
   supports printing).


Logger API
----------

The central piece of the Logger API is the :c:type:`SUNLogger` type:

.. c:type:: struct SUNLogger_* SUNLogger

When SUNDIALS is built with logging enabled, a default logging object is stored
in the :c:type:`SUNContext` object and can be accessed with a call to
:c:func:`SUNContext_GetLogger`.

The enumerated type :c:enum:`SUNLogLevel` is used by some of the logging
functions to identify the output level or file.

.. c:enum:: SUNLogLevel

   The SUNDIALS logging level

.. c:enumerator:: SUN_LOGLEVEL_ALL

   Represents all output levels

.. c:enumerator:: SUN_LOGLEVEL_NONE

   Represents none of the output levels

.. c:enumerator:: SUN_LOGLEVEL_ERROR

   Represents error-level logging messages

.. c:enumerator:: SUN_LOGLEVEL_WARNING

   Represents warning-level logging messages

.. c:enumerator:: SUN_LOGLEVEL_INFO

   Represents info-level logging messages

.. c:enumerator:: SUN_LOGLEVEL_DEBUG

   Represents deubg-level logging messages


The :c:type:`SUNLogger` class provides the following methods.


.. c:function:: int SUNLogger_Create(void* comm, int output_rank, SUNLogger* logger)

   Creates a new :c:type:`SUNLogger` object.

   **Arguments:**
      * ``comm`` -- a pointer to the MPI communicator if MPI is enabled,
        otherwise can be ``NULL``.
      * ``output_rank`` -- the MPI rank used for output (can be ``-1`` to print
        to all ranks).
      * ``logger`` -- [in,out] On input this is a pointer to a
         :c:type:`SUNLogger`, on output it will point to a new
         :c:type:`SUNLogger` instance.

   **Returns:**
      * Returns zero if successful, or non-zero if an error occurred.


.. c:function:: int SUNLogger_CreateFromEnv(void* comm, SUNLogger* logger)

   Creates a new :c:type:`SUNLogger` object and opens the output streams/files
   from the environment variables:

   .. code-block::

      SUNLOGGER_ERROR_FILENAME
      SUNLOGGER_WARNING_FILENAME
      SUNLOGGER_INFO_FILENAME
      SUNLOGGER_DEBUG_FILENAME

   **Arguments:**
      * ``comm`` -- a pointer to the MPI communicator if MPI is enabled,
        otherwise can be ``NULL``.
      * ``logger`` -- [in,out] On input this is a pointer to a
         :c:type:`SUNLogger`, on output it will point to a new
         :c:type:`SUNLogger` instance.

   **Returns:**
      * Returns zero if successful, or non-zero if an error occurred.


.. c:function:: int SUNLogger_SetErrorFilename(SUNLogger logger, const char* error_filename)

   Sets the filename for error output.

   **Arguments:**
      * ``logger`` -- a :c:type:`SUNLogger` object.
      * ``error_filename`` -- the name of the file to use for error output.

   **Returns:**
      * Returns zero if successful, or non-zero if an error occurred.


.. c:function:: int SUNLogger_SetWarningFilename(SUNLogger logger, const char* warning_filename)

   Sets the filename for warning output.

   **Arguments:**
      * ``logger`` -- a :c:type:`SUNLogger` object.
      * ``warning_filename`` -- the name of the file to use for warning output.

   **Returns:**
      * Returns zero if successful, or non-zero if an error occurred.


.. c:function:: int SUNLogger_SetInfoFilename(SUNLogger logger, const char* info_filename)

   Sets the filename for info output.

   **Arguments:**
      * ``logger`` -- a :c:type:`SUNLogger` object.
      * ``info_filename`` -- the name of the file to use for info output.

   **Returns:**
      * Returns zero if successful, or non-zero if an error occurred.


.. c:function:: int SUNLogger_SetDebugFilename(SUNLogger logger, const char* debug_filename)

   Sets the filename for debug output.

   **Arguments:**
      * ``logger`` -- a :c:type:`SUNLogger` object.
      * ``debug_filename`` -- the name of the file to use for debug output.

   **Returns:**
      * Returns zero if successful, or non-zero if an error occurred.


.. c:function:: int SUNLogger_QueueMsg(SUNLogger logger, SUNLogLevel lvl, const char* scope, const char* label, const char* msg_txt, ...)

   Queues a message to the output log level.

   **Arguments:**
      * ``logger`` -- a :c:type:`SUNLogger` object.
      * ``lvl`` -- the message log level (i.e. error, warning, info, debug).
      * ``scope`` -- the message scope (e.g. the function name).
      * ``label`` -- the message label.
      * ``msg_txt`` -- the message text itself.
      * ``...`` -- the format string arguments

   **Returns:**
      * Returns zero if successful, or non-zero if an error occurred.

   .. warning::

      When compiling for ANSI C / C89 / C90 (and without compiler extensions),
      it is dangerous to pass any user input to this function because it falls
      back to using ``sprintf`` with a fixed buffer size.

      It is **highly recommended** to compile with C99 or newer if your compiler
      does not support ``snprintf`` through extensions.


.. c:function:: int SUNLogger_Flush(SUNLogger logger, SUNLogLevel lvl)

   Flush the message queue(s).

   **Arguments:**
      * ``logger`` -- a :c:type:`SUNLogger` object.
      * ``lvl`` -- the message log level (i.e. error, warning, info, debug or
        all).

   **Returns:**
      * Returns zero if successful, or non-zero if an error occurred.


.. c:function:: int SUNLogger_GetOutputRank(SUNLogger logger, int* output_rank)

   Get the output MPI rank for the logger.

   **Arguments:**
      * ``logger`` -- a :c:type:`SUNLogger` object.
      * ``output_rank`` -- [in,out] On input this is a pointer to an int, on
        output it points to the int holding the output rank.

   **Returns:**
      * Returns zero if successful, or non-zero if an error occurred.


.. c:function:: int SUNLogger_Destroy(SUNLogger* logger)

   Free the memory for the :c:type:`SUNLogger` object.

   **Arguments:**
      * ``logger`` -- a pointer to the :c:type:`SUNLogger` object.

   **Returns:**
      * Returns zero if successful, or non-zero if an error occur.


.. _SUNDIALS.Logging.Example:

Example Usage
-------------

As previously mentioned, if it is enabled at build time, there is a default
:c:type:`SUNLogger` attached to a :c:type:`SUNContext` instance when it is
created. This logger can be configured using the environment variables, e.g.,

.. code-block::

   SUNDIALS_INFO_FILENAME=stdout ./examples/cvode/serial/cvKrylovDemo_ls

SUNDIALS also includes several example codes that demonstrate how to use the
logging interface via the C API.

.. code-block::

   examples/arkode/CXX_serial/ark_analytic_sys.cpp
   examples/cvode/serial/cvAdvDiff_bnd.c
   examples/cvode/parallel/cvAdvDiff_diag_p.c
   examples/kinsol/CXX_parallel/kin_em_p.cpp
   examples/kinsol/CUDA_mpi/kin_em_mpicuda.cpp
