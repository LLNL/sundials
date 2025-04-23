..
   Author(s): David J. Gardner, Cody J. Balos @ LLNL
   -----------------------------------------------------------------------------
   SUNDIALS Copyright Start
   Copyright (c) 2002-2025, Lawrence Livermore National Security
   and Southern Methodist University.
   All rights reserved.

   See the top-level LICENSE and NOTICE files for details.

   SPDX-License-Identifier: BSD-3-Clause
   SUNDIALS Copyright End
   -----------------------------------------------------------------------------

.. _SourceCode.Style:

Style
=====

In this section we describe the style conventions and guidelines for SUNDIALS
source code.

Formatting
----------

All new code added to SUNDIALS should be formatted with `clang-format
<https://clang.llvm.org/docs/ClangFormat.html>`_ for C/C++, `fprettify
<https://github.com/fortran-lang/fprettify>`_ for Fortran, `cmake-format
<https://cmake-format.readthedocs.io>`_ for CMake, and `black
<https://black.readthedocs.io>`_ for Python. The ``.clang-format`` file in the
root of the project defines our configuration for clang-format. We use the
default fprettify settings, except we use 2-space indentation. The
``.cmake-format.py`` file in the root of the project defines our configuration
for cmake-format. We also use the default black settings.


To apply ``clang-format``, ``fprettify``, ``cmake-format``, and ``black`` you
can run:

.. code-block:: shell

   ./scripts/format.sh <path to directories or files to format>

.. warning::

   The output of ``clang-format`` is sensitive to the ``clang-format`` version. We recommend
   that you use version ``17.0.4``, which can be installed from source or with Spack. Alternatively,
   when you open a pull request on GitHub, an action will run ``clang-format`` on the code. If any
   formatting is required, the action will fail. Commenting with the magic keyword ``/autofix`` will
   kick off a GitHub action which will automatically apply the formatting changes needed.

If clang-format breaks lines in a way that is unreadable, use ``//`` to break the line. For example,
sometimes (mostly in C++ code) you may have code like this:

.. code-block:: cpp

   MyClass::callAFunctionOfSorts::doSomething().doAnotherThing().doSomethingElse();

That you would like to format as (for readability):

.. code-block:: cpp

   MyObject::callAFunctionOfSorts()
         .doSomething()
         .doAnotherThing()
         .doSomethingElse();

Clang-format might produce something like:

.. code-block:: cpp

   MyObject::callAFunctionOfSorts().doSomething().doAnotherThing()
         .doSomethingElse();


unless you add the ``//``

.. code-block:: cpp

   MyObject::callAFunctionOfSorts()
         .doSomething()       //
         .doAnotherThing()    //
         .doSomethingElse();  //

There are other scenarios (e.g., a function call with a lot of parameters) where
doing this type of line break is useful for readability too.

.. It may be necessary to override clang-tidy at times. This can be done with the
.. ``NOLINT`` magic comments e.g.,

.. .. code-block:: cpp

..   template<class GkoSolverType, class GkoMatrixType>
..   int SUNLinSolFree_Ginkgo(SUNLinearSolver S)
..   {
..     auto solver{static_cast<LinearSolver<GkoSolverType, GkoMatrixType>*>(S->content)};
..     delete solver; // NOLINT
..     return SUNLS_SUCCESS;
..   }

..   class BaseObject {
..   protected:
..     // NOLINTNEXTLINE(cppcoreguidelines-non-private-member-variables-in-classes)
..     SUNContext sunctx_{};
..   };

.. See the clang-tidy documentation for more details.


.. _Style.Output:

Output
------

For consistent formatting of :c:type:`sunrealtype`, the following macros are
available.

.. c:macro:: SUN_FORMAT_E

   A format specifier for scientific notation. This should be used when
   displaying arrays, matrices, and tables where fixed width alignment aids with
   readability.

   **Example usage:**

   .. code-block:: C

      for (i = 0; i < N; i++) {
         fprintf(outfile, SUN_FORMAT_E "\n", xd[i]);
      }

.. c:macro:: SUN_FORMAT_G

   A format specifier for scientific or standard notation, whichever is more
   compact. It is more reader-friendly than :c:macro:`SUN_FORMAT_E` and should
   be used in all cases not covered by that macro.

   **Example usage:**

   .. code-block:: C

      SUNLogInfo(sunctx->logger, "label", "x = " SUN_FORMAT_G, x);

.. c:macro:: SUN_FORMAT_SG

   Like :c:macro:`SUN_FORMAT_G` but with a leading plus or minus sign.


To aid in printing statistics in functions like :c:func:`CVodePrintAllStats`,
the following utility functions are available.

.. c:function:: void sunfprintf_real(FILE* fp, SUNOutputFormat fmt, sunbooleantype start, const char* name, sunrealtype value)

   Writes a :c:type:`sunrealtype` value to a file pointer using the specified
   format.

   :param fp: The output file pointer.
   :param fmt: The output format.
   :param start: :c:macro:`SUNTRUE` if the value is the first in a series of
                 statistics, and :c:macro:`SUNFALSE` otherwise.
   :param name: The name of the statistic.
   :param value: The value of the statistic.

.. c:function:: void sunfprintf_long(FILE* fp, SUNOutputFormat fmt, sunbooleantype start, const char* name, long value)

   Writes a long value to a file pointer using the specified format.

   :param fp: The output file pointer.
   :param fmt: The output format.
   :param start: :c:macro:`SUNTRUE` if the value is the first in a series of
                 statistics, and :c:macro:`SUNFALSE` otherwise.
   :param name: The name of the statistic.
   :param value: The value of the statistic.

.. c:function:: void sunfprintf_long_array(FILE* fp, SUNOutputFormat fmt, sunbooleantype start, const char* name, long* value, size_t count)

   Writes an array of long values to a file pointer using the specified format.

   :param fp: The output file pointer.
   :param fmt: The output format.
   :param start: :c:macro:`SUNTRUE` if the value is the first in a series of
                 statistics, and :c:macro:`SUNFALSE` otherwise.
   :param name: The name of the statistic.
   :param value: Pointer to the array.
   :param count: The number of elements in the array.

.. _Style.Logging:

Logging
-------

Use the macros below to add informational and debugging messages to SUNDIALS
code rather than adding ``#ifdef SUNDIALS_LOGGING_<level>`` / ``#endif`` blocks
containing calls to :c:func:`SUNLogger_QueueMsg`. Error and warning messages are
handled through package-specific ``ProcessError`` functions or the ``SUNAssert``
and ``SUNCheck`` macros.

The logging macros help ensure messages follow the required format presented in
:numref:`SUNDIALS.Logging.Enabling` and used by the ``suntools`` Python module
for parsing logging output. For informational and debugging output the log
message payload (the part after the brackets) must be either be a
comma-separated list of key-value pairs with the key and value separated by an
equals sign with a space on either side e.g.,

.. code-block:: C

   /* log an informational message */
   SUNLogInfo(sunctx->logger, "begin-step", "t = " SUN_FORMAT_G ", h = " SUN_FORMAT_G, t, h);

   /* log a debugging message */
   SUNLogDebug(sunctx->logger, "error-estimates",
               "eqm1 = " SUN_FORMAT_G ", eq = " SUN_FORMAT_G ", eqp1 = " SUN_FORMAT_G,
               eqm1, eq, eqp1);

or the name of a vector/array followed by ``(:) =`` with each vector/array entry
written to a separate line e.g., a vector may be logged with

.. code-block:: C

   SUNLogExtraDebugVec(sunctx->logger, "new-solution", ynew, "ynew(:) =");

where the message can contain format specifiers e.g., if ``Fe`` is an array of
vectors you may use

.. code-block:: C

   SUNLogExtraDebugVec(sunctx->logger, "new-solution", Fe[i], "Fe_%d(:) =", i);

To assist in parsing logging messages, ``begin-`` and ``end-`` markers are used
in the log message ``label`` field to denote where particular regions begin and
end. When adding a new ``begin-`` / ``end-`` label the ``logs.py`` script will
need to be updated accordingly. The region markers currently supported by the
Python module for parsing log files are as follows:

* ``begin-step-attempt`` / ``end-step-attempt``

* ``begin-nonlinear-solve`` / ``end-nonlinear-solve``

* ``begin-nonlinear-iterate`` / ``end-nonlinear-iterate``

* ``begin-linear-solve`` / ``end-linear-solve``

* ``begin-linear-iterate`` / ``end-linear-iterate``

* ``begin-group`` / ``end-group``

* ``begin-stage`` / ``end-stage``

* ``begin-fast-steps`` / ``end-fast-steps``

* ``begin-mass-linear-solve`` / ``end-mass-linear-solve``

* ``begin-compute-solution`` / ``end-compute-solution``

* ``begin-compute-embedding`` / ``end-compute-embedding``

Logging Macros
^^^^^^^^^^^^^^

.. versionadded:: 7.2.0

To log informational messages use the following macros:

.. c:macro:: SUNLogInfo(logger, label, msg_txt, ...)

   When information logging is enabled this macro expands to a call to
   :c:func:`SUNLogger_QueueMsg` to log an informational message. Otherwise, this
   expands to nothing.

   :param logger: the :c:type:`SUNLogger` to handle the message.
   :param label: the ``const char*`` message label.
   :param msg_txt: the ``const char*`` message text, may contain format
                   specifiers.
   :param ...: the arguments for format specifiers in ``msg_txt``.

.. c:macro:: SUNLogInfoIf(condition, logger, label, msg_txt, ...)

   When information logging is enabled this macro expands to a conditional call
   to :c:func:`SUNLogger_QueueMsg` to log an informational message. Otherwise,
   this expands to nothing.

   :param condition: a boolean expression that determines if the log message
                     should be queued.
   :param logger: the :c:type:`SUNLogger` to handle the message.
   :param label: the ``const char*`` message label.
   :param msg_txt: the ``const char*`` message text, may contain format.
                   specifiers.
   :param ...: the arguments for format specifiers in ``msg_txt``.

To log debugging messages use the following macros:

.. c:macro:: SUNLogDebug(logger, label, msg_txt, ...)

   When debugging logging is enabled this macro expands to a call to
   :c:func:`SUNLogger_QueueMsg` to log a debug message. Otherwise, this expands
   to nothing.

   :param logger: the :c:type:`SUNLogger` to handle the message.
   :param label: the ``const char*`` message label.
   :param msg_txt: the ``const char*`` message text, may contain format.
                   specifiers.
   :param ...: the arguments for format specifiers in ``msg_txt``.

.. c:macro:: SUNLogDebugIf(condition, logger, label, msg_txt, ...)

   When debugging logging is enabled this macro expands to a conditional call to
   :c:func:`SUNLogger_QueueMsg` to log a debug message. Otherwise, this expands
   to nothing.

   :param condition: a boolean expression that determines if the log message
                     should be queued.
   :param logger: the :c:type:`SUNLogger` to handle the message.
   :param label: the ``const char*`` message label.
   :param msg_txt: the ``const char*`` message text, may contain format.
                   specifiers.
   :param ...: the arguments for format specifiers in ``msg_txt``.

To log extra debugging messages use the following macros:

.. c:macro:: SUNLogExtraDebug(logger, label, msg_txt, ...)

   When extra debugging logging is enabled, this macro expands to a call to
   :c:func:`SUNLogger_QueueMsg` to log an extra debug message. Otherwise, this expands
   to nothing.

   :param logger: the :c:type:`SUNLogger` to handle the message.
   :param label: the ``const char*`` message label.
   :param msg_txt: the ``const char*`` message text, may contain format
                   specifiers.
   :param ...: the arguments for format specifiers in ``msg_txt``.

.. c:macro:: SUNLogExtraDebugIf(condition, logger, label, msg_txt, ...)

   When extra debugging logging is enabled, this macro expands to a conditional
   call to :c:func:`SUNLogger_QueueMsg` to log an extra debug message. Otherwise, this
   expands to nothing.

   :param condition: a boolean expression that determines if the log message
                     should be queued.
   :param logger: the :c:type:`SUNLogger` to handle the message.
   :param label: the ``const char*`` message label.
   :param msg_txt: the ``const char*`` message text, may contain format
                   specifiers.
   :param ...: the arguments for format specifiers in ``msg_txt``.

.. c:macro:: SUNLogExtraDebugVec(logger, label, vec, msg_txt, ...)

   When extra debugging logging is enabled, this macro expands to a call to
   :c:func:`SUNLogger_QueueMsg` and :c:func:`N_VPrintFile` to log an extra
   debug message and output the vector data. Otherwise, this expands to nothing.

   :param logger: the :c:type:`SUNLogger` to handle the message.
   :param label: the ``const char*`` message label.
   :param vec: the ``N_Vector`` to print.
   :param msg_txt: the ``const char*`` message text, may contain format
                   specifiers.
   :param ...: the arguments for format specifiers in ``msg_txt``.

.. c:macro:: SUNLogExtraDebugVecIf(condition, logger, label, vec, msg_txt, ...)

   When extra debugging logging is enabled, this macro expands to a conditional
   call to :c:func:`SUNLogger_QueueMsg` and :c:func:`N_VPrintFile` to log an extra
   debug message and output the vector data. Otherwise, this expands to nothing.

   :param condition: a boolean expression that determines if the log message
                     should be queued.
   :param logger: the :c:type:`SUNLogger` to handle the message.
   :param label: the ``const char*`` message label.
   :param vec: the ``N_Vector`` to print.
   :param msg_txt: the ``const char*`` message text, may contain format
                   specifiers.
   :param ...: the arguments for format specifiers in ``msg_txt``.

.. c:macro:: SUNLogExtraDebugVecArray(logger, label, nvecs, vecs, msg_txt)

   When extra debugging logging is enabled, this macro expands to a loop calling
   :c:func:`SUNLogger_QueueMsg` and :c:func:`N_VPrintFile` for each vector in
   the vector array to log an extra debug message and output the vector data.
   Otherwise, this expands to nothing.

   :param logger: the :c:type:`SUNLogger` to handle the message.
   :param label: the ``const char*`` message label.
   :param nvecs: the ``int`` number of vectors to print.
   :param vecs: the ``N_Vector*`` (vector array) to print.
   :param msg_txt: the ``const char*`` message text, must contain a format
                   specifier for the vector array index.

   .. warning::

      The input parameter ``msg_txt`` **must** include a format specifier for
      the vector array index (of type ``int``) **only** e.g.,

      .. code-block:: C

         SUNLogExtraDebugVecArray(logger, "YS-vector-array", "YS[%d](:) =", YS, 5);


Struct Accessor Macros
----------------------

Since many SUNDIALS structs use a type-erased (i.e., `void*`) "content" pointer, 
a common idiom occurring in SUNDIALS code is extracting the content, casting it to its original
type, and then accessing the struct member of interest. To ensure readability, it is 
recommended to use locally (to the source file in question) defined macros `GET_CONTENT`
and `IMPL_MEMBER` like the following example:

.. code-block:: c

   #define GET_CONTENT(S)       ((SUNAdjointCheckpointScheme_Fixed_Content)S->content)
   #define IMPL_MEMBER(S, prop) (GET_CONTENT(S)->prop)

   SUNAdjointCheckpointScheme self;
   IMPL_MEMBER(self, current_insert_step_node)   = step_data_node;
   IMPL_MEMBER(self, step_num_of_current_insert) = step_num;
