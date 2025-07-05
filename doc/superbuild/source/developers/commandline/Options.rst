..
   Author(s): Daniel R. Reynolds @ UMBC
   -----------------------------------------------------------------------------
   SUNDIALS Copyright Start
   Copyright (c) 2002-2025, Lawrence Livermore National Security
   and Southern Methodist University.
   All rights reserved.

   See the top-level LICENSE and NOTICE files for details.

   SPDX-License-Identifier: BSD-3-Clause
   SUNDIALS Copyright End
   -----------------------------------------------------------------------------

.. _Infrastructure:

Infrastructure
==============

SUNDIALS provides support for parsing command-line arguments to control SUNDIALS "Set"
routines through a set of data structures and utility routines that compare a given
command-line argument against a list of reserved keywords.  If the given argument is
matched to a relevant key, then the associated "Set" routine is called with the
argument(s) that follow the key on the command line.  This allows users to set SUNDIALS
options at runtime, without needing to modify the source code or recompile the library.

Prototypes for the SUNDIALS-provided infrastructure routines and corresponding data
structures for command-line support are in the ``src/sundials/sundials_cli.h``
header file. Their implementations are found in ``src/sundials/sundials_cli.c``.
These routines assume that each SUNDIALS module provides a suite of
"Set" routines of the form ``<SetFn>(void* mem, arg1, arg2, ...)``, where
``<SetFn>`` is the name of the set routine (e.g., :c:func:`ARKodeSetOrder`), ``mem``
is an opaque pointer to the SUNDIALS module (e.g., the pointer returned by
:c:func:`CVodeCreate`), and ``arg1``, ..., are the input arguments for the set routine.
The SUNDIALS-provided command-line processing routines and structures are grouped
according to the number and type of these arguments, as elaborated in the
following sections:

* :ref:`IntArg` - single ``int`` argument
* :ref:`TwoIntArgs` - two ``int`` arguments
* :ref:`LongArg` - single ``long`` argument
* :ref:`IntRealArgs` - a single ``int`` followed by a single ``sunrealtype`` argument
* :ref:`IntRealRealArgs` - a single ``int`` followed by two ``sunrealtype`` arguments
* :ref:`IntLongArgs` - a single ``int`` followed by a single ``long`` argument
* :ref:`RealArg` - single ``sunrealtype`` argument
* :ref:`TwoRealArgs` - two ``sunrealtype`` arguments
* :ref:`CharArg` - single ``char*`` argument
* :ref:`TwoCharArgs` - two ``char*`` arguments
* :ref:`ActionArg` - perform a given action

.. _IntArg:

One ``int`` argument
--------------------

.. c:type:: int (*sunIntSetFn)(void* mem, int arg1)

   These functions set a single integer option for a SUNDIALS module.

   :param mem: an opaque pointer to the SUNDIALS module.
   :param arg1: the integer value to set.

   :return: An *sunIntSetFn* should return 0 if successful, or a nonzero value on failure.

   .. versionadded:: x.y.z

.. c:struct:: sunKeyIntPair

   This is a structure that contains

   .. c:member:: const char* key;

         The command-line key to match.

   .. c:member:: sunIntSetFn set;

         The function to call if the key is matched.

   .. versionadded:: x.y.z

.. c:function:: SUNErrCode sunCheckAndSetIntArgs(void* mem, int* argidx, char* argv[], const size_t offset, const struct sunKeyIntPair* testpairs, int numpairs, sunbooleantype* arg_used, int *failedarg)

   This function loops over an array of potential key/function pairs to check whether any match ``argv[*argidx]``, and if so it calls the corresponding set routine.

   :param mem: an opaque pointer to the SUNDIALS module.
   :param argidx: a pointer to the index of the current command-line argument.  If the argument is found and set, this will be incremented by the number of arguments consumed (e.g., 1 for a single argument, 2 for two arguments, etc.).
   :param argv: the command-line argument vector.
   :param offset: the offset width to ignore (stores a module-specific prefix for the key).
   :param testpairs: an array of key-value pairs to test against.
   :param numpairs: the number of key-value pairs in ``testpairs``.
   :param arg_used: a pointer to a boolean indicating if the argument was used.
   :param failedarg: a pointer to an integer indicating the failed argument index (if any).

   :return: SUN_SUCCESS if either the argument was not found, or if it was matched and set correctly.  If it was found but the set routine failed, then this returns the value emanating from the module-specific set routine.

   .. versionadded:: x.y.z

.. _TwoIntArgs:

Two ``int`` arguments
---------------------

.. c:type:: int (*sunTwoIntSetFn)(void* mem, int arg1, int arg2)

   These functions set two integer options for a SUNDIALS module.

   :param mem: an opaque pointer to the SUNDIALS module.
   :param arg1: the first integer value to set.
   :param arg2: the second integer value to set.

   :return: An *sunTwoIntSetFn* should return 0 if successful, or a nonzero value on failure.

   .. versionadded:: x.y.z

.. c:struct:: sunKeyTwoIntPair

   This is a structure that contains

   .. c:member:: const char* key;

         The command-line key to match.

   .. c:member:: sunTwoIntSetFn set;

         The function to call if the key is matched.

   .. versionadded:: x.y.z

.. c:function:: SUNErrCode sunCheckAndSetTwoIntArgs(void* mem, int* argidx, char* argv[], const size_t offset, const struct sunKeyTwoIntPair* testpairs, int numpairs, sunbooleantype* arg_used, int *failedarg)

   This function loops over an array of potential key/function pairs to check whether any match ``argv[*argidx]``, and if so it calls the corresponding set routine.

   :param mem: an opaque pointer to the SUNDIALS module.
   :param argidx: a pointer to the index of the current command-line argument.  If the argument is found and set, this will be incremented by the number of arguments consumed (e.g., 1 for a single argument, 2 for two arguments, etc.).
   :param argv: the command-line argument vector.
   :param offset: the offset width to ignore (stores a module-specific prefix for the key).
   :param testpairs: an array of key-value pairs to test against.
   :param numpairs: the number of key-value pairs in ``testpairs``.
   :param arg_used: a pointer to a boolean indicating if the argument was used.
   :param failedarg: a pointer to an integer indicating the failed argument index (if any).

   :return: SUN_SUCCESS if either the argument was not found, or if it was matched and set correctly.  If it was found but the set routine failed, then this returns the value emanating from the module-specific set routine.

   .. versionadded:: x.y.z

.. _LongArg:

One ``long int`` argument
-------------------------

.. c:type:: int (*sunLongSetFn)(void* mem, long int arg1)

   These functions set a single long integer option for a SUNDIALS module.

   :param mem: an opaque pointer to the SUNDIALS module.
   :param arg1: the long integer value to set.

   :return: An *sunLongSetFn* should return 0 if successful, or a nonzero value on failure.

   .. versionadded:: x.y.z

.. c:struct:: sunKeyLongPair

   This is a structure that contains

   .. c:member:: const char* key;

         The command-line key to match.

   .. c:member:: sunLongSetFn set;

         The function to call if the key is matched.

   .. versionadded:: x.y.z

.. c:function:: SUNErrCode sunCheckAndSetLongArgs(void* mem, int* argidx, char* argv[], const size_t offset, const struct sunKeyLongPair* testpairs, int numpairs, sunbooleantype* arg_used, int *failedarg)

   This function loops over an array of potential key/function pairs to check whether any match ``argv[*argidx]``, and if so it calls the corresponding set routine.

   :param mem: an opaque pointer to the SUNDIALS module.
   :param argidx: a pointer to the index of the current command-line argument.  If the argument is found and set, this will be incremented by the number of arguments consumed (e.g., 1 for a single argument, 2 for two arguments, etc.).
   :param argv: the command-line argument vector.
   :param offset: the offset width to ignore (stores a module-specific prefix for the key).
   :param testpairs: an array of key-value pairs to test against.
   :param numpairs: the number of key-value pairs in ``testpairs``.
   :param arg_used: a pointer to a boolean indicating if the argument was used.
   :param failedarg: a pointer to an integer indicating the failed argument index (if any).

   :return: SUN_SUCCESS if either the argument was not found, or if it was matched and set correctly.  If it was found but the set routine failed, then this returns the value emanating from the module-specific set routine.

   .. versionadded:: x.y.z

.. _IntRealArgs:

One ``int`` and one ``sunrealtype`` argument
--------------------------------------------

.. c:type:: int (*sunIntRealSetFn)(void* mem, int arg1, sunrealtype arg2)

   These functions set a single integer option and a single real option for a SUNDIALS module.

   :param mem: an opaque pointer to the SUNDIALS module.
   :param arg1: the integer value to set.
   :param arg2: the real value to set.

   :return: An *sunIntRealSetFn* should return 0 if successful, or a nonzero value on failure.

   .. versionadded:: x.y.z

.. c:struct:: sunKeyIntRealPair

   This is a structure that contains

   .. c:member:: const char* key;

         The command-line key to match.

   .. c:member:: sunIntRealSetFn set;

         The function to call if the key is matched.

   .. versionadded:: x.y.z

.. c:function:: SUNErrCode sunCheckAndSetIntRealArgs(void* mem, int* argidx, char* argv[], const size_t offset, const struct sunKeyIntRealPair* testpairs, int numpairs, sunbooleantype* arg_used, int *failedarg)

   This function loops over an array of potential key/function pairs to check whether any match ``argv[*argidx]``, and if so it calls the corresponding set routine.

   :param mem: an opaque pointer to the SUNDIALS module.
   :param argidx: a pointer to the index of the current command-line argument.  If the argument is found and set, this will be incremented by the number of arguments consumed (e.g., 1 for a single argument, 2 for two arguments, etc.).
   :param argv: the command-line argument vector.
   :param offset: the offset width to ignore (stores a module-specific prefix for the key).
   :param testpairs: an array of key-value pairs to test against.
   :param numpairs: the number of key-value pairs in ``testpairs``.
   :param arg_used: a pointer to a boolean indicating if the argument was used.
   :param failedarg: a pointer to an integer indicating the failed argument index (if any).

   :return: SUN_SUCCESS if either the argument was not found, or if it was matched and set correctly.  If it was found but the set routine failed, then this returns the value emanating from the module-specific set routine.

   .. versionadded:: x.y.z

.. _IntRealRealArgs:

One ``int`` and two ``sunrealtype`` arguments
---------------------------------------------

.. c:type:: int (*sunIntRealRealSetFn)(void* mem, int arg1, sunrealtype arg2, sunrealtype arg3)

   These functions set a single integer option and two real options for a SUNDIALS module.

   :param mem: an opaque pointer to the SUNDIALS module.
   :param arg1: the integer value to set.
   :param arg2: the first real value to set.
   :param arg3: the second real value to set.

   :return: An *sunIntRealRealSetFn* should return 0 if successful, or a nonzero value on failure.

   .. versionadded:: x.y.z

.. c:struct:: sunKeyIntRealRealPair

   This is a structure that contains

   .. c:member:: const char* key;

         The command-line key to match.

   .. c:member:: sunIntRealRealSetFn set;

         The function to call if the key is matched.

   .. versionadded:: x.y.z

.. c:function:: SUNErrCode sunCheckAndSetIntRealRealArgs(void* mem, int* argidx, char* argv[], const size_t offset, const struct sunKeyIntRealRealPair* testpairs, int numpairs, sunbooleantype* arg_used, int *failedarg)

   This function loops over an array of potential key/function pairs to check whether any match ``argv[*argidx]``, and if so it calls the corresponding set routine.

   :param mem: an opaque pointer to the SUNDIALS module.
   :param argidx: a pointer to the index of the current command-line argument.  If the argument is found and set, this will be incremented by the number of arguments consumed (e.g., 1 for a single argument, 2 for two arguments, etc.).
   :param argv: the command-line argument vector.
   :param offset: the offset width to ignore (stores a module-specific prefix for the key).
   :param testpairs: an array of key-value pairs to test against.
   :param numpairs: the number of key-value pairs in ``testpairs``.
   :param arg_used: a pointer to a boolean indicating if the argument was used.
   :param failedarg: a pointer to an integer indicating the failed argument index (if any).

   :return: SUN_SUCCESS if either the argument was not found, or if it was matched and set correctly.  If it was found but the set routine failed, then this returns the value emanating from the module-specific set routine.

   .. versionadded:: x.y.z

.. _IntLongArgs:

One ``int`` and one ``long int`` argument
-----------------------------------------

.. c:type:: int (*sunIntLongSetFn)(void* mem, int arg1, long int arg2)

   These functions set a single integer option and a long integer option for a SUNDIALS module.

   :param mem: an opaque pointer to the SUNDIALS module.
   :param arg1: the integer value to set.
   :param arg2: the long integer value to set.

   :return: An *sunIntLongSetFn* should return 0 if successful, or a nonzero value on failure.

   .. versionadded:: x.y.z

.. c:struct:: sunKeyIntLongPair

   This is a structure that contains

   .. c:member:: const char* key;

         The command-line key to match.

   .. c:member:: sunIntLongSetFn set;

         The function to call if the key is matched.

   .. versionadded:: x.y.z

.. c:function:: SUNErrCode sunCheckAndSetIntLongArgs(void* mem, int* argidx, char* argv[], const size_t offset, const struct sunKeyIntLongPair* testpairs, int numpairs, sunbooleantype* arg_used, int *failedarg)

   This function loops over an array of potential key/function pairs to check whether any match ``argv[*argidx]``, and if so it calls the corresponding set routine.

   :param mem: an opaque pointer to the SUNDIALS module.
   :param argidx: a pointer to the index of the current command-line argument.  If the argument is found and set, this will be incremented by the number of arguments consumed (e.g., 1 for a single argument, 2 for two arguments, etc.).
   :param argv: the command-line argument vector.
   :param offset: the offset width to ignore (stores a module-specific prefix for the key).
   :param testpairs: an array of key-value pairs to test against.
   :param numpairs: the number of key-value pairs in ``testpairs``.
   :param arg_used: a pointer to a boolean indicating if the argument was used.
   :param failedarg: a pointer to an integer indicating the failed argument index (if any).

   :return: SUN_SUCCESS if either the argument was not found, or if it was matched and set correctly.  If it was found but the set routine failed, then this returns the value emanating from the module-specific set routine.

   .. versionadded:: x.y.z

.. _RealArg:

One ``sunrealtype`` argument
----------------------------

.. c:type:: int (*sunRealSetFn)(void* mem, sunrealtype arg1)

   These functions set a single real option for a SUNDIALS module.

   :param mem: an opaque pointer to the SUNDIALS module.
   :param arg1: the real value to set.

   :return: An *sunRealSetFn* should return 0 if successful, or a nonzero value on failure.

   .. versionadded:: x.y.z

.. c:struct:: sunKeyRealPair

   This is a structure that contains

   .. c:member:: const char* key;

         The command-line key to match.

   .. c:member:: sunRealSetFn set;

         The function to call if the key is matched.

   .. versionadded:: x.y.z

.. c:function:: SUNErrCode sunCheckAndSetRealArgs(void* mem, int* argidx, char* argv[], const size_t offset, const struct sunKeyRealPair* testpairs, int numpairs, sunbooleantype* arg_used, int *failedarg)

   This function loops over an array of potential key/function pairs to check whether any match ``argv[*argidx]``, and if so it calls the corresponding set routine.

   :param mem: an opaque pointer to the SUNDIALS module.
   :param argidx: a pointer to the index of the current command-line argument.  If the argument is found and set, this will be incremented by the number of arguments consumed (e.g., 1 for a single argument, 2 for two arguments, etc.).
   :param argv: the command-line argument vector.
   :param offset: the offset width to ignore (stores a module-specific prefix for the key).
   :param testpairs: an array of key-value pairs to test against.
   :param numpairs: the number of key-value pairs in ``testpairs``.
   :param arg_used: a pointer to a boolean indicating if the argument was used.
   :param failedarg: a pointer to an integer indicating the failed argument index (if any).

   :return: SUN_SUCCESS if either the argument was not found, or if it was matched and set correctly.  If it was found but the set routine failed, then this returns the value emanating from the module-specific set routine.

   .. versionadded:: x.y.z

.. _TwoRealArgs:

Two ``sunrealtype`` arguments
-----------------------------

.. c:type:: int (*sunTwoRealSetFn)(void* mem, sunrealtype arg1, sunrealtype arg2)

   These functions set two real options for a SUNDIALS module.

   :param mem: an opaque pointer to the SUNDIALS module.
   :param arg1: the first real value to set.
   :param arg2: the second real value to set.

   :return: An *sunTwoRealSetFn* should return 0 if successful, or a nonzero value on failure.

   .. versionadded:: x.y.z

.. c:struct:: sunKeyTwoRealPair

   This is a structure that contains

   .. c:member:: const char* key;

         The command-line key to match.

   .. c:member:: sunTwoRealSetFn set;

         The function to call if the key is matched.

   .. versionadded:: x.y.z

.. c:function:: SUNErrCode sunCheckAndSetTwoRealArgs(void* mem, int* argidx, char* argv[], const size_t offset, const struct sunKeyTwoRealPair* testpairs, int numpairs, sunbooleantype* arg_used, int *failedarg)

   This function loops over an array of potential key/function pairs to check whether any match ``argv[*argidx]``, and if so it calls the corresponding set routine.

   :param mem: an opaque pointer to the SUNDIALS module.
   :param argidx: a pointer to the index of the current command-line argument.  If the argument is found and set, this will be incremented by the number of arguments consumed (e.g., 1 for a single argument, 2 for two arguments, etc.).
   :param argv: the command-line argument vector.
   :param offset: the offset width to ignore (stores a module-specific prefix for the key).
   :param testpairs: an array of key-value pairs to test against.
   :param numpairs: the number of key-value pairs in ``testpairs``.
   :param arg_used: a pointer to a boolean indicating if the argument was used.
   :param failedarg: a pointer to an integer indicating the failed argument index (if any).

   :return: SUN_SUCCESS if either the argument was not found, or if it was matched and set correctly.  If it was found but the set routine failed, then this returns the value emanating from the module-specific set routine.

   .. versionadded:: x.y.z

.. _CharArg:

One ``char*`` argument
----------------------

.. c:type:: int (*sunCharSetFn)(void* mem, const char* arg1)

   These functions set a single string option for a SUNDIALS module.

   :param mem: an opaque pointer to the SUNDIALS module.
   :param arg1: the string value to set.

   :return: An *sunCharSetFn* should return 0 if successful, or a nonzero value on failure.

   .. versionadded:: x.y.z

.. c:struct:: sunKeyCharPair

   This is a structure that contains

   .. c:member:: const char* key;

         The command-line key to match.

   .. c:member:: sunCharSetFn set;

         The function to call if the key is matched.

   .. versionadded:: x.y.z

.. c:function:: SUNErrCode sunCheckAndSetCharArgs(void* mem, int* argidx, char* argv[], const size_t offset, const struct sunKeyCharPair* testpairs, int numpairs, sunbooleantype* arg_used, int *failedarg)

   This function loops over an array of potential key/function pairs to check whether any match ``argv[*argidx]``, and if so it calls the corresponding set routine.

   :param mem: an opaque pointer to the SUNDIALS module.
   :param argidx: a pointer to the index of the current command-line argument.  If the argument is found and set, this will be incremented by the number of arguments consumed (e.g., 1 for a single argument, 2 for two arguments, etc.).
   :param argv: the command-line argument vector.
   :param offset: the offset width to ignore (stores a module-specific prefix for the key).
   :param testpairs: an array of key-value pairs to test against.
   :param numpairs: the number of key-value pairs in ``testpairs``.
   :param arg_used: a pointer to a boolean indicating if the argument was used.
   :param failedarg: a pointer to an integer indicating the failed argument index (if any).

   :return: SUN_SUCCESS if either the argument was not found, or if it was matched and set correctly.  If it was found but the set routine failed, then this returns the value emanating from the module-specific set routine.

   .. versionadded:: x.y.z

.. _TwoCharArgs:

Two ``char*`` arguments
-----------------------

.. c:type:: int (*sunTwoCharSetFn)(void* mem, const char* arg1, const char* arg2)

   These functions set two string options for a SUNDIALS module.

   :param mem: an opaque pointer to the SUNDIALS module.
   :param arg1: the first string value to set.
   :param arg2: the second string value to set.

   :return: An *sunTwoCharSetFn* should return 0 if successful, or a nonzero value on failure.

   .. versionadded:: x.y.z

.. c:struct:: sunKeyTwoCharPair

   This is a structure that contains

   .. c:member:: const char* key;

         The command-line key to match.

   .. c:member:: sunTwoCharSetFn set;

         The function to call if the key is matched.

   .. versionadded:: x.y.z

.. c:function:: SUNErrCode sunCheckAndSetTwoCharArgs(void* mem, int* argidx, char* argv[], const size_t offset, const struct sunKeyTwoCharPair* testpairs, int numpairs, sunbooleantype* arg_used, int *failedarg)

   This function loops over an array of potential key/function pairs to check whether any match ``argv[*argidx]``, and if so it calls the corresponding set routine.

   :param mem: an opaque pointer to the SUNDIALS module.
   :param argidx: a pointer to the index of the current command-line argument.  If the argument is found and set, this will be incremented by the number of arguments consumed (e.g., 1 for a single argument, 2 for two arguments, etc.).
   :param argv: the command-line argument vector.
   :param offset: the offset width to ignore (stores a module-specific prefix for the key).
   :param testpairs: an array of key-value pairs to test against.
   :param numpairs: the number of key-value pairs in ``testpairs``.
   :param arg_used: a pointer to a boolean indicating if the argument was used.
   :param failedarg: a pointer to an integer indicating the failed argument index (if any).

   :return: SUN_SUCCESS if either the argument was not found, or if it was matched and set correctly.  If it was found but the set routine failed, then this returns the value emanating from the module-specific set routine.

   .. versionadded:: x.y.z

.. _ActionArg:

No arguments (action only)
--------------------------

.. c:type:: int (*sunActionSetFn)(void* mem)

   These functions set a single integer option for a SUNDIALS module.

   :param mem: an opaque pointer to the SUNDIALS module.

   :return: An *sunActionSetFn* should return 0 if successful, or a nonzero value on failure.

   .. versionadded:: x.y.z

.. c:struct:: sunKeyActionPair

   This is a structure that contains

   .. c:member:: const char* key;

         The command-line key to match.

   .. c:member:: sunActionSetFn set;

         The function to call if the key is matched.

   .. versionadded:: x.y.z

.. c:function:: SUNErrCode sunCheckAndSetActionArgs(void* mem, int* argidx, char* argv[], const size_t offset, const struct sunKeyActionPair* testpairs, int numpairs, sunbooleantype* arg_used, int *failedarg)

   This function loops over an array of potential key/function pairs to check whether any match ``argv[*argidx]``, and if so it calls the corresponding set routine.

   :param mem: an opaque pointer to the SUNDIALS module.
   :param argidx: a pointer to the index of the current command-line argument.  If the argument is found and set, this will be incremented by the number of arguments consumed (e.g., 1 for a single argument, 2 for two arguments, etc.).
   :param argv: the command-line argument vector.
   :param offset: the offset width to ignore (stores a module-specific prefix for the key).
   :param testpairs: an array of key-value pairs to test against.
   :param numpairs: the number of key-value pairs in ``testpairs``.
   :param arg_used: a pointer to a boolean indicating if the argument was used.
   :param failedarg: a pointer to an integer indicating the failed argument index (if any).

   :return: SUN_SUCCESS if either the argument was not found, or if it was matched and set correctly.  If it was found but the set routine failed, then this returns the value emanating from the module-specific set routine.

   .. versionadded:: x.y.z



Package-specific Command-line Support
=====================================

Each SUNDIALS module that wishes to support command-line options should provide a
routine of the form
``<module>SetFromCommandLine(void* mem, const char* moduleid, int argc, char* argv[])``.
This routine can then be called by users to indicate that they wish to use
command-line options to control the corresponding SUNDIALS module.  The arguments to
this function are:

* ``mem``: an opaque pointer to the SUNDIALS module (e.g., the pointer returned by
  :c:func:`CVodeCreate`).
* ``moduleid``: a desired string identifier for arguments to that module (e.g., "arkode").
  Note that each module should specify a default string identifier, that would be
  used if the user specifies ``NULL`` for this argument.  However, users can supply
  non-default identifiers so that they can control multiple instances of the same module
  to be independently (e.g., when using multiple ARKode integrators in the
  same program).
* ``argc``: the number of command-line arguments.
* ``argv``: the command-line argument vector.

Within this module-provided routine, arrays of key-value pairs having the correct type
for the corresponding "Set" routine should be defined (e.g., see the file
``src/arkode/arkode_cli.c``).

.. note::
   When adding new "Set" routines to an existing SUNDIALS module, developers should
   attempt to add a corresponding entry in the appropriate key-value pair array, and note
   the new key in the module's documentation.

After defining the allowable command-line arguments (and their corresponding "Set"
routines), the module-provided routine should loop over all ``argc`` command-line
arguments, and perform the following steps:

#. Check whether the prefix for the current command-line argument matches the module's
   identifier (e.g., "arkode").  If it does not match, then skip to the next argument.
#. If the prefix matches, then call each of the SUNDIALS-provided command-line processing
   routines above (e.g., :c:func:`sunCheckAndSetActionArgs`) to attempt processing of that
   command-line argument.  If that routine indicates that the argument was used, then
   continue to the next command-line argument; else the next SUNDIALS-provided
   command-line processing routine should be called.
#. If no SUNDIALS-provided command-line processing routine indicates that the argument
   was used, then the module routine can process additional arguments that fall outside
   the expertise of the SUNDIALS-provided routines.
#. By the end of the loop body, if a given argument that has the correct prefix has
   still not been processed, then the routine should print a warning that the argument
   was not handled.