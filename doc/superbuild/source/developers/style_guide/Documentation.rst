..
   Author(s): David J. Gardner @ LLNL
   -----------------------------------------------------------------------------
   SUNDIALS Copyright Start
   Copyright (c) 2002-2024, Lawrence Livermore National Security
   and Southern Methodist University.
   All rights reserved.

   See the top-level LICENSE and NOTICE files for details.

   SPDX-License-Identifier: BSD-3-Clause
   SUNDIALS Copyright End
   -----------------------------------------------------------------------------

.. _Style.Documentation:

*******************
Documentation Style
*******************

Style guide for reStructuredText with Sphinx. For the most part, we attempt
to follow the Python developer's guide commentary on documentation (where relevant):
https://devguide.python.org/documenting/#style-guide.

File Structure
==============

Each package in SUNDIALS has a directory within the `doc` folder,
because each package documentation can be compiled separately.
Shared documentation goes into the `doc/shared` directory and
can be included from the package documentation using the
Sphinx ``.. include::`` directive. Generally speaking, the use
of symlinks to the `doc/shared` directory is discouraged.
However, for the `figs` folder, we allow it.

Figures and other assets
------------------------

Figures and other assets go under `doc/shared/figs`. Package specific
figures go in the `doc/shared/figs/<package>` directory.

Documentation superbuild
------------------------

An HTML build of the documentation can be generated that includes
all of the package documentations (with shared parts unrepeated).
This is what is built by readthedocs and what can be accessed
at readthedocs.org/project/sundials.

Headings
========

Follow the Python documentation convention:

.. code-block::

   # with overline, for parts

   * with overline, for chapters

   =, for sections

   -, for subsections

   ^, for subsubsections

   ", for paragraphs

Note that by following this convention, the `#` headers should **never**
occur in the package documentation directories. This is reserved for
the documentation superbuild.

Colors
======

+---------------------+--------------------+-----------------------+-----------------------------+
| Color Roles         | Text               | Color Italic Roles    | Text                        |
+---------------------+--------------------+-----------------------+-----------------------------+
| ``:red:`text```     | :red:`red`         | ``:redit:`text```     | :redit:`red italic`         |
+---------------------+--------------------+-----------------------+-----------------------------+
| ``:black:`text```   | :black:`black`     | ``:blackit:`text```   | :blackit:`black italic`     |
+---------------------+--------------------+-----------------------+-----------------------------+
| ``:gray:`text```    | :gray:`gray`       | ``:grayit:`text```    | :grayit:`gray italic`       |
+---------------------+--------------------+-----------------------+-----------------------------+
| ``:silver:`text```  | :silver:`silver`   | ``:silverit:`text```  | :silverit:`silver italic`   |
+---------------------+--------------------+-----------------------+-----------------------------+
| ``:magenta:`text``` | :magenta:`magenta` | ``:magentait:`text``` | :magentait:`magenta italic` |
+---------------------+--------------------+-----------------------+-----------------------------+
| ``:pink:`text```    | :pink:`pink`       | ``:pinkit:`text```    | :pinkit:`pink italic`       |
+---------------------+--------------------+-----------------------+-----------------------------+
| ``:orange:`text```  | :orange:`orange`   | ``:orangeit:`text```  | :orangeit:`orange italic`   |
+---------------------+--------------------+-----------------------+-----------------------------+
| ``:yellow:`text```  | :yellow:`yellow`   | ``:yellowit:`text```  | :yellowit:`yellow italic`   |
+---------------------+--------------------+-----------------------+-----------------------------+
| ``:lime:`text```    | :lime:`lime`       | ``:limeit:`text```    | :limeit:`lime italic`       |
+---------------------+--------------------+-----------------------+-----------------------------+
| ``:green:`text```   | :green:`green`     | ``:greenit:`text```   | :greenit:`green italic`     |
+---------------------+--------------------+-----------------------+-----------------------------+
| ``:olive:`text```   | :olive:`olive`     | ``:oliveit:`text```   | :oliveit:`olive italic`     |
+---------------------+--------------------+-----------------------+-----------------------------+
| ``:teal:`text```    | :teal:`teal`       | ``:tealit:`text```    | :tealit:`teal italic`       |
+---------------------+--------------------+-----------------------+-----------------------------+
| ``:cyan:`text```    | :cyan:`cyan`       | ``:cyanit:`text```    | :cyanit:`cyan italic`       |
+---------------------+--------------------+-----------------------+-----------------------------+
| ``:blue:`text```    | :blue:`blue`       | ``:blueit:`text```    | :blueit:`blue italic`       |
+---------------------+--------------------+-----------------------+-----------------------------+
| ``:purple:`text```  | :purple:`purple`   | ``:purpleit:`text```  | :purpleit:`purple italic`   |
+---------------------+--------------------+-----------------------+-----------------------------+

Capitalization
==============

Special terms in the SUNDIALS documentation that should be capitalized:
TODO: enumerate them

Footnotes
=========

Sphinx footnotes do not compile when generating the PDF from Latex,
therefore the use of footnotes is entirely banned. Restructure the
text, use ``.. notes::``, or ``.. warning::`` directives instead.

References
==========

All citations go into `doc/shared/sundials.bib`.
TODO: add citation and reference key style.


.. _Style.Documentation.UserCallable:

User-Callable Functions
=======================

Document user-callable functions with the :external+sphinx:rst:dir:`c:function`
or :external+sphinx:rst:dir:`cpp:function` directives, as appropriate. The
directive is followed by the C/C++ function signature. Under the signature
(skipping a line and indenting 3 spaces) provide a brief description of the
function followed by any information on its usage. When describing the function
parameters use ``:param <name>:``. If the function returns a specific set of
SUNDIALS error codes, describe the possible return values using ``:retval
<value>:`` for each value. Otherwise, use a single ``:returns:`` item to
describe the result of the function. If the function returns ``void``, a return
entry should not be included. Below we give two examples, the first returns an
error code (int) while the second is a constructor that returns an object
(pointer).

.. code-block:: rst

   .. c:function:: int SetFoo(param_type p1, param_type p2)

      Brief description of what the user-callable function does.

      Additional information about the function and its usage.

      :param p1: description of the first parameter.
      :param p2: description of the second parameter.

      :retval SUCCESS: under some conditions.
      :retval FAILURE_1: under some other conditions.
      :retval FAILURE_2: under yet some other conditions.

.. code-block:: rst

   .. c:function:: FooObject CreateFooObject(param_type p1, param_type p2)

      Brief description of what the user-callable function does.

      Additional information about the function and its usage.

      :param p1: description of the first parameter.
      :param p2: description of the second parameter.

      :returns: If successful some object, otherwise ``NULL``.

When adding, updating, or deprecating a function, use the
:external+sphinx:rst:dir:`versionadded`,
:external+sphinx:rst:dir:`versionchanged`, or
:external+sphinx:rst:dir:`deprecated` directives with the placeholder version
number ``x.y.z`` after the return description. The release script will find and
replace all instances of ``x.y.z`` in the documentation with the actual release
number. When altering the behavior of a function or deprecating a function
include a description for the change under the directive (skipping a line and
indenting 3 spaces). For example,

.. code-block:: rst

   .. versionadded:: x.y.z

.. code-block:: rst

   .. versionchanged:: x.y.z

      Describe how the function behavior has changed from before.

.. code-block:: rst

   .. deprecated:: x.y.z

      If a replacement function/procedure is available, describe what users
      should do to replace the deprecated function e.g., cross reference the
      function superseding this one or list the new steps to follow. Otherwise,
      note the feature/capability is no longer supported/provided and, if
      possible, state why this function was removed.

If special attention needs to be drawn to some behavior, consideration, or
limitation of a function that could be overlooked in the description, use the
:external+sphinx:rst:dir:`note` or :external+sphinx:rst:dir:`warning` directives
as appropriate. These should be used sparingly to avoid diluting their impact.
For example,

.. code-block:: rst

   .. note::

      Something users should not over look e.g., a feature is only compatible
      with a subset of methods.

.. code-block:: rst

   .. warning::

      Something critical users should be aware of e.g., performance impacts.

Finally, at the end of the function documentation, you may include (a
non-trivial) example usage of the function and/or a list of example programs
that utilize the function. For example,

.. code-block:: rst

   **Example usage:**

   .. code-block:: C

      /* Short code block demonstrating typical usage */

      /* Create the object */
      FooObject foo_obj = CreateFooObject(p1, p2);
      if (foo_obj == NULL) { return 1; }

      /* Attach the object to mem */
      int retval = SetFoo(mem, foo_obj);
      if (retval != 0) { return 1; }

      /* Perform some actions */
      ...

      /* Destroy the object */
      retval = DestroyFooObject(&foo_obj);
      if (retval != 0) { return 1; }

.. code-block:: rst

   **Examples codes:**

   * ``examples/package/subdir/pkg_some_code.c``

Putting it all together, the rendered documentation should look like the
following.

.. c:function:: int FooSetBar(void* foo_obj, int bar_value)
   :nocontentsentry:
   :noindexentry:

   This function sets the value of Bar in a FooObject.

   The default value for Bar is :math:`10`. An input value :math:`< 0` will
   reset Bar to the default value.

   :param foo_obj: the FooObject.
   :param bar_value: the value of Bar.

   :retval SUCCESS: if the value was successfully set.
   :retval NULL_OBJ: if the ``foo_obj`` was ``NULL``.

   .. versionadded:: 1.1.0

   .. versionchanged:: 2.0.0

      The type of p1 was changed from ``unsigned int`` to ``int``

   .. note::

      Utilizing this capability requires building with Bar enabled.

   .. warning::

      Setting values greater than 100 may degrade performance.

   **Example usage:**

   .. code-block:: C

      /* Create the object */
      void* foo_obj = CreateFooObject(p1, p2);
      if (foo_obj == NULL) { return 1; }

      /* Update the value of Bar */
      int retval = FooSetBar(foo_obj, 50);
      if (retval != 0) { return 1; }

      /* Perform some actions */
      ...

   **Examples codes:**

   * ``examples/package/subdir/pkg_foo_demo.c``


.. _Style.Documentation.UserSupplied:

User-Supplied Functions
=======================

Document user-supplied functions with the :external+sphinx:rst:dir:`c:type`
directive. The directive is followed by the ``typedef`` for the function
pointer. The description of the function type mirrors the style used for
user-callable functions (see :ref:`Style.Documentation.UserCallable`) with one
exception. As :external+sphinx:rst:dir:`c:type` does not currently support the
``param``, ``retval``, and ``returns`` fields, these sections must be manually
created. The style that follows is chosen to reflect that of ``param``,
``retval``, and ``returns`` fields as much as possible. Function parameters
should be listed under a boldface "Parameters:" section with the parameters in
boldface and separated from their description by an en-dash. As user-supplied
functions typically return a ``int``, but specific values are not required, a
description of how the return value is interpreted should be given under a
boldface "Returns:" section (skipping a line and indenting 2 spaces). If
specific return values are required, these should be documented similarly to the
function parameters and listed under a boldface "Return values:" section. If the
function returns ``void``, a return section should not be included. Below we
give, two examples describing user-supplied functions.

.. code-block:: rst

   .. c:type:: int (*FooFn)(param_type p1, param_type p2)

      Brief description of what the user-provided function should do.

      Additional information about the function and its usage.

      **Parameters:**

      * **p1** -- description of the first parameter.
      * **p2** -- description of the second parameter.

      **Returns:**

        A :c:type:`FooFn` function should return 0 if successful, a positive
        value if a recoverable error occurred, or a negative value if an
        unrecoverable error occurred.

.. code-block:: rst

   .. c:type:: int (*BarFn)(param_type p1, param_type p2)

      Brief description of what the user-provided function should do.

      Additional information about the function and its usage.

      **Parameters:**

      * **p1** -- description of the first parameter.
      * **p2** -- description of the second parameter.

      **Return values:**

      * **VALUE_1** -- under some circumstances.
      * **VALUE_2** -- under some other circumstances.

Other than the difference in the function parameter and return value sections
the remaining guidelines from the user-callable function documentation are the
same. Putting it all together, the rendered documentation should look like the
following.

.. c:type:: int (*FooFn)(double* p1, double* p2)
   :nocontentsentry:
   :noindexentry:

   Brief description of what the user-provided function should do.

   Additional information about the function and its usage.

   **Parameters:**

   * **p1** -- the input array of values.
   * **p2** -- the output array of values.

   **Returns:**

     A :c:type:`FooFn` function should return 0 if successful, a positive value
     if a recoverable error occurred, or a negative value if an unrecoverable
     error occurred.

   .. versionadded:: 2.2.0

   .. note::

      This function is required when using the Foo option.

   **Examples codes:**

   * ``examples/package/subdir/pkg_bar_demo.c``


.. _Style.Documentation.MacroFunction:

Function-like Macros
====================

Document function-like macros with the :external+sphinx:rst:dir:`c:macro`
directive followed by the macro. The guidelines for documenting function-like
macros are the same as those used for documenting user-callable functions (see
:ref:`Style.Documentation.UserCallable`) with one exception. As
:external+sphinx:rst:dir:`c:macro` does not include the parameter types, the
types should be included in the parameter descriptions when relevant i.e., when
the macro is a wrapper to function (see :c:macro:`SUNLogError`). For example,

.. code-block:: rst

   .. c:macro:: FnLikeMacro(p1, p2)

      Brief description of what the function-like macro does.

      Additional information about the macro and its usage.

      :param p1: the :c:type:`p1_type` parameter.
      :param p2: the :c:type:`p1_type` parameter.

      :retval retval1: under some conditions.
      :retval retval2: under some other conditions.

      .. versionadded:: x.y.z
