..
   Author(s): David J. Gardner @ LLNL
   -----------------------------------------------------------------------------
   SUNDIALS Copyright Start
   Copyright (c) 2002-2023, Lawrence Livermore National Security
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


Documenting Functions
=====================

Adding New Functions
--------------------

The documentation for new functions should include the ``.. versionadded::``
directive at the end of the documentation text noting the *package version*
number in which the function was added.

Changes to Existing Functions
-----------------------------

If the signature or behavior of a function changes in any release the
``.. versionchanged::`` directive should be added to the function documentation
noting the *package version* number in which the change happened and describing
the change.

Deprecating Functions
---------------------

When a function is deprecated the ``.. deprecated::`` directive should be added
to the function documentation noting the *package version* number in which the
function was deprecated and describing what function should be used instead
if appropriate.
