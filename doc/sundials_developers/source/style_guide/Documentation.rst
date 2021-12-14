..
   Author(s): David J. Gardner @ LLNL
   -----------------------------------------------------------------------------
   SUNDIALS Copyright Start
   Copyright (c) 2002-2021, Lawrence Livermore National Security
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

Style guide for RestructuredText with Sphinx. For the most part, we attempt
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
This is what is built by readthedocs and what can be acessed
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

Capatilization
==============

Special terms in the SUNDIALS documentation that should be captialized:
TODO: enumerate them

Footnotes
=========

Sphinx footnotes do not compile when generating the PDF from Latex,
therefore the use of footnotes is entirely banned. Restructure the
text, use ``.. notes::``, or ``.. warning::`` directives instead.

References
----------

All citations go into `doc/shared/sundials.bib`.
TODO: add citation and reference key style.
