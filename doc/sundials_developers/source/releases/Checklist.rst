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

.. _Versioning:

Versioning
==========

SUNDIALS follows the `semantic versioning <https://semver.org/>`_ scheme and
each release is given a version number ``x.y.z`` where

* ``x`` is the major version number and is incremented when the release includes
  an incompatible API change

* ``y`` is the minor version number and is incremented when new functionality
  has been added and is backwards compatible with the prior version

* ``z`` is the patch version and is incremented when the release only includes
  backwards compatible bug fixes

Additionally a pre-release versions will have ``-l.w`` appended to the version
number i.e., ``x.y.z-l.w`` where

* ``l`` is a pre-release label e.g., ``alpha``, ``beta``, or ``rc``.

* ``w`` is the pre-release version e.g, ``alpha.0``, ``alpha.1``, etc.

The convention for shared library versioning is ``libfoo.so.x.y.z``. Note that
there was no library versioning before sundials-2.2.0. The version information
is specified by setting the ``VERSION`` and ``SOVERSION`` properties to the
library target

.. code-block:: cmake

   add_library(foo file1.c file2.c)
   set_target_properties(foo PROPERTIES VERSION x.y.z SOVERSION x)

and build system is responsible for creating the symlinks

.. code-block:: none

   libfoo.so   -> libfoo.so.x.y.z
   libfoo.so.x -> libfoo.so.x.y.z

Note that we force ``SOVERSION`` to be the same as the major number for
compatibility with the libtool versioning scheme (at least under Linux).

.. _ReleaseCheckList:

Pre-Release Tasks
=================

This is a list of tasks that need to be done before a SUNDIALS release.
The order is bottom-up starting with solver source files and ending with
web pages.

#. Create release branch from trunk in the Git repository

#. If this is a major release, search the SUNDIALS code for
   'DEPRECATION NOTICE' and 'SUNDIALS_DEPRECATED'. All deprecated
   functions should be removed.

#. Regenerate the Fortran 2003 interfaces. It is possible nothing will be updated.

#. Update the "Changes in ..." sections in all user guides. The changes should be
   sorted so that major new features are above bug fixes.

#. Update version numbers and release date information using the ``updateVerson.sh``
   script. This will update the following files:

   * ``CMakeLists.txt``
   * ``README.md``
   * ``src/arkode/README``
   * ``src/cvode/README``
   * ``src/cvodes/README``
   * ``src/ida/README``
   * ``src/idas/README``
   * ``src/kinsol/README``
   * ``doc/arkode/SphinxDocs/examples/source/conf.py``
   * ``doc/arkode/SphinxDocs/examples/source/References.rst``
   * ``doc/arkode/SphinxDocs/guide/source/conf.py``
   * ``doc/arkode/SphinxDocs/guide/source/History.rst``
   * ``doc/arkode/SphinxDocs/guide/source/References.rst``
   * ``doc/sundials/biblio.bib``
   * ``doc/sundials/sundials_release_history.tex``
   * ``doc/sundials/ug.tex``
   * ``scripts/tarscript``

   The ARKode LaTeX files ``doc/arkode/ARKode_example.tex`` and
   ``doc/arkode/ARKode.tex`` need to be regenerated and updated manually.

   The following files are no longer maintaianed:

   * ``html/main.html`` (This is no longer maintained as of at least 2016)
   * ``sundialsTB/install_STB.m`` (This is no longer maintained as of 2016)

#. Update version numbers of third party libraries in the Install Guide
   in doc directory.

#. Create the new distribution tarballs using the ``tarscript`` shell script
   under the ``scripts`` directory. This also compiles the documents (user
   guides and example docs) and creates all tarballs in their final form,
   appropriate for uploading to the website.

#. Update Internal Drupal Web pages for SUNDIALS:
   https://computation-external.llnl.gov/user

   a) Modify content (save edits on each page as you go)

      * Edit Main Page:
        https://computation-external.llnl.gov/projects/sundials

        * Update Release History page:
          https://computation-external.llnl.gov/projects/sundials/release-history
        * Add list of primary changes and child page containing
          full list of changes (this information is shared with
          "Changes in ..." sections of user guides).

      * Edit Download Page:
        https://computation-external.llnl.gov/projects/sundials/sundials-software

        * Update main download table with new versions of solvers and
          new sizes
        * Upload new documentation (pdf) files.
        * Update Previous releases table with new entry for previous
          release of full SUNDIALS suite.

      * Edit Usage Notes Example Page if major release:

        * Upload new ``cvRoberts_dns_negsol.c``, ``cvDisc_dns.c``, and
          ``cvsRoberts_FSA_dns_Switch.c`` examples if necessary.

      * Edit FAQ if necessary:
        https://computation-external.llnl.gov/projects/sundials/faq

   b) Once each sub page is complete, ask for team review of draft pages:
      https://computation-external.llnl.gov/projects/sundials

   c) After team comments are included and saved, select the
      "Publishing options" button in the bottom left group of buttons on the
      draft page. Select the Moderation state reflecting the amount of
      required review then Save. This must be done for each page and is the
      final action before pages are uploaded for external release.

#. After final push, ensure web content and behavior is as expected on the main
   page: http://computation.llnl.gov/projects/sundials

#. Tag the release

Note as of 28 Aug 2019 the addresses web address above are still vaild but may
change from "computation" to "computing" in the future.

**Old steps for maintaianed code:**

#. Create PDF files for SundialsTB:

   a) Create the PDF doc for SundialsTB by running the Matlab program
      ``texdoc.m`` available in ``sundialsTB/doc``.

   b) The program uses the m2html toolbox, freely available. It creates doc
      files in PS and PDF formats as ``sundialsTB.ps`` and ``sundialsTB.pdf``.

   c) Follow Radu's instructions in ``sundials/sundialsTB/doc/README_texdoc``.
