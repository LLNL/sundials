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

.. _OpenPR:

Opening a Pull Request
======================

When a branch is ready to be reviewed for integration into the ``master`` or
``develop`` branches follow the steps below to open a pull request:

#. Browse to `https://github.com/LLNL/sundials <https://github.com/LLNL/sundials>`_

#. Click on the branch icon on the left side of screen - you'll see a list
   of all available branches

#. Click on your branch - you'll see a 'Compare' screen that lets you pick a
   branch (source on top) to merge with another branch (target on bottom)

#. Select the desired branches and click 'Create pull request'

#. Edit the title of the pull request (defaults to branch name), add a
   description, and select reviewers that can approve the request

#. Click 'Create'

The selected reviewers will go over the changes in the pull request and may
ask for additional changes before merging the branch. After the pull request is
merged, delete the local copy the branch:

.. code-block:: none

   $ git checkout PARENT
   $ git branch -D <branch-name>
