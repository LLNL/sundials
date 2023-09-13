.. ----------------------------------------------------------------
   SUNDIALS Copyright Start
   Copyright (c) 2002-2023, Lawrence Livermore National Security
   and Southern Methodist University.
   All rights reserved.

   See the top-level LICENSE and NOTICE files for details.

   SPDX-License-Identifier: BSD-3-Clause
   SUNDIALS Copyright End
   ----------------------------------------------------------------

.. _SUNHeuristics:

#####################################
Heuristic Time Step Constraints
#####################################

The SUNDIALS library comes packaged with a variety of :c:type:`SUNHeuristics`
implementations, designed to support various forms of heuristic time step
constraints within SUNDIALS time integrators.  To support applications that may
want to adjust or disable these heuristic controls, SUNDIALS provides a
:c:type:`SUNHeuristics` base class, along with two default implementations: a
"default" heuristics control module that encodes standard SUNDIALS heuristic
controls, and an "unconstrained" module that disables these.

.. toctree::
   :maxdepth: 1

   SUNHeuristics_links.rst
