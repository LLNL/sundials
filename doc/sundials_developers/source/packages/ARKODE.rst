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

.. _ARKODE.ERKFullRHS:

ERKStep Full RHS Flow
=====================

Dotted lines signify execution paths that should not be possible.

.. digraph:: fullrhs

   node [shape=box]

   start [label="Initialize fn_current to False"];
   h0    [label="Is h0 provided?", target="_top"];
   f0q1  [label="Has f(t0, y0) been computed?"];
   eval1 [label="Evaluate f(t0, y0)\lStore in F[0]\lCopy F[0] to fn\lSet fn_current to True"];
   step1 [label="Start Step 1"];

   start -> h0
   h0 -> step1 [label="Yes"]
   h0 -> f0q1 [label="No"]
   f0q1 -> step1 [label="Yes", style="dotted"]
   f0q1 -> eval1 [label="No"]
   eval1 -> step1

   f0q2     [label="Has f(t0, y0) been computed?"];
   eval2    [label="Evaluate f(t0, y0)\lStore in F[0]\lCopy F[0] to fn\lSet fn_current to True"];
   stages1  [label="Compute Stages"]
   complete [label="Complete Step"]

   step1 -> f0q2
   f0q2 -> stages1 [label="Yes"]
   f0q2 -> eval2 [label="No"]
   eval2 -> stages1
   stages1 -> complete

   interp_update [label="Update interpolation module?"];
   hermite1      [label="Hermite interpolation?"];
   lagrange1     [label="Add ycur to history"]
   f0q3          [label="Has f(t0, y0) been computed?"];
   eval3         [label="Evaluate f(t0, y0)\lStore in F[0]\lCopy F[0] to fn\lSet fn_current to True"];
   update_yn     [label="Copy tcur to tn\lCopy ycur to yn\lSet fn_current to False"];

   complete -> interp_update
   interp_update -> update_yn [label="No"]
   interp_update -> hermite1 [label="Yes"]
   hermite1 -> lagrange1 [label="No"]
   hermite1 -> f0q3 [label="Yes"]
   lagrange1 -> update_yn
   f0q3 -> update_yn [label="Yes"]
   f0q3 -> eval3 [label="No", style="dotted"]
   eval3 -> update_yn [style="dotted"]

   interp      [label="Interpolate output?"];
   interp_eval [label="Evaluate interpolant"]
   step2       [label="Start Step 2"];
   f1q1        [label="Has f(t1, y1) been computed?"];
   hermite2    [label="Hermite interpolation?"];
   fsal1       [label="Is the method FSAL?"]
   eval4       [label="Evaluate f(t1, y1)\lStore in F[0]\lCopy F[0] to fn\lSet fn_current to True"];
   copy1       [label="Copy F[s-1] to F[0]\lCopy F[0] to fn\lSet fn_current to True"];
   interp_yout [label="Compute yout"]

   update_yn -> interp
   interp -> interp_yout [label="No"]
   interp -> interp_eval [label="Yes"]
   interp_eval -> hermite2
   hermite2 -> interp_yout [label="No"]
   hermite2 -> f1q1 [label="Yes"]
   f1q1 -> interp_yout [label="Yes", style="dotted"]
   f1q1 -> fsal1 [label="No"]
   fsal1 -> eval4 [label="No"]
   fsal1 -> copy1 [label="Yes"]
   copy1 -> interp_yout
   eval4 -> interp_yout
   interp_yout -> step2

   f1q2    [label="Has f(t1, y1) been computed?"];
   fsal2   [label="Is the method FSAL?"]
   eval5   [label="Evaluate f(t1, y1)\lStore in F[0]\lCopy F[0] to fn\lSet fn_current to True"];
   copy2   [label="Copy F[s-1] to F[0]\lCopy F[0] to fn\lSet fn_current to True"];
   stages2 [label="Compute Stages"]

   step2 -> f1q2
   f1q2 -> stages2 [label="Yes"]
   f1q2 -> fsal2 [label="No"]
   fsal2 -> copy2 [label="Yes"]
   fsal2 -> eval5 [label="No"]
   copy2 -> stages2
   eval5 -> stages2
   stages2 -> "Complete Step"
