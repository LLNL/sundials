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

ERKStep Full RHS
================

Dotted lines signify execution paths that should not be possible.

.. digraph:: erk_fullrhs

   node [shape=box]

   start [label="Initialize fn_current to False"];
   h0    [label="Is h0 provided?", target="_top"];
   f0q1  [label="Has f(t0, y0) been computed?"];
   eval1 [label="Evaluate f(t0, y0)\lStore in F[0]\lCopy F[0] to fn\lSet fn_current to True"];
   step1 [label="Start Step 1"];

   start -> h0
   h0 -> step1 [label="Yes"]
   h0 -> f0q1 [label="No"]
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
   lagrange1     [label="Add ycur to history"]
   f0q3          [label="Has f(t0, y0) been computed?"];
   update_yn     [label="Copy tcur to tn\lCopy ycur to yn\lSet fn_current to False"];

   complete -> interp_update
   interp_update -> update_yn [label="No"]
   interp_update -> lagrange1 [label="Lagrange"]
   lagrange1 -> update_yn
   interp_update -> f0q3 [label="Hermite"]
   f0q3 -> update_yn [label="Yes"]

   interp_eval [label="Evaluate interpolant?"]
   step2       [label="Start Step 2"];
   f1q1        [label="Has f(t1, y1) been computed?"];
   fsal1       [label="Is the method FSAL?"]
   eval3       [label="Evaluate f(t1, y1)\lStore in F[0]"];
   copy1       [label="Copy F[s-1] to F[0]"];
   fcur1       [label="Copy F[0] to fn\lSet fn_current to True"];
   interp_yout [label="Compute yout"]

   update_yn -> interp_eval
   interp_eval -> step2 [label="No"]
   interp_eval -> interp_yout [label="Lagrange"]
   interp_eval -> f1q1 [label="Hermite"]
   f1q1 -> fsal1 [label="No"]
   fsal1 -> eval3 [label="No"]
   eval3 -> fcur1
   fsal1 -> copy1 [label="Yes"]
   copy1 -> fcur1
   fcur1 -> interp_yout
   interp_yout -> step2

   f1q2    [label="Has f(t1, y1) been computed?"];
   fsal2   [label="Is the method FSAL?"]
   eval4   [label="Evaluate f(t1, y1)\lStore in F[0]"];
   copy2   [label="Copy F[s-1] to F[0]"];
   fcur2   [label="Copy F[0] to fn\lSet fn_current to True"];
   stages2 [label="Compute Stages"]

   step2 -> f1q2
   f1q2 -> stages2 [label="Yes"]
   f1q2 -> fsal2 [label="No"]
   fsal2 -> copy2 [label="Yes"]
   copy2 -> fcur2
   fsal2 -> eval4 [label="No"]
   eval4 -> fcur2
   fcur2 -> stages2
   stages2 -> "Complete Step"


ARKStep Full RHS
================

Dotted lines signify execution paths that should not be possible.

.. digraph:: ark_fullrhs

   node [shape=box]

   start  [label="Initialize fn_current to False"];
   h0     [label="Is h0 provided?", target="_top"];
   f0q1   [label="Has f(t0, y0) been computed?"];
   eval1  [label="Evaluate fe(t0, y0), fi(t0,y0)\lStore in Fe[0], Fi[0]"];
   mass1a [label="Is there a mass matrix?"];
   mass1b [label="Solve M(t) u = Fe[0], M(t) v = Fi[0]\lStore u, v in Fe[0], Fi[0]\lCopy Fe[0] + Fi[0] to fn"];
   mass1c [label="Copy Fe[0] + Fi[0] to fn\lSolve M x = fn\lCopy x to fn"];
   fcur1  [label="Set fn_current to True"];
   step1  [label="Start Step 1"];

   start -> h0
   h0 -> step1 [label="Yes"]
   h0 -> f0q1 [label="No"]
   f0q1 -> step1 [label="Yes", style="dotted"]
   f0q1 -> eval1 [label="No"]
   eval1 -> mass1a
   mass1a -> mass1b [label="Yes\lM(t)"]
   mass1a -> mass1c [label="Yes\lM"]
   mass1a -> fcur1 [label="No"]
   mass1b -> fcur1
   mass1c -> fcur1
   fcur1 -> step1

   f0q2a    [label="Is the first stage explicit?\nor\nIs the method stiffly accurate and Hermite interpolation is used?"];
   f0q2b    [label="Has f(t0, y0) been computed?"];
   eval2    [label="Evaluate fe(t0, y0), fi(t0,y0)\lStore in Fe[0], Fi[0]"];
   mass2a   [label="Is there a mass matrix?"];
   mass2b   [label="Solve M(t) u = Fe[0], M(t) v = Fi[0]\lStore u, v in Fe[0], Fi[0]\lCopy Fe[0] + Fi[0] to fn"];
   mass2c   [label="Copy Fe[0] + Fi[0] to fn\lSolve M x = fn\lCopy x to fn"];
   fcur2    [label="Set fn_current to True"];
   stages1  [label="Compute Stages"]
   complete [label="Complete Step"]

   step1 -> f0q2a
   f0q2a -> f0q2b [label="Yes"]
   f0q2a -> stages1 [label="No"]
   f0q2b -> stages1 [label="Yes"]
   f0q2b -> eval2 [label="No"]
   eval2 -> mass2a
   mass2a -> mass2b [label="Yes\lM(t)"]
   mass2a -> mass2c [label="Yes\lM"]
   mass2a -> fcur2 [label="No"]
   mass2b -> fcur2
   mass2c -> fcur2
   fcur2 -> stages1
   stages1 -> complete

   interp_update [label="Update interpolation module?"];
   hermite1      [label="Hermite interpolation?"];
   lagrange1     [label="Add ycur to history"]
   f0q3          [label="Has f(t0, y0) been computed?"];
   eval3         [label="Evaluate fe(t0, y0), fi(t0,y0)\lStore in Fe[0], Fi[0]"];
   mass3a        [label="Is there a mass matrix?"];
   mass3b        [label="Solve M(t) u = Fe[0], M(t) v = Fi[0]\lStore u, v in Fe[0], Fi[0]\lCopy Fe[0] + Fi[0] to fn"];
   mass3c        [label="Copy Fe[0] + Fi[0] to fn\lSolve M x = fn\lCopy x to fn"];
   fcur3         [label="Set fn_current to True"];
   update_yn     [label="Copy tcur to tn\lCopy ycur to yn\lSet fn_current to False"];

   complete -> interp_update
   interp_update -> update_yn [label="No"]
   interp_update -> hermite1 [label="Yes"]
   hermite1 -> lagrange1 [label="No"]
   hermite1 -> f0q3 [label="Yes"]
   lagrange1 -> update_yn
   f0q3 -> update_yn [label="Yes"]
   f0q3 -> eval3 [label="No"]
   eval3 -> mass3a
   mass3a -> mass3b [label="Yes\lM(t)"]
   mass3a -> mass3c [label="Yes\lM"]
   mass3a -> fcur3 [label="No"]
   mass3b -> fcur3
   mass3c -> fcur3
   fcur3 -> update_yn

   interp      [label="Interpolate output?"];
   interp_eval [label="Evaluate interpolant"]
   step2       [label="Start Step 2"];
   f1q1        [label="Has f(t1, y1) been computed?"];
   hermite2    [label="Hermite interpolation?"];
   sa          [label="Is the method stiffly accurate?"]
   eval4       [label="Evaluate fe(t1, y1), fi(t1, y1)\lStore in Fe[0], Fi[0]"];
   mass4a      [label="Is there a mass matrix?"];
   mass4b      [label="Solve M(t) u = Fe[0], M(t) v = Fi[0]\lStore u, v in Fe[0], Fi[0]"];
   mass4c      [label="Solve M x = fn\lCopy x to fn"];
   copy1       [label="Copy Fe[s-1], Fi[s-1] to Fe[0], Fi[0]"];
   copy2       [label="Copy Fe[0] + Fi[0] to fn"];
   fcur4       [label="Set fn_current to True"];
   interp_yout [label="Compute yout"]

   update_yn -> interp
   interp -> step2 [label="No"]
   interp -> interp_eval [label="Yes"]
   interp_eval -> hermite2
   hermite2 -> interp_yout [label="No"]
   hermite2 -> f1q1 [label="Yes"]
   f1q1 -> interp_yout [label="Yes", style="dotted"]
   f1q1 -> sa [label="No"]
   sa -> copy1 [label="Yes"]
   sa -> eval4 [label="No"]
   copy1 -> copy2
   eval4 -> copy2
   copy2 -> mass4a
   mass4a -> mass4b [label="Yes\lM(t)"]
   mass4a -> mass4c [label="Yes\lM"]
   mass4a -> fcur4 [label="No"]
   mass4c -> fcur4
   fcur4 -> interp_yout
   interp_yout -> step2

   f1q2    [label="Has f(t1, y1) been computed?"];
   fsal    [label="Is the method FSAL?"]
   eval5   [label="Evaluate f(t1, y1)\lStore in F[0]\lCopy F[0] to fn\lSet fn_current to True"];
   copy3   [label="Copy F[s-1] to F[0]\lCopy F[0] to fn\lSet fn_current to True"];
   stages2 [label="Compute Stages"]

   step2 -> f1q2
   f1q2 -> stages2 [label="Yes"]
   f1q2 -> fsal [label="No"]
   fsal -> copy3 [label="Yes"]
   fsal -> eval5 [label="No"]
   copy3 -> stages2
   eval5 -> stages2
   stages2 -> "Complete Step"
