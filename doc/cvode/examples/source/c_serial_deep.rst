..
   Programmer(s): Daniel M. Margolis @ SMU
   ----------------------------------------------------------------
   SUNDIALS Copyright Start
   Copyright (c) 2002-2023, Lawrence Livermore National Security
   and Southern Methodist University.
   All rights reserved.

   See the top-level LICENSE and NOTICE files for details.

   SPDX-License-Identifier: BSD-3-Clause
   SUNDIALS Copyright End
   ----------------------------------------------------------------

:tocdepth: 3


.. _serial_deep_c:

============================================
Serial C example problems -- Deep dive
============================================



.. _deep_dive.cvRoberts_dns:

A dense example: cvRoberts_dns [DD]
============================================

As an initial illustration of the use of the CVODE package for the
integration of IVP ODEs, we give a sample program called :literal:`cvRoberts_dns.c`.
It uses the CVODE linear solver interface CVode with
dense matrix and linear solver modules (SUNMAT_DENSE and SUNLINSOL_DENSE)
and the NVECS module (which provides a serial implementation of NVECTOR)
in the solution of a 3-species chemical kinetics problem.

The problem consists of the following three rate equations:

.. math::

    \dot{y}_1 &= -0.04 \cdot y_1 + 10^4 \cdot y_2 \cdot y_3 \\
    \dot{y}_2 &=  0.04 \cdot y_1 - 10^4 \cdot y_2 \cdot y_3
                                  - 3 \cdot 10^7 \cdot y_2^2 \\
    \dot{y}_3 &=  3 \cdot 10^7 \cdot y_2^2

on the interval :math:`t \in [0, ~4 \cdot 10^{10}]`, with initial conditions
:math:`y_1(0) = 1.0`, :math:`y_2(0) = y_3(0) = 0.0`.
While integrating the system, we also use the root-finding
feature to find the points at which :math:`y_1 = 10^{-4}` or at which
:math:`y_3 = 0.01`.

For the source we give a rather detailed explanation of the parts of the program and
their interaction with CVODE.

Following the initial comment block, this program has a number
of ``#include`` lines, which allow access to useful items in CVODE
header files.  The :literal:`sundials_types.h` file provides the definition of the
type :literal:`realtype` (see `Types<https://sundials.readthedocs.io/en/latest/arkode/Usage/General.html#data-types>`_
for details).  For now, it suffices to read :literal:`realtype` as :literal:`double`.
The :literal:`cvode.h` file provides prototypes for the CVODE
functions to be called (excluding the linear solver selection
function), and also a number of constants that are to be used in
setting input arguments and testing the return value of CVode.
The :literal:`sunlinsol_dense.h`
file is the header file for the dense implementation of the
SUNLINSOL module and includes definitions of the
:literal:`SUNLinearSolver` type.  Similarly, the :literal:`sunmatrix_dense.h`
file is the header file for the dense implementation of the
SUNMATRIX module, including definitions of the :literal:`SUNMatrix` type
as well as macros and functions to access matrix components.  
We have explicitly included \id{sunmatrix\_dense.h}, but this is not
necessary because it is included by \id{sunlinsol\_dense.h}.
The \id{nvector\_serial.h} file is the header file for the serial
implementation of the {\nvector} module and includes definitions of the 
\id{N\_Vector} type, a macro to access vector components, and prototypes 
for the serial implementation specific machine environment memory allocation
and freeing functions.

This program includes two user-defined accessor macros, \id{Ith} and
\id{IJth}, that are useful in writing the problem functions in a form
closely matching the mathematical description of the ODE system,
i.e. with components numbered from 1 instead of from 0.  The \id{Ith}
macro is used to access components of a vector of type \id{N\_Vector}
with a serial implementation.  It is defined using the {\nvecs}
accessor macro \id{NV\_Ith\_S} which numbers components starting with
0. The \id{IJth} macro is used to access elements of a dense matrix of
type \id{SUNMatrix}.  It is similarly defined using the {\sunmatdense}
accessor macro \id{SM\_ELEMENT\_D} which numbers matrix rows and
columns starting with 
0. 
The macro \id{NV\_Ith\_S} is fully described in \ugref{ss:nvec_ser}.
The macro \id{SM\_ELEMENT\_D} is fully described in \ugref{ss:sunmat_dense}.

Next, the program includes some problem-specific constants, which are
isolated to this early location to make it easy to change them as needed.  
The program prologue ends with prototypes of four private helper
functions and the three user-supplied functions that are called by {\cvode}.

The \id{main} program begins with some dimensions and type declarations,
including use of the generic types \id{N\_Vector}, \id{SUNMatrix} and
\id{SUNLinearSolver}.  The next several lines 
allocate memory for the \id{y} and \id{abstol} vectors using
\id{N\_VNew\_Serial} with a length argument of \id{NEQ} ($= 3$). The
lines following that load the initial values of the dependendent
variable vector into \id{y} and the absolute tolerances into \id{abstol}
using the \id{Ith} macro.

The calls to \id{N\_VNew\_Serial}, and also later calls to \id{CVode***}
functions, make use of a private function, \id{check\_flag}, which examines
the return value and prints a message if there was a failure.  The
\id{check\_flag} function was written to be used for any serial {\sundials}
application.

The call to \id{CVodeCreate} creates the {\cvode} solver memory block,
specifying the \id{CV\_BDF} integration method with \id{CV\_NEWTON} iteration.
Its return value is a pointer to that memory block for this
problem.  In the case of failure, the return value is \id{NULL}.  This
pointer must be passed in the remaining calls to {\cvode} functions.

The call to \id{CVodeInit} allocates and initializes the solver memory block.
Its arguments include the name of the {\CC} function \id{f} defining the
right-hand side function $f(t,y)$, and the initial values of $t$ and $y$.
The call to \id{CVodeSVtolerances} specifies a vector of absolute tolerances,
and includes the value of the relative tolerance \id{reltol} and the absolute 
tolerance vector \id{abstol}.  See \ugref{sss:cvodemalloc} and
\ugref{sss:cvtolerances} for full details of these calls.

The call to \id{CVodeRootInit} specifies that a rootfinding problem
is to be solved along with the integration of the ODE system, that the
root functions are specified in the function \id{g}, and that there are
two such functions.  Specifically, they are set to $y_1 - 0.0001$ and 
$y_3 - 0.01$, respectively.
See \ugref{ss:cvrootinit} for a detailed description of this call.

The call to \id{SUNDenseMatrix} (see \ugref{ss:sunmat_dense}) creates
a \id{NEQ}$\times$\id{NEQ} dense {\sunmatrix} object to use within the
Newton solve in {\cvode}.  The following call to
\id{SUNLinSol\_Dense} (see \ugref{ss:sunlinsol_dense}) creates the 
dense {\sunlinsol} object that will perform the linear solves within
the Newton method.  These are attached to the {\cvls} linear
solver interface with the  call to \id{CVodeSetLinearSolver} (see
\ugref{sss:lin_solv_init}), and the subsequent call to
\id{CVodeSetJacFn} (see \ugref{ss:optional_input}) specifies the
analytic Jacobian supplied by the user-supplied function \id{Jac}.

The actual solution of the ODE initial value problem is accomplished in
the loop over values of the output time \id{tout}.  In each pass of the
loop, the program calls \id{CVode} in the \id{CV\_NORMAL} mode, meaning that
the integrator is to take steps until it overshoots \id{tout} and then
interpolate to $t = $\id{tout}, putting the computed value of $y$(\id{tout})
into \id{y}, with \id{t} = \id{tout}.  The return value in this case is
\id{CV\_SUCCESS}.  However, if \id{CVode} finds a root before reaching the next
value of \id{tout}, it returns \id{CV\_ROOT\_RETURN} and stores the root
location in \id{t} and the solution there in \id{y}.  In either case, the
program prints \id{t} and \id{y}.  In the case of a root, it calls
\id{CVodeGetRootInfo} to get a length-2 array \id{rootsfound} of bits showing
which root function was found to have a root.  If \id{CVode} returned any
negative value (indicating a failure), the program breaks out of the loop.  
In the case of a \id{CV\_SUCCESS} return, the value of \id{tout} is
advanced (multiplied by 10) and a counter (\id{iout}) is advanced, so
that the loop can be ended when that counter reaches the preset number
of output times, \id{NOUT} = 12.  See \ugref{sss:cvode} for full
details of the call to \id{CVode}.

Finally, the main program calls \id{PrintFinalStats} to get and print
all of the relevant statistical quantities.  It then calls \id{NV\_Destroy}
to free the vectors \id{y} and \id{abstol}, \id{CVodeFree} to free the 
{\cvode} memory block, \id{SUNLinSolFree} to free the linear solver
memory, and \id{SUNMatDestroy} to free the matrix \id{A}.

The function \id{PrintFinalStats} used here is actually suitable for
general use in applications of {\cvode} to any problem with a direct
linear solver.  It calls various \id{CVodeGet***} 
functions to obtain the relevant counters, and then prints them.
Specifically, these are: the cumulative number of steps (\id{nst}), 
the number of \id{f} evaluations (\id{nfe}) (excluding those for
difference-quotient Jacobian evaluations),
the number of matrix factorizations (\id{nsetups}),
the number of \id{f} evaluations for Jacobian evaluations (\id{nfeLS}
= 0 here),
the number of Jacobian evaluations (\id{nje}),
the number of nonlinear (Newton) iterations (\id{nni}),
the number of nonlinear convergence failures (\id{ncfn}),
the number of local error test failures (\id{netf}), and
the number of \id{g} (root function) evaluations (\id{nge}).
These optional outputs are described in \ugref{ss:optional_output}.

The function \id{f} is a straightforward expression of the ODEs. 
It uses the user-defined macro \id{Ith} to extract the components of \id{y}
and to load the components of \id{ydot}.
See \ugref{ss:rhsFn} for a detailed specification of \id{f}.

Similarly, the function \id{g} defines the two functions, $g_0$ and $g_1$,
whose roots are to be found.  See \ugref{ss:rootFn} for a detailed description
of the \id{g} function.

The function \id{Jac} sets the nonzero elements of the Jacobian as a
dense matrix.  (Zero elements need not be set because \id{J} is preset
to zero.)  It uses the user-defined macro \id{IJth} to reference the
elements of a dense matrix of type {\sunmatrix}.  Here the problem
size is small, so we need not worry about the inefficiency of using
\id{NV\_Ith\_S} and \id{SM\_ELEMENT\_D} to access \id{N\_Vector} and
{\sunmatdense} elements.  Note that in this example, \id{Jac}
only accesses the \id{y} and \id{J} arguments.  See \ugref{ss:jacFn}
for a detailed description of the \id{Jac} function.

The output generated by \id{cvRoberts\_dns} is shown below.  It shows the output
values at the 12 preset values of \id{tout}.  It also shows the two root
locations found, first at a root of $g_1$, and then at a root of $g_0$.


.. _deep_dive.cvAdvDiff_bnd:

A banded example: cvAdvDiff_bnd [DD]
============================================




.. _deep_dive.cvDiurnal_kry:

A Krylov example: cvDiurnal_kry [DD]
============================================



