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

Description -- step-by-step
----------------------------

As an initial illustration of the use of the **CVODE** package for the
integration of IVP ODEs, we give a sample program called :literal:`cvRoberts_dns.c`.
It uses the **CVODE** *linear solver interface* **CVode** with
dense matrix and *linear solver* modules (**SUNMATRIX_DENSE** and **SUNLinSol_Dense**)
and the **NVECTOR_SERIAL** module (which provides a *serial* implementation of **NVECTOR**)
in the solution of a :math:`3`-species chemical kinetics problem.

The problem consists of the following three rate equations:

.. math::

   \dot{y}_1 &= -0.04 \cdot y_1 + 10^4 \cdot y_2 \cdot y_3 \\
   \dot{y}_2 &=  0.04 \cdot y_1 - 10^4 \cdot y_2 \cdot y_3
                                - 3 \cdot 10^7 \cdot y_2^2 \\
   \dot{y}_3 &=  3 \cdot 10^7 \cdot y_2^2

on the interval :math:`t \in [0, ~4 \cdot 10^{10}]`, with initial conditions
:math:`y_1(0) = 1.0`, :math:`y_2(0) = y_3(0) = 0.0`.
While integrating the system, we also use the *root-finding
feature* to find the points at which :math:`y_1 = 10^{-4}` or at which
:math:`y_3 = 0.01`.

For the source we give a rather detailed explanation of the parts of the program and
their interaction with **CVODE**.

Following the initial comment block, this program has a number
of ``#include`` lines, which allow access to useful items in **CVODE**
header files.  The :literal:`sundials_types.h` file provides the definition of the
type :literal:`realtype` (see `Types
<https://sundials.readthedocs.io/en/latest/arkode/Usage/General.html#data-types>`_
for details).  For now, it suffices to read :literal:`realtype` as :literal:`double`.
The :literal:`cvode.h` file provides prototypes for the **CVODE**
functions to be called (excluding the *linear solver* selection
function), and also a number of constants that are to be used in
setting input arguments and testing the return value of **CVode**.
The :literal:`sunlinsol_dense.h`
file is the header file for the *dense* implementation of the
**SUNLINSOL** module and includes definitions of the
:literal:`SUNLinearSolver` type.  Similarly, the :literal:`sunmatrix_dense.h`
file is the header file for the *dense* implementation of the
**SUNMATRIX** module, including definitions of the :literal:`SUNMatrix` type
as well as macros and functions to access *matrix components*.  
We have explicitly included :literal:`sunmatrix_dense.h`, but this is not
necessary because it is included by :literal:`sunlinsol_dense.h`.
The :literal:`nvector_serial.h` file is the header file for the *serial*
implementation of the **NVECTOR** module and includes definitions of the 
:literal:`N_Vector` type, a macro to access vector components, and prototypes 
for the *serial* implementation specific machine environment memory allocation
and freeing functions.

This program includes two *user-defined accessor macros*, :literal:`Ith` and
:literal:`IJth`, that are useful in writing the problem functions in a form
closely matching the mathematical description of the ODE system,
i.e. with components numbered from :math:`1` instead of from :math:`0`.  The :literal:`Ith`
*macro* is used to access components of a vector of type :literal:`N_Vector`
with a *serial* implementation.  It is defined using the **NVECTOR_SERIAL**
*accessor macro* :literal:`NV_Ith_S` which numbers components starting with
:math:`0`. The :literal:`IJth` *macro* is used to access elements of a dense matrix of
type :literal:`SUNMatrix`.  It is similarly defined using the **SUNMATRIX_DENSE**
*accessor macro* :literal:`SM_ELEMENT_D` which numbers matrix rows and
columns starting with :math:`0`.

The *macro* :literal:`NV_Ith_S` is fully described in `NVECTOR_SERIAL
<https://sundials.readthedocs.io/en/latest/nvectors/NVector_links.html>`_.
The *macro* :literal:`SM_ELEMENT_D` is fully described in `SUNMATRIX_DENSE
<https://sundials.readthedocs.io/en/latest/sunmatrix/SUNMatrix_links.html>`_.

Next, the program includes some problem-specific constants, which are
isolated to this early location to make it easy to change them as needed.  
The program prologue ends with prototypes of four private helper
functions and the three *user-supplied functions* that are called by **CVODE**.

The :literal:`main` program begins with some *dimensions* and *type* declarations,
including use of the generic types :literal:`N_Vector`, :literal:`SUNMatrix` and
:literal:`SUNLinearSolver`.  The next several lines 
allocate memory for the :math:`y` and :math:`abstol` vectors using
:literal:`N_VNew_Serial` with a length argument of :math:`NEQ (= 3)`. The
lines following that load the *initial values* of the dependendent
variable vector into :math:`y` and the *absolute tolerances* into :math:`abstol`
using the :literal:`Ith` macro.

The calls to :literal:`N_VNew_Serial`, and also later calls to :literal:`CVode***`
functions, make use of a private function, :literal:`check_flag`, which examines
the return value and prints a message if there was a failure.  The
:literal:`check_flag` function was written to be used for any *serial* **SUNDIALS**
application.

The call to :literal:`CVodeCreate` creates the **CVODE** solver memory block,
specifying the :literal:`CV_BDF` integration method with built-in ``CV_NEWTON`` iteration.
Its return value is a pointer to that memory block for this
problem.  In the case of failure, the return value is ``NULL``.  This
pointer must be passed in the remaining calls to **CVODE** functions.

The call to :literal:`CVodeInit` *allocates* and *initializes* the solver memory block.
Its arguments include the name of the :literal:`C` function :math:`f` defining the
right-hand side function :math:`f(t,y)`, and the *initial values* of :math:`t` and :math:`y`.
The call to :literal:`CVodeSVtolerances` specifies a *vector* of *absolute tolerances*,
and includes the *value* of the *relative tolerance* :math:`reltol` and the *absolute 
tolerance vector* :math:`abstol`.  See
`CVodeCreate and CVodeInit
<https://sundials.readthedocs.io/en/latest/cvode/Usage/index.html#cvode-initialization-and-deallocation-functions>`_
and `CVodeSStolerances and CVodeSVtolerances
<https://sundials.readthedocs.io/en/latest/cvode/Usage/index.html#cvode-tolerance-specification-functions>`_
for full details of these calls.

The call to :literal:`CVodeRootInit` specifies that a *root-finding* problem
is to be solved along with the integration of the ODE system, that the
root functions are specified in the function :math:`g`, and that there are
two such functions.  Specifically, they are set to :math:`y_1 - 0.0001` and 
:math:`y_3 - 0.01`, respectively. See `CVodeRootInit
<https://sundials.readthedocs.io/en/latest/cvode/Usage/index.html#rootfinding-initialization-function>`_
for a detailed description of this call.

The call to :literal:`SUNDenseMatrix` (see `SUNMATRIX_DENSE
<https://sundials.readthedocs.io/en/latest/sunmatrix/SUNMatrix_links.html#>`_)
creates a :math:`NEQ \times NEQ` *dense* **SUNMATRIX** object to use within the
*Newton* solve in **CVODE**.  The following call to :literal:`SUNLinSol_Dense` (see `SUNLinSol_Dense 
<https://sundials.readthedocs.io/en/latest/sunlinsol/SUNLinSol_links.html#the-sunlinsol-dense-module>`_)
creates the *dense* **SUNLINSOL** object that will perform the linear solves within
the *Newton* method.  These are attached to the **CVODE** *linear
solver interface* with the call to :literal:`CVodeSetLinearSolver` (see `CVodeSetLinearSolver interface 
<https://sundials.readthedocs.io/en/latest/cvode/Usage/index.html#linear-solver-interface-functions>`_),
and the subsequent call to :literal:`CVodeSetJacFn` (see `CVode Linear Solver Optional Input Functions 
<https://sundials.readthedocs.io/en/latest/cvode/Usage/index.html#linear-solver-interface-optional-input-functions>`_)
specifies the analytic *Jacobian* supplied by the *user-supplied function* :math:`Jac`.

The actual solution of the ODE initial value problem is accomplished in
the loop over values of the output time :math:`t_{out}`.  In each pass of the
loop, the program calls :literal:`CVode` in the ``CV_NORMAL`` mode, meaning that
the integrator is to take steps until it overshoots :math:`t_{out}` and then
interpolate to :math:`t = t_{out}`, putting the computed value of :math:`y(t_{out})`
into :math:`y`, with :math:`t = t_{out}`.  The return value in this case is
``CV_SUCCESS``.  However, if :literal:`CVode` finds a root before reaching the next
value of :math:`t_{out}`, it returns ``CV_ROOT_RETURN`` and stores the root
location in :math:`t` and the solution at that point in :math:`y`.  In either case, the
program prints :math:`t` and :math:`y`.  In the case of a root, it calls
:literal:`CVodeGetRootInfo` to get a length-:math:`2` array :math:`rootsfound` of bits showing
which root function was found to have a root.  If :literal:`CVode` returned any
negative value (*indicating a failure*), the program breaks out of the loop.  
In the case of a ``CV_SUCCESS`` return, the value of :math:`t_{out}` is
advanced (multiplied by :math:`10`) and a counter (:math:`i_{out}`) is advanced, so
that the loop can be ended when that counter reaches the preset number
of output times, :math:`N_{OUT} = 12`.  See `CVode Solver Function 
<https://sundials.readthedocs.io/en/latest/cvode/Usage/index.html#cvode-solver-function>`_
for full details of the call to :literal:`CVode`.

Finally, the main program calls :literal:`PrintFinalStats` to get and print
all of the relevant statistical quantities.  It then calls :literal:`NV_Destroy`
to free the vectors :math:`y` and :math:`abstol`, :literal:`CVodeFree` to free the 
**CVODE** *memory block*, :literal:`SUNLinSolFree` to free the *linear solver
memory*, and :literal:`SUNMatDestroy` to *free the matrix* :math:`A`.

The function :literal:`PrintFinalStats` used here is actually suitable for
general use in applications of **CVODE** to any problem with a direct
*linear solver*.  It calls various :literal:`CVodeGet***` 
functions to obtain the relevant counters, and then prints them.
Specifically, these are:

* the cumulative number of steps (``nst``), 
* the number of :math:`f` evaluations (``nfe``) (excluding those for *difference-quotient Jacobian* evaluations),
* the number of matrix factorizations (``nsetups``),
* the number of :math:`f` evaluations for *Jacobian* evaluations (``nfeLS`` :math:`= 0` here),
* the number of *Jacobian* evaluations (``nje``),
* the number of nonlinear (*Newton*) iterations (``nni``),
* the number of nonlinear convergence failures (``ncfn``),
* the number of local error test failures (``netf``), and
* the number of :math:`g` (*root function*) evaluations (``nge``).

These optional outputs are described in `CVode Optional Output Functions 
<https://sundials.readthedocs.io/en/latest/cvode/Usage/index.html#optional-output-functions>`_ .

The function :literal:`f` is a straightforward expression of the ODEs. 
It uses the *user-defined macro* :literal:`Ith` to extract the components of :math:`y`
and to load the components of :math:`\dot y =` ``ydot``.  See `ODE RHS Function
<https://sundials.readthedocs.io/en/latest/cvode/Usage/index.html#ode-right-hand-side>`_
for a detailed specification of :math:`f`.

Similarly, the function :literal:`g` defines the two functions, :math:`g_0` and :math:`g_1`,
whose roots are to be found.  See `Rootfinding Function 
<https://sundials.readthedocs.io/en/latest/cvode/Usage/index.html#rootfinding-function>`_
for a detailed description of the :math:`g` function.

The function :literal:`Jac` sets the nonzero elements of the Jacobian as a
dense matrix.  (Zero elements need not be set because :math:`J` is preset
to zero.)  It uses the *user-defined macro* :literal:`IJth` to reference the
elements of a dense matrix of type **SUNMATRIX**.  Here the problem
size is small, so we need not worry about the inefficiency of using
:literal:`NV_Ith_S` and :literal:`SM_ELEMENT_D` to access :literal:`N_Vector` and
**SUNMATRIX_DENSE** elements.  Note that in this example, :literal:`Jac`
only accesses the :math:`y` and :math:`J` arguments.  See `Jacobian Construction 
<https://sundials.readthedocs.io/en/latest/cvode/Usage/index.html#jacobian-construction-matrix-based-linear-solvers>`_
for a detailed description of the :math:`Jac` function.

Program output
---------------

The output generated by :literal:`cvRoberts_dns` is shown again below.  It shows the output
values at the :math:`12` preset values of :math:`t_{out}`.  It also shows the two root
locations found, first at a root of :math:`g_1`, and then at a root of :math:`g_0`.

.. include:: ../../../../examples/cvode/serial/cvRoberts_dns.out
   :literal:


.. _deep_dive.cvAdvDiff_bnd:

A banded example: cvAdvDiff_bnd [DD]
============================================

Description -- step-by-step
----------------------------

The example program :literal:`cvAdvDiff_bnd.c` solves the semi-discretized form of
the :math:`2`-D advection-diffusion equation

.. math::
   
   \frac{\partial v}{\partial t} = \frac{\partial^2 v}{\partial x^2}
                                 + 0.5 \frac{\partial v}{\partial x}
                                 + \frac{\partial^2 v}{\partial y^2}

on a rectangle, with zero Dirichlet boundary conditions. The PDE is 
discretized with standard central finite differences on a 
:math:`(M_X + 2) \times (M_Y + 2)` mesh, giving an ODE system of size
:math:`M_X * M_Y`.  The discrete value :math:`v_{ij}` approximates :math:`v` at
:math:`x = i \Delta x`, :math:`y = j \Delta y`. The ODEs are

.. math::
   
   \frac{dv_{ij}}{dt} = f_{ij} =
         \frac{v_{i-1,j} - 2 v_{ij} + v_{i+1,j}}{(\Delta x)^2}
       + 0.5  \frac{v_{i+1,j} - v_{i-1,j}}{2 \Delta x}
       + \frac{v_{i,j-1} - 2 v_{ij} + v_{i,j+1}}{(\Delta y)^2} \, ,

where :math:`1 \leq i \leq M_X` and :math:`1 \leq j \leq M_Y`.  The boundary
conditions are imposed by taking :math:`v_{ij} = 0` above if :math:`i = 0`
or :math:`M_X + 1`, or if :math:`j = 0` or :math:`M_Y + 1`. 
If we set :math:`u_{(j-1)+(i-1)*M_Y} = v_{ij}`, so that the ODE system is
:math:`\dot u = f(u)`, then the system *Jacobian* :math:`J = \frac{\partial f}{\partial u}` is
a *band matrix* with *upper* and *lower half-bandwidths* both equal to :math:`M_Y`.
In the example, we take :math:`M_X = 10` and :math:`M_Y = 5`.

The :literal:`cvAdvDiff_bnd.c` program includes files
:literal:`sunmatrix_band.h` and :literal:`sunlinsol_band.h` in order to use the
**SUNLinSol_Band** *linear solver*. The :literal:`sunmatrix_band.h` file
contains the definition of the *banded* **SUNMATRIX** type, and the
:literal:`SM_COLUMN_B` and :literal:`SM_COLUMN_ELEMENT_B` *macros* for
accessing *banded matrix elements* (see `SUNMATRIX_BAND 
<https://sundials.readthedocs.io/en/latest/sunmatrix/SUNMatrix_links.html#the-sunmatrix-band-module>`_).
The :literal:`sunlinsol_band.h` file contains the definition of the *banded*
**SUNLINSOL** type.  We note that have explicitly included
:literal:`sunmatrix_band.h`, but this is not necessary because it is
included by :literal:`sunlinsol_band.h`.  The file :literal:`nvector_serial.h`
is included for the definition of the *serial* :literal:`N_Vector` type. 

The *include* lines at the top of the file are followed by *definitions* of
*problem constants* which include the :math:`x` and :math:`y` mesh dimensions, :math:`M_X` and
:math:`M_Y`, the *number of equations* :math:`NEQ`, the *scalar absolute tolerance*
:math:`ATOL`, the *initial time* :math:`T_0`, and the *initial output time* :math:`T_1`.

Spatial discretization of the PDE naturally produces an ODE system in
which equations are numbered by mesh coordinates :math:`(i,j)`. The
*user-defined macro* :literal:`IJth` isolates the *translation* for the
mathematical *two-dimensional* index to the *one-dimensional*
:literal:`N_Vector` index and allows the user to write clean, readable code
to access components of the dependent variable.  The :literal:`NV_DATA_S`
macro returns the component array for a given :literal:`N_Vector`, and this
array is passed to :literal:`IJth` in order to do the actual :literal:`N_Vector`
access.

The type :literal:`UserData` is a pointer to a structure containing problem
data used in the :literal:`f`` and :literal:`Jac` functions.  This structure is
allocated and initialized at the beginning of :literal:`main`. The pointer
to it, called :literal:`data`, is passed to :literal:`CVodeSetUserData`, and as a
result it will be passed back to the :literal:`f` and :literal:`Jac` functions
each time they are called.  The use of the :literal:`data` pointer
eliminates the need for global program data.

The :literal:`main` program is straightforward.  The :literal:`CVodeCreate` call specifies
the ``CV_BDF`` method with a built-in ``CV_NEWTON`` iteration. Following the
:literal:`CVodeInit` call, the call to :literal:`CVodeSStolerances` indicates *scalar
relative and absolute tolerances*, and *values* :math:`reltol` and :math:`abstol` are passed.
The call to :literal:`SUNBandMatrix` (see `SUNMATRIX_BAND 
<https://sundials.readthedocs.io/en/latest/sunmatrix/SUNMatrix_links.html#the-sunmatrix-band-module>`_)
creates a *banded* **SUNMATRIX** *Jacobian* template, and specifies that both
*half-bandwidths* of the *Jacobian* are equal to :math:`M_Y`.  The calls to
:literal:`SUNBandLinearSolver` (see `SUNLinSol_Band 
<https://sundials.readthedocs.io/en/latest/sunlinsol/SUNLinSol_links.html>`_)
and :literal:`CVodeSetLinearSolver` (see `CVode Linear Solver Functions 
<https://sundials.readthedocs.io/en/latest/cvode/Usage/index.html#linear-solver-interface-functions>`_)
specifies the **SUNLinSol_Band** *linear solver* to the **CVode** *interface*.
The call to :literal:`CVodeSetJacFn` (see `CVode Linear Solver Optional Input Functions
<https://sundials.readthedocs.io/en/latest/cvode/Usage/index.html#linear-solver-interface-optional-input-functions>`_)
specifies that a *user-supplied Jacobian* function :literal:`Jac` is to be used.

The actual solution of the problem is performed by
the call to :literal:`CVode` within the loop over the output times :math:`t_{out}`.
The max-norm of the solution vector (from a call to :literal:`N_VMaxNorm`) and
the cumulative number of time steps (from a call to :literal:`CVodeGetNumSteps`) are
printed at each output time. Finally, the calls to :literal:`PrintFinalStats`,
:literal:`N_VDestroy`, and :literal:`CVodeFree` print statistics and free problem memory.

Following the :literal:`main` program in the :literal:`cvAdvDiff_bnd.c` file are definitions of
five functions: :literal:`f`, :literal:`Jac`, :literal:`SetIC`, :literal:`PrintHeader`, :literal:`PrintOutput`,
:literal:`PrintFinalStats`, and :literal:`check_flag`.   The last five functions are called
only from within the :literal:`cvAdvDiff_bnd.c` file.
The :literal:`SetIC` function sets the initial dependent variable vector;
:literal:`PrintHeader` prints the heading of the output page;
:literal:`PrintOutput` prints a line of solution output;
:literal:`PrintFinalStats` gets and prints statistics at the end of the run;
and :literal:`check_flag` aids in checking return values.  The statistics
printed include counters such as the total number of steps (``nst``), 
:math:`f` evaluations (excluding those for *Jaobian* evaluations) (``nfe``),
:math:`LU` decompositions (``nsetups``), :math:`f` evaluations for
*difference-quotient Jacobians* (:math:`nfeLS = 0` here),
*Jacobian* evaluations (``nje``), and nonlinear iterations (``nni``).
These optional outputs are described in `CVode Optional Output Functions
<https://sundials.readthedocs.io/en/latest/cvode/Usage/index.html#optional-output-functions>`_.
Note that :literal:`PrintFinalStats` is suitable for general use in applications
of **CVODE** to any problem with a direct *linear solver*.

.. role:: C(code)
   :language: C

The :literal:`f` function implements the *central difference approximation* (`IVP Solutions
Mathematics <https://sundials.readthedocs.io/en/latest/cvode/Mathematics_link.html>`_)
with :math:`u` identically zero on the boundary. 
The constant coefficients :math:`(\Delta x)^{-2}`, :math:`0.5 (2 \Delta x)^{-1}`, and 
:math:`(\Delta y)^{-2}` are computed only once at the beginning of :literal:`main`, 
and stored in the locations :C:`data->hdcoef`, :C:`data->hacoef`, and
:C:`data->vdcoef`, respectively.   When :literal:`f` receives the :literal:`data`
pointer (renamed :literal:`user_data` here), it pulls out these values from storage
in the local variables ``hordc``, ``horac``, and ``verdc``.  It then
uses these to construct the diffusion and advection terms, which are
combined to form :math:`\dot u =` ``udot``.  Note the extra lines setting out-of-bounds
values of :math:`u` to zero.

The :literal:`Jac` function is an expression of the derivatives

.. math::

   \frac{\partial f_{ij}}{\partial v_{ij}} &=&
         -2 \left[ (\Delta x)^{-2} + (\Delta y)^{-2} \right] \\
   \frac{\partial f_{ij}}{\partial v_{i \pm 1,j}} &=& (\Delta x)^{-2} 
                  \pm 0.5 (2 \Delta x)^{-1}, ~~~~
   \frac{\partial f_{ij}}{\partial v_{i,j \pm 1}}  =  (\Delta y)^{-2} ~.

This function loads the *Jacobian* by *columns*, and like :literal:`f` it
makes use of the preset coefficients in :literal:`data`.
It loops over the mesh points :math:`(i,j)`.  For each such mesh
point, the one-dimensional index :math:`k = j - 1 + (i - 1) * M_Y` is computed
and the :math:`k^{th}` *column* of the *Jacobian matrix* :math:`J` is set.
The row index :math:`k'` of each component :math:`f_{i',j'}` that depends on
:math:`v_{i,j}` must be identified in order to load the corresponding element.
The elements are loaded with the :literal:`SM_COLUMN_ELEMENT_B` macro.
Note that the formula for the global index :math:`k` implies that decreasing 
(increasing) :math:`i` by :math:`1` corresponds to decreasing (increasing) 
:math:`k` by :math:`M_Y`, while decreasing (increasing) :math:`j` by :math:`1` 
corresponds of decreasing (increasing) :math:`k` by :math:`1`. 
These statements are reflected in the arguments to
:literal:`SM_COLUMN_ELEMENT_B`.  The first argument passed to the
:literal:`SM_COLUMN_ELEMENT_B` *macro* is a pointer to the diagonal element
in the *column* to be accessed.  This pointer is obtained via a call to
the :literal:`SM_COLUMN_B` *macro* and is stored in :literal:`kthCol` in 
the :literal:`Jac` function.  When setting the components of :math:`J` we must be
careful not to index out of bounds. The guards :C:`(i != 1)` etc.
in front of the calls to :literal:`SM_COLUMN_ELEMENT_B` prevent illegal indexing.
See `CVode Jacobian Construction 
<https://sundials.readthedocs.io/en/latest/cvode/Usage/index.html#jacobian-construction-matrix-based-linear-solvers>`_
for a detailed description of the :literal:`Jac` function.

Program output
---------------

The output generated by :literal:`cvAdvDiff_bnd` is again shown below.

.. include:: ../../../../examples/cvode/serial/cvAdvDiff_bnd.out
   :literal:


.. _deep_dive.cvDiurnal_kry:

A Krylov example: cvDiurnal_kry [DD]
============================================

Description -- step-by-step
----------------------------

We give here an example that illustrates the use of **CVODE** with the *Krylov*
method **SPGMR**, in the **SUNLinSol_SPGMR** module, as the linear system
solver through the **CVode** *interface*.

This program solves the semi-discretized form of a pair of
kinetics-advection-diffusion partial differential equations, which
represent a simplified model for the transport, production, and loss
of ozone and the oxygen singlet in the upper atmosphere.  The problem
includes nonlinear diurnal kinetics, horizontal advection and diffusion, 
and nonuniform vertical diffusion.  The PDEs can be written as

.. math::
   :label: cvDiurnalpde

   \frac{\partial c^i}{\partial t} = K_h \frac{\partial^2 c^i}{\partial x^2}
    + V \frac{\partial c^i}{\partial x}
    + \frac{\partial} {\partial y} K_v(y) \frac{\partial c^i}{\partial y}
    + R^i(c^1,c^2,t) \quad (i = 1, 2)~,

where the superscripts :math:`i` are used to distinguish the two chemical
species, and where the reaction terms are given by

.. math::
   :label: eq:cvDiurnal:r

   R^1(c^1,c^2,t) & = -q_1 c^1 c^3 - q_2 c^1 c^2 + 2 q_3(t) c^3 + q_4(t) c^2 ~, \\
   R^2(c^1,c^2,t) & = q_1 c^1 c^3 - q_2 c^1 c^2 - q_4(t) c^2 ~.

The spatial domain is :math:`0 \leq x \leq 20,\;30 \leq y \leq 50` (in *km*). 
The various constants and parameters are: 
:math:`K_h = 4.0 \cdot 10^{-6},` :math:`V = 10^{-3},` :math:`K_v = 10^{-8} e^{\frac{y}{5}},`
:math:`q_1 = 1.63 \cdot 10^{-16},` :math:`q_2 = 4.66 \cdot 10^{-16},` :math:`c^3 = 3.7 \cdot 10^{16},`
and the diurnal rate constants are defined as:

.. math::

   q_i(t) = 
   \left\{ \begin{array}{ll}
     \exp \left[ \frac{-a_i}{\sin (\omega t)} \right], & \mbox{for } \sin (\omega t) > 0 \\
     0, & \mbox{for } \sin (\omega t) \leq 0
     \end{array} \right\} ~~~(i = 3, 4) ~,

where :math:`\omega = \frac{\pi}{43200},` :math:`a_3 = 22.62,` :math:`a_4 = 7.601.`  The time interval of
integration is :math:`[0, 86400]`, representing 24 *hours* measured in *seconds*.

Homogeneous Neumann boundary conditions are imposed on each boundary, and the
initial conditions are 

.. math::
   :label: cvDiurnalic

   c^1 (x,y,0) &= 10^6 \alpha (x) \beta (y) ~,~~~ 
                   c^2 (x,y,0) = 10^{12} \alpha (x) \beta (y) ~, \\
   \alpha (x) &= 1 - (0.1 x - 1)^2 + \frac{1}{2} (0.1 x - 1)^4 ~, \\
   \beta (y) &= 1 - (0.1 y - 4)^2 + \frac{1}{2} (0.1 y - 4)^4 ~.

For this example, the equations (:math:numref:`cvDiurnalpde`) are discretized spatially
with standard central finite differences on a :math:`10 \times 10` mesh,
giving an ODE system of size :math:`200`.

Among the initial ``#include`` lines in this case are lines to
include :literal:`sunlinsol_spgmr.h` and :literal:`sundials_math.h`.  The first
contains constants and function prototypes associated with the
**SUNLinSol_SPGMR** module, including the values of the ``pretype``
argument to :literal:`SUNLinSol_SPGMR`. The inclusion of :literal:`sundials_math.h` is
done to access the :literal:`SUNSQR` macro for the square of a :literal:`realtype`
number. 

The :literal:`main` program calls :literal:`CVodeCreate` specifying the ``CV_BDF`` method
and built-in ``CV_NEWTON`` iteration, and then calls :literal:`CVodeInit`, and
:literal:`CVodeSetSStolerances` specifies the *scalar tolerances*.
It calls :literal:`SUNLinSol_SPGMR` to create the **SPGMR** *linear solver* with *left
preconditioning*, and the default value (indicated by a zero argument)
for ``maxl``.  It then calls :literal:`CVodeSetLinearSolver` (see `CVode Linear Solver Functions 
<https://sundials.readthedocs.io/en/latest/cvode/Usage/index.html#linear-solver-interface-functions>`_)
to attach this *linear solver* to the
**CVode** *interface*.  The call to :literal:`CVodeSetJacTimes` specifies a
*user-supplied function* for *Jacobian-vector products* (the ``NULL``
argument specifies that no *Jacobian-vector setup routine* is needed). 
Next, *user-supplied preconditioner setup* and *solve* functions,
:literal:`Precond` and :literal:`PSolve`, are specified. See `CVode Linear Solver Optional Input Functions 
<https://sundials.readthedocs.io/en/latest/cvode/Usage/index.html#linear-solver-interface-optional-input-functions>`_
for details on the :literal:`CVodeSetPreconditioner` function.

For a sequence of :math:`t_{out}` *values*, :literal:`CVode` is called in the
``CV_NORMAL`` mode, sampled output is printed, and the return value is
tested for error conditions.  After that, :literal:`PrintFinalStats` is called
to get and print final statistics, and memory is freed by calls to
:literal:`N_VDestroy`, :literal:`FreeUserData`, and :literal:`CVodeFree`.  The printed
statistics include various counters, such as the total numbers of steps
(``nst``), of :math:`f` evaluations (excluding those for :math:`J \cdot v` product
evaluations) (``nfe``), of :math:`f` evaluations for :math:`J \cdot v` evaluations (``nfeLS``),
of nonlinear iterations (``nni``), of linear (*Krylov*) iterations (``nli``),
of *preconditioner* setups (``nsetups``), of *preconditioner* evaluations
(``npe``), and of *preconditioner* solves (``nps``), among others.  
Also printed are the lengths of the problem-dependent *real and integer
workspaces* used by the main integrator :literal:`CVode`, denoted ``lenrw`` and
``leniw``, and those used by **CVode**, denoted ``lenrwLS`` and ``leniwLS``.
All of these optional outputs are described in `CVode Optional Output Functions 
<https://sundials.readthedocs.io/en/latest/cvode/Usage/index.html#optional-output-functions>`_.
The :literal:`PrintFinalStats` function is suitable for general use in applications
of **CVODE** to any problem with an iterative *linear solver*.

Mathematically, the dependent variable has three dimensions: species
number, :math:`x` mesh point, and :math:`y` mesh point.  But in **NVECTOR_SERIAL**, a vector of
type :literal:`N_Vector` works with a one-dimensional contiguous array of
data components. The *macro* :literal:`IJKth` isolates the translation from
three dimensions to one. Its use results in clearer code and makes it
easy to change the underlying layout of the three-dimensional data. 
Here the problem size is :math:`200`, so we use the :literal:`NV_DATA_S` *macro*
for efficient :literal:`N_Vector` access.  The :literal:`NV_DATA_S` *macro* gives
a pointer to the first component of an :literal:`N_Vector` which we pass to
the :literal:`IJKth` *macro* to do an :literal:`N_Vector` access.

The *preconditioner* used here is the *block-diagonal part* of the true *Newton
matrix*.  It is generated and factored in the :literal:`Precond` routine (see `CVode Preconditioner Setup Function
<https://sundials.readthedocs.io/en/latest/cvode/Usage/index.html#preconditioner-setup-iterative-linear-solvers>`_)
and backsolved in the :literal:`PSolve` routine (see `CVode Preconditioner Solve Function
<https://sundials.readthedocs.io/en/latest/cvode/Usage/index.html#preconditioner-solve-iterative-linear-solvers>`_).
Its diagonal blocks are :math:`2 \times 2` *matrices* that include the interaction
*Jacobian elements* and the *diagonal contribution* of the *diffusion Jacobian
elements*.  The *block-diagonal* part of the *Jacobian* itself, :math:`J_{bd}`, is
saved in separate storage each time it is generated, on calls to :literal:`Precond`
with :C:`jok == SUNFALSE`.  On calls with :C:`jok == SUNTRUE`, signifying that
saved *Jacobian data* can be reused, the *preconditioner* :math:`P = I - \gamma J_{bd}`
is formed from the *saved matrix* :math:`J_{bd}` and factored.  (A call to
:literal:`Precond` with :C:`jok == SUNTRUE` can only occur after a prior call with
:C:`jok == SUNFALSE`.)  The :literal:`Precond` routine must also set the value
of :C:`jcur`, i.e. :C:`*jcurPtr`, to ``SUNTRUE`` when :math:`J_{bd}` is re-evaluated,
and ``SUNFALSE`` otherwise, to inform **CVode** of the status of *Jacobian data*.

We need to take a brief detour to explain one last important aspect of
this program.  While the generic **SUNLinSol_Dense** *linear solver*
module serves as the *interface* to *dense matrix* solves for the main
**SUNDIALS** solvers, the underlying algebraic operations operate on
*dense matrices* with :literal:`realtype **` as the underlying *dense matrix*
type.  To avoid the extra layer of function calls and *dense matrix* and
*linear solver data structures*, :literal:`cvDiurnal_kry.c` uses underlying small
dense functions for all operations on the :math:`2 \times 2` *preconditioner blocks*.  
Thus it includes :literal:`sundials_dense.h`, and calls the *small dense
matrix* functions :literal:`newDenseMat`, :literal:`newIndexArray`, :literal:`denseCopy`,
:literal:`denseScale`, :literal:`denseAddIdentity`, :literal:`denseGETRF`, and
:literal:`denseGETRS`. The *macro* :literal:`IJth` defined near the top of the file
is used to access individual elements in each *preconditioner block*,
numbered from :math:`1`.  The underlying *dense algebra functions* are
available for **CVODE** user programs generally.

In addition to the functions called by **CVODE**, :literal:`cvDiurnal_kry.c` includes
definitions of several private functions.  These are: :literal:`AllocUserData`
to allocate space for :math:`J_{bd}`, :math:`P`, and the pivot arrays; :literal:`InitUserData`
to load problem constants in the ``data`` block; :literal:`FreeUserData` to free
that block; :literal:`SetInitialProfiles` to load the *initial values* in :math:`y`; 
:literal:`PrintOutput` to *retreive* and *print* selected solution values and
statistics; :literal:`PrintFinalStats` to *print* statistics; and :literal:`check_flag`
to *check return values* for error conditions.

Program output
---------------

The output generated by :literal:`cvDiurnal_kry.c` is shown below.  Note that the
number of *preconditioner* evaluations, ``npe``, is much smaller than
the number of *preconditioner* setups, ``nsetups``, as a result of the
*Jacobian* re-use scheme.

.. include:: ../../../../examples/cvode/serial/cvDiurnal_kry.out
   :literal:

