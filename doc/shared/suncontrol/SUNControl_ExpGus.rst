..
   Programmer(s): Daniel R. Reynolds @ SMU
   ----------------------------------------------------------------
   SUNDIALS Copyright Start
   Copyright (c) 2002-2023, Lawrence Livermore National Security
   and Southern Methodist University.
   All rights reserved.

   See the top-level LICENSE and NOTICE files for details.

   SPDX-License-Identifier: BSD-3-Clause
   SUNDIALS Copyright End
   ----------------------------------------------------------------

.. _SUNControl.ExpGus:

The SUNControl_ExpGus Module
======================================

The explicit Gustafsson implementation of the SUNControl class, SUNControl_ExpGus, implements a controller proposed by K. Gustafsson for explicit Runge--Kutta methods in :cite:p:`Gust:91`.  This controller has the form

.. math::
   h' \;=\; \begin{cases}
      h_1\; \varepsilon_1^{-1/p}, &\quad\text{on the first step}, \\
      h_n\; \varepsilon_n^{-k_1/p}\;
        \left(\dfrac{\varepsilon_n}{\varepsilon_{n-1}}\right)^{k_2/p}, &
      \quad\text{on subsequent steps},
   \end{cases}

with default values :math:`k_1=0.367` and :math:`k_2=0.268`.

The SUNControl_ExpGus controller is implemented as a derived SUNControl class, and defines its *content* field as:

.. code-block:: c

   struct _SUNControlContent_ExpGus {
     realtype k1;
     realtype k2;
     realtype bias;
     realtype ep;
     int p;
     sunbooleantype pq;
     sunbooleantype firststep;
   };

These entries of the *content* field contain the following information:

* ``k1, k2`` - controller parameters above.

* ``bias`` - error bias factor, that converts from an input temporal error estimate via :math:`\varepsilon = \text{bias}*\text{dsm}`.

* ``ep`` - storage for the previous error estimate, :math:`\varepsilon_{n-1}`.

* ``p`` - asymptotic order to use in error control.

* ``pq`` - flag indicating whether ``p`` corresponds to the order of accuracy for the time integration method (``SUNTRUE``) or the embedding (``SUNFALSE``).


The header file to be included when using this module is
``suncontrol/suncontrol_expgus.h``.


The SUNControl_ExpGus class provides implementations of all controller operations listed in :numref:`SUNControl.Description.operations`. The SUNControl_ExpGus class also provides the following additional user-callable routines:


.. c:function:: SUNControl SUNControlExpGus(SUNContext sunctx)

   This constructor function creates and allocates memory for a SUNControl_ExpGus object, and inserts its default parameters.  The only argument is the SUNDIALS context object.  Upon successful completion it will return a :c:type:`SUNControl` object; otherwise it will return ``NULL``.


.. c:function:: int SUNControlExpGus_SetParams(SUNControl C, sunbooleantype pq, realtype k1, realtype k2)

   This user-callable function provides control over the relevant parameters above.  The *pq* input is stored directly.  The *k1* and *k2* parameters are only stored if the corresponding input is non-negative.  Upon completion, this returns ``SUNCONTROL_SUCCESS``.
