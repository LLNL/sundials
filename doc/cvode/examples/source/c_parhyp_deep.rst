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


.. _parhyp_deep_c:

==================================================
Parallel HYPRE C example problems -- Deep dive
==================================================



.. _deep_dive.cvAdvDiff_non_ph:

A nonstiff example: cvAdvDiff_non_ph [DD]
======================================================

This example is same as ``cvAdvDiff_non_p``, except that it
uses the *hypre* vector type instead of the SUNDIALS native
*parallel vector* implementation.

The outputs from the two examples are identical. In the following, we will point
out only the differences between the two. Familiarity with `HYPRE <https://github.com/hypre-space/hypre>`_
library at `HYPRE User Manual <https://hypre.readthedocs.io/en/latest/>`_ is helpful.

We use the *hypre* IJ vector interface to allocate the template vector and
create parallel partitioning:

.. code-block:: c

   HYPRE_IJVectorCreate(comm, my_base, my_base + local_N - 1, &Uij);
   HYPRE_IJVectorSetObjectType(Uij, HYPRE_PARCSR);
   HYPRE_IJVectorInitialize(Uij);

The initialize call means that vector elements are ready to be set using
the IJ interface. We choose an *initial condition vector* :math:`x_0 = x(t_0)` as the
template vector and we set its values in the :code:`SetIC(...)` function. We
complete the *hypre* vector assembly by:

.. code-block:: c

   HYPRE_IJVectorAssemble(Uij);
   HYPRE_IJVectorGetObject(Uij, (void**) &Upar);

The assemble call is collective and it makes the *hypre* vector ready to use.
This sets the handle :code:`Upar` to the actual *hypre* vector.
The handle is then passed to the :code:`N_VMake` function, which creates
the template :code:`N_Vector` as a *wrapper* around the *hypre* vector.
All other vectors in the computation are created by cloning the template
vector. The template vector does not own the underlying *hypre* vector,
and it is the user's responsibility to destroy it using a
:code:`HYPRE_IJVectorDestroy(Uij)` call after the template vector has been
destroyed. This function will destroy both the *hypre* vector and its IJ
interface.

To access individual elements of solution vectors :math:`u =` ``u`` and :math:`\dot u =` ``udot``
in the *residual function*, the user needs to extract the *hypre* vector first
by calling :code:`N_VGetVector_ParHyp`, and then use *hypre* methods from
that point on.

.. warning::

   At this point interfaces to *hypre* solvers and preconditioners are
   not available. They will be provided in subsequent SUNDIALS releases.
   The interface to the *hypre* vector is included in this release mainly for
   testing purposes and as a preview of functionality to come.
