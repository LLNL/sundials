.. ----------------------------------------------------------------
   Programmer(s): Steven B. Roberts @ LLNL
   ----------------------------------------------------------------
   SUNDIALS Copyright Start
   Copyright (c) 2002-2024, Lawrence Livermore National Security
   and Southern Methodist University.
   All rights reserved.

   See the top-level LICENSE and NOTICE files for details.

   SPDX-License-Identifier: BSD-3-Clause
   SUNDIALS Copyright End
   ----------------------------------------------------------------

.. _ARKODE.Usage.SplittingStep.SplittingStepCoefficients:

Operator Splitting Coefficients Data Structure
----------------------------------------------

SplittingStep supplies several functions to construct operator splitting
coefficients of various orders and partitions. There are also a number of 
built-in methods of fixed orders and partitions (see
:numref:`ARKODE.Usage.SplittingStep.SplittingStepCoefficients.Coefficients`).
Finally, a user may construct a custom set of coefficients and attach it with
:c:func:`SplittingStepSetCoefficients`. The operator splitting coefficients are
stored in a :c:func:`SplittingStepCoefficients` object which is a pointer to a
:c:struct:`SplittingStepCoefficientsMem` structure:

.. c:type:: SplittingStepCoefficientsMem *SplittingStepCoefficients

.. c:struct::  SplittingStepCoefficientsMem

   Structure for storing the coefficients defining an operator splitting method.

   As described in :numref:`ARKODE.Mathematics.SplittingStep`, an operator
   splitting method is defined by a vector :math:`\alpha \in \mathbb{R}^{r}` and
   a tensor :math:`\beta \in \mathbb{R}^{r \times (s + 1) \times P}`.

   .. c:member:: sunrealtype *alpha

      An array of length ``[sequential_methods]`` containing the weight of each
      sequential method used to produce the overall operator splitting solution

   .. c:member:: sunrealtype ***beta

      A three-dimensional array with dimensions
      ``[sequential_methods][stages+1][partitions]`` containing the time nodes
      of the inner integrations.

   .. c:member:: int sequential_methods

      The number :math:`r` of sequential methods combined to produce the overall
      operator splitting solution

   .. c:member:: int stages

      The number :math:`s` of stages

   .. c:member:: int partitions

      The number :math:`P` of partitions in the IVP

   .. c:member:: int order

      The method order of accuracy


.. _ARKODE.Usage.SplittingStep.SplittingStepCoefficients.Functions:

SplittingStepCoefficients Functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


This section describes the functions for creating and interacting with operator
splitting coefficients. The function prototypes and as well as the relevant
integer constants are defined ``arkode/arkode_splitting_coefficients.h``.

.. _ARKODE.Usage.SplittingStep.SplittingStepCoefficients.Functions.Table:
.. table:: SplittingStepCoefficients functions

   +--------------------------------------------------------------+-------------------------------------------------------------+
   | Function name                                                | Description                                                 |
   +--------------------------------------------------------------+-------------------------------------------------------------+
   | :c:func:`SplittingStepCoefficients_LoadCoefficients()`       | Load a pre-defined SplittingStepCoefficients by ID          |
   +--------------------------------------------------------------+-------------------------------------------------------------+
   | :c:func:`SplittingStepCoefficients_LoadCoefficientsByName()` | Load a pre-defined SplittingStepCoefficients by name        |
   +--------------------------------------------------------------+-------------------------------------------------------------+
   | :c:func:`SplittingStepCoefficients_IDToName()`               | Convert a pre-defined SplittingStepCoefficients to its name |
   +--------------------------------------------------------------+-------------------------------------------------------------+
   | :c:func:`SplittingStepCoefficients_LieTrotter()`             | Create a Lie-Trotter splitting method                       |
   +--------------------------------------------------------------+-------------------------------------------------------------+
   | :c:func:`SplittingStepCoefficients_Strang()`                 | Create a Strang splitting method                            |
   +--------------------------------------------------------------+-------------------------------------------------------------+
   | :c:func:`SplittingStepCoefficients_SymmetricParallel()`      | Create a symmetrization of the Lie-Trotter splitting method |
   +--------------------------------------------------------------+-------------------------------------------------------------+
   | :c:func:`SplittingStepCoefficients_ThirdOrderSuzuki()`       | Create a third order composition method of Suzuki           |
   +--------------------------------------------------------------+-------------------------------------------------------------+
   | :c:func:`SplittingStepCoefficients_TripleJump()`             | Create an arbitrary order, three-jump composition method    |
   +--------------------------------------------------------------+-------------------------------------------------------------+
   | :c:func:`SplittingStepCoefficients_SuzukiFractal()`          | Create an arbitrary order, five-jump composition method     |
   +--------------------------------------------------------------+-------------------------------------------------------------+
   | :c:func:`SplittingStepCoefficients_Alloc()`                  | Allocate an empty SplittingStepCoefficient                  |
   +--------------------------------------------------------------+-------------------------------------------------------------+
   | :c:func:`SplittingStepCoefficients_Create()`                 | Create a new SplittingStepCoefficient from coefficients     |
   +--------------------------------------------------------------+-------------------------------------------------------------+
   | :c:func:`SplittingStepCoefficients_Copy()`                   | Create a copy of a SplittingStepCoefficients                |
   +--------------------------------------------------------------+-------------------------------------------------------------+
   | :c:func:`SplittingStepCoefficients_Free()`                   | Deallocate a SplittingStepCoefficients                      |
   +--------------------------------------------------------------+-------------------------------------------------------------+
   | :c:func:`SplittingStepCoefficients_Write()`                  | Write the SplittingStepCoefficients to an output file       |
   +--------------------------------------------------------------+-------------------------------------------------------------+


.. c:function:: SplittingStepCoefficients SplittingStepCoefficients_LoadCoefficients(ARKODE_SplittingCoefficientsID method)

   Retrieves specified splitting coefficients. For further information on the
   current set of splitting coefficients and their corresponding identifiers,
   see
   :numref:`ARKODE.Usage.SplittingStep.SplittingStepCoefficients.Coefficients`.


   **Arguments:**
      * ``method`` -- the splitting coefficients identifier.

   **Return value:**
      * A :c:type:`SplittingStepCoefficients` structure if successful.
      * A ``NULL`` pointer if *method* was invalid or an allocation error occurred.
   
   .. versionadded:: x.y.z



.. c:function:: SplittingStepCoefficients SplittingStepCoefficients_LoadCoefficientsByName(const char *method)

   Retrieves specified splitting coefficients. For further information on the
   current set of splitting coefficients and their corresponding name, see
   :numref:`ARKODE.Usage.SplittingStep.SplittingStepCoefficients.Coefficients`.


   **Arguments:**
      * ``method`` -- the splitting coefficients identifier.

   **Return value:**
      * A :c:type:`SplittingStepCoefficients` structure if successful.
      * A ``NULL`` pointer if *method* was invalid, *method* was
        ``"ARKODE_SPLITTING_NONE"``, or an allocation error occurred.

   .. note::

      This function is case sensitive.
   
   .. versionadded:: x.y.z


.. c:function:: const char* SplittingStepCoefficients_LoadCoefficientsByName(ARKODE_SplittingCoefficientsID method)

   Converts specified splitting coefficients ID to a string of the same name.
   For further information on the current set of splitting coefficients and
   their corresponding name, see
   :numref:`ARKODE.Usage.SplittingStep.SplittingStepCoefficients.Coefficients`.

   **Arguments:**
      * *method* -- the splitting coefficients identifier.

   **Return value:**
      * The name associated with *method*.
      * ``NULL`` pointer if *method* was invalid.
   
   .. versionadded:: x.y.z


.. c:function:: SplittingStepCoefficients SplittingStepCoefficients_LieTrotter(int partitions)

   Create the coefficients for the first order Lie-Trotter splitting

   .. math::
      y_n = L_h(y_{n-1}) = \left( \phi^P_{h} \circ \phi^{P-1}_{h}
      \circ \dots \circ \phi^1_{h} \right) (y_{n-1})
      :label: ARKODE_Lie-Trotter

   **Arguments:**
      * *partitions* -- The number :math:`P > 1` of partitions in the IVP.

   **Return value:**
      * A :c:type:`SplittingStepCoefficients` structure if successful.
      * ``NULL`` if ``partitions`` was invalid or an allocation error occurred.
   
   .. versionadded:: x.y.z


.. c:function:: SplittingStepCoefficients SplittingStepCoefficients_Strang(int partitions)

   Create the coefficients for the second order Strang splitting

   .. math::
      y_n = S_h(y_{n-1}) = \left( L^*_{h/2} \circ L_{h/2} \right) (y_{n-1}),
      :label: ARKODE_Strang

   where :math:`L` is the Lie-Trotter splitting :eq:`ARKODE_Lie` and
   :math:`L*_h = L^{-1}_{-h}` is its adjoint.

   **Arguments:**
      * *partitions* -- The number :math:`P > 1` of partitions in the IVP.

   **Return value:**
      * A :c:type:`SplittingStepCoefficients` structure if successful.
      * ``NULL`` if ``partitions`` was invalid or an allocation error occurred.
   
   .. versionadded:: x.y.z


.. c:function:: SplittingStepCoefficients SplittingStepCoefficients_Parallel(int partitions)

   Create the coefficients for the first order splitting method

   .. math::
      y_n = \phi^1_h(y_{n-1}) + \phi^2_h(y_{n-1}) + \dots + \phi^P(y_{n-1}) +
      (1 - p) y_{n-1},

   **Arguments:**
      * *partitions* -- The number :math:`P > 1` of partitions in the IVP.

   **Return value:**
      * A :c:type:`SplittingStepCoefficients` structure if successful.
      * ``NULL`` if ``partitions`` was invalid or an allocation error occurred.
   
   .. versionadded:: x.y.z


.. c:function:: SplittingStepCoefficients SplittingStepCoefficients_SymmetricParallel(int partitions)

   Create the coefficients for the second order, symmetrized Lie-Trotter
   splitting

   .. math::
      y_n = \frac{1}{2} \left( L_h(y_{n-1}) + L^*_h(y_{n-1}) \right),

   where :math:`L` is the Lie-Trotter splitting :eq:`ARKODE_Lie` and
   :math:`L^*_h = L^{-1}_{-h}` is its adjoint.
   
   **Arguments:**
      * *partitions* -- The number :math:`P > 1` of partitions in the IVP.

   **Return value:**
      * A :c:type:`SplittingStepCoefficients` structure if successful.
      * ``NULL`` if ``partitions`` was invalid or an allocation error occurred.
   
   .. versionadded:: x.y.z


.. c:function:: SplittingStepCoefficients SplittingStepCoefficients_ThirdOrderSuzuki(int partitions)

   Create the coefficients for a splitting method of Suzuki :cite:p:`Suzuki:92`

   .. math::
      y_n = \left( L_{p_1 h} \circ L^*_{p_2 h} \circ L_{p_3 h} \circ L^*_{p_4 h}
      \circ L_{p_5 h} \right) (y_{n-1}),

   where :math:`L` is the Lie-Trotter splitting :eq:`ARKODE_Lie` and
   :math:`L^*_h = L^{-1}_{-h}` is its adjoint. The parameters
   :math:`p_1, \dots, p_5` are selected to give third order.
   
   **Arguments:**
      * *partitions* -- The number :math:`P > 1` of partitions in the IVP.

   **Return value:**
      * A :c:type:`SplittingStepCoefficients` structure if successful.
      * ``NULL`` if ``partitions`` was invalid or an allocation error occurred.
   
   .. versionadded:: x.y.z


.. c:function:: SplittingStepCoefficients SplittingStepCoefficients_TripleJump(int partitions, int order)

   Create the coefficients for the triple jump splitting method
   :cite:p:`CrGo:89`

   .. math::
      \begin{align*}
      T_h^{[2]} &= S, \\
      T_h^{[i+2]} &= T_{\gamma_1 h}^{[i]} \circ T_{(1 - 2 \gamma_1) h}^{[i]}
      \circ T_{\gamma_1 h}^{[i]}, \\
      y_n &= T_h^{[order]}(y_{n-1}),
      \end{align*}
   
   where :math:`S` is the Strang splitting :eq:`ARKODE_Stang` and parameter
   :math:`\gamma_1` selected to increase the order by two each recursion.
   
   **Arguments:**
      * *partitions* -- The number :math:`P > 1` of partitions in the IVP.
      * *order* -- A positive even number for the method order of accuracy.

   **Return value:**
      * A :c:type:`SplittingStepCoefficients` structure if successful.
      * ``NULL`` if ``partitions`` or ``order`` was invalid or an allocation
        error occurred.
   
   .. versionadded:: x.y.z


.. c:function:: SplittingStepCoefficients SplittingStepCoefficients_SuzukiFractal(int partitions, int order)

   Create the coefficients for the quintuple jump splitting method
   :cite:p:`Suzuki:90`

   .. math::
      \begin{align*}
      Q_h^{[2]} &= S, \\
      Q_h^{[i+2]} &= Q_{\gamma_1 h}^{[i]} \circ Q_{\gamma_1 h}^{[i]} \circ
      Q_{(1 - 4 \gamma_1) h}^i \circ Q_{\gamma_1 h}^{[i]} \circ
      Q_{\gamma_1 h}^{[i]}, \\
      y_n &= Q_h^{[order]}(y_{n-1}),
      \end{align*}
   
   where :math:`S` is the Strang splitting :eq:`ARKODE_Stang` and parameter
   :math:`\gamma_1` selected to increase the order by two each recursion.
   
   **Arguments:**
      * *partitions* -- The number :math:`P > 1` of partitions in the IVP.
      * *order* -- A positive even number for the method order of accuracy.

   **Return value:**
      * A :c:type:`SplittingStepCoefficients` structure if successful.
      * ``NULL`` if ``partitions`` or ``order`` was invalid or an allocation
        error occurred.
   
   .. versionadded:: x.y.z


.. c:function:: SplittingStepCoefficients SplittingStepCoefficients_Alloc(int sequential_methods, int stages, int partitions)

   Allocates an empty SplittingStepCoefficients.

   **Arguments:**
      * *sequential_methods* -- The number :math:`r` of sequential methods
        combined to produce the overall operator splitting solution.
      * *stages* -- The number :math:`s` of stages.
      * *partitions* -- The number :math:`P` of partitions in the IVP.

   **Return value:**
      * An :c:type:`SplittingStepCoefficients` structure if successful.
      * A ``NULL`` pointer if *sequential_methods*, *stages* or *partitions* was
        invalid or an allocation error occurred.
   
   .. versionadded:: x.y.z


.. c:function:: SplittingStepCoefficients SplittingStepCoefficients_Create(int sequential_methods, int stages, int partitions, int order, sunrealtype* alpha, sunrealtype* beta)

   Allocates a SplittingStepCoefficients and fills it with the given values.

   **Arguments:**
      * *sequential_methods* -- The number :math:`r` of sequential methods
        combined to produce the overall operator splitting solution.
      * *stages* -- The number :math:`s` of stages.
      * *partitions* -- The number :math:`P` of partitions in the IVP.
      * *order* -- The method order of accuracy.
      * *alpha* -- An array of length ``sequential_methods`` containing the
        weight of each sequential method used to produce the overall operator
        splitting solution.
      * *beta* -- An array of length
        ``sequential_methods * (stages+1) * partitions`` containing the time nodes
        of the inner integrations in the order

        .. math::
           \beta_{1,0,1}, \dots, \beta_{1,0,P}, \beta_{1,1,1}, \dots, \beta_{1,1,P}, \dots, \beta_{2,0,1}, \dots, \beta_{r,s,P}.

   **Return value:**
      * An :c:type:`SplittingStepCoefficients` structure if successful.
      * A ``NULL`` pointer if an argument was invalid or an allocation error
        occurred.


.. c:function:: SplittingStepCoefficients SplittingStepCoefficients_Copy(SplittingStepCoefficients coefficients)

   Creates copy of the given splitting coefficients.

   **Arguments:**
      * ``coefficients`` -- The splitting coefficients to copy.

   **Return value:**
      * An :c:type:`SplittingStepCoefficients` structure if successful.
      * A ``NULL`` pointer if an allocation error occurred.


.. c:function:: void SplittingStepCoefficients_Free(SplittingStepCoefficients coefficients)

   Deallocate the splitting coefficients memory.

   **Arguments:**
      * ``coefficients`` -- The splitting coefficients.


.. c:function:: void SplittingStepCoefficients_Write(SplittingStepCoefficients coefficients, FILE* outfile)

   Write the splitting coefficients to the provided file pointer.

   **Arguments:**
      * ``coefficients`` -- The splitting coefficients.
      * ``outfile`` -- Pointer to use for printing the splitting coefficients.

   .. note::

      The *outfile* argument can be ``stdout`` or ``stderr``, or it may point to
      a specific file created using ``fopen``.


.. _ARKODE.Usage.SplittingStep.SplittingStepCoefficients.Coefficients:

Operator Splitting Coefficients
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

SplittingStep currently provides several pre-defined coefficients for problems
with two partitions.  We list the identifiers, order of accuracy, and relevant
references for each in the table below. We use the naming convention

.. code-block:: text

   <name>_<stages>_<order>_<partitions>

Each of the splitting coefficients that are packaged with SplittingStep are
specified by a unique ID having type:

.. c:type:: int ARKODE_SplittingCoefficientsID

with values specified for each method below (e.g.,
``ARKODE_SPLITTING_LIE_TROTTER_1_1_2``).

.. table:: Operator splitting coefficients.

   ======================================  =====  =====================
   Table name                              Order        Reference
   ======================================  =====  =====================
   ``ARKODE_SPLITTING_LIE_TROTTER_1_1_2``  1      
   ``ARKODE_SPLITTING_STRANG_2_2_2``       2      :cite:p:`Strang:68`
   ``ARKODE_SPLITTING_BEST_2_2_2``         2      :cite:p:`AuHoKeKo:16`
   ``ARKODE_SPLITTING_SUZUKI_3_3_2``       3      :cite:p:`Suzuki:92`
   ``ARKODE_SPLITTING_RUTH_3_3_2``         3      :cite:p:`Ruth:93`
   ``ARKODE_SPLITTING_YOSHIDA_4_4_2``      4      :cite:p:`Yoshida:90`
   ``ARKODE_SPLITTING_YOSHIDA_8_6_2``      6      :cite:p:`Yoshida:90`
   ======================================  =====  =====================
