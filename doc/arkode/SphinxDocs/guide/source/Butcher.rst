..
   Programmer(s): Daniel R. Reynolds @ SMU
   ----------------------------------------------------------------
   Copyright (c) 2013, Southern Methodist University.
   All rights reserved.
   For details, see the LICENSE file.
   ----------------------------------------------------------------

:tocdepth: 3


.. _Butcher:

=========================
Appendix: Butcher tables
=========================

Here we catalog the full set of Butcher tables included in ARKode.
We group these into three categories: *explicit*, *implicit* and
*additive*.  However, since the methods that comprise an additive
Runge Kutta method are themselves explicit and implicit, their
component Butcher tables are listed within their separate
sections, but are referenced together in the additive section.

In each of the following tables, we use the following notation (shown
for a 3-stage method):

.. math::

   \begin{array}{r|ccc}
     c_1 & a_{1,1} & a_{1,2} & a_{1,3} \\
     c_2 & a_{2,1} & a_{2,2} & a_{2,3} \\
     c_3 & a_{3,1} & a_{3,2} & a_{3,3} \\
     \hline
     q & b_1 & b_2 & b_3 \\
     p & \tilde{b}_1 & \tilde{b}_2 & \tilde{b}_3
   \end{array}

where here the method and embedding share stage :math:`A` and
:math:`c` values, but use their stages :math:`z_i` differently through
the coefficients :math:`b` and :math:`\tilde{b}` to generate methods
of orders :math:`q` (the main method) and :math:`p` (the embedding,
typically :math:`q = p+1`, though sometimes this is reversed).

Method authors often use different naming conventions to categorize
their methods.  For each of the methods below, we follow a uniform
naming convention:

.. code-block:: text

   NAME-S-P-Q

where here

* ``NAME`` is the author or the name provided by the author (if applicable),
* ``S`` is the number of stages in the method,
* ``P`` is the global order of accuracy for the embedding,
* ``Q`` is the global order of accuracy for the method.

In the code, unique integer IDs are defined inside ``arkode_butcher.h`` for
each method, which may be used by calling routines to specify the
desired method.  These names are specified in ``fixed width font`` at
the start of each method's section below.

Additionally, for each method we provide a plot of the linear
stability region in the complex plane.  These have been computed via
the following approach.  For any Runge Kutta method as defined above,
we may define the stability function

.. math::

   R(\eta) = 1 + \eta b [I - \eta A]^{-1} e,

where :math:`e\in\mathbb{R}^s` is a column vector of all ones, :math:`\eta =
h\lambda` and :math:`h` is the time step size.  If the stability
function satisfies :math:`|R(\eta)| \le 1` for all eigenvalues,
:math:`\lambda`, of :math:`\frac{\partial }{\partial y}f(t,y)` for a
given IVP, then the method will be linearly stable for that problem
and step size.  The stability region 

.. math::

   S = \{ \eta\in\mathbb{C}\; :\; \left| R(\eta) \right| \le 1\}

is typically given by an enclosed region of the complex plane, so it
is standard to search for the border of that region in order to
understand the method.  Since all complex numbers with unit magnitude
may be written as :math:`e^{i\theta}` for some value of :math:`\theta`,
we perform the following algorithm to trace out this boundary.

1. Define an array of values ``Theta``.  Since we wish for a
   smooth curve, and since we wish to trace out the entire boundary,
   we choose 10,000 linearly-spaced points from 0 to :math:`16\pi`.
   Since some angles will correspond to multiple locations on the
   stability boundary, by going beyond :math:`2\pi` we ensure that all
   boundary locations are plotted, and by using such a fine
   discretization the Newton method (next step) is more likely to
   converge to the root closest to the previous boundary point,
   ensuring a smooth plot.

2. For each value :math:`\theta \in` ``Theta``, we solve the nonlinear
   equation 

   .. math::

      0 = f(\eta) = R(\eta) - e^{i\theta}

   using a finite-difference Newton iteration, using tolerance
   :math:`10^{-7}`, and differencing parameter
   :math:`\sqrt{\varepsilon}` (:math:`\approx 10^{-8}`).

   In this iteration, we use as initial guess the solution from the
   previous value of :math:`\theta`, starting with an initial-initial
   guess of :math:`\eta=0` for :math:`\theta=0`.  

3. We then plot the resulting :math:`\eta` values that trace the
   stability region boundary.

We note that for any stable IVP method, the value :math:`\eta_0 =
-\varepsilon + 0i` is always within the stability region.  So in each
of the following pictures, the interior of the stability region is the
connected region that includes :math:`\eta_0`.  Resultingly, methods
whose linear stability boundary is located entirely in the right
half-plane indicate an `A-stable` method.



.. _Butcher.explicit:

Explicit Butcher tables
---------------------------

In the category of explicit Runge-Kutta methods, ARKode includes
methods that have orders 2 through 6, with embeddings that are of
orders 1 through 5.


.. _Butcher.Heun_Euler:

Heun-Euler-2-1-2 
^^^^^^^^^^^^^^^^^^^^

.. index:: Heun-Euler-2-1-2 ERK method

Accessible via the constant ``HEUN_EULER_2_1_2`` to
:c:func:`ARKStepSetERKTableNum()`, :c:func:`ERKStepSetERKTableNum()` or
:c:func:`ARKodeLoadButcherTable()`.  This is the default 2nd order
explicit method.  

.. math::

   \begin{array}{r|cc}
     0 & 0 & 0 \\
     1 & 1 & 0 \\
     \hline
     2 & \frac{1}{2} & \frac{1}{2} \\
     1 & 1 & 0
   \end{array}

.. figure:: figs/stab_region_0.png
   :scale: 50 %
   :align: center
   
   Linear stability region for the Heun-Euler method.  The method's
   region is outlined in blue; the embedding's region is in red.



.. _Butcher.Bogacki_Shampine:

Bogacki-Shampine-4-2-3 
^^^^^^^^^^^^^^^^^^^^^^^^^^

.. index:: Bogacki-Shampine-4-2-3 ERK method

Accessible via the constant ``BOGACKI_SHAMPINE_4_2_3`` to
:c:func:`ARKStepSetERKTableNum()`, :c:func:`ERKStepSetERKTableNum()` or
:c:func:`ARKodeLoadButcherTable()`.  This is the default 3rd order 
explicit method (from [BS1989]_). 

.. math::

   \begin{array}{r|cccc}
     0 &   0 & 0 & 0 & 0 \\
     \frac{1}{2} & \frac{1}{2} & 0 & 0 & 0 \\
     \frac{3}{4} & 0 & \frac{3}{4} & 0 & 0 \\
     1   & \frac{2}{9} & \frac{1}{3} & \frac{4}{9} & 0 \\
     \hline
     3 & \frac{2}{9} & \frac{1}{3} & \frac{4}{9} \\
     2 & \frac{7}{24} & \frac{1}{4} & \frac{1}{3} & \frac{1}{8}
   \end{array}

.. figure:: figs/stab_region_1.png
   :scale: 50 %
   :align: center
   
   Linear stability region for the Bogacki-Shampine method.  The method's
   region is outlined in blue; the embedding's region is in red.




.. _Butcher.ARK_4_2_3_E:

ARK-4-2-3 (explicit) 
^^^^^^^^^^^^^^^^^^^^^^^^

.. index:: ARK-4-2-3 ERK method

Accessible via the constant ``ARK324L2SA_ERK_4_2_3`` to
:c:func:`ARKStepSetERKTableNum()`, :c:func:`ERKStepSetERKTableNum()` or
:c:func:`ARKodeLoadButcherTable()`. This is the explicit portion of
the default 3rd order additive method (from [KC2003]_). 

.. math::

   \begin{array}{r|cccc}
     0 & 0 & 0 & 0 & 0 \\
     \frac{1767732205903}{2027836641118} & \frac{1767732205903}{2027836641118} & 0 & 0 & 0 \\
     \frac{3}{5} & \frac{5535828885825}{10492691773637} & \frac{788022342437}{10882634858940} & 0 & 0 \\
     1 & \frac{6485989280629}{16251701735622} & -\frac{4246266847089}{9704473918619} & \frac{10755448449292}{10357097424841} & 0 \\
     \hline
     3 & \frac{1471266399579}{7840856788654} & -\frac{4482444167858}{7529755066697} & \frac{11266239266428}{11593286722821} & \frac{1767732205903}{4055673282236} \\
     2 & \frac{2756255671327}{12835298489170} & -\frac{10771552573575}{22201958757719} & \frac{9247589265047}{10645013368117} & \frac{2193209047091}{5459859503100}
   \end{array}

.. figure:: figs/stab_region_2.png
   :scale: 50 %
   :align: center
   
   Linear stability region for the explicit ARK-4-2-3 method.  The method's
   region is outlined in blue; the embedding's region is in red.





.. _Butcher.Zonneveld:

Zonneveld-5-3-4 
^^^^^^^^^^^^^^^^^^

.. index:: Zonneveld-5-3-4 ERK method

Accessible via the constant ``ZONNEVELD_5_3_4`` to
:c:func:`ARKStepSetERKTableNum()`, :c:func:`ERKStepSetERKTableNum()`
or :c:func:`ARKodeLoadButcherTable()`.  This is  
the default 4th order explicit method (from [Z1963]_).

.. math::

   \begin{array}{r|ccccc}
       0 & 0 & 0 & 0 & 0 & 0 \\
     \frac{1}{2} & \frac{1}{2} & 0 & 0 & 0 & 0 \\
     \frac{1}{2} & 0 & \frac{1}{2} & 0 & 0 & 0 \\
       1 & 0 & 0 & 1 & 0 & 0 \\
     \frac{3}{4} & \frac{5}{32} & \frac{7}{32} & \frac{13}{32} & -\frac{1}{32} & 0 \\
     \hline
     4 & \frac{1}{6} & \frac{1}{3} & \frac{1}{3} & \frac{1}{6} & 0 \\
     3 & -\frac{1}{2} & \frac{7}{3} & \frac{7}{3} & \frac{13}{6} & -\frac{16}{3}
   \end{array}

.. figure:: figs/stab_region_3.png
   :scale: 50 %
   :align: center
   
   Linear stability region for the Zonneveld method.  The method's
   region is outlined in blue; the embedding's region is in red.





.. _Butcher.ARK_6_3_4_E:

ARK-6-3-4 (explicit) 
^^^^^^^^^^^^^^^^^^^^^^^^

.. index:: ARK-6-3-4 ERK method

Accessible via the constant ``ARK436L2SA_ERK_6_3_4`` to
:c:func:`ARKStepSetERKTableNum()`, :c:func:`ERKStepSetERKTableNum()`
or :c:func:`ARKodeLoadButcherTable()`.  This is the explicit portion
of the default 4th order additive method (from [KC2003]_). 

.. math::

   \begin{array}{r|cccccc}
     0 & 0 & 0 & 0 & 0 & 0 & 0 \\
     \frac12 & \frac12 & 0 & 0 & 0 & 0 & 0 \\
     \frac{83}{250} & \frac{13861}{62500} & \frac{6889}{62500} & 0 & 0 & 0 & 0 \\
     \frac{31}{50} & -\frac{116923316275}{2393684061468} & -\frac{2731218467317}{15368042101831} & \frac{9408046702089}{11113171139209} & 0 & 0 & 0 \\
     \frac{17}{20} & -\frac{451086348788}{2902428689909} & -\frac{2682348792572}{7519795681897} & \frac{12662868775082}{11960479115383} & \frac{3355817975965}{11060851509271} & 0 & 0 \\
     1 & \frac{647845179188}{3216320057751} & \frac{73281519250}{8382639484533} & \frac{552539513391}{3454668386233} & \frac{3354512671639}{8306763924573} & \frac{4040}{17871} & 0 \\
     \hline
     4 & \frac{82889}{524892} & 0 & \frac{15625}{83664} & \frac{69875}{102672} & -\frac{2260}{8211} & \frac14 \\
     3 & \frac{4586570599}{29645900160} & 0 & \frac{178811875}{945068544} & \frac{814220225}{1159782912} & -\frac{3700637}{11593932} & \frac{61727}{225920}
   \end{array}

.. figure:: figs/stab_region_4.png
   :scale: 50 %
   :align: center
   
   Linear stability region for the explicit ARK-6-3-4 method.  The method's
   region is outlined in blue; the embedding's region is in red.




.. _Butcher.Sayfy_Aburub:

Sayfy-Aburub-6-3-4 
^^^^^^^^^^^^^^^^^^^^^

.. index:: Sayfy-Aburub-6-3-4 ERK method

Accessible via the constant ``SAYFY_ABURUB_6_3_4`` to
:c:func:`ARKStepSetERKTableNum()`, :c:func:`ERKStepSetERKTableNum()`
or :c:func:`ARKodeLoadButcherTable()` (from [SA2002]_). 

.. math::

   \begin{array}{r|cccccc}
     0 & 0 & 0 & 0 & 0 & 0 & 0 \\
     \frac{1}{2} & \frac{1}{2} & 0 & 0 & 0 & 0 & 0 \\
     1 & -1 & 2 & 0 & 0 & 0 & 0 \\
     1 & \frac{1}{6} & \frac{2}{3} & \frac{1}{6} & 0 & 0 & 0 \\
     \frac{1}{2} & 0.137 & 0.226 & 0.137 & 0 & 0 & 0 \\
     1 & 0.452 & -0.904 & -0.548 & 0 & 2 & 0 \\
     \hline
     4 & \frac{1}{6} & \frac{1}{3} & \frac{1}{12} & 0 & \frac{1}{3} & \frac{1}{12} \\
     3 & \frac{1}{6} & \frac{2}{3} & \frac{1}{6} & 0 & 0 & 0
   \end{array}

.. figure:: figs/stab_region_5.png
   :scale: 50 %
   :align: center
   
   Linear stability region for the Sayfy-Aburub-6-3-4 method.  The method's
   region is outlined in blue; the embedding's region is in red.





.. _Butcher.Cash-Karp:

Cash-Karp-6-4-5 
^^^^^^^^^^^^^^^^^^

.. index:: Cash-Karp-6-4-5 ERK method

Accessible via the constant ``CASH_KARP_6_4_5`` to
:c:func:`ARKStepSetERKTableNum()`, :c:func:`ERKStepSetERKTableNum()`
or :c:func:`ARKodeLoadButcherTable()`.  This is the default 5th order 
explicit method (from [CK1990]_). 

.. math::

   \begin{array}{r|cccccc}
     0 & 0 & 0 & 0 & 0 & 0 & 0 \\
     \frac{1}{5} & \frac{1}{5} & 0 & 0 & 0 & 0 & 0 \\
     \frac{3}{10} & \frac{3}{40} & \frac{9}{40} & 0 & 0 & 0 & 0 \\
     \frac{3}{5} & \frac{3}{10} & -\frac{9}{10} & \frac{6}{5} & 0 & 0 & 0 \\
     1 & -\frac{11}{54} & \frac{5}{2} & -\frac{70}{27} & \frac{35}{27} & 0 & 0 \\
     \frac{7}{8} & \frac{1631}{55296} & \frac{175}{512} & \frac{575}{13824} & \frac{44275}{110592} & \frac{253}{4096} & 0 \\
     \hline
     5 & \frac{37}{378} & 0 & \frac{250}{621} & \frac{125}{594} & 0 & \frac{512}{1771} \\
     4 & \frac{2825}{27648} & 0 & \frac{18575}{48384} & \frac{13525}{55296} & \frac{277}{14336} & \frac{1}{4}
   \end{array}

.. figure:: figs/stab_region_6.png
   :scale: 50 %
   :align: center
   
   Linear stability region for the Cash-Karp method.  The method's
   region is outlined in blue; the embedding's region is in red.







.. _Butcher.Fehlberg:

Fehlberg-6-4-5 
^^^^^^^^^^^^^^^^^

.. index:: Fehlberg-6-4-5 ERK method

Accessible via the constant ``FEHLBERG_6_4_5`` to
:c:func:`ARKStepSetERKTableNum()`, :c:func:`ERKStepSetERKTableNum()`
or :c:func:`ARKodeLoadButcherTable()` (from [F1969]_). 

.. math::

   \begin{array}{r|cccccc}
     0 & 0 & 0 & 0 & 0 & 0 & 0 \\
     \frac{1}{4} & \frac{1}{4} & 0 & 0 & 0 & 0 & 0 \\
     \frac{3}{8} & \frac{3}{32} & \frac{9}{32} & 0 & 0 & 0 & 0 \\
     \frac{12}{13} & \frac{1932}{2197} & -\frac{7200}{2197} & \frac{7296}{2197} & 0 & 0 & 0 \\
     1 & \frac{439}{216} & -8 & \frac{3680}{513} & -\frac{845}{4104} & 0 & 0 \\
     \frac{1}{2} & -\frac{8}{27} & 2 & -\frac{3544}{2565} & \frac{1859}{4104} & -\frac{11}{40} & 0 \\
     \hline
     5 & \frac{16}{135} & 0 & \frac{6656}{12825} & \frac{28561}{56430} & -\frac{9}{50} & \frac{2}{55} \\
     4 & \frac{25}{216} & 0 & \frac{1408}{2565} & \frac{2197}{4104} & -\frac{1}{5} & 0
   \end{array}

.. figure:: figs/stab_region_7.png
   :scale: 50 %
   :align: center
   
   Linear stability region for the Fehlberg method.  The method's
   region is outlined in blue; the embedding's region is in red.





.. _Butcher.Dormand_Prince:

Dormand-Prince-7-4-5 
^^^^^^^^^^^^^^^^^^^^^^^

.. index:: Dormand-Prince-7-4-5 ERK method

Accessible via the constant ``DORMAND_PRINCE_7_4_5`` to
:c:func:`ARKStepSetERKTableNum()`, :c:func:`ERKStepSetERKTableNum()`
or :c:func:`ARKodeLoadButcherTable()` (from [DP1980]_).  

.. math::

   \begin{array}{r|ccccccc}
     0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
     \frac{1}{5} & \frac{1}{5} & 0 & 0 & 0 & 0 & 0 & 0 \\
     \frac{3}{10} & \frac{3}{40} & \frac{9}{40} & 0 & 0 & 0 & 0 & 0 \\
     \frac{4}{5} & \frac{44}{45} & -\frac{56}{15} & \frac{32}{9} & 0 & 0 & 0 & 0 \\
     \frac{8}{9} & \frac{19372}{6561} & -\frac{25360}{2187} & \frac{64448}{6561} & -\frac{212}{729} & 0 & 0 & 0 \\
     1 & \frac{9017}{3168} & -\frac{355}{33} & \frac{46732}{5247} & \frac{49}{176} & -\frac{5103}{18656} & 0 & 0 \\
     1 & \frac{35}{384} & 0 & \frac{500}{1113} & \frac{125}{192} & -\frac{2187}{6784} & \frac{11}{84} & 0 \\
     \hline
     5 & \frac{35}{384} & 0 & \frac{500}{1113} & \frac{125}{192} & -\frac{2187}{6784} & \frac{11}{84} & 0 \\
     4 & \frac{5179}{57600} & 0 & \frac{7571}{16695} & \frac{393}{640} & -\frac{92097}{339200} & \frac{187}{2100} & \frac{1}{40}
   \end{array}

.. figure:: figs/stab_region_8.png
   :scale: 50 %
   :align: center
   
   Linear stability region for the Dormand-Prince method.  The method's
   region is outlined in blue; the embedding's region is in red.





.. _Butcher.ARK_8_4_5_E:

ARK-8-4-5 (explicit) 
^^^^^^^^^^^^^^^^^^^^^^

.. index:: ARK-8-4-5 ERK method

Accessible via the constant ``ARK548L2SA_ERK_8_4_5`` to
:c:func:`ARKStepSetERKTableNum()`, :c:func:`ERKStepSetERKTableNum()`
or :c:func:`ARKodeLoadButcherTable()`.  This is the explicit portion
of the default 5th order additive method (from [KC2003]_). 

.. math::

   \begin{array}{r|cccccccc}
     0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
     \frac{41}{100} & \frac{41}{100} & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
     \frac{2935347310677}{11292855782101} & \frac{367902744464}{2072280473677} & \frac{677623207551}{8224143866563} & 0 & 0 & 0 & 0 & 0 & 0 \\
     \frac{1426016391358}{7196633302097} & \frac{1268023523408}{10340822734521} & 0 & \frac{1029933939417}{13636558850479} & 0 & 0 & 0 & 0 & 0 \\
     \frac{92}{100} & \frac{14463281900351}{6315353703477} & 0 & \frac{66114435211212}{5879490589093} & -\frac{54053170152839}{4284798021562} & 0 & 0 & 0 & 0 \\
     \frac{24}{100} & \frac{14090043504691}{34967701212078} & 0 & \frac{15191511035443}{11219624916014} & -\frac{18461159152457}{12425892160975} & -\frac{281667163811}{9011619295870} & 0 & 0 & 0 \\
     \frac{3}{5} & \frac{19230459214898}{13134317526959} & 0 & \frac{21275331358303}{2942455364971} & -\frac{38145345988419}{4862620318723} & -\frac{1}{8} & -\frac{1}{8} & 0 & 0 \\
     1 & -\frac{19977161125411}{11928030595625} & 0 & -\frac{40795976796054}{6384907823539} & \frac{177454434618887}{12078138498510} & \frac{782672205425}{8267701900261} & -\frac{69563011059811}{9646580694205} & \frac{7356628210526}{4942186776405} & 0 \\
     \hline
     5 & -\frac{872700587467}{9133579230613} & 0 & 0 & \frac{22348218063261}{9555858737531} & -\frac{1143369518992}{8141816002931} & -\frac{39379526789629}{19018526304540} & \frac{32727382324388}{42900044865799} & \frac{41}{200} \\
     4 & -\frac{975461918565}{9796059967033} & 0 & 0 & \frac{78070527104295}{32432590147079} & -\frac{548382580838}{3424219808633} & -\frac{33438840321285}{15594753105479} & \frac{3629800801594}{4656183773603} & \frac{4035322873751}{18575991585200}
   \end{array}

.. figure:: figs/stab_region_9.png
   :scale: 50 %
   :align: center
   
   Linear stability region for the explicit ARK-8-4-5 method.  The method's
   region is outlined in blue; the embedding's region is in red.






.. _Butcher.Verner-6-5:

Verner-8-5-6 
^^^^^^^^^^^^^^

.. index:: Verner-8-5-6 ERK method

Accessible via the constant ``VERNER_8_5_6`` to
:c:func:`ARKStepSetERKTableNum()`, :c:func:`ERKStepSetERKTableNum()`
or :c:func:`ARKodeLoadButcherTable()`.  This is the default 6th order 
explicit method (from [V1978]_).

.. math::

   \begin{array}{r|cccccccc}
     0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
     \frac{1}{6} & \frac{1}{6} & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
     \frac{4}{15} & \frac{4}{75} & \frac{16}{75} & 0 & 0 & 0 & 0 & 0 & 0 \\
     \frac{2}{3} & \frac{5}{6} & -\frac{8}{3} & \frac{5}{2} & 0 & 0 & 0 & 0 & 0 \\
     \frac{5}{6} & -\frac{165}{64} & \frac{55}{6} & -\frac{425}{64} & \frac{85}{96} & 0 & 0 & 0 & 0 \\
     1 & \frac{12}{5} & -8 & \frac{4015}{612} & -\frac{11}{36} & \frac{88}{255} & 0 & 0 & 0 \\
     \frac{1}{15} & -\frac{8263}{15000} & \frac{124}{75} & -\frac{643}{680} & -\frac{81}{250} & \frac{2484}{10625} & 0 & 0 & 0 \\
     1 & \frac{3501}{1720} & -\frac{300}{43} & \frac{297275}{52632} & -\frac{319}{2322} & \frac{24068}{84065} & 0 & \frac{3850}{26703} & 0 \\
     \hline
     6 & \frac{3}{40} & 0 & \frac{875}{2244} & \frac{23}{72} & \frac{264}{1955} & 0 & \frac{125}{11592} & \frac{43}{616} \\
     5 & \frac{13}{160} & 0 & \frac{2375}{5984} & \frac{5}{16} & \frac{12}{85} & \frac{3}{44} & 0 & 0
   \end{array}

.. figure:: figs/stab_region_10.png
   :scale: 50 %
   :align: center
   
   Linear stability region for the Verner-8-5-6 method.  The method's
   region is outlined in blue; the embedding's region is in red.



.. _Butcher.Fehlberg-8-7:

Fehlberg-13-7-8 
^^^^^^^^^^^^^^^^^^

.. index:: Fehlberg-13-7-8 ERK method

Accessible via the constant ``FEHLBERG_13_7_8`` to
:c:func:`ARKStepSetERKTableNum()`, :c:func:`ERKStepSetERKTableNum()`
or :c:func:`ARKodeLoadButcherTable()`.  This is the default 8th order 
explicit method (from [B2008]_). 

.. math::

   \begin{array}{r|ccccccccccccc}
     0&   0& 0& 0& 0& 0& 0& 0& 0& 0& 0& 0& 0& 0\\
     \frac{2}{27}&   \frac{2}{27}& 0& 0& 0& 0& 0& 0& 0& 0& 0& 0& 0& 0\\
     \frac{1}{9}&   \frac{1}{36}& \frac{1}{12}& 0& 0& 0& 0& 0& 0& 0& 0& 0& 0& 0\\
     \frac{1}{6}&   \frac{1}{24}& 0& \frac{1}{8}& 0& 0& 0& 0& 0& 0& 0& 0& 0& 0\\
     \frac{5}{12}&   \frac{5}{12}& 0& -\frac{25}{16}& \frac{25}{16}& 0& 0& 0& 0& 0& 0& 0& 0& 0\\
     \frac{1}{2}&   \frac{1}{20}& 0& 0& \frac{1}{4}& \frac{1}{5}& 0& 0& 0& 0& 0& 0& 0& 0\\
     \frac{5}{6}&   -\frac{25}{108}& 0& 0& \frac{125}{108}& -\frac{65}{27}& \frac{125}{54}& 0& 0& 0& 0& 0& 0& 0\\
     \frac{1}{6}&   \frac{31}{300}& 0& 0& 0& \frac{61}{225}& -\frac{2}{9}& \frac{13}{900}& 0& 0& 0& 0& 0& 0\\
     \frac{2}{3}&   2& 0& 0& -\frac{53}{6}& \frac{704}{45}& -\frac{107}{9}& \frac{67}{90}& 3& 0& 0& 0& 0& 0\\
     \frac{1}{3}&   -\frac{91}{108}& 0& 0& \frac{23}{108}& -\frac{976}{135}& \frac{311}{54}& -\frac{19}{60}& \frac{17}{6}& -\frac{1}{12}& 0& 0& 0& 0\\
     1&   \frac{2383}{4100}& 0& 0& -\frac{341}{164}& \frac{4496}{1025}& -\frac{301}{82}& \frac{2133}{4100}& \frac{45}{82}& \frac{45}{164}& \frac{18}{41}& 0& 0& 0\\
     0&   \frac{3}{205}& 0& 0& 0& 0& -\frac{6}{41}& -\frac{3}{205}& -\frac{3}{41}& \frac{3}{41}& \frac{6}{41}& 0& 0& 0\\
     1&   -\frac{1777}{4100}& 0& 0& -\frac{341}{164}& \frac{4496}{1025}& -\frac{289}{82}& \frac{2193}{4100}& \frac{51}{82}& \frac{33}{164}& \frac{12}{41}& 0& 1& 0\\
     \hline
     8& 0& 0& 0& 0& 0& \frac{34}{105}& \frac{9}{35}& \frac{9}{35}& \frac{9}{280}& \frac{9}{280}& 0& \frac{41}{840}& \frac{41}{840} \\
     7& \frac{41}{840}& 0& 0& 0& 0& \frac{34}{105}& \frac{9}{35}& \frac{9}{35}& \frac{9}{280}& \frac{9}{280}& \frac{41}{840}& 0& 0
   \end{array}


.. figure:: figs/stab_region_23.png
   :scale: 50 %
   :align: center
   
   Linear stability region for the Fehlberg-13-7-8 method.  The method's
   region is outlined in blue; the embedding's region is in red.





.. _Butcher.implicit:

Implicit Butcher tables
---------------------------


In the category of diagonally implicit Runge-Kutta methods, ARKode
includes methods that have orders 2 through 5, with embeddings that are of
orders 1 through 4.



.. _Butcher.SDIRK-2-1:

SDIRK-2-1-2 
^^^^^^^^^^^^^^

.. index:: SDIRK-2-1-2 method

Accessible via the constant ``SDIRK_2_1_2`` to
:c:func:`ARKStepSetIRKTableNum()` or
:c:func:`ARKodeLoadButcherTable()`.  This is the default 2nd order 
implicit method.  Both the method and embedding are A- and B-stable.

.. math::

   \begin{array}{r|cc}
     1 & 1 & 0 \\
     0 & -1 & 1 \\
     \hline
     2 & \frac{1}{2} & \frac{1}{2} \\
     1 & 1 & 0
   \end{array}

.. figure:: figs/stab_region_11.png
   :scale: 50 %
   :align: center
   
   Linear stability region for the SDIRK-2-1-2 method.  The method's
   region is outlined in blue; the embedding's region is in red.






.. _Butcher.Billington:

Billington-3-3-2 
^^^^^^^^^^^^^^^^^^^

.. index:: Billington-3-3-2 SDIRK method

Accessible via the constant ``BILLINGTON_3_3_2`` to
:c:func:`ARKStepSetIRKTableNum()` or
:c:func:`ARKodeLoadButcherTable()`.  Here, the higher-order
embedding is less stable than the lower-order method (from [B1983]_).

.. math::

   \begin{array}{r|ccc}
     0.292893218813 & 0.292893218813 & 0 & 0 \\
     1.091883092037 & 0.798989873223 & 0.292893218813 & 0 \\
     1.292893218813 & 0.740789228841 & 0.259210771159 & 0.292893218813 \\
     \hline
     2 & 0.740789228840 & 0.259210771159 & 0 \\
     3 & 0.691665115992 & 0.503597029883 & -0.195262145876
   \end{array}

.. figure:: figs/stab_region_12.png
   :scale: 50 %
   :align: center
   
   Linear stability region for the Billington method.  The method's
   region is outlined in blue; the embedding's region is in red.






.. _Butcher.TRBDF2:

TRBDF2-3-3-2 
^^^^^^^^^^^^^^^

.. index:: TRBDF2-3-3-2 ESDIRK method

Accessible via the constant ``TRBDF2_3_3_2`` to
:c:func:`ARKStepSetIRKTableNum()` or
:c:func:`ARKodeLoadButcherTable()`.  As with Billington, here the 
higher-order embedding is less stable than the lower-order method
(from [B1985]_). 

.. math::

   \begin{array}{r|ccc}
     0 & 0 & 0 & 0 \\
     2-\sqrt{2} & \frac{2-\sqrt{2}}{2} & \frac{2-\sqrt{2}}{2} & 0 \\
     1 & \frac{\sqrt{2}}{4} & \frac{\sqrt{2}}{4} & \frac{2-\sqrt{2}}{2} \\
     \hline
     2 & \frac{\sqrt{2}}{4} & \frac{\sqrt{2}}{4} & \frac{2-\sqrt{2}}{2} \\
     3 & \frac{1-\frac{\sqrt{2}}{4}}{3} & \frac{\frac{3\sqrt{2}}{4}+1}{3} & \frac{2-\sqrt{2}}{6}
   \end{array}

.. figure:: figs/stab_region_13.png
   :scale: 50 %
   :align: center
   
   Linear stability region for the TRBDF2 method.  The method's
   region is outlined in blue; the embedding's region is in red.






.. _Butcher.Kvaerno_4_2_3:

Kvaerno-4-2-3 
^^^^^^^^^^^^^^^^

.. index:: Kvaerno-4-2-3 ESDIRK method

Accessible via the constant ``KVAERNO_4_2_3`` to
:c:func:`ARKStepSetIRKTableNum()` or
:c:func:`ARKodeLoadButcherTable()`.  Both the method and embedding are 
A-stable; additionally the method is L-stable (from [K2004]_).

.. math::

   \begin{array}{r|cccc}
     0 & 0 & 0 & 0 & 0 \\
     0.871733043 & 0.4358665215 & 0.4358665215 & 0 & 0 \\
     1 & 0.490563388419108 & 0.073570090080892 & 0.4358665215 & 0 \\
     1 & 0.308809969973036 & 1.490563388254106 & -1.235239879727145 & 0.4358665215 \\
     \hline
     3 & 0.308809969973036 & 1.490563388254106 & -1.235239879727145 & 0.4358665215 \\
     2 & 0.490563388419108 & 0.073570090080892 & 0.4358665215 & 0
   \end{array}

.. figure:: figs/stab_region_14.png
   :scale: 50 %
   :align: center
   
   Linear stability region for the Kvaerno-4-2-3 method.  The method's
   region is outlined in blue; the embedding's region is in red.






.. _Butcher.ARK_4_2_3_I:

ARK-4-2-3 (implicit) 
^^^^^^^^^^^^^^^^^^^^^^^

.. index:: ARK-4-2-3 ESDIRK method

Accessible via the constant ``ARK324L2SA_DIRK_4_2_3`` to
:c:func:`ARKStepSetIRKTableNum()` or
:c:func:`ARKodeLoadButcherTable()`.  This is the default 3rd order 
implicit method, and the implicit portion of the default 3rd order
additive method.  Both the method and embedding are A-stable;
additionally the method is L-stable (from [KC2003]_). 

.. math::

   \begin{array}{r|cccc}
     0 & 0 & 0 & 0 & 0 \\
     \frac{1767732205903}{2027836641118} & \frac{1767732205903}{4055673282236} & \frac{1767732205903}{4055673282236} & 0 & 0 \\
     \frac{3}{5} & \frac{2746238789719}{10658868560708} & -\frac{640167445237}{6845629431997} & \frac{1767732205903}{4055673282236} & 0 \\
     1 & \frac{1471266399579}{7840856788654} & -\frac{4482444167858}{7529755066697} & \frac{11266239266428}{11593286722821} & \frac{1767732205903}{4055673282236} \\
     \hline
     3 & \frac{1471266399579}{7840856788654} & -\frac{4482444167858}{7529755066697} & \frac{11266239266428}{11593286722821} & \frac{1767732205903}{4055673282236} \\
     2 & \frac{2756255671327}{12835298489170} & -\frac{10771552573575}{22201958757719} & \frac{9247589265047}{10645013368117} & \frac{2193209047091}{5459859503100}
   \end{array}

.. figure:: figs/stab_region_15.png
   :scale: 50 %
   :align: center
   
   Linear stability region for the implicit ARK-4-2-3 method.  The method's
   region is outlined in blue; the embedding's region is in red.







.. _Butcher.Cash_5_2_4:

Cash-5-2-4 
^^^^^^^^^^^^^^

.. index:: Cash-5-2-4 SDIRK method

Accessible via the constant ``CASH_5_2_4`` to
:c:func:`ARKStepSetIRKTableNum()` or
:c:func:`ARKodeLoadButcherTable()`.  Both the method and embedding are
A-stable; additionally the method is L-stable (from [C1979]_).  

.. math::

   \begin{array}{r|ccccc}
     0.435866521508 & 0.435866521508 & 0 & 0 & 0 & 0 \\
     -0.7 & -1.13586652150 & 0.435866521508 & 0 & 0 & 0 \\
     0.8 & 1.08543330679 & -0.721299828287 & 0.435866521508 & 0 & 0 \\
     0.924556761814 & 0.416349501547 & 0.190984004184 & -0.118643265417 & 0.435866521508 & 0 \\
     1 & 0.896869652944 & 0.0182725272734 & -0.0845900310706 & -0.266418670647 & 0.435866521508 \\
     \hline
     4 & 0.896869652944 & 0.0182725272734 & -0.0845900310706 & -0.266418670647 & 0.435866521508 \\
     2 & 1.05646216107052 & -0.0564621610705236 & 0 & 0 & 0
   \end{array}

.. figure:: figs/stab_region_16.png
   :scale: 50 %
   :align: center
   
   Linear stability region for the Cash-5-2-4 method.  The method's
   region is outlined in blue; the embedding's region is in red.






.. _Butcher.Cash_5_3_4:

Cash-5-3-4
^^^^^^^^^^^

.. index:: Cash-5-3-4 SDIRK method

Accessible via the constant ``CASH_5_3_4`` to
:c:func:`ARKStepSetIRKTableNum()` or
:c:func:`ARKodeLoadButcherTable()`.  Both the method and embedding are
A-stable; additionally the method is L-stable (from [C1979]_).

.. math::

   \begin{array}{r|ccccc}
     0.435866521508 & 0.435866521508 & 0 & 0 & 0 & 0 \\
     -0.7 & -1.13586652150 & 0.435866521508 & 0 & 0 & 0 \\
     0.8 & 1.08543330679 & -0.721299828287 & 0.435866521508 & 0 & 0 \\
     0.924556761814 & 0.416349501547 & 0.190984004184 & -0.118643265417 & 0.435866521508 & 0 \\
     1 & 0.896869652944 & 0.0182725272734 & -0.0845900310706 & -0.266418670647 & 0.435866521508 \\
     \hline
     4 & 0.896869652944 & 0.0182725272734 & -0.0845900310706 & -0.266418670647 & 0.435866521508 \\
     3 & 0.776691932910 & 0.0297472791484 & -0.0267440239074 & 0.220304811849 & 0
   \end{array}

.. figure:: figs/stab_region_17.png
   :scale: 50 %
   :align: center
   
   Linear stability region for the Cash-5-3-4 method.  The method's
   region is outlined in blue; the embedding's region is in red.





.. _Butcher.SDIRK-5-4:

SDIRK-5-3-4 
^^^^^^^^^^^^^^

.. index:: SDIRK-5-3-4 method

Accessible via the constant ``SDIRK_5_3_4`` to
:c:func:`ARKStepSetIRKTableNum()` or
:c:func:`ARKodeLoadButcherTable()`.  This is the default 4th order 
implicit method.  Here, the method is both A- and L-stable, although
the embedding has reduced stability (from [HW1996]_). 

.. math::

   \begin{array}{r|ccccc}
     \frac{1}{4} & \frac{1}{4} & 0 & 0 & 0 & 0 \\
     \frac{3}{4} & \frac{1}{2} & \frac{1}{4} & 0 & 0 & 0 \\
     \frac{11}{20} & \frac{17}{50} & -\frac{1}{25} & \frac{1}{4} & 0 & 0 \\
     \frac{1}{2} & \frac{371}{1360} & -\frac{137}{2720} & \frac{15}{544} & \frac{1}{4} & 0 \\
     1 & \frac{25}{24} & -\frac{49}{48} & \frac{125}{16} & -\frac{85}{12} & \frac{1}{4} \\
     \hline
     4 & \frac{25}{24} & -\frac{49}{48} & \frac{125}{16} & -\frac{85}{12} & \frac{1}{4} \\
     3 & \frac{59}{48} & -\frac{17}{96} & \frac{225}{32} & -\frac{85}{12} & 0
   \end{array}

.. figure:: figs/stab_region_18.png
   :scale: 50 %
   :align: center
   
   Linear stability region for the SDIRK-5-3-4 method.  The method's
   region is outlined in blue; the embedding's region is in red.








.. _Butcher.Kvaerno_5_3_4:

Kvaerno-5-3-4 
^^^^^^^^^^^^^^^^

.. index:: Kvaerno-5-3-4 ESDIRK method

Accessible via the constant ``KVAERNO_5_3_4`` to
:c:func:`ARKStepSetIRKTableNum()` or
:c:func:`ARKodeLoadButcherTable()`.  Both the method and embedding are 
A-stable (from [K2004]_).

.. math::

   \begin{array}{r|ccccc}
     0 & 0 & 0 & 0 & 0 & 0 \\
     0.871733043 & 0.4358665215  & 0.4358665215  & 0 & 0 & 0 \\
     0.468238744853136 & 0.140737774731968 & -0.108365551378832 & 0.4358665215 & 0 & 0 \\
     1 & 0.102399400616089 & -0.376878452267324 & 0.838612530151233 & 0.4358665215 & 0 \\
     1 & 0.157024897860995 & 0.117330441357768 & 0.61667803039168 & -0.326899891110444 & 0.4358665215 \\
     \hline
     4 & 0.157024897860995 & 0.117330441357768 & 0.61667803039168 & -0.326899891110444 & 0.4358665215 \\
     3 & 0.102399400616089 & -0.376878452267324 & 0.838612530151233 & 0.4358665215 & 0
   \end{array}

.. figure:: figs/stab_region_19.png
   :scale: 50 %
   :align: center
   
   Linear stability region for the Kvaerno-5-3-4 method.  The method's
   region is outlined in blue; the embedding's region is in red.






.. _Butcher.ARK_6_3_4_I:

ARK-6-3-4 (implicit) 
^^^^^^^^^^^^^^^^^^^^^^^

.. index:: ARK-6-3-4 ESDIRK method

Accessible via the constant ``ARK436L2SA_DIRK_6_3_4`` to
:c:func:`ARKStepSetIRKTableNum()` or
:c:func:`ARKodeLoadButcherTable()`.  This is the implicit portion
of the default 4th order additive method. Both the method and
embedding are A-stable; additionally the method is L-stable (from [KC2003]_). 

.. math::

   \begin{array}{r|cccccc}
     0 & 0 & 0 & 0 & 0 & 0 & 0 \\
     \frac{1}{2} & \frac{1}{4} & \frac{1}{4} & 0 & 0 & 0 & 0 \\
     \frac{83}{250} & \frac{8611}{62500} & -\frac{1743}{31250} & \frac{1}{4} & 0 & 0 & 0 \\
     \frac{31}{50} & \frac{5012029}{34652500} & -\frac{654441}{2922500} & \frac{174375}{388108} & \frac{1}{4} & 0 & 0 \\
     \frac{17}{20} & \frac{15267082809}{155376265600} & -\frac{71443401}{120774400} & \frac{730878875}{902184768} & \frac{2285395}{8070912} & \frac{1}{4} & 0 \\
     1 & \frac{82889}{524892} & 0 & \frac{15625}{83664} & \frac{69875}{102672} & -\frac{2260}{8211} & \frac{1}{4} \\
     \hline
     4 & \frac{82889}{524892} & 0 & \frac{15625}{83664} & \frac{69875}{102672} & -\frac{2260}{8211} & \frac{1}{4} \\
     3 & \frac{4586570599}{29645900160} & 0 & \frac{178811875}{945068544} & \frac{814220225}{1159782912} & -\frac{3700637}{11593932} & \frac{61727}{225920}
   \end{array}

.. figure:: figs/stab_region_20.png
   :scale: 50 %
   :align: center
   
   Linear stability region for the implicit ARK-6-3-4 method.  The method's
   region is outlined in blue; the embedding's region is in red.






.. _Butcher.Kvaerno_7_4_5:

Kvaerno-7-4-5 
^^^^^^^^^^^^^^^^^

.. index:: Kvaerno-7-4-5 ESDIRK method

Accessible via the constant ``KVAERNO_7_4_5`` to
:c:func:`ARKStepSetIRKTableNum()` or
:c:func:`ARKodeLoadButcherTable()`.  Both the method and embedding are 
A-stable; additionally the method is L-stable (from [K2004]_).

.. math::

   \begin{array}{r|ccccccc}
     0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
     0.52 & 0.26 & 0.26 & 0 & 0 & 0 & 0 & 0 \\
     1.230333209967908 & 0.13 & 0.84033320996790809 & 0.26 & 0 & 0 & 0 & 0 \\
     0.895765984350076 & 0.22371961478320505 & 0.47675532319799699 & -0.06470895363112615 & 0.26 & 0 & 0 & 0 \\
     0.436393609858648 & 0.16648564323248321 & 0.10450018841591720 & 0.03631482272098715 & -0.13090704451073998 & 0.26 & 0 & 0 \\
     1 & 0.13855640231268224 & 0 & -0.04245337201752043 & 0.02446657898003141 & 0.61943039072480676 & 0.26 & 0 \\
     1 & 0.13659751177640291 & 0 & -0.05496908796538376 & -0.04118626728321046 & 0.62993304899016403 & 0.06962479448202728 & 0.26 \\
     \hline
     5 & 0.13659751177640291 & 0 & -0.05496908796538376 & -0.04118626728321046 & 0.62993304899016403 & 0.06962479448202728 & 0.26 \\
     4 & 0.13855640231268224 & 0 & -0.04245337201752043 & 0.02446657898003141 & 0.61943039072480676 & 0.26 & 0
   \end{array}

.. figure:: figs/stab_region_21.png
   :scale: 50 %
   :align: center
   
   Linear stability region for the Kvaerno-7-4-5 method.  The method's
   region is outlined in blue; the embedding's region is in red.







.. _Butcher.ARK_8_4_5_I:

ARK-8-4-5 (implicit)
^^^^^^^^^^^^^^^^^^^^^^

.. index:: ARK-8-4-5 ESDIRK method

Accessible via the constant ``ARK548L2SA_DIRK_8_4_5`` for
:c:func:`ARKStepSetIRKTableNum()` or
:c:func:`ARKodeLoadButcherTable()`.  This is the default 5th order 
implicit method, and the implicit portion of the default 5th order
additive method.  Both the method and embedding are A-stable;
additionally the method is L-stable (from [KC2003]_).

.. math::

   \begin{array}{r|cccccccc}
     0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
     \frac{41}{100} & \frac{41}{200} & \frac{41}{200} & 0 & 0 & 0 & 0 & 0 & 0 \\
     \frac{2935347310677}{11292855782101} & \frac{41}{400} & -\frac{567603406766}{11931857230679} & \frac{41}{200} & 0 & 0 & 0 & 0 & 0 \\
     \frac{1426016391358}{7196633302097} & \frac{683785636431}{9252920307686} & 0 & -\frac{110385047103}{1367015193373} & \frac{41}{200} & 0 & 0 & 0 & 0 \\
     \frac{92}{100} & \frac{3016520224154}{10081342136671} & 0 & \frac{30586259806659}{12414158314087} & -\frac{22760509404356}{11113319521817} & \frac{41}{200} & 0 & 0 & 0 \\
     \frac{24}{100} & \frac{218866479029}{1489978393911} & 0 & \frac{638256894668}{5436446318841} & -\frac{1179710474555}{5321154724896} & -\frac{60928119172}{8023461067671} & \frac{41}{200} & 0 & 0 \\
     \frac{3}{5} & \frac{1020004230633}{5715676835656} & 0 & \frac{25762820946817}{25263940353407} & -\frac{2161375909145}{9755907335909} & -\frac{211217309593}{5846859502534} & -\frac{4269925059573}{7827059040749} & \frac{41}{200} & 0 \\
     1 & -\frac{872700587467}{9133579230613} & 0 & 0 & \frac{22348218063261}{9555858737531} & -\frac{1143369518992}{8141816002931} & -\frac{39379526789629}{19018526304540} & \frac{32727382324388}{42900044865799} & \frac{41}{200} \\
     \hline
     5 & -\frac{872700587467}{9133579230613} & 0 & 0 & \frac{22348218063261}{9555858737531} & -\frac{1143369518992}{8141816002931} & -\frac{39379526789629}{19018526304540} & \frac{32727382324388}{42900044865799} & \frac{41}{200} \\
     4 & -\frac{975461918565}{9796059967033} & 0 & 0 & \frac{78070527104295}{32432590147079} & -\frac{548382580838}{3424219808633} & -\frac{33438840321285}{15594753105479} & \frac{3629800801594}{4656183773603} & \frac{4035322873751}{18575991585200}
   \end{array}

.. figure:: figs/stab_region_22.png
   :scale: 50 %
   :align: center
   
   Linear stability region for the implicit ARK-8-4-5 method.  The method's
   region is outlined in blue; the embedding's region is in red.







.. _Butcher.additive:

Additive Butcher tables
---------------------------

In the category of additive Runge-Kutta methods for split implicit and
explicit calculations, ARKode includes methods that have orders 3
through 5, with embeddings that are of orders 2 through 4.  These
Butcher table pairs are as follows:

* :index:`3rd-order pair <ARK-4-2-3 ARK method>`:
  :ref:`Butcher.ARK_4_2_3_E` with :ref:`Butcher.ARK_4_2_3_I`,
  corresponding to Butcher tables ``ARK324L2SA_ERK_4_2_3`` and
  ``ARK324L2SA_DIRK_4_2_3`` for :c:func:`ARKStepSetARKTableNum()`.

* :index:`4th-order pair <ARK-6-3-4 ARK method>`:  
  :ref:`Butcher.ARK_6_3_4_E` with :ref:`Butcher.ARK_6_3_4_I`,
  corresponding to Butcher tables ``ARK436L2SA_ERK_6_3_4`` and
  ``ARK436L2SA_DIRK_6_3_4`` for :c:func:`ARKStepSetARKTableNum()`.  

* :index:`5th-order pair <ARK-8-4-5 ARK method>`:  
  :ref:`Butcher.ARK_8_4_5_E` with :ref:`Butcher.ARK_8_4_5_I`,
  corresponding to Butcher tables ``ARK548L2SA_ERK_8_4_5`` and 
  ``ARK548L2SA_ERK_8_4_5`` for :c:func:`ARKStepSetARKTableNum()`.

