..
   Programmer(s): Cody J. Balos @ LLNL
   ----------------------------------------------------------------
   SUNDIALS Copyright Start
   Copyright (c) 2002-2023, Lawrence Livermore National Security
   and Southern Methodist University.
   All rights reserved.

   See the top-level LICENSE and NOTICE files for details.

   SPDX-License-Identifier: BSD-3-Clause
   SUNDIALS Copyright End
   ----------------------------------------------------------------

######################
SUNDIALS Documentation
######################

This is the documentation for the `SUNDIALS
<https://computing.llnl.gov/projects/sundials>`_ suite of
nonlinear and differential/algebraic equation solvers.

SUNDIALS is developed on `GitHub <https://github.com/LLNL/sundials>`_.

.. toctree::
   :maxdepth: 1
   :numbered:
   :hidden:

   Organization_link.rst
   sundials/index.rst
   arkode/index.rst
   cvode/index.rst
   cvodes/index.rst
   ida/index.rst
   idas/index.rst
   kinsol/index.rst
   nvectors/index.rst
   sunmatrix/index.rst
   sunlinsol/index.rst
   sunnonlinsol/index.rst
   sunadaptcontroller/index.rst
   sunmemory/index.rst
   Install_link.rst
   History_link.rst
   developers/index.rst
   References


Contributors
============

The SUNDIALS library has been developed over many years by a number of
contributors. The current SUNDIALS team consists of Cody J.  Balos,
David J. Gardner, Alan C. Hindmarsh, Daniel R. Reynolds, and Carol S.
Woodward. We thank Radu Serban for significant and critical past contributions.

Other contributors to SUNDIALS include: James Almgren-Bell, Lawrence E. Banks,
Peter N. Brown, George Byrne, Rujeko Chinomona, Scott D. Cohen, Aaron Collier,
Keith E. Grant, Steven L. Lee, Shelby L. Lockhart, John Loffeld, Daniel McGreer,
Yu Pan, Slaven Peles, Cosmin Petra, Steven B. Roberts, H. Hunter Schwartz,
Jean M. Sexton, Dan Shumaker, Steve G. Smith, Shahbaj Sohal, Allan G. Taylor,
Hilari C. Tiedeman, Chris White, Ting Yan, and Ulrike M. Yang.


Citing
======

.. include:: ../../shared/cite_sundials.rst

When using the ARKODE package from SUNDIALS, please also cite:

.. code-block:: latex

   @article{reynolds2023arkode,
     title   = {{ARKODE: A flexible IVP solver infrastructure for one-step methods}},
     author  = {Reynolds, Daniel R and Gardner, David J and Woodward, Carol S and Chinomona, Rujeko},
     journal = {ACM Transactions on Mathematical Software},
     volume  = {49},
     number  = {2},
     pages   = {1--26},
     year    = {2023},
     doi     = {10.1145/3594632}
   }

We also ask that users cite the documentation for the package and version that
they are using rather than the combined SUNDIALS online guide:

.. parsed-literal::

   @Misc{arkodeDocumentation,
      author       = {Daniel R. Reynolds and David J. Gardner and Carol S. Woodward, Rujeko Chinomona and Cody J. Balos},
      title        = {User Documentation for ARKODE},
      year         = {|YEAR|},
      note         = {|ARKODE_VERSION|},
      url          = {https://sundials.readthedocs.io/en/latest/arkode},
      howpublished = {\url{https://sundials.readthedocs.io/en/latest/arkode}}
   }

.. parsed-literal::

   @Misc{cvodeDocumentation,
      author       = {Alan C. Hindmarsh and Radu Serban and Cody J. Balos and David J. Gardner and Daniel R. Reynolds and Carol S. Woodward},
      title        = {User Documentation for CVODE},
      year         = {|YEAR|},
      note         = {|CVODE_VERSION|},
      url          = {https://sundials.readthedocs.io/en/latest/cvode},
      howpublished = {\url{https://sundials.readthedocs.io/en/latest/cvode}}
   }

.. parsed-literal::

   @Misc{cvodesDocumentation,
      author       = {Alan C. Hindmarsh and Radu Serban and Cody J. Balos and David J. Gardner and Daniel R. Reynolds and Carol S. Woodward},
      title        = {User Documentation for CVODES},
      year         = {|YEAR|},
      note         = {|CVODES_VERSION|},
      url          = {https://sundials.readthedocs.io/en/latest/cvodes},
      howpublished = {\url{https://sundials.readthedocs.io/en/latest/cvodes}}
   }

.. parsed-literal::

   @Misc{idaDocumentation,
      author       = {Alan C. Hindmarsh and Radu Serban and Cody J. Balos and David J. Gardner and Daniel R. Reynolds and Carol S. Woodward},
      title        = {User Documentation for IDA},
      year         = {|YEAR|},
      note         = {|IDA_VERSION|},
      url          = {https://sundials.readthedocs.io/en/latest/ida},
      howpublished = {\url{https://sundials.readthedocs.io/en/latest/ida}}
   }

.. parsed-literal::

   @Misc{idasDocumentation,
      author       = {Radu Serban and Cosmin Petra and Alan C. Hindmarsh and Cody J. Balos and David J. Gardner and Daniel R. Reynolds and Carol S. Woodward},
      title        = {User Documentation for IDAS},
      year         = {|YEAR|},
      note         = {|IDAS_VERSION|},
      url          = {https://sundials.readthedocs.io/en/latest/idas},
      howpublished = {\url{https://sundials.readthedocs.io/en/latest/idas}}
   }

.. parsed-literal::

   @Misc{kinsolDocumentation,
      author       = {Alan C. Hindmarsh and Radu Serban and Cody J. Balos and David J. Gardner and Daniel R. Reynolds and Carol S. Woodward},
      title        = {User Documentation for KINSOL},
      year         = {|YEAR|},
      note         = {|KINSOL_VERSION|},
      url          = {https://sundials.readthedocs.io/en/latest/kinsol},
      howpublished = {\url{https://sundials.readthedocs.io/en/latest/kinsol}}
   }


SUNDIALS License and Notices
============================

.. include:: ../../shared/LicenseReleaseNumbers.rst
