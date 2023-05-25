# CVODE
### Version 6.5.1 (Mar 2023)

**Alan C. Hindmarsh, Radu Serban, Cody J. Balos, David J. Gardner,
  and Carol S. Woodward, Center for Applied Scientific Computing, LLNL**

**Daniel R. Reynolds, Department of Mathematics, Southern Methodist University**

CVODE is a package for the solution of stiff and nonstiff ordinary differential
equation (ODE) systems (initial value problems) given in explicit form
```
dy/dt = f(t,y), y(t0) = y0.
```
CVODE provides a choice of two variable-order, variable-coefficient multistep
methods, Adams-Moulton methods for non-stiff problems or BDF (Backward
Differentiation Formula) methods in fixed-leading-coefficient form for stiff
problems.

CVODE is part of the SUNDIALS Suite of Nonlinear and Differential/Algebraic
equation Solvers which consists of ARKode, CVODE, CVODES, IDA, IDAS and KINSOL.
It is written in ANSI standard C and can be used in a variety of computing
environments including serial, shared memory, distributed memory, and
accelerator-based (e.g., GPU) systems. This flexibility is obtained from a
modular design that leverages the shared vector, matrix, linear solver, and
nonlinear solver APIs used across SUNDIALS packages.

For use with Fortran applications, a set of Fortran/C interface routines, called
FCVODE, is also supplied. These are written in C, but assume that the user
calling program and all user-supplied routines are in Fortran.

## Documentation

See the CVODE documentation at [Read the Docs](https://sundials.readthedocs.io/en/latest/cvode)
for more information about CVODE usage.

## Installation

For installation instructions see the
[SUNDIALS Installation Guide](https://sundials.readthedocs.io/en/latest/Install_link.html).

## Release History

Information on recent changes to CVODE can be found in the "Introduction"
chapter of the CVODE User Guide and a complete release history is available in
the "SUNDIALS Release History" appendix of the CVODE User Guide.

## References

* A. C. Hindmarsh, R. Serban, C. J. Balos, D. J. Gardner, D. R. Reynolds
  and C. S. Woodward, "User Documentation for CVODE v6.5.1,"
  LLNL technical report UCRL-SM-208108, Mar 2023.

* A. C. Hindmarsh and R. Serban, "Example Programs for CVODE v6.5.1,"
  LLNL technical report UCRL-SM-208110, Mar 2023.

* S.D. Cohen and A.C. Hindmarsh, "CVODE, a Stiff/nonstiff ODE Solver in C,"
  Computers in Physics, 10(2), pp. 138-143, 1996.

* A. C. Hindmarsh, P. N. Brown, K. E. Grant, S. L. Lee, R. Serban,
  D. E. Shumaker, and C. S. Woodward, "SUNDIALS, Suite of Nonlinear and
  Differential/Algebraic Equation Solvers," ACM Trans. Math. Softw.,
  31(3), pp. 363-396, 2005.
