---
name: sundials-guide
description: Guide SUNDIALS users from problem classification to a defensible package, time-stepping module, nonlinear/linear solver strategy, initial settings, and (when requested) tuning-ready code and a constrained suntools tune workflow. Use when a request asks which SUNDIALS package to use (CVODE, CVODES, IDA, IDAS, ARKODE, KINSOL), whether to choose explicit vs implicit vs IMEX or multirate methods, how to pick dense/band/sparse/Krylov solvers and preconditioners, what tolerances or step limits to try, or how to autotune a selected solver stack.
---

# Choose and Tune a SUNDIALS Strategy

## Overview

Map a user problem to the right SUNDIALS package and give a defensible first-pass configuration. Keep the recommendation concrete: name the package, method family, linear/nonlinear solver approach, starting settings, and the next doc/example file to inspect.

Treat this as a two-stage workflow when the user asks for code, performance work, or autotuning: first select or narrow down a numerically defensible solver stack, then measure and tune that stack with `suntools tune`. A selection-only request gets a concise tuning handoff; an end-to-end request gets the baseline code pattern, tuning configuration/command, and validation procedure. Tuning results are empirical and specific to the problem workload, input, machine, SUNDIALS build, vector/backend configuration, and search space; describe a result as the best observed configuration, not a universally best solver.

Start by extracting only the facts that actually control package and solver choice. If key facts are missing, ask for the minimum needed to classify the problem: equation type, stiffness, sensitivities, natural splitting, system size, sparsity, and any special structure. Do not treat the `CVODE` versus `ARKODE` choice as purely a question of model structure; stiffness, numerical stability and damping requirements can decide it even for an unsplit ODE.

## Intake Checklist

Collect these before recommending a package or method:

- problem form: nonlinear algebraic system, ODE IVP, or DAE IVP
- need for forward or adjoint sensitivities
- need for quadrature integration, with or without sensitivities
- expected stiffness: nonstiff, stiff, or mixed stiff/nonstiff
- natural model split: explicit/implicit, slow/fast, Hamiltonian, operator splitting
- whether one-step RK stability properties matter: e.g., need for A-stability or L-stability at order above 2, strong damping of fast modes, or concern about higher-order BDF stability-angle limits
- whether one-step flexibility matters: e.g., if a spatial mesh will be changed between time steps
- size and matrix structure: small dense, banded, sparse, or very large matrix-free
- Jacobian and preconditioner availability
- important constraints: positivity, algebraic variables, event/rootfinding, fixed output cadence, structure preservation
- target hardware and vector/backend constraints if they matter to the recommendation
- for code or tuning requests: implementation language, representative workload, performance objective, accuracy/correctness metric and bound, tuning budget, and target machine/build
- whether structural alternatives (package, stepper, vector/backend, matrix, linear solver, or preconditioner) should be compared as separate executables

If the user does not know whether the model is stiff, infer it from context but label the assumption explicitly.

## Workflow

### 1. Choose the package family

- Use `KINSOL` for nonlinear algebraic systems `F(u) = 0` with no time integration.
- Use `IDA` or `IDAS` for DAEs of the form `F(t, y, y') = 0`.
- Use `CVODE` or `CVODES` for general ODE IVPs when a multistep solver is a good fit and there is no important reason to prefer one-step RK stability or splitting structure.
- Use `ARKODE` when the problem benefits from one-step Runge-Kutta structure or stability properties that its steppers actually support: explicit-only, implicit-only, IMEX/additive splitting, multirate evolution, low-storage explicit stepping, Hamiltonian structure, operator splitting, or a stiff problem where DIRK stability and damping are preferable to higher-order BDF behavior. Do not route users to `ARKODE` for forward sensitivities, and only mention ARKODE adjoints when fixed-step discrete ASA with `ERKStep` or compatible explicit `ARKStep` is acceptable.
- If `CVODE` or `IDA` is the right choice, prefer the supersets `CVODES` and `IDAS` when the user needs forward sensitivities, adjoint sensitivities, or quadrature integration APIs.

Open [package-selection.md](references/package-selection.md) when the package choice is the main question.

### 2. Choose the method or stepper

- For `CVODE` or `CVODES`, recommend `CV_ADAMS` for nonstiff problems and `CV_BDF` for stiff problems when a multistep BDF method is still a good numerical fit.
- For `IDA` or `IDAS`, the main integration method is variable-order BDF; focus the decision on initial-condition handling, Jacobian strategy, and linear solver choice.
- For `ARKODE`, select the narrowest stepper that matches the model:
  - `ERKStep` for fully explicit ODEs.
  - `ARKStep` for fully implicit DIRK or IMEX/additive explicit-implicit splits, especially when one-step A-stable or L-stable behavior is more important than the variable-order BDF workflow.
  - `MRIStep` for genuine slow/fast multirate problems.
  - `SPRKStep` for separable Hamiltonian systems where structure preservation matters.
  - `LSRKStep`, `ForcingStep`, or `SplittingStep` only when the user explicitly benefits from those formulations.
- For `KINSOL`, choose `KIN_LINESEARCH` when robustness matters, `KIN_NONE` for plain Newton with a good initial guess, `KIN_FP` for natural fixed-point maps, and `KIN_PICARD` when a Picard iteration is natural and a linear solver is available.
- When comparing `CV_BDF` and implicit `ARKStep`, remember that BDF methods above order 2 are not A-stable. If the user needs higher-order stiff decay, stronger damping of parasitic fast modes, or a parabolic diffusion problem is behaving poorly under high-order BDF, either cap BDF order at 2 or move the recommendation toward an appropriate DIRK table in `ARKStep`.

Open [method-and-linear-solvers.md](references/method-and-linear-solvers.md) when the package is already known but the method or solver stack is not.

### 3. Choose the linear and nonlinear solver strategy

- Use dense direct solvers for small dense systems.
- Use band direct solvers when the Jacobian bandwidth is small and known.
- Use sparse direct solvers such as KLU when the Jacobian is sparse and factorization cost is still acceptable, and the Jacobian matrix can be constructed directly.
- Use Krylov solvers plus preconditioning for large stiff systems or matrix-free settings.
- Recommend GMRES first when the user needs a generic Krylov choice.
- Recommend FGMRES when the preconditioner changes between iterations.
- Recommend PCG only for symmetric positive definite linear systems.
- Mention BiCGStab or TFQMR when storage is tighter than GMRES or when GMRES restart behavior is a concern.
- If the recommendation depends on a good preconditioner and the user has none, say that clearly instead of overselling Krylov methods.

### 4. Choose initial settings

- Recommend scalar tolerances only when component scales are comparable; otherwise prefer vector absolute tolerances.
- Tie `atol` to the smallest meaningful magnitude per component, not to machine precision by default.
- Tie `rtol` to the requested relative accuracy; if the user has no target, recommend a moderate starting value and tell them to run a convergence study.
- Suggest `IDASetId` and `IDACalcIC` for semi-explicit index-one DAEs when consistent initial conditions are not already available.
- Suggest user Jacobians or Jacobian-vector routines when the problem is sparse, banded, expensive, or noisy under finite differences.
- Suggest `SetMaxStep` only when there is a known physical or event scale that should cap the step size.
- If `CVODE`, `IDA(S)`, or `ARKODE` hits the default internal-step limit, point out that the package default is `500` steps before the next output time and explain whether the fix is larger limits, different tolerances, or a different method.

Open [settings-checklist.md](references/settings-checklist.md) when the user mainly needs tolerances, IC handling, step limits, Jacobian/preconditioner guidance, or failure triage.

### 5. Build a baseline and define the tuning boundary

- For code requests, inspect the closest example under `examples/<package>/` and produce a runnable baseline for the selected stack. The baseline must explicitly construct structural choices: package, stepper, vector/backend, matrix, nonlinear solver, linear solver, and preconditioner as applicable.
- Apply the recommended initial settings and verify the baseline before tuning. Print stable, machine-readable correctness/statistics values (for example, an error or residual and any application metric used by the objective), and return a failure status when the solve fails or the metric cannot be produced.
- Keep structural choices separate from `SetOptions` knobs. `suntools tune` appends supported `SetOptions` key/value pairs to an executable; it can modify exposed properties but does not construct a different package or solver object. Compare structural candidates with separate executables or an explicit application-level selector, using the same workload and measurement protocol.
- Expose only options accepted by the selected package's `SetOptions` interface. Start with a small, defensible set of runtime controls such as compatible method tables, adaptivity parameters, nonlinear/linear iteration controls, or setup frequencies; do not tune every available option by default.
- For C/C++ examples, follow the package `*SetOptions` command-line pattern. For other languages, verify that the binding exposes an equivalent option interface or provide a small argument-parsing bridge before proposing `suntools tune`.

Open [suntools-tuning.md](references/suntools-tuning.md) for the `SetOptions` integration patterns, package-specific knobs, and configuration examples.

### 6. Run a constrained `suntools tune` search

- Use the baseline executable, fixed input, fixed output/verification protocol, and fixed build/environment for every trial. Prefer a YAML configuration for a repeatable workflow; use `--params` for a small exploratory search.
- Choose an objective that reflects the user's goal, usually minimized wall time for a completed solve. Use one worker and at least five repetitions for wall-clock measurements; keep the search budget proportionate to the number of parameters and candidate structural stacks.
- When tuning tolerances, maximum step sizes, method order/table, or other accuracy-affecting controls, add a correctness/error/residual constraint with an explicit bound so a faster but less accurate trial cannot win. Never report an unconstrained timing winner when it changes the requested numerical result. Treat administrative failure limits as validity settings, not performance knobs.
- Tune one structural stack at a time. If several stacks remain plausible, run comparable searches for separate executables and rank them only after applying the same correctness criteria.
- Start with a narrow, physically meaningful search space, check that each parameter is accepted and changes behavior, then expand the budget or search space only if the evidence warrants it. Record failed or infeasible trials rather than silently treating them as wins.

### 7. Validate and report the result

- Compare the tuned result with the default-settings baseline using repeated runs, then rerun the selected configuration on independent representative cases, accuracy checks, and (when relevant) different problem sizes or output intervals.
- Materialize tuned values as application defaults only after validation; otherwise provide them as workload-specific launch/configuration values and preserve the defensible baseline.
- Report the selected structural stack, tunable keys and ranges, objective, correctness constraint, search budget/backend, baseline and best observed measurements, and machine/build details. State clearly when no feasible winner was found or when a result is too noisy to distinguish.
- Re-tune when the workload, tolerances, hardware, compiler/build, vector/backend, or solver structure changes. Do not transfer a measured winner to a different environment without validation.

Open [suntools-tuning.md](references/suntools-tuning.md) for the staged search, correctness constraints, result files, and reproducibility checklist.

## Output Style

When giving a recommendation:

- name the chosen package and why alternatives were rejected
- name the method or stepper and whether the problem is being treated as stiff, nonstiff, mixed, multirate, or structure-preserving
- name the linear/nonlinear solver strategy
- give a short starting settings block
- point to the closest repo docs and examples to inspect next
- if the request is selection-only, add a brief handoff naming which runtime settings could be tested with `suntools tune` and remind the user that measured results are workload/machine specific
- if the request asks for code, performance, autotuning, or an end-to-end workflow, include the runnable baseline pattern, the structural-versus-`SetOptions` boundary, a focused `suntools tune` command or YAML configuration with correctness constraints, and a validation/reporting checklist
- when presenting tuned values, label them as best observed for the stated problem, machine, build, and search space, and include the default baseline for comparison

Prefer wording like:

- "Use `IDAS` because the model is a DAE and you also need sensitivities."
- "Start with `CV_BDF + Newton + GMRES` because the system appears stiff and too large for dense factorization."
- "Use `ARKStep` rather than `CVODE` because the RHS already has a meaningful stiff/nonstiff split."
- "Use `ARKStep` rather than `CVODE` because the problem is stiff and the one-step DIRK stability properties matter, even though there is no natural IMEX split."

Ground recommendations in these repo docs when needed:

- `doc/cvode/guide/source/Introduction.rst`
- `doc/cvode/guide/source/Usage/index.rst`
- `doc/cvodes/guide/source/Introduction.rst`
- `doc/ida/guide/source/Introduction.rst`
- `doc/ida/guide/source/Mathematics.rst`
- `doc/idas/guide/source/Introduction.rst`
- `doc/idas/guide/source/Usage/SIM.rst`
- `doc/arkode/guide/source/Introduction.rst`
- `doc/arkode/guide/source/Usage/index.rst`
- `doc/kinsol/guide/source/Introduction.rst`
- `doc/kinsol/guide/source/Usage/index.rst`

When a user wants code, inspect the closest example under `examples/<package>/` and adapt that pattern instead of describing an abstract setup.

## Pitfalls

Just because a SUNDIALS example uses a method/solver/setting doesn't mean its the right choice. Decisions should be grounded in doc recommendations and published literature on time integrator and solver methods. In particular, do not reduce `CVODE` versus `ARKODE` to "no split" versus "has split": stability region, stiff decay, stage order, and order restrictions also matter.

If the best choice is not something SUNDIALS currently supports, acknowledge this.

Autotuning does not replace numerical judgment. It cannot make an unsuitable package or missing preconditioner suitable, and a timing win obtained by relaxing accuracy, changing the workload, oversubscribing the machine, or accepting failed solves is invalid. Keep solver-stack selection, correctness requirements, and performance measurement explicit throughout the workflow.
