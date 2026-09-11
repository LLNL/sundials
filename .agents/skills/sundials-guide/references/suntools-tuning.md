# `suntools tune` handoff

Use this reference after selecting a SUNDIALS package and a plausible solver
stack. Tuning is an experiment around a runnable application; it is not a
replacement for classifying the equations or choosing a viable vector,
matrix, nonlinear solver, linear solver, and preconditioner.

The useful handoff is:

1. establish a baseline with the recommended stack and representative data;
2. expose a small set of supported `SetOptions` parameters;
3. search a focused space with an objective and an accuracy constraint; and
4. validate the best feasible trial before making it a new application default.

The result is the best observed configuration for the measured problem,
machine, SUNDIALS build, and search space. It is not a generally optimal
SUNDIALS configuration.

## Define the experiment

Before writing a tuning file, record:

- the selected package, stepper, vector implementation, matrix/linear solver,
  nonlinear solver, and preconditioner;
- a representative problem size, input, time interval, output schedule, and
  number of runs; include the hard cases that matter in production;
- the objective, such as `wall_time`, solver work, memory, or a reported
  application metric, and whether it should be minimized or maximized;
- a correctness metric and an acceptable upper bound. If changing a setting
  can change the numerical answer, use this as a constraint instead of letting
  speed alone select the result;
- the compiler, optimization flags, SUNDIALS configuration, accelerator and
  CPU placement, thread counts, MPI layout, and relevant environment
  variables; and
- a search budget large enough to compare the proposed space, but small
  enough that repeated confirmation runs remain affordable.

Make the executable return a nonzero status for solver failure. Print stable,
machine-readable metrics, for example:

```text
error_linf=2.1e-8
solver_work=1842
```

`wall_time` is built in and needs no application output or regular expression.
Other objectives and constraints are extracted from `stdout`, `stderr`, or a
file with a regular expression containing the numeric capture group. A list
of regular expressions can form a composite metric using `sum` or `mean`.
The configuration supports one upper-bound constraint. If correctness needs
several checks, have the application return nonzero when any check fails or
emit one aggregate violation metric with a meaningful upper bound.

## Make a C or C++ application tuning-ready

The C/C++ `SetOptions` interface consumes an argument vector containing
key/value pairs. The option key includes its object prefix, such as
`arkode.lsetup_frequency` or `cvode.max_order`; it is not a flag beginning with
`-`. Call the relevant routine after constructing and configuring the solver
stack so command-line values override the defaults:

```c
/* Construct the selected SUNDIALS stack and set application defaults. */
void* arkode_mem = ARKStepCreate(fe, fi, t0, y, ctx);
/* ... tolerances, matrix, linear solver, nonlinear solver, callbacks ... */

int flag = ARKodeSetOptions(arkode_mem, NULL, NULL, argc, argv);
if (flag != ARK_SUCCESS) {
  return 1;
}

/* Evolve the representative workload. Return nonzero on any failure. */
/* ... */
printf("error_linf=%.17g\n", error_linf);
printf("solver_work=%ld\n", solver_work);
```

Use the corresponding routine for the selected package (`CVodeSetOptions`,
`IDASetOptions`, `KINSetOptions`, or `ARKodeSetOptions`). The default package
prefixes are `cvode`, `cvodes`, `ida`, `idas`, `kinsol`, and `arkode`, as
applicable. Options are applied in argument order; if a key is repeated, its
last value wins. For more than one ARKODE integrator, call
`ARKodeSetOptions` for each object with a distinct identifier such as `fast`
and `slow`, then tune `fast.order` and `slow.order` separately. A linear or
nonlinear solver object can likewise be configured with its
`SUNLinSolSetOptions` or `SUNNonlinSolSetOptions` entry point when that object
implements one.

The Fortran interface does not provide the package `*SetOptions` routines.
For a Fortran application, parse the arguments in the application and call
the ordinary Fortran setter interfaces, or add a small C/C interoperability
bridge. `suntools tune` can run the resulting executable, but it cannot make
Fortran consume the appended key/value tokens automatically.

Confirm each key against the guide for the exact package and object. A
successful process is not evidence that a key was used if the application
does not check the `SetOptions` return code or if the prefix belongs to a
different object.

## Separate structural choices from runtime tuning

`suntools tune` appends each sampled parameter to the executable command. It
does not reconstruct the application or replace objects created in C/C++.
Keep these choices in the executable or an application-level selector:

- SUNDIALS package and stepper (`CVODE` versus `ARKStep`, `MRIStep`, and so
  on);
- `N_Vector` implementation and CPU/GPU/MPI backend;
- matrix type, linear solver family, nonlinear solver family, and
  preconditioner; and
- callbacks, explicit/implicit splitting, and data layout.

Compare those alternatives with separate executables built from the same
workload, or with a deliberate application-level mode argument. Do not put
incompatible stack choices in one `parameters` list and assume a
`SetOptions` key will change construction. Some runtime options are
method- or object-specific: for example, ARKStep `table_names` can be a
choice only when the executable already has the compatible ARKStep setup.

Once the stack is fixed, useful tuning candidates usually include method
order/table, meaningful step-size controls, nonlinear convergence controls,
linear-solver iteration limits, setup/Jacobian frequencies, adaptivity
controls, and fused-kernel switches. Administrative failure limits such as
`max_num_steps` should be large enough for valid runs, not optimized for
speed. Change one family at a time when the interaction is not understood.

## Search-space and runner rules

The YAML model accepts three parameter kinds:

- `float` with inclusive `bounds`, optionally `scale: log` for a positive
  multiplicative range;
- `int` with inclusive `bounds`; and
- `choice` with an explicit list of allowed strings.

A choice value is split on whitespace before execution. This supports a
single key that takes multiple arguments, such as an ARKStep table pair or a
two-value `scalar_tolerances` option:

```yaml
parameters:
  - name: arkode.table_names
    type: choice
    values:
      - "ARKODE_ARK324L2SA_DIRK_4_2_3 ARKODE_ERK_NONE"
      - "ARKODE_ARK436L2SA_DIRK_6_3_4 ARKODE_ERK_NONE"
```

Encode every multi-argument value as one categorical string containing the
whitespace-separated tuple. Do not use a continuous parameter for it:
numeric `float` and `int` parameters each emit exactly one command-line
token, while a `choice` emits one token for each whitespace-separated word.

Use only values supported by the selected object. Keep spaces narrow and
meaningful; a large Cartesian product spends evaluations on combinations
that are invalid, unstable, or irrelevant. A practical sequence is:

1. compare a few structurally compatible methods or tables;
2. hold the winner fixed and tune one or two runtime families; then
3. narrow the ranges around promising values and repeat with more
   repetitions.

Every run first executes the executable with no tune parameters as the
baseline; this baseline does not consume `max_evals`. Each sampled
configuration is then executed with the sampled key/value pairs. With
`search.repetitions > 1`, objective, constraint, and wall-time values are
averaged, and every repetition must succeed and satisfy the constraint for
the trial to be feasible.

Use `workers: 1` and `repetitions: 5` or greater for `wall_time`. Concurrent
trials contend for CPU, memory, I/O, and accelerator resources, so their
timings are not directly comparable. Repetitions reduce noise, but they do
not replace validation on independent workloads.

Each trial runs in an isolated temporary current working directory.
`executable.cwd` is the base used to resolve a relative executable command
(and is itself resolved relative to the configuration file); it is not the
subprocess current directory during a trial. Resolve application input files
with absolute paths or make the application locate them robustly from its
executable/configuration path. Relative output files are trial-local and are
not a durable result store. A relative file used as a metric source is read
from that trial directory.

## Generic YAML starting point

This is a valid configuration for an ARKODE executable that accepts the
shown keys and prints `error_linf`. Replace the executable, workload
arguments, and parameter ranges with values justified by the selected stack.

```yaml
backend:
  name: deephyper

search:
  max_evals: 24
  workers: 1
  repetitions: 5
  output_dir: tune-results

executable:
  command: ./arkode_app
  args: [--input, /absolute/path/to/representative-input]
  cwd: .
  env: {}

parameters:
  - name: arkode.lsetup_frequency
    type: int
    bounds: [1, 20]
  - name: arkode.nonlin_conv_coef
    type: float
    bounds: [0.001, 0.3]
    scale: log

objective:
  metric: wall_time
  direction: minimize

constraint:
  metric: error_linf
  source: stdout
  regex: "error_linf=([0-9.eE+-]+)"
  group: 1
  upper_bound: 1.0e-6
```

Run it with:

```bash
suntools tune --config tune.yaml
```

The configuration is validated before a backend starts. The backend must be
installed for the selected name (`deephyper`, `gptune`, or `ytopt`).

## Package-specific candidates

Use the package guide to verify availability and interactions before putting
an option in a space. These are candidate families, not defaults or a reason
to tune every setting:

- **CVODE/CVODES:** `cvode.max_order` or `cvodes.max_order`, initial/max/min
  step, nonlinear convergence controls, linear setup and Jacobian
  frequencies, `eps_lin`, and integrator fused-kernel controls where
  supported.
- **IDA/IDAS:** `ida.max_order` or `idas.max_order`, initial/max/min step,
  consistent-initial-condition limits, nonlinear convergence controls,
  and `eps_lin`. Keep the differential/algebraic partition and consistent-IC
  strategy fixed while measuring solver settings.
- **ARKODE:** `arkode.order` or compatible `arkode.table_names`,
  initial/max/min step, nonlinear convergence and
  adaptivity controls, Jacobian/setup frequencies, and stepper-specific
  options. For `MRIStep` or nested integrators use the distinct prefixes
  supplied when each object is configured.
- **KINSOL:** `kinsol.num_max_iters`, damping and nonlinear stopping
  tolerances, `kinsol.max_newton_step`, and Anderson-acceleration controls
  when a fixed-point or Picard strategy is already appropriate.
- **SUNLinearSolver/SUNNonlinearSolver:** solver-object iteration limits,
  restart or setup controls, and solver-specific tolerances only when the
  selected object exposes `SetOptions`; the identifier must match the prefix
  passed to that object's options routine.

Do not tune an accuracy tolerance as an unconstrained speed knob. If the
requested numerical accuracy is not represented by a constraint or an
independent comparison, a looser and faster but wrong trial can win.

## Read, validate, and materialize a result

The result directory contains machine-readable summaries including:
`baseline.json`, `best.json`, `worst.json`, `results.csv`, and
`trials.jsonl` (plus backend-specific logs). `best.json` records only a
successful feasible trial, including sampled parameters and a reproducible
command. If no successful feasible trial exists, fix the executable, metric
extraction, constraint, or search space before interpreting the search.

For a candidate winner:

1. rerun the command from `best.json` with the same build and environment;
2. compare it with the baseline over repeated runs and the representative
   workload;
3. run independent problem sizes, parameter sets, output schedules, and
   tighter reference tolerances; and
4. check solver return codes, solution norms, conservation/positivity or
   other application invariants, and relevant solver statistics.

Only then copy the selected values into application defaults or a checked-in
configuration. Record the machine, compiler/build, thread/MPI layout,
workload, search space, repetitions, and result files beside the materialized
settings. Retune after changing any of those materially, after changing the
SUNDIALS version, or when the production workload changes.
