# SUNDIALS Testing

SUNDIALS utilizes CTest to build and run the example programs distributed with
SUNDIALS as a test suite. In order to run the tests SUNDIALS examples must be
enabled by setting the CMake options `EXAMPLES_ENABLE_<language>` to `ON`.

To run the standard C tests use the following commands:
```
mkdir <build directory>
cd <build directory>
cmake <sundials directory> -DCMAKE_INSTALL_PREFIX=<prefix> -DEXAMPLES_ENABLE_C=ON
make
make test
```
This will create an out of source build, build sundials and all the C examples,
and run the tests (which are the examples) that return a pass/fail status. All
output is saved to the file `<build directory>/Testing/Temporary/LastTest.log`.

### Development Test

Examples that do not return a pass/fail status are considered "development" only
tests and are excluded by default when running `make test`. To enable all
examples when testing set the CMake option `SUNDIALS_DEVTESTS` to `ON`. With
this option enabled tests are run using a Python test runner `test/testRunner`.
The test runner determines if a test passes or fails by comparing the test
output against saved output files that are deemed correct by the SUNDIALS team.

To run the development C tests use the following commands:
```
mkdir <build directory>
cd <build directory>
cmake <sundials directory> -DCMAKE_INSTALL_PREFIX=<prefix> -DEXAMPLES_ENABLE_C=ON -DSUNDIALS_DEVTESTS=ON
make
make test
```
This will create an out of source build, build sundials and all the C examples,
and run all the tests (which are the examples). The output from each test is
saved to `<build directory>/Testing/output`.

## Testing Scripts

Several shell scripts are provided for setting up and running SUNDIALS tests and
are used with continuous integration tools to test each push and pull-request to
the SUNDIALS repository. The primary driver script is `test_driver.sh` which
sets up the testing environment and selects which test type to run: BRANCH, PR,
or RELEASE. The typical usage of `test_driver.sh` is as follows:
```
./test_driver.sh <number of build threads> <test type>
```

The main test driver configures the testing environment by sourcing one of the
following scripts (listed in the order checked):

1. Script specificed by the environment variable `SUNDIALS_ENV`
2. A user's local environment script: `sunrepo/test/env.sh`
3. A user's global environment script: `~/.sundials_config/env.sh`
4. The default SUNDIALS environment script: `sunrepo/test/env.default.sh`

Environment scripts should set the following environment variables that are
used when configuring SUNDIALS for testing:
```
CC  = C compiler
FC  = Fortran compiler
CXX = C++ compiler

MPIDIR  = path to MPI installation
MPIEXEC = executble for launching MPI runs

BLAS_LIBRARIES   = full path to BLAS library
LAPACK_LIBRARIES = full path to LAPACK library

KLUDIR = path to KLU installation

SLUMTDIR_32 = path to 32-bit int SuperLU_MT installation
SLUMTDIR_64 = path to 64-bit int SuperLU_MT installation

SLUDISTDIR_32 = path to 32-bit int SuperLU_DIST installation
SLUDISTDIR_64 = path to 64-bit int SuperLU_DIST installation

HYPREDIR_32 = path to 32-bit int hypre installation
HYPREDIR_64 = path to 64-bit int hypre installation

PETSCDIR_32 = path to 32-bit int PETSc installation
PETSCDIR_64 = path to 64-bit int PETSc installation
```
Note: At this time the testing scripts only run development tests when SUNDIALS
is configured with real type double (either index size can be used).

Once the environment is setup the main driver calls one of three sub-driver
scripts depending on the test type:
* BRANCH tests (`test_branch.sh`) are run on each push to a branch in the
SUNDIALS repository.
* PR tests (`test_pr.sh`) are run for each pull-request issued to the SUNDIALS
repository.
* RELEASE tests (`test_release.sh`) are run for each push to a release branch or
release pull-request.

The BRANCH, PR, and RELEASE scripts call one or more of the following test
scripts to configure, build, and test SUNDIALS:
* `suntest.sh` -- tests configured using standard options
* `suntest_xsdk.sh` -- tests configured using xSDK options
* `suntest_tarscript.sh` -- create tarballs and test with standard options

Note: Any of the scripts may be called independently of the main driver script
for manual testing. See the comment block at the top of each file for more
information on running scripts.