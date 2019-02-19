%module fsunlinsol_mod

// Load the typedefs and generate a "use" statement in the module
%import "fsundials_types.i"
%import "fsunmatrix_mod.i"
%import "fnvector_mod.i"

// Process and wrap functions in the following files
%include "sundials/sundials_iterative.h"
%include "sundials/sundials_linearsolver.h"

