%module fsunmatrix_mod

// Load the typedefs and generate a "use" statement in the module
%import "../sundials/fsundials_types.i"

// Process and wrap functions in the following files
%include "sundials/sundials_matrix.h"

