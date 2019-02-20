%module farkode_erkstep_mod

%include "fsundials.i"

// Load the typedefs and generate a "use" statements in the module
%import "fsundials_types.i"
%import "farkode_mod.i"
%import "fnvector_mod.i"
%import "fsunlinsol_mod.i"
%import "fsunnonlinsol_mod.i"
%import "fsunmatrix_mod.i"

// Process definitions from these files
%include "arkode/arkode_erkstep.h"

