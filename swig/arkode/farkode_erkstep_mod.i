%module farkode_erkstep_mod

%include "../sundials/fsundials.i"

// Load the typedefs and generate a "use" statements in the module
%import "farkode_mod.i"

// Process definitions from these files
%include "arkode/arkode_erkstep.h"

