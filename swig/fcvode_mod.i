%module fcvode_mod

%include "fsundials.i"

// Load the typedefs and generate a "use" statements in the module
%import "fnvector.i"
%import "fsunlinsol.i"
%import "fsunnonlinsol.i"
%import "fsunmatrix.i"

// Process definitions from these files
%include "cvode/cvode.h"
%include "cvode/cvode_ls.h"
