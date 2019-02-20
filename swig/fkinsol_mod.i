%module fkinsol_mod

%include "fsundials.i"

// Load the typedefs and generate a "use" statements in the module
%import "fnvector_mod.i"
%import "fsunlinsol_mod.i"
%import "fsunmatrix_mod.i"

// Process definitions from these files
%include "kinsol/kinsol.h"
%include "kinsol/kinsol_bbdpre.h"
%include "kinsol/kinsol_ls.h"

