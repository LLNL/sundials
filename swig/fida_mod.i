%module fida_mod

%include "fsundials.i"

// Load the typedefs and generate a "use" statements in the module
%import "fnvector_mod.i"
%import "fsunlinsol_mod.i"
%import "fsunnonlinsol_mod.i"
%import "fsunmatrix_mod.i"

// Process definitions from these files
%include "ida/ida.h"
%include "ida/ida_bbdpre.h"
%include "ida/ida_ls.h"

