%module farkode_mod

%include "fsundials.i"

%apply void * { ARKodeButcherTable };

// Load the typedefs and generate a "use" statements in the module
%import "fsundials_types.i"
%import "fnvector_mod.i"
%import "fsunlinsol_mod.i"
%import "fsunnonlinsol_mod.i"
%import "fsunmatrix_mod.i"

// Process definitions from these files
%include "arkode/arkode.h"
%include "arkode/arkode_bandpre.h"
%include "arkode/arkode_bbdpre.h"
%include "arkode/arkode_butcher.h"
%include "arkode/arkode_butcher_dirk.h"
%include "arkode/arkode_butcher_erk.h"
%include "arkode/arkode_ls.h"

