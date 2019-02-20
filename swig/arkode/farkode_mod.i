%module farkode_mod

%include "../sundials/fsundials.i"

// Load the typedefs and generate a "use" statements in the module
%import "../sundials/fsundials_types.i"

%apply void * { ARKodeButcherTable };

// Process definitions from these files
%include "arkode/arkode.h"
%include "arkode/arkode_bandpre.h"
%include "arkode/arkode_bbdpre.h"
%include "arkode/arkode_butcher.h"
%include "arkode/arkode_butcher_dirk.h"
%include "arkode/arkode_butcher_erk.h"
%include "arkode/arkode_ls.h"

