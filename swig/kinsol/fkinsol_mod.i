%module fkinsol_mod

%include "../sundials/fsundials.i"

// Load the typedefs and generate a "use" statements in the module
%import "../sundials/fsundials_types.i"

// Process definitions from these files
%include "kinsol/kinsol.h"
%include "kinsol/kinsol_bbdpre.h"
%include "kinsol/kinsol_ls.h"

