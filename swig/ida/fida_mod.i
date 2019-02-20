%module fida_mod

%include "../sundials/fsundials.i"

// Load the typedefs and generate a "use" statements in the module
%import "../sundials/fsundials_types.i"

// Process definitions from these files
%include "ida/ida.h"
%include "ida/ida_bbdpre.h"
%include "ida/ida_ls.h"

