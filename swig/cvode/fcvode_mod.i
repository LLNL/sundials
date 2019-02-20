%module fcvode_mod

%include "../sundials/fsundials.i"

// Load the typedefs and generate a "use" statements in the module
%import "../sundials/fnvector_mod.i"

// Process definitions from these files
%include "cvode/cvode.h"
%include "cvode/cvode_bandpre.h"
%include "cvode/cvode_bbdpre.h"
%include "cvode/cvode_diag.h"
%include "cvode/cvode_ls.h"

