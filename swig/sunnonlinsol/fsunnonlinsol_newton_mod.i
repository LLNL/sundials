%module fsunnonlinsol_newton_mod

// include code common to all implementations
%include "fsunnonlinsol.i"

%sunnonlinsol_impl(Newton)

// Process and wrap functions in the following files
%include "sunnonlinsol/sunnonlinsol_newton.h"

