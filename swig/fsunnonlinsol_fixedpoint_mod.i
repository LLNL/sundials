%module fsunnonlinsol_fixedpoint_mod

// include code common to all implementations
%include "fsunnonlinsol.i"

%sunnonlinsol_impl(FixedPoint)

// Process and wrap functions in the following files
%include "sunnonlinsol/sunnonlinsol_fixedpoint.h"

