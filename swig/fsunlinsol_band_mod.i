%module fsunlinsol_band_mod

// include code common to all nvector implementations
%include "fsunlinsol.i"

// sunlinsol_impl macro defines some ignore and inserts with the linear solver name appended
%sunlinsol_impl(Band)

// Process and wrap functions in the following files
%include "sunlinsol/sunlinsol_band.h"

