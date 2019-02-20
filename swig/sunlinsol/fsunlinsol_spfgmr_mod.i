%module fsunlinsol_spfgmr_mod

// include code common to all nvector implementations
%include "fsunlinsol.i"

// sunlinsol_impl macro defines some ignore and inserts with the linear solver name appended
%sunlinsol_impl(SPFGMR)

// Process and wrap functions in the following files
%include "sunlinsol/sunlinsol_spfgmr.h"

