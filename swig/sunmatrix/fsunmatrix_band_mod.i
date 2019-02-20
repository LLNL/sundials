%module fsunmatrix_band_mod

// include code common to all nvector implementations
%include "fsunmatrix.i"

// sunmatrix_impl macro defines some ignore and inserts with the matrix name appended
%sunmatrix_impl(Band)

// Process and wrap functions in the following files
%include "sunmatrix/sunmatrix_band.h"

