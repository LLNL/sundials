%module fsunmatrix_dense_mod

// include code common to all nvector implementations
%include "fsunmatrix.i"

// sunmatrix_impl macro defines some ignore and inserts with the matrix name appended
%sunmatrix_impl(Dense)

// Process and wrap functions in the following files
%include "sunmatrix/sunmatrix_dense.h"

