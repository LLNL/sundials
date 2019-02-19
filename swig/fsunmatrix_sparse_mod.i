%module fsunmatrix_sparse_mod

// include code common to all nvector implementations
%include "fsunmatrix.i"

// sunmatrix_impl macro defines some ignore and inserts with the matrix name appended
%sunmatrix_impl(Sparse)

// Process and wrap functions in the following files
%include "sunmatrix/sunmatrix_sparse.h"

