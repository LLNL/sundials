%module fnvector_openmp_mod

// include code common to all nvector implementations
%include "fnvector.i"

// nvector_impl macro defines some ignore and inserts with the vector name appended
%nvector_impl(OpenMP)

// Process and wrap functions in the following files
%include "nvector/nvector_openmp.h"

