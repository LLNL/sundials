%module fsundials_types

%include "fsundials.i"
%include <stdint.i>

// Inform SWIG of the configure-provided types
#define SUNDIALS_EXPORT
#define SUNDIALS_INDEX_TYPE int64_t
#define SUNDIALS_DOUBLE_PRECISION
#define booleantype bool

// Insert code into the C wrapper to check that the sizes match
%{
#include "sundials/sundials.h"

#ifndef SUNDIALS_DOUBLE_PRECISION
#error "The Fortran bindings are only targeted at double-precision"
#endif
const char _your_int_size_is_wrong = ""[sizeof(SUNDIALS_INDEX_TYPE) != 8 ? -1 : 0];
%}

// Process and wrap functions in the following files
%include "sundials/sundials_types.h"

