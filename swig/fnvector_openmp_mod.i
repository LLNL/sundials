%module fnvector_openmp_mod

%include "fsundials.i"

// Load the typedefs and generate "use" statements
%import "fsundials_types.i"
%import "fnvector.i"

// Ignore content structure typedef etc.
%nvector_impl(OpenMP)

// Assume double* is an array of doubles
%typemap(bindc, in="real(C_DOUBLE), dimension(*)") double* "real(C_DOUBLE)"

%insert(fwrapper) %{
subroutine FN_VGetData_OpenMP(vec, varray)

    use, intrinsic :: iso_c_binding
    implicit none

    type(c_ptr)     :: vec
    integer(c_long) :: length
    real(c_double)  :: vptr
    real(c_double), dimension(:) :: varray

    varray = FN_VGetArrayPointer_OpenMP(vec)

end subroutine FN_VGetData_OpenMP
%}

// Process and wrap functions in the following files
%include "nvector/nvector_openmp.h"

