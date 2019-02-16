%module fnvector_serial_mod

%include "fsundials.i"

// Load the typedefs and generate "use" statements
%import "fsundials_types.i"
%import "fnvector.i"

// Ignore content structure typedef etc.
%nvector_impl(Serial)

// Assume double* is an array of doubles
%typemap(bindc, in="type(C_PTR)") double* "type(C_PTR)"

// Add function to make it easier to get data out
%insert(fpublic) %{
public :: FN_VGetData_Serial
%}

%insert(fwrapper) %{
subroutine FN_VGetData_Serial(vec, vdata)

    use, intrinsic :: iso_c_binding
    implicit none

    type(C_PTR)        :: vec
    integer(C_INT64_T) :: len
    type(C_PTR)        :: cptr
    real(C_DOUBLE), dimension(:), pointer :: vdata

    len = FN_VGetLength_Serial(vec)
    cptr = FN_VGetArrayPointer_Serial(vec)

    call c_f_pointer(cptr, vdata, (/len/)) 

end subroutine FN_VGetData_Serial
%}

// Process and wrap functions in the following files
%include "nvector/nvector_serial.h"

