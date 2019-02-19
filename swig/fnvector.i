// Include shared configuration
%include "fsundials.i"

// Load the typedefs and generate "use" statements
%import "fsundials_types.i"
%import "fnvector_mod.i"

// Assume double* is an array of doubles
// NOTE: F2003 does not allow assumed-shape, but does allow assumed-size
%typemap(bindc, in="real(C_DOUBLE), dimension(*)") realtype* v_data "type(C_PTR)"

// Macro for creating an interface to an N_Vector
%define %nvector_impl(TYPE)
  %ignore _N_VectorContent_## TYPE ##;

  %insert(fpublic) %{
  public :: FN_VGetData_## TYPE ##
  %}

  %insert(fwrapper) %{
  subroutine FN_VGetData_## TYPE ##(vec, vdata)

      use, intrinsic :: iso_c_binding
      implicit none

      type(C_PTR)        :: vec
      integer(C_INT64_T) :: len
      type(C_PTR)        :: cptr
      real(C_DOUBLE), dimension(:), pointer :: vdata

      len = FN_VGetLength_## TYPE ##(vec)
      cptr = FN_VGetArrayPointer_## TYPE ##(vec)

      call c_f_pointer(cptr, vdata, (/len/))

  end subroutine FN_VGetData_## TYPE ##
  %}
%enddef

