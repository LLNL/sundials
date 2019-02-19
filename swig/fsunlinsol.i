// Include shared configuration
%include "fsundials.i"

// Load the typedefs and generate "use" statements
%import "fsundials_types.i"
%import "fnvector_mod.i"
%import "fsunlinsol_mod.i"
%import "fsunmatrix_mod.i"

// Macro for creating an interface to an N_Vector
%define %sunlinsol_impl(TYPE)
  %ignore _SUNLinearSolverContent_## TYPE ##;
%enddef

