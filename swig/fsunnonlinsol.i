// Include shared configuration
%include "fsundials.i"

// Load the typedefs and generate "use" statements
%import "fsundials_types.i"
%import "fsunnonlinsol_mod.i"
%import "fnvector_mod.i"

// Macro for creating an interface to an N_Vector
%define %sunnonlinsol_impl(TYPE)
  %ignore _SUNNonlinearSolverContent_## TYPE ##;
%enddef

