// Include shared configuration
%include "../sundials/fsundials.i"

// Load the typedefs and generate "use" statements
%import "../sundials/fsunlinsol_mod.i"

// Macro for creating an interface to an N_Vector
%define %sunlinsol_impl(TYPE)
  %ignore _SUNLinearSolverContent_## TYPE ##;
%enddef

