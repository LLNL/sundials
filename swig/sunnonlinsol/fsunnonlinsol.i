// Include shared configuration
%include "../sundials/fsundials.i"

// Load the typedefs and generate "use" statements
%import "../sundials/fsunnonlinsol_mod.i"

// Macro for creating an interface to an N_Vector
%define %sunnonlinsol_impl(TYPE)
  %ignore _SUNNonlinearSolverContent_## TYPE ##;
%enddef

