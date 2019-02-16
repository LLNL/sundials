// Provide direct BIND(C) interfaces to *all* functions instead of generating wrapper functions
%fortranbindc;
// By default, wrap all constants as native fortran PARAMETERs
%fortranconst;

// Treat FILE* as an opaque pointer
%apply void * { FILE * };

// Prefix all functions with F
// E.g. CVodeCreate -> FCVodeCreate
%rename("F%s", %$isfunction) "";

// Macro for treating the opaque types that SUNDIALS defines as void pointers
%define %sundials_generic(TYPE)
  %ignore _generic_ ## TYPE ## _Ops;
  %ignore _generic_ ## TYPE;
  %apply void * { _generic_ ## TYPE *, TYPE };
%enddef

%sundials_generic(N_Vector)
%sundials_generic(SUNLinearSolver)
%sundials_generic(SUNNonlinearSolver)
%sundials_generic(SUNMatrix)

// Macro for ignoring the content structures
%define %nvector_impl(TYPE)
  %ignore _N_VectorContent_## TYPE ##;
%enddef

