// Provide direct BIND(C) interfaces to *all* functions instead of generating wrapper functions
%fortranbindc;
// By default, wrap all constants as native fortran PARAMETERs
%fortranconst;

// Prefix all functions with F
// E.g. CVodeCreate -> FCVodeCreate
%rename("F%s", %$isfunction) "";

// Macro for treating the opaque types that SUNDIALS defines as void pointers
%define %sundials_generic(TYPE)
  %ignore _generic_ ## TYPE ## _Ops;
  %ignore _generic_ ## TYPE;
  %apply void * { _generic_ ## TYPE *, TYPE };
%enddef

// Treat sundials generics as void pointers
%sundials_generic(N_Vector)
%sundials_generic(SUNLinearSolver)
%sundials_generic(SUNNonlinearSolver)
%sundials_generic(SUNMatrix)

// Treat FILE* as an opaque pointer
%apply void * { FILE * };

// Don't use value attribute on void**
%typemap(bindc,in="type(C_PTR)") void** "type(C_PTR)"

