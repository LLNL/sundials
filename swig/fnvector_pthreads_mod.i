%module fnvector_pthreads_mod

// include code common to all nvector implementations
%include "fnvector.i"

// nvector_impl macro defines some ignore and inserts with the vector name appended
%nvector_impl(Pthreads)

// ignore the Pthreads_Data struct  
%ignore _Pthreads_Data;
  
// Process and wrap functions in the following files
%include "nvector/nvector_pthreads.h"

