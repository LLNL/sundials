% SUNDIALS NVECTOR operations
%
% Functions:
%
%   N_VMax      - returns the largest element of x
%   N_VMaxNorm  - returns the maximum norm of x
%   N_VMin      - returns the smallest element of x
%   N_VDotProd  - returns the dot product of two vectors
%   N_VWrmsNorm - returns the weighted root mean square norm of x
%   N_VWL2Norm  - returns the weighted Euclidean L2 norm of x
%   N_VL1Norm   - returns the L1 norm of x
%
% NOTE For serial vectors, all of the above operations default to
%   the corresponding MATLAB functions. For parallel vectors, they
%   can be used either on the local portion of the distributed vector
%   or on the global vector (in which case they will trigger an MPI
%   allreduce operation).
