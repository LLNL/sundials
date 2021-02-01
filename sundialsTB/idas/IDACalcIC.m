function [status, varargout] = IDACalcIC(tout,icmeth)
%IDACalcIC computes consistent initial conditions
%
%   Usage: STATUS = IDACalcIC ( TOUT, ICMETH )
%          [STATUS, YY0, YP0] = IDACalcIC ( TOUT, ICMETH )
%
%   IDACalcIC corrects the guess for initial conditions passed
%   to IDAInit or IDAReInit so that the algebraic constraints
%   are satisfied. 
%
%   The argument TOUT is the first value of t at which a soluton will be      
%   requested (from IDASolve). This is needed here to determine the 
%   direction of integration and rough scale in the independent variable. 
%
%   If ICMETH is 'FindAlgebraic', then IDACalcIC attempts to compute 
%   the algebraic components of y and differential components of y', 
%   given the differential components of y.  
%   This option requires that the vector id was set through IDASetOptions  
%   specifying the differential and algebraic components.
%   If ICMETH is 'FindAll', then IDACalcIC attempts to compute all  
%   components of y, given y'.  In this case, id is not required. 
%
%   On return, STATUS is one of the following:
% SUCCESS             IDACalcIC was successful.  The corrected   
%                     initial value vectors are in y0 and yp0.
% IDA_MEM_NULL        The argument ida_mem was NULL.             
% IDA_ILL_INPUT       One of the input arguments was illegal.    
%                     See printed message.                       
% IDA_LINIT_FAIL      The linear solver's init routine failed.   
% IDA_BAD_EWT         Some component of the error weight vector  
%                     is zero (illegal), either for the input    
%                     value of y0 or a corrected value.          
% IDA_RES_FAIL        The user's residual routine returned 
%                     a non-recoverable error flag.              
% IDA_FIRST_RES_FAIL  The user's residual routine returned 
%                     a recoverable error flag on the first call,
%                     but IDACalcIC was unable to recover.       
% IDA_LSETUP_FAIL     The linear solver's setup routine had a    
%                     non-recoverable error.                     
% IDA_LSOLVE_FAIL     The linear solver's solve routine had a    
%                     non-recoverable error.                     
% IDA_NO_RECOVERY     The user's residual routine, or the linear 
%                     solver's setup or solve routine had a      
%                     recoverable error, but IDACalcIC was       
%                     unable to recover.                         
% IDA_CONSTR_FAIL     IDACalcIC was unable to find a solution    
%                     satisfying the inequality constraints.     
% IDA_LINESEARCH_FAIL The Linesearch algorithm failed to find a  
%                     solution with a step larger than steptol   
%                     in weighted RMS norm.
% IDA_CONV_FAIL       IDACalcIC failed to get convergence of the 
%                     Newton iterations.  
%
%   If the output arguments YY0 and YP0 are present, they will
%   contain the consistent initial conditions.
%
%  See also: IDASetOptions, IDAInit, IDAReInit

% Radu Serban <radu@llnl.gov>
% SUNDIALS Copyright Start
% Copyright (c) 2002-2021, Lawrence Livermore National Security
% and Southern Methodist University.
% All rights reserved.
%
% See the top-level LICENSE and NOTICE files for details.
%
% SPDX-License-Identifier: BSD-3-Clause
% SUNDIALS Copyright End
% $Revision$Date: 2007/02/05 20:23:46 $

mode = 25;

if nargout == 1
  status = idm(mode, tout, icmeth);
elseif nargout == 3
  [status, yy, yp] = idm(mode, tout, icmeth);
  varargout(1) = {yy};
  varargout(2) = {yp};
else
  disp('IDACalcIC:: wrong number of output arguments');
end
