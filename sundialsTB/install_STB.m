function [] = install_STB
%
% INSTALL_STB Compilation of sundialsTB MEX files

% Radu Serban <radu@llnl.gov>
% Copyright (c) 2005, The Regents of the University of California.
% $Revision: 1.6 $Date: 2006/02/13 23:01:18 $

% Turn on debugging?

dbg = true;

% Location of sundialsTB

stb = pwd;

% Location of SUNDIALS

sun = input('SUNDIALS home directory: ','s');

if ~exist(sun, 'dir')
  error('SUNDIALS not found');
end

% Should we enable parallel support?

par = true;
if isempty(getenv('LAMHOME'))
  par = false;
end
if isempty(getenv('MPITB_ROOT'))
  par = false;
end
q = [sun '/nvec_par'];
if ~exist(q, 'dir')
  par = false;
end

% Compile cvm MEX file

q = [sun '/cvodes'];
if exist(q,'dir')
  install_CVM(stb,sun,par,dbg);
end

% Compile kim MEX file

q = [sun '/kinsol'];
if exist(q,'dir')
  install_KIM(stb,sun,par,dbg);
end

%---------------------------------------------------------------------------------
% compilation of cvm MEX file
%---------------------------------------------------------------------------------

function [] = install_CVM(stb,sun,par,dbg)

% Move to cvm MEX directory
q = [stb '/src/cvm'];
cd(q);

includes = '-I.. -I../nvm';
if par
  sources = 'cvm.c cvmWrap.c cvmOpts.c ../nvm/nvm_parallel.c ../nvm/nvm_ops.c';
else
  sources = 'cvm.c cvmWrap.c cvmOpts.c ../nvm/nvm_serial.c ../nvm/nvm_ops.c';
end
libraries = '';

% Add CVODES sources and header files

cvodes_src = {'cvodes_band.c'
              'cvodes_bandpre.c'
              'cvodes_bbdpre.c'
              'cvodes_dense.c'
              'cvodes_diag.c'
              'cvodea.c'
              'cvodes.c'
              'cvodes_io.c'
              'cvodea_io.c'
              'cvodes_spils.c'
              'cvodes_spbcgs.c'
              'cvodes_spgmr.c'
              'cvodes_sptfqmr.c'};

for i=1:length(cvodes_src)
  tmp = strcat(sun,'/cvodes/source/',cvodes_src{i});
  sources = sprintf('%s %s',sources,tmp);
end

includes = sprintf('%s -I%s/cvodes/include -I%s/cvodes/source',includes,sun,sun);

% Add SHARED sources and header files

shared_src = {'sundials_band.c'
              'sundials_dense.c'
              'sundials_iterative.c'
              'sundials_nvector.c'
              'sundials_smalldense.c'
              'sundials_spbcgs.c'
              'sundials_spgmr.c'
              'sundials_sptfqmr.c'
              'sundials_math.c'};

for i=1:length(shared_src)
  tmp = strcat(sun,'/shared/source/',shared_src{i});
  sources = sprintf('%s %s',sources,tmp);
end

includes = sprintf('%s -I%s/shared/include',includes,sun);

% Add NVEC_SER sources and header files

nvecser_src = {'nvector_serial.c'};

for i=1:length(nvecser_src)
  tmp = strcat(sun,'/nvec_ser/',nvecser_src{i});
  sources = sprintf('%s %s',sources,tmp);
end

includes = sprintf('%s -I%s/nvec_ser',includes,sun);


if par
  
% Add NVEC_PAR sources and header files

  nvecpar_src = {'nvector_parallel.c'};

  for i=1:length(nvecpar_src)
    tmp = strcat(sun,'/nvec_par/',nvecpar_src{i});
    sources = sprintf('%s %s',sources,tmp);
  end
  
  includes = sprintf('%s -I%s/nvec_par',includes,sun);

% Add LAM headers and libraries
  
  lam = getenv('LAMHOME');
  includes = sprintf('%s -I%s/include',includes,lam);
  libraries = sprintf('%s -L%s/lib -lmpi -llam -lutil',libraries,lam);
  
end

% Create MEX file

cvm = [stb '/cvodes/cvm'];
if dbg
  mex_cmd = 'mex -g';
else
  mex_cmd = 'mex';
end
mex_cmd = sprintf('%s -v -outdir %s %s %s %s', mex_cmd, cvm, includes, sources, libraries);
eval(mex_cmd);

% Move back to sundialsTB

cd(stb)


%---------------------------------------------------------------------------------
% compilation of KINSOL MEX file
%---------------------------------------------------------------------------------

function [] = install_KIM(stb,sun,par,dbg)

% Move to kim MEX directory
q = [stb '/src/kim'];
cd(q);

includes = '-I.. -I../nvm';
if par
  sources = 'kim.c kimWrap.c kimOpts.c ../nvm/nvm_parallel.c ../nvm/nvm_ops.c';
else
  sources = 'kim.c kimWrap.c kimOpts.c ../nvm/nvm_serial.c ../nvm/nvm_ops.c';
end
libraries = '';

% Add KINSOL sources and header files

kinsol_src = {'kinsol_band.c'
              'kinsol_bbdpre.c'
              'kinsol_dense.c'
              'kinsol.c'
              'kinsol_io.c'
              'kinsol_spils.c'
              'kinsol_spbcgs.c'
              'kinsol_spgmr.c'
              'kinsol_sptfqmr.c'};

for i=1:length(kinsol_src)
  tmp = strcat(sun,'/kinsol/source/',kinsol_src{i});
  sources = sprintf('%s %s',sources,tmp);
end

includes = sprintf('%s -I%s/kinsol/include -I%s/kinsol/source',includes,sun,sun);

% Add SHARED sources and header files

shared_src = {'sundials_band.c'
              'sundials_dense.c'
              'sundials_iterative.c'
              'sundials_nvector.c'
              'sundials_smalldense.c'
              'sundials_spbcgs.c'
              'sundials_spgmr.c'
              'sundials_sptfqmr.c'
              'sundials_math.c'};

for i=1:length(shared_src)
  tmp = strcat(sun,'/shared/source/',shared_src{i});
  sources = sprintf('%s %s',sources,tmp);
end

includes = sprintf('%s -I%s/shared/include',includes,sun);

% Add NVEC_SER sources and header files

nvecser_src = {'nvector_serial.c'};

for i=1:length(nvecser_src)
  tmp = strcat(sun,'/nvec_ser/',nvecser_src{i});
  sources = sprintf('%s %s',sources,tmp);
end

includes = sprintf('%s -I%s/nvec_ser',includes,sun);


if par
  
% Add NVEC_PAR sources and header files

  nvecpar_src = {'nvector_parallel.c'};

  for i=1:length(nvecpar_src)
    tmp = strcat(sun,'/nvec_par/',nvecpar_src{i});
    sources = sprintf('%s %s',sources,tmp);
  end
  
  includes = sprintf('%s -I%s/nvec_par',includes,sun);

% Add LAM headers and libraries
  
  lam = getenv('LAMHOME');
  includes = sprintf('%s -I%s/include',includes,lam);
  libraries = sprintf('%s -L%s/lib -lmpi -llam -lutil',libraries,lam);
    
end

% Create MEX file

kim = [stb '/kinsol/kim'];
if dbg
  mex_cmd = 'mex -g';
else
  mex_cmd = 'mex';
end
mex_cmd = sprintf('%s -v -outdir %s %s %s %s', mex_cmd, kim, includes, sources, libraries);
eval(mex_cmd);

% Move back to sundialsTB

cd(stb)


