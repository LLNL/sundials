function [] = install_STB
%
% INSTALL_STB Interactive compilation and installtion of sundialsTB

% Radu Serban <radu@llnl.gov>
% Copyright (c) 2005, The Regents of the University of California.
% $Revision: 1.17 $Date: 2007/05/11 18:52:47 $

% MEX compiler command
% --------------------

mexcompiler = 'mex';

% Location of sundialsTB and top of sundials source tree
% ------------------------------------------------------

stb = pwd;
cd('..');
sun = pwd;
cd(stb);

% Should we enable parallel support?
% ----------------------------------

par = true;
if isempty(getenv('LAMHOME'))
    par = false;
end
if isempty(getenv('MPITB_ROOT'))
    par = false;
end
q = fullfile(sun,'src','nvec_par');
if ~exist(q, 'dir')
    par = false;
end

% Create sundials_config.h
% ------------------------

mkdir('sundials');
fi = fopen(fullfile('sundials','sundials_config.h'),'w');
fprintf(fi,'#define SUNDIALS_PACKAGE_VERSION "2.3.0"\n');
fprintf(fi,'#define SUNDIALS_DOUBLE_PRECISION 1\n');
fprintf(fi,'#define SUNDIALS_USE_GENERIC_MATH 1\n');
fprintf(fi,'#define SUNDIALS_EXPORT\n');
fclose(fi);

% Compile MEX file
% ----------------

cvm_ok = false;
q = fullfile(sun,'src','cvodes');
if exist(q, 'dir')
  answ = input('    Compile CVODES interface? (y/n) ','s');
  if answ == 'y'
    compile_CVM(mexcompiler,stb,sun,par);
    cvm_ok = true;
  end
end
  
idm_ok = false;
q = fullfile(sun,'src','ida');
if exist(q, 'dir')
  answ = input('    Compile IDAS interface? (y/n) ','s');
  if answ == 'y'
    compile_IDM(mexcompiler,stb,sun,par);
    idm_ok = false;
  end
end
  
kim_ok = false;
q = fullfile(sun,'src','kinsol');
if exist(q, 'dir')
  answ = input('    Compile KINSOL interface? (y/n) ','s');
  if answ == 'y'
    compile_KIM(mexcompiler,stb,sun,par);
    kim_ok = true;
  end
end

% Remove sundials_config.h
% ------------------------

rmdir('sundials','s');

% Install sundialsTB
% ------------------

fprintf('\nMEX files were successfully created.\n');

answ = input('    Install toolbox? (y/n) ','s');
if answ ~= 'y'
  fprintf('\nOK. All done.\n');
  return
end

while true
  fprintf('\nSpecify the location where you wish to install the toolbox.\n');
  fprintf('The toolbox will be installed in a subdirectory "sundialsTB".\n');
  fprintf('Enter return to cancel the installation.\n');
  where = input('    Installation directory: ','s');
  if isempty(where)
    go = 0;
    break;
  end
  if exist(where,'dir')
    go = 1;
    break
  end
  fprintf('\n%s does not exist!\n', where);
end

if ~go
  fprintf('\nOK. All done.\n');
  return
end

stbi = fullfile(where,'sundialsTB');

go = 1;
if exist(stbi,'dir')
  fprintf('\nDirectory %s exists!\n',stbi);
  answ = input('    Replace? (y/n) ','s');
  if answ == 'y'
    rmdir(stbi,'s');
    go = 1;
  else
    go = 0;
  end
end

if ~go
  fprintf('\nOK. All done.\n');
  return
end

mkdir(where,'sundialsTB');
mkdir(fullfile(where,'sundialsTB'),'nvector');
if par
  mkdir(fullfile(where,'sundialsTB'),'putils');
end

instSTB(stb, where, par);

if cvm_ok
  mkdir(fullfile(where,'sundialsTB'),'cvodes');
  mkdir(fullfile(where,'sundialsTB','cvodes'),'cvm');
  mkdir(fullfile(where,'sundialsTB','cvodes'),'function_types');
  mkdir(fullfile(where,'sundialsTB','cvodes'),'examples_ser');
  if par
    mkdir(fullfile(where,'sundialsTB','cvodes'),'examples_par');
  end
  instCVM(stb, where, par);
end

if idm_ok
  mkdir(fullfile(where,'sundialsTB'),'idas');
  mkdir(fullfile(where,'sundialsTB','idas'),'idm');
  mkdir(fullfile(where,'sundialsTB','idas'),'function_types');
  mkdir(fullfile(where,'sundialsTB','idas'),'examples_ser');
  if par
    mkdir(fullfile(where,'sundialsTB','idas'),'examples_par');
  end
  instIDM(stb, where, par);
end

if kim_ok
  mkdir(fullfile(where,'sundialsTB'),'kinsol');
  mkdir(fullfile(where,'sundialsTB','kinsol'),'kim');
  mkdir(fullfile(where,'sundialsTB','kinsol'),'function_types');
  mkdir(fullfile(where,'sundialsTB','kinsol'),'examples_ser');
  if par
    mkdir(fullfile(where,'sundialsTB','kinsol'),'examples_par');
  end
  instKIM(stb, where, par);
end
  
fprintf('\nThe sundialsTB toolbox was installed in %s\n',stbi);
fprintf('\nA startup file, "startup_STB.m" was created in %s.\n',stbi);
fprintf('Use it as your Matlab startup file, or, if you already have a startup.m file,\n');
fprintf('add a call to %s\n',fullfile(stbi,'startup_STB.m'));
fprintf('\nEnjoy!\n\n');

%---------------------------------------------------------------------------------
% compilation of cvm MEX file
%---------------------------------------------------------------------------------

function [] = compile_CVM(mexcompiler,stb,sun,par)

cvm_sources = {
    fullfile(stb,'cvodes','cvm','src','cvm.c')
    fullfile(stb,'cvodes','cvm','src','cvmWrap.c')
    fullfile(stb,'cvodes','cvm','src','cvmOpts.c')
              };
if par
    nvm_sources = {
        fullfile(stb,'nvector','src','nvm_parallel.c')
        fullfile(stb,'nvector','src','nvm_ops.c')
                  };
else
    nvm_sources = {
        fullfile(stb,'nvector','src','nvm_serial.c')
        fullfile(stb,'nvector','src','nvm_ops.c')
                  };
end
sources = '';
for i=1:length(cvm_sources)
    sources = sprintf('%s "%s"',sources,cvm_sources{i});
end
for i=1:length(nvm_sources)
    sources = sprintf('%s "%s"',sources,nvm_sources{i});
end

cvm_incdir = fullfile(stb,'cvodes','cvm','src'); % for cvm.h
nvm_incdir = fullfile(stb,'nvector','src');      % for nvm.h
includes = sprintf('-I"%s" -I"%s" -I"%s"',stb,cvm_incdir,nvm_incdir);

libraries = '';

% Add CVODES sources and header files

cvs_sources = {
    fullfile(sun,'src','cvodes','cvodes_band.c')
    fullfile(sun,'src','cvodes','cvodes_bandpre.c')
    fullfile(sun,'src','cvodes','cvodes_bbdpre.c')
    fullfile(sun,'src','cvodes','cvodes_direct.c')
    fullfile(sun,'src','cvodes','cvodes_dense.c')
    fullfile(sun,'src','cvodes','cvodes_diag.c')
    fullfile(sun,'src','cvodes','cvodea.c')
    fullfile(sun,'src','cvodes','cvodes.c')
    fullfile(sun,'src','cvodes','cvodes_io.c')
    fullfile(sun,'src','cvodes','cvodea_io.c')
    fullfile(sun,'src','cvodes','cvodes_spils.c')
    fullfile(sun,'src','cvodes','cvodes_spbcgs.c')
    fullfile(sun,'src','cvodes','cvodes_spgmr.c')
    fullfile(sun,'src','cvodes','cvodes_sptfqmr.c')
              };
shr_sources = {
    fullfile(sun,'src','sundials','sundials_band.c')
    fullfile(sun,'src','sundials','sundials_dense.c')
    fullfile(sun,'src','sundials','sundials_iterative.c')
    fullfile(sun,'src','sundials','sundials_nvector.c')
    fullfile(sun,'src','sundials','sundials_direct.c')
    fullfile(sun,'src','sundials','sundials_spbcgs.c')
    fullfile(sun,'src','sundials','sundials_spgmr.c')
    fullfile(sun,'src','sundials','sundials_sptfqmr.c')
    fullfile(sun,'src','sundials','sundials_math.c')
              };
for i=1:length(cvs_sources)
    sources = sprintf('%s "%s"',sources,cvs_sources{i});
end
for i=1:length(shr_sources)
    sources = sprintf('%s "%s"',sources,shr_sources{i});
end

sun_incdir = fullfile(sun,'include');      % for SUNDIALS exported headers
cvs_srcdir = fullfile(sun,'src','cvodes'); % for cvodes_impl.h
includes = sprintf('%s -I"%s" -I"%s"',includes,sun_incdir,cvs_srcdir);

% Add NVEC_SER sources and header files

nvs_sources = fullfile(sun,'src','nvec_ser','nvector_serial.c');
sources = sprintf('%s "%s"',sources, nvs_sources);

if par
  
  % Add NVEC_PAR sources and header files

  nvp_sources = fullfile(sun,'src','nvec_par','nvector_parallel.c');
  sources = sprintf('%s "%s"',sources, nvp_sources);

% Add LAM headers and libraries

  lam = getenv('LAMHOME');
  lam_incdir = fullfile(lam, 'include');
  lam_libdir = fullfile(lam, 'lib');
  includes = sprintf('%s -I"%s"',includes,lam_incdir);
  libraries = sprintf('%s -L"%s" -lmpi -llam -lutil',libraries,lam_libdir);

end

% Create MEX file

cvm_dir = fullfile(stb,'cvodes','cvm');
cd(cvm_dir)
mex_cmd = sprintf('%s -v %s %s %s', mexcompiler, includes, sources, libraries);
disp(mex_cmd);
eval(mex_cmd);

% Move back to sundialsTB

cd(stb)

%---------------------------------------------------------------------------------
% compilation of idm MEX file
%---------------------------------------------------------------------------------

function [] = compile_IDM(mexcompiler,stb,sun,par)

idm_sources = {
    fullfile(stb,'idas','idm','src','idm.c')
    fullfile(stb,'idas','idm','src','idmWrap.c')
    fullfile(stb,'idas','idm','src','idmOpts.c')
              };
if par
    nvm_sources = {
        fullfile(stb,'nvector','src','nvm_parallel.c')
        fullfile(stb,'nvector','src','nvm_ops.c')
                  };
else
    nvm_sources = {
        fullfile(stb,'nvector','src','nvm_serial.c')
        fullfile(stb,'nvector','src','nvm_ops.c')
                  };
end
sources = '';
for i=1:length(idm_sources)
    sources = sprintf('%s "%s"',sources,idm_sources{i});
end
for i=1:length(nvm_sources)
    sources = sprintf('%s "%s"',sources,nvm_sources{i});
end

idm_incdir = fullfile(stb,'idas','idm','src');   % for idm.h
nvm_incdir = fullfile(stb,'nvector','src');      % for nvm.h
includes = sprintf('-I"%s" -I"%s" -I"%s"',stb,idm_incdir,nvm_incdir);

libraries = '';

% Add IDAS sources and header files

ids_sources = {
    fullfile(sun,'src','ida','ida_band.c')
    fullfile(sun,'src','ida','ida_bbdpre.c')
    fullfile(sun,'src','ida','ida_dense.c')
    fullfile(sun,'src','ida','ida_direct.c')
%    fullfile(sun,'src','idas','idaa.c')
    fullfile(sun,'src','ida','ida.c')
    fullfile(sun,'src','ida','ida_ic.c')
    fullfile(sun,'src','ida','ida_io.c')
%    fullfile(sun,'src','idas','idaa_io.c')
    fullfile(sun,'src','ida','ida_spils.c')
    fullfile(sun,'src','ida','ida_spbcgs.c')
    fullfile(sun,'src','ida','ida_spgmr.c')
    fullfile(sun,'src','ida','ida_sptfqmr.c')
              };
shr_sources = {
    fullfile(sun,'src','sundials','sundials_band.c')
    fullfile(sun,'src','sundials','sundials_dense.c')
    fullfile(sun,'src','sundials','sundials_iterative.c')
    fullfile(sun,'src','sundials','sundials_nvector.c')
    fullfile(sun,'src','sundials','sundials_direct.c')
    fullfile(sun,'src','sundials','sundials_spbcgs.c')
    fullfile(sun,'src','sundials','sundials_spgmr.c')
    fullfile(sun,'src','sundials','sundials_sptfqmr.c')
    fullfile(sun,'src','sundials','sundials_math.c')
              };
for i=1:length(ids_sources)
    sources = sprintf('%s "%s"',sources,ids_sources{i});
end
for i=1:length(shr_sources)
    sources = sprintf('%s "%s"',sources,shr_sources{i});
end

sun_incdir = fullfile(sun,'include');     % for SUNDIALS exported headers
ids_srcdir = fullfile(sun,'src','ida');   % for idas_impl.h
includes = sprintf('%s -I"%s" -I"%s"',includes,sun_incdir,ids_srcdir);

% Add NVEC_SER sources and header files

nvs_sources = fullfile(sun,'src','nvec_ser','nvector_serial.c');
sources = sprintf('%s "%s"',sources, nvs_sources);

if par
  
  % Add NVEC_PAR sources and header files

  nvp_sources = fullfile(sun,'src','nvec_par','nvector_parallel.c');
  sources = sprintf('%s "%s"',sources, nvp_sources);

% Add LAM headers and libraries

  lam = getenv('LAMHOME');
  lam_incdir = fullfile(lam, 'include');
  lam_libdir = fullfile(lam, 'lib');
  includes = sprintf('%s -I"%s"',includes,lam_incdir);
  libraries = sprintf('%s -L"%s" -lmpi -llam -lutil',libraries,lam_libdir);

end

% Create MEX file

idm_dir = fullfile(stb,'idas','idm');
cd(idm_dir)
mex_cmd = sprintf('%s -v %s %s %s', mexcompiler, includes, sources, libraries);
disp(mex_cmd);
eval(mex_cmd);

% Move back to sundialsTB

cd(stb)


%---------------------------------------------------------------------------------
% compilation of KINSOL MEX file
%---------------------------------------------------------------------------------

function [] = compile_KIM(mexcompiler,stb,sun,par)

kim_sources = {
    fullfile(stb,'kinsol','kim','src','kim.c')
    fullfile(stb,'kinsol','kim','src','kimWrap.c')
    fullfile(stb,'kinsol','kim','src','kimOpts.c')
              };
if par
    nvm_sources = {
        fullfile(stb,'nvector','src','nvm_parallel.c')
        fullfile(stb,'nvector','src','nvm_ops.c')
                  };
else
    nvm_sources = {
        fullfile(stb,'nvector','src','nvm_serial.c')
        fullfile(stb,'nvector','src','nvm_ops.c')
                  };
end
sources = '';
for i=1:length(kim_sources)
    sources = sprintf('%s "%s"',sources,kim_sources{i});
end
for i=1:length(nvm_sources)
    sources = sprintf('%s "%s"',sources,nvm_sources{i});
end

kim_incdir = fullfile(stb,'kinsol','kim','src'); % for kim.h
nvm_incdir = fullfile(stb,'nvector','src');      % for nvm.h
includes = sprintf('-I"%s" -I"%s" -I"%s"',stb,kim_incdir,nvm_incdir);

libraries = '';

% Add KINSOL sources and header files

kin_sources = {
    fullfile(sun,'src','kinsol','kinsol_band.c')
    fullfile(sun,'src','kinsol','kinsol_bbdpre.c')
    fullfile(sun,'src','kinsol','kinsol_dense.c')
    fullfile(sun,'src','kinsol','kinsol_direct.c')
    fullfile(sun,'src','kinsol','kinsol.c')
    fullfile(sun,'src','kinsol','kinsol_io.c')
    fullfile(sun,'src','kinsol','kinsol_spils.c')
    fullfile(sun,'src','kinsol','kinsol_spbcgs.c')
    fullfile(sun,'src','kinsol','kinsol_spgmr.c')
    fullfile(sun,'src','kinsol','kinsol_sptfqmr.c')
              };
shr_sources = {
    fullfile(sun,'src','sundials','sundials_band.c')
    fullfile(sun,'src','sundials','sundials_dense.c')
    fullfile(sun,'src','sundials','sundials_iterative.c')
    fullfile(sun,'src','sundials','sundials_nvector.c')
    fullfile(sun,'src','sundials','sundials_direct.c')
    fullfile(sun,'src','sundials','sundials_spbcgs.c')
    fullfile(sun,'src','sundials','sundials_spgmr.c')
    fullfile(sun,'src','sundials','sundials_sptfqmr.c')
    fullfile(sun,'src','sundials','sundials_math.c')
              };

for i=1:length(kin_sources)
    sources = sprintf('%s "%s"',sources,kin_sources{i});
end
for i=1:length(shr_sources)
    sources = sprintf('%s "%s"',sources,shr_sources{i});
end

sun_incdir = fullfile(sun,'include');      % for SUNDIALS exported headers
kin_srcdir = fullfile(sun,'src','kinsol'); % for kinsol_impl.h
includes = sprintf('%s -I"%s" -I"%s"',includes,sun_incdir,kin_srcdir);

% Add NVEC_SER sources and header files

nvs_sources = fullfile(sun,'src','nvec_ser','nvector_serial.c');
sources = sprintf('%s "%s"',sources, nvs_sources);

if par
  
  % Add NVEC_PAR sources and header files

  nvp_sources = fullfile(sun,'src','nvec_par','nvector_parallel.c');
  sources = sprintf('%s "%s"',sources, nvp_sources);

% Add LAM headers and libraries

  lam = getenv('LAMHOME');
  lam_incdir = fullfile(lam, 'include');
  lam_libdir = fullfile(lam, 'lib');
  includes = sprintf('%s -I"%s"',includes,lam_incdir);
  libraries = sprintf('%s -L"%s" -lmpi -llam -lutil',libraries,lam_libdir);

end

% Create MEX file

kim_dir = fullfile(stb, 'kinsol', 'kim');
cd(kim_dir)
mex_cmd = sprintf('%s -v %s %s %s', mexcompiler, includes, sources, libraries);
disp(mex_cmd);
eval(mex_cmd);

% Move back to sundialsTB

cd(stb)

%---------------------------------------------------------------------------------
% Installation of common sundialsTB files
%---------------------------------------------------------------------------------

function [] = instSTB(stb, where, par)

stbi = fullfile(where,'sundialsTB');

% Create startup_STB.m (use the template startup_STB.in)

in_file = fullfile(stb,'startup_STB.in');
fi = fopen(in_file,'r');
out_file = fullfile(stbi,'startup_STB.m');
fo = fopen(out_file,'w');
while(~feof(fi))
  l = fgets(fi);
  i = strfind(l,'@STB_PATH@');
  if ~isempty(i)
    l = sprintf('stb_path = ''%s'';\n',where);
  end
  fprintf(fo,'%s',l);
end
fclose(fo);
fclose(fi);

nvm_files = {
    fullfile('nvector','Contents.m')
    fullfile('nvector','N_VDotProd.m')
    fullfile('nvector','N_VL1Norm.m')
    fullfile('nvector','N_VMax.m')
    fullfile('nvector','N_VMaxNorm.m')
    fullfile('nvector','N_VMin.m')
    fullfile('nvector','N_VWL2Norm.m')
    fullfile('nvector','N_VWrmsNorm.m')
            };

put_files = {
    fullfile('putils','Contents.m')
    fullfile('putils','mpistart.m')
    fullfile('putils','mpirun.m')
    fullfile('putils','mpiruns.m')    
            };

stb_files = [nvm_files];
if par
  stb_files = [stb_files; put_files];
end

fprintf('\n\n');
for i=1:length(stb_files)
  src = fullfile(stb,stb_files{i});
  dest = fullfile(stbi,stb_files{i});
  fprintf('Install %s\n',dest);
  [success,msg,msgid] = copyfile(src,dest);
  if ~success
    disp(msg);
    break;
  end
end

%---------------------------------------------------------------------------------
% Installation of CVODES files
%---------------------------------------------------------------------------------

function [] = instCVM(stb, where, par)

stbi = fullfile(where,'sundialsTB');

% Copy files to installation directory

cvmmex = ['cvm.' mexext];

cvm_files = {
    'LICENSE'
    'Contents.m'
    fullfile('cvodes','Contents.m')
    fullfile('cvodes','CVodeSetOptions.m')
    fullfile('cvodes','CVodeQuadSetOptions.m')
    fullfile('cvodes','CVodeSensSetOptions.m')
    fullfile('cvodes','CVodeInit.m')
    fullfile('cvodes','CVodeReInit.m')
    fullfile('cvodes','CVodeQuadInit.m')
    fullfile('cvodes','CVodeQuadReInit.m')
    fullfile('cvodes','CVodeSensInit.m')
    fullfile('cvodes','CVodeSensReInit.m')
    fullfile('cvodes','CVodeAdjInit.m')
    fullfile('cvodes','CVodeAdjReInit.m')
    fullfile('cvodes','CVodeInitB.m')
    fullfile('cvodes','CVodeReInitB.m')
    fullfile('cvodes','CVodeQuadInitB.m')
    fullfile('cvodes','CVodeQuadReInitB.m')
    fullfile('cvodes','CVodeSensToggleOff.m')
    fullfile('cvodes','CVode.m')
    fullfile('cvodes','CVodeB.m')
    fullfile('cvodes','CVodeGet.m')
    fullfile('cvodes','CVodeGetStatsB.m')
    fullfile('cvodes','CVodeGetStats.m')
    fullfile('cvodes','CVodeFree.m')
    fullfile('cvodes','CVodeMonitor.m')
    fullfile('cvodes','CVodeMonitorB.m')
    fullfile('cvodes','cvm','Contents.m')
    fullfile('cvodes','cvm','cvm_rhs.m')
    fullfile('cvodes','cvm','cvm_rhsQ.m')
    fullfile('cvodes','cvm','cvm_rhsS.m')
    fullfile('cvodes','cvm','cvm_root.m')
    fullfile('cvodes','cvm','cvm_bjac.m')
    fullfile('cvodes','cvm','cvm_djac.m')
    fullfile('cvodes','cvm','cvm_jtv.m')
    fullfile('cvodes','cvm','cvm_pset.m')
    fullfile('cvodes','cvm','cvm_psol.m')
    fullfile('cvodes','cvm','cvm_gcom.m')
    fullfile('cvodes','cvm','cvm_gloc.m')
    fullfile('cvodes','cvm','cvm_rhsB.m')
    fullfile('cvodes','cvm','cvm_rhsQB.m')
    fullfile('cvodes','cvm','cvm_bjacB.m')
    fullfile('cvodes','cvm','cvm_djacB.m')
    fullfile('cvodes','cvm','cvm_jtvB.m')
    fullfile('cvodes','cvm','cvm_psetB.m')
    fullfile('cvodes','cvm','cvm_psolB.m')
    fullfile('cvodes','cvm','cvm_gcomB.m')
    fullfile('cvodes','cvm','cvm_glocB.m')
    fullfile('cvodes','cvm','cvm_monitor.m')
    fullfile('cvodes','cvm','cvm_monitorB.m')
    fullfile('cvodes','cvm','cvm_options.m')    
    fullfile('cvodes','cvm',cvmmex)    
            };

cvm_ftypes = {
    fullfile('cvodes','function_types','CVRhsFn.m')
    fullfile('cvodes','function_types','CVQuadRhsFn.m')
    fullfile('cvodes','function_types','CVSensRhsFn.m')
    fullfile('cvodes','function_types','CVRootFn.m')
    fullfile('cvodes','function_types','CVBandJacFn.m')
    fullfile('cvodes','function_types','CVDenseJacFn.m')
    fullfile('cvodes','function_types','CVJacTimesVecFn.m')
    fullfile('cvodes','function_types','CVPrecSetupFn.m')
    fullfile('cvodes','function_types','CVPrecSolveFn.m')
    fullfile('cvodes','function_types','CVGcommFn.m')
    fullfile('cvodes','function_types','CVGlocalFn.m')
    fullfile('cvodes','function_types','CVMonitorFn.m')
    fullfile('cvodes','function_types','CVRhsFnB.m')
    fullfile('cvodes','function_types','CVQuadRhsFnB.m')
    fullfile('cvodes','function_types','CVBandJacFnB.m')
    fullfile('cvodes','function_types','CVDenseJacFnB.m')
    fullfile('cvodes','function_types','CVJacTimesVecFnB.m')
    fullfile('cvodes','function_types','CVPrecSetupFnB.m')
    fullfile('cvodes','function_types','CVPrecSolveFnB.m')
    fullfile('cvodes','function_types','CVGcommFnB.m')
    fullfile('cvodes','function_types','CVGlocalFnB.m')
    fullfile('cvodes','function_types','CVMonitorFnB.m')
             };    

cvm_exs = {
    fullfile('cvodes','examples_ser','cvadx.m')
    fullfile('cvodes','examples_ser','cvbx_f.m')
    fullfile('cvodes','examples_ser','cvbx_J.m')
    fullfile('cvodes','examples_ser','cvbx.m')
    fullfile('cvodes','examples_ser','cvbx_q.m')
    fullfile('cvodes','examples_ser','cvdiscx.m')
    fullfile('cvodes','examples_ser','cvdx_fB.m')
    fullfile('cvodes','examples_ser','cvdx_f.m')
    fullfile('cvodes','examples_ser','cvdx_fS.m')
    fullfile('cvodes','examples_ser','cvdx_g.m')
    fullfile('cvodes','examples_ser','cvdx_JB.m')
    fullfile('cvodes','examples_ser','cvdx_J.m')
    fullfile('cvodes','examples_ser','cvdx.m')
    fullfile('cvodes','examples_ser','cvdx_qB.m')
    fullfile('cvodes','examples_ser','cvdx_q.m')
    fullfile('cvodes','examples_ser','cvfdx.m')
    fullfile('cvodes','examples_ser','cvkxb.m')
    fullfile('cvodes','examples_ser','cvkx_f.m')
    fullfile('cvodes','examples_ser','cvkx.m')
    fullfile('cvodes','examples_ser','cvkx_pset.m')
    fullfile('cvodes','examples_ser','cvkx_psol.m')
    fullfile('cvodes','examples_ser','pleiades_f.m')
    fullfile('cvodes','examples_ser','pleiades_J.m')
    fullfile('cvodes','examples_ser','pleiades.m')
    fullfile('cvodes','examples_ser','vdp_f.m')
    fullfile('cvodes','examples_ser','vdp_J.m')
    fullfile('cvodes','examples_ser','vdp.m')
    fullfile('cvodes','examples_ser','cvhess.m')
    fullfile('cvodes','examples_ser','cvhess_f.m')
    fullfile('cvodes','examples_ser','cvhess_q.m')
    fullfile('cvodes','examples_ser','cvhess_fB.m')
    fullfile('cvodes','examples_ser','cvhess_qB.m')          
          };

cvm_exp = {
    fullfile('cvodes','examples_par','pvfnx_f.m')
    fullfile('cvodes','examples_par','pvfnx.m')
    fullfile('cvodes','examples_par','pvkx_fl.m')
    fullfile('cvodes','examples_par','pvkx_f.m')
    fullfile('cvodes','examples_par','pvkx.m')
    fullfile('cvodes','examples_par','pvnx_f.m')
    fullfile('cvodes','examples_par','pvnx.m')    
          };

stb_files = [cvm_files ; cvm_ftypes ; cvm_exs];
if par
  stb_files = [stb_files ; cvm_exp];
end

fprintf('\n\n');
for i=1:length(stb_files)
  src = fullfile(stb,stb_files{i});
  dest = fullfile(stbi,stb_files{i});
  fprintf('Install %s\n',dest);
  [success,msg,msgid] = copyfile(src,dest);
  if ~success
    disp(msg);
    break;
  end
end

%---------------------------------------------------------------------------------
% Installation of IDAS files
%---------------------------------------------------------------------------------

function [] = instIDM(stb, where, par)

stbi = fullfile(where,'sundialsTB');

% Copy files to installation directory

idmmex = ['idm.' mexext];

idm_files = {
    fullfile('idas','Contents.m')
    fullfile('idas','IDACalcIC.m')
    fullfile('idas','IDAFree.m')
    fullfile('idas','IDAGet.m')
    fullfile('idas','IDAGetStats.m')
    fullfile('idas','IDASolve.m')
    fullfile('idas','IDAMalloc.m')
    fullfile('idas','IDAMonitor.m')
    fullfile('idas','IDASetOptions.m')
    fullfile('idas','IDAReInit.m')
    fullfile('idas','idm','Contents.m')
    fullfile('idas','idm','idm_bjac.m')
    fullfile('idas','idm','idm_djac.m')
    fullfile('idas','idm','idm_gcom.m')
    fullfile('idas','idm','idm_gloc.m')
    fullfile('idas','idm','idm_jtv.m')
    fullfile('idas','idm','idm_monitor.m')
    fullfile('idas','idm','idm_pset.m')
    fullfile('idas','idm','idm_psol.m')
    fullfile('idas','idm','idm_res.m')
    fullfile('idas','idm','idm_root.m')
    fullfile('idas','idm',idmmex)    
            };

idm_ftypes = {
    fullfile('idas','function_types','IDAResFn.m')
    fullfile('idas','function_types','IDARootFn.m')
    fullfile('idas','function_types','IDADenseJacFn.m')
    fullfile('idas','function_types','IDABandJacFn.m')
    fullfile('idas','function_types','IDAJacTimesVecFn.m')
    fullfile('idas','function_types','IDAGcommFn.m')
    fullfile('idas','function_types','IDAGlocalFn.m')
    fullfile('idas','function_types','IDAPrecSetupFn.m')
    fullfile('idas','function_types','IDAPrecSolveFn.m')
    fullfile('idas','function_types','IDAMonitorFn.m')
             };

idm_exs = {
    fullfile('idas','examples_ser','idabanx.m')
    fullfile('idas','examples_ser','idabanx_ic.m')
    fullfile('idas','examples_ser','idabanx_f.m')
    fullfile('idas','examples_ser','idadenx.m')
    fullfile('idas','examples_ser','idadenx_f.m')
    fullfile('idas','examples_ser','idadenx_J.m')
    fullfile('idas','examples_ser','idadenx_g.m')
    fullfile('idas','examples_ser','pend.m')
    fullfile('idas','examples_ser','pendGGL.m')
          };

stb_files = [idm_files ; idm_ftypes ; idm_exs];
%if par
%  stb_files = [stb_files ; idm_exp];
%end

fprintf('\n\n');
for i=1:length(stb_files)
  src = fullfile(stb,stb_files{i});
  dest = fullfile(stbi,stb_files{i});
  fprintf('Install %s\n',dest);
  [success,msg,msgid] = copyfile(src,dest);
  if ~success
    disp(msg);
    break;
  end
end

%---------------------------------------------------------------------------------
% Installation of KINSOL files
%---------------------------------------------------------------------------------

function [] = instKIM(stb, where, par)

stbi = fullfile(where,'sundialsTB');

% Copy files to installation directory

kimmex = ['kim.' mexext];

kim_files = {
    fullfile('kinsol','Contents.m')
    fullfile('kinsol','KINFree.m')
    fullfile('kinsol','KINGetStats.m')
    fullfile('kinsol','KINMalloc.m')
    fullfile('kinsol','KINSetOptions.m')
    fullfile('kinsol','KINSol.m')
    fullfile('kinsol','kim','Contents.m')
    fullfile('kinsol','kim','kim_bjac.m')
    fullfile('kinsol','kim','kim_djac.m')
    fullfile('kinsol','kim','kim_gcom.m')
    fullfile('kinsol','kim','kim_gloc.m')
    fullfile('kinsol','kim','kim_info.m')
    fullfile('kinsol','kim','kim_jtv.m')
    fullfile('kinsol','kim','kim_pset.m')
    fullfile('kinsol','kim','kim_psol.m')
    fullfile('kinsol','kim','kim_sys.m')
    fullfile('kinsol','kim',kimmex)
            };

kim_ftypes = {
    fullfile('kinsol','function_types','KINSysFn.m')
    fullfile('kinsol','function_types','KINBandJacFn.m')
    fullfile('kinsol','function_types','KINDenseJacFn.m')
    fullfile('kinsol','function_types','KINJacTimesVecFn.m')
    fullfile('kinsol','function_types','KINPrecSetupFn.m')
    fullfile('kinsol','function_types','KINPrecSolveFn.m')
    fullfile('kinsol','function_types','KINGcommFn.m')
    fullfile('kinsol','function_types','KINGlocalFn.m')
             };

kim_exs = {
    fullfile('kinsol','examples_ser','kindiag.m')
    fullfile('kinsol','examples_ser','kindiag_pset.m')
    fullfile('kinsol','examples_ser','kindiag_psol.m')
    fullfile('kinsol','examples_ser','kindiag_sys.m')
    fullfile('kinsol','examples_ser','kindx.m')
    fullfile('kinsol','examples_ser','kindx_sys.m')
          };

kim_exp = {
    fullfile('kinsol','examples_par','kindiagp.m')
    fullfile('kinsol','examples_par','kindiagp_pset.m')
    fullfile('kinsol','examples_par','kindiagp_psol.m')
    fullfile('kinsol','examples_par','kindiagp_sys.m') 
          };

stb_files = [kim_files ; kim_ftypes ; kim_exs];
if par
  stb_files = [stb_files ; kim_exp];
end

fprintf('\n\n');
for i=1:length(stb_files)
  src = fullfile(stb,stb_files{i});
  dest = fullfile(stbi,stb_files{i});
  fprintf('Install %s\n',dest);
  [success,msg,msgid] = copyfile(src,dest);
  if ~success
    disp(msg);
    break;
  end
end


