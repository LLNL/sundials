function [] = texdoc
% TEXDOC - Creates LaTeX documentation for sundialsTB
%

% Radu Serban <radu@llnl.gov>
% Copyright (c) 2005, The Regents of the University of California.
% $Revision: 1.4 $Date: 2006/03/14 00:47:29 $


global cv_fct cv_ftp ces_drv ces_fct cep_drv cep_fct
global kin_fct kin_ftp kes_drv kes_fct kep_drv kep_fct
global nvec_fct putils_fct

% If the output directory does not exist, create it
system('mkdir -p tex');

% Set output dir
cd('..');
doc_dir = 'doc/tex';

% Set list of files to be processed
set_files;

%==================================================================      

file = sprintf('%s/sundialsTB.tex',doc_dir);
fid = fopen(file,'w');

fprintf(fid,'\\input{front}\n\n');

%==================================================================      

fprintf(fid,'\\input{cvodes_top}\n\n');

% CVODES interface functions

fprintf(fid,'\\newpage\n\\subsection{Interface functions}\n');
fprintf(fid,'\\rule{0in}{0.25in}\n\n');
for i = 1:length(cv_fct)
  process_no_source(cv_fct{i}, doc_dir, fid);
end
process_with_source('cvodes/CVodeMonitor.m', doc_dir, fid);

% CVODES function types

fprintf(fid,'\n\\newpage\n\\subsection{Function types}\n'); 
fprintf(fid,'\\rule{0in}{0.25in}\n\n');
for i = 1:length(cv_ftp)
  process_no_source(cv_ftp{i}, doc_dir, fid);
end

%==================================================================

fprintf(fid,'\\input{kinsol_top}\n\n');

% KINSOL interface functions

fprintf(fid,'\\newpage\n\\subsection{Interface functions}\n');
fprintf(fid,'\\rule{0in}{0.25in}\n\n');
for i = 1:length(kin_fct)
  process_no_source(kin_fct{i}, doc_dir, fid);
end

% KINSOL function types

fprintf(fid,'\n\\newpage\n\\subsection{Function types}\n'); 
fprintf(fid,'\\rule{0in}{0.25in}\n\n');
for i = 1:length(kin_ftp)
  process_no_source(kin_ftp{i}, doc_dir, fid);
end

%==================================================================

fprintf(fid,'\n\n\\input{other_top}\n\n');

fprintf(fid,'\n\n\\newpage\n\\subsection{{\\nvector} functions}\n');
fprintf(fid,'\\rule{0in}{0.25in}\n\n');
for i = 1:length(nvec_fct)
  process_with_source(nvec_fct{i}, doc_dir, fid);
end

fprintf(fid,'\n\n\\newpage\n\\subsection{Parallel utilities}\n');
fprintf(fid,'\\rule{0in}{0.25in}\n\n');
for i = 1:length(putils_fct)
  process_with_source(putils_fct{i}, doc_dir, fid);
end

%==================================================================   

fprintf(fid,'\n\\include{references}\n\n\\end{document}');
fclose(fid);

system('cp doc/tex_files/*.tex doc/tex/');
system('cp doc/tex_files/biblio.bib doc/tex/');
system('cp doc/tex_files/doc_logo.eps doc/tex/');


cd(doc_dir);

fprintf('Compile LaTeX files...\n');
system('latex sundialsTB');
system('bibtex sundialsTB');
system('latex sundialsTB');
system('latex sundialsTB');

fprintf('Generate PS file...\n');
system('dvips -o sundialsTB.ps sundialsTB');

fprintf('Generate PDF file...\n');
system('dvipdf sundialsTB');

system('cp sundialsTB.ps ..');
system('cp sundialsTB.pdf ..');

cd('..');

%==================================================================   
% Process matlab file. Do not generate source
%==================================================================   

function process_no_source(name, doc_dir, fid)

m2html('mfiles',name,...
       'verbose','off',...
       'syntaxHighlighting','off',...
       'source','off',...
       'htmldir',doc_dir,...
       'extension','.tex',...
       'indexFile','foo',...
       'template','latex');

fprintf('%s\n',name);

ii = strfind(name,'.m');     
fprintf(fid,'\\input{%s}\n',name(1:ii-1)); 


%==================================================================   
% Process matlab file. Generated highlighted source code.
%==================================================================   

function process_with_source(name, doc_dir, fid)

m2html('mfiles',name,...
       'verbose','off',...
       'syntaxHighlighting','off',...
       'htmldir',doc_dir,...
       'extension','.tex',...
       'source','on', ...
       'indexFile','foo',...
       'template','latex');
  
ii = strfind(name,'.m');
fi = fopen(name,'r');
fo = fopen(sprintf('%s/%ssrc.tex',doc_dir, name(1:ii-1)),'w');

start = 0;
finish = 0;

l = fgets(fi);
if l(1) == '%'
  start = 1;
end
i = 2;
while ~feof(fi)
  l = fgets(fi);
  if l(1) ~= '%'
    if finish == 0, finish = i-1; end
  else
    if start == 0,  start = i;    end
  end  
  i = i + 1;
end
last = i;
frewind(fi);

fprintf('%s   %d  %d  %d\n',name,start,finish,i);


if start == 1
  lines = sprintf('linerange={%d-%d}',finish+1,last);
else
  lines = sprintf('linerange={1-%d,%d-%d}',start-1,finish+1,last);  
end

fprintf(fo,'\\begin{lstlisting}[%s]\n',lines);
for i = 1: last-1
  fprintf(fo,'%s',fgets(fi));
end
fprintf(fo,'\\end{lstlisting}\n');

fclose(fo);
fclose(fi);

fprintf(fid,'\\input{%s}\n',name(1:ii-1)); 

%==================================================================   
% Set files
%==================================================================   

function [] = set_files()

global cv_fct cv_ftp ces_drv ces_fct cep_drv cep_fct
global kin_fct kin_ftp kes_drv kes_fct kep_drv kep_fct
global nvec_fct putils_fct

cv_fct = {
    'cvodes/CVodeSetOptions.m'...
    'cvodes/CVodeSetFSAOptions.m'...
    'cvodes/CVodeMalloc.m'...
    'cvodes/CVodeSensMalloc.m'...
    'cvodes/CVadjMalloc.m'...
    'cvodes/CVodeMallocB.m'...
    'cvodes/CVode.m'...
    'cvodes/CVodeB.m'...
    'cvodes/CVodeGetStats.m'...
    'cvodes/CVodeGetStatsB.m'...
    'cvodes/CVodeGet.m'...
    'cvodes/CVodeFree.m'...
         };

cv_ftp = {
    'cvodes/CVBandJacFn.m'...
    'cvodes/CVDenseJacFn.m'...
    'cvodes/CVGcommFn.m'...
    'cvodes/CVGlocalFn.m'...
    'cvodes/CVMonitorFn.m'...
    'cvodes/CVQuadRhsFn.m'...
    'cvodes/CVRhsFn.m'...
    'cvodes/CVRootFn.m'...
    'cvodes/CVSensRhsFn.m'...
    'cvodes/CVJacTimesVecFn.m'...
    'cvodes/CVPrecSetupFn.m'...
    'cvodes/CVPrecSolveFn.m'...
         };

kin_fct = {
    'kinsol/KINSetOptions.m'...
    'kinsol/KINMalloc.m'...
    'kinsol/KINSol.m'...
    'kinsol/KINGetStats.m'...
    'kinsol/KINFree.m'...
          };

kin_ftp = {
    'kinsol/KINDenseJacFn.m'...
    'kinsol/KINBandJacFn.m'...
    'kinsol/KINGcommFn.m'...
    'kinsol/KINGlocalFn.m'...
    'kinsol/KINJacTimesVecFn.m'...
    'kinsol/KINPrecSetupFn.m'...
    'kinsol/KINPrecSolveFn.m'...
    'kinsol/KINSysFn.m'...
    };

nvec_fct = {
    'nvector/N_VDotProd.m'...
    'nvector/N_VL1Norm.m'...
    'nvector/N_VMax.m'...
    'nvector/N_VMaxNorm.m'...
    'nvector/N_VMin.m'...
    'nvector/N_VWL2Norm.m'...
    'nvector/N_VWrmsNorm.m'...
           };
putils_fct = {
    'putils/mpirun.m'...
    'putils/mpiruns.m'...
    'putils/mpistart.m'...
             };
ces_drv = {
    'cvodes/examples_ser/cvdx.m'...
    'cvodes/examples_ser/cvbx.m'...
    'cvodes/examples_ser/cvkx.m'...
    'cvodes/examples_ser/cvkxb.m'...
    'cvodes/examples_ser/pleiades.m'...
    'cvodes/examples_ser/vdp.m'...
    'cvodes/examples_ser/cvfdx.m'...
    'cvodes/examples_ser/cvadx.m'...
         };

ces_fct ={
    'cvodes/examples_ser/cvdx_f.m'...
    'cvodes/examples_ser/cvdx_fB.m'...
    'cvodes/examples_ser/cvdx_J.m'...
    'cvodes/examples_ser/cvdx_JB.m'...
    'cvodes/examples_ser/cvdx_q.m'...
    'cvodes/examples_ser/cvdx_qB.m'...
    'cvodes/examples_ser/cvdx_g.m'...
    'cvodes/examples_ser/cvdx_fS.m'...
    'cvodes/examples_ser/cvbx_f.m'...
    'cvodes/examples_ser/cvbx_J.m'...
    'cvodes/examples_ser/cvbx_q.m'...
    'cvodes/examples_ser/cvkx_f.m'...
    'cvodes/examples_ser/cvkx_pset.m'...
    'cvodes/examples_ser/cvkx_psol.m'...
    'cvodes/examples_ser/pleiades_f.m'...
    'cvodes/examples_ser/pleiades_J.m'...
    'cvodes/examples_ser/vdp_f.m'...
    'cvodes/examples_ser/vdp_J.m'...
        };

cep_drv = {
    'cvodes/examples_par/pvnx.m'...
    'cvodes/examples_par/pvfnx.m'...
         };

cep_drv ={
    'cvodes/examples_par/pvnx_f.m'...
    'cvodes/examples_par/pvfnx_f.m'...
         };


