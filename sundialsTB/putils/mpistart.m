function mpistart(nslaves, rpi, hosts)
%MPISTART invokes lamboot (if required) and MPI_Init (if required).
%
%   Usage: MPISTART [ ( NSLAVES [, RPI [, HOSTS] ] ) ]
%
%   MPISTART boots LAM and initializes MPI to match a given number of slave 
%   hosts (and rpi) from a given list of hosts. All three args optional.
%
%   If they are not defined, HOSTS are taken from a builtin HOSTS list
%   (edit HOSTS at the beginning of this file to match your cluster)
%   or from the bhost file if defined through LAMBHOST (in this order).
%
%   If not defined, RPI is taken from the builtin variable RPI (edit it
%   to suit your needs) or from the LAM_MPI_SSI_rpi environment variable
%   (in this order).

% Heavily based on the LAM_Init function in MPITB.

%------------------------------
% ARGCHECK
%------------------------------

% List of hosts

if nargin>2

% Hosts passed as an argument...

  if ~iscell(hosts)
    error('MPISTART: 3rd arg is not a cell');
  end        
  for i=1:length(hosts)
    if ~ischar(hosts{i})
      error('MPISTART: 3rd arg is not cell-of-strings');
    end
  end

else

% Get hosts from file specified in env. var. LAMBHOST

  bfile = getenv('LAMBHOST');
  if isempty(bfile)
    error('MPISTART: cannot find list of hosts');
  end
  hosts = readHosts(bfile);

end

% RPI

if nargin>1

% RPI passed as an argument

  if ~ischar(rpi)
    error('MPISTART: 2nd arg is not a string')
  end
% Get full rpi name, if single letter used
  rpi = rpi_str(rpi);
  if isempty(rpi)
    error('MPISTART: 2nd arg is not a known RPI')
  end

else

% Get RPI from env. var. LAM_MPI_SSI_rpi

  RPI = getenv('LAM_MPI_SSI_rpi');
  if isempty(RPI)
%   If LAM_MPI_SSI_rpi not defined, use RPI='tcp'
    RPI = 'tcp';
  end
  rpi = rpi_str(RPI);

end

% Number of slaves

if nargin>0
  if ~isreal(nslaves) || fix(nslaves)~=nslaves  || nslaves>=length(hosts)
    error('MPISTART: 1st arg is not a valid #slaves')
  end
else
  nslaves = length(hosts)-1;
end

%------------------------------
% LAMHALT %
%-------------------------------------------------------------
% reasons to lamhalt:
% - not enough nodes  (nslv+1) % NHL < NSLAVES+1
% - localhost not      in list % weird - just lamboot (NHL=0)
% - localhost not last in list % weird - just lamboot (NHL=0)
%-------------------------------------------------------------

% Lam Nodes Output
[stat, LNO] = system('lamnodes');
if ~stat                                      % already lambooted
  
  emptyflag = false;
  if isempty(LNO)
    % this shouldn't happen
    emptyflag=true;
    % it's MATLAB's fault I think
    fprintf('pushing stubborn MATLAB "system" call (lamnodes): ');
  end
  
  while isempty(LNO) || stat
    fprintf('.');
    [stat, LNO] = system('lamnodes');
  end
  if emptyflag
    fprintf('\n'); 
  end 

  LF = char(10);
  LNO = split(LNO,LF);                        % split lines in rows at \n
  
  [stat, NHL] = system('lamnodes|wc -l');     % Number of Hosts in Lamnodes
        
  emptyflag = false;                          % again,
  if isempty(NHL)                             % this shouldn't happen
    emptyflag=true;                           % it's MATLAB's fault I think
    fprintf('pushing stubborn MATLAB "system" call (lamnodes|wc): ');
  end
  while isempty(NHL) || stat
    fprintf('.');
    [stat, NHL] = system('lamnodes|wc -l');
  end
  if emptyflag
    fprintf('\n'); 
  end
  
  NHL = str2num(NHL);
  if NHL ~= size(LNO,1) || ~ NHL>0            % Oh my, logic error
    NHL= 0;                                   % pretend there are no nodes
    disp('MPISTART: internal logic error: lamboot')
  end                                         % to force lamboot w/o lamhalt
  if isempty(findstr(LNO(end,:),'this_node')) % master computer last in list
    disp('MPISTART: local host is not last in nodelist, hope that''s right')
    beforeflag=0;
    for i=1:size(LNO,1)
      if ~isempty(findstr(LNO(i,:),'this_node'))
        beforeflag=1; 
        break;                                % well, not 1st but it's there
      end
    end                                       % we already warned the user
    if ~beforeflag                            % Oh my, incredible, not there
      NHL= 0;                                 % pretend there are no nodes
      disp('MPISTART: local host not in LAM? lamboot')
    end
  end                                         % to force lamboot w/o lamhalt

  if NHL > 0                                  % accurately account multiprocessors
    NCL = 0;                                  % number of CPUs in lamnodes
    for i=1:size(LNO,1)                       % add the 2nd ":"-separated
      fields=split(LNO(i,:),':');             % field, ie, #CPUs
      NCL = NCL + str2num(fields(2,:));
    end
    if NCL<NHL                                % Oh my, logic error
      NHL= 0;                                 % pretend there are no nodes
      disp('MPISTART: internal logic error: lamboot')
    else
      % update count
      NHL=NCL;                                
    end                                       % can't get count from MPI, 
  end                                         % since might be not _Init'ed

  if NHL < nslaves+1                          % we have to lamboot

    % but avoid getting caught
    [infI flgI]=MPI_Initialized;              % Init?
    [infF flgF]=MPI_Finalized;                % Finalize?
    if infI ||  infF
      error('MPISTART: error calling _Initialized/_Finalized?')
    end
    if flgI && ~flgF                          % avoid hangup due to
      MPI_Finalize;                           % inminent lamhalt
      clear MPI_*                             % force MPI_Init in Mast/Ping
      disp('MPISTART: MPI already used- clearing before lamboot')
    end                                       % by pretending "not _Init"
    if NHL > 0                                % avoid lamhalt in weird cases
      disp('MPISTART: halting LAM')
      system('lamhalt');                      % won't get caught on this
    end
  end
end

%------------------------------
% LAMBOOT
%-------------------------------------------------------------
% reasons to lamboot:          %
% - not lambooted yet          % stat~=0
% - lamhalted above (or weird) % NHL < NSLAVES+1 (0 _is_ <)
%-------------------------------------------------------------

if stat || NHL<nslaves+1

  HNAMS=hosts{end};
  for i=nslaves:-1:1
    HNAMS=strvcat(hosts{i},HNAMS); 
  end
  HNAMS = HNAMS';                             % transpose for "for"

  fid=fopen('bhost','wt');
  for h = HNAMS
    fprintf(fid,'%s\n',h');                   % write slaves' hostnames
  end        
  fclose(fid);
  disp  ('MPISTART: booting LAM')

  stat = system('lamboot -s -v bhost');

  if stat                                     % again, this shouldn't happen
    fprintf('pushing stubborn MATLAB "system" call (lamboot): ');
    while stat
      fprintf('.'); stat = system('lamboot -s -v bhost');
    end
    fprintf('\n');
  end

  system('rm -f bhost');                      % don't need bhost anymore
end                                           % won't wipe on exit/could lamhalt

%------------------------------
% RPI CHECK
%------------------------------

[infI flgI] = MPI_Initialized;                % Init?
[infF flgF] = MPI_Finalized;                  % Finalize?

if infI || infF
  error('MPISTART: error calling _Initialized/_Finalized?')
end
 
if  flgI && ~flgF                             % Perfect, ready to start
else                                          % something we could fix?
  if flgI ||  flgF                            % MPI used, will break
    clear MPI_*                               % unless we clear MPITB
    disp('MPISTART: MPI already used- clearing')  % must start over
  end
                                                                                
  MPI_Init;
end

%------------------------------
% NSLAVES CHECK
%------------------------------

[info attr flag] = MPI_Attr_get(MPI_COMM_WORLD,MPI_UNIVERSE_SIZE);
if info | ~flag
  error('MPISTART: attribute MPI_UNIVERSE_SIZE does not exist?')
end
if attr<2
  error('MPISTART: required 2 computers in LAM')
end

%====================================================================

function hosts = readHosts(bfile)

hosts = [];

fid = fopen(bfile);
if fid == -1
  fprintf('Cannot open bhost file %s\n',bfile);
  return;
end

i = 0;
while ~feof(fid)
% get a line
  l = fgetl(fid);
% Discard comments
  ic = min(strfind(l,'#'));
  if ~isempty(ic), l = l(1:ic-1); end
% Test if there is anything left :-)
  if isempty(l), continue; end
% Got a new host  
  i = i + 1;
% Stop at first blank or tab=char(9)
  indx = find((l==' ') | (l==char(9)));
  if isempty(indx)
    hosts{i} = l;
  else
    hosts{i} = l(1:min(indx));
  end
end

fclose(fid);


%====================================================================

function rpi = rpi_str(c)
%RPI_STR Full LAM SSI RPI string given initial letter(s)
%
%  rpi = rpi_str (c)
%
%  c    initial char(s) of rpi name: t,l,u,s
%  rpi  full rpi name, one of: tcp, lamd, usysv, sysv
%       Use '' if c doesn't match to any supported rpi
%
 
flag = nargin~=1 || isempty(c) || ~ischar(c);
if flag
  return
end
                                                                                
c=lower(c(1));
rpis={'tcp','lamd','usysv','sysv','none'};    % 'none' is sentinel

for i=1:length(rpis)
  if rpis{i}(1)==c
    break
  end
end

if i<length(rpis)
  rpi=rpis{i};                                % normal cases
else
  rpi='';                                     % no way, unknown rpi
end
