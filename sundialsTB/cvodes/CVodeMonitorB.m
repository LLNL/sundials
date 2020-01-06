function [new_data] = CVodeMonitorB(call, idxB, T, Y, YQ, data)
%CVodeMonitorB is the default CVODES monitoring function for backward problems.
%   To use it, set the Monitor property in CVodeSetOptions to
%   'CVodeMonitorB' or to @CVodeMonitorB and 'MonitorData' to mondata
%   (defined as a structure).
%  
%   With default settings, this function plots the evolution of the step 
%   size, method order, and various counters.
%   
%   Various properties can be changed from their default values by passing
%   to CVodeSetOptions, through the property 'MonitorData', a structure
%   MONDATA with any of the following fields. If a field is not defined, 
%   the corresponding default value is used.
%
%   Fields in MONDATA structure:
%     o stats [ {true} | false ]
%         If true, report the evolution of the step size and method order.
%     o cntr [ {true} | false ]
%         If true, report the evolution of the following counters:
%         nst, nfe, nni, netf, ncfn (see CVodeGetStats)
%     o mode [ {'graphical'} | 'text' | 'both' ] 
%         In graphical mode, plot the evolutions of the above quantities.
%         In text mode, print a table.
%     o sol  [ true | {false} ]
%         If true, plot solution components.
%     o select [ array of integers ]
%         To plot only particular solution components, specify their indeces in
%         the field select. If not defined, but sol=true, all components are plotted.
%     o updt [ integer | {50} ]
%         Update frequency. Data is posted in blocks of dimension n.
%     o skip [ integer | {0} ]
%         Number of integrations steps to skip in collecting data to post.
%     o post [ {true} | false ]
%         If false, disable all posting. This option is necessary to disable
%         monitoring on some processors when running in parallel.
%
%   See also CVodeSetOptions, CVMonitorFnB
%
%   NOTES:
%     1. The argument mondata is REQUIRED. Even if only the default options
%        are desired, set mondata=struct; and pass it to CVodeSetOptions.
%     2. The yQ argument is currently ignored.     

% Radu Serban <radu@llnl.gov>
% SUNDIALS Copyright Start
% Copyright (c) 2002-2020, Lawrence Livermore National Security
% and Southern Methodist University.
% All rights reserved.
%
% See the top-level LICENSE and NOTICE files for details.
%
% SPDX-License-Identifier: BSD-3-Clause
% SUNDIALS Copyright End
% $Revision$Date: 2006/10/05 22:12:20 $


if (nargin ~= 6) 
  error('Monitor data not defined.');
end


new_data = [];

if call == 0

% Initialize unspecified fields to default values.
  data = initialize_data(data);  

% Open figure windows
  if data.post

    if data.grph
      if data.stats | data.cntr
        data.hfg = figure;
      end
%     Number of subplots in figure hfg
      if data.stats
        data.npg = data.npg + 2;
      end
      if data.cntr
        data.npg = data.npg + 1;
      end
    end
    
    if data.text
      if data.cntr | data.stats
        data.hft = figure;
      end
    end

    if data.sol
      data.hfs = figure; 
    end
  
  end
  
% Initialize other private data
  data.i = 0;
  data.n = 1;  
  data.t = zeros(1,data.updt);
  if data.stats
    data.h = zeros(1,data.updt);
    data.q = zeros(1,data.updt);
  end
  if data.cntr
    data.nst = zeros(1,data.updt);
    data.nfe = zeros(1,data.updt);
    data.nni = zeros(1,data.updt);
    data.netf = zeros(1,data.updt);
    data.ncfn = zeros(1,data.updt);
  end

  data.first = true;        % the next one will be the first call = 1
  data.initialized = false; % the graphical windows were not initalized
  
  new_data = data;
  
  return;

else

% If this is the first call ~= 0, 
% use Y for additional initializations
  
  if data.first

    if data.sol
      
      if isempty(data.select)
      
        data.N = length(Y);
        data.select = [1:data.N];
        
      else
        
        data.N = length(data.select);
        
      end

      if data.sol
        data.y = zeros(data.N,data.updt);
        data.nps = data.nps + 1;
      end
        
    end
    
    data.first = false;
  
  end
  
% Extract variables from data

  hfg  = data.hfg;
  hft  = data.hft;
  hfs  = data.hfs;
  npg  = data.npg;
  nps  = data.nps;
  i    = data.i;
  n    = data.n;
  t    = data.t;
  N    = data.N;
  y    = data.y;
  h    = data.h;
  q    = data.q;
  nst  = data.nst;
  nfe  = data.nfe;
  nni  = data.nni;
  netf = data.netf;
  ncfn = data.ncfn;
  
end


% Load current statistics?

if call == 1

  if i ~= 0
    i = i-1;
    data.i = i;
    new_data = data;
    return;
  end

  si = CVodeGetStatsB(idxB);

  t(n) = si.tcur;
  
  if data.stats
    h(n) = si.hlast;
    q(n) = si.qlast;
  end
  
  if data.cntr
    nst(n) = si.nst;
    nfe(n) = si.nfe;
    nni(n) = si.nni;
    netf(n) = si.netf;
    ncfn(n) = si.ncfn;
  end

  if data.sol
    for j = 1:N
      y(j,n) = Y(data.select(j));
    end
  end

end

% Is it time to post?

if data.post & (n == data.updt | call==2)

  if call == 2
    n = n-1;
  end
  
  if ~data.initialized

    if (data.stats | data.cntr) & data.grph
      graphical_init(n, hfg, npg, data.stats, data.cntr, ...
                     t, h, q, nst, nfe, nni, netf, ncfn);
    end
    
    if (data.stats | data.cntr) & data.text
      text_init(n, hft, data.stats, data.cntr, ...
                t, h, q, nst, nfe, nni, netf, ncfn);
    end

    if data.sol
      sol_init(n, hfs, nps, data.sol, ...
               N, t, y);
    end
    
    data.initialized = true;
  
  else

    if (data.stats | data.cntr) & data.grph
      graphical_update(n, hfg, npg, data.stats, data.cntr, ...
                       t, h, q, nst, nfe, nni, netf, ncfn);
    end

    if (data.stats | data.cntr) & data.text
      text_update(n, hft, data.stats, data.cntr, ...
                  t, h, q, nst, nfe, nni, netf, ncfn);
    end
    
    if data.sol
      sol_update(n, hfs, nps, data.sol, N, t, y);
    end
      
  end

  if call == 2
    
    if (data.stats | data.cntr) & data.grph
      graphical_final(hfg, npg, data.cntr, data.stats);
    end
    
    if data.sol
      sol_final(hfs, nps, data.sol, N);
    end

    return;
  
  end
  
  n = 1;

else

  n = n + 1;

end


% Save updated values in data

data.i    = data.skip;
data.n    = n;
data.npg  = npg;
data.t    = t;
data.y    = y;
data.h    = h;
data.q    = q;
data.nst  = nst;
data.nfe  = nfe;
data.nni  = nni;
data.netf = netf;
data.ncfn = ncfn;
  
new_data = data;

return;

%-------------------------------------------------------------------------

function data = initialize_data(data)

if ~isfield(data,'mode')
  data.mode = 'graphical';
end
if ~isfield(data,'updt')
  data.updt = 50;
end
if ~isfield(data,'skip')
  data.skip = 0;
end
if ~isfield(data,'stats')
  data.stats = true;
end
if ~isfield(data,'cntr')
  data.cntr = true;
end
if ~isfield(data,'sol')
  data.sol = false;
end
if ~isfield(data,'select')
  data.select = [];
end
if ~isfield(data,'post')
  data.post = true;
end

data.grph = true;
data.text = true;
if strcmp(data.mode,'graphical')
  data.text = false;
end
if strcmp(data.mode,'text')
  data.grph = false;
end

if ~data.sol
  data.select = [];
end
  
% Other initializations
data.npg = 0;
data.nps = 0;
data.hfg = 0;
data.hft = 0;
data.hfs = 0;
data.h = 0;
data.q = 0;
data.nst = 0;
data.nfe = 0;
data.nni = 0;
data.netf = 0;
data.ncfn = 0;
data.N = 0;
data.y = 0;

%-------------------------------------------------------------------------

function [] = graphical_init(n, hfg, npg, stats, cntr, ...
                             t, h, q, nst, nfe, nni, netf, ncfn)

fig_name = 'CVODES run statistics';

% If this is a parallel job, look for the MPI rank in the global
% workspace and append it to the figure name

global sundials_MPI_rank

if ~isempty(sundials_MPI_rank)
  fig_name = sprintf('%s (PE %d)',fig_name,sundials_MPI_rank);
end

figure(hfg);
set(hfg,'Name',fig_name);
set(hfg,'color',[1 1 1]);
pl = 0;

% Time label and figure title

tlab = '\leftarrow   t   \leftarrow';

% Step size and order
if stats
  pl = pl+1;
  subplot(npg,1,pl)
  semilogy(t(1:n),abs(h(1:n)),'-');
  hold on;
  box on;
  grid on;
  xlabel(tlab);
  ylabel('|Step size|');
  
  pl = pl+1;
  subplot(npg,1,pl)
  plot(t(1:n),q(1:n),'-');
  hold on;
  box on;
  grid on;
  xlabel(tlab);
  ylabel('Order');
end
  
% Counters
if cntr
  pl = pl+1;
  subplot(npg,1,pl)
  plot(t(1:n),nst(1:n),'k-');
  hold on;
  plot(t(1:n),nfe(1:n),'b-');
  plot(t(1:n),nni(1:n),'r-');
  plot(t(1:n),netf(1:n),'g-');
  plot(t(1:n),ncfn(1:n),'c-');
  box on;
  grid on;
  xlabel(tlab);
  ylabel('Counters');
end

drawnow;

%-------------------------------------------------------------------------

function [] = graphical_update(n, hfg, npg, stats, cntr, ...
                               t, h, q, nst, nfe, nni, netf, ncfn)

figure(hfg);
pl = 0;
  
% Step size and order
if stats
  pl = pl+1;
  subplot(npg,1,pl)
  hc = get(gca,'Children');
  xd = [get(hc,'XData') t(1:n)];
  yd = [get(hc,'YData') abs(h(1:n))];
  set(hc, 'XData', xd, 'YData', yd);
  
  pl = pl+1;
  subplot(npg,1,pl)
  hc = get(gca,'Children');
  xd = [get(hc,'XData') t(1:n)];
  yd = [get(hc,'YData') q(1:n)];
  set(hc, 'XData', xd, 'YData', yd);
end

% Counters
if cntr
  pl = pl+1;
  subplot(npg,1,pl)
  hc = get(gca,'Children');
% Attention: Children are loaded in reverse order!
  xd = [get(hc(1),'XData') t(1:n)];
  yd = [get(hc(1),'YData') ncfn(1:n)];
  set(hc(1), 'XData', xd, 'YData', yd);
  yd = [get(hc(2),'YData') netf(1:n)];
  set(hc(2), 'XData', xd, 'YData', yd);
  yd = [get(hc(3),'YData') nni(1:n)];
  set(hc(3), 'XData', xd, 'YData', yd);
  yd = [get(hc(4),'YData') nfe(1:n)];
  set(hc(4), 'XData', xd, 'YData', yd);
  yd = [get(hc(5),'YData') nst(1:n)];
  set(hc(5), 'XData', xd, 'YData', yd);
end

drawnow;

%-------------------------------------------------------------------------

function [] = graphical_final(hfg,npg,stats,cntr)

figure(hfg);
pl = 0;

if stats
  pl = pl+1;
  subplot(npg,1,pl)
  hc = get(gca,'Children');
  xd = get(hc,'XData');
  set(gca,'XLim',sort([xd(1) xd(end)]));
  
  pl = pl+1;
  subplot(npg,1,pl)
  ylim = get(gca,'YLim');
  ylim(1) = ylim(1) - 1;
  ylim(2) = ylim(2) + 1;
  set(gca,'YLim',ylim);
  set(gca,'XLim',sort([xd(1) xd(end)]));
end

if cntr
  pl = pl+1;
  subplot(npg,1,pl)
  hc = get(gca,'Children');
  xd = get(hc(1),'XData');
  set(gca,'XLim',sort([xd(1) xd(end)]));
  legend('nst','nfe','nni','netf','ncfn',2);
end

%-------------------------------------------------------------------------

function [] = text_init(n,hft,stats,cntr,t,h,q,nst,nfe,nni,netf,ncfn)

fig_name = 'CVODES run statistics';

% If this is a parallel job, look for the MPI rank in the global
% workspace and append it to the figure name

global sundials_MPI_rank

if ~isempty(sundials_MPI_rank)
  fig_name = sprintf('%s (PE %d)',fig_name,sundials_MPI_rank);
end

figure(hft);
set(hft,'Name',fig_name);
set(hft,'color',[1 1 1]);
set(hft,'MenuBar','none');
set(hft,'Resize','off');

% Create text box

margins=[10 10 50 50]; % left, right, top, bottom
pos=get(hft,'position');
tbpos=[margins(1) margins(4) pos(3)-margins(1)-margins(2) ...
       pos(4)-margins(3)-margins(4)];
tbpos(tbpos<1)=1;

htb=uicontrol(hft,'style','listbox','position',tbpos,'tag','textbox');
set(htb,'BackgroundColor',[1 1 1]);
set(htb,'SelectionHighlight','off');
set(htb,'FontName','courier');

% Create table head

tpos = [tbpos(1) tbpos(2)+tbpos(4)+10 tbpos(3) 20];
ht=uicontrol(hft,'style','text','position',tpos,'tag','text');
set(ht,'BackgroundColor',[1 1 1]);
set(ht,'HorizontalAlignment','left');
set(ht,'FontName','courier');
newline = '   time         step      order  |    nst   nfe   nni  netf  ncfn';
set(ht,'String',newline);

% Create OK button
  
bsize=[60,28];
badjustpos=[0,25];
bpos=[pos(3)/2-bsize(1)/2+badjustpos(1) -bsize(2)/2+badjustpos(2)...
      bsize(1) bsize(2)];
bpos=round(bpos);
bpos(bpos<1)=1;
hb=uicontrol(hft,'style','pushbutton','position',bpos,...
             'string','Close','tag','okaybutton');
set(hb,'callback','close');

% Save handles

handles=guihandles(hft);
guidata(hft,handles);

for i = 1:n
  newline = '';
  if stats
    newline = sprintf('%10.3e   %10.3e     %1d    |',t(i),h(i),q(i));
  end
  if cntr
    newline = sprintf('%s %5d %5d %5d %5d %5d',...
                      newline,nst(i),nfe(i),nni(i),netf(i),ncfn(i));
  end
  string = get(handles.textbox,'String');
  string{end+1}=newline;
  set(handles.textbox,'String',string);
end

drawnow

%-------------------------------------------------------------------------

function [] = text_update(n,hft,stats,cntr,t,h,q,nst,nfe,nni,netf,ncfn)

figure(hft);

handles=guidata(hft);

for i = 1:n
   if stats
    newline = sprintf('%10.3e   %10.3e     %1d    |',t(i),h(i),q(i));
  end
  if cntr
    newline = sprintf('%s %5d %5d %5d %5d %5d',...
                      newline,nst(i),nfe(i),nni(i),netf(i),ncfn(i));
  end
  string = get(handles.textbox,'String');
  string{end+1}=newline;
  set(handles.textbox,'String',string);
end

drawnow

%-------------------------------------------------------------------------

function [] = sol_init(n, hfs, nps, sol, N, t, y)

fig_name = 'CVODES solution';

% If this is a parallel job, look for the MPI rank in the global
% workspace and append it to the figure name

global sundials_MPI_rank

if ~isempty(sundials_MPI_rank)
  fig_name = sprintf('%s (PE %d)',fig_name,sundials_MPI_rank);
end


figure(hfs);
set(hfs,'Name',fig_name);
set(hfs,'color',[1 1 1]);

% Time label

tlab = '\leftarrow   t   \leftarrow';

% Get number of colors in colormap
map = colormap;
ncols = size(map,1);

% Initialize current subplot counter
pl = 0;

if sol

  pl = pl+1;
  subplot(nps,1,pl);
  hold on;

  for i = 1:N
    hp = plot(t(1:n),y(i,1:n),'-');
    ic = 1+(i-1)*floor(ncols/N);
    set(hp,'Color',map(ic,:));
  end
  box on;
  grid on;
  xlabel(tlab);
  ylabel('y');
  title('Solution');

end

drawnow;

%-------------------------------------------------------------------------

function [] = sol_update(n, hfs, nps, sol, N, t, y)

figure(hfs);

pl = 0;

if sol
  
  pl = pl+1;
  subplot(nps,1,pl);
  
  hc = get(gca,'Children');
  xd = [get(hc(1),'XData') t(1:n)];
% Attention: Children are loaded in reverse order!
  for i = 1:N
    yd = [get(hc(i),'YData') y(N-i+1,1:n)];
    set(hc(i), 'XData', xd, 'YData', yd);
  end

end

drawnow;


%-------------------------------------------------------------------------

function [] = sol_final(hfs, nps, sol, N)

figure(hfs);

pl = 0;

if sol

  pl = pl +1;
  subplot(nps,1,pl);
  
  hc = get(gca,'Children');
  xd = get(hc(1),'XData');
  set(gca,'XLim',sort([xd(1) xd(end)]));

  ylim = get(gca,'YLim');
  addon = 0.1*abs(ylim(2)-ylim(1));
  ylim(1) = ylim(1) + sign(ylim(1))*addon;
  ylim(2) = ylim(2) + sign(ylim(2))*addon;
  set(gca,'YLim',ylim);
  
  for i = 1:N
    cstring{i} = sprintf('y_{%d}',i);
  end
  legend(cstring);
  
end

drawnow
