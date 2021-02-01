function [new_data] = CVodeMonitor(call, T, Y, YQ, YS, data)
%CVodeMonitor is the default CVODES monitoring function.
%   To use it, set the Monitor property in CVodeSetOptions to
%   'CVodeMonitor' or to @CVodeMonitor and 'MonitorData' to mondata
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
%     o sol  [ true | {false} ]
%         If true, plot solution components.
%     o sensi [ true | {false} ]
%         If true and if FSA is enabled, plot sensitivity components.
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
%   See also CVodeSetOptions, CVMonitorFn
%
%   NOTES:
%     1. The argument mondata is REQUIRED. Even if only the default options
%        are desired, set mondata=struct; and pass it to CVodeSetOptions.
%     2. The yQ argument is currently ignored.     

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
% $Revision$Date: 2007/05/11 18:51:32 $

% NOTES:
%   - Unlike Matlab, Octave loads children in the normal order
%   - Unlike Matlab, Octave stores 'XData' and 'YData' as column vectors

if (nargin ~= 6) 
  error('Monitor data not defined.');
end

new_data = [];

if call == 0

% Initialize unspecified fields to default values.
  data = initialize_data(data);  

% Open figure windows
  if data.post

    if data.stats | data.cntr
      data.hfg = figure;
    end
%   Number of subplots in figure hfg
    if data.stats
      data.npg = data.npg + 2;
    end
    if data.cntr
      data.npg = data.npg + 1;
    end
    
    if data.sol | data.sensi
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
% use Y and YS for additional initializations
  
  if data.first

    if isempty(YS)
      data.sensi = false;
    end
    
    if data.sol | data.sensi
      
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
        
      if data.sensi
        data.Ns = size(YS,2);
        data.ys = zeros(data.N, data.Ns, data.updt);
        data.nps = data.nps + data.Ns;
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
  Ns   = data.Ns;
  y    = data.y;
  ys   = data.ys;
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

  si = CVodeGetStats;

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

  if data.sensi
    for k = 1:Ns
      for j = 1:N
        ys(j,k,n) = YS(data.select(j),k);
      end
    end
  end
  
end

% Is it time to post?

if data.post & (n == data.updt | call==2)

  if call == 2
    n = n-1;
  end
  
  if ~data.initialized

    if (data.stats | data.cntr)
      graphical_init(n, hfg, npg, data.stats, data.cntr, ...
                     t, h, q, nst, nfe, nni, netf, ncfn);
    end
    
    if data.sol | data.sensi
      sol_init(n, hfs, nps, data.sol, data.sensi, ...
               N, Ns, t, y, ys);
    end
    
    data.initialized = true;
  
  else

    if (data.stats | data.cntr)
      graphical_update(n, hfg, npg, data.stats, data.cntr, ...
                       t, h, q, nst, nfe, nni, netf, ncfn);
    end

    if data.sol
      sol_update(n, hfs, nps, data.sol, data.sensi, N, Ns, t, y, ys);
    end
      
  end

  if call == 2
    
    if (data.stats | data.cntr)
      graphical_final(hfg, npg, data.cntr, data.stats);
    end
    
    if data.sol | data.sensi
      sol_final(hfs, nps, data.sol, data.sensi, N, Ns);
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
data.ys   = ys;
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
if ~isfield(data,'sensi')
  data.sensi = false;
end
if ~isfield(data,'select')
  data.select = [];
end
if ~isfield(data,'post')
  data.post = true;
end

if ~data.sol & ~data.sensi
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
data.Ns = 0;
data.y = 0;
data.ys = 0;

%-------------------------------------------------------------------------

function [] = graphical_init(n, hfg, npg, stats, cntr, ...
                             t, h, q, nst, nfe, nni, netf, ncfn)

figure(hfg);
pl = 0;

% Time label and figure title

tlab = '->   t   ->';

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
  xd = [get(hc,'XData') ; t(1:n)'];
  yd = [get(hc,'YData') ; abs(h(1:n)')];
  set(hc, 'XData', xd, 'YData', yd);
  
  pl = pl+1;
  subplot(npg,1,pl)
  hc = get(gca,'Children');
  xd = [get(hc,'XData') ; t(1:n)'];
  yd = [get(hc,'YData') ; q(1:n)'];
  set(hc, 'XData', xd, 'YData', yd);
end

% Counters
if cntr
  pl = pl+1;
  subplot(npg,1,pl)
  hc = get(gca,'Children');
  xd = [get(hc(1),'XData') ; t(1:n)'];
  yd = [get(hc(1),'YData') ; ncfn(1:n)'];
  set(hc(1), 'XData', xd, 'YData', yd);
  yd = [get(hc(2),'YData') ; netf(1:n)'];
  set(hc(2), 'XData', xd, 'YData', yd);
  yd = [get(hc(3),'YData') ; nni(1:n)'];
  set(hc(3), 'XData', xd, 'YData', yd);
  yd = [get(hc(4),'YData') ; nfe(1:n)'];
  set(hc(4), 'XData', xd, 'YData', yd);
  yd = [get(hc(5),'YData') ; nst(1:n)'];
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

function [] = sol_init(n, hfs, nps, sol, sensi, N, Ns, t, y, ys)

figure(hfs);

% Time label

tlab = '->   t   ->';

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

if sensi
  
  for is = 1:Ns
    
    pl = pl+1;
    subplot(nps,1,pl);
    hold on;
      
    ys_crt = ys(:,is,1:n);
    for i = 1:N
      hp = plot(t(1:n),ys_crt(i,1:n),'-');
      ic = 1+(i-1)*floor(ncols/N);
      set(hp,'Color',map(ic,:));
    end
    box on;
    grid on;
    xlabel(tlab);
    str = sprintf('s_{%d}',is); ylabel(str);
    str = sprintf('Sensitivity %d',is); title(str);
    
  end
  
end


drawnow;

%-------------------------------------------------------------------------

function [] = sol_update(n, hfs, nps, sol, sensi, N, Ns, t, y, ys)

figure(hfs);

pl = 0;

if sol
  
  pl = pl+1;
  subplot(nps,1,pl);
  
  hc = get(gca,'Children');
  xd = [get(hc(1),'XData') ; t(1:n)'];
  for i = 1:N
    yd = [get(hc(i),'YData') ; y(i,1:n)'];
    set(hc(i), 'XData', xd, 'YData', yd);
  end

end
  
if sensi
  
  for is = 1:Ns
    
    pl = pl+1;
    subplot(nps,1,pl);

    ys_crt = ys(:,is,:);
    
    hc = get(gca,'Children');
    xd = [get(hc(1),'XData') ; t(1:n)'];
    for i = 1:N
      yd = [get(hc(i),'YData') ; ys_crt(i,1:n)'];
      set(hc(i), 'XData', xd, 'YData', yd);
    end
    
  end
  
end


drawnow;


%-------------------------------------------------------------------------

function [] = sol_final(hfs, nps, sol, sensi, N, Ns)

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

if sensi
  
  for is = 1:Ns
    
    pl = pl+1;
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
      cstring{i} = sprintf('s%d_{%d}',is,i);
    end
    legend(cstring);
    
  end
  
end

drawnow
