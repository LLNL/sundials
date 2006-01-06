function [] = CVodeMonitor(call, time, sol, yQ, yS, varargin)
%CVodeMonitor is a simple monitoring function example.
%   To use it, set the Monitor property in CVodeSetOptions to
%   'CVodeMonitor' or to @CVodeMonitor.
%  
%   With default settings, this function plots the evolution of the step 
%   size, method order, and various counters.
%   
%   Various properties can be changed from their default values by passing
%   to CvodeSetOptions, through the property 'MonitorData', a structure
%   MONDATA with any of the following fields. If a field is not defined, 
%   the corresponding default value is used.
%
%   Fields in MONDATA structure:
%     o stats [ {true} | false ]
%         If true, CVodeMonitor reports the evolution of the step size and
%         method order.
%     o cntr [ {true} | false ]
%         If true, CVodeMonitor reports the evolution of the following counters:
%         nst, nfe, nni, netf, ncfn (see CVodeGetStats)
%     o sol  [ true | {false} ]
%         If true, CvodeMonitor plots all solution components (graphical mode only).
%     o grph [ {true} | false ] 
%         If true, CvodeMonitor plots the evolutions of the above quantities.
%         Otherwise, it prints to the screen.
%     o updt [ integer | {50} ]
%         CvodeMonitor update frequency.
%     o select [ array of integers ]
%         To plot only particular solution components, specify their indeces in
%         the field select. If not defined, but sol=true, CVodeMonitor plots all 
%         components (graphical mode only).
%     o xaxis [ linear | {log} ]
%         Type of the time axis for the stepsize, order, and counter plots 
%         (graphical mode only).
%     o dir [ {1} | -1 ]
%         Specifies forward or backward integration.
%
%   See also CVodeSetOptions, CVMonitorFn
%
%   NOTE:
%     CVodeMonitor ignores the yQ and yS arguments.

% Radu Serban <radu@llnl.gov>
% Copyright (c) 2005, The Regents of the University of California.
% $Revision: 1.1 $Date$

persistent data
persistent first
persistent hf1 hf2 npl
persistent i
persistent t y h q nst nfe nni netf ncfn


% If this is the first call, initialize static variables and return

if call == 0

  if nargin > 5
    data = varargin{1};
  end

  data = initialize_data(data, length(sol));
  
  first = true;
  if data.grph
    npl = 0;
    if data.stats
      npl = npl + 2;
    end
    if data.cntr
      npl = npl + 1;
    end
    if npl ~= 0
      hf1 = figure;
    end
  end
  if data.sol
    hf2 = figure; 
    colormap(data.map);
  end

  i = 1;
  t = zeros(1,data.updt);
  if data.stats
    h = zeros(1,data.updt);
    q = zeros(1,data.updt);
  end
  if data.cntr
    nst = zeros(1,data.updt);
    nfe = zeros(1,data.updt);
    nni = zeros(1,data.updt);
    netf = zeros(1,data.updt);
    ncfn = zeros(1,data.updt);
  end
  if data.sol
    N = length(data.select);
    y = zeros(N,data.updt);
  end

  return;
  
end

% Load current statistics

if data.dir == 1
  si = CVodeGetStats;
else
  si = CVodeGetStatsB;
end

t(i) = si.tcur;

if data.stats
  h(i) = si.hlast;
  q(i) = si.qlast;
end
  
if data.cntr
  nst(i) = si.nst;
  nfe(i) = si.nfe;
  nni(i) = si.nni;
  netf(i) = si.netf;
  ncfn(i) = si.ncfn;
end

if data.sol
  N = length(data.select);
  for j = 1:N
    y(j,i) = sol(data.select(j));
  end
else
  N = 0;
end

% Finalize post

if call == 2
  if data.grph
    graphical_final(i,...
                    hf1, npl, data.stats, data.cntr, data.sol, data.dir,...
                    t, h, q, nst, nfe, nni, netf, ncfn,...
                    hf2, y, N, data.select);
  else
    text_final(i,data.stats,data.cntr,t,h,q,nst,nfe,nni,netf,ncfn);
  end
  return
end

% Is it time to post?

if i == data.updt

  if first
    if data.grph
      graphical_init(hf1, npl, data.stats, data.cntr, data.sol, data.dir,...
                     t, h, q, nst, nfe, nni, netf, ncfn,...
                     hf2, y, N, data.xaxis);
    else
      text_update(data.stats,data.cntr,t,h,q,nst,nfe,nni,netf,ncfn);
    end
    first = false;
  else
    if data.grph
      graphical_update(hf1, npl, data.stats, data.cntr, data.sol, data.dir,...
                       t, h, q, nst, nfe, nni, netf, ncfn,...
                       hf2, y, N);
    else
      text_update(data.stats,data.cntr,t,h,q,nst,nfe,nni,netf,ncfn);
    end
  end
  i = 1;
else
  i = i + 1;
end

% If this was the last call, reset data

if call == 2
  data = [];
end


%-------------------------------------------------------------------------

function data = initialize_data(data, N)

if ~isfield(data,'grph')
  data.grph = true;
end
if ~isfield(data,'updt')
  data.updt = 50;
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
if ~isfield(data,'map')
  data.map = 'default';
end
if ~isfield(data,'select')
  data.select = [1:N];
end
if ~isfield(data,'xaxis')
  data.xaxis = 'log';
end
if ~isfield(data,'dir')
  data.dir = 1;
end

if ~data.grph
  data.sol = false;
end

%-------------------------------------------------------------------------

function [] = graphical_init(hf1, npl, stats, cntr, sol, dir,...
                             t, h, q, nst, nfe, nni, netf, ncfn,...
                             hf2, y, N, xaxis)

if npl ~= 0
  figure(hf1);
  pl = 0;
end

% Step size and order
if stats
  pl = pl+1;
  subplot(npl,1,pl)
  semilogy(t,abs(h),'-');
  if strcmp(xaxis,'log')
    set(gca,'XScale','log');
  end
  hold on;
  box on;
  grid on;
  xlabel('t');
  ylabel('|Step size|');
  
  pl = pl+1;
  subplot(npl,1,pl)
  plot(t,q,'-');
  if strcmp(xaxis,'log')
    set(gca,'XScale','log');
  end
  hold on;
  box on;
  grid on;
  xlabel('t');
  ylabel('Order');
end
  
% Counters
if cntr
  pl = pl+1;
  subplot(npl,1,pl)
  semilogy(t,nst,'k-');
  hold on;
  semilogy(t,nfe,'b-');
  semilogy(t,nni,'r-');
  semilogy(t,netf,'g-');
  semilogy(t,ncfn,'c-');
  if strcmp(xaxis,'log')
    set(gca,'XScale','log');
  end
  box on;
  grid on;
  xlabel('t');
  ylabel('Counters');
end

% Solution components
if sol
  figure(hf2);
  map = colormap;
  ncols = size(map,1);
  hold on;
  for i = 1:N
    hp = plot(t,y(i,:),'-');
    ic = 1+(i-1)*floor(ncols/N);
    set(hp,'Color',map(ic,:));
  end
  if strcmp(xaxis,'log')
    set(gca,'XScale','log');
  end
  box on;
  grid on;
  xlabel('t');
  ylabel('y');
  title('Solution');
end

drawnow;

%-------------------------------------------------------------------------

function [] = graphical_update(hf1, npl, stats, cntr, sol, dir,...
                               t, h, q, nst, nfe, nni, netf, ncfn,...
                               hf2, y, N)

if npl ~= 0
  figure(hf1);
  pl = 0;
end
  
% Step size and order
if stats
  pl = pl+1;
  subplot(npl,1,pl)
  hc = get(gca,'Children');
  xd = [get(hc,'XData') t];
  yd = [get(hc,'YData') abs(h)];
  if length(xd) ~= length(yd)
    disp('h');
  end
  set(hc, 'XData', xd, 'YData', yd);

  pl = pl+1;
  subplot(npl,1,pl)
  hc = get(gca,'Children');
  xd = [get(hc,'XData') t];
  yd = [get(hc,'YData') q];
  if length(xd) ~= length(yd)
    disp('q');
  end
  set(hc, 'XData', xd, 'YData', yd);
end
  
% Counters
if cntr
  pl = pl+1;
  subplot(npl,1,pl)
  hc = get(gca,'Children');
  % Attention: Children are loaded in reverse order!
  xd = [get(hc(1),'XData') t];
  yd = [get(hc(1),'YData') ncfn];
  set(hc(1), 'XData', xd, 'YData', yd);
  yd = [get(hc(2),'YData') netf];
  set(hc(2), 'XData', xd, 'YData', yd);
  yd = [get(hc(3),'YData') nni];
  set(hc(3), 'XData', xd, 'YData', yd);
  yd = [get(hc(4),'YData') nfe];
  set(hc(4), 'XData', xd, 'YData', yd);
  yd = [get(hc(5),'YData') nst];
  set(hc(5), 'XData', xd, 'YData', yd);
end

% Solution components
if sol
  figure(hf2);
  hc = get(gca,'Children');
  xd = [get(hc(1),'XData') t];
  % Attention: Children are loaded in reverse order!
  for i = 1:N
    yd = [get(hc(i),'YData') y(N-i+1,:)];
    set(hc(i), 'XData', xd, 'YData', yd);
  end
end

drawnow;

%-------------------------------------------------------------------------

function [] = graphical_final(n, hf1, npl, stats, cntr, sol, dir,...
                              t, h, q, nst, nfe, nni, netf, ncfn,...
                              hf2, y, N, select)

if npl ~= 0
  figure(hf1);
  pl = 0;
end
  
% Step size and order
if stats
  pl = pl+1;
  subplot(npl,1,pl)
  hc = get(gca,'Children');
  xd = [get(hc,'XData') t(1:n-1)];
  yd = [get(hc,'YData') abs(h(1:n-1))];
  set(hc, 'XData', xd, 'YData', yd);
%  xlim = get(gca,'XLim');
%  set(gca,'XLim',[xlim(1) t(n-1)]);
  
  pl = pl+1;
  subplot(npl,1,pl)
  hc = get(gca,'Children');
  xd = [get(hc,'XData') t(1:n-1)];
  yd = [get(hc,'YData') q(1:n-1)];
  set(hc, 'XData', xd, 'YData', yd);
%  xlim = get(gca,'XLim');
%  set(gca,'XLim',[xlim(1) t(n-1)]);
  ylim = get(gca,'YLim');
  set(gca,'YLim',[ylim(1)-1 ylim(2)+1]);
end

% Counters
if cntr
  pl = pl+1;
  subplot(npl,1,pl)
  hc = get(gca,'Children');
  xd = [get(hc(1),'XData') t(1:n-1)];
  yd = [get(hc(1),'YData') ncfn(1:n-1)];
  set(hc(1), 'XData', xd, 'YData', yd);
  yd = [get(hc(2),'YData') netf(1:n-1)];
  set(hc(2), 'XData', xd, 'YData', yd);
  yd = [get(hc(3),'YData') nni(1:n-1)];
  set(hc(3), 'XData', xd, 'YData', yd);
  yd = [get(hc(4),'YData') nfe(1:n-1)];
  set(hc(4), 'XData', xd, 'YData', yd);
  yd = [get(hc(5),'YData') nst(1:n-1)];
  set(hc(5), 'XData', xd, 'YData', yd);
%  xlim = get(gca,'XLim');
%  set(gca,'XLim',[xlim(1) t(n-1)]);
  legend('nst','nfe','nni','netf','ncfn',2);
end

% Solution components
if sol
  figure(hf2);
  hc = get(gca,'Children');
  xd = [get(hc(1),'XData') t(1:n-1)];
  % Attention: Children are loaded in reverse order!
  for i = 1:N
    yd = [get(hc(i),'YData') y(N-i+1,1:n-1)];
    set(hc(i), 'XData', xd, 'YData', yd);
    cstring{i} = sprintf('y_{%d}',i);
  end
  legend(cstring);
end

drawnow;

%-------------------------------------------------------------------------

function [] = text_init(stats,cntr,t,h,q,nst,nfe,nni,netf,ncfn)

%-------------------------------------------------------------------------

function [] = text_update(stats,cntr,t,h,q,nst,nfe,nni,netf,ncfn)

n = length(t);
for i = 1:n
  if stats
    fprintf('%8.3e %12.6e %1d | ',t(i),h(i),q(i));
  end
  if cntr
    fprintf('%5d %5d %5d %5d %5d\n',nst(i),nfe(i),nni(i),netf(i),ncfn(i));
  else
    fprintf('\n');
  end
end
fprintf('-----\n');
  
%-------------------------------------------------------------------------

function [] = text_final(n,stats,cntr,t,h,q,nst,nfe,nni,netf,ncfn)

for i = 1:n-1
  if stats
    fprintf('%8.3e %12.6e %1d | ',t(i),h(i),q(i));
  end
  if cntr
    fprintf('%5d %5d %5d %5d %5d\n',nst(i),nfe(i),nni(i),netf(i),ncfn(i));
  else
    fprintf('\n');
  end
end
fprintf('-----\n');
