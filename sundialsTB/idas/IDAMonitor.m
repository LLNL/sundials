function [new_data] = IDAMonitor(call, T, YY, YP, YQ, YYS, YPS, data)
%IDAMonitor is the default IDAS monitoring function.
%   To use it, set the Monitor property in IDASetOptions to
%   'IDAMonitor' or to @IDAMonitor and 'MonitorData' to mondata
%   (defined as a structure).
%  
%   With default settings, this function plots the evolution of the step 
%   size, method order, and various counters.
%   
%   Various properties can be changed from their default values by passing
%   to IDASetOptions, through the property 'MonitorData', a structure
%   MONDATA with any of the following fields. If a field is not defined, 
%   the corresponding default value is used.
%
%   Fields in MONDATA structure:
%     o stats [ {true} | false ]
%         If true, report the evolution of the step size and method order.
%     o cntr [ {true} | false ]
%         If true, report the evolution of the following counters:
%         nst, nre, nni, netf, ncfn (see IDAGetStats)
%     o mode [ {'graphical'} | 'text' | 'both' ] 
%         In graphical mode, plot the evolutions of the above quantities.
%         In text mode, print a table.
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
%     o reinit [ {'continue'} | 'scratch' ]
%         Specifies whether monitoring after a solver reinitialization continues
%         previous monitoring or it starts from scratch.
%     o dir [ {'forward'} | 'backward' ]
%         Specifies forward or backward integration.
%     o post [ {true} | false ]
%         If false, disable all posting. This option is necessary to disable
%         monitoring on some processors when running in parallel.
%
%   See also IDASetOptions, IDAMonitorFn
%
%   NOTES:
%     1. The argument mondata is REQUIRED. Even if only the default options
%        are desired, set mondata=struct; and pass it to IDASetOptions.
%     2. The arguments YP, YQ, and YPS are currently ignored.     

% Radu Serban <radu@llnl.gov>
% Copyright (c) 2005, The Regents of the University of California.
% $Revision: 1.2 $Date: 2006/07/17 16:49:50 $

new_data = [];

% The sequence of call values is:
% 1 0 0 .... 0 0 [ 2 1 0 0 ... 0 0 [ 2 1 0 0 ... 0 0 ] ...] 3

% ------------------------------------------------------------------------
% Initialize data structure with user options and/or default values
% ONLY the first time IDAMonitor is called (the first call with call=1)
% The field 'init_figs' does not exist initially. 
% ------------------------------------------------------------------------

if ~isfield(data,'init_figs')

% Set to their default values all fields in the data structure that were
% not set by the user
  data = initialize_data(data);

% Initialize arrays for statistics and counters
  data.t = zeros(1,data.updt);
  if data.stats
    data.h = zeros(1,data.updt);
    data.q = zeros(1,data.updt);
  end
  if data.cntr
    data.nst = zeros(1,data.updt);
    data.nre = zeros(1,data.updt);
    data.nni = zeros(1,data.updt);
    data.netf = zeros(1,data.updt);
    data.ncfn = zeros(1,data.updt);
  end

% Test if FSA was actually enabled
  if isempty(YYS)
    data.sensi = false;
  end

% Initialize arrays for solution and sensitivities
  if data.sol || data.sensi
    
    if isempty(data.select)
      data.N = length(YY);
      data.select = [1:data.N];
    else
      data.N = length(data.select);
    end
    
    if data.sol
      data.y = zeros(data.N,data.updt);
    end
    
    if data.sensi
      data.Ns = size(YYS,2);
      data.ys = zeros(data.N, data.Ns, data.updt);
    end
    
  end
    
% Set the field init_figs to true so that new figures will be created 
% at the following call=1
  data.init_figs = true;

% No solver reinitialization happened yet.
  data.slv_reinit = false;
  
end
  
% ------------------------------------------------------------------------
% If posting is disabled, return now
% ------------------------------------------------------------------------

if ~data.post
  new_data = data;
  return;
end

% ------------------------------------------------------------------------
% Right before a solver reinitialization, if we will restart monitoring
% from scratch, set all counter offsets to 0; otherwise, get the solver 
% statistics and set offsets for counters.
% ------------------------------------------------------------------------

if call == 2 

  if data.scratch

    data.nst_offset = 0;
    data.nre_offset = 0;
    data.nni_offset = 0;
    data.netf_offset = 0;
    data.ncfn_offset = 0;

  else
    
    data.nst_offset =  data.nst(data.n);
    data.nre_offset = data.nre(data.n);
    data.nni_offset = data.nni(data.n);
    data.netf_offset = data.netf(data.n);
    data.ncfn_offset = data.ncfn(data.n);
    
  end
    
  data.slv_reinit = true;
  
end

% ------------------------------------------------------------------------
% At the initial time (call=1) perform the following operations:
% 1. if needed, initialize the figures
% 2. initialize the skip and buffer counters
% 3. load information from the initial time
% ------------------------------------------------------------------------

if (call == 1)

% If needed, initialize figures
  
  if data.init_figs
  
%   If graphical display is enabled, open graphical figure
%   and count number of subplots
    if data.grph
      data.hfg = figure;
      data.npg = 0;
      if data.stats
        data.npg = data.npg + 2;
      end
      if data.cntr
        data.npg = data.npg + 1;
      end
    end
    
%   If text display is enabled, open text figure
    if data.text
      data.hft = figure;
    end
  
%   If solution and/or sensitivity plotting is enabled, open
%   solution figure and count number of subplots.
    if data.sol || data.sensi
      data.hfs = figure; 
      data.nps = 0;
      if data.sol
        data.nps = data.nps + 1;
      end
      if data.sensi
        data.nps = data.nps + data.Ns;
      end
    end
    
%   Flag that figures have been initialized
    data.init_figs = false;

%   Flag that the next one will be the first data posting
    data.first_post = true;
  
  end

% Data at the initial time is posted if this call is not right after a 
% reinitialization OR if after a reintialization we start from scratch  
  if ~data.slv_reinit || data.scratch
  
%   Initialize buffer counter
    data.n = 1;

%   Load stats at initial time (the negative step size will be treated later)
    data.t(data.n) = T;
    if data.stats
      data.h(data.n) = -1.0;
      data.q(data.n) = 1;
    end
  
%   Load counters at initial time
    if data.cntr
      data.nst(data.n)  = data.nst_offset;
      data.nre(data.n)  = data.nre_offset;
      data.nni(data.n)  = data.nni_offset;
      data.netf(data.n) = data.netf_offset;
      data.ncfn(data.n) = data.ncfn_offset;
    end
    
%   Load solution and sensitivities at initial time
    if data.sol
      for j = 1:data.N
        data.y(j,data.n) = YY(data.select(j));
      end
    end
    if data.sensi
      for k = 1:data.Ns
        for j = 1:data.N
          data.ys(j,k,data.n) = YYS(data.select(j),k);
        end
      end
    end
    
  end
  
% Initialize skip counter  
  data.i = data.skip;

  
end

% ------------------------------------------------------------------------
% At any intermediate call (call=0), test if we should append new info.
% into the buffers (i.e. test whether this step should be skipped or not)
% ------------------------------------------------------------------------

if call == 0

  if data.i ~= 0
    data.i = data.i-1;
    new_data = data;
    return;
  end

  data.n = data.n + 1;
  
  if data.fwd
    si = IDAGetStats;
  else
    si = IDAGetStatsB;
  end

  data.t(data.n) = T;
  
  if data.stats
    data.h(data.n) = si.hlast;
    data.q(data.n) = si.qlast;
  end
  
  if data.cntr
    data.nst(data.n) = si.nst + data.nst_offset;
    data.nre(data.n) = si.nre + data.nre_offset;
    data.nni(data.n) = si.nni + data.nni_offset;
    data.netf(data.n) = si.netf + data.netf_offset;
    data.ncfn(data.n) = si.ncfn + data.ncfn_offset;
  end

  if data.sol
    for j = 1:data.N
      data.y(j,data.n) = YY(data.select(j));
    end
  end

  if data.sensi
    for k = 1:data.Ns
      for j = 1:data.N
        data.ys(j,k,data.n) = YYS(data.select(j),k);
      end
    end
  end
  
  data.i = data.skip;
  
end

% ------------------------------------------------------------------------
% Post information if any of the following conditions is met:
% 1. data.n reached the value data.updt
% 2. the solver will be reinitialized
% 3. this is the final call
% ------------------------------------------------------------------------

if (data.n == data.updt) || (call==2) || (call==3)

%  if call == 3
%    data.n = data.n-1;
%  end

% Post new data
  if data.grph
    graphical_post(data);
  end
  if data.text
    text_post(data);
  end
  if data.sol || data.sensi
    solution_post(data);
  end

  if data.first_post
    data.first_post = false;
  end
    
% Finalize plots if this is the last call OR if the solver will
% be reinitialized and we will start monitoring from scratch
  if (call==3) || (call==2 && data.scratch)
    if data.grph
      graphical_final(data);
    end
    if data.sol || data.sensi
      solution_final(data);
    end
    data.init_figs = true;
  end

% If the solver will be reinitialized AND if we continue monitoring
% add a line in the text figure
  if(call==2 && ~data.scratch) 
    if data.text
      text_final(data);
    end
  end
  
% Reset buffer counter
  data.n = 0;

end

new_data = data;
return;

% ==========================================================================
%
% ==========================================================================

function data = initialize_data(data)

if ~isfield(data,'post')
  data.post = true;
end

data.scratch=false;
if isfield(data,'reinit') && strcmp(data.reinit,'scratch')
  data.scratch=true;
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

if ~isfield(data,'mode')
  data.grph=true;
  data.text=false;
else
  data.grph=true;
  data.text=true;
  if strcmp(data.mode,'graphical')
    data.text=false;
  end
  if strcmp(data.mode,'text')
    data.grph=false;
  end
end

if ~data.stats && ~data.cntr
  data.grph=false;
  data.text=false;
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
if ~data.sol && ~data.sensi
  data.select = [];
end

data.fwd=true;
if isfield(data,'dir') && strcmp(data.dir,'backward')
  data.fwd=false;
end

% Set everything else to zero
data.i = 0;
data.n = 0;
data.npg = 0;
data.nps = 0;
data.hfg = 0;
data.hft = 0;
data.hfs = 0;
data.h = 0;
data.q = 0;
data.nst = 0;
data.nre = 0;
data.nni = 0;
data.netf = 0;
data.ncfn = 0;
data.nst_offset = 0;
data.nre_offset = 0;
data.nni_offset = 0;
data.netf_offset = 0;
data.ncfn_offset = 0;
data.N = 0;
data.Ns = 0;
data.y = 0;
data.ys = 0;

% ==========================================================================
%
% ==========================================================================

function [] = graphical_post(data)

% Focus the appropriate figure 
figure(data.hfg);

% 'Fix' the negative value of stepsize at iintial time
if data.h(1) < 0
  data.h(1) = data.h(2);
end

if data.first_post

% -------------------------------------------
% The first time we post, create the plots
% -------------------------------------------
  
% Set figure name. If this is a parallel job, look for the MPI rank in 
% the global workspace and append it to the figure name
  global sundials_MPI_rank
  fig_name = 'IDAS run statistics';
  if ~isempty(sundials_MPI_rank)
    fig_name = sprintf('%s (PE %d)',fig_name,sundials_MPI_rank);
  end

  set(data.hfg,'Name',fig_name);
  set(data.hfg,'color',[1 1 1]);

% Time label
  if data.fwd
    tlab = '\rightarrow   t   \rightarrow';
  else
    tlab = '\leftarrow   t   \leftarrow';
  end

  pl = 0;

% Step size and order
  if data.stats
    pl = pl+1;
    subplot(data.npg,1,pl)
    semilogy(data.t(1:data.n),abs(data.h(1:data.n)),'-');
    hold on;  box on;  grid on;
    xlabel(tlab);
    ylabel('|Step size|');
    
    pl = pl+1;
    subplot(data.npg,1,pl)
    plot(data.t(1:data.n),data.q(1:data.n),'-');
    hold on;  box on;  grid on;
    xlabel(tlab);
    ylabel('Order');
  end
  
% Counters
  if data.cntr
    pl = pl+1;
    subplot(data.npg,1,pl)
    plot(data.t(1:data.n),data.nst(1:data.n),'k-');
    hold on;
    plot(data.t(1:data.n),data.nre(1:data.n),'b-');
    plot(data.t(1:data.n),data.nni(1:data.n),'r-');
    plot(data.t(1:data.n),data.netf(1:data.n),'g-');
    plot(data.t(1:data.n),data.ncfn(1:data.n),'c-');
    box on;  grid on;
    xlabel(tlab);
    ylabel('Counters');
  end

else

% -------------------------------------------
% At subsequent posts, update plot data
% -------------------------------------------

  pl = 0;
  
% Step size and order
  if data.stats
    pl = pl+1;
    subplot(data.npg,1,pl)
    hc = get(gca,'Children');
    xd = [get(hc,'XData') data.t(1:data.n)];
    yd = [get(hc,'YData') abs(data.h(1:data.n))];
    set(hc, 'XData', xd, 'YData', yd);
    
    pl = pl+1;
    subplot(data.npg,1,pl)
    hc = get(gca,'Children');
    xd = [get(hc,'XData') data.t(1:data.n)];
    yd = [get(hc,'YData') data.q(1:data.n)];
    set(hc, 'XData', xd, 'YData', yd);
  end
  
% Counters
  if data.cntr
    pl = pl+1;
    subplot(data.npg,1,pl)
    hc = get(gca,'Children');
%   Attention: Children are loaded in reverse order!
    xd = [get(hc(1),'XData') data.t(1:data.n)];
    yd = [get(hc(1),'YData') data.ncfn(1:data.n)];
    set(hc(1), 'XData', xd, 'YData', yd);
    yd = [get(hc(2),'YData') data.netf(1:data.n)];
    set(hc(2), 'XData', xd, 'YData', yd);
    yd = [get(hc(3),'YData') data.nni(1:data.n)];
    set(hc(3), 'XData', xd, 'YData', yd);
    yd = [get(hc(4),'YData') data.nre(1:data.n)];
    set(hc(4), 'XData', xd, 'YData', yd);
    yd = [get(hc(5),'YData') data.nst(1:data.n)];
    set(hc(5), 'XData', xd, 'YData', yd);
  end

end
  
drawnow;

% ==========================================================================
%
% ==========================================================================

function [] = graphical_final(data)

% Focus the appropriate figure 
figure(data.hfg);

pl = 0;

if data.stats
  pl = pl+1;
  subplot(data.npg,1,pl)
  hc = get(gca,'Children');
  xd = get(hc,'XData');
  set(gca,'XLim',sort([xd(1) xd(end)]));
  
  pl = pl+1;
  subplot(data.npg,1,pl)
  ylim = get(gca,'YLim');
  ylim(1) = ylim(1) - 1;
  ylim(2) = ylim(2) + 1;
  set(gca,'YLim',ylim);
  set(gca,'XLim',sort([xd(1) xd(end)]));
end

if data.cntr
  pl = pl+1;
  subplot(data.npg,1,pl)
  hc = get(gca,'Children');
  xd = get(hc(1),'XData');
  set(gca,'XLim',sort([xd(1) xd(end)]));
  legend('nst','nre','nni','netf','ncfn',2);
end

% ==========================================================================
%
% ==========================================================================

function [] = text_post(data)

% Focus the appropriate figure 
figure(data.hft);

if data.first_post

% -------------------------------------------
% The first time we post, create table
% -------------------------------------------

% Set figure name. If this is a parallel job, look for the MPI rank in 
% the global workspace and append it to the figure name
  global sundials_MPI_rank
  fig_name = 'IDAS run statistics';
  if ~isempty(sundials_MPI_rank)
    fig_name = sprintf('%s (PE %d)',fig_name,sundials_MPI_rank);
  end

  set(data.hft,'Name',fig_name);
  set(data.hft,'color',[1 1 1]);
  set(data.hft,'MenuBar','none');
  set(data.hft,'Resize','off');
  
% Create text box
  margins=[10 10 50 50]; % left, right, top, bottom
  pos=get(data.hft,'position');
  tbpos=[margins(1) margins(4) pos(3)-margins(1)-margins(2) ...
         pos(4)-margins(3)-margins(4)];
  tbpos(tbpos<1)=1;
  
  htb=uicontrol(data.hft,'style','listbox','position',tbpos,'tag','textbox');
  set(htb,'BackgroundColor',[1 1 1]);
  set(htb,'SelectionHighlight','off');
  set(htb,'FontName','courier');
  
% Create table head
  tpos = [tbpos(1) tbpos(2)+tbpos(4)+10 tbpos(3) 20];
  ht=uicontrol(data.hft,'style','text','position',tpos,'tag','text');
  set(ht,'BackgroundColor',[1 1 1]);
  set(ht,'HorizontalAlignment','left');
  set(ht,'FontName','courier');
  newline = '   time         step      order  |    nst   nre   nni  netf  ncfn';
  set(ht,'String',newline);
  
% Create 'Close' button
  bsize=[60,28];
  badjustpos=[0,25];
  bpos=[pos(3)/2-bsize(1)/2+badjustpos(1) -bsize(2)/2+badjustpos(2)...
        bsize(1) bsize(2)];
  bpos=round(bpos);
  bpos(bpos<1)=1;
  hb=uicontrol(data.hft,'style','pushbutton','position',bpos,...
               'string','Close','tag','okaybutton');
  set(hb,'callback','close');
  
% Save handles
  handles=guihandles(data.hft);
  guidata(data.hft,handles);

  for i = 1:data.n
    newline = '';
    if data.stats
      if(data.h(i) < 0)
        newline = sprintf('%10.3e                       |',data.t(i));
      else
        newline = sprintf('%10.3e   %10.3e     %1d    |',data.t(i),data.h(i),data.q(i));
      end
    end
    if data.cntr
      newline = sprintf('%s %5d %5d %5d %5d %5d',...
                        newline,data.nst(i),data.nre(i),data.nni(i),data.netf(i),data.ncfn(i));
    end
    string = get(handles.textbox,'String');
    string{end+1}=newline;
    set(handles.textbox,'String',string);
  end

else
  
% -------------------------------------------
% At subsequent posts, update table data
% -------------------------------------------
  
  handles=guidata(data.hft);
  
  for i = 1:data.n
    newline = '';
    if data.stats
      if(data.h(i) < 0)
        newline = sprintf('%10.3e                       |',data.t(i));
      else
        newline = sprintf('%10.3e   %10.3e     %1d    |',data.t(i),data.h(i),data.q(i));
      end
    end
    if data.cntr
      newline = sprintf('%s %5d %5d %5d %5d %5d',...
                        newline,data.nst(i),data.nre(i),data.nni(i),data.netf(i),data.ncfn(i));
    end
    string = get(handles.textbox,'String');
    string{end+1}=newline;
    set(handles.textbox,'String',string);
  end

end
  
drawnow

% ==========================================================================
%
% ==========================================================================

function [] = text_final(data)

% Focus the appropriate figure 
figure(data.hft);

handles=guidata(data.hft);
newline = '---------------------------------+-------------------------------';
string = get(handles.textbox,'String');
string{end+1}=newline;
set(handles.textbox,'String',string);


% ==========================================================================
%
% ==========================================================================

function [] = solution_post(data)

% Focus the appropriate figure 
figure(data.hfs);

if data.first_post

% -------------------------------------------
% The first time we post, create the plots
% -------------------------------------------

% Set figure name. If this is a parallel job, look for the MPI rank in 
% the global workspace and append it to the figure name
  global sundials_MPI_rank
  fig_name = 'IDAS solution';
  if ~isempty(sundials_MPI_rank)
    fig_name = sprintf('%s (PE %d)',fig_name,sundials_MPI_rank);
  end

  set(data.hfs,'Name',fig_name);
  set(data.hfs,'color',[1 1 1]);

% Time label
  if data.fwd
    tlab = '\rightarrow   t   \rightarrow';
  else
    tlab = '\leftarrow   t   \leftarrow';
  end

% Get number of colors in colormap
  map = colormap;
  ncols = size(map,1);

% Initialize current subplot counter
  pl = 0;

  if data.sol

    pl = pl+1;
    subplot(data.nps,1,pl);
    hold on;

    for i = 1:data.N
      hp = plot(data.t(1:data.n),data.y(i,1:data.n),'-');
      ic = 1+(i-1)*floor(ncols/data.N);
      set(hp,'Color',map(ic,:));
    end
    
    box on;  grid on;
    xlabel(tlab);
    ylabel('y');
    title('Solution');
    
  end

  if data.sensi
    
    for is = 1:data.Ns
      
      pl = pl+1;
      subplot(data.nps,1,pl);
      hold on;
      
      ys_crt = data.ys(:,is,1:data.n);
      for i = 1:data.N
        hp = plot(data.t(1:data.n),ys_crt(i,1:data.n),'-');
        ic = 1+(i-1)*floor(ncols/data.N);
        set(hp,'Color',map(ic,:));
      end
      box on;  grid on;
      xlabel(tlab);
      str = sprintf('s_{%d}',is); ylabel(str);
      str = sprintf('Sensitivity %d',is); title(str);
      
    end
    
  end

else

% -------------------------------------------
% At subsequent posts, update plot data
% -------------------------------------------

  pl = 0;

  if data.sol
  
    pl = pl+1;
    subplot(data.nps,1,pl);
    
    hc = get(gca,'Children');
    xd = [get(hc(1),'XData') data.t(1:data.n)];
%   Attention: Children are loaded in reverse order!
    for i = 1:data.N
      yd = [get(hc(i),'YData') data.y(data.N-i+1,1:data.n)];
      set(hc(i), 'XData', xd, 'YData', yd);
    end
    
  end
  
  if data.sensi
  
    for is = 1:data.Ns
    
      pl = pl+1;
      subplot(data.nps,1,pl);

      ys_crt = data.ys(:,is,:);
    
      hc = get(gca,'Children');
      xd = [get(hc(1),'XData') data.t(1:data.n)];
%     Attention: Children are loaded in reverse order!
      for i = 1:data.N
        yd = [get(hc(i),'YData') ys_crt(data.N-i+1,1:data.n)];
        set(hc(i), 'XData', xd, 'YData', yd);
      end
      
    end
    
  end

end

drawnow;

% ==========================================================================
%
% ==========================================================================

function [] = solution_final(data)

figure(data.hfs);

pl = 0;

if data.sol

  pl = pl +1;
  subplot(data.nps,1,pl);
  
  hc = get(gca,'Children');
  xd = get(hc(1),'XData');
  set(gca,'XLim',sort([xd(1) xd(end)]));

  ylim = get(gca,'YLim');
  addon = 0.1*abs(ylim(2)-ylim(1));
  ylim(1) = ylim(1) + sign(ylim(1))*addon;
  ylim(2) = ylim(2) + sign(ylim(2))*addon;
  set(gca,'YLim',ylim);
  
  for i = 1:data.N
    cstring{i} = sprintf('y_{%d}',i);
  end
  legend(cstring);
  
end

if data.sensi
  
  for is = 1:data.Ns
    
    pl = pl+1;
    subplot(data.nps,1,pl);

    hc = get(gca,'Children');
    xd = get(hc(1),'XData');
    set(gca,'XLim',sort([xd(1) xd(end)]));

    ylim = get(gca,'YLim');
    addon = 0.1*abs(ylim(2)-ylim(1));
    ylim(1) = ylim(1) + sign(ylim(1))*addon;
    ylim(2) = ylim(2) + sign(ylim(2))*addon;
    set(gca,'YLim',ylim);
  
    for i = 1:data.N
      cstring{i} = sprintf('s%d_{%d}',is,i);
    end
    legend(cstring);
    
  end
  
end

drawnow
