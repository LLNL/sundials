function varargout = kim_info(action, fin, message);

switch action
  
 case 0  % initialize 

  % Create figure
  
  f = figure;

  set(f,'resizefcn','kim_info(2,0,0)');
  set(f,'name','KINSOL info','numbertitle','off');
  set(f,'menubar','none','tag','figure');
  
  % Create text box
   
  tbpos=getTBPos(f);
  h=uicontrol(f,'style','listbox','position',tbpos,'tag','textbox');
  set(h,'BackgroundColor',[1 1 1]);
  set(h,'SelectionHighlight','off');

  % Create OK button
  
  bpos=getOKPos(f);
  h=uicontrol(f,'style','pushbutton','position',bpos,'string','Close','tag','okaybutton');
  set(h,'callback','kim_info(3,0,0)');

  % Save handles
  
  handles=guihandles(f);
  guidata(f,handles);
 
  varargout{1} = f;
 
 case 1 % append text

  f = fin;
  new_str = message;

  handles=guidata(f);
  string = get(handles.textbox,'String');
  string{end+1}=new_str;
  set(handles.textbox,'String',string);
  
 case 2 % resize

  handles=guidata(gcbo);
  tbpos=getTBPos(handles.figure);
  bpos=getOKPos(handles.figure);
  set(handles.okaybutton,'position',bpos);
  set(handles.textbox,'position',tbpos);

 case 3 % close
  
  handles=guidata(gcbo);
  close(handles.figure);
  
end

%------------------------------------
function tbpos=getTBPos(f)

margins=[10 10 10 50]; % left, right, top, bottom
pos=get(f,'position');
tbpos=[margins(1) margins(4) pos(3)-margins(1)-margins(2) ...
    pos(4)-margins(3)-margins(4)];
tbpos(tbpos<1)=1;


%------------------------------------
function tbpos=getOKPos(f)

bsize=[60,30];
badjustpos=[0,25];

pos=get(f,'position');

tbpos=[pos(3)/2-bsize(1)/2+badjustpos(1) -bsize(2)/2+badjustpos(2)...
    bsize(1) bsize(2)];
tbpos=round(tbpos);
tbpos(tbpos<1)=1;
