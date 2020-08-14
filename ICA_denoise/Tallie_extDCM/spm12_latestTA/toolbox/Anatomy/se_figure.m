
function varargout=se_figure(varargin)

global st;
global MAP;
global SPM;
global VOL;
global CLUSTER;
global PMap;
global displayType;


%-Condition arguments
%-----------------------------------------------------------------------
if (nargin==0), Action = 'Create'; else, Action = varargin{1}; end

switch lower(Action), case 'create'
%=======================================================================
% F = se_figure('Create',Tag,Name,Visible)
%-Condition arguments
if nargin<4, Visible='on'; else, Visible=varargin{4}; end
if nargin<3, Name=''; else, Name=varargin{3}; end
if nargin<2, Tag=''; else, Tag=varargin{2}; end

F = se_figure('CreateWin',Tag,Name,Visible);
se_figure('CreateBar',F)
varargout = {F};


case 'findwin'
%=======================================================================
% F=se_figure('FindWin',F)
% F=se_figure('FindWin',Tag)
%-Find window: Find window with FigureNumber# / 'Tag' attribute
%-Returns empty if window cannot be found - deletes multiple tagged figs.

if nargin<2, F='Graphics'; else, F=varargin{2}; end

if isempty(F)
	% Leave F empty
elseif ischar(F)
	% Finds Graphics window with 'Tag' string - delete multiples
	Tag=F;
	F = findobj(get(0,'Children'),'Flat','Tag',Tag);
	if length(F) > 1
		% Multiple Graphics windows - close all but most recent
		close(F(2:end))
		F = F(1);
	end
else
	% F is supposed to be a figure number - check it
	if ~any(F==get(0,'Children')), F=[]; end
end
varargout = {F};

case 'getwin'
%=======================================================================
% F=se_figure('GetWin',Tag)

if nargin<2, Tag='Graphics'; else, Tag=varargin{2}; end
F = se_figure('FindWin',Tag);

if isempty(F)
	if ischar(Tag)
		switch Tag, case 'Graphics'
			F = se_figure('Create','Graphics','Graphics');
		case 'Interactive'
			F = spm('CreateIntWin');
		end
	end
else
	set(0,'CurrentFigure',F);
end
varargout = {F};

case 'parentfig'
%=======================================================================
% F=se_figure('ParentFig',h)
if nargin<2, error('No object specified'), else, h=varargin{2}; end
F = get(h(1),'Parent');
while ~strcmp(get(F,'Type'),'figure'), F=get(F,'Parent'); end
varargout = {F};


case 'clear'
%=======================================================================
% se_figure('Clear',F,Tags)

%-Sort out arguments
%-----------------------------------------------------------------------
if nargin<3, Tags=[]; else, Tags=varargin{3}; end
if nargin<2, F=get(0,'CurrentFigure'); else, F=varargin{2}; end
F = se_figure('FindWin',F);
if isempty(F), return, end

%-Clear figure
%-----------------------------------------------------------------------
if isempty(Tags)
	%-Clear figure of objects with 'HandleVisibility' 'on'
	delete(findobj(get(F,'Children'),'flat','HandleVisibility','on'));
	%-Reset figures callback functions
	set(F,'KeyPressFcn','',...
		'WindowButtonDownFcn','',...
		'WindowButtonMotionFcn','',...
		'WindowButtonUpFcn','')
	%-If this is the 'Interactive' window, reset name & UserData
	if strcmp(get(F,'Tag'),'Interactive')
		set(F,'Name','','UserData',[]), end
else
	%-Clear specified objects from figure
	cSHH = get(0,'ShowHiddenHandles');
	set(0,'ShowHiddenHandles','on')
	if ischar(Tags); Tags=cellstr(Tags); end
	if any(strcmp(Tags(:),'!all'))
		delete(get(F,'Children'))
	else
	    for tag = Tags(:)'
		delete(findobj(get(F,'Children'),'flat','Tag',tag{:}));
	    end
	end	
	set(0,'ShowHiddenHandles',cSHH)
end
set(F,'Pointer','Arrow')


case 'defprintcmd'
%=======================================================================
% se_figure('DefPrintCmd')
varargout = {'print -dpsc2 -painters -append -noui'};


case 'print'
%=======================================================================
%fg = get(0,'CurrentFigure');
fg = se_figure('Findwin','Graphics');
figure(fg);

if strcmp(displayType,'OL')

    FS        = spm('FontSizes');			%-Scaled font sizes
    PF        = spm_platform('fonts');		%-Font names (for this platform)
    hA2   = axes('Position',[0.001 0.05 0.45 0.15],...
    	'DefaultTextInterpreter','Tex',...
    	'DefaultTextVerticalAlignment','Baseline',...
    	'Units','points',...
    	'Visible','off');
        AxPos = get(hA2,'Position'); set(hA2,'YLim',[0,AxPos(4)]);
    
    y     = floor(AxPos(4));
    xyzmm = se_orthviews('pos');

    text(0,y-10,[sprintf('Crosshair position: x=%1.0f, y=%1.0f, z=%1.0f,',xyzmm)],...
         'FontSize',FS(9))

     
     hAx   = axes('Position',[0.48 0.05 0.51 0.56],...
    	'DefaultTextInterpreter','Tex',...
    	'DefaultTextVerticalAlignment','Baseline',...
    	'Units','points',...
    	'Visible','off');
        AxPos = get(hAx,'Position'); set(hAx,'YLim',[0,AxPos(4)]);
    
    y     = floor(AxPos(4));

    text(0,y-10,['Maximum probability map: ' spm_str_manip(MAP(1).MaxMap.fname,'rt')],...
         'FontSize',FS(10))
    text(0,y-30,get(st.titel,'String'),'FontWeight','Bold','FontSize',FS(9));
    text(0,y-75,get(st.vAssign,'String'),'FontSize',FS(7));
    text(0,y-125,get(st.mAssign,'String'),'FontSize',FS(7));
    
elseif strcmp(displayType,'PX')
    FS        = spm('FontSizes');			%-Scaled font sizes
    PF        = spm_platform('fonts');		%-Font names (for this platform)
    hA2   = axes('Position',[0.05 0.05 0.45 0.15],...
    	'DefaultTextInterpreter','Tex',...
    	'DefaultTextVerticalAlignment','Baseline',...
    	'Units','points',...
    	'Visible','off');
        AxPos = get(hA2,'Position'); set(hA2,'YLim',[0,AxPos(4)]);
    
    y     = floor(AxPos(4));
    xyzmm = se_orthviews('pos');

    haha = text(0,y-10,[sprintf('Crosshair position: x=%1.0f, y=%1.0f, z=%1.0f,',xyzmm)],...
         'FontSize',FS(12));

     
     hAx   = axes('Position',[0.51 0.05 0.45 0.56],...
    	'DefaultTextInterpreter','Tex',...
    	'DefaultTextVerticalAlignment','Baseline',...
    	'Units','points',...
    	'Visible','off');
        AxPos = get(hAx,'Position'); set(hAx,'YLim',[0,AxPos(4)]);
    
    y     = floor(AxPos(4));

    text(0,y-10,['Maximum probability map: ' spm_str_manip(MAP(1).MaxMap.fname,'rt')],...
         'FontSize',FS(10));
    haha2 = text(0,y-40,get(st.mAssign,'String'),'FontSize',FS(10));
 
elseif strcmp(displayType,'GA')
    if isfield(st,'hAx')
    else
    FS        = spm('FontSizes');			%-Scaled font sizes
    PF        = spm_platform('fonts');		%-Font names (for this platform)
    
    
     hA2   = axes('Position',[0.05 0.05 0.45 0.15],...
    	'DefaultTextInterpreter','Tex',...
    	'DefaultTextVerticalAlignment','Baseline',...
    	'Units','points',...
    	'Visible','off');
        AxPos = get(hA2,'Position'); set(hA2,'YLim',[0,AxPos(4)]);
    
    y     = floor(AxPos(4));
    xyzmm = se_orthviews('pos');

    text(0,y-10,[sprintf('Crosshair position: x=%1.0f, y=%1.0f, z=%1.0f,',xyzmm)],...
         'FontSize',FS(13))
   
    
    hAx   = axes('Position',[0.51 0.05 0.45 0.56],...
    	'DefaultTextInterpreter','Tex',...
    	'DefaultTextVerticalAlignment','Baseline',...
    	'Units','points',...
    	'Visible','off');
        AxPos = get(hAx,'Position'); set(hAx,'YLim',[0,AxPos(4)]);
    
    y     = floor(AxPos(4));
    
    text(0,y-10,['Maximum probability map: ' spm_str_manip(MAP(1).MaxMap.fname,'rt')],...
         'FontSize',FS(11))
    text(0,y-30,get(st.titel,'String'),'FontWeight','Bold','FontSize',FS(13));
    text(0,y-55,get(st.vAssign,'String'),'FontSize',FS(9));
    if get(st.maxima,'Value') > 1
        text(0,y-80,['Maximum ' sprintf('%2.0f', get(st.maxima,'Value')-1) ' selected'],'FontSize',FS(9));
    end
    text(0,y-100,get(st.mAssign,'String'),'FontSize',FS(9));
        
        
        
        
    end
    
elseif strcmp(displayType,'PM')
    FS        = spm('FontSizes');			%-Scaled font sizes
    PF        = spm_platform('fonts');		%-Font names (for this platform)
    pos = round(se_orthviews('pos'));

    pAx   = axes('Position',[0.5 0.58 0.5 0.06],...
    	'DefaultTextInterpreter','Tex',...
    	'DefaultTextVerticalAlignment','Baseline',...   
    	'DefaultTextHorizontalAlignment','center',...
        	'Units','points',...
    	'Visible','off');
        AxPos = get(pAx,'Position'); 
        set(pAx,'YLim',[0,AxPos(4)]); set(pAx,'XLim',[0,AxPos(3)]);

    text(20,20,['Crosshair position:    x = ' int2str(pos(1))...
             '    y = ' int2str(pos(2)) '    z = ' int2str(pos(3))],'FontSize',FS(8))

end    

FNote = sprintf('%s%s: %s',spm('ver'),spm('GetUser',' (%s)'),spm('time'));
%-Delete old tag lines, and print new one
delete(findobj(fg,'Tag','SPMprintFootnote'));
fAx = axes('Position',[0.005,0.005,0.1,0.1],...
	'Visible','off',...
	'Tag','SPMprintFootnote');
text(0,0,FNote,'FontSize',6);

H  = findobj(get(fg,'Children'),'flat','Type','axes');
un = cellstr(get(H,'Units'));
set(H,'Units','normalized')

err = 0;
%try, print -dpsc2 -painters -append -noui Anatomy.ps, catch, err=1; end

try
    
    try
        pd = defaults.ui.print;
    catch
        pd = struct('opt',{{'-dpsc2'  '-append'}},'append',true,'ext','.ps');
    end
    opts = {'Anatomy.ps','-noui','-painters',pd.opt{:}};
    fg = spm_figure('FindWin','Graphics');
    print(fg,opts{:});
catch
    print ('-dpng',  'Anatomy.png','-noui');
    spm('alert!',{'Postscript output failed';'';'Image printed as Anatomy.png' },sqrt(-1));
end





%print('-dmeta','test.emf','-noui')

if err
	errstr = lasterr;
	tmp = [find(abs(errstr)==10),length(errstr)+1];
	str = {errstr(1:tmp(1)-1)};
	for i = 1:length(tmp)-1
		if tmp(i)+1 < tmp(i+1) 
			str = [str, {errstr(tmp(i)+1:tmp(i+1)-1)}];
		end
	end
	str = {str{:},	'','- print command is:',['    print -dpsc2 -painters -append -noui Anatomy.ps'],...
			'','- current directory is:',['    ',pwd],...
			'','            * nothing has been printed *'};
	spm('alert!',str,'printing problem...',sqrt(-1));
end

if strcmp(displayType,'OL')
    delete(hAx);
    delete(hA2);
elseif strcmp(displayType,'GA')
    if isfield(st,'hAx')
    else
        delete(hAx);
        delete(hA2);
    end
elseif strcmp(displayType,'PM')
    delete(pAx);
end    
delete(fAx);

try; delete(haha);  end
try; delete(hAx);   end
try; delete(hA2);   end    
try; delete(haha2); end









case 'printflip'
%=======================================================================
%fg = get(0,'CurrentFigure');
fg = se_figure('Findwin','Graphics');
figure(fg);


    FS        = spm('FontSizes');			%-Scaled font sizes
    PF        = spm_platform('fonts');		%-Font names (for this platform)
    hA2   = axes('Position',[0.001 0.05 0.45 0.15],...
    	'DefaultTextInterpreter','Tex',...
    	'DefaultTextVerticalAlignment','Baseline',...
    	'Units','points',...
    	'Visible','off');
        AxPos = get(hA2,'Position'); set(hA2,'YLim',[0,AxPos(4)]);
    
    y     = floor(AxPos(4));
    xyzmm = se_orthviews('pos');

    text(0,y-10,[sprintf('Crosshair position: x=%1.0f, y=%1.0f, z=%1.0f,',xyzmm)],...
         'FontSize',FS(9),'Fontweight','bold')
    text(0,y-25,['Maximum (' int2str(get(st.maxima,'Value')-1) ')'],...
         'FontSize',FS(9),'Fontweight','bold')

     
     hAx   = axes('Position',[0.48 0.05 0.51 0.56],...
    	'DefaultTextInterpreter','Tex',...
    	'DefaultTextVerticalAlignment','Baseline',...
    	'Units','points',...
    	'Visible','off');
        AxPos = get(hAx,'Position'); set(hAx,'YLim',[0,AxPos(4)]);
    
    y     = floor(AxPos(4));

    text(0,y-10,['Maximum probability map: ' spm_str_manip(MAP(1).MaxMap.fname,'rt')],...
         'FontSize',FS(10))
    text(0,y-30,get(st.titel,'String'),'FontWeight','Bold','FontSize',FS(9));
    text(0,y-75,get(st.vAssign,'String'),'FontSize',FS(7));
    text(0,y-125,get(st.mAssign,'String'),'FontSize',FS(7));

    
    
    
FNote = sprintf('%s%s: %s',spm('ver'),spm('GetUser',' (%s)'),spm('time'));
%-Delete old tag lines, and print new one
delete(findobj(fg,'Tag','SPMprintFootnote'));
fAx = axes('Position',[0.005,0.005,0.1,0.1],...
	'Visible','off',...
	'Tag','SPMprintFootnote');
text(0,0,FNote,'FontSize',6);

H  = findobj(get(fg,'Children'),'flat','Type','axes');
un = cellstr(get(H,'Units'));
set(H,'Units','normalized')

err = 0;
%try, print -dpsc2 -painters -append -noui Anatomy.ps, catch, err=1; end

try
    
    try
        pd = defaults.ui.print;
    catch
        pd = struct('opt',{{'-dpsc2'  '-append'}},'append',true,'ext','.ps');
    end
    opts = {'Anatomy.ps','-noui','-painters',pd.opt{:}};
    fg = spm_figure('FindWin','Graphics');
    print(fg,opts{:});
catch
    print ('-dpng',  'Anatomy.png','-noui');
    spm('alert!',{'Postscript output failed';'';'Image printed as Anatomy.png' },sqrt(-1));
end





%print('-dmeta','test.emf','-noui')

if err
	errstr = lasterr;
	tmp = [find(abs(errstr)==10),length(errstr)+1];
	str = {errstr(1:tmp(1)-1)};
	for i = 1:length(tmp)-1
		if tmp(i)+1 < tmp(i+1) 
			str = [str, {errstr(tmp(i)+1:tmp(i+1)-1)}];
		end
	end
	str = {str{:},	'','- print command is:',['    print -dpsc2 -painters -append -noui Anatomy.ps'],...
			'','- current directory is:',['    ',pwd],...
			'','            * nothing has been printed *'};
	spm('alert!',str,'printing problem...',sqrt(-1));
end

if strcmp(displayType,'OL')
    delete(hAx);
    delete(hA2);
elseif strcmp(displayType,'GA')
    if isfield(st,'hAx')
    else
        delete(hAx);
        delete(hA2);
    end
elseif strcmp(displayType,'PM')
    delete(pAx);
end    
delete(fAx);

try; delete(haha);  end
try; delete(hAx);   end
try; delete(hA2);   end    
try; delete(haha2); end









case 'newpage'
%=======================================================================
% [hNextPage, hPrevPage, hPageNo] = se_figure('NewPage',h)
if nargin<2 | isempty(varargin{2}), error('No handles to paginate')
	else, h=varargin{2}(:)'; end

%-Work out which figure we're in
F = se_figure('ParentFig',h(1));

hNextPage = findobj(F,'Tag','NextPage');
hPrevPage = findobj(F,'Tag','PrevPage');
hPageNo   = findobj(F,'Tag','PageNo');

%-Create pagination widgets if required
%-----------------------------------------------------------------------
if isempty(hNextPage)
	WS = spm('WinScale');
	FS = spm('FontSizes');
	hNextPage = uicontrol(F,'Style','Pushbutton',...
		'HandleVisibility','on',...
		'String','>','FontSize',FS(10),...
		'ToolTipString','next page',...
		'Callback','se_figure(''TurnPage'',''+1'',gcbf)',...
		'Position',[580 020 015 015].*WS,...
		'ForegroundColor',[0 0 0],...
		'Tag','NextPage','UserData',[]);
	hPrevPage = uicontrol(F,'Style','Pushbutton',...
		'HandleVisibility','on',...
		'String','<','FontSize',FS(10),...
		'ToolTipString','previous page',...
		'Callback','se_figure(''TurnPage'',''-1'',gcbf)',...
		'Position',[565 020 015 015].*WS,...
		'Visible','on',...
		'Enable','off',...
		'Tag','PrevPage');
	hPageNo = uicontrol(F,'Style','Text',...
		'HandleVisibility','on',...
		'String','1',...
		'FontSize',FS(6),...
		'HorizontalAlignment','center',...
		'BackgroundColor','w',...
		'Position',[550 005 060 015].*WS,...
		'Visible','on',...
		'UserData',1,...
		'Tag','PageNo','UserData',1);
end

%-Add handles for this page to UserData of hNextPage
%-Make handles for this page invisible if PageNo>1
%-----------------------------------------------------------------------
mVis    = strcmp('on',get(h,'Visible'));
hPg     = get(hNextPage,'UserData');
if isempty(hPg)
	hPg = {h(mVis), h(~mVis)};
else
	hPg = [hPg; {h(mVis), h(~mVis)}];
	set(h(mVis),'Visible','off')
end
set(hNextPage,'UserData',hPg)

%-Return handles to pagination controls if requested
if nargout>0, varargout = {[hNextPage, hPrevPage, hPageNo]}; end


case 'turnpage'
%=======================================================================
% se_figure('TurnPage',move,F)
if nargin<3, F='Graphics'; else, F=varargin{3}; end
if nargin<2, move=1; else, move=varargin{2}; end
F = se_figure('FindWin',F);
if isempty(F), error('No Graphics window'), end

hNextPage = findobj(F,'Tag','NextPage');
hPrevPage = findobj(F,'Tag','PrevPage');
hPageNo   = findobj(F,'Tag','PageNo');
if isempty(hNextPage), return, end
hPg       = get(hNextPage,'UserData');
Cpage     = get(hPageNo,  'UserData');
nPages    = size(hPg,1);

%-Sort out new page number
if ischar(move), Npage = Cpage+eval(move); else, Npage = move; end
Npage = max(min(Npage,nPages),1);

%-Make current page invisible, new page visible, set page number string
set(hPg{Cpage,1},'Visible','off')
set(hPg{Npage,1},'Visible','on')
set(hPageNo,'UserData',Npage,'String',sprintf('%d / %d',Npage,nPages))

%-Disable appropriate page turning control if on first/last page (for neatness)
if Npage==1, set(hPrevPage,'Enable','off')
else, set(hPrevPage,'Enable','on'), end
if Npage==nPages, set(hNextPage,'Enable','off')
else, set(hNextPage,'Enable','on'), end



case 'deletepagecontrols'
%=======================================================================
% se_figure('DeletePageControls',F)
if nargin<2, F='Graphics'; else, F=varargin{2}; end
F = se_figure('FindWin',F);
if isempty(F), error('No Graphics window'), end

hNextPage = findobj(F,'Tag','NextPage');
hPrevPage = findobj(F,'Tag','PrevPage');
hPageNo   = findobj(F,'Tag','PageNo');

delete([hNextPage hPrevPage hPageNo])


case '#page'
%=======================================================================
% n = se_figure('#Page',F)
if nargin<2, F='Graphics'; else, F=varargin{2}; end
F = se_figure('FindWin',F);
if isempty(F), error('No Graphics window'), end

hNextPage = findobj(F,'Tag','NextPage');
if isempty(hNextPage)
	n = 1;
else
	n = size(get(hNextPage,'UserData'),1)+1;
end
varargout = {n};


case 'watermark'
%=======================================================================
% se_figure('WaterMark',F,str,Tag,Angle,Perm)
if nargin<6, HVis='on'; else, HVis='off'; end
if nargin<5, Angle=-45; else, Angle=varargin{5}; end
if nargin<4 | isempty(varargin{4}), Tag = 'WaterMark'; else, Tag=varargin{4}; end
if nargin<3 | isempty(varargin{3}), str = 'SPM';       else, str=varargin{3}; end
if nargin<2, if any(get(0,'Children')), F=gcf; else, F=''; end
	else, F=varargin{2}; end
F = se_figure('FindWin',F);
if isempty(F), return, end

%-Specify watermark color from background colour
%-----------------------------------------------------------------------
Colour = get(F,'Color');
%-Only mess with grayscale backgrounds
if ~all(Colour==Colour(1)), return, end
%-Work out colour - lighter unless grey value > 0.9
Colour = Colour+(2*(Colour(1)<0.9)-1)*0.02;

cF = get(0,'CurrentFigure');
set(0,'CurrentFigure',F)
Units=get(F,'Units');
set(F,'Units','normalized');
h = axes('Position',[0.45,0.5,0.1,0.1],...
	'Units','normalized',...
	'Visible','off',...
	'Tag',Tag);
set(F,'Units',Units)
text(0.5,0.5,str,...
	'FontSize',spm('FontSize',80),...
	'FontWeight','Bold',...
	'FontName',spm_platform('Font','times'),...
	'Rotation',Angle,...
	'HorizontalAlignment','Center',...
	'VerticalAlignment','middle',...
	'EraseMode','normal',...
	'Color',Colour,...
	'ButtonDownFcn',[...
		'if strcmp(get(gcbf,''SelectionType''),''open''),',...
			'delete(get(gcbo,''Parent'')),',...
		'end'])
set(h,'HandleVisibility',HVis)
set(0,'CurrentFigure',cF)


case 'createwin'
%=======================================================================
% F=se_figure('CreateWin',Tag,Name,Visible)

%-Condition arguments
%-----------------------------------------------------------------------
if nargin<4 | isempty(varargin{4}), Visible='on'; else, Visible=varargin{4}; end
if nargin<3, Name=''; else, Name = varargin{3}; end
if nargin<2, Tag='';  else, Tag  = varargin{2}; end

WS   = spm('WinScale');				%-Window scaling factors
FS   = spm('FontSizes');			%-Scaled font sizes
PF   = spm_platform('fonts');			%-Font names (for this platform)
Rect = spm('WinSize','Graphics','raw').*WS;	%-Graphics window rectangle
Rect = Rect+[100 50 0 -50];

F      = figure(...
	'Tag',Tag,...
	'Position',Rect,...
	'Resize','off',...
	'Color','w',...
	'ColorMap',gray(64),...
	'DefaultTextColor','k',...
	'DefaultTextInterpreter','none',...
	'DefaultTextFontName',PF.helvetica,...
	'DefaultTextFontSize',FS(10),...
	'DefaultAxesColor','w',...
	'DefaultAxesXColor','k',...
	'DefaultAxesYColor','k',...
	'DefaultAxesZColor','k',...
	'DefaultAxesFontName',PF.helvetica,...
	'DefaultPatchFaceColor','k',...
	'DefaultPatchEdgeColor','k',...
	'DefaultSurfaceEdgeColor','k',...
	'DefaultLineColor','k',...
	'DefaultUicontrolFontName',PF.helvetica,...
	'DefaultUicontrolFontSize',FS(10),...
	'DefaultUicontrolInterruptible','on',...
	'PaperType','A4',...
	'PaperUnits','normalized',...
	'PaperPosition',[.0726 .0644 .854 .870],...
	'InvertHardcopy','off',...
	'Renderer','zbuffer',...
	'Visible','off');
if ~isempty(Name)
	set(F,'Name',sprintf('%s%s',[spm('ver') ' Anatomy Toolbox'],...
		spm('GetUser',' (%s)')),'NumberTitle','off')
end
se_figure('FigContextMenu',F);
set(F,'Visible',Visible)
varargout = {F};


case 'getwinscale'
%=======================================================================
% WS = se_figure('GetWinScale')
warning('se_figure(''GetWinScale''... is Grandfathered: use spm(''WinScale''')
varargout = {spm('WinScale')};


case 'fontsizes'
%=======================================================================
% FS = se_figure('FontSizes',FS)
warning('se_figure(''FontSizes''... is Grandfathered: use spm(''FontSizes''')
if nargin<2, FS=[08,09,11,13,14,6:36]; else, FS=varargin{2}; end
varargout = {round(FS*min(spm('WinScale')))};


case 'createbar'
%=======================================================================
% se_figure('CreateBar',F)
if nargin<2, if any(get(0,'Children')), F=gcf; else, F=''; end
	else, F=varargin{2}; end
F = se_figure('FindWin',F);
if isempty(F), return, end

%-Get position and size parameters
%-----------------------------------------------------------------------
cUnits = get(F,'Units');
set(F,'Units','Pixels');
P     = get(F,'Position'); P  = P(3:4);		% Figure dimensions {pixels}
S_Gra = P./[600, 865];				% x & y scaling coefs

nBut  = 12;
nGap  = 2;
sx    = floor(P(1)./(nBut+(nGap+2)/6));		% uicontrol object width
dx    = floor(2*sx/6);				% inter-uicontrol gap
sy    = floor(20*S_Gra(1));			% uicontrol object height
x0    = dx;					% initial x position
x     = dx/2;					% uicontrol x position
y     = P(2) - sy;				% uicontrol y position
y2    = P(2) - 2.25*sy;				% uicontrol y position
FS    = round(10*min(S_Gra));			% uicontrol font size

%-Delete any existing 'ToolBar' 'Tag'ged objects
%-----------------------------------------------------------------------
cSHH = get(0,'ShowHiddenHandles');
set(0,'ShowHiddenHandles','on')
delete(findobj(F,'Tag','ToolBar'));
set(0,'ShowHiddenHandles',cSHH)

%-Create Frame for controls
%-----------------------------------------------------------------------
uicontrol(F,'Style', 'Frame',...
	'Position',[-4 (P(2) - 1.25*sy) P(1)+8 1.25*sy+4],...
	'Tag','ToolBar','HandleVisibility','callback');

%-Create uicontrol objects
%-----------------------------------------------------------------------
uicontrol(F,'String','Print','ToolTipString','print figure',...
	'Position',[x y sx sy],...
	'CallBack','se_figure(''Print'',gcf)',...
	'FontSize',FS,...
	'Interruptible','off','BusyAction','cancel',...
	'Tag','ToolBar','HandleVisibility','callback',...
	'ForegroundColor','b'); x = x+sx;

h = uicontrol(F,'String','Clear','ToolTipString','clear figure',...
	'Position',[x y sx sy],...
	'CallBack','se_figure(''Clear'',gcbf)',...
	'FontSize',FS,...
	'Interruptible','off','BusyAction','cancel',...
        'Tag','ToolBar','HandleVisibility','callback',...
	'ForegroundColor','b'); x = x+sx+dx;
if strcmp(get(F,'Tag'),'Graphics')	%-Do a full SPM interface clear
	set(h,'CallBack','spm(''Clear'',''Interactive'',gcbf)',...
		'ToolTipString','clear figure & reset SPM GUI')
end

uicontrol(F,'Style','PopUp','String',...
	'ColorMap|gray|hot|pink|gray-hot|gray-pink',...
	'ToolTipString','change colormap',...
	'Position',[x,y,2*sx,sy],...
	'CallBack',['if (get(gco,''Value'') > 1),',...
			'set(gco,''UserData'',get(gco,''Value'')),',...
			'set(gco,''Value'',1),',...
			'se_figure(''ColorMap'',get(gco,''UserData'')),',...
			'end'],...
	'FontSize',FS,...
	'Tag','ToolBar','HandleVisibility','callback',...
	'UserData','Pop'); x = x + 2*sx;

uicontrol(F,'Style','PopUp','String','Effects|invert|brighten|darken',...
	'ToolTipString','colormap effects',...
	'Position',[x,y,2*sx,sy],...
	'CallBack',['if (get(gco,''Value'') > 1),',...
			'set(gco,''UserData'',get(gco,''Value'')),',...
			'set(gco,''Value'',1),',...
			'se_figure(''ColorMap'',get(gco,''UserData'')),',...
			'end'],...
	'FontSize',FS,...
	'Tag','ToolBar','HandleVisibility','callback',...
	'UserData','Pop'); x = x + 2*sx + dx;

uicontrol(F,'String','cut','ToolTipString','delete a graphics object',...
	'Position',[x y sx sy],...
	'CallBack','se_figure(''GraphicsCut'')','FontSize',FS,...
	'Interruptible','off','BusyAction','cancel',...
	'Tag','ToolBar','HandleVisibility','callback'); x = x+sx;
uicontrol(F,'String','move','ToolTipString','move a graphics object',...
	'Position',[x y sx sy],...
	'CallBack','se_figure(''GraphicsMove'')','FontSize',FS,...
	'Interruptible','off','BusyAction','cancel',...
	'Tag','ToolBar','HandleVisibility','callback'); x = x+sx;
uicontrol(F,'String','resize','ToolTipString','resize a graphics object',...
	'Position',[x y sx sy],...
	'CallBack','se_figure(''GraphicsReSize'')','FontSize',FS,...
	'Interruptible','off','BusyAction','cancel',...
	'Tag','ToolBar','HandleVisibility','callback'); x = x+sx;
uicontrol(F,'String','text','ToolTipString','create text annotation',...
	'Position',[x y sx sy],...
	'CallBack','se_figure(''GraphicsText'')','FontSize',FS,...
	'Interruptible','off','BusyAction','cancel',...
	'Tag','ToolBar','HandleVisibility','callback'); x = x+sx;
uicontrol(F,'String','edit','ToolTipString','edit a text string',...
	'Position',[x y sx sy],...
	'CallBack','se_figure(''GraphicsTextEdit'')','FontSize',FS,...
	'Interruptible','off','BusyAction','cancel',...
	'Tag','ToolBar','HandleVisibility','callback'); x = x+sx+dx;

uicontrol(F,'String','?','ToolTipString','se_figure help',...
	'Position',[x y dx sy],...
	'CallBack','spm_help(''se_figure.m'')','FontSize',FS,...
	'Interruptible','off','BusyAction','queue',...
	'Tag','ToolBar','HandleVisibility','callback',...
	'ForegroundColor','g'); x = x+2*dx;

set(F,'Units',cUnits)


case 'figcontextmenu'
%=======================================================================
% h = se_figure('FigContextMenu',F)
if nargin<2
	F = get(0,'CurrentFigure');
	if isempty(F), error('no figure'), end
else
	F = se_figure('FindWin',varargin{2});
	if isempty(F), error('no such figure'), end
end

h = uicontextmenu('Parent',F,'HandleVisibility','CallBack');
uimenu(h,'Label','Print','HandleVisibility','CallBack',...
	'CallBack','se_figure(''Print'',gcf)')
uimenu(h,'Label','Clear','HandleVisibility','CallBack',...
	    'CallBack','se_figure(''Clear'',gcbf)');
%uimenu(h,'Label','Clear','HandleVisibility','CallBack',...
%	'CallBack','spm(''Clear'',''Interactive'',gcbf)');

hC = uimenu(h,'Label','Colormap','HandleVisibility','CallBack','Separator','on');
uimenu(hC,'Label','gray','HandleVisibility','CallBack',...
	'CallBack','se_figure(''ColorMap'',''gray'')')
uimenu(hC,'Label','hot','HandleVisibility','CallBack',...
	'CallBack','se_figure(''ColorMap'',''hot'')')
uimenu(hC,'Label','pink','HandleVisibility','CallBack',...
	'CallBack','se_figure(''ColorMap'',''pink'')')
uimenu(hC,'Label','gray-hot','HandleVisibility','CallBack',...
	'CallBack','se_figure(''ColorMap'',''gray-hot'')')
uimenu(hC,'Label','gray-pink','HandleVisibility','CallBack',...
	'CallBack','se_figure(''ColorMap'',''gray-pink'')')

hE = uimenu(h,'Label','Effects','HandleVisibility','CallBack');
uimenu(hE,'Label','invert','HandleVisibility','CallBack',...
	'CallBack','se_figure(''ColorMap'',''invert'')')
uimenu(hE,'Label','brighten','HandleVisibility','CallBack',...
	'CallBack','se_figure(''ColorMap'',''brighten'')')
uimenu(hE,'Label','darken','HandleVisibility','CallBack',...
	'CallBack','se_figure(''ColorMap'',''darken'')')

uimenu(h,'Label','cut','HandleVisibility','CallBack','Separator','on',...
	'CallBack','se_figure(''GraphicsCut'')')
uimenu(h,'Label','move','HandleVisibility','CallBack',...
	'CallBack','se_figure(''GraphicsMove'')')
uimenu(h,'Label','resize','HandleVisibility','CallBack',...
	'CallBack','se_figure(''GraphicsReSize'')')
uimenu(h,'Label','text','HandleVisibility','CallBack',...
	'CallBack','se_figure(''GraphicsText'')')
uimenu(h,'Label','edit','HandleVisibility','CallBack',...
	'CallBack','se_figure(''GraphicsTextEdit'')')

uimenu(h,'Label','help','HandleVisibility','CallBack','Separator','on',...
	'CallBack','spm_help(''se_figure.m'')')

%uimenu(h,'Label','Identify handle','HandleVisibility','CallBack','Separator','on',...
%	'CallBack','se_figure(''GraphicsHandle'')')

set(F,'UIContextMenu',h)
varargout = {h};


case 'txtcontextmenu'
%=======================================================================
% h = se_figure('TxtContextMenu',t)
if nargin<2
	F = get(0,'CurrentFigure'); t=[];
	if isempty(F), error('no figure'), end
elseif ischar(varargin{2})		%-Figure 'Tag'
	F = se_figure('FindWin',varargin{2}); t=[];
	if isempty(F), error(sprintf('no figure ''Tag''ged ''%s''',...
		varargin{2})), end
elseif ishandle(varargin{2})	%-handle: figure or text object
	switch get(varargin{2},'Type')
	case 'figure'
		F = varargin{2};
		t = [];
	case 'text'
		t = varargin{2};
		F = se_figure('ParentFig',t);
	otherwise
		error('handle not a figure or text object')
	end
else
	error('invalid handle')
end



h = uicontextmenu('Parent',F);
uimenu(h,'Label','Delete','CallBack','delete(gco)')
uimenu(h,'Label','Move','CallBack','se_figure(''GraphicsMove'')')

tmp = uimenu(h,'Label','Font','Separator','on');
uimenu(tmp,	'Label','normal','CallBack','set(gco,''FontAngle'',''normal'')')
uimenu(tmp,	'Label','italic','CallBack','set(gco,''FontAngle'',''italic'')')
uimenu(tmp,	'Label','oblique',...
		'CallBack','set(gco,''FontAngle'',''oblique'')')

uimenu(tmp,'Separator','on',...
		'Label','Helvetica','CallBack',...
		'set(gco,''FontName'',spm_platform(''font'',''helvetica''))')
uimenu(tmp,	'Label','Times','CallBack',...
		'set(gco,''FontName'',spm_platform(''font'',''times''))')
uimenu(tmp,	'Label','Courier','CallBack',...
		'set(gco,''FontName'',spm_platform(''font'',''courier''))')
uimenu(tmp,	'Label','Symbol','CallBack',...
		'set(gco,''FontName'',spm_platform(''font'',''symbol''))')

uimenu(tmp,'Separator','on',...
		'Label','light','CallBack','set(gco,''FontWeight'',''light'')')
uimenu(tmp,	'Label','normal','CallBack','set(gco,''FontWeight'',''normal'')')
uimenu(tmp,	'Label','demi','CallBack','set(gco,''FontWeight'',''demi'')')
uimenu(tmp,	'Label','bold','CallBack','set(gco,''FontWeight'',''bold'')')

tmp = uimenu(h,'Label','FontSize');
uimenu(tmp,	'Label', '6','CallBack','set(gco,''FontSize'', 6)')
uimenu(tmp,	'Label', '8','CallBack','set(gco,''FontSize'', 8)')
uimenu(tmp,	'Label','10','CallBack','set(gco,''FontSize'',10)')
uimenu(tmp,	'Label','12','CallBack','set(gco,''FontSize'',12)')
uimenu(tmp,	'Label','14','CallBack','set(gco,''FontSize'',14)')
uimenu(tmp,	'Label','16','CallBack','set(gco,''FontSize'',16)')
uimenu(tmp,	'Label','18','CallBack','set(gco,''FontSize'',18)')
uimenu(tmp,	'Label','20','CallBack','set(gco,''FontSize'',20)')
uimenu(tmp,	'Label','24','CallBack','set(gco,''FontSize'',24)')
uimenu(tmp,	'Label','28','CallBack','set(gco,''FontSize'',28)')
uimenu(tmp,	'Label','36','CallBack','set(gco,''FontSize'',36)')
uimenu(tmp,	'Label','48','CallBack','set(gco,''FontSize'',48)')
uimenu(tmp,	'Label','64','CallBack','set(gco,''FontSize'',64)')
uimenu(tmp,	'Label','96','CallBack','set(gco,''FontSize'',96)')

tmp = uimenu(h,'Label','Rotation');
uimenu(tmp,	'Label',  '0','CallBack','set(gco,''Rotation'',0)')
uimenu(tmp,	'Label', '90','CallBack','set(gco,''Rotation'',90)')
uimenu(tmp,	'Label','180','CallBack','set(gco,''Rotation'',180)')
uimenu(tmp,	'Label','-90','CallBack','set(gco,''Rotation'',270)')

uimenu(h,	'Label','Edit','CallBack','set(gco,''Editing'',''on'')')

tmp = uimenu(h,	'Label','Interpreter');
uimenu(tmp,	'Label','none','CallBack','set(gco,''Interpreter'',''none'')')
uimenu(tmp,	'Label','TeX','CallBack','set(gco,''Interpreter'',''tex'')')

uimenu(h,'Separator','on','Label','Get handle','CallBack','gco')

if ~isempty(t), set(t,'UIContextMenu',h), end
varargout = {h};


case 'colormap'
%=======================================================================
% se_figure('Colormap',ColAction,h)
if nargin<3, h=[]; else, h=varargin{3}; end
if nargin<2, ColAction='gray'; else, ColAction=varargin{2}; end
if ~ischar(ColAction)
	if ColAction==1, return, end
	Actions   = get(gcbo,'String');
	ColAction = deblank(Actions(ColAction,:));
end

switch lower(ColAction), case 'gray'
	colormap(gray(64))
case 'hot'
	colormap(hot(64))
case 'pink'
	colormap(pink(64))
case 'gray-hot'
	tmp = hot(64 + 16);  tmp = tmp([1:64] + 16,:);
	colormap([gray(64); tmp])
case 'gray-jet'
    
    txx = [0 0 .5; 
           0 0 1;
           0 1 1;
           0 1 0;           
           .5 1 0;
           1 1 0;
           1 .5 0;
           1 0 0;
           .5 0 0];
      
      
 	colormap([gray(64); interp2(1:3,[1:size(txx,1)],txx,1:3,linspace(1,size(txx,1),100)')])
case 'gray-pink'
	tmp = pink(64 + 16); tmp = tmp([1:64] + 16,:);
	colormap([gray(64); tmp])
case 'invert'
	colormap(flipud(colormap))
case 'brighten'
	colormap(brighten(colormap, 0.2))
case 'darken'
	colormap(brighten(colormap, -0.2))
otherwise
	error('Illegal ColAction specification')
end

case 'graphicshandle'
%=======================================================================
% h = se_figure('GraphicsHandle',F)
if nargin<2, F=gcbf; else, F=se_figure('FindWin',varargin{2}); end
if isempty(F), return, end

tmp = get(F,'Name');
set(F,'Name',...
	'Handle: Select item to identify, MiddleMouse=parent, RightMouse=cancel...');
set(F,'Pointer','CrossHair')
waitforbuttonpress;
h        = gco(F);
hType    = get(h,'Type');
SelnType = get(gcf,'SelectionType');
set(F,'Pointer','Arrow','Name',tmp)

if ~strcmp(SelnType,'alt') & ~isempty(h) & gcf==F
	str = sprintf('Selected (%s) object',get(h,'Type'));
	if strcmp(SelnType,'normal')
		str = sprintf('%s: handle',str);
	else
		h = get(h,'Parent');
		str = sprintf('%s: handle of parent (%s) object',str,get(h,'Type'));
	end
	if nargout==0
		assignin('base','ans',h)
		fprintf('\n%s: \n',str)
		ans = h
	else
		varargout={h};
	end
else
	varargout={[]};
end



case 'graphicscut'
%=======================================================================
% se_figure('GraphicsCut',F)
% Delete next object clicked, provided it's deletable, in this gcbf & not UIcon
% "normal" mouse button deletes object (or parent if image, line, patch, surface)
% "extend" mouse button deletes parent (if text), or object (if image, line,
%          patch, surface)
% "alt"    mouse button cancels operation

if nargin<2, F=gcbf; else, F=se_figure('FindWin',varargin{2}); end
if isempty(F), return, end

hBut  = gcbo;
set(hBut,'ForegroundColor','r')
tmp = get(F,'Name');
set(F,'Name','Cut: Select object to delete. RightMouse=cancel...');
set(F,'Pointer','Circle')
waitforbuttonpress;
h        = gco(F);
hType    = get(h,'Type');
SelnType = get(gcf,'SelectionType');
set(F,'Pointer','Arrow','Name',tmp)

if strcmp(get(h,'HandleVisibility'),'on') & ...
		~strcmp(SelnType,'alt') & ...
		~strcmp(hType,{'root','figure','uimenu','uicontrol'}) & ...
		gcf==F
	if strcmp(hType,'axes')
		delete(h)
	elseif strcmp(hType,'text')
		if strcmp(SelnType,'extend'), delete(get(h,'Parent'))
			else, delete(h), end
	elseif any(strcmp(hType,{'image','line','patch','surface'}))
		if strcmp(SelnType,'extend'), delete(h)
			else, delete(get(h,'Parent')), end
	end
end
set(hBut,'ForegroundColor','k')


case 'graphicsmove'
%=======================================================================
% se_figure('GraphicsMove',F)
% Move the next object clicked, provided it's movable and in this figure
% "normal" mouse button moves object (or parent if image, line, patch, surface)
% "extend" mouse button moves parent (if text)
% "alt"    mouse button cancels operation

if nargin<2, F=gcbf; else, F=se_figure('FindWin',varargin{2}); end
if isempty(F), return, end

hMoveBut   = gcbo;
set(hMoveBut,'ForegroundColor','r')
tmp        = get(F,'Name');
set(F,'Name','Move: MiddleMouse moves parent. RightMouse=cancel');
set(F,'Pointer','CrossHair')
waitforbuttonpress;
hPress     = gco(F);
hPressType = get(hPress,'Type');
SelnType   = get(gcf,'SelectionType');
set(F,'Pointer','Fleur','Name',tmp)

if ~strcmp(get(hPress,'HandleVisibility'),'off') & ...
		~strcmp(SelnType,'alt') & ...
		~strcmp(hPressType,{'root','figure','uimenu','uicontrol'})&...
		gcf==F
	MS.cFUnits = get(F,'Units');
	set(F,'Units','Pixels')
	MS.OPt  = get(F,'CurrentPoint');

	if ( strcmp(SelnType,'extend') & strcmp(hPressType,'text') ) | ...
			any(strcmp(hPressType,{'image','line','patch','surface'}))
		hMove = get(hPress,'Parent');
	elseif any(strcmp(hPressType,{'axes','text'}))
		hMove = hPress;
	else
		set(hMoveBut,'ForegroundColor','k')
		set(F,'Pointer','Arrow')
		return
	end

	%-Store info in UserData of hPress
	MS.hMove    = hMove;
	MS.hMoveBut = hMoveBut;
	MS.chMUnits = get(hMove,'Units');
	set(hMove,'Units','Pixels');
	MS.OPos     = get(hMove,'Position');
	MS.UserData = get(F,'UserData');
	set(hPress,'UserData',MS)

	%-Set Motion & ButtonUp functions
	set(F,'WindowButtonMotionFcn','se_figure(''GraphicsMoveMotion'')')
	set(F,'WindowButtonUpFcn','se_figure(''GraphicsMoveEnd'')')
else
	set(hMoveBut,'ForegroundColor','k')
	set(F,'Pointer','Arrow')
end


case 'graphicsmovemotion'
%=======================================================================
% se_figure('GraphicsMoveMotion')
MS = get(gco,'UserData');
set(MS.hMove,'Units','Pixels',...
	'Position',...
	MS.OPos + [get(gcf,'CurrentPoint')-MS.OPt(1:2),0*MS.OPos(3:end)])


case 'graphicsmoveend'
%=======================================================================
% se_figure('GraphicsMoveEnd')
MS = get(gco,'UserData');
hType = get(gco,'Type');
if any(strcmp(hType,{'image','line','patch','surface'}))
	set(get(gco,'Parent'),'Units',MS.chMUnits);
	set(gco,'UserData',MS.UserData);
else
	set(gco,'Units',MS.chMUnits,...
		'UserData',MS.UserData);
end

set(gcf,'Units',MS.cFUnits,...
	'WindowButtonMotionFcn','',...
	'WindowButtonUpFcn','',...
	'Pointer','Arrow')
set(MS.hMoveBut,'ForegroundColor','k')


case 'graphicsresize'
%=======================================================================
% se_figure('GraphicsReSize',F)
% Change size of next object clicked, provided it's editable and in figure
% "normal" mouse button decreases size
% "extend" mouse button increases size
% "alt"    mouse button cancels operation

if nargin<2, F=gcbf; else, F=se_figure('FindWin',varargin{2}); end
if isempty(F), return, end

hBut  = gcbo;
set(hBut,'ForegroundColor','r')
tmp   = get(F,'Name');
set(F,'Name',...
	'Resize: LeftMouse=shrink, MiddleMouse=grow, RightMouse=cancel...');
set(F,'Pointer','Circle')
waitforbuttonpress;
h        = gco(F);
hType    = get(h,'Type');
SelnType = get(gcf,'SelectionType');
u = 2*strcmp(SelnType,'extend')-1;
set(F,'Pointer','Arrow','Name',tmp)

if ~strcmp(get(h,'HandleVisibility'),'off') & ...
		~strcmp(SelnType,'alt') & ...
		~strcmp(hType,{'root','figure','uimenu','uicontrol'}) & ...
		gcf==F
	if any(strcmp(hType,{'image','line','patch','surface'}))
		h = get(h,'Parent'); hType = get(h,'Type'); end
	if strcmp(hType,'text')
		set(h,'Fontsize',(get(h,'FontSize')+2*u))
	elseif strcmp(hType,'axes')
		if u==1; u = 1.24; else, u = 1/1.24; end
		P = get(h,'Position');
		set(h,'Position',[P(1:2)-P(3:4)*(u-1)/2, P(3:4)*u])
	end
end
set(hBut,'ForegroundColor','k')



case 'graphicstext'
%=======================================================================
% se_figure('GraphicsText',F)
% Add text annotation to a figure

if nargin<2, F=gcbf; else, F=se_figure('FindWin',varargin{2}); end
if isempty(F), return, end

hBut  = gcbo;
set(hBut,'ForegroundColor','r')
tmp = get(F,'Name');
set(F,'Name',...
	'Select starting position, edit text widget. RightMouse=cancel..');
set(F,'Pointer','BotL')
waitforbuttonpress;
set(F,'Pointer','Arrow','Name',tmp)

if ~strcmp(get(gcf,'SelectionType'),'alt') & gcf==F
	cUnits = get(F,'Units');
	set(F,'Units','Normalized')
	CPt = get(F,'CurrentPoint');
	
	%-Set up axes for text widget - draw empty text object
	axes('Position',[CPt, 0.1, 0.1],'Visible','off');
	h = text(0,0,'','UIContextMenu',se_figure('TxtContextMenu'));
	
	%-Set uicontrol object and Callback
	set(F,'Units','Pixels')
	P = get(F,'Position');
	uicontrol(F,'Style','Edit',...
		'ToolTipString','enter text here',...
		'Position',[CPt(1:2).*P(3:4), (1-CPt(1))*P(3), 22],...
		'BackGroundColor',[0.8,0.8,1.0],...
		'HorizontalAlignment','Left',...
		'UserData',h,...
		'Callback',['set(get(gcbo,''UserData''),'...
			'''String'',get(gcbo,''String'')), delete(gcbo)']);
	set(F,'Units',cUnits)
end
set(hBut,'ForegroundColor','k')


case 'graphicstextedit'
%=======================================================================
% se_figure('GraphicsTextEdit',F)
% Edit text annotation to a figure

if nargin<2, F=gcbf; else, F=se_figure('FindWin',varargin{2}); end
if isempty(F), return, end

hBut     = gcbo;
set(hBut,'ForegroundColor','r')
tmp      = get(F,'Name');
set(F,'Name','Select text to edit. RightMouse=cancel..');
set(F,'Pointer','Circle')
waitforbuttonpress;
h        = gco(F);
SelnType = get(gcf,'SelectionType');
set(F,'Pointer','Arrow','Name',tmp)

if strcmp(get(h,'HandleVisibility'),'off') | ...
		strcmp(SelnType,'alt') | ...
		~strcmp(get(h,'Type'),'text') | ...
		gcf ~= F
	set(hBut,'ForegroundColor','k')
	return
end

if strcmp(SelnType,'extend')
	set(hBut,'ForegroundColor','k')
	set(h,'UIContextMenu',se_figure('TxtCOntextMenu'))
	return
end

%-Save units of various objects
cFUnits     = get(F,'Units');
cAUnits     = get(gca,'Units');
chUnits     = get(h,'Units');
chFontUnits = get(h,'FontUnits');

%-Get locations
set(F,'Units','Pixels')
set(gca,'Units','Pixels')
set(h,'Units','Pixels')
set(h,'FontUnits','Points')
tExtent = get(h,'Extent');
tmp = [1,1,0,0].*get(gca,'Position') ...
	+ [1,1,0,0].*tExtent...
	+ [0,0,1.2*max([0,0,1,1].*tExtent),ceil(22*get(h,'FontSize')/12)];

%-Create editable text widget to adjust text string
uicontrol(F,'Style','Edit',...
	'ToolTipString','edit text & press return',...
	'String',get(h,'String'),...
	'FontAngle',get(h,'FontAngle'),...
	'FontName',get(h,'FontName'),...
	'FontSize',get(h,'FontSize'),...
	'Position',tmp,...
	'BackGroundColor',[0.8,0.8,1.0],...
	'HorizontalAlignment',get(h,'HorizontalAlignment'),...
	'UserData',h,...
	'Callback',['set(get(gcbo,''UserData''),'...
		'''String'',get(gcbo,''String'')), delete(gcbo)']);

%-Reset stuff
set(hBut,'ForegroundColor','k')
set(F,'Units',cFUnits)
set(gca,'Units',cAUnits)
set(h,'Units',chUnits)
set(h,'FontUnits',chFontUnits)

otherwise
%=======================================================================
warning(['Illegal Action string: ',Action])


%=======================================================================
end
