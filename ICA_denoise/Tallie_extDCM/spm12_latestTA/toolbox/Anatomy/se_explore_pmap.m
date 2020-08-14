function se_explore_pmap(op,varargin)

global st
global PMap
global displayType
global PMap
global defaults

displayType = 'PM';

if nargin==0,
    op = 'init';
end    
   
switch op
    
    case 'init0'
        try, PMap = rmfield(PMap,'BG'); end
        se_explore_pmap('init');
    case 'initM'
        PMap.mat = NaN;
        se_explore_pmap('init')
    case 'init'
        
        try
            try,
                st.flip = defaults.analyze.flip; defaults.analyze.flip = 0;
            catch;
                try
                    spm('ChMod','FMRI')
                catch
                    spm_defaults;
                end
                st.flip = defaults.analyze.flip; defaults.analyze.flip = 0;
            end
        end
        
spm_figure('GetWin','Interactive');


try 
    PMap.mat;
    if isnan(PMap.mat)
        MAP = spm_vol(spm_select(1,'image','Select PMap to analyze',[],[spm('Dir','se_anatomy') filesep 'PMaps']));
    end
catch
    MAP = spm_vol(spm_select(1,'image','Select PMap to analyze',[],[spm('Dir','se_anatomy') filesep 'PMaps']));
end

try
    BG = PMap.BG;
catch
    BG    = spm_vol(spm_select(1,'image','Select Background (reference) image',[],spm('Dir','se_anatomy')));
    BG.orig  = spm_input('Background origin AC ?',-1,'y/n',[1,2],1);
end

try
    if BG.orig == 2
        st.move = 1;
    else
        st.move = 0;
    end
catch
    BG.orig  = spm_input('Background origin corrected to AC ?',-1,'y/n',[1,2],1);
    if BG.orig == 2
        st.move = 1;
    else
        st.move = 0;
    end
end
    
try, PMap = MAP; end
PMap.BG = BG;

AM = spm_vol(fullfile(spm('Dir','se_anatomy'), 'AnatMask.nii'));

PMap.XYZ = [1;1;1]; PMap.Z = [0];
%spm_progress_bar('Init',PMap.dim(3),'Preparing data');
for p = 1:PMap.dim(3)
    [i,j,z] = find(spm_slice_vol(PMap,[1 0 0 0; 0 1 0 0; 0 0 1 p; 0 0 0 1],PMap.dim(1:2),[0, NaN]));
    if any(i)
            PMap.XYZ = [PMap.XYZ [i'; j'; p*ones(1,length(i))]]; PMap.Z = [PMap.Z z'*250];
    end
%    spm_progress_bar('Set',p);
end
%spm_progress_bar('Clear')

LR = round(spm_sample_vol(AM,PMap.XYZ(1,:),PMap.XYZ(2,:),PMap.XYZ(3,:),0));


PMap.XYZ_li = PMap.XYZ(:,LR == 2);
PMap.XYZ_re = PMap.XYZ(:,LR == 1);

PMap.Z_li = PMap.Z(:,LR == 2);
PMap.Z_re = PMap.Z(:,LR == 1);


PMap.XYZmm    = PMap.mat * [PMap.XYZ; ones(1,size(PMap.XYZ,2))]; 
PMap.XYZmm_li = PMap.mat * [PMap.XYZ_li; ones(1,size(PMap.XYZ_li,2))];
PMap.XYZmm_re = PMap.mat * [PMap.XYZ_re; ones(1,size(PMap.XYZ_re,2))];

PMap.cog_li = round(PMap.XYZmm_li * [PMap.Z_li'])/sum(PMap.Z_li);
PMap.cog_re = round(PMap.XYZmm_re * [PMap.Z_re'])/sum(PMap.Z_re);

fg = se_figure('GetWin','Graphics');
if isempty(fg), error('Can''t create graphics window'); end
se_figure('Clear','Graphics');

WS = spm('WinScale');

se_orthviews('Reset');
se_orthviews('Image', BG.fname, [0.1 -0.04 .7 .7]);     % Image appears
se_orthviews('AddBlobs',1,PMap.XYZ,(PMap.Z/25)*10,PMap.mat);
pos = PMap.XYZmm(:,find(PMap.Z==max(PMap.Z))); pos = pos(:,1);
se_orthviews('Reposition',pos(1:3))

st.callback = 'se_explore_pmap(''shopos'');';

hAx   = axes('Position',[0.03 0.56 0.97 0.44],...
	'DefaultTextInterpreter','Tex',...
	'DefaultTextVerticalAlignment','Baseline',...
	'DefaultTextHorizontalAlignment','center',...
	'Units','points',...
	'Visible','off');
    AxPos = get(hAx,'Position'); 
    set(hAx,'YLim',[0,AxPos(4)]); set(hAx,'XLim',[0,AxPos(3)]);


text(10*WS(3), 220*WS(4),strrep(spm_str_manip(PMap.fname,'rt'),'_',':  '),...
             'FontSize',25,'FontWeight','bold','HorizontalAlignment','left');

text(40*WS(3), 190*WS(4),'Probability','FontWeight','bold','FontSize',12);
 
text(120*WS(3), 190*WS(4),'Voxel = mm3','FontWeight','bold','FontSize',12);
 
text(98*WS(3), 170*WS(4),'left','FontWeight','bold','FontSize',10);

text(142*WS(3), 170*WS(4),'right','FontWeight','bold','FontSize',10);


text(50*WS(3),(167-0*16-8)*WS(4),[' 0%'],'FontSize',10);
for i=1:10
    text(50*WS(3),(167-i*16-8)*WS(4),[int2str(i*10) '%'],'FontSize',10);
        
    text(98*WS(3), (167-i*16)*WS(4),[int2str( sum(PMap.Z_li>((i-1)*25) & PMap.Z_li<(i*25+1))  )],'FontSize',10);
        
    text(142*WS(3), (167-i*16)*WS(4),[int2str( sum(PMap.Z_re>((i-1)*25) & PMap.Z_re<(i*25+1)) )],'FontSize',10);
end

text(215*WS(3), 160*WS(4),...
             {['Center: '];[]; 
              ['Minimum: '];
              ['Maximum: '];[];
              ['50% min.'];
              ['50% max.']},...
             'FontWeight','bold','FontSize',10);
         

text(290*WS(3), 190*WS(4),'left','FontWeight','bold','FontSize',12);

text(375*WS(3), 190*WS(4),'right','FontWeight','bold','FontSize',12);

text(265*WS(3), 175*WS(4),'x','FontWeight','bold','FontSize',10);         
text(290*WS(3), 175*WS(4),'y','FontWeight','bold','FontSize',10); 
text(315*WS(3), 175*WS(4),'z','FontWeight','bold','FontSize',10);

text(350*WS(3), 175*WS(4),'x','FontWeight','bold','FontSize',10);      
text(375*WS(3), 175*WS(4),'y','FontWeight','bold','FontSize',10);        
text(400*WS(3), 175*WS(4),'z','FontWeight','bold','FontSize',10);       

for i=1:3
  text((265+((i-1)*25))*WS(3), 160*WS(4),...
    {int2str(PMap.cog_li(i));[]; 
     int2str(min(PMap.XYZmm_li(i,:)));
    int2str(max(PMap.XYZmm_li(i,:)));[];
    int2str(min(PMap.XYZmm_li(i,PMap.Z_li>=125)));
    int2str(max(PMap.XYZmm_li(i,PMap.Z_li>=125)))},...
    'FontSize',10);         
  text((350+((i-1)*25))*WS(3), 160*WS(4),...
    {int2str(PMap.cog_re(i));[]; 
     int2str(min(PMap.XYZmm_re(i,:)));
    int2str(max(PMap.XYZmm_re(i,:)));[];
    int2str(min(PMap.XYZmm_re(i,PMap.Z_re>=125)));
    int2str(max(PMap.XYZmm_re(i,PMap.Z_re>=125)))},...
    'FontSize',10);         
end

uicontrol(fg,'Style','PushButton','Position',[5 750 200 35].*WS,'String','Change Map','Callback','se_explore_pmap(''initM'')',...
             'ToolTipString','Change the displayed map');

uicontrol(fg,'Style','PushButton','Position',[215 750 260 35].*WS,'String','Change Background','Callback','se_explore_pmap(''init0'')',...
             'ToolTipString','Change the background (reference) image');

uicontrol(fg,'Style','PushButton','Position',[495 750 100 35].*WS,'String','EXIT','Callback','Anatomy(1)',...
             'ToolTipString','Back to main menue','ForegroundColor','r','FontWeight','bold');

%x_shift = 265; y_shift = 353;
x_shift = 330; y_shift = -30;
uicontrol(fg,'Style','Frame','Position',[45+x_shift 125+y_shift 190 82].*WS);
uicontrol(fg,'Style','Frame','Position',[50+x_shift  130+y_shift 180 72].*WS);

uicontrol(fg,'Style','Text', 'Position',[55+x_shift  183+y_shift 170 018].*WS,'String','Crosshair Position','ToolTipString','In anatomical MNI space');
uicontrol(fg,'Style','PushButton', 'Position',[55+x_shift  174+y_shift 170 006].*WS,...
	'Callback','se_orthviews(''Reposition'',[0 0 0]);','ToolTipString','move crosshairs to origin');
uicontrol(fg,'Style','Text', 'Position',[55+x_shift  155+y_shift 35 020].*WS,'String','mm:','ToolTipString','In anatomical MNI space');
uicontrol(fg,'Style','Text', 'Position',[55+x_shift 132+y_shift 65 020].*WS,'String','Probability:');

st.mp = uicontrol(fg,'Style','edit', 'Position',[90+x_shift 155+y_shift 135 020].*WS,'String','','Callback','se_explore_pmap(''setposmm'')','ToolTipString','move crosshairs to mm coordinates (in anatomical MNI space)');
st.in = uicontrol(fg,'Style','Text', 'Position',[120+x_shift 132+y_shift  85 020].*WS,'String','');


se_explore_pmap('shopos')



case 'shopos'
    
	% The position of the crosshairs has been moved.
	%-----------------------------------------------------------------------
	if isfield(st,'mp'),
		fg = se_figure('Findwin','Graphics');
		if any(findobj(fg) == st.mp),
		set(st.mp,'String',sprintf('%.0f  %.0f  %.0f',se_orthviews('pos')));
		pos = se_orthviews('pos');
        pos = inv(PMap.mat )* [pos;1];
        
        
        intX = spm_sample_vol(PMap,pos(1),pos(2),pos(3),0)*100;
        if intX>0
            set(st.in,'String',[int2str(intX) '%']);
        else
            set(st.in,'String',['0%']);
        end
		else,
			st.Callback = ';';
			rmfield(st,{'mp','vp','in'});
		end;
	else,
		st.Callback = ';';
	end;
	return;
end;

if strcmp(op,'setposmm'),
	% Move the crosshairs to the specified position
	%-----------------------------------------------------------------------
	if isfield(st,'mp'),
		fg = se_figure('Findwin','Graphics');
		if any(findobj(fg) == st.mp),
			pos = sscanf(get(st.mp,'String'), '%g %g %g');
			if length(pos)~=3,
				pos = se_orthviews('pos');
			end;
			se_orthviews('Reposition',pos);
		end;
	end;
	return;
end;
