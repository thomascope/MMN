function se_Achart(varargin)

global group
global st
global MAP
global ArResp
global index
global displayType;

displayType = 'AC';

if (nargin==0), Action = 'init'; else, Action = varargin{1}; end

switch lower(Action), 
    case 'checked'
        for i=1:numel(st.area)
            a(i) = get(st.area(i),'Value');
        end
        if prod(size(find(a)))>20;
            spm('alert!','No more than 20 areas may be selected',sqrt(-1));
            set(gcbo,'Value',0);
        end
    case 'clear'
        if isfield(st,'figs'); delete(st.figs); st = rmfield(st,'figs');; end
        set(st.area(:),'Value',0)
    case 'display'
        PSCmin = []; PSCmax = [];
        if isfield(st,'figs'); try, delete(st.figs); st = rmfield(st,'figs'); end, end
        a = get(st.area,'Value'); for fn = 1:size(a,1); b(fn) = a{fn}; end
        if find(b)
            todo = find(b);
        	mn = size(todo,2);
        	n  = 2; if size(todo,2) < 3; n = 1; end
        	m  = ceil(mn/n); if size(todo,2) < 3; m = 2; end
        	w  = .94/n;
        	h  = .94/m;
        	ds = (w+h)*0.02;
            for i=1:prod(size(todo))
                if 1 % ArResp((todo(i))).nvox==0
                    [ArResp((todo(i))).PSC, ArResp((todo(i))).locNr, ArResp((todo(i))).nvox, ArResp((todo(i))).side] =...
                        getPSC((todo(i)),get(st.mode,'Value'));
                end
                x = .03 + w*rem(i-1,n);
                y = .96 - h*(1+(floor((i-1)/n)));
                [PSCmin(i), PSCmax(i)] = displayPSC(ArResp((todo(i))).PSC,size(ArResp((todo(i))).PSC,1),...
                    [x+ds/2 y+ds/2 w-ds (h-ds)*.94],i);
            end
            fg = se_figure('GetWin','Graphics');
            for i=1:size(todo,2)
                set(fg,'CurrentAxes',st.figs(i));
                set(st.figs(i),'XLim',[min(PSCmin),max(PSCmax)],'Box','on')
                
                if ArResp((todo(i))).side == -1;
                    Tstr = ['Left ' MAP(ArResp((todo(i))).locNr).name '  (' int2str(ArResp((todo(i))).nvox) ' voxel)'];
                else
                    Tstr = ['Right ' MAP(ArResp((todo(i))).locNr).name '  (' int2str(ArResp((todo(i))).nvox) ' voxel)'];
                end
                ArResp((todo(i))).Title = Tstr;
                title(Tstr,'VerticalAlignment','middle','FontWeight','bold');
                 if (floor((i-1)/n))<m-1 & size(todo,2)> 2
                     set(st.figs(i),'XTickLabel','');
                 end
            end
            PSC = ArResp(todo); PSC = rmfield(PSC,'side'); PSC = rmfield(PSC,'locNr');
            
            fileName = spm_input('Filename for PSC output','1','s','PSC results');
            save([fileName '.mat'],'PSC');
        else
            spm('alert!','No area selected',sqrt(-1));
        end
    case 'init'
        ArResp = struct('PSC',{},'locNr',{},'nvox',{},'side',{});
        fg = se_figure('GetWin','Graphics');
        Finter = spm('CreateIntWin','on');	
        se_figure('Clear','Graphics');
        load(spm_select(1,'mat',['Select Map'],[],spm('Dir','se_anatomy'),'MPM.',1));
		FS        = spm('FontSizes');
        hFS = FS(get(gcf,'DefaultUicontrolFontSize'));
        [B,index] = sortrows(char(MAP.name));


        st.FX = findobj(get(0,'Children'),'Flat','Tag','SATB');
        if ~isempty(st.FX)
            try
                set(0,'CurrentFigure',st.FX);
                clf
            catch
                st.FX = [];
            end
        end

        if isempty(st.FX)
            FS   = spm('FontSizes');                 PF   = spm_platform('fonts');
            a = get(0,'ScreenSize');
            st.FX = figure(...
                'Tag','SATB',                             'Position',[10 40 a(3)*.6 a(4)*.89],...
                'Resize','off',                           'MenuBar','figure',...
                'Color','w',                              'ColorMap',gray(64),...
                'DefaultTextColor','k',                   'DefaultTextInterpreter','tex',...
                'DefaultTextFontName',PF.helvetica,       'DefaultTextFontSize',FS(12),...
                'DefaultAxesColor','w',                   'DefaultAxesXColor','k',...
                'DefaultAxesYColor','k',                  'DefaultAxesZColor','k',...
                'DefaultAxesFontName',PF.helvetica,       'DefaultPatchFaceColor','k',...
                'DefaultPatchEdgeColor','k',              'DefaultSurfaceEdgeColor','k',...
                'DefaultLineColor','k',                   'DefaultUicontrolFontName',PF.helvetica,...
                'DefaultUicontrolFontSize',FS(12),        'DefaultUicontrolInterruptible','on',...
                'PaperType','A4',                         'PaperUnits','normalized',...
                'PaperPosition',[.0726 .0644 .854 .870],  'InvertHardcopy','off',...
                'Renderer','zbuffer',                     'Visible','on');
                set(st.FX,'Name',[spm('ver') ' Anatomy Toolbox Controller'],'NumberTitle','off')
        end


        uicontrol(st.FX,'Style','Frame','Units','normalized','Position',[.005 .01 .65 .98]);
        uicontrol(st.FX,'Style','Pushbutton','Units','normalized','Position',[0.02 .94 .46 .04],'ForegroundColor','r',...
            'FontWeight','bold','String','DISPLAY','Callback','se_Achart(''display'')','ToolTipString','Calculate the response (% signal change) for the selected areas');

        
        scalY = .75 / ceil(size(MAP,2)/2);
        
        st.area = [];

        for i = 1:size(MAP,2)
            if i<=ceil(size(MAP,2)/3)
                st.area((index(i)*2)-1) = uicontrol(st.FX,'Style','checkbox','Units','normalized','Position',[0.01 .903-.02*i .02 .02],'Callback','se_createROI(''checked'')'); 
                st.area((index(i)*2)) =   uicontrol(st.FX,'Style','checkbox','Units','normalized','Position',[0.03 .903-.02*i .02 .02],'Callback','se_createROI(''checked'')');
                uicontrol(st.FX,'Style','text','Units','normalized','Position',[0.05 .90-.02*i .34 .025],'String',{MAP(index(i)).name},'HorizontalAlignment','left','FontSize',hFS-5);
            elseif i<=(ceil(size(MAP,2)/3)*2)
                st.area((index(i)*2)-1) = uicontrol(st.FX,'Style','checkbox','Units','normalized','Position',[0.21 .903-.02*(i-ceil(size(MAP,2)/3)) .02 .02],'Callback','se_createROI(''checked'')'); 
                st.area((index(i)*2)) =   uicontrol(st.FX,'Style','checkbox','Units','normalized','Position',[0.23 .903-.02*(i-ceil(size(MAP,2)/3)) .02 .02],'Callback','se_createROI(''checked'')');
                uicontrol(st.FX,'Style','text','Units','normalized','Position',[0.25 .90-.02*(i-ceil(size(MAP,2)/3)) .24 .025],'String',{MAP(index(i)).name},'HorizontalAlignment','left','FontSize',hFS-5);
            else
                st.area((index(i)*2)-1) = uicontrol(st.FX,'Style','checkbox','Units','normalized','Position',[0.41 .903-.02*(i-ceil(size(MAP,2)/3)*2) .02 .02],'Callback','se_createROI(''checked'')'); 
                st.area((index(i)*2)) =   uicontrol(st.FX,'Style','checkbox','Units','normalized','Position',[0.43 .903-.02*(i-ceil(size(MAP,2)/3)*2) .02 .02],'Callback','se_createROI(''checked'')');
                uicontrol(st.FX,'Style','text','Units','normalized','Position',[0.45 .90-.02*(i-ceil(size(MAP,2)/3)*2) .2 .025],'String',{MAP(index(i)).name},'HorizontalAlignment','left','FontSize',hFS-5);
            end
            ArResp(((index(i)*2)-1)).nvox = 0; ArResp(((index(i)*2))).nvox = 0;
        end

        uicontrol(st.FX,'Style','text','Units','normalized','Position',[0.008 .905 .02 .02],'FontWeight','bold','String','L','HorizontalAlignment','center');
        uicontrol(st.FX,'Style','text','Units','normalized','Position',[0.028 .905 .02 .02],'FontWeight','bold','String','R','HorizontalAlignment','center');
        
        uicontrol(st.FX,'Style','text','Units','normalized','Position',[0.208 .905 .02 .02],'FontWeight','bold','String','L','HorizontalAlignment','center');
        uicontrol(st.FX,'Style','text','Units','normalized','Position',[0.228 .905 .02 .02],'FontWeight','bold','String','R','HorizontalAlignment','center');

        uicontrol(st.FX,'Style','text','Units','normalized','Position',[0.408 .905 .02 .02],'FontWeight','bold','String','L','HorizontalAlignment','center');
        uicontrol(st.FX,'Style','text','Units','normalized','Position',[0.428 .905 .02 .02],'FontWeight','bold','String','R','HorizontalAlignment','center');
       
        
        
        
        
        
        uicontrol(st.FX,'Style','Pushbutton','Units','normalized','Position',[.05 .02 .11 .04],...
            'String','CLEAR','Callback','se_Achart(''clear'')');
        uicontrol(st.FX,'Style','PushButton','Units','normalized','Position',[.34 .02 .11 .04],'ForegroundColor','r','FontWeight','bold',...
          	'String','EXIT','Callback','se_Achart(''exit'');','ToolTipString','quit');
        
        st.mode = uicontrol(st.FX,'Style','popupmenu','Units','normalized','Position',[0.03 .075 .44 .03],...
            'String',str2mat('Highest prob.','All assigned'),'FontSize',FS(11),'Value',1);
        


        se_getMap('groupStat');

    regressors = size(group(1).xSPM.Sess(1).Fc,2);
    
    for i=1:regressors
        switch st.SPM
            case 'SPM99'; strin = group(1).xSPM.Sess{1}.name{i};
            otherwise, strin = group(1).xSPM.Sess(1).Fc(i).name;
        end
        string = strrep(strin,'0',''); string = strrep(string,'1',''); string = strrep(string,'_',' '); string = strrep(string,'  ',' ');
        conString{i} = [int2str(i) ':' string];
    end

   uicontrol(st.FX,'Style','text','Units','normalized','Position',[.68 .85 .3 .1],'FontWeight','bold','FontSize',FS(13),'String','Regressors -> Conditions','HorizontalAlignment','center');
   uicontrol(st.FX,'Style','text','Units','normalized','Position',[.68 .85-regressors*.05 .3 regressors*.05],'FontWeight','demi','FontSize',FS(12),'String',conString,'HorizontalAlignment','center');
    
    
    uicontrol(st.FX,'Style','text','Units','normalized','Position',[.68 .35 .3 .1],'FontWeight','bold','FontSize',FS(13),'String','Significance (T - Test)','HorizontalAlignment','center');
    uicontrol(st.FX,'Style','text','Units','normalized','Position',[.68 .35-3*.05 .3 3*.05],'FontWeight','demi','FontSize',FS(12),'String',{'* p < 0.01  (unc)';'** p < 0.001  (unc)';'*** p < 0.0001  (unc)'},'HorizontalAlignment','center');
        
        
    case 'exit'
        try,  delete(st.FX), end
        se_figure('Clear','Graphics');
        clear all;
        Anatomy('select');
end



function [PSC, locNr, nvox, side] = getPSC(area,mode)
	global group
	global st
	global MAP
    global ArResp
    global SPM
    
    if rem(area,2); xyzmm = -1; locNr = (area+1)/2;
    else; xyzmm = 1; locNr = area/2;
    end; side = xyzmm.*MAP(1).orient;

    if mode == 1

        total = sum(MAP(locNr).allZ(MAP(locNr).allLR==side));
        [B I] = sort(MAP(locNr).allZ(MAP(locNr).allLR==side),'descend');
        
        XYZ = MAP(locNr).allXYZ(:,MAP(locNr).allLR==side & MAP(locNr).allZ>B(find(cumsum(B)>total/10,1,'first')));
        XYZmm = MAP(1).MaxMap.mat * [XYZ; ones(1,size(XYZ,2))];

    else
        XYZmm = MAP(locNr).XYZmm(:,sign(MAP(locNr).LR) == sign(xyzmm));
    end
    targets = XYZmm(1:3,:);

    PSC = [];
    spm_progress_bar('Init',0,'Calculating response',[MAP(locNr).name ' Subject']);
    if side == -1
        spm_progress_bar('Init',prod(size(group)),['Left ' MAP(locNr).name ': Subjects']);
    else
        spm_progress_bar('Init',prod(size(group)),['Right ' MAP(locNr).name ': Subjects']);
    end
	for i=1:prod(size(group))
        psc(i) = se_PSC_area(group(i).xSPM,targets);

        
        
        PSC = [PSC nanmean(psc(i).PSC')'];
        spm_progress_bar('Set',i);
    end
    
    
%    PSC = PSC(2:5,:);
    
    
	spm_progress_bar('Clear')
    nvox = psc(i).number;
    
function [PSCmin, PSCmax] = displayPSC(PSC,regressors,position,handleNr)
    global st
    global group
    
    fg = se_figure('GetWin','Graphics');
    st.figs(handleNr)   = axes('Position',position,'DefaultTextInterpreter','Tex',...
	'DefaultTextVerticalAlignment','Baseline','Units','points','Visible','on');
    sem = sqrt(var(PSC')/size(PSC,2));
    PSCmin = min([0 min(mean(PSC')-sem*1.2)]); PSCmax = max([0 max(mean(PSC')+sem*1.2)]); 
    barh(mean(PSC')); colormap white;
    line([0 0],[0.5 regressors+0.5],'Color','k','LineWidth',1);
    for k = 1:regressors,
       line([mean(PSC(k,:))-sem(k) mean(PSC(k,:))+sem(k)],[k k],'Color',[0 0 0],'LineWidth',1)
       line([mean(PSC(k,:))-sem(k) mean(PSC(k,:))-sem(k)],[k-.2 k+.2],'Color',[0 0 0],'LineWidth',1)
       line([mean(PSC(k,:))+sem(k) mean(PSC(k,:))+sem(k)],[k-.2 k+.2],'Color',[0 0 0],'LineWidth',1)
       a = get(st.area,'Value'); for fn = 1:size(a,1); b(fn) = a{fn}; end
       as = ceil(prod(size(find(b)))/2)*1.75;
       if abs(mean(PSC(k,:))/(var(PSC(k,:))/sqrt(size(PSC,2))))>spm_invTcdf(1 - 0.0001,size(PSC,2)-1)
            text(0,k,[int2str(k) '***'],'Color','r','FontWeight','bold','FontSize',max([18-as 1]),...
                'HorizontalAlignment','Center','VerticalAlignment','Middle') 
        elseif abs(mean(PSC(k,:))/(var(PSC(k,:))/sqrt(size(PSC,2))))>spm_invTcdf(1 - 0.001,size(PSC,2)-1)
            text(0,k,[int2str(k) '**'],'Color','r','FontWeight','bold','FontSize',max([18-as 1]),...
                'HorizontalAlignment','Center','VerticalAlignment','Middle') 
        elseif abs(mean(PSC(k,:))/(var(PSC(k,:))/sqrt(size(PSC,2))))>spm_invTcdf(1 - 0.01,size(PSC,2)-1)
            text(0,k,[int2str(k) '*'],'Color','r','FontWeight','bold','FontSize',max([18-as 1]),...
                'HorizontalAlignment','Center','VerticalAlignment','Middle') 
        else
            text(0,k,int2str(k),'Color','r','FontWeight','bold','FontSize',max([18-as 1]),...
                'HorizontalAlignment','Center','VerticalAlignment','Middle') 
        end
    end
    set(st.figs(handleNr),'YLim',[0.5 regressors+0.5],'YDir','reverse','YTick',[],'YTickLabel',''); 
    set(st.figs(handleNr),'XLim',[PSCmin,PSCmax]);
 
