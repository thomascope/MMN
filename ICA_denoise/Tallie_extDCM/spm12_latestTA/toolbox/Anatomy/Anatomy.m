function Anatomy(varargin)

    
global defaults
global st

    
    try
        if numel(defaults)==0
            try
                spm_defaults;
            catch
                try
                    spm('ChMod','FMRI')
                end
            end
        end
        defaults.oldDefaults = defaults;
    end



spmpath=spm('dir');
rmpath(spmpath)
addpath(spmpath)


fg = se_figure('GetWin','Graphics');
delete(fg);
fg = se_figure('GetWin','Graphics');
if isempty(fg), error('Can''t create graphics window'); end
se_figure('Clear','Graphics');
set(gcf,'DefaultUicontrolFontSize',spm('FontSizes',get(gcf,'DefaultUicontrolFontSize')));

WS = spm('WinScale');

    st.SPM = spm('FnBanner');

    uicontrol(fg,'Style','Text','Position',[25 680 550 55].*WS,'String','SPM ANATOMY TOOLBOX  v2.1',...
             'FontSize',30,'FontWeight','bold','BackgroundColor',[1 1 1],'HorizontalAlignment','Center');

   
   


uicontrol(fg,'Style','PushButton','Units','normalized','Position',[.05 .70 .9 .055],'Callback','se_explore_pmap;',...
    'String','Visualisation and statistics of cytoarchitectonic probabilistic maps','FontSize',spm('FontSizes',12));

%uicontrol(fg,'Style','PushButton','Units','normalized','Position',[.05 .70-1*0.075 .9 0.055],'Callback','se_create_MPM;',...
%    'String','Calculation of maximum probability maps','FontSize',spm('FontSizes',12));

uicontrol(fg,'Style','PushButton','Units','normalized','Position',[.05 .70-1*0.075 .9 0.055],'Callback','se_anatomy(''initP'');',...
    'String',{'Cytoarchitectonic probabilities at defined MNI coordinates'},'FontSize',spm('FontSizes',12));

uicontrol(fg,'Style','PushButton','Units','normalized','Position',[.05 .70-2*0.075 .9 0.055],'Callback','se_TabList;',...
    'String',{'Cytoarchitectonic probabilities at defined MNI coordinates (batch)'},'FontSize',spm('FontSizes',12));

uicontrol(fg,'Style','PushButton','Units','normalized','Position',[.05 .70-3*0.075 .9 0.055],'Callback','se_anatomy(''initO'');',...
    'String',{'Overlap between structure and function (SPM/images)'},'FontSize',spm('FontSizes',12));

uicontrol(fg,'Style','PushButton','Units','normalized','Position',[.05 .70-4*0.075 .9 0.055],'Callback','se_anatomy(''initG'');',...
    'String',{'Mean response (group analysis)'},'FontSize',spm('FontSizes',12));
         
uicontrol(fg,'Style','PushButton','Units','normalized','Position',[.05 .70-5*0.075 .9 0.055],'Callback','se_anatomy(''initA'');',...
    'String',{'Display functional response of anatomical areas'},'FontSize',spm('FontSizes',12));

uicontrol(fg,'Style','PushButton','Units','normalized','Position',[.05 .70-6*0.075 .9 0.055],'Callback','se_Achart;',...
    'String',{'Functional response of anatomical areas (summary)'},'FontSize',spm('FontSizes',12));

uicontrol(fg,'Style','PushButton','Units','normalized','Position',[.05 .70-7*0.075 .9 0.055],'Callback','se_createROI;',...
    'String',{'Create anatomical ROIs'},'FontSize',spm('FontSizes',12));

uicontrol(fg,'Style','PushButton','Units','normalized','Position',[.05 .70-8*0.075 .9 0.055],'Callback','se_ROIextract;',...
    'String',{'Calculate image means within anatomical ROIs'},'FontSize',spm('FontSizes',12));


uicontrol(fg,'Style','PushButton','Units','normalized','Position',[.9 .01 .05 .03],'Callback','se_exit;',...
    'String',{'Quit'},'FontSize',spm('FontSizes',10));


axis off

if nargin==0
    I = imread('InfoScreen.png');
    a = get(0,'ScreenSize');
    figure(99), image(I), axis off, set(gca,'OuterPosition',[0 0 1 1],'Position',[0 0 1 1]),
    set(gcf,'MenuBar','none','Position',[round((a(3)-round((a(4)-60)/size(I,1)*size(I,2)))/2) 60 round((a(4)-60)/size(I,1)*size(I,2)*.95) a(4)-80])
end