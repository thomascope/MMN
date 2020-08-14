fg = se_figure('GetWin','Graphics');
se_figure('Clear','Graphics');
WS = spm('WinScale');
uicontrol(fg,'Style','Text','Position',[10 780 590 55].*WS,'String','SPM ANATOMY TOOLBOX',...
             'FontSize',25,'FontWeight','bold','BackgroundColor',[1 1 1],'HorizontalAlignment','Center');

uicontrol(fg,'Style','Text','Position',[210 760 160 25].*WS,'String','Version 1.5',...
             'FontSize',14,'BackgroundColor',[1 1 1]);
         
         
         
uicontrol(fg,'Style','Text','Position',[80 715 160 25].*WS,'String','written by:',...
             'FontSize',10,'BackgroundColor',[1 1 1]);

uicontrol(fg,'Style','Text','Position',[200 600 400 150].*WS,'String',...
    {'Simon Eickhoff  (s.eickhoff@fz-juelich.de)';
     'Institut for Medicine (IME), Research Center Jülich'},'FontSize',10,'BackgroundColor',[1 1 1],...
     'HorizontalAlignment','left');
         
uicontrol(fg,'Style','Text','Position',[50 595 490 100].*WS,'String','Primary references:','FontWeight','bold','FontSize',13,'BackgroundColor',[1 1 1]);
uicontrol(fg,'Style','Text','Position',[50 520 500 150].*WS,'String',...
{'Eickhoff SB et al.:  A new SPM toolbox for combining probabilistic cytoarchitectonic maps and functional imaging data. (2005) NeuroImage 25(4): 1325-1335';...
'Eickhoff SB et al.:  Testing anatomically specified hypotheses in functional imaging using cytoarchitectonic maps. (2006) NeuroImage 32(2): 570-82';...
'Eickhoff SB et al.:  Assignment of functional activations to probabilistic cytoarchitectonic areas revisited. (2007) NeuroImage 36(3): 511-521 '},...
  'FontSize',10,'BackgroundColor',[1 1 1],...
     'HorizontalAlignment','left');


 
 uicontrol(fg,'Style','Text','Position',[10 480 590 55].*WS,'String','Publications describing included cytoarchitectonic maps:',...
             'FontSize',13,'BackgroundColor',[1 1 1],'FontWeight','bold');

uicontrol(fg,'Style','Text','Position',[20 60 200 445].*WS,'String',...
    {'Auditory cortex';
     '';
     'Broca''s area';
     '';
     'Motor cortex';
     '';
     '';
     'Somatosensory cortex';
     ''
     ''
     'Parietal operculum / SII'
     '';
     'Amygdala'
     'Hippocampus';
     '';
     'anterior intraparietal sulcus'
     '';
     'Visual cortex';
     ''; '';
     'Fiber tracts'},'FontSize',10,'BackgroundColor',[1 1 1],...
     'HorizontalAlignment','left');
 
uicontrol(fg,'Style','Text','Position',[180 60 200 445].*WS,'String',...
    {'TE 1.0, TE 1.1, TE 1.2';
     '';
     'BA 44, BA 45';
     '';
     'BA 4a, BA 4p';
     'BA 6';
     '';
     'BA 3a, BA 3b, BA 1';
     'BA 2';
     ''
     'OP 1, OP 2, OP 3, OP 4'
     '';
     'CM/LB/SF'
     'FD/CA/SUB/EC/HATA';
     '';
     'hIP1, hIP2'
     '';
     'BA 17, BA 18';
     'hOC5'; '';
     'ar, cb, cing, ct, forn, iof, lgb';'mb, mgb, or, slf, sof, uf'},'FontSize',10,'BackgroundColor',[1 1 1],...
     'HorizontalAlignment','left');
 
uicontrol(fg,'Style','Text','Position',[350 60 300 445].*WS,'String',...
    {'Morosan et al., NeuroImage 2001';
     '';
     'Amunts et al., J Comp Neurol 1999';
     '';
     'Geyer et al., Nature 1996';
     'S. Geyer, Springer press 2003';
     '';
     'Geyer et al., NeuroImage, 1999, 2000';
     'Grefkes et al., NeuroImage 2001';
     '';
     'Eickhoff et al., Cerebral Cortex 2006a,b'
     '';
     'Amunts et al., Anat Embryol 2005'
     'Amunts et al., Anat Embryol 2005'
     '';
     'Choi et al., J Comp Neurol 2006'
     '';
     'Amunts et al., NeuroImage 2000';
     'Malikovic et al., Cerebral Cortex 2006';
     '';'Bürgel et al., NeuroImage 1999, 2006'},'FontSize',10,'BackgroundColor',[1 1 1],...
     'HorizontalAlignment','left');

 uicontrol(fg,'Style','Text','Position',[40 12 500 55].*WS,'String','Other areas may only be used with authors'' permission !',...
             'FontSize',12,'BackgroundColor',[1 1 1],'FontWeight','bold');

uicontrol(fg,'Style','PushButton', 'Position',[210 10 160 28].*WS,'String','Start','FontSize',20,...
	'Callback','Anatomy(''select'')','FontWeight','bold','FontSize',14,'ForegroundColor',[0 .5 0]);
