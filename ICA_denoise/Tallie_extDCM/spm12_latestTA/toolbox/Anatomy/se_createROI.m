function se_createROI(varargin)

try
    warning off MATLAB:divideByZero
end

global group
global st
global MAP
global ArResp
global index
global displayType;
global defaults

st.SPM = spm('ver');

if (nargin==0), Action = 'init'; else, Action = varargin{1}; end

switch lower(Action), 
    case 'clear'
        if isfield(st,'figs'); delete(st.figs); st = rmfield(st,'figs');; end
        set(st.area(:),'Value',0)
    case 'calculate'
        a = get(st.area,'Value'); for fn = 1:size(a,1); b(fn) = a{fn}; end
        if find(b)
            se_figure('getWin','Interactive');
            ROIname = spm_input('ROI file name',1,'s','');
            fwhm = 0; %spm_input('smoothing {FWHM in mm, 0 for none}','+1','i',0);
            
            Vo = MAP(1).MaxMap;  
            Vo.dim = Vo.dim(1:3);
            Vo.dt  = [spm_type('uint8') spm_platform('bigend')];
                
            cMode = spm_input('Output space','!+1','b','Anatomical MNI | MNI',[1 2],2,'Select mode');
           
            if cMode == 2
                Vo.mat(2,4) = Vo.mat(2,4)+4;
                Vo.mat(3,4) = Vo.mat(3,4)-5;
                ROIname = [ROIname '_MNI'];
            else
                ROIname = [ROIname '_aMNI'];
            end
            
            Vo.fname = ['ROI_' ROIname '.img']; Vo.M = Vo.mat; Vo.descrip = 'Anatomy toolbox ROI';
            
            Vtmp= zeros(MAP(1).MaxMap.dim(1:3));
            todo = find(b);
            for i=1:size(todo,2)
                if rem(todo(i),2); 
                    vox = MAP((todo(i)+1)/2).XYZ(:,MAP((todo(i)+1)/2).LR == -1);
                else; 
                    vox = MAP(todo(i)/2).XYZ(:,MAP(todo(i)/2).LR == 1);
                end;
                
                A     = spm_clusters(vox); Q     = [];
                for i = 1:max(A); j = find(A == i); if length(j) >= 25; Q = [Q j]; end; end
                vox   = vox(:,Q);

                for vx=1:size(vox,2); 
                    Vtmp(vox(1,vx),vox(2,vx),vox(3,vx)) = 1; 
                end
            end
            Vo = spm_create_vol(Vo);
            for p=1:Vo.dim(3),
                Vo = spm_write_plane(Vo,Vtmp(:,:,p),p);
            end;

            if strcmpi(st.SPM,'spm5') | strcmpi(st.SPM,'spm8')
            else, try, Vo= spm_close_vol(Vo); end, end

            mat = Vo.mat; M = Vo.M;
            save(['ROI_' ROIname '.mat'],'mat','M');
            se_createROI('clear')
            se_figure('Clear','Interactive');
        else
            spm('alert!','No area selected',sqrt(-1));
        end

        
    case 'calculati'
        a = get(st.area,'Value'); for fn = 1:size(a,1); b(fn) = a{fn}; end
        if find(b)
            se_figure('getWin','Interactive');
            fwhm = 0; %spm_input('smoothing {FWHM in mm, 0 for none}','+1','i',0);
            
            Vo = MAP(1).MaxMap;  
            Vo.dim = Vo.dim(1:3);
            Vo.dt  = [spm_type('uint8') spm_platform('bigend')];
                
            cMode = spm_input('Output space','!+1','b','Anatomical | MNI',[1 2],2,'Select mode');
            
            if cMode == 2
                Vo.mat(2,4) = Vo.mat(2,4)+4;
                Vo.mat(3,4) = Vo.mat(3,4)-5;
            end

            Vtmp= zeros(MAP(1).MaxMap.dim(1:3));
            todo = find(b);
            for i=1:size(todo,2)
                Vtmp= zeros(MAP(1).MaxMap.dim(1:3));
                if rem(todo(i),2);
                    vox = MAP((todo(i)+1)/2).XYZ(:,MAP((todo(i)+1)/2).LR == -1);
                    ROIname = [spm_str_manip(strrep(strrep(MAP((todo(i)+1)/2).ref,'/',filesep'),'\',filesep),'rt') '_L'];
                else; 
                    vox = MAP(todo(i)/2).XYZ(:,MAP(todo(i)/2).LR == 1);
                    ROIname = [spm_str_manip(strrep(strrep(MAP((todo(i))/2).ref,'/',filesep'),'\',filesep),'rt') '_R'];
                end;
                
                
                
                if cMode == 2
                    ROIname = [ROIname '_MNI'];
                else
                    ROIname = [ROIname '_aMNI'];
                end
                Vo.fname = ['ROI_' ROIname '.nii']; Vo.M = Vo.mat; Vo.descrip = 'Anatomy toolbox ROI';


                A     = spm_clusters(vox); Q     = [];
                for i = 1:max(A); j = find(A == i); if length(j) >= 25; Q = [Q j]; end; end
                vox   = vox(:,Q);

                for vx=1:size(vox,2);
                    Vtmp(vox(1,vx),vox(2,vx),vox(3,vx)) = 1;
                end
                
                try
                    Vo = spm_write_vol(Vo,Vtmp);
                catch
                    Vo = spm_create_vol(Vo);
                for p=1:Vo.dim(3),
                    Vo = spm_write_plane(Vo,Vtmp(:,:,p),p);
                end;
                end
            end



            if strcmpi(st.SPM,'spm5') | strcmpi(st.SPM,'spm8')
            else, try, Vo= spm_close_vol(Vo); end, end



            se_createROI('clear')
            se_figure('Clear','Interactive');
        else
            spm('alert!','No area selected',sqrt(-1));
        end
        
        
    case 'init'
        fg = se_figure('GetWin','Graphics');
        Finter = spm('CreateIntWin','on');	
        se_figure('Clear','Graphics');
        load(spm_select(1,'mat',['Select Map'],[],spm('Dir','se_anatomy'),'MPM.',1));
		FS        = spm('FontSizes');
		hFS = FS(get(gcf,'DefaultUicontrolFontSize'));
        [B,index] = sortrows(char(MAP.name));
        
        uicontrol(fg,'Style','Frame','Units','normalized','Position',[0 .85-.02*ceil(size(MAP,2)/4) 1 .115+.02*(ceil(size(MAP,2)/4))]);
        uicontrol(fg,'Style','Pushbutton','Units','normalized','Position',[0.01 .92 .31 .04],'ForegroundColor','r',...
            'FontWeight','bold','String','Create composite','Callback','se_createROI(''calculate'')','ToolTipString','Creates an ROI containing all selected areas');
        uicontrol(fg,'Style','Pushbutton','Units','normalized','Position',[0.33 .92 .31 .04],'ForegroundColor','r',...
            'FontWeight','bold','String','Create individual','Callback','se_createROI(''calculati'')','ToolTipString','Creates ROI or each if the selected areas');

        uicontrol(fg,'Style','Pushbutton','Units','normalized','Position',[.64 .92 .15 .04],...
            'String','CLEAR','Callback','se_createROI(''clear'')');
        uicontrol(fg,'Style','PushButton','Units','normalized','Position',[.79 .92 .15 .04],'ForegroundColor','r','FontWeight','bold',...
          	'String','EXIT','Callback','se_createROI(''exit'');','ToolTipString','quit');

        
        for i = 1:size(MAP,2)
            if i<=ceil(size(MAP,2)/4)
                st.area((index(i)*2)-1) = uicontrol(fg,'Style','checkbox','Units','normalized','Position',[0.01 .895-.02*i .03 .02],'Callback','se_createROI(''checked'')'); 
                st.area((index(i)*2)) =   uicontrol(fg,'Style','checkbox','Units','normalized','Position',[0.04 .895-.02*i .03 .02],'Callback','se_createROI(''checked'')');
                uicontrol(fg,'Style','text','Units','normalized','Position',[0.08 .89-.02*i .17 .025],'String',{MAP(index(i)).name},'HorizontalAlignment','left','FontSize',hFS-1);
            elseif i<=(ceil(size(MAP,2)/4)*2)
                st.area((index(i)*2)-1) = uicontrol(fg,'Style','checkbox','Units','normalized','Position',[0.26 .895-.02*(i-ceil(size(MAP,2)/4)) .03 .02],'Callback','se_createROI(''checked'')'); 
                st.area((index(i)*2)) =   uicontrol(fg,'Style','checkbox','Units','normalized','Position',[0.29 .895-.02*(i-ceil(size(MAP,2)/4)) .03 .02],'Callback','se_createROI(''checked'')');
                uicontrol(fg,'Style','text','Units','normalized','Position',[0.33 .89-.02*(i-ceil(size(MAP,2)/4)) .17 .025],'String',{MAP(index(i)).name},'HorizontalAlignment','left','FontSize',hFS-1);
            elseif i<=(ceil(size(MAP,2)/4)*3)
                st.area((index(i)*2)-1) = uicontrol(fg,'Style','checkbox','Units','normalized','Position',[0.51 .895-.02*(i-ceil(size(MAP,2)/4)*2) .03 .02],'Callback','se_createROI(''checked'')'); 
                st.area((index(i)*2)) =   uicontrol(fg,'Style','checkbox','Units','normalized','Position',[0.54 .895-.02*(i-ceil(size(MAP,2)/4)*2) .03 .02],'Callback','se_createROI(''checked'')');
                uicontrol(fg,'Style','text','Units','normalized','Position',[0.58 .89-.02*(i-ceil(size(MAP,2)/4)*2) .17 .025],'String',{MAP(index(i)).name},'HorizontalAlignment','left','FontSize',hFS-1);
            else
                st.area((index(i)*2)-1) = uicontrol(fg,'Style','checkbox','Units','normalized','Position',[0.76 .895-.02*(i-ceil(size(MAP,2)/4)*3) .03 .02],'Callback','se_createROI(''checked'')'); 
                st.area((index(i)*2)) =   uicontrol(fg,'Style','checkbox','Units','normalized','Position',[0.79 .895-.02*(i-ceil(size(MAP,2)/4)*3) .03 .02],'Callback','se_createROI(''checked'')');
                uicontrol(fg,'Style','text','Units','normalized','Position',[0.83 .89-.02*(i-ceil(size(MAP,2)/4)*3) .16 .025],'String',{MAP(index(i)).name},'HorizontalAlignment','left','FontSize',hFS-1);
            end
            ArResp(((index(i)*2)-1)).nvox = 0; ArResp(((index(i)*2))).nvox = 0;
        end
        
        uicontrol(fg,'Style','text','Units','normalized','Position',[0.01 .895 .029 .02],'FontWeight','bold','String','L','HorizontalAlignment','center');
        uicontrol(fg,'Style','text','Units','normalized','Position',[0.03 .895 .029 .02],'FontWeight','bold','String','R','HorizontalAlignment','center');
        
        uicontrol(fg,'Style','text','Units','normalized','Position',[0.26 .895 .029 .02],'FontWeight','bold','String','L','HorizontalAlignment','center');
        uicontrol(fg,'Style','text','Units','normalized','Position',[0.29 .895 .029 .02],'FontWeight','bold','String','R','HorizontalAlignment','center');

        uicontrol(fg,'Style','text','Units','normalized','Position',[0.51 .895 .029 .02],'FontWeight','bold','String','L','HorizontalAlignment','center');
        uicontrol(fg,'Style','text','Units','normalized','Position',[0.53 .895 .029 .02],'FontWeight','bold','String','R','HorizontalAlignment','center');

        uicontrol(fg,'Style','text','Units','normalized','Position',[0.76 .895 .029 .02],'FontWeight','bold','String','L','HorizontalAlignment','center');
        uicontrol(fg,'Style','text','Units','normalized','Position',[0.79 .895 .029 .02],'FontWeight','bold','String','R','HorizontalAlignment','center');

        
        
       
    case 'exit'
        se_figure('Clear','Graphics');
        clear all;
        Anatomy('select');
end


