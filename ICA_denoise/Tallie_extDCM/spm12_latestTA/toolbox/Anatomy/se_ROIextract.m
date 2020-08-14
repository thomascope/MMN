function se_ROIextract(varargin)

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
        a = get(st.area,'Value'); for fn = 1:size(a,1); b(fn) = a{fn}; end
    case 'clear'
        try
        if isfield(st,'figs'); delete(st.figs); st = rmfield(st,'figs');; end
        end
        set(st.area(:),'Value',0)
    case 'all'
        if isfield(st,'figs'); delete(st.figs); st = rmfield(st,'figs');; end
        set(st.area(:),'Value',1)
    case 'display'

        OUT = struct('Y',{},'YMean',{},'Name',{});
        Vi = spm_vol(spm_select(Inf,'image','Select images to process'));

        orig = spm_input('Are the origins already corrected to AC ?','+1','y/n',[1,2],2);
        if orig == 2
            xtr = 0; ytr = -4; ztr = +5;
            for im = 1:size(Vi,1)
                Vi(im).mat(1,4) = Vi(im).mat(1,4) + xtr;
                Vi(im).mat(2,4) = Vi(im).mat(2,4) + ytr;
                Vi(im).mat(3,4) = Vi(im).mat(3,4) + ztr;
            end
        end

        
        [p,YPos] = spm_input('Data adjustment','+1','b',{'Scalar','Global','None'},[1 2 0]);
        
        if p == 1
                pmult  = spm_input(['Premultiply data by '],'+0','r',1);
                while numel(pmult) ~= 1 & numel(pmult) ~= numel(Vi)
                    spm('alert',['Please enter either a scalar or a ' int2str(numel(Vi)) ' element vector '])
                    pmult  = spm_input(['Premultiply data by '],'+0','r',1);
                end
                if numel(pmult) == 1
                    pmult = repmat(pmult,1,numel(Vi));
                end
                
        elseif p == 2
            fprintf('%-40s: %30s','Calculating globals',' ')                     %-#
            for i = 1:numel(Vi)
                fprintf('%s%30s',repmat(sprintf('\b'),1,30),sprintf('%4d/%-4d',i,numel(Vi))) %-#
                g(i) = spm_global(Vi(i));
            end
            GM      = 1; % mean(g)
            pmult   = GM./g;
            
        else
        	pmult = repmat(1,1,numel(Vi));
            
        end

        OUTPUTname = spm_input('Name of group for file ouptut','!+1','s','Group_1');

        a = get(st.area,'Value'); for fn = 1:size(a,1); b(fn) = a{fn}; end
        if find(b)
            todo = find(b);
            cnt = 1;
            se_figure('GetWin','Interactive');
            se_figure('Clear','Interactive');
            spm_progress_bar('Init',prod(size(todo)),'Calculating ROI means');
            for area = todo
                if rem(area,2); xyzmm = -1; locNr = (area+1)/2; 
                else; xyzmm = 1; locNr = area/2;
                end; side = xyzmm;
                
                modus = get(st.mode,'Value');
                if modus == 1

                    total = sum(MAP(locNr).allZ(MAP(locNr).allLR==side));
                    [B I] = sort(MAP(locNr).allZ(MAP(locNr).allLR==side),'descend');
                    
                    XYZ = MAP(locNr).allXYZ(:,MAP(locNr).allLR==side & MAP(locNr).allZ>B(find(cumsum(B)>total/10,1,'first')));
                    XYZmm = MAP(1).MaxMap.mat * [XYZ; ones(1,size(XYZ,2))];

                else
                    XYZmm = MAP(locNr).XYZmm(:,sign(MAP(locNr).LR) == sign(xyzmm));
                end


                Ytmp1 = []; Ytmp2 = []; XYZmm = XYZmm(1:3,:);
                for im = 1:size(Vi,1)
                    XYZ = inv(Vi(im).mat) * [XYZmm; ones(1,size(XYZmm,2))];
                    Y = spm_sample_vol(Vi(im),XYZ(1,:),XYZ(2,:),XYZ(3,:),1)*pmult(im);
                    Ytmp1 = [Ytmp1; mean(Y(~isnan(Y)))];
                    Ytmp2 = [Ytmp2; median(Y(~isnan(Y)))];
                end


                OUT(cnt).means   = Ytmp1;
                OUT(cnt).medians = Ytmp2;
                OUT(cnt).YMean = mean(Ytmp1);
                if xyzmm == -1;
                    OUT(cnt).Name = [MAP(locNr).name ' L'];
                else
                    OUT(cnt).Name = [MAP(locNr).name ' R'];
                end
                spm_progress_bar('Set',cnt);
                cnt = cnt+1;
            end
            spm_progress_bar('Clear');

            if ispc
                fid = fopen([OUTPUTname '.txt'],'wt');
            else
                fid = fopen([OUTPUTname '.txt'],'w+');
            end

            if modus == 1
                fprintf(fid,'%s\n\n\n\r','ROIs defined by to 10-percentile of the probability maps');
            else
                fprintf(fid,'%s\n\n\n\r','ROIs defined by the MPM');
            end
            for roi = 1:size(OUT,2)
                fprintf(fid,'%s\t',OUT(roi).Name);
                fprintf(fid,[repmat('%-10.4f\t',1,size(OUT(roi).means,1)) '\n'],OUT(roi).means');
%                fprintf(fid,'%s\n\n\n\r',['Mean value: ' num2str(OUT(roi).YMean,'%10.4f')]);
            end

            fprintf(fid,'\n\n\n\r%s\n\r','Files used:');
            for roi = 1:size(Vi,1)
                fprintf(fid,'%s\n\r',Vi(roi).fname);
            end


            status = fclose(fid);
            se_ROIextract('clear')
            
            save([OUTPUTname '.mat'],'OUT');
            msgbox(['File ' OUTPUTname '.txt written to current directory (' pwd ')'],'SPM Anatomy Toolbox','warn');
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



        uicontrol(fg,'Style','Frame','Units','normalized','Position',[0.005 .76-.02*ceil(size(MAP,2)/2) .99 .205+.02*(ceil(size(MAP,2)/2))]);
        uicontrol(fg,'Style','Pushbutton','Units','normalized','Position',[0.01 .93 .2 .03],'ForegroundColor','r',...
            'FontWeight','bold','String','Compute','Callback','se_ROIextract(''display'')','ToolTipString','Calculate the mean signal for the selected areas within a series of images');

     
        uicontrol(fg,'Style','Pushbutton','Units','normalized','Position',[.525 .93 .15 .03],...
            'String','ALL','FontWeight','bold','Callback','se_ROIextract(''all'')');
        uicontrol(fg,'Style','Pushbutton','Units','normalized','Position',[.675 .93 .18 .03],...
            'String','CLEAR','FontWeight','bold','Callback','se_ROIextract(''clear'')');
        uicontrol(fg,'Style','PushButton','Units','normalized','Position',[.86 .93 .12 .03],'ForegroundColor','r','FontWeight','bold',...
          	'String','EXIT','Callback','se_ROIextract(''exit'');','ToolTipString','quit');

        st.mode = uicontrol(fg,'Style','popupmenu','Units','normalized','Position',[0.215 .93 .3 .03],...
            'String',str2mat('Highest prob.','All assigned'),'Value',2,'FontWeight','bold');

        
         for i = 1:size(MAP,2)
            if i<=ceil(size(MAP,2)/3)
                st.area((index(i)*2)-1) = uicontrol(fg,'Style','checkbox','Units','normalized','Position',[0.01 .893-.02*i .02 .02],'Callback','se_createROI(''checked'')'); 
                st.area((index(i)*2)) =   uicontrol(fg,'Style','checkbox','Units','normalized','Position',[0.04 .893-.02*i .02 .02],'Callback','se_createROI(''checked'')');
                uicontrol(fg,'Style','text','Units','normalized','Position',[0.075 .89-.02*i .34 .025],'String',{MAP(index(i)).name},'HorizontalAlignment','left','FontSize',hFS-4);
            elseif i<=(ceil(size(MAP,2)/3)*2)
                st.area((index(i)*2)-1) = uicontrol(fg,'Style','checkbox','Units','normalized','Position',[0.31 .893-.02*(i-ceil(size(MAP,2)/3)) .02 .02],'Callback','se_createROI(''checked'')'); 
                st.area((index(i)*2)) =   uicontrol(fg,'Style','checkbox','Units','normalized','Position',[0.34 .893-.02*(i-ceil(size(MAP,2)/3)) .02 .02],'Callback','se_createROI(''checked'')');
                uicontrol(fg,'Style','text','Units','normalized','Position',[0.375 .89-.02*(i-ceil(size(MAP,2)/3)) .24 .025],'String',{MAP(index(i)).name},'HorizontalAlignment','left','FontSize',hFS-4);
            else
                st.area((index(i)*2)-1) = uicontrol(fg,'Style','checkbox','Units','normalized','Position',[0.61 .893-.02*(i-ceil(size(MAP,2)/3)*2) .02 .02],'Callback','se_createROI(''checked'')'); 
                st.area((index(i)*2)) =   uicontrol(fg,'Style','checkbox','Units','normalized','Position',[0.64 .893-.02*(i-ceil(size(MAP,2)/3)*2) .02 .02],'Callback','se_createROI(''checked'')');
                uicontrol(fg,'Style','text','Units','normalized','Position',[0.675 .89-.02*(i-ceil(size(MAP,2)/3)*2) .24 .025],'String',{MAP(index(i)).name},'HorizontalAlignment','left','FontSize',hFS-4);
            end
            ArResp(((index(i)*2)-1)).nvox = 0; ArResp(((index(i)*2))).nvox = 0;
        end
        
        uicontrol(fg,'Style','text','Units','normalized','Position',[0.01 .895 .029 .02],'FontWeight','bold','String','L','HorizontalAlignment','center');
        uicontrol(fg,'Style','text','Units','normalized','Position',[0.04 .895 .029 .02],'FontWeight','bold','String','R','HorizontalAlignment','center');
        
        uicontrol(fg,'Style','text','Units','normalized','Position',[0.305 .895 .029 .02],'FontWeight','bold','String','L','HorizontalAlignment','center');
        uicontrol(fg,'Style','text','Units','normalized','Position',[0.335 .895 .029 .02],'FontWeight','bold','String','R','HorizontalAlignment','center');

        uicontrol(fg,'Style','text','Units','normalized','Position',[0.605 .895 .029 .02],'FontWeight','bold','String','L','HorizontalAlignment','center');
        uicontrol(fg,'Style','text','Units','normalized','Position',[0.635 .895 .029 .02],'FontWeight','bold','String','R','HorizontalAlignment','center');


        





    case 'exit'
        se_figure('Clear','Graphics');
        clear all;
        Anatomy('select');
end
