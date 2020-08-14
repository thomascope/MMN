function se_anatomy(op,varargin)

global st;
global MAP;
global SPM;
global xSPM;
global CLUSTER;
global displayType;
global group
global defaults;

if strncmp(op,'init',4) & prod(size(op))>4;

    clear global MAP SPM xSPM CLUSTER group
    fg = se_figure('getWin','Graphics');
    delete(fg)
    [Finter,Fgraph,CmdLine] = spm('FnUIsetup','Display',0);
    SPMid = spm('FnBanner',mfilename,'2.15');
    spm_help('!ContextHelp',[mfilename,'.m']);

    % get the image's filename {P}
    %-----------------------------------------------------------------------
    MapName = spm_select(1,'mat',['Select Map'],[],spm('Dir','se_anatomy'),'MPM.',1);

    se_getMap('anat',MapName);
    if numel(dir([MapName(1:end-8) '.img']))>0
        P      = deblank([MapName(1:end-8) '.img']);
    else
        P      = deblank([MapName(1:end-8) '.nii']);
    end

    fg = se_figure('GetWin','Graphics');
    if isempty(fg), error('Can''t create graphics window'); end
    se_figure('Clear','Graphics');

    se_orthviews('Reset');

    switch spm('FnBanner')
        case 'SPM99'
            se_orthviews('Image', P, [0.0 0.22 1 .8]);     % Image appears
        otherwise
            se_orthviews('Image', P, [0.0 0.2 1 .8]);     % Image appears
    end


    if strcmp(op,'initO');
        st.prog = 'OV';
    elseif strcmp(op,'initG');
        st.prog = 'GR';
    elseif strcmp(op,'initP');
        st.prog = 'PX';
    elseif strcmp(op,'initA');
        st.prog = 'AR';
    end
    se_anatomy('init',fg);
    return;
end;

if ~strcmp(op,'init') & ~strcmp(op,'reset') & isempty(st.vols{1}), my_reset; warning('Lost all the image information'); return; end;

if strcmp(op,'shopos'),
    % The position of the crosshairs has been moved.
    %-----------------------------------------------------------------------
    if strcmp(st.prog,'PX')
        fg = se_figure('Findwin','Graphics');
        if xSPM.xVOL.orig == 2
            set(st.mp,'String',sprintf('%.1f %.1f %.1f',se_orthviews('pos')+[0;+4;-5]));
        else
            set(st.mp,'String',sprintf('%.1f %.1f %.1f',se_orthviews('pos')));
        end
        [loc, aVox] = hier_area;
        set(st.in,'String',loc);

        CLUSTER.XYZmm = se_orthviews('pos');
        CLUSTER.XYZp = se_orthviews('pos',1);

        [MLev, VLev, m_anz] = se_overlap(1);
        set(st.mAssign,'String',MLev);
    else
        if isfield(st,'mp'),
            fg = se_figure('Findwin','Graphics');
            if any(findobj(fg) == st.mp),
                set(st.mp,'String',sprintf('%.1f %.1f %.1f',se_orthviews('pos')));
                pos = se_orthviews('pos',1);
                set(st.vp,'String',sprintf('%.1f %.1f %.1f',pos));
                [loc, aVox] = hier_area;
                set(st.in,'String',loc);
                if prod(size(CLUSTER))>0
                    fDa = 0;
                    try ,fDa = CLUSTER(1).VOL(round(pos(1)),round(pos(2)),round(pos(3))); end
                    if fDa;
                        if CLUSTER(1).VOL(round(pos(1)),round(pos(2)),round(pos(3))) ~= get(st.cluster,'Value')-1;
                            set(st.cluster,'Value',1+CLUSTER(1).VOL(round(pos(1)),round(pos(2)),round(pos(3))));
                            clu_title;
                        end
                    else
                        try; delete(st.FX);             end
                        try; se_anatomy('clearRender'); end
                        set(st.titel,'String','');
                        if ~strcmp(st.prog,'AR')
                            set(st.vAssign,'String',''); set(st.mAssign,'String','');
                        end
                        set(st.cluster,'Value',1)
                    end
                end
            else,
                st.Callback = ';'; rmfield(st,{'mp','vp','in'});
            end;
        else,
            st.Callback = ';';
        end;
    end
    return;
end;

if strcmp(op,'clearRender')
    st.Fgraph = findobj(get(0,'Children'),'Flat','Tag','xSATB');
    if ~isempty(st.Fgraph)
        try
            set(0,'CurrentFigure',st.Fgraph);
            clf
        catch
            st.Fgraph = [];
        end
    end
    ax=axes('Parent',st.Fgraph,'units','normalized','Position',[0, 0, 1, 1],'Visible','off'); image(0,'Parent',ax); set(ax,'YTick',[],'XTick',[]);
end


if strcmp(op,'setposmm'),
    % Move the crosshairs to the specified position
    %-----------------------------------------------------------------------
    if isfield(st,'mp'),
        fg = se_figure('Findwin','Graphics');
        if any(findobj(fg) == st.mp),
            pos = sscanf(get(st.mp,'String'), '%g %g %g');
            if strcmp(st.prog,'PX') & xSPM.xVOL.orig == 2
                pos = pos+[0;-4;5];
                se_orthviews('Reposition',pos);
            else
                se_orthviews('Reposition',pos);
            end
        end;
    end;
    return;
end;

if strcmp(op,'setposclust'),
    % Move the crosshairs to the specified cluster
    %--------------------------------------------------------------------------------------------------------------------------------
    fg = se_figure('Findwin','Graphics');
    poscl = get(st.cluster,'Value')-1;
    if poscl > 0
        pos =  [round(mean(CLUSTER(poscl).XYZmm(1,:)));
            round(mean(CLUSTER(poscl).XYZmm(2,:)));
            round(mean(CLUSTER(poscl).XYZmm(3,:)))];
        set(st.cluster,'Value',poscl+1)
        clu_title;
        set(st.maxima,'Value',2)
        se_anatomy('setposmax');
    end;
    if get(st.zoomer,'Value')>1
        se_anatomy('zoom_in')
    end
    return;
end;

if strcmp(op,'setposmax'),
    % Move the crosshairs to the specified position
    %--------------------------------------------------------------------------------------------------------------------------------
    fg = se_figure('Findwin','Graphics');
    posmx = get(st.maxima,'Value')-1;
    poscl = get(st.cluster,'Value')-1;
    if posmx > 0
        [B,INDEX] = sort(CLUSTER(poscl).maxZ);
        posmx = size(B,2)+1-posmx;
        if CLUSTER(poscl).maxZ == 0
            pos   = xSPM.xVOL.M * [mean(CLUSTER(poscl).XYZ(1,:)');...
                mean(CLUSTER(poscl).XYZ(2,:)');...
                mean(CLUSTER(poscl).XYZ(3,:)'); 1];		%-Reorientation matrix
            pos = pos(1:3);
        else
            pos   = xSPM.xVOL.M * [CLUSTER(poscl).maxM(1,INDEX(posmx));...
                CLUSTER(poscl).maxM(2,INDEX(posmx));...
                CLUSTER(poscl).maxM(3,INDEX(posmx)); 1];		%-Reorientation matrix
            pos = pos(1:3);
        end
        if size(CLUSTER(poscl).XYZ,2) ~= 1
            se_orthviews('Reposition',pos);
        end
    end
end;

if strcmp(op,'setposarea'),
    poscl = get(st.cluster,'Value')-1;
    if prod(size(group))<2; se_getMap('groupStat'); end
    se_figure('Clear','Interactive');
    if isfield(st,'hAx');
        delete(st.hAx); st = rmfield(st,'hAx');
    end
    se_orthviews('rmblobs',1);
    PSC = []; GV = [MAP.GV];

    locNr = get(st.cluster,'Value')-1;
    check = spm_input(['Process area ' MAP(locNr).name '  ?'],'+1','y/n',[1,2],1);
    if check == 1
        side = MAP(1).orient*spm_input('Side','!+1','b','Left|Right',[-1 1],1,'Select side');
        mode = 1; % spm_input('Voxel','!+1','b','All assigned|Highest probability',[2 1],2,'Select mode');
        [nvox, QQ, targets] = se_defineVoxel(side,locNr,mode);

        XYZ = MAP(locNr).XYZ(:,QQ);
        Z = MAP(locNr).Z(QQ);

        se_orthviews('addcolouredblobs',1,XYZ,Z,MAP(1).MaxMap.mat,[1 0 0]);
        se_orthviews('Reposition',mean(targets'));

        spm_progress_bar('Init',prod(size(group)),'Subjects completed');
        for i=1:prod(size(group))
            psc(i) = se_PSC_area(group(i).xSPM,targets);
            PSC = [PSC mean(psc(i).PSC')'];
            spm_progress_bar('Set',i);
        end
        spm_progress_bar('Clear')

        nvox = psc(i).number;
        regressors = size(PSC,1)
        [PSCmin PSCmax Ay] = displayPSC(PSC,regressors,[0.51 0.17 0.46 0.335],1);
        text(PSCmin,Ay-10,sprintf('Area: %s overlap %s (%s voxel)',MAP(locNr).name,int2str(min(Z)/25),int2str(nvox)),'FontSize',12);
    else
        se_figure('Clear','Interactive');
        return
    end

end;

if strcmp(op,'setposvx'),
    % Move the crosshairs to the specified position
    %-----------------------------------------------------------------------
    if isfield(st,'mp'),
        fg = se_figure('Findwin','Graphics');
        if any(findobj(fg) == st.vp),
            pos = sscanf(get(st.vp,'String'), '%g %g %g');
            if length(pos)~=3,
                pos = se_orthviews('pos',1);
            end;
            tmp = st.vols{1}.premul*st.vols{1}.mat;
            pos = tmp(1:3,:)*[pos ; 1];
            se_orthviews('Reposition',pos);
        end;
    end;
    return;
end;


if strcmp(op,'addblobs'),
    if st.prog == 'AR'
    else
        % Add blobs to the image - in full colour
        se_figure('Clear','Interactive');
        if st.prog == 'GR'
            try
                group(1);
                spm('alert"',{'Please specify the secound level statistic'},'Next Step',0,1);
                se_getMap('SPM','');
            catch
                se_getMap('group','');
            end
        elseif st.prog == 'PX'
            spm('alert!','No activation map in quick-check mode',sqrt(-1));
        elseif st.prog == 'OV'
            OLwhat = spm_input('Overlay','!+1','b','SPM|Image',[1 2],1,'Select mode');
            if OLwhat == 1
                se_getMap('SPM','');
            else
                se_getMap('IMAGE','');
            end
        else
            se_getMap('SPM','');
        end

        if numel(SPM.XYZ)>0
            se_getMap('cluster','')
            
            if min(SPM.Z) == max(SPM.Z)
                c = spm_input('Colour',1,'m','Red blobs|Yellow blobs|Green blobs|Cyan blobs|Blue blobs|Magenta blobs',[1 2 3 4 5 6],1);
                colours = [1 0 0;1 1 0;0 1 0;0 1 1;0 0 1;1 0 1];
                se_orthviews('addcolouredblobs',1,SPM.XYZ,SPM.Z,xSPM.xVOL.M,colours(c,:));
            else
                c = spm_input('Colour',1,'m','Red -> White colormap|Red blobs|Yellow blobs|Green blobs|Cyan blobs|Blue blobs|Magenta blobs',[1 2 3 4 5 6 7],1);
                if c > 1
                    colours = [1 0 0;1 1 0;0 1 0;0 1 1;0 0 1;1 0 1];
                    se_orthviews('addcolouredblobs',1,SPM.XYZ,SPM.Z,xSPM.xVOL.M,colours(c-1,:));
                else
                    se_orthviews('AddBlobs',1,SPM.XYZ,SPM.Z-min(SPM.Z),xSPM.xVOL.M);
                end
            end
            set(st.blobber,'String','Remove Blobs','Callback','se_anatomy(''rmblobs'');');
            set(st.cluster,'String',str2mat(getCluster(1)));
            se_orthviews('Redraw');
        else
            spm('alert"',{'No superthreshold voxel !'},'',0,1);
        end
    end
end;

if strcmp(op,'rmblobs'),
    % Remove all blobs from the images
    se_orthviews('rmblobs',1);
    set(st.blobber,'String','Add Blobs','Callback','se_anatomy(''addblobs'');');
    set(st.cluster,'String','No activation cluster found','Value',1);
    set(st.titel,'String',''); 
    try; delete(st.FX);             end
    try; se_anatomy('clearRender'); end

    if ~(st.prog == 'AR')
        set(st.vAssign,'String','');
        set(st.mAssign,'String','');
    end
end;

if strcmp(op,'window'),
    op = get(st.win,'Value');
    if op == 1,
        se_orthviews('window',1);
    else,
        se_orthviews('window',1,spm_input('Range','+1','e','',2));
    end;
end;


if strcmp(op,'reset'),
    my_reset;
end;


if strcmp(op,'exit'),
    my_reset;
    Anatomy('select');
end;


if strcmp(op,'zoom_in'),
    op = get(st.zoomer,'Value');
    if op==1,
        se_orthviews('resolution',1);
        se_orthviews('MaxBB');
    else,
        vx = sqrt(sum(st.Space(1:3,1:3).^2));
        vx = vx.^(-1);
        pos = se_orthviews('pos');
        pos = st.Space\[pos ; 1];
        pos = pos(1:3)';
        if     op == 2, st.bb = [pos-80*vx ; pos+80*vx] ; se_orthviews('resolution',1);
        elseif op == 3, st.bb = [pos-40*vx ; pos+40*vx] ; se_orthviews('resolution',.5);
        elseif op == 4, st.bb = [pos-20*vx ; pos+20*vx] ; se_orthviews('resolution',.25);
        elseif op == 5, st.bb = [pos-10*vx ; pos+10*vx] ; se_orthviews('resolution',.125);
        else          , st.bb = [pos- 5*vx ; pos+ 5*vx] ; se_orthviews('resolution',.125);
        end;
    end;
    return;
end;

if strcmp(op,'Info'),
    st.mode = 1;
    if isfield(st,'hAx');
        delete(st.hAx);
        st = rmfield(st,'hAx');
    end
    set(st.titel,'Visible','on');
    set(st.mAssign,'Visible','on');
    set(st.vAssign,'Visible','on');
    set(st.clustBox,'Visible','on');
    set(st.maxBox,'Visible','on');
    set(st.cluster,'Visible','on');
    set(st.maxima,'Visible','on');
    se_figure('Clear','Interactive');
end


if strcmp(op,'SVC'),

    if strcmp(spm('FnBanner'),'SPM99')
        h = spm('alert*','Only available for SPM2 onwards');
    else
        if strcmp(get(st.blobber,'String'),'Remove Blobs') & isfield(SPM,'STAT')

            cMode = spm_input('ROI definition','!+1','b','(anat.) Mask | Sphere',[1 2],1,'Select ROI definition');
            
            if cMode == 1
                Msk   = spm_select(1,'image','Image defining search volume',[],pwd,'ROI');
                D     = spm_vol(Msk);
                str   = sprintf('image mask: %s',spm_str_manip(Msk,'a30'));
                VOX   = sqrt(sum(D.mat(1:3,1:3).^2));
                FWHM  = SPM.FWHM.*(SPM.VOX./VOX);
                XYZm   = inv(D.mat) * (xSPM.xVOL.M*[xSPM.xVOL.XYZ; ones(1, xSPM.xVOL.S)]);
                k     = find(spm_sample_vol(D, XYZm(1,:), XYZm(2,:), XYZm(3,:),0) > 0);

                SPM.XYZ = xSPM.xVOL.XYZ(:,k);

                SPM.S     = length(k);
                SPM.R     = spm_resels(FWHM,D,'I');

            else
                XYZmm    = (xSPM.xVOL.M*[xSPM.xVOL.XYZ; ones(1, xSPM.xVOL.S)]);
                FWHM     = SPM.FWHM;
                
                space    = spm_input('Centre in','!+1','b','Anat. MNI | MNI',[1 2],2,'Select space');
                xyzmm    = spm_input('Centre of VOI {mm}',-2);
                if space == 2; xyzmm = xyzmm + [0 -4 5]; end
                D          = spm_input('radius of VOI {mm}',-2);
                
                str        = sprintf('%0.1fmm sphere',D);
                k          = find(sum((     XYZmm(1:3,:) - xyzmm' * ones(1,size(XYZmm,2))).^2) <= D^2);
                D          = D./SPM.VOX;

                SPM.XYZ = xSPM.xVOL.XYZ(:,k);

                SPM.S     = length(k);
                SPM.R     = spm_resels(FWHM,D,'S');
           
            end


            SPM.Z         = Inf;
            for i     = SPM.Ic
                SPM.Z = min(SPM.Z,spm_get_data(xSPM.xCon(i).Vspm,SPM.XYZ));
            end

            switch SPM.STAT
                case 'T'
                    SPM.Ps = (1 - spm_Tcdf(SPM.Z,SPM.df(2))).^SPM.n;
                case 'P'
                    SPM.Ps = (1 - SPM.Z).^SPM.n;
                case 'F'
                    SPM.Ps = (1 - spm_Fcdf(SPM.Z,SPM.df)).^SPM.n;
            end

            SPM.u      = -Inf;  SPM.k      = 0;

            if SPM.STAT ~= 'P'
                str = 'FWE|none';
                switch spm_input('p value adjustment to control','+1','b',str,[],1)
                    case 'FWE' % family-wise false positive rate
                        u  = spm_input('p value (family-wise error)','+0','r',0.05,1,[0,1]);
                        SPM.u  = spm_uc(u,SPM.df,SPM.STAT,SPM.R,SPM.n,SPM.S);
                    otherwise  %-NB: no adjustment
                        SPM.u  = spm_input(['threshold {',SPM.STAT,' or p value}'],'+0','r',0.001,1);
                        if SPM.u <= 1; SPM.u = spm_u(SPM.u^(1/SPM.n),SPM.df,SPM.STAT); end
                end
            elseif STAT == 'P'
                SPM.u  = spm_input(['p value threshold for PPM'],'+0','r',.95,1);
            end % (if STAT)

            Q      = find(SPM.Z > SPM.u);
            SPM.Z      = SPM.Z(:,Q);
            SPM.XYZ    = SPM.XYZ(:,Q);
            if isempty(Q)
                warning(sprintf('No voxels survive height threshold u=%0.2g',u))
            end

            if ~isempty(SPM.XYZ)
                SPM.k     = spm_input('& extent threshold {voxels}','+1','r',0,1,[0,Inf]);
                A     = spm_clusters(SPM.XYZ); Q     = [];
                for i = 1:max(A); j = find(A == i); if length(j) >= SPM.k; Q = [Q j]; end; end

                SPM.Z     = SPM.Z(:,Q);
                SPM.XYZ   = SPM.XYZ(:,Q);
                if isempty(Q); warning(sprintf('No voxels survive extent threshold k=%0.2g',k));  end
            end
            SPM.XYZmm = xSPM.xVOL.M * [SPM.XYZ; ones(1,size(SPM.XYZ,2))]; SPM.XYZmm = SPM.XYZmm(1:3,:);
            SPM.XYZp = [round(SPM.XYZmm(1,:)-MAP(1).MaxMap.mat(1,4)); round(SPM.XYZmm(2,:)-MAP(1).MaxMap.mat(2,4)); round(SPM.XYZmm(3,:)-MAP(1).MaxMap.mat(3,4))];

            if numel(SPM.XYZ)>0
                se_getMap('cluster','')
                c = spm_input('Colour',1,'m','Red blobs|Yellow blobs|Green blobs|Cyan blobs|Blue blobs|Magenta blobs',[1 2 3 4 5 6],1);
                colours = [1 0 0;1 1 0;0 1 0;0 1 1;0 0 1;1 0 1];

                se_orthviews('rmblobs',1);
                se_orthviews('addcolouredblobs',1,SPM.XYZ,SPM.Z,xSPM.xVOL.M,colours(c,:));
                set(st.cluster,'String',str2mat(getCluster(1)));
            else
                spm('alert"',{'No superthreshold voxel !'},'',0,1);
            end

        else
            h = spm('alert*','No activation map loaded','SVC not possible');
        end
    end
end

if strcmp(op,'Plot'),
%     se_figure('clear','Interactive')
%     if strcmp(get(st.blobber,'String'),'Remove Blobs') & isfield(SPM,'STAT')
%         try; delete(st.FX);             end
%         try; se_anatomy('clearRender'); end
%         xyzp = round(se_orthviews('pos',1));
%         fDa = 0;, try ,fDa = CLUSTER(1).VOL((xyzp(1)),(xyzp(2)),(xyzp(3))); end
%         if 0; fDa
%             definition = spm_input('Process whole cluster ?','+1','y/n',[1,2],2);
%             if definition == 1
%                 xyz =  xSPM.xVOL.M*[CLUSTER(fDa).XYZ;ones(1,size(CLUSTER(fDa).XYZ,2))];
%             else
%                 xyz = [se_orthviews('pos'); 1];
%             end
%         else
%             definition = 2;
%             xyz = [se_orthviews('pos'); 1];
%         end
%         se_graph(SPM,xSPM,xyz,definition)
%     else
%         h = spm('alert*','No activation map loaded');
%     end
    se_figure('clear','Interactive')
    if strcmp(get(st.blobber,'String'),'Remove Blobs') & isfield(SPM,'STAT')
         try; delete(st.FX);             end
        xyz = [se_orthviews('pos'); 1];
        se_graph(SPM,xSPM,xyz);
    else
        h = spm('alert*','No activation map loaded');
    end
end



if strcmp(op,'Tab'),
    if strcmp(get(st.blobber,'String'),'Remove Blobs') & strcmp(st.prog,'OV')

        try
            thres = strrep(SPM.thresDef,'0.','');, thres = strrep(thres,'(','-'); thres = strrep(thres,')',''); thres = strrep(thres,'p<','_p'); thres = strrep(thres,'P>','_P');
        catch
            thres = '';
        end

        filename = [strrep(SPM.title,' ','_') strrep(thres,' ','') '.txt'];
        filename = strrep(filename,'>','vs');    filename = strrep(filename,'<','vs');
        filename = strrep(filename,')','-');     filename = strrep(filename,'(','-');
        filename = strrep(filename,'/','-');     filename = strrep(filename,'\','-');
        filename = fullfile(pwd,filename);
        if strcmp(computer,'PCWIN')
            fid = fopen(filename,'wt');
        else
            fid = fopen(filename,'w+');
        end

        for poscl = 1:size(CLUSTER,2)
            se_tabelate(poscl,fid)
            fprintf(fid,repmat('\n',1,5));
        end

        status = fclose(fid);
        spm('alert"',{'Results saved to file'; spm_str_manip(filename,'t'); 'in path'; pwd});

    else
        h = spm('alert*','No activation map loaded');
    end
end


if strcmp(op,'PSCG'),
    try
        group(1);
    catch
        se_anatomy('addblobs')
    end
    PSC = [];
    xyzp = round(se_orthviews('pos',1));

    fDa = 0;, try ,fDa = CLUSTER(1).VOL((xyzp(1)),(xyzp(2)),(xyzp(3))); end
    %     if fDa;
    %         PSCmode = spm_input('Which mode ?', '!+1', 'm','Single voxel|Cluster by area|Cluster',[1 2 4],1);
    %     else
    PSCmode = 1;
    %     end

    switch PSCmode
        case 1 % Single voxel
            targets  = round(se_orthviews('pos'));

        case 2  % Cluster by area
            cl  = (get(st.cluster,'Value')-1); GV = [MAP.GV];
            tmp   = CLUSTER(cl).XYZp; tms = sign(CLUSTER(cl).XYZmm); clear gtmp;

            for xbx = 0:xSPM.xVOL.VOX(1)-1, for ybx = 0:xSPM.xVOL.VOX(2)-1, for zbx = 0:xSPM.xVOL.VOX(3)-1
                        if exist('gtmp')
                            gtmp = [gtmp(1,:) tmp(1,:)-(tms(1,:)*xbx);...
                                gtmp(2,:) tmp(2,:)-(tms(2,:)*ybx);...
                                gtmp(3,:) tmp(3,:)-(tms(3,:)*zbx)];
                        else
                            gtmp = [tmp(1,:)-(tms(1,:)*xbx);...
                                tmp(2,:)-(tms(2,:)*ybx);...
                                tmp(3,:)-(tms(3,:)*zbx)];
                        end
                    end, end, end

            Msk = spm_sample_vol(MAP(1).AnatMask,gtmp(1,:),gtmp(2,:),gtmp(3,:),0);
            Volume = spm_sample_vol(MAP(1).MaxMap,gtmp(1,:),gtmp(2,:),gtmp(3,:),0); % Corresponding voxel + next towards origin in all directions

            VLev = {}; vp = [];
            for v=1:size(MAP,2)
                if any(Volume(Msk == 2) == GV(v)), vp = [vp; size(find(Volume(Msk == 2) == GV(v)),2) v 2]; end
                if any(Volume(Msk == 1) == GV(v)), vp = [vp; size(find(Volume(Msk == 1) == GV(v)),2) v 1]; end
            end, vp = sortrows(vp);

            if ~any(vp); spm('alert!','No areas found at the specified positions',sqrt(-1)); return; end
            selec = []; sel = []; vp = sortrows(vp); lr = [];
            for v=(size(vp,1)):-1:1
                if (100*vp(v,1)/(CLUSTER(cl).voxel(1)*prod(xSPM.xVOL.VOX))) >= 5
                    if  vp(v,3) == -1;
                        selec = [selec '|' ['left ' MAP(vp(v,2)).name] ': ' int2str((100*vp(v,1)/(CLUSTER(cl).voxel(1)*prod(xSPM.xVOL.VOX)))) '%']; sel = [sel vp(v,2)]; lr = [lr vp(v,3)];
                    else
                        selec = [selec '|' ['right ' MAP(vp(v,2)).name] ': ' int2str((100*vp(v,1)/(CLUSTER(cl).voxel(1)*prod(xSPM.xVOL.VOX)))) '%']; sel = [sel vp(v,2)]; lr = [lr vp(v,3)];
                    end
                end
            end

            choice = spm_input('Select area ', '!+1', 'm',selec,[1:size(sel,2)],1); locNr = sel(choice);
            pxyz = round(se_orthviews('pos')); onXYZ = find(SPM.XYZmm(1,:) == pxyz(1) & SPM.XYZmm(2,:) == pxyz(2) & SPM.XYZmm(3,:) == pxyz(3));
            A     = spm_clusters(SPM.XYZ); Q     = find(A == A(onXYZ));
            xx = round(SPM.XYZmm(1,Q)-MAP(1).MaxMap.mat(1,4)); yy = round(SPM.XYZmm(2,Q)-MAP(1).MaxMap.mat(2,4)); zz = round(SPM.XYZmm(3,Q)-MAP(1).MaxMap.mat(3,4));
            Q = Q(spm_sample_vol(MAP(1).MaxMap,xx,yy,zz,0) == MAP(locNr).GV & ...
                spm_sample_vol(MAP(1).AnatMask,xx,yy,zz,0) == lr(choice));
        case 4
            A     = spm_clusters(SPM.XYZ);
            pxyz = round(se_orthviews('pos')); onXYZ = find(SPM.XYZmm(1,:) == pxyz(1) & SPM.XYZmm(2,:) == pxyz(2) & SPM.XYZmm(3,:) == pxyz(3));
            Q     = find(A == A(onXYZ));
    end

    if PSCmode > 1
        targets = SPM.XYZmm(:,Q);
    end

    spm_progress_bar('Init',prod(size(group)),'Subjects completed');
    for i=1:prod(size(group))
        psc(i) = se_PSC_area(group(i).xSPM,targets);
%         if size(targets,2)>1
%             PSC = [PSC mean(psc(i).PSC')'];
%         else
            PSC = [PSC psc(i).PSC];
%         end
        spm_progress_bar('Set',i);
    end
    spm_progress_bar('Clear')

    if psc(i).Sessions == 1
        descrip = ['Rows: Conditions;   Columns: Sub1, Sub2, Sub3... Sub' int2str(prod(size(group)))];
    else
        descrip = ['Rows: Conditions;   Columns: Sub1 Ses1 - Ses' int2str(psc(i).Sessions ) ',  Sub2 Ses1 - Ses' int2str(psc(i).Sessions ) '... Sub' int2str(prod(size(group))) ' Ses1 - Ses' int2str(psc(i).Sessions )]; 
    end
    
    save PSC.mat PSC descrip;

    
    if PSCmode > 1
        nvox = size(targets,2);
    else
        nvox = 1;
    end

    regressors = size(PSC,1);
    fg = se_figure('GetWin','Graphics'); WS = spm('WinScale');
    set(st.mAssign,'Visible','off'); set(st.vAssign,'Visible','off');
    set(st.clustBox,'Visible','off'); set(st.maxBox,'Visible','off'); set(st.cluster,'Visible','off'); set(st.maxima,'Visible','off');


    st.mode = 3;
    se_figure('Clear','Interactive');
    if isfield(st,'hAx');
        delete(st.hAx);
        st = rmfield(st,'hAx');
    end

    if any(~isnan(sum(PSC)))
    [PSCmin, PSCmax, Ay] = displayPSC(PSC(:,~isnan(sum(PSC))),regressors,[0.51 0.21 0.46 0.335],0);
    switch PSCmode
        case 1
            text(PSCmin,Ay,sprintf('Position: x=%1.0f, y=%1.0f, z=%1.0f,',se_orthviews('pos')),'FontSize',12);
            [Area, code] = hier_area;
            text(PSCmin,Ay-17*WS(2),sprintf('Area: %s',Area),'FontSize',12);
        case 2
            text(PSCmin,Ay,['Activated voxels in area ' MAP(locNr).name],'FontSize',12);
        case 4
            text(PSCmin,Ay,['All ' int2str(CLUSTER(CLUSTER(1).VOL(round(xyzp(1)),round(xyzp(2)),round(xyzp(3)))).voxel(1)) ' voxels in cluster ' sprintf('%2.0f', get(st.cluster,'Value')-1)],'FontSize',12);
    end

    else
        spm('alert"',{'No data  available !'},'',0,1);
        se_anatomy('Info')    
    end
end




% ----------------------- INITIALIZE ----------------------------

if strcmp(op,'init'),
    st.pwd = pwd;
    st.SPM = spm('FnBanner');
    fg = se_figure('GetWin','Graphics');
    se_orthviews('Interp',0)
    if isempty(st.vols{1}), return; end;

    se_orthviews('MaxBB');
    st.callback = 'se_anatomy(''shopos'');';

    st.B = [0 0 0  0 0 0  1 1 1  0 0 0];

    WS = spm('WinScale');

    xp = -20;
    switch st.SPM
        case 'SPM99'
        otherwise
            set(gcf,'DefaultUicontrolFontSize',spm('FontSizes',8));
            xp = -30;
    end

    % Crosshair position
    %-----------------------------------------------------------------------
    %	uicontrol(fg,'Style','Frame','Position',[45 110+xp 190 100].*WS,'DeleteFcn','se_anatomy(''reset'');');
    uicontrol(fg,'Style','Frame','Position',[50+10 105-31 180 82].*WS);

    uicontrol(fg,'Style','Text', 'Position',[55+10 166-29 170 016].*WS,'String','Crosshair Position');
    uicontrol(fg,'Style','PushButton', 'Position',[55+10 161-29 170 006].*WS,...
        'Callback','se_orthviews(''Reposition'',[0 0 0]);','ToolTipString','move crosshairs to origin');
    uicontrol(fg,'Style','Text', 'Position',[55+10 135-29 35 020].*WS,'String','mm:');
    uicontrol(fg,'Style','Text', 'Position',[55+10 115-29 35 020].*WS,'String','vx:');
    uicontrol(fg,'Style','Text', 'Position',[55+10 106-29 65 014].*WS,'String','Area:');

    st.mp = uicontrol(fg,'Style','edit', 'Position',[90+10 140-25 135 020].*WS,'String','','Callback','se_anatomy(''setposmm'')','ToolTipString','move crosshairs to mm coordinates');
    st.vp = uicontrol(fg,'Style','edit', 'Position',[90+10 120-25 135 020].*WS,'String','','Callback','se_anatomy(''setposvx'')','ToolTipString','move crosshairs to voxel coordinates');
    st.in = uicontrol(fg,'Style','Text', 'Position',[120+10 106-25  85 014].*WS,'String','');

    % Assorted other buttons.
    %-----------------------------------------------------------------------
    uicontrol(fg,'Style','Frame','Position',[5 30+xp 270 70].*WS);
    st.zoomer = uicontrol(fg,'Style','popupmenu' ,'Position',[10 75+xp 125 20].*WS,...
        'String',str2mat('Full Volume','160x160x160mm','80x80x80mm','40x40x40mm','20x20x20mm','10x10x10mm'),...
        'Callback','se_anatomy(''zoom_in'')','ToolTipString','zoom in by different amounts');
    c = 'if get(gco,''Value'')==1, se_orthviews(''Space''), else, se_orthviews(''Space'', 1);end;se_anatomy(''zoom_in'')';
    uicontrol(fg,'Style','popupmenu' ,'Position',[10 55+xp 125 20].*WS,...
        'String',str2mat('World Space','Voxel Space'),...
        'Callback',c,'ToolTipString','display in aquired/world orientation');
    c = 'if get(gco,''Value'')==1, se_orthviews(''Xhairs'',''off''), else, se_orthviews(''Xhairs'',''on''); end;';
    uicontrol(fg,'Style','togglebutton','Position',[145 75+xp 125 20].*WS,...
        'String','Hide Crosshairs','Callback',c,'ToolTipString','show/hide crosshairs');
    uicontrol(fg,'Style','popupmenu' ,'Position',[145 55+xp 125 20].*WS,...
        'String',str2mat('NN interp','bilin interp','sinc interp'),...
        'Callback','tmp_ = [0 1 -4];se_orthviews(''Interp'',tmp_(get(gco,''Value'')))',...
        'Value',1,'ToolTipString','interpolation method for displaying images');
    uicontrol(fg,'Style','PushButton','Position',[10 35+xp 125 20].*WS,'ForegroundColor','r','FontWeight','bold',...
        'String','EXIT','Callback','se_anatomy(''exit'');','ToolTipString','quit');

    if ~strcmp(st.prog,'AR')
        st.blobber = uicontrol(fg,'Style','pushbutton','Position',[145 35+xp 125 20].*WS,...
            'String','Add Blobs','Callback','se_anatomy(''addblobs'');','ToolTipString','superimpose activations');
    end


    FS  = spm('FontSizes');
    hFS = FS(8);

    if (~strcmp(st.prog,'PX') & ~strcmp(st.prog,'AR'))
        if ~isfield(st,'mAssignNr');
            st.mAssignNr = 0;
        end
        if st.mAssignNr
            Smax = 'No maximum selected'; for cl=1:st.mAssignNr; Smax = [Smax '|' 'Maximum (' sprintf('%2.0f', cl) ')']; end
        else
            Smax = 'No maxima found';
        end

        st.maxBox = uicontrol(fg,'Style','Frame','Position',[285 322 197 35].*WS);
        st.maxima = uicontrol(fg,'Style','popupmenu','Position',[290 332 190 20].*WS,...
            'String',str2mat(Smax),'Callback','se_anatomy(''setposmax'')','ToolTipString','move crosshairs to maximum','FontSize',hFS+1);


        st.clustBox = uicontrol(fg,'Style','Frame','Position',[285 459 260 37].*WS);
        st.cluster = uicontrol(fg,'Style','popupmenu','Position',[290 469 250 20].*WS,...
            'String',str2mat(getCluster(2)),'Callback','se_anatomy(''setposclust'')','ToolTipString','move crosshairs to cluster','FontSize',hFS+1);
        st.titel = uicontrol(fg,'Style','text','Position',[290 435 308 20].*WS,'BackgroundColor',[1 1 1],'FontWeight','bold','HorizontalAlignment','left',...
            'String','','FontSize',hFS+1);

        st.vAssign = uicontrol(fg,'Style','text','Position',[290 358 300 80].*WS,'BackgroundColor',[1 1 1],'HorizontalAlignment','left',...
            'String',{},'FontSize',hFS);

        if strcmp(st.prog,'OV')
            st.mAssign = uicontrol(fg,'Style','text','Position',[290 10 300 305].*WS,'BackgroundColor',[1 1 1],'HorizontalAlignment','left',...
                'String',{},'FontSize',FS(7));
        else
            st.mAssign = uicontrol(fg,'Style','text','Position',[290 20+xp 300 300-xp].*WS,'BackgroundColor',[1 1 1],'HorizontalAlignment','left',...
                'String',{},'FontSize',hFS);
        end
    end

    if strcmp(st.prog,'AR');
        displayType = 'AR';
        st.blobber = uicontrol(fg,'Style','pushbutton','Position',[145+20 35+xp 125 20].*WS,...
            'String','Read statistics','Callback','se_getMap(''groupStat'')','ToolTipString','load statistical data');
        clust = 'No area selected';
        for cl=1:size(MAP,2); clust = [clust '|' MAP(cl).name]; end
        st.clustBox = uicontrol(fg,'Style','Frame','Position',[285 490 260 35].*WS);
        st.cluster = uicontrol(fg,'Style','popupmenu','Position',[290 500 250 20].*WS,...
            'String',str2mat(clust),'Callback','se_anatomy(''setposarea'')','ToolTipString','move crosshairs to area','FontSize',hFS+1);
        st.titel = uicontrol(fg,'Style','text','Position',[290 460 300 20].*WS,'BackgroundColor',[1 1 1],'FontWeight','bold','HorizontalAlignment','left',...
            'String','','FontSize',hFS+1);

    elseif strcmp(st.prog,'OV');
        displayType = 'OL';
        st.SVC = uicontrol(fg,'Style','PushButton','Position',[300 35+xp 75 20].*WS,'FontWeight','bold',...
            'String','S. V. C.','Callback','se_anatomy(''SVC'');','ToolTipString','Small volume correction');
        st.Tab = uicontrol(fg,'Style','PushButton','Position',[400 35+xp 75 20].*WS,'FontWeight','bold',...
            'String','Tab','Callback','se_anatomy(''Tab'');','ToolTipString','Tabelated output to file');
        st.Plot = uicontrol(fg,'Style','PushButton','Position',[500 35+xp 75 20].*WS,'FontWeight','bold',...
            'String','Plot','Callback','se_anatomy(''Plot'');','ToolTipString','Plot contrast estimates');

    elseif strcmp(st.prog,'GR')
        displayType = 'GA';
        st.mode = 1;
        st.Info = uicontrol(fg,'Style','PushButton','Position',[300 35+xp 80 20].*WS,'FontWeight','bold',...
            'String','Info','Callback','se_anatomy(''Info'');','ToolTipString','Localization of volumes and maxima');
        if strcmp(st.prog,'GR')
            st.PSC = uicontrol(fg,'Style','PushButton','Position',[480 35+xp 100 20].*WS,'FontWeight','bold',...
                'String','% signal change','Callback','se_anatomy(''PSCG'');','ToolTipString','Calculate % signal change');
        end

    elseif strcmp(st.prog,'PX')
        displayType = 'PX';
        hFS = FS(9);
        st.mAssign = uicontrol(fg,'Style','text','Position',[290 390 300 135].*WS,'BackgroundColor',[1 1 1],'HorizontalAlignment','left',...
            'String',{},'FontSize',hFS);
        if ~exist(xSPM)
            xSPM = struct('xVOL',{});
        end
        xSPM(1).xVOL.M = MAP(1).MaxMap.mat;
        xF = spm_figure('GetWin','Graphics'); xI = spm_figure('GetWin','Interactive');
        set(xF,'visible','off'), figure(xI)
        xSPM(1).xVOL.orig = spm_input('Coordinates','!+1','b','Anatomical|MNI',[1 2],2,'Select mode'); set(xF,'visible','on'),
    end

end;
return;




function my_reset
st.FX = findobj(get(0,'Children'),'Flat','Tag','SATB');
st.Fgraph = findobj(get(0,'Children'),'Flat','Tag','xSATB');

try; delete(st.FX); delete(st.Fgraph); end
se_figure('Clear','Graphics');
return;

function [PSCmin, PSCmax, Ay] = displayPSC(PSC,regressors,position,dotext)
global st
global group

FS        = spm('FontSizes');			%-Scaled font sizes
WS = spm('WinScale');
set(st.titel,'Visible','off'); fg = se_figure('GetWin','Graphics');
st.hAx   = axes('Position',position,'DefaultTextInterpreter','Tex',...
    'DefaultTextVerticalAlignment','Baseline','Units','points','Visible','on');
AxPos = get(st.hAx,'Position'); Ay     = floor(AxPos(4))-10;
sem = sqrt(var(PSC')/size(PSC,2)); PSCmin = min([0 min(mean(PSC')-sem*1.2)]); PSCmax = max([0 max(mean(PSC')+sem*1.2)]);
boxMax = Ay-25-regressors*35;
set(st.hAx,'YLim',[boxMax,AxPos(4)-28]); set(st.hAx,'XLim',[PSCmin,PSCmax]);
set(st.hAx,'YTick',[]); set(st.hAx,'YTickLabel',''); pp = 0;
line([PSCmin PSCmin],[Ay-19 boxMax],'Color','k');
line([PSCmax PSCmax],[Ay-19 boxMax],'Color','k');
line([PSCmin PSCmax],[boxMax boxMax],'Color','k');
line([PSCmin PSCmax],[Ay-19 Ay-19],'Color','k');
line([0 0],[Ay-19 boxMax],'Color','k','LineWidth',1);
lab = {'% signal change','(Mean and SEM indicated)'};
off = (max(max(PSCmax))-min(min(PSCmin)))*0.01;
for k = 1:regressors,
    line([0 mean(PSC(k,:))],[Ay-40-pp*35 Ay-40-pp*35],'Color',[0 0 0],'LineWidth',25)
    if mean(PSC(k,:))>0
        line([off/2 mean(PSC(k,:))-off],[Ay-40-pp*35 Ay-40-pp*35],'Color',[0.5 0.5 0.5],'LineWidth',22)
    else
        line([-off/2 mean(PSC(k,:))+off],[Ay-40-pp*35 Ay-40-pp*35],'Color',[0.5 0.5 0.5],'LineWidth',22)
    end
    line([mean(PSC(k,:))-sem(k) mean(PSC(k,:))+sem(k)],[Ay-40-pp*35 Ay-40-pp*35],'Color',[0 0 0],'LineWidth',1)
    line([mean(PSC(k,:))-sem(k) mean(PSC(k,:))-sem(k)],[Ay-45-pp*35 Ay-35-pp*35],'Color',[0 0 0],'LineWidth',1)
    line([mean(PSC(k,:))+sem(k) mean(PSC(k,:))+sem(k)],[Ay-45-pp*35 Ay-35-pp*35],'Color',[0 0 0],'LineWidth',1)
    text(0,Ay-40-pp*35,int2str(k),'Color','r','FontWeight','bold','FontSize',17,...
        'HorizontalAlignment','Center','VerticalAlignment','Middle')
    pp = pp+1;
end
xlabel(lab);
ys = get(st.hAx,'YLim');
Ay2 = boxMax-(ys(2)-ys(1))*.27;
for i=1:ceil(regressors/3)
    text(PSCmin,Ay2-((i-1)*(10+regressors)),[int2str(i*3-2) ': ' getCondLabel(i*3-2) ' (' num2str(mean(PSC(i*3-2,:)),'%5.3f') ')'],'VerticalAlignment','middle','FontSize',FS(10));
    if (i-1)*3+2 <= regressors
        text(PSCmin+((PSCmax-PSCmin)/3)*1,Ay2-((i-1)*(10+regressors)),[int2str(i*3-1) ': ' getCondLabel(i*3-1) ' (' num2str(mean(PSC(i*3-1,:)),'%5.3f') ')'],'VerticalAlignment','middle','FontSize',FS(10))
    end
    if (i-1)*3+3 <= regressors
        text(PSCmin+((PSCmax-PSCmin)/3)*2,Ay2-((i-1)*(10+regressors)),[int2str(i*3) ': ' getCondLabel(i*3) ' (' num2str(mean(PSC(i*3,:)),'%5.3f') ')'],'VerticalAlignment','middle','FontSize',FS(10))
    end
end


function string = getCondLabel(i)
global st
global group
switch st.SPM
    case 'SPM99'; strin = group(1).xSPM.Sess{1}.name{i};
    otherwise, strin = group(1).xSPM.Sess(1).Fc(i).name;
end
string = strrep(strin,'0',''); string = strrep(string,'1',''); string = strrep(string,'_',' '); string = strrep(string,'  ',' ');

function [Area, code] = hier_area
global MAP
global st
xyzv = se_orthviews('pos',1); locNr = 0;
aVox = find([MAP.GV] == spm_sample_vol(st.vols{1},xyzv(1),xyzv(2),xyzv(3),st.hld));
if aVox; Area = MAP(aVox).name; code = aVox;
else;
    ML = round(spm_sample_vol(MAP(1).Macro,xyzv(1),xyzv(2),xyzv(3),0)');
    code = 0;
    if ML > 0;
        Area = MAP(1).MLabels.Labels{ML};
    else
        Area = 'Unknown area';
    end
end


function clu_title
global st
global SPM
global CLUSTER


if isfield(SPM,'STAT')
    try
        set(st.titel,'String',['Cluster ' int2str(get(st.cluster,'Value')-1) ' (' int2str(CLUSTER(get(st.cluster,'Value')-1).voxel(1)) ' vox)'...
            ': ' SPM.title ' '  SPM.thresDef]);
    catch
        set(st.titel,'String',['Cluster ' int2str(get(st.cluster,'Value')-1) ' (' int2str(CLUSTER(get(st.cluster,'Value')-1).voxel(1)) ' vox)'...
            ': ' SPM.title '  ( ' SPM.STAT '> ' sprintf('%4.2f', SPM.u) ')']);
    end
else
    titl = SPM.title;
    if length(titl)>15; titl = titl(1:15); end
    set(st.titel,'String',['Cluster ' int2str(get(st.cluster,'Value')-1) ' (' int2str(CLUSTER(get(st.cluster,'Value')-1).voxel(1)) ' vox)'...
        ': ' titl ]); % '  (u > ' sprintf('%4.2f', SPM.u) ', k > ' int2str(SPM.k)  ')'
end
[MLev, VLev, m_anz] = se_overlap(get(st.cluster,'Value')-1);
set(st.vAssign,'String',VLev); set(st.mAssign,'String',MLev); st.mAssignNr = m_anz;
Smax = 'No maximum selected';
for cl=1:st.mAssignNr;
    Smax = [Smax '|' 'Maximum (' sprintf('%2.0f', cl) ')'];
end
set(st.maxima,'String',Smax,'Value',1);


function clust = getCluster(mode)
global st;
global MAP;
global SPM;
global xSPM;
global CLUSTER;
global displayType;
global group
global defaults;

if mode == 1
    clust = 'No cluster selected';

    for cl=1:size(CLUSTER,2)
        clust = [clust '|Cluster ' sprintf('%2.0f', cl) ':  x= ' sprintf('%+3.0f', mean(CLUSTER(cl).XYZmm(1,:)))...
            '  y= ' sprintf('%+2.0f', mean(CLUSTER(cl).XYZmm(2,:))) '  z= ' sprintf('%+2.0f', mean(CLUSTER(cl).XYZmm(3,:)))];
        if xSPM.xVOL.orig == 2
            clust = [clust ' (MNI: ' sprintf('%+3.0f', mean(CLUSTER(cl).XYZmm(1,:))) '/'...
                sprintf('%+3.0f', mean(CLUSTER(cl).XYZmm(2,:))+4) '/' sprintf('%+3.0f', mean(CLUSTER(cl).XYZmm(3,:))-4) ')'];
        end
    end
else
    clear clust
    if prod(size(CLUSTER))>0
        clust = 'No cluster selected';
        for cl=1:size(CLUSTER,2)
            clust = [clust '|Cluster ' sprintf('%2.0f', cl) ':  x= ' sprintf('%+3.0f', mean(CLUSTER(cl).XYZmm(1,:)))...
                '  y= ' sprintf('%+2.0f', mean(CLUSTER(cl).XYZmm(2,:))) '  z= ' sprintf('%+2.0f', mean(CLUSTER(cl).XYZmm(3,:)))];
            if xSPM.xVOL.orig == 2
                clust = [clust ' (orig. MNI: ' sprintf('%+3.0f', mean(CLUSTER(cl).XYZmm(1,:))) '/'...
                    sprintf('%+3.0f', mean(CLUSTER(cl).XYZmm(2,:))+4) '/' sprintf('%+3.0f', mean(CLUSTER(cl).XYZmm(3,:))-4) ')'];
            end
        end
    end
    if ~exist('clust')
        clust = 'No activation cluster found';
    end
end

