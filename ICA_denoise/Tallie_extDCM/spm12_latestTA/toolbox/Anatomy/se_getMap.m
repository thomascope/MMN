function getMap(welche,MapName)
global st;
global MAP;
global SPM;
global xSPM;
global CLUSTER;
global group;
global defaults

xtr = 0;
ytr = -4;
ztr = +5;
switch welche
    case 'anat'
        load(MapName);
        flip = -1;
        try, flip = defaults.analyze.flip; defaults.analyze.flip = 0; end
        if numel(dir([MapName(1:end-8) '.nii']))>0
            MAP(1).MaxMap = spm_vol([MapName(1:end-8) '.nii']);
        else
            MAP(1).MaxMap = spm_vol([MapName(1:end-8) '.nii']);
        end
        MAP(1).AnatMask = spm_vol(fullfile(spm('Dir','se_anatomy'),'AnatMask.nii'));
        MAP(1).Macro = spm_vol(fullfile(spm('Dir','se_anatomy'),'MacroLabels.nii'));
        MAP(1).MLabels = load('Macro.mat');
        for i=1:size(MAP,2)
            MAP(i).ref = strrep(MAP(i).ref,'\',filesep);
            MAP(i).ref = strrep(MAP(i).ref,'/',filesep);
            MAP(i).PMap = spm_vol([spm('Dir','se_anatomy') filesep 'PMaps' filesep spm_str_manip(MAP(i).ref,'t')]);
            for side = [-1 1]
                MAP(i).AM(side+2)   = mean(MAP(i).allZ(MAP(i).allLR == side));
                tmp                 = histc(MAP(i).allZ(MAP(i).allLR == side),[0.01:.1:1.01]); tmp = tmp(1:end-1);
                MAP(i).nA{side+2}   = tmp(1:end)/sum(tmp)*100;
            end
        end
        if flip ~= -1; defaults.analyze.flip = flip; end

        try
            MAP(1).cmap;
        catch
            load(fullfile(spm('Dir','se_anatomy'),'cmap2.mat'));
            MAP(1).cmap = cmap;
        end
        
    try
        MAP(1).rend;
    catch
        load(fullfile(spm('Dir','spm'),'rend','render_single_subj.mat'));
        if (exist('rend') ~= 1), % Assume old format...
            rend = cell(size(Matrixes,1),1);
            for i=1:size(Matrixes,1),
                rend{i}=struct('M',eval(Matrixes(i,:)),...
                    'ren',eval(Rens(i,:)),...
                    'dep',eval(Depths(i,:)));
                rend{i}.ren = rend{i}.ren/max(max(rend{i}.ren));
            end;
        end;

        for i=1:length(rend),
            rend{i}.max=0;
            rend{i}.data = cell(1,1);
            if issparse(rend{i}.ren),
                d = size(rend{i}.ren);
                B1 = spm_dctmtx(d(1),d(1));
                B2 = spm_dctmtx(d(2),d(2));
                rend{i}.ren = B1*rend{i}.ren*B2';
                rend{i}.dep = exp(B1*rend{i}.dep*B2')-1;
            end;
            msk = find(rend{i}.ren>1);rend{i}.ren(msk)=1;
            msk = find(rend{i}.ren<0);rend{i}.ren(msk)=0;
        end;
        MAP(1).rend = rend;    
    end
        
    MAP(1).GM = spm_vol(fullfile(spm('Dir','se_anatomy'),'wgrey.nii'));  MAP(1).GM.pinfo(1) = 0.00392156885936856270;    
    MAP(1).WM = spm_vol(fullfile(spm('Dir','se_anatomy'),'wwhite.nii')); MAP(1).WM.pinfo(1) = 0.00392156885936856270;       

        
        
        
        
    case 'IMAGE'
        OLimg = spm_vol(spm_select(1,'image','Select overlay image'));
        xSPM = struct('xVOL',{});
        xSPM(1).xVOL.M = OLimg.mat; xSPM.xVOL.iM = inv(OLimg.mat);
        xSPM.xVOL.DIM = OLimg.dim(1:3); xSPM.xVOL.VOX = [abs(OLimg.mat(1,1)); abs(OLimg.mat(2,2)); abs(OLimg.mat(3,3))];
        xSPM.xVOL.orig = spm_input('Origin corrected to AC ?','+1','y/n',[1,2],2);

        pmult  = spm_input(['Premultiply data by '],'+0','r',1,1);

        u  = spm_input(['Hight threshold (0 for none)'],'+0','r',0,1);
        k  = spm_input(['Extend threshold (0 for none)'],'+0','r',0,1);

        if xSPM.xVOL.orig == 2
            xSPM.xVOL.M(1,4) = xSPM.xVOL.M(1,4) + xtr;
            xSPM.xVOL.M(2,4) = xSPM.xVOL.M(2,4) + ytr;
            xSPM.xVOL.M(3,4) = xSPM.xVOL.M(3,4) + ztr;
            try, set(st.mAssign,'ToolTipString','Corrdinates given in anatomical MNI space, orig. MNI coordinates in brackets'); end
        end

        SPM = struct('title','');
        
        SPM.title = spm_input('Title','+1','s',spm_str_manip(OLimg.fname,'rt'));


        SPM.XYZ = []; SPM.Z = [];
        spm_progress_bar('Init',OLimg.dim(3),'Preparing data');
        for p = 1:OLimg.dim(3)
            d = spm_slice_vol(OLimg,inv(spm_matrix([0 0 -p 0 0 0 1 1 1])),OLimg.dim(1:2),[0, NaN])*pmult;
            if any(any(d>=u))
                [i,j] = find(d>u);
                SPM.XYZ = [SPM.XYZ [i'; j'; p*ones(1,size(j,1))]];
                SPM.Z = [SPM.Z d(find(d>u))'];
            end
            spm_progress_bar('Set',p);
        end
        if any(isinf(SPM.Z))
            SPM.Z(isinf(SPM.Z)) = max(SPM.Z(~isinf(SPM.Z)));
            spm('alert!','Infinitive values replace by maximum finite value');
        end
        spm_progress_bar('Clear')

        k = max(k,2);
        
        if numel(SPM.XYZ)>0


            A     = spm_clusters(SPM.XYZ);
            Q     = [];
            for i = 1:max(A)
                j = find(A == i);
                if length(j) >= k; Q = [Q j]; end
            end

            SPM.Z     = SPM.Z(:,Q);
            SPM.XYZ   = SPM.XYZ(:,Q);

            SPM.XYZmm = xSPM.xVOL.M * [SPM.XYZ; ones(1,size(SPM.XYZ,2))]; SPM.XYZmm = SPM.XYZmm(1:3,:);
            SPM.XYZp = [round(SPM.XYZmm(1,:)-MAP(1).MaxMap.mat(1,4)); round(SPM.XYZmm(2,:)-MAP(1).MaxMap.mat(2,4)); round(SPM.XYZmm(3,:)-MAP(1).MaxMap.mat(3,4))];
            SPM.u = u; SPM.k = k;
        end

    case 'SPM'
        st.SPM = spm('FnBanner');
        switch st.SPM;
            case 'SPM99'
                [SPM,VOL,xX,xCon,xSDM] = spm_getSPM;
                xSPM = struct('xVOL',VOL,'xX',xX,'xCon',xCon,'Vbeta',xSDM.Vbeta);
                load(fullfile(SPM.swd,'SPMcfg.mat'));
                SPM.FWHM = VOL.FWHM;
                SPM.VOX = VOL.VOX';
                xSPM.xY.VY = VY;
                xSPM.xsDes = xsDes;
                SPM.S = xSDM.S;
            otherwise
                [xSPM,SPM] = spm_getSPM;
                if isfield(xSPM.xVol,'mat');
                    xSPM.xVol.M = xSPM.xVol.mat
                end
                xSPM.xVol.VOX = [abs(xSPM.xVol.M(1,1)) abs(xSPM.xVol.M(2,2)) abs(xSPM.xVol.M(3,3))];
                xSPM.xVOL = xSPM.xVol;
                xSPM = rmfield(xSPM,'xVol');
        end
        if numel(SPM.XYZ)>0
            xSPM.xVOL.orig = spm_input('Origin corrected to AC ?','+1','y/n',[1,2],2);
            if xSPM.xVOL.orig == 2
                xSPM.xVOL.M(1,4) = xSPM.xVOL.M(1,4) + xtr;
                xSPM.xVOL.M(2,4) = xSPM.xVOL.M(2,4) + ytr;
                xSPM.xVOL.M(3,4) = xSPM.xVOL.M(3,4) + ztr;
                SPM.XYZmm = xSPM.xVOL.M * [SPM.XYZ; ones(1,size(SPM.XYZ,2))]; SPM.XYZmm = SPM.XYZmm(1:3,:);
                try, set(st.mAssign,'ToolTipString','Corrdinates given in anatomical MNI space, orig. MNI coordinates in brackets'); end
            end

            SPM.XYZp = [round(SPM.XYZmm(1,:)-MAP(1).MaxMap.mat(1,4)); round(SPM.XYZmm(2,:)-MAP(1).MaxMap.mat(2,4)); round(SPM.XYZmm(3,:)-MAP(1).MaxMap.mat(3,4))];
        end

    case 'cluster'
        [A] = spm_clusters(SPM.XYZ);
        CLUSTER = struct('Z',{},'XYZ',{},'XYZmm',{},'voxel',{},'nrMax',{},'maxZ',{},'maxM',{});
        Buffer = struct('Z',{},'XYZ',{},'XYZmm',{},'voxel',{},'nrMax',{},'maxZ',{},'maxM',{});
        mx = zeros(1,size(unique(A),2));
        for cl=1: size(unique(A),2)
            CLUSTER(cl).Z = SPM.Z(:,A == cl);
            CLUSTER(cl).XYZ = SPM.XYZ(:,A == cl);
            CLUSTER(cl).XYZmm = SPM.XYZmm(:,A == cl);
            CLUSTER(cl).XYZp = [CLUSTER(cl).XYZmm(1,:)-MAP(1).MaxMap.mat(1,4); CLUSTER(cl).XYZmm(2,:)-MAP(1).MaxMap.mat(2,4); CLUSTER(cl).XYZmm(3,:)-MAP(1).MaxMap.mat(3,4)];
            [CLUSTER(cl).voxel CLUSTER(cl).maxZ CLUSTER(cl).maxM dummy] = spm_max(SPM.Z(:,A == cl),SPM.XYZ(:,A == cl));
            
            [CLUSTER(cl).maxZ I] = sort(CLUSTER(cl).maxZ,'descend');
             CLUSTER(cl).maxM    = CLUSTER(cl).maxM(:,I);
            
             
             
             Q    = zeros(size(CLUSTER(cl).maxZ));
             Q(1) = 1;
             for d = 2:numel(CLUSTER(cl).maxZ)
                 if CLUSTER(cl).maxZ(d) ~= CLUSTER(cl).maxZ(d-1)
 
                     distX = sqrt(sum(((CLUSTER(cl).maxM(:,d)-CLUSTER(cl).maxM(:,1)) .* abs([xSPM.xVOL.M(1,1) xSPM.xVOL.M(2,2) xSPM.xVOL.M(3,3)]')).^2));

                     kif = find(Q(2:d-1)>0);
                     for dd = kif
                             distX = min(distX,sqrt(sum(((CLUSTER(cl).maxM(:,d)-CLUSTER(cl).maxM(:,dd)) .* abs([xSPM.xVOL.M(1,1) xSPM.xVOL.M(2,2) xSPM.xVOL.M(3,3)]')).^2)));
                             if distX<4
                                 break
                             end
                     end
                     if distX>4;
                         % sqrt(xSPM.xVOL.M(1,1)^2+xSPM.xVOL.M(2,2)^2+xSPM.xVOL.M(3,3)^2)
                         Q(d) = 1;
                     end
                 end
             end
            
            CLUSTER(cl).maxZ = CLUSTER(cl).maxZ(Q>0);
            CLUSTER(cl).maxM = CLUSTER(cl).maxM(:,Q>0);
            
            if isempty(max(CLUSTER(cl).maxZ)); % binary image
                CLUSTER(cl).voxel = prod(size(find(A == cl)));
            end
            CLUSTER(cl).nrMax = size(CLUSTER(cl).maxZ,2);
            mx(cl) = CLUSTER(cl).voxel(1);
        end
        CLUSTER(1).VOL = zeros(MAP(1).MaxMap.dim(1:3));
        [B,INDEX] = sort(mx); INDEX = fliplr(INDEX); CLUSTER = CLUSTER(INDEX);
        for i=1:size(unique(A),2); RIX(i) = find(INDEX == i); end

        for i=1:size(A,2)
            try
                CLUSTER(1).VOL([SPM.XYZp(1,i)-1:SPM.XYZp(1,i)+1],...
                    [SPM.XYZp(2,i)-1:SPM.XYZp(2,i)+1],...
                    [SPM.XYZp(3,i)-1:SPM.XYZp(3,i)+1]) = RIX(A(i));
            end
        end


    case 'groupStat'
        st.SPM = spm('FnBanner');
        group = [];
        def = spm_input('Load an existing group ?','+1','y/n',[1,0],1);
        if def == 0 % Create new group
            cwd = pwd;
            nsubjects = spm_input('Number of subjects','+1','i',1);
            group = struct('xSPM',{});
            for i=1:nsubjects
                a = spm_select(1,'^SPM\.mat$',['Select SPM.mat for subject ' num2str(i)]);
                load(a);
                st.SPM = spm('FnBanner');
                switch st.SPM;
                    case 'SPM99'
                        for ii=1:size(Vbeta,1)
                            try
                                f = cell2struct(Vbeta(ii,1),'name',1);
                                group(i).xSPM.Vbeta(ii) = spm_vol([spm_str_manip(a,'H') filesep spm_str_manip(f.name,'t')]);

                            catch
                                f = struct('name',Vbeta{ii,1})
                                group(i).xSPM.Vbeta(ii) = spm_vol([spm_str_manip(a,'H') filesep spm_str_manip(f.name,'t')]);
                            end
                        end
                        group(i).xSPM.xX = xX;
                        group(i).xSPM.Sess = Sess;
                    otherwise
                        for ii=1:size(SPM.Vbeta,2)
                            group(i).xSPM.Vbeta(ii)  = spm_vol([spm_str_manip(a,'H') '/' SPM.Vbeta(1,ii).fname]);
                        end
                        group(i).xSPM.xX = SPM.xX;
                        group(i).xSPM.VY = SPM.xY.VY;
                        group(i).xSPM.XYZ = SPM.xVol.XYZ;
                        group(i).xSPM.xVOL.M = SPM.xVol.M;
                        group(i).xSPM.xVOL.iM = inv(SPM.xVol.M);
                        group(i).xSPM.xCon = SPM.xCon;
                        group(i).xSPM.Sess = SPM.Sess;
                        group(i).xSPM.swd = SPM.swd;
                end
            end
            load(spm_select(1,'^SPM\.mat$',['Select SPM.mat for the 2nd level analysis']));
            group(1).xSPM.xVOL.orig = spm_input('Origin corrected to AC ?','+1','y/n',[1,2],2);

            if group(1).xSPM.xVOL.orig == 2
                for i=1:nsubjects
                    group(i).xSPM.xVOL.M(1,4) = group(i).xSPM.xVOL.M(1,4) + xtr;
                    group(i).xSPM.xVOL.M(2,4) = group(i).xSPM.xVOL.M(2,4) + ytr;
                    group(i).xSPM.xVOL.M(3,4) = group(i).xSPM.xVOL.M(3,4) + ztr;
                    try, set(st.mAssign,'ToolTipString','Corrdinates given in anatomical MNI space, orig. MNI coordinates in brackets'); end
                    group(i).xSPM.XYZmm = group(i).xSPM.xVOL.M * [group(i).xSPM.XYZ; ones(1,size(group(i).xSPM.XYZ,2))];
                    group(i).xSPM.XYZmm = group(i).xSPM.XYZmm(1:3,:);
                end
            end

            name = spm_input('Name of group','!+1','s','GroupAnalysis');
            cd(cwd);
            save(deblank(fullfile(pwd,name)),'group','xSPM','SPM')

        else    % load existing group
            load(spm_select(1,'mat',['Select information about the subjects']));
            xSPM.xVOL = group(1).xSPM.xVOL;
        end
        spm('alert!','Statistical information loaded');



    case 'group'
        st.SPM = spm('FnBanner');
        def = spm_input('Load an existing group ?','+1','y/n',[1,0],1);
        if def == 0 % Create new group
            cwd = pwd;
            nsubjects = spm_input('Number of subjects','+1','i',1);
            group = struct('xSPM',{});
            for i=1:nsubjects
                a = spm_select(1,'^SPM\.mat$',['Select SPM.mat for subject ' num2str(i)]);
                load(a);
                switch st.SPM;
                    case 'SPM99'
                        for ii=1:size(Vbeta,1)
                            try
                                f = cell2struct(Vbeta(ii,1),'name',1);
                                group(i).xSPM.Vbeta(ii) = spm_vol([spm_str_manip(a,'H') filesep spm_str_manip(f.name,'t')]);

                            catch
                                f = struct('name',Vbeta{ii,1})
                                group(i).xSPM.Vbeta(ii) = spm_vol([spm_str_manip(a,'H') filesep spm_str_manip(f.name,'t')]);
                            end
                        end
                        group(i).xSPM.xX = xX;
                        group(i).xSPM.Sess = Sess;
                    otherwise
                        for ii=1:size(SPM.Vbeta,2)
                            group(i).xSPM.Vbeta(ii)  = spm_vol([spm_str_manip(a,'H') '/' SPM.Vbeta(1,ii).fname]);
                        end
                        group(i).xSPM.xX = SPM.xX;
                        group(i).xSPM.VY = SPM.xY.VY;
                        group(i).xSPM.XYZ = SPM.xVol.XYZ;
                        group(i).xSPM.xVOL.M = SPM.xVol.M;
                        group(i).xSPM.xVOL.iM = inv(SPM.xVol.M);
                        group(i).xSPM.xCon = SPM.xCon;
                        group(i).xSPM.Sess = SPM.Sess;
                        group(i).xSPM.swd = SPM.swd;
                end
            end
            spm('alert"',{'Please specify now the second level statistic'},'Next Step',0,1);
            se_getMap('SPM','');
            if xSPM.xVOL.orig == 2
                for i=1:nsubjects
                    group(i).xSPM.xVOL.M(1,4) = group(i).xSPM.xVOL.M(1,4) + xtr;
                    group(i).xSPM.xVOL.M(2,4) = group(i).xSPM.xVOL.M(2,4) + ytr;
                    group(i).xSPM.xVOL.M(3,4) = group(i).xSPM.xVOL.M(3,4) + ztr;
                    try, set(st.mAssign,'ToolTipString','Corrdinates given in anatomical MNI space, orig. MNI coordinates in brackets'); end
                    group(i).xSPM.XYZmm = group(i).xSPM.xVOL.M * [group(i).xSPM.XYZ; ones(1,size(group(i).xSPM.XYZ,2))];
                    group(i).xSPM.XYZmm = group(i).xSPM.XYZmm(1:3,:);
                end
            end
            name = spm_input('Name of group','!+1','s','GroupAnalysis');
            cd(cwd);
            save(deblank(fullfile(pwd,name)),'group','xSPM','SPM')
        else    % load existing group
            load(spm_select(1,'mat',['Select information about the subjects']));
            spm('alert"',{'Please specify now the second level statistic'},'Next Step',0,1);
            se_getMap('SPM','');
        end
end

