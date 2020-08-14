function [MLev, VLev, m_anz] = se_overlap(cl)

global st;
global MAP;
global SPM;
global xSPM;
global CLUSTER;

GV = [MAP.GV];

doPlot = 0;

if strcmp(st.prog,'PX')
    MLev = {};
    xyzTMP = [CLUSTER.XYZmm; 1];
    xyz = inv(MAP(1).MaxMap.mat) * xyzTMP;
    wXYZ = ['x= ' sprintf('%+3.0f', xyzTMP(1)) '  y= ' sprintf('%+2.0f', xyzTMP(2)) '  z= ' sprintf('%+2.0f', xyzTMP(3))];
    if xSPM.xVOL.orig == 2
        wXYZ = [wXYZ '  (MNI: ' sprintf('%+3.0f', xyzTMP(1)) '/' sprintf('%+2.0f', xyzTMP(2)+4) '/' sprintf('%+2.0f', xyzTMP(3)-5) ')'];
    end
    MLev = {wXYZ};
    
        ML = round(spm_sample_vol(MAP(1).Macro,xyz(1),xyz(2),xyz(3),0)');
        if ML > 0;
            MLl = MAP(1).MLabels.Labels{ML};
            MLev(size(MLev,1)+1,1) = {[MAP(1).MLabels.Labels{ML}]};
        end
    
        ProbMax(1) = spm_sample_vol(MAP(1).MaxMap,xyz(1),xyz(2),xyz(3),0);
        for PM = 1:size(MAP,2)
            ProbMax(PM+1) = 100*spm_sample_vol(MAP(PM).PMap,xyz(1),xyz(2),xyz(3),0);
        end
        ProbMax = round(ProbMax);
        
        if any(ProbMax(2:end))
            Probs = find(ProbMax(2:end)>0); [value sortP]= sort(ProbMax(Probs+1));
            
            if numel(find(GV== ProbMax(1)))
                try
                    posMPM = find(Probs(sortP) == find(GV== ProbMax(1)) );
                    if posMPM~=numel(sortP)
                        sortP = [sortP(find(Probs(sortP) ~= find(GV== ProbMax(1)) )) sortP(posMPM)];
                    end
                end
            end

            for getPr = size(Probs,2):-1:max([size(Probs,2)-2 1])
                    if getPr == size(Probs,2)
                        [Ploc, Pmin, Pmax] = MinMax(MAP(Probs(sortP(getPr))).PMap,xyz);
                        if ProbMax(1,1)>99 
                             MLev(size(MLev,1)+1,1) = {['-> Assigned to ' MAP(Probs(sortP(getPr))).name...
                                  ',  Probability: ' int2str(Ploc) '%  [' int2str(Pmin) '-' int2str(Pmax) '%]']} ;
                        else
                            ML = round(spm_sample_vol(MAP(1).Macro,xyz(1),xyz(2),xyz(3),0)');
                            
                            
                            if ML > 0;
                                MLl = MAP(1).MLabels.Labels{ML};
                                MLev(size(MLev,1)+1,1) = {MAP(1).MLabels.Labels{ML}};
                                MLev(size(MLev,1)+1,1) = {['Probability for ' MAP(Probs(sortP(getPr))).name ': ' int2str(Ploc) '%  [' int2str(Pmin) '-' int2str(Pmax) '%]']} ;
                            else
                                MLev(size(MLev,1)+1,1) = {['-> Not assigned'...
                                  ',  Probability for ' MAP(Probs(sortP(getPr))).name ': ' int2str(Ploc) '%  [' int2str(Pmin) '-' int2str(Pmax) '%]']} ;
                            end
                        end
                    else
                        if (ProbMax(1+Probs(sortP(getPr))))>=5
                            [Ploc, Pmin, Pmax] = MinMax(MAP(Probs(sortP(getPr))).PMap,xyz);
                            MLev(size(MLev,1)+1,1) = {['  Probability for '...
                              MAP(Probs(sortP(getPr))).name ': ' int2str(Ploc) '%  [' int2str(Pmin) '-' int2str(Pmax) '%]']} ;
                        end
                    end
            end
            MLev(size(MLev,1)+1,1) = {''};
        end

        if ~any(ProbMax(2:end))
            if ProbMax(1)>0
                MLev(size(MLev,1)+1,1) = {['-> Assigned to ' MAP(find([MAP.GV]==ProbMax(1))).name]};
            else
                MLev(size(MLev,1)+1,1) = {['-> Not found in any probability map']};
            end
        end
        
    VLev = ''; m_anz = 0;
else
    
    
    
    % Maxima
    MLev = {}; m_anz = 0;
    if (min(CLUSTER(cl).Z) == max(CLUSTER(cl).Z))   % Binary image
        CLUSTER(cl).maxM = mean(CLUSTER(cl).XYZ,2); m_anz = 1; CLUSTER(cl).maxZ = 0;
		xyzTMP   = xSPM.xVOL.M * [CLUSTER(cl).maxM(1); CLUSTER(cl).maxM(2); CLUSTER(cl).maxM(3); 1];		%-Reorientation matrix
        xyz(1:4,1) = inv(MAP(1).MaxMap.mat) * xyzTMP(1:4,1);
        for PM = 1:size(MAP,2)
            ProbMax(1:size(xyz,2),PM+1) = 100*spm_sample_vol(MAP(PM).PMap,xyz(1,:),xyz(2,:),xyz(3,:),0)';
        end
        ProbMax(:,1) = spm_sample_vol(MAP(1).MaxMap,xyz(1,:),xyz(2,:),xyz(3,:),0)'; ProbMax(:,1) = ProbMax(:,1) .* (ProbMax(:,1) > 99);
        ProbMax = round(ProbMax);
        wXYZ = ['Center of gravity:  x= ' sprintf('%+3.0f', xyzTMP(1)) '  y= ' sprintf('%+2.0f', xyzTMP(2)) '  z= ' sprintf('%+2.0f', xyzTMP(3))];
        if xSPM.xVOL.orig == 2
            wXYZ = [wXYZ '  (MNI: ' sprintf('%+3.0f', xyzTMP(1)) '/' sprintf('%+2.0f', xyzTMP(2)+4) '/' sprintf('%+2.0f', xyzTMP(3)-5) ')'];
        end
        MLev = {wXYZ};

        indxx = 1;
        if any(ProbMax(:))
            ML = round(spm_sample_vol(MAP(1).Macro,xyz(1,indxx),xyz(2,indxx),xyz(3,indxx),0)');
            if ML > 0;
                MLl = MAP(1).MLabels.Labels{ML};
                wXYZ = [wXYZ '   (' MLl ')'];
            end
            Probs = find(ProbMax(2:end)>0); [value sortP]= sort(ProbMax(Probs+1));
                posMPM = find(Probs(sortP) == find(GV== ProbMax(1)) );

                if posMPM~=numel(sortP)
                    sortP = [sortP(find(Probs(sortP) ~= find(GV== ProbMax(1)) )) sortP(posMPM)];
                end

                for getPr = size(Probs,2):-1:max([size(Probs,2)-2 1])
                if getPr == size(Probs,2)
                    [Ploc, Pmin, Pmax] = MinMax(MAP(Probs(sortP(getPr))).PMap,xyz(:));
                    if ProbMax(1) 
                         MLev(size(MLev,1)+1,1) = {['-> Assigned to ' MAP(Probs(sortP(getPr))).name...
                              ',  Probability: ' int2str(Ploc) '%  [' int2str(Pmin) '-' int2str(Pmax) '%]']} ;
                    else
                        ML = round(spm_sample_vol(MAP(1).Macro,xyz(1),xyz(2),xyz(3),0)');
                        
                        if ML > 0;
                            MLl = MAP(1).MLabels.Labels{ML};
                            MLev(size(MLev,1)+1,1) = {['-> ' MAP(1).MLabels.Labels{ML}...
                                  ',  Probability for ' MAP(Probs(sortP(getPr))).name ': ' int2str(Ploc) '%  [' int2str(Pmin) '-' int2str(Pmax) '%]']} ;
                        else
                            MLev(size(MLev,1)+1,1) = {['-> Not assigned'...
                                  ',  Probability for ' MAP(Probs(sortP(getPr))).name ': ' int2str(Ploc) '%  [' int2str(Pmin) '-' int2str(Pmax) '%]']} ;
                        end
                    end
                else
                    if sum(ProbMax(1+Probs(sortP(getPr+1:size(Probs,2))))) < 100 & (ProbMax(1+Probs(sortP(getPr))))>5
                        [Ploc, Pmin, Pmax] = MinMax(MAP(Probs(sortP(getPr))).PMap,xyz(:));
                        MLev(size(MLev,1)+1,1) = {['  Probability for '...
                          MAP(Probs(sortP(getPr))).name ': ' int2str(Ploc) '%  [' int2str(Pmin) '-' int2str(Pmax) '%]']} ;
                    end
                end
            end
        else
            ML = round(spm_sample_vol(MAP(1).Macro,xyz(1),xyz(2),xyz(3),0)');
            if ML > 0;
                MLl = MAP(1).MLabels.Labels{ML};
                MLev(size(MLev,1)+1,1) = {['  ' MAP(1).MLabels.Labels{ML}]};
            else
                MLev(size(MLev,1)+1,1) = {['-> Not found in any probability map']};
            end
        end
        MLev(size(MLev,1)+1,1) = {''}; 
    
    
    
    else    % Kein binärimage!
        if isempty(max(CLUSTER(cl).maxZ))
            CLUSTER(cl).maxZ = CLUSTER(cl).Z(:,find(CLUSTER(cl).Z == max(CLUSTER(cl).Z)));
            CLUSTER(cl).maxM = CLUSTER(cl).XYZ(:,find(CLUSTER(cl).Z == max(CLUSTER(cl).Z)));
        end
        [B,INDEX] = sort(CLUSTER(cl).maxZ); xyzTMP = []; xyz = [];  ProbMax = [];
        
        for indxx  = INDEX(size(B,2):-1:max([size(B,2)-6 1])); 
		    xyzTMP(1:4,indxx)   = xSPM.xVOL.M * [CLUSTER(cl).maxM(1,indxx); CLUSTER(cl).maxM(2,indxx); CLUSTER(cl).maxM(3,indxx); 1];		%-Reorientation matrix
            xyz(1:4,indxx) = inv(MAP(1).MaxMap.mat) * xyzTMP(1:4,indxx);
        end
        
        for PM = 1:size(MAP,2)
            ProbMax(1:size(xyz,2),PM+1) = 100*spm_sample_vol(MAP(PM).PMap,xyz(1,:),xyz(2,:),xyz(3,:),0)';
        end
        ProbMax(:,1) = spm_sample_vol(MAP(1).MaxMap,xyz(1,:),xyz(2,:),xyz(3,:),0)'; ProbMax(:,1) = ProbMax(:,1) .* (ProbMax(:,1) > 99);
        ProbMax = round(ProbMax);
        for indxx  = INDEX(size(B,2):-1:max([size(B,2)-5 1])); 
            m_anz = m_anz+1;
            wXYZ = ['(' int2str(m_anz) ')  ' 'x= ' sprintf('%+3.0f', xyzTMP(1,indxx)) '  y= ' sprintf('%+2.0f', xyzTMP(2,indxx)) '  z= ' sprintf('%+2.0f', xyzTMP(3,indxx))...
                    ';  T = ' sprintf('%+2.2f', CLUSTER(cl).maxZ(indxx))];
            if xSPM.xVOL.orig == 2
                wXYZ = [wXYZ '  (MNI: ' sprintf('%+3.0f', xyzTMP(1,indxx)) '/' sprintf('%+2.0f', xyzTMP(2,indxx)+4) '/' sprintf('%+2.0f', xyzTMP(3,indxx)-5) ')'];
            end
            if prod(size(MLev)); MLev(size(MLev,1)+1,1) = {wXYZ}; else MLev = {wXYZ}; end                       
            
                ML = round(spm_sample_vol(MAP(1).Macro,xyz(1,indxx),xyz(2,indxx),xyz(3,indxx),0)');
                            GM = spm_sample_vol(MAP(1).GM,xyz(1,indxx),xyz(2,indxx),xyz(3,indxx),0)';
                            WM = spm_sample_vol(MAP(1).WM,xyz(1,indxx),xyz(2,indxx),xyz(3,indxx),0)';
                            
                            if GM>WM
                                was = ['Probably GM [' int2str(GM*100/128) '/' int2str(WM*100/128) ']'];
                            else
                                was = ['Probably WM [' int2str(WM*100/128) '/' int2str(GM*100/128) ']'];
                            end
                if ML > 0;
                    MLl = MAP(1).MLabels.Labels{ML};
                    MLev(size(MLev,1)+1,1) = {['(' MLl ' - ' was ')']};
                else
                    MLev(size(MLev,1)+1,1) = {['(' was ')']};
                end
            
            
            if any(ProbMax(indxx,:))
                Probs = find(ProbMax(indxx,2:end)>0); [value sortP]= sort(ProbMax(indxx,Probs+1));
                
                if numel(find(GV== ProbMax(indxx,1)))==1
                    try
                        posMPM = find(Probs(sortP) == find(GV== ProbMax(indxx,1)) );
                        if posMPM~=numel(sortP)
                            sortP = [sortP(find(Probs(sortP) ~= find(GV== ProbMax(indxx,1)) )) sortP(posMPM)];
                        end
                    end
                end

                for getPr = size(Probs,2):-1:max([size(Probs,2)-2 1])
                        if getPr == size(Probs,2)
                            [Ploc, Pmin, Pmax] = MinMax(MAP(Probs(sortP(getPr))).PMap,xyz(:,indxx));
                            if ProbMax(indxx,1) 
                                 MLev(size(MLev,1)+1,1) = {['-> Assigned to ' MAP(Probs(sortP(getPr))).name...
                                      ',  Probability: ' int2str(Ploc) '%  [' int2str(Pmin) '-' int2str(Pmax) '%]']} ;
                            else
                                MLev(size(MLev,1)+1,1) = {['-> Not assigned'...
                                          ',  Probability for ' MAP(Probs(sortP(getPr))).name ': ' int2str(Ploc) '%  [' int2str(Pmin) '-' int2str(Pmax) '%]']} ;
                            end
                        else
                            if sum(ProbMax(indxx,1+Probs(sortP(getPr+1:size(Probs,2))))) < 100 & (ProbMax(indxx,1+Probs(sortP(getPr))))>5
                                [Ploc, Pmin, Pmax] = MinMax(MAP(Probs(sortP(getPr))).PMap,xyz(:,indxx));
                                MLev(size(MLev,1)+1,1) = {['   Probability for '...
                                  MAP(Probs(sortP(getPr))).name ': ' int2str(Ploc) '%  [' int2str(Pmin) '-' int2str(Pmax) '%]']} ;
                            end
                        end
                end
            else
                ML = round(spm_sample_vol(MAP(1).Macro,xyz(1,indxx),xyz(2,indxx),xyz(3,indxx),0)');
                MLev(size(MLev,1)+1,1) = {['-> Not found in any probability map']};
            end
            MLev(size(MLev,1)+1,1) = {''};
        end
    end
    
    
    
    
    
    
    
    
    % Volume
        
    if ~isfield(CLUSTER,'gtmp') || prod(size(CLUSTER(cl).gtmp)) < 1 | ~isfield(CLUSTER,'gZ') || prod(size(CLUSTER(cl).gZ)) < 1
        tmp   = CLUSTER(cl).XYZp;
        tmpZ  = CLUSTER(cl).Z;
        tms = sign(CLUSTER(cl).XYZmm); 
        clear gtmp; gZ = []; mult = 1;
        for xbx = 0:.5:xSPM.xVOL.VOX(1)-1
            for ybx = 0:.5:xSPM.xVOL.VOX(2)-1
                for zbx = 0:.5:xSPM.xVOL.VOX(3)-1
                    if exist('gtmp')
                        gtmp = [gtmp(1,:) tmp(1,:)-(tms(1,:)*xbx);...
                                gtmp(2,:) tmp(2,:)-(tms(2,:)*ybx);...
                                gtmp(3,:) tmp(3,:)-(tms(3,:)*zbx)];
                        gZ   = [gZ tmpZ];
                    else
                        gtmp = [tmp(1,:)-(tms(1,:)*xbx);...
                                tmp(2,:)-(tms(2,:)*ybx);...
                                tmp(3,:)-(tms(3,:)*zbx)];
                            
                        gZ   = [tmpZ];
                    end
                end
            end
        end
        gtmp      = round(gtmp);
        [b, m, n] = unique(gtmp','rows');
        
        gtmp = gtmp(:,m);
        gZ   = gZ(:,m);
        
        CLUSTER(cl).gtmp = gtmp;
        CLUSTER(cl).gZ   = gZ;
        CLUSTER(cl).mult = length(gZ)/length(tmpZ);
    end
    
    Msk    = round(spm_sample_vol(MAP(1).AnatMask,CLUSTER(cl).gtmp(1,:),CLUSTER(cl).gtmp(2,:),CLUSTER(cl).gtmp(3,:),0));
    Volume = spm_sample_vol(MAP(1).MaxMap,  CLUSTER(cl).gtmp(1,:),CLUSTER(cl).gtmp(2,:),CLUSTER(cl).gtmp(3,:),0); % Corresponding voxel + next towards origin in all directions
       
    
    VLev = {}; vp = [];
    for v=1:size(MAP,2)
        if any(Volume(Msk == 2) == GV(v))
            vp = [vp; size(find(Volume(Msk == 2) == GV(v)),2) v -1];
        end
        if any(Volume(Msk == 1) == GV(v))
            vp = [vp; size(find(Volume(Msk == 1) == GV(v)),2) v 1];
        end
    end
    vp = sortrows(vp);
    if any(any(vp))
        doPlot = 1;
    end
    for v=size(vp,1):-1:max([size(vp,1)-5 1])
            if 1% (100*vp(v,1)/(CLUSTER(cl).voxel(1)*CLUSTER(cl).mult)) >1
                
                if vp(v,3) == -1; 
                    side = 'left'; vol = MAP(vp(v,2)).VOL(1);
                else; 
                    side = 'right'; vol = MAP(vp(v,2)).VOL(2);
                end
                if v==size(vp,1)
                    VLev = {[sprintf('%2.1f', (100*vp(v,1)/(CLUSTER(cl).voxel(1)*CLUSTER(cl).mult))) '% in ' side ' ' MAP(vp(v,2)).name...
                        ' (' sprintf('%3.1f', 100*vp(v,1)/vol)  '% activated)']};
                else
                    VLev(size(VLev,1)+1,1) = {[sprintf('%2.1f', (100*vp(v,1)/(CLUSTER(cl).voxel(1)*CLUSTER(cl).mult))) '%  in ' side ' ' MAP(vp(v,2)).name...
                        ' (' sprintf('%3.1f', 100*vp(v,1)/vol)  '% activated)']};
                end
            end
    end
    
    
    if ~prod(size(VLev))
        VLev = {'Volume not matching any probability map'};
    end


    if size(VLev,1)>4; VLev = VLev(1:4,:); end
    
    if numel(CLUSTER(cl).Z) == 1
        doRender = 0;
        doPlot   = 0;
    else
        doRender = 1;
    end
    
se_plotCorrespond2(CLUSTER(cl));    

%     if  doPlot
%         try
% %            notNow
%             [PQ, Z, nP, index] = se_plotCorrespond2(CLUSTER(cl));
%             CLUSTER(cl).PQ = PQ;
%             CLUSTER(cl).xZ = Z;
%             CLUSTER(cl).nP = nP;
%             CLUSTER(cl).index = index;
%         catch
%            try; delete(st.FX); end
%         end
%     else
%         try; delete(st.FX); end
%     end
    
    if doRender
        try
            se_render(cl)
        end
    end
    
end


function [Ploc, Pmin, Pmax] = MinMax(map,tmp)
                               sample = spm_sample_vol(map,...
                                            [tmp(1)-1 tmp(1) tmp(1)+1 tmp(1)-1 tmp(1) tmp(1)+1 tmp(1)-1 tmp(1) tmp(1)+1 ...
                                             tmp(1)-1 tmp(1) tmp(1)+1 tmp(1)-1 tmp(1) tmp(1)+1 tmp(1)-1 tmp(1) tmp(1)+1 ...
                                             tmp(1)-1 tmp(1) tmp(1)+1 tmp(1)-1 tmp(1) tmp(1)+1 tmp(1)-1 tmp(1) tmp(1)+1], ....
                                            [tmp(2)-1 tmp(2)-1 tmp(2)-1 tmp(2) tmp(2) tmp(2) tmp(2)+1 tmp(2)+1 tmp(2)+1 ...
                                             tmp(2)-1 tmp(2)-1 tmp(2)-1 tmp(2) tmp(2) tmp(2) tmp(2)+1 tmp(2)+1 tmp(2)+1 ...
                                             tmp(2)-1 tmp(2)-1 tmp(2)-1 tmp(2) tmp(2) tmp(2) tmp(2)+1 tmp(2)+1 tmp(2)+1], ...
                                            [tmp(3)-1 tmp(3)-1 tmp(3)-1 tmp(3)-1 tmp(3)-1 tmp(3)-1 tmp(3)-1 tmp(3)-1 tmp(3)-1 ...
                                             tmp(3) tmp(3) tmp(3) tmp(3) tmp(3) tmp(3) tmp(3) tmp(3) tmp(3) ...
                                             tmp(3)+1 tmp(3)+1 tmp(3)+1 tmp(3)+1 tmp(3)+1 tmp(3)+1 tmp(3)+1 tmp(3)+1 tmp(3)+1],...
                                            0);
                               Pmin = min(sample)*100;
                               Pmax = max(sample)*100;
                               Ploc =  sample(14)*100;
