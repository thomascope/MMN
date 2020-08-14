function se_tabelate(cl,fid)

global st;
global MAP;
global SPM;
global xSPM;
global CLUSTER;

GV = [MAP.GV];





% Volume
if isfield(SPM,'STAT')
    try
        titl = ['Cluster ' int2str(cl) ' (' int2str(CLUSTER(cl).voxel(1)) ' vox)'...
            ': ' SPM.title ' '  SPM.thresDef];
    catch
        titl = ['Cluster ' int2str(cl) ' (' int2str(CLUSTER(cl).voxel(1)) ' vox)'...
            ': ' SPM.title '  ( ' SPM.STAT '> ' sprintf('%4.2f', SPM.u) ')'];
    end
else
    titl = SPM.title;
    if length(titl)>15; titl = titl(1:15); end
    titl = ['Cluster ' int2str(cl) ' (' int2str(CLUSTER(cl).voxel(1)) ' vox)'...
        ': ' SPM.title '  (u > ' sprintf('%4.2f', SPM.u) ', k > ' int2str(SPM.k)  ')'];
    SPM.STAT = '';
end

fprintf(fid,'%s\n\n',titl);
if ~isfield(CLUSTER,'gtmp') || prod(size(CLUSTER(cl).gtmp)) < 1
    tmp   = CLUSTER(cl).XYZp;
    tms = sign(CLUSTER(cl).XYZmm);
    clear gtmp; gZ = []; mult = 1;
    for xbx = 0:xSPM.xVOL.VOX(1)-1
        for ybx = 0:xSPM.xVOL.VOX(2)-1
            for zbx = 0:xSPM.xVOL.VOX(3)-1
                if exist('gtmp')
                    gtmp = [gtmp(1,:) tmp(1,:)-(tms(1,:)*xbx);...
                        gtmp(2,:) tmp(2,:)-(tms(2,:)*ybx);...
                        gtmp(3,:) tmp(3,:)-(tms(3,:)*zbx)];
                    mult = mult+1;
                else
                    gtmp = [tmp(1,:)-(tms(1,:)*xbx);...
                        tmp(2,:)-(tms(2,:)*ybx);...
                        tmp(3,:)-(tms(3,:)*zbx)];
                end
            end
        end
    end
    CLUSTER(cl).gtmp = gtmp;
    CLUSTER(cl).mult = mult;
end

Msk    = round(spm_sample_vol(MAP(1).AnatMask,CLUSTER(cl).gtmp(1,:),CLUSTER(cl).gtmp(2,:),CLUSTER(cl).gtmp(3,:),0));
Volume = round(spm_sample_vol(MAP(1).MaxMap,  CLUSTER(cl).gtmp(1,:),CLUSTER(cl).gtmp(2,:),CLUSTER(cl).gtmp(3,:),0)); % Corresponding voxel + next towards origin in all directions

vp = [];
for v=1:numel(MAP)
    if any(Volume(Msk == 2) == GV(v))
        vp = [vp; size(find(Volume(Msk == 2) == GV(v)),2) v -1];
    end
    if any(Volume(Msk == 1) == GV(v))
        vp = [vp; size(find(Volume(Msk == 1) == GV(v)),2) v 1];
    end
end
vp = sortrows(vp);
for v=size(vp,1):-1:1
    if (100*vp(v,1)/(CLUSTER(cl).voxel(1)*CLUSTER(cl).mult))>=0.01
        if vp(v,3) == -1;
            side = 'left'; vol = MAP(vp(v,2)).VOL(1);
        else;
            side = 'right'; vol = MAP(vp(v,2)).VOL(2);
        end
        fprintf(fid,'%2.1f \t %s \t %2.1f \t      %s \t   %s \t  %s \t   %2.1f \t %s \n',...
            vp(v,1)/CLUSTER(cl).mult,...
            ' voxel = ',...
            (100*vp(v,1)/(CLUSTER(cl).voxel(1)*CLUSTER(cl).mult)),...
            '% in ',...
            side, MAP(vp(v,2)).name,...
            100*vp(v,1)/vol,...
            '% of this area activated');
    end
end

if ~any(any(vp))
    fprintf(fid,'%s\n','Volume not matching any probability map');
else
    fprintf(fid,'\n%2.1f \t %s \t %2.1f \t      %s\n',...
        (sum(vp(:,1))/CLUSTER(cl).mult),...
        ' voxel = ',...
        100*sum(vp(:,1))/(CLUSTER(cl).voxel(1)*CLUSTER(cl).mult),...
        '% of the cluster volume are assigned in total');

end

fprintf(fid,'\n\n\n\n');






% Maxima
m_anz = 0;
if (isempty(max(CLUSTER(cl).maxZ)) & min(CLUSTER(cl).Z) == max(CLUSTER(cl).Z)) | CLUSTER(cl).maxZ == 0
else

    if isempty(max(CLUSTER(cl).maxZ))
        CLUSTER(cl).maxZ = CLUSTER(cl).Z(:,find(CLUSTER(cl).Z == max(CLUSTER(cl).Z)));
        CLUSTER(cl).maxM = CLUSTER(cl).XYZ(:,find(CLUSTER(cl).Z == max(CLUSTER(cl).Z)));
    end
    [B,INDEX] = sort(CLUSTER(cl).maxZ); xyzTMP = []; xyz = [];  ProbMax = [];

    for indxx  = INDEX(size(B,2):-1:1);
        xyzTMP(1:4,indxx)   = xSPM.xVOL.M * [CLUSTER(cl).maxM(1,indxx); CLUSTER(cl).maxM(2,indxx); CLUSTER(cl).maxM(3,indxx); 1];		%-Reorientation matrix
        xyz(1:4,indxx) = inv(MAP(1).MaxMap.mat) * xyzTMP(1:4,indxx);
    end

    for PM = 1:size(MAP,2)
        ProbMax(1:size(xyz,2),PM+1) = spm_sample_vol(MAP(PM).PMap,xyz(1,:),xyz(2,:),xyz(3,:),0)';
    end
    ProbMax(:,1) = spm_sample_vol(MAP(1).MaxMap,xyz(1,:),xyz(2,:),xyz(3,:),0)'; ProbMax(:,1) = ProbMax(:,1) .* (ProbMax(:,1) > 99);


    if size(B,2)>10
        genug = size(B,2)-10;
    else
        genug = 1;
    end

    for indxx  = INDEX(size(B,2):-1:genug);
        m_anz = m_anz+1;

        if numel(SPM.STAT) == 0
            tellstat = 'STAT';
        else
            tellstat = SPM.STAT;
        end
        if m_anz < 10; pad = '0'; else; pad = ''; end
        fprintf(fid,'%s\t %s\t %s \t %3.0f\t %3.0f\t %3.0f ',...
            ['Maximum ' pad int2str(m_anz) ],...
            [tellstat ' = ' num2str(CLUSTER(cl).maxZ(indxx),'%2.2f')],...
            'X / Y / Z = ',...
            xyzTMP(1,indxx),xyzTMP(2,indxx),xyzTMP(3,indxx));

        if xSPM.xVOL.orig == 2
            fprintf(fid,'\t\t%s \t %3.0f\t %3.0f\t %3.0f',...
                'MNI:',...
                xyzTMP(1,indxx),xyzTMP(2,indxx)+4,xyzTMP(3,indxx)-5);
        end

        ML = round(spm_sample_vol(MAP(1).Macro,xyz(1,indxx),xyz(2,indxx),xyz(3,indxx),0)');
        if ML > 0;
            MLl = MAP(1).MLabels.Labels{ML};
        else
            MLl = 'N/A';
        end
        fprintf(fid,'\t\t%s',[MLl]);

        %        fprintf(fid,'\n');

        if any(ProbMax(indxx,:))
            Probs = find(ProbMax(indxx,2:end)>0); [value sortP]= sort(ProbMax(indxx,Probs+1));
            for getPr = size(Probs,2):-1:1
                [Ploc, Pmin, Pmax] = MinMax(MAP(Probs(sortP(getPr))).PMap,xyz(:,indxx));
                if getPr == size(Probs,2)

                    if spm_sample_vol(MAP(1).AnatMask,xyz(1,indxx),xyz(2,indxx),xyz(3,indxx),0)>1.5
                        side = 'left'; else; side = 'right'; end
                    woAssign = repmat(' ',10,1);
                    if ProbMax(indxx,1)
                        woAssign(1:length(MAP(Probs(sortP(getPr))).name)) = MAP(Probs(sortP(getPr))).name;
                        fprintf(fid,'\t %s\t %s\t %s\t \n',...
                            '-> Assigned to ', side, woAssign);
                    else
                        fprintf(fid,'\n');
                    end
                end
                if (ProbMax(indxx,1+Probs(sortP(getPr)))*100)>1
                    woAssign = repmat(' ',1,15);
                    woAssign(1:length(MAP(Probs(sortP(getPr))).name)) = MAP(Probs(sortP(getPr))).name;
                    fprintf(fid,['\t%s\t %s\t %2.0f\t %s\t %2.0f\t %2.0f\t %s \n'],...
                        'Probability for  ',...
                        woAssign,...
                        Ploc,...
                        '%           [', Pmin, Pmax, '%]');
                end
            end
        else
            fprintf(fid,'\n');
        end
        fprintf(fid,'\n');
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