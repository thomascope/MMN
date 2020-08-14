function se_TabList

global st;
global MAP;
global SPM;
global xSPM;
global CLUSTER;



MapName = spm_select(1,'mat',['Select Map'],[],spm('Dir','se_anatomy'),'MPM.',1);
se_getMap('anat',MapName);

TxtName = spm_select(1,'any',['Select coordinate list'],[],pwd,'.txt',1);


try
    [X Y Z] = textread(TxtName,'%f %f %f');
    Orig = spm_input('Coordinate space','!+1','b','Anatomical|MNI',[1 2],2,'Select mode');
    if Orig == 1
        XYZmm = [X Y Z]';
    else
        XYZmm = [X Y-4 Z+5]';
    end
    descDa = 0;
    
catch
    
    try
        [X Y Z description] = textread(TxtName,'%f %f %f %s');
        Orig = spm_input('Coordinate space','!+1','b','Anatomical|MNI',[1 2],2,'Select mode');
        if Orig == 1
            XYZmm = [X Y Z]';
        else
            XYZmm = [X Y-4 Z+5]';
        end
        descDa = 1;
        
    catch
        spm('alert"',{'Only text files with three numeric entries per line can be read'});
        return
        
    end
end



ResName = spm_input('Output filename',1,'s','');
if ~any(ResName == '.')
    ResName = [ResName '.txt'];
end

try
    if strcmp(computer,'PCWIN')
        fid = fopen(ResName,'wt');
    else
        fid = fopen(ResName,'w+');
    end
catch
    spm('alert"',{'Could not open output file'});
end

xyz = inv(MAP(1).MaxMap.mat) * [XYZmm; ones(1,size(XYZmm,2))] ;


for PM = 1:size(MAP,2)
    ProbMax(1:size(xyz,2),PM+1) = spm_sample_vol(MAP(PM).PMap,xyz(1,:),xyz(2,:),xyz(3,:),0)' * 100;
end
ProbMax(:,1) = spm_sample_vol(MAP(1).MaxMap,xyz(1,:),xyz(2,:),xyz(3,:),0)'; ProbMax(:,1) = ProbMax(:,1) .* (ProbMax(:,1) > 99);



for indxx  = 1:size(XYZmm,2)
    if descDa
        fprintf(fid,'%s \t %s \t \n',['(' int2str(indxx) ')'],description{indxx});
    else
        fprintf(fid,'%s \t \n',['(' int2str(indxx) ')']);
    end
    
    fprintf(fid,'%s \t %3.0f\t %3.0f\t %3.0f ',...
        'X / Y / Z = ',...
        XYZmm(1,indxx),XYZmm(2,indxx),XYZmm(3,indxx));

    if Orig == 2
        fprintf(fid,'\t\t%s \t %3.0f\t %3.0f\t %3.0f',...
            'MNI:',...
            XYZmm(1,indxx),XYZmm(2,indxx)+4,XYZmm(3,indxx)-5);
    end

    ML = round(spm_sample_vol(MAP(1).Macro,xyz(1,indxx),xyz(2,indxx),xyz(3,indxx),0)');
    if ML > 0;
        MLl = MAP(1).MLabels.Labels{ML};
        fprintf(fid,'\t\t%s\n',['(' MLl ')']);
    end
    %        fprintf(fid,'\n');

    if any(ProbMax(indxx,:))
        Probs = find(ProbMax(indxx,2:end)>0); [value sortP]= sort(ProbMax(indxx,Probs+1));
        for getPr = size(Probs,2):-1:1
            [Ploc, Pmin, Pmax] = MinMax(MAP(Probs(sortP(getPr))).PMap,xyz(:,indxx));
            if getPr == size(Probs,2)
                woAssign = repmat(' ',10,1);
                if ProbMax(indxx,1)
                    woAssign(1:length(MAP(Probs(sortP(getPr))).name)) = MAP(Probs(sortP(getPr))).name;
                    fprintf(fid,'\t %s\t %s\t \n',...
                        '-> Assigned to ', woAssign);
                else
                    fprintf(fid,'\n');
                end
            end
            if (ProbMax(indxx,1+Probs(sortP(getPr))))>0
                woAssign = repmat(' ',1,15);
                woAssign(1:length(MAP(Probs(sortP(getPr))).name)) = MAP(Probs(sortP(getPr))).name;
                fprintf(fid,['%s\t %s\t %2.1f\t %s\t %2.1f\t %2.1f\t %s \n'],...
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


status = fclose(fid)
spm('alert"',{'Results saved to file'; spm_str_manip(ResName,'t'); 'in path'; pwd});





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
    0)*100;
Pmin = min(sample);
Pmax = max(sample);
Ploc =  sample(14);