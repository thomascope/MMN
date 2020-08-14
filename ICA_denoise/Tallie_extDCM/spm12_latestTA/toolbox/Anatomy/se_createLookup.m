clear,clc

MapName = spm_select(1,'mat',['Select Map'],[],pwd,'MPM.',1); 

load(MapName);

if strcmp(computer,'PCWIN')
    fid = fopen([spm_str_manip(MapName,'rt') '.txt'],'wt');
else
    fid = fopen([spm_str_manip(MapName,'rt') '.txt'],'w+');
end

for i=1:size(MAP,2)
    name{i} = MAP(i).name;
    GV(i)   = MAP(i).GV;
    nr(i)   = i;
end

[names idx] = sort(name);
GVs         = GV(idx);
nrs         = nr(idx);


fprintf(fid,[repmat('%s\t',1,23) '\n'],'Area','Index','GV',...
                   'MPM Vol-L', 'MPM Vol-R', 'MPM CoG-L', '','', 'MPM CoG-R', '','', 'MPM mean-L', 'MPM mean-R',...
                   'PMap Vol-L','PMap Vol-R','PMap CoG-L','','','PMap CoG-R', '','','PMap mean-L','PMap mean-R');


for i=1:size(MAP,2)
    namen = repmat(' ',1,15);
    namen(1:length(names{i})) = names{i};
    index = nrs(i);
    
    MAP(index).allXYZmm = MAP(1).MaxMap.mat * [MAP(index).allXYZ; ones(1,size(MAP(index).allXYZ,2))];
    MAP(index).allXYZmm = MAP(index).allXYZmm(1:3,:);
    
    fprintf(fid,[repmat('%s\t',1,3) repmat('%s\t',1,2) repmat('%3.1f\t',1,8) repmat('%s\t',1,2)  repmat('%3.1f\t',1,8) '\n'],...
       namen,int2str(index),int2str(MAP(index).GV),...
       int2str(MAP(index).VOL(1)),int2str(MAP(index).VOL(2)),...
       mean(MAP(index).XYZmm(:,MAP(index).LR == -1),2),mean(MAP(index).XYZmm(:,MAP(index).LR == 1),2),...
       mean(MAP(index).Z(:,MAP(index).LR == -1),2)*100,mean(MAP(index).Z(:,MAP(index).LR == 1),2)*100,...
       int2str(sum(MAP(index).allZ(:,MAP(index).allLR == -1))),int2str(sum(MAP(index).allZ(:,MAP(index).allLR == 1))),...
       mean(MAP(index).allXYZmm(:,MAP(index).allLR == -1),2),mean(MAP(index).allXYZmm(:,MAP(index).allLR == 1),2),...
       mean(MAP(index).allZ(:,MAP(index).allLR == -1),2)*100,mean(MAP(index).allZ(:,MAP(index).allLR == 1),2)*100);
       
    
     if round(i/3) == i/3
         fprintf(fid,'\n');
     end
end
status = fclose(fid);