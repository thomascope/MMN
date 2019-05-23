function D = Fullpipeline_extraction(D,folder,outfolder, val)


Lpos = [-42, -22, 7;
        -61, -32, 8;  
        -46, 20, 8;
        -49, -38, 38;
        46, -14, 8;
        59, -25, 8;
        46, 20, 8;
        57, -38, 42];
Sname = {'left A1';
         'left STG';
         'left IFG';
         'left IPC';
         'right A1';
         'right STG';
         'right IFG';
         'right IPC'};
        
    
[f1,f2,f3] = fileparts(D);
fn = sprintf('%s/%s/%dLFP_%s%s',outfolder,folder, length(Sname), f2, f3);
[destfolder,~,~] = fileparts(fn);
if ~exist(destfolder)
    mkdir(destfolder)
end

% specify:
m{1}.spm.meeg.source.extract.D(1) = {D};
m{1}.spm.meeg.source.extract.val = val;
m{1}.spm.meeg.source.extract.source(1).label = Sname{1};
m{1}.spm.meeg.source.extract.source(1).xyz = Lpos(1,:);
m{1}.spm.meeg.source.extract.source(2).label = Sname{2};
m{1}.spm.meeg.source.extract.source(2).xyz = Lpos(2,:);
m{1}.spm.meeg.source.extract.source(3).label = Sname{3};
m{1}.spm.meeg.source.extract.source(3).xyz = Lpos(3,:);
m{1}.spm.meeg.source.extract.source(4).label = Sname{4};
m{1}.spm.meeg.source.extract.source(4).xyz = Lpos(4,:);
m{1}.spm.meeg.source.extract.source(5).label = Sname{5};
m{1}.spm.meeg.source.extract.source(5).xyz = Lpos(5,:);
m{1}.spm.meeg.source.extract.source(6).label = Sname{6};
m{1}.spm.meeg.source.extract.source(6).xyz = Lpos(6,:);
m{1}.spm.meeg.source.extract.source(7).label = Sname{7};
m{1}.spm.meeg.source.extract.source(7).xyz = Lpos(7,:);
m{1}.spm.meeg.source.extract.source(8).label = Sname{8};
m{1}.spm.meeg.source.extract.source(8).xyz = Lpos(8,:);
m{1}.spm.meeg.source.extract.rad = 5;
m{1}.spm.meeg.source.extract.type = 'trials';
m{1}.spm.meeg.source.extract.fname = fn;

spm('defaults','EEG');
spm_jobman('run',m);

clear D
load(fn);
x = num2cell(Lpos(:,1));
y = num2cell(Lpos(:,2));
[D.channels.X_plot2D] = deal(x{:});
[D.channels.Y_plot2D] = deal(y{:});
save(fn,'D')


end

