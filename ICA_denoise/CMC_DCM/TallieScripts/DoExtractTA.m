function D = DoExtractTA(D,varargin)


%addpath(genpath('/home/as08/old_spm12/'));
Lpos = [-46    20     8;  % These from Garrido 2008 / Phillips 2015
        -61   -32     8;
        %-42   -14     7;
        46    20     8;
        59   -25     8;
        46   -14     8];
Sname = {'LIFG';
        'LSTG';
        %'LAud';
        'RIFG';
        'RSTG';
        'RAud'};

% specify:
m{1}.spm.meeg.source.extract.D(1) = {D};
m{1}.spm.meeg.source.extract.val = 1;
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
%m{1}.spm.meeg.source.extract.source(6).label = Sname{6};
%m{1}.spm.meeg.source.extract.source(6).xyz = Lpos(6,:);
m{1}.spm.meeg.source.extract.rad = 5;
m{1}.spm.meeg.source.extract.type = 'trials';
m{1}.spm.meeg.source.extract.fname = ['LFP' num2str(length(Sname))];

spm('defaults','EEG');
spm_jobman('run',m);

clear D
load(['LFP' num2str(length(Sname))]);
x = num2cell(Lpos(:,1));
y = num2cell(Lpos(:,2));
[D.channels.X_plot2D] = deal(x{:});
[D.channels.Y_plot2D] = deal(y{:});


end

