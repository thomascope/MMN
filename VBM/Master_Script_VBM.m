%function Master_Script_VBM

rmpath(genpath('/imaging/local/software/spm_cbu_svn/releases/spm12_latest/'))
%addpath /imaging/local/software/spm_cbu_svn/releases/spm12_fil_r6906
addpath /group/language/data/thomascope/spm12_fil_r6906/
% addpath('/group/language/data/thomascope/SD_Wordending/') % For mask scripts
addpath('/group/language/data/thomascope/7T_full_paradigm_pilot_analysis_scripts') % For VBM scripts

groupstodo = {'matched_HCs' 'pca' 'bvFTD' 'pnfa'};
dirnames_inv = {'matched_HCs' 'pca' 'bvftd' 'vespa'};
all_subjs = [];

pathstem_structurals = '/imaging/hp02/pnfa_mmn/preprocessed/For_Thomas_dvts_sep/mri_scans/';
outdir = '/imaging/tc02/Holly_MMN/VBM/';
workingdir = pwd;

%% Which stages to run

segment_this = 0;
tiv_this = 0;
make_dartel_this = 0;
normalise_this = 0;
make_average_brain_this = 1;
stats_this = 1;


%% Data needed to work out ages for covariate - use lookup from Holly's spreadsheet
meg_numbers = {
    '13_0220'
    '13_0225'
    '13_0236'
    '13_0277'
    '13_0284'
    '13_0300'
    '13_0303'
    '13_0324'
    '13_0356'
    '13_0437'
    '13_0526'
    '14_0061'
    '14_0199'
    '14_0287'
    '14_0327'
    '14_0333'
    '11_0267'
    '12_0031'
    '12_0033'
    '12_0036'
    '12_0040'
    '12_0056'
    '12_0064'
    '12_0075'
    '12_0520'
    '14_0094'
    '14_0107'
    '08_0247'
    '08_0255'
    '08_0273'
    '08_0274'
    '08_0357'
    '08_359'
    '08_0365'
    '08_0368'
    '08_0378'
    '08_0379'
    '08_384'
    '09_0087'
    '11_0179'
    '11_0249'
    '11_0270'
    '11_0238'
    '12_0072'
    '12_0060'
    '12_0092'
    '12_0161'
    '12_0228'
    '12_0389'
    '12_0366'
    '12_0504'
    '12_0496'
    '12_0519'
    '13_0016'
    '13_0162'
    '13_0279'
    '13_0315'
    '13_0410'
    '12_0314'
    '13_0454'
    '14_0143'
    '14_0315'
    '14_0559'
    '15_0061'
    'vp1'
    'vp2'
    'vp4'
    'vp5'
    'vp6'
    'vp7'
    'vp8'
    'vp9'
    'vp10'
    'vp11'
    'vp12'
    'vc4'
    '410179'
    '510639'
    '520065'
    '420261'
    '412021'
    '08_0368'
    '12_0040'
    '520745'
    '520253'
    '08_0379'
    '420198'
    'vc3'
    '520097'
    '610496'
    '520127'
    '510629'
    '420383'
    '410084'
    '14_0094'
    '510648'
    '420100'
    '420157'
    '11_0267'
    '710446'
    '420486'
    '410390'
    '12_0031'
    '610405'
    '510355'
    '09_0087'
    '08_0365'
    '12_0520'
    '12_0036'
    '410284'
    '420435'
    '510259'
    '512003'
    '08_0378'
    '620572'
    '420162'
    '520055'
    '620592'
    '08_0273'
    '610061'
    '610288'
    '721291'
    '610932'
    
    };


all_ages = [
    68
    60
    59
    57
    63
    53
    55
    64
    65
    59
    63
    87
    59
    58
    55
    74
    58
    55
    47
    63
    63
    45
    52
    45
    57
    61
    67
    68
    67
    63
    60
    52
    66
    58
    57
    70
    61
    72
    60
    61
    61
    50
    68
    62
    59
    63
    59
    60
    64
    58
    55
    62
    78
    56
    60
    63
    63
    64
    62
    61
    59
    60
    64
    73
    80
    64
    76
    63
    75
    63
    79
    72
    82
    78
    70
    61
    58
    64
    62
    56
    56
    57
    63
    64
    61
    61
    59
    60
    65
    73
    68
    61
    50
    59
    61
    63
    59
    60
    58
    88
    59
    58
    55
    74
    68
    60
    58
    57
    63
    53
    55
    64
    65
    70
    80
    54
    63
    75
    63
    79
    72
    82
    78
    ];

assert(numel(meg_numbers)==numel(all_ages),'There is a problem with the lookup tables not being the same length')

%% First organise files
cd(pathstem_structurals)
filenames = cell(1,length(dirnames_inv));
these_ages = cell(1,length(dirnames_inv));
all_mrilist = {};
all_agelist = [];
for i = 1:length(dirnames_inv)
    cd(dirnames_inv{i});
    thesedirs = dir;
    for j = 3:length(thesedirs) %Assume that first two entries are . and ..
        if exist(thesedirs(j).name,'dir')
            cd(thesedirs(j).name)
            thesefiles = dir('*.nii');
            if isempty(thesefiles) || strcmp(thesefiles(1).name,'avg152T1.nii')
                disp(['No structural found in directory ' thesedirs(j).name])
                
            elseif length(thesefiles) == 1
                disp(['Structural filename ' thesefiles(1).name ' found in ' thesedirs(j).name])
                filenames{i}{end+1} = [pathstem_structurals dirnames_inv{i} '/' thesedirs(j).name '/' thesefiles(1).name];
                this_id = strsplit(thesedirs(j).name,'meg');
                this_lookup_loc = strncmp(this_id{2},meg_numbers,min(length(this_id{2}),7));
                if sum(this_lookup_loc) == 1
                    these_ages{i}(end+1) = all_ages(this_lookup_loc);
                elseif sum(this_lookup_loc) > 1
                    if range(all_ages(this_lookup_loc)) == 0
                        all_these_ages = all_ages(this_lookup_loc);
                        these_ages{i}(end+1) = all_these_ages(1);
                    else
                        disp(['ERROR, MORE THAN ONE AGE MATCH FOUND FOR ' thesedirs(j).name])
                    end
                else
                    disp(['Multiple images found in ' thesedirs(j).name ', using ' thesefiles(loc_shortest_filename).name])
                    filenames{i}{end+1} = [pathstem_structurals dirnames_inv{i} '/' thesedirs(j).name '/' thesefiles(loc_shortest_filename).name];
                    this_id = strsplit(thesedirs(j).name,'meg'); % For those where unique ID is a MEG number
                    this_lookup_loc = strncmp(this_id{2},meg_numbers,min(length(this_id{2}),7));
                    if sum(this_lookup_loc) == 1
                        these_ages{i}(end+1) = all_ages(this_lookup_loc);
                    else
                        this_id = strsplit(thesedirs(j).name,'cc'); % For those where unique ID is a camcan number
                        if numel(this_id) == 2
                            this_lookup_loc = strncmp(this_id{2},meg_numbers,min(length(this_id{2}),7));
                        else
                            this_id = strsplit(thesedirs(j).name,'_'); % For those where unique ID is a VESPA number
                            this_lookup_loc = strncmp(this_id{end},meg_numbers,length(this_id{end}));
                        end
                        if sum(this_lookup_loc) == 1
                            these_ages{i}(end+1) = all_ages(this_lookup_loc);
                        elseif sum(this_lookup_loc) > 1
                            if range(all_ages(this_lookup_loc)) == 0
                                all_these_ages = all_ages(this_lookup_loc);
                                these_ages{i}(end+1) = all_these_ages(1);
                            else
                                disp(['ERROR, MORE THAN ONE AGE MATCH FOUND FOR ' thesedirs(j).name])
                            end
                        else
                            disp(['ERROR, NO AGE MATCH FOUND FOR ' thesedirs(j).name])
                        end
                    end
                end
            else
                all_filelengths = zeros(1,length(thesefiles));
                for k = 1:length(thesefiles)
                    if strncmp(thesefiles(k).name,'Template',8) %Exclude DARTEL template scans
                        all_filelengths(k) = 999;
                    else
                        all_filelengths(k) = length(thesefiles(k).name);
                    end
                end
                [~, loc_shortest_filename] = min(all_filelengths);
                if strncmp(thesefiles(loc_shortest_filename).name,'vc',2)
                    disp(['VESPA control in ' thesedirs(j).name ', moving on'])
                else
                    disp(['Multiple images found in ' thesedirs(j).name ', using ' thesefiles(loc_shortest_filename).name])
                    filenames{i}{end+1} = [pathstem_structurals dirnames_inv{i} '/' thesedirs(j).name '/' thesefiles(loc_shortest_filename).name];
                    this_id = strsplit(thesedirs(j).name,'meg'); % For those where unique ID is a MEG number
                    this_lookup_loc = strncmp(this_id{2},meg_numbers,min(length(this_id{2}),7));
                    if sum(this_lookup_loc) == 1
                        these_ages{i}(end+1) = all_ages(this_lookup_loc);
                    elseif sum(this_lookup_loc) > 1
                        if range(all_ages(this_lookup_loc)) == 0
                            all_these_ages = all_ages(this_lookup_loc);
                            these_ages{i}(end+1) = all_these_ages(1);
                        else
                            disp(['ERROR, MORE THAN ONE AGE MATCH FOUND FOR ' thesedirs(j).name])
                        end
                    else
                        this_id = strsplit(thesedirs(j).name,'cc'); % For those where unique ID is a camcan number
                        if numel(this_id) == 2
                            this_lookup_loc = strncmp(this_id{2},meg_numbers,min(length(this_id{2}),7));
                        else
                            this_id = strsplit(thesedirs(j).name,'_'); % For those where unique ID is a VESPA number
                            this_lookup_loc = strcmp(this_id{end},meg_numbers);
                        end
                        if sum(this_lookup_loc) == 1
                            these_ages{i}(end+1) = all_ages(this_lookup_loc);
                        elseif sum(this_lookup_loc) > 1
                            if range(all_ages(this_lookup_loc)) == 0
                                all_these_ages = all_ages(this_lookup_loc);
                                these_ages{i}(end+1) = all_these_ages(1);
                            else
                                disp(['ERROR, MORE THAN ONE AGE MATCH FOUND FOR ' thesedirs(j).name])
                            end
                        else
                            disp(['ERROR, NO AGE MATCH FOUND FOR ' thesedirs(j).name])
                        end
                    end
                end
            end
            
            cd ..
        end
    end
    cd(pathstem_structurals)
    
    assert(numel(filenames{i})==numel(these_ages{i}),['There is a problem with the number of ages and scans not being the same length for ' dirnames_inv{i}])
    disp([num2str(numel(filenames{i})) ' scans found for ' dirnames_inv{i} ' each with a matching age in lookup table '])
    all_mrilist = [all_mrilist filenames{i}];
    all_agelist = [all_agelist these_ages{i}];
end
cd(workingdir)


%% Next segment all the images

if segment_this == 1
    nrun = length(all_mrilist);
    jobfile = {'/imaging/tc02/vespa/scans/PNFA_VBM/tom/VBM_batch_segment.m'};
    jobs = repmat(jobfile, 1, nrun);
    inputs = cell(1, nrun);
    
    for crun = 1:nrun
        inputs{1, crun} = cellstr([all_mrilist{crun}]); % for dartel templating
        debug_split = strsplit(all_mrilist{crun},'/'); %Debug step
        failedhere(crun) = isempty(dir(['/' fullfile(debug_split{1:end-1}) '/mwc3*']));
    end
    
    segmentworkedcorrectly = zeros(1,nrun);
    jobs = repmat(jobfile, 1, 1);
    
    if isempty(gcp('nocreate'))
        cbupool(nrun)
    end
    
    parfor crun = 1:nrun
        spm('defaults', 'PET');
        spm_jobman('initcfg')
        try
            spm_jobman('run', jobs, inputs{:,crun});
            segmentworkedcorrectly(crun) = 1;
        catch
            segmentworkedcorrectly(crun) = 0;
        end
    end
end

%%  Now calculate the TIV and then rename all of the mc files as will be overwritten by DARTEL (if not smoothed)
if tiv_this == 1
    nrun = 1;
    jobfile = {'/imaging/tc02/vespa/scans/PNFA_VBM/tom/VBM_batch_TIV.m'};
    jobs = repmat(jobfile, 1, nrun);
    inputs = cell(2, nrun);
    
    for crun = 1:nrun
        inputs{1, crun} = cell(length(all_mrilist),1);
        for i = 1:length(all_mrilist)
            inputs{1, crun}(i) = cellstr([all_mrilist{i}(1:end-4) '_seg8.mat']); % for dartel templating
        end
    end
    
    inputs{2,1} = [outdir 'total_intractranial_volumes_VBM'];
    tiv_filename = [inputs{2,1} '.csv'];
    
    TIVworkedcorrectly = zeros(1,nrun);
    jobs = repmat(jobfile, 1, 1);
    
    parfor crun = 1:nrun %Even though single shot, submit to parfor to avoid overloading login node
        spm('defaults', 'PET');
        spm_jobman('initcfg')
        try
            spm_jobman('run', jobs, inputs{:,crun});
            TIVworkedcorrectly(crun) = 1;
        catch
            TIVworkedcorrectly(crun) = 0;
        end
    end
    
    split_stem = regexp(all_mrilist, '/', 'split');
    old_imagepaths = cell(length(all_mrilist),3);
    new_imagepaths = cell(length(all_mrilist),3);
    for i = 1:length(all_mrilist)
        for j = 1:3
            old_imagepaths(i,j) = cellstr(['/' fullfile(split_stem{i}{1:end-1}) '/mwc' num2str(j) split_stem{i}{end}]);
            new_imagepaths(i,j) = cellstr(['/' fullfile(split_stem{i}{1:end-1}) '/mwc' num2str(j) '_ns_' split_stem{i}{end}]);
            movefile(char(old_imagepaths(i,j)),char(new_imagepaths(i,j)))
        end
    end
end

%% Then make a DARTEL template based on all images (impossible to age match and do two-stage approach)
if make_dartel_this == 1
    
    nrun = 1; % enter the number of runs here
    jobfile = {'/imaging/tc02/vespa/scans/PNFA_VBM/tom/VBM_batch_dartel.m'};
    jobs = repmat(jobfile, 1, nrun);
    inputs = cell(2,1);
    inputs{1,1} = cell(length(all_mrilist),1);
    inputs{2,1} = cell(length(all_mrilist),1);
    split_stem = regexp(all_mrilist, '/', 'split');
    
    for i = 1:length(all_mrilist)
        inputs{1,1}(i) = cellstr(['/' fullfile(split_stem{i}{1:end-1}) '/rc1' split_stem{i}{end}]);
        inputs{2,1}(i) = cellstr(['/' fullfile(split_stem{i}{1:end-1}) '/rc2' split_stem{i}{end}]);
    end
    
    dartelworkedcorrectly = zeros(1,nrun);
    jobs = repmat(jobfile, 1, 1);
    
    parfor crun = 1:nrun %Still submit as a parfor to avoid overloading a login node
        spm('defaults', 'PET');
        spm_jobman('initcfg')
        try
            spm_jobman('run', jobs, inputs{:,crun});
            dartelworkedcorrectly(crun) = 1;
        catch
            dartelworkedcorrectly(crun) = 0;
        end
    end
end

%% Now normalise all scans

if normalise_this == 1
    nrun = length(all_mrilist);
    jobfile = {'/imaging/tc02/vespa/scans/PNFA_VBM/tom/VBM_batch_normalise.m'};
    jobs = repmat(jobfile, 1, nrun);
    inputs = cell(3, nrun);
    
    split_stem = regexp(all_mrilist, '/', 'split');
    split_stem_template = regexp(all_mrilist, '/', 'split');
    path_to_template_6 = cellstr(['/' fullfile(split_stem_template{1}{1:end-1}) '/Template_6.nii']);
    
    for crun = 1:nrun
        inputs{1, crun} = path_to_template_6; % Normalise to MNI Space: Dartel Template - cfg_files
        inputs{2, crun} = cellstr(['/' fullfile(split_stem{crun}{1:end-1}) '/u_rc1' split_stem{crun}{end}(1:end-4) '_Template.nii']);
        if ~exist(char(inputs{2,crun}),'file')
            inputs{2, crun} = cellstr(['/' fullfile(split_stem{crun}{1:end-1}) '/u_rc1' split_stem{crun}{end}]); %Deal with the difference in naming convention depending if part of the dartel templating or not
        end
        inputs{3, crun} = cellstr(['/' fullfile(split_stem{crun}{1:end-1}) '/c1' split_stem{crun}{end}]);
    end
    
    normaliseworkedcorrectly = zeros(1,nrun);
    jobs = repmat(jobfile, 1, 1);
    
    parfor crun = 1:nrun
        spm('defaults', 'PET');
        spm_jobman('initcfg')
        try
            spm_jobman('run', jobs, inputs{:,crun});
            normaliseworkedcorrectly(crun) = 1;
        catch
            normaliseworkedcorrectly(crun) = 0;
        end
    end
end

%% Now read in TIV file and do group stats with TIV and age file as covariates in the ANOVA
if stats_this ==1 ;
    if exist(tiv_filename)
        filename =tiv_filename;
    else
        error('no tiv_file found')
    end
    delimiter = ',';
    startRow = 2;
    endRow = inf;
    formatSpec = '%s%f%f%f%[^\n\r]';
    fileID = fopen(filename,'r');
    dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'HeaderLines', startRow(1)-1, 'ReturnOnError', false);
    fclose(fileID);
    tiv= dataArray{2}+dataArray{3}+dataArray{4};
    
    nrun = 1;
    jobfile = {'/imaging/tc02/vespa/scans/PNFA_VBM/tom/VBM_batch_multigroup_TIV_age.m'};
    jobs = repmat(jobfile, 1, nrun);
    inputs = cell(8, nrun);
    
    stats_folder = {[outdir 'VBM_stats/factorial_full_group_vbm_TIVnormalised_agecovaried_unsmoothedmask']};
    
    split_stem = cell(1,4);
    for i = 1:4
        split_stem{i} = regexp(filenames{i}, '/', 'split');
    end
    
    inputs{1, 1} = stats_folder;
    
    for crun = 1:nrun
        for i = 1:4
            inputs{1+i, crun} = cell(length(filenames{i}),1);
            for j = 1:length(filenames{i})
                inputs{1+i,crun}(j) = cellstr(['/' fullfile(split_stem{i}{j}{1:end-1}) '/smwc1' split_stem{i}{j}{end}]);
            end
        end
    end
    
    try
        inputs{6, 1} = tiv;
    catch
        filename = tiv_filename;
        delimiter = ',';
        startRow = 2;
        endRow = inf;
        formatSpec = '%s%f%f%f%[^\n\r]';
        fileID = fopen(filename,'r');
        dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'HeaderLines', startRow(1)-1, 'ReturnOnError', false);
        fclose(fileID);
        tiv= dataArray{2}+dataArray{3}+dataArray{4};
        inputs{6, 1} = tiv;
    end
    inputs{7, 1} = all_agelist;
    inputs{8, 1} = {'control_majority_unsmoothed_mask_c1_thr0.05_cons0.8.img'};
    
    if ~exist(char(inputs{8, 1}),'file')
        split_stem_template = regexp(all_mrilist, '/', 'split');
        path_to_template_6 = cellstr(['/' fullfile(split_stem_template{1}{1:end-1}) '/Template_6.nii']);
        make_VBM_explicit_mask(filenames{1}, path_to_template_6, 'control')
    end
    
    spm_jobman('run', jobs, inputs{:});
    
    inputs = cell(1, nrun);
    inputs{1, 1} =  {[char(stats_folder) '/SPM.mat']};
    
    jobfile = {'/imaging/tc02/vespa/scans/PNFA_VBM/tom/VBM_batch_estimate.m'};
    jobs = repmat(jobfile, 1, nrun);
    
    spm_jobman('run', jobs, inputs{:});
    
    jobfile = {'/imaging/tc02/vespa/scans/PNFA_VBM/tom/VBM_batch_contrast_MMN.m'};
    jobs = repmat(jobfile, 1, nrun);
    
    spm_jobman('run', jobs, inputs{:});
    
    %     jobfile = {'/imaging/tc02/vespa/scans/PNFA_VBM/tom/VBM_batch_results.m'};
    %     jobs = repmat(jobfile, 1, nrun);
    %
    %     spm_jobman('run', jobs, inputs{:});
    
end

%% Make average brain by normalising MPRAGES to template and then averaging them
if make_average_brain_this == 1
    
    nrun = length(all_mrilist);
    jobfile = {'/imaging/tc02/vespa/scans/PNFA_VBM/tom/VBM_batch_normalise_unmodulated_unsmoothed.m'};
    jobs = repmat(jobfile, 1, nrun);
    inputs = cell(3, nrun);
    
    split_stem = regexp(all_mrilist, '/', 'split');
    split_stem_template = regexp(all_mrilist, '/', 'split');
    path_to_template_6 = cellstr(['/' fullfile(split_stem_template{1}{1:end-1}) '/Template_6.nii']);
    
    for crun = 1:nrun
        inputs{1, crun} = path_to_template_6; % Normalise to MNI Space: Dartel Template - cfg_files
        inputs{2, crun} = cellstr(['/' fullfile(split_stem{crun}{1:end-1}) '/u_rc1' split_stem{crun}{end}(1:end-4) '_Template.nii']);
        if ~exist(char(inputs{2,crun}),'file')
            error(['Cannot find template at /' fullfile(split_stem{crun}{1:end-1}) '/u_rc1' split_stem{crun}{end}(1:end-4) '_Template.nii']);
        end
        inputs{3, crun} = cellstr(['/' fullfile(split_stem{crun}{1:end-1}) '/' split_stem{crun}{end}]);
    end
    
    normaliseforaverageworkedcorrectly = zeros(1,nrun);
    jobs = repmat(jobfile, 1, 1);
    
    parfor crun = 1:nrun
        spm('defaults', 'PET');
        spm_jobman('initcfg')
        try
            spm_jobman('run', jobs, inputs{:,crun});
            normaliseforaverageworkedcorrectly(crun) = 1;
        catch
            normaliseforaverageworkedcorrectly(crun) = 0;
        end
    end
    
    nrun = 1;
    jobfile = {'/imaging/tc02/vespa/scans/PNFA_VBM/tom/VBM_batch_imcalc_average.m'};
    inputs = cell(2, nrun);
    
    split_stem = regexp(all_mrilist, '/', 'split');
    inputs{1,1} = cell(length(all_mrilist),1);
    
    for i = 1:length(all_mrilist)
        inputs{1,1}(i) = cellstr(['/' fullfile(split_stem{i}{1:end-1}) '/w' split_stem{i}{end}]);
    end
    
    inputs{2,1} = ['average_T1head'];
    
    imcalcworkedcorrectly = zeros(1,nrun);
    jobs = repmat(jobfile, 1, 1);
    
    parfor crun = 1:nrun %Still submit as a parfor to avoid overloading a login node
        spm('defaults', 'PET');
        spm_jobman('initcfg')
        try
            spm_jobman('run', jobs, inputs{:,crun});
            imcalcworkedcorrectly(crun) = 1;
        catch
            imcalcworkedcorrectly(crun) = 0;
        end
    end
end

%pause