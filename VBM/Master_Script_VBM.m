function Master_Script_VBM

rmpath(genpath('/imaging/local/software/spm_cbu_svn/releases/spm12_latest/'))
%addpath /imaging/local/software/spm_cbu_svn/releases/spm12_fil_r6906
addpath /group/language/data/thomascope/spm12_fil_r6906/

groupstodo = {'matched_HCs' 'pca' 'bvFTD' 'pnfa'};
dirnames_inv = {'matched_HCs' 'pca' 'bvftd' 'vespa'};
all_subjs = [];

pathstem_structurals = '/imaging/hp02/pnfa_mmn/preprocessed/For_Thomas_dvts_sep/mri_scans/';
outdir = '/imaging/tc02/Holly_MMN/VBM/';
workingdir = pwd;

% Data needed to work out ages for covariate - use lookup from Holly's spreadsheet
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

% First organise files
cd(pathstem_structurals)
filenames = cell(1,length(dirnames_inv));
these_ages = cell(1,length(dirnames_inv));
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
                filenames{i}{end+1} = [dirnames_inv{i} '/' thesedirs(j).name '/' thesefiles(1).name];
                this_id = strsplit(thesedirs(j).name,'meg');
                this_lookup_loc = strncmp(this_id{2},meg_numbers,min(length(this_id{2}),7));
                if sum(this_lookup_loc) == 1
                    these_ages{i}{end+1} = all_ages(this_lookup_loc);
                elseif sum(this_lookup_loc) > 1
                    if range(all_ages(this_lookup_loc)) == 0
                        all_these_ages = all_ages(this_lookup_loc);
                        these_ages{i}{end+1} = all_these_ages(1);
                    else
                        disp(['ERROR, MORE THAN ONE AGE MATCH FOUND FOR ' thesedirs(j).name])
                    end
                else
                    disp(['Multiple images found in ' thesedirs(j).name ', using ' thesefiles(loc_shortest_filename).name])
                    filenames{i}{end+1} = [dirnames_inv{i} '/' thesedirs(j).name '/' thesefiles(loc_shortest_filename).name];
                    this_id = strsplit(thesedirs(j).name,'meg'); % For those where unique ID is a MEG number
                    this_lookup_loc = strncmp(this_id{2},meg_numbers,min(length(this_id{2}),7));
                    if sum(this_lookup_loc) == 1
                        these_ages{i}{end+1} = all_ages(this_lookup_loc);
                    else
                        this_id = strsplit(thesedirs(j).name,'cc'); % For those where unique ID is a camcan number
                        if numel(this_id) == 2
                            this_lookup_loc = strncmp(this_id{2},meg_numbers,min(length(this_id{2}),7));
                        else
                            this_id = strsplit(thesedirs(j).name,'_'); % For those where unique ID is a VESPA number
                            this_lookup_loc = strncmp(this_id{end},meg_numbers,length(this_id{end}));
                        end
                        if sum(this_lookup_loc) == 1
                            these_ages{i}{end+1} = all_ages(this_lookup_loc);
                        elseif sum(this_lookup_loc) > 1
                            if range(all_ages(this_lookup_loc)) == 0
                                all_these_ages = all_ages(this_lookup_loc);
                                these_ages{i}{end+1} = all_these_ages(1);
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
                    all_filelengths(k) = length(thesefiles(k).name);
                end
                [~, loc_shortest_filename] = min(all_filelengths);
                if strncmp(thesefiles(loc_shortest_filename).name,'vc',2)
                    disp(['VESPA control in ' thesedirs(j).name ', moving on'])
                else
                    disp(['Multiple images found in ' thesedirs(j).name ', using ' thesefiles(loc_shortest_filename).name])
                    filenames{i}{end+1} = [dirnames_inv{i} '/' thesedirs(j).name '/' thesefiles(loc_shortest_filename).name];
                    this_id = strsplit(thesedirs(j).name,'meg'); % For those where unique ID is a MEG number
                    this_lookup_loc = strncmp(this_id{2},meg_numbers,min(length(this_id{2}),7));
                    if sum(this_lookup_loc) == 1
                        these_ages{i}{end+1} = all_ages(this_lookup_loc);
                    elseif sum(this_lookup_loc) > 1
                        if range(all_ages(this_lookup_loc)) == 0
                            all_these_ages = all_ages(this_lookup_loc);
                            these_ages{i}{end+1} = all_these_ages(1);
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
                            these_ages{i}{end+1} = all_ages(this_lookup_loc);
                        elseif sum(this_lookup_loc) > 1
                            if range(all_ages(this_lookup_loc)) == 0
                                all_these_ages = all_ages(this_lookup_loc);
                                these_ages{i}{end+1} = all_these_ages(1);
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
end
cd(workingdir)



pause