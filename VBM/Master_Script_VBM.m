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

cd(pathstem_structurals)
filenames = cell(1,length(dirnames_inv));

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
                end
            end
            cd ..
        end
    end
    cd(pathstem_structurals)
end
cd(workingdir)
pause