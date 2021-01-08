% Makes an explicit mask for input data, smoothed and unsmoothed but unmodulated
function make_VBM_explicit_mask(scan_paths, path_to_template_6, prefix)

split_stem = regexp(scan_paths, '/', 'split');

oldfilenames = cell(length(scan_paths),1);
newfilenames = cell(length(scan_paths),1);
for i = 1:length(scan_paths)
    oldfilenames{i}= cellstr(['/' fullfile(split_stem{i}{1:end-1}) '/mwc1_ns_' split_stem{i}{end}]);
    newfilenames{i}= cellstr(['/' fullfile(split_stem{i}{1:end-1}) '/s_0_mwc1_ns_' split_stem{i}{end}]);
    try
        movefile(char(oldfilenames{i}),char(newfilenames{i}));
    catch
    end
end

addpath('/group/language/data/thomascope/SD_Wordending/Masking/')
%  make_majority_mask(thresholds, consensus, outputfilname, files)
for i = 1:length(scan_paths)
    newfilenames{i}= char(newfilenames{i});
end
make_majority_mask([0.2 0.1 0.05 0.001], 0.8, [prefix '_majority_unsmoothed_mask_c1'], char(newfilenames))

