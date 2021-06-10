Overlap_fnames = {'./VBM_bycond/Positive_Correlates_001unc.nii';
    './VBM/all.nii';
    './Language_overlap_n220.img';
    './MD_overlap_n63.img'};
spm_reslice(Overlap_fnames)
for i = 1:length(Overlap_fnames)
    [PATHSTR,NAME,EXT] = fileparts(Overlap_fnames{i});
    Overlap_fnames{i} = [PATHSTR filesep 'r' NAME EXT];
end
all_fstructs = spm_vol(Overlap_fnames);

%First calculate proportion of voxels that overlap with the probablistic maps thresholded at 5%
disp('Image 1 is VBM MMN correlates, 0.001unc, Image 2 is Language probablistic map 5% thresholded')
[totvox1,totvox2,totvox3] = propvox(all_fstructs{1}, 0.1, all_fstructs{3}, 0.05); %NB: DO NOT THRESHOLD AT ZERO because of interpolation in the reslicing.

disp('Image 1 is VBM MMN correlates, 0.001unc, Image 2 is Multiple Demand probablistic map 5% thresholded')
[totvox1,totvox2,totvox3] = propvox(all_fstructs{1}, 0.1, all_fstructs{4}, 0.05);

%Now check specificiity of this, by looking at overall atrophy map thresholded FWE 0.05, (t = 4.75)
disp('Image 1 is VBM overall, FWE 0.05, Image 2 is Language probablistic map 5% thresholded')
[totvox1,totvox2,totvox3] = propvox(all_fstructs{2}, 4.75, all_fstructs{3}, 0.05);

disp('Image 1 is VBM oveall, FWE 0.05, Image 2 is Multiple Demand probablistic map 5% thresholded')
[totvox1,totvox2,totvox3] = propvox(all_fstructs{2}, 4.75, all_fstructs{4}, 0.05);

%Now do weighted multiplications to get overall association strength scores
mask = './control_majority_unsmoothed_mask_c1_thr0.05_cons0.8.img';

disp('Image 1 is VBM MMN correlates, 0.001unc, Image 2 is Language probablistic map 5% thresholded')
[total1,total2,total3] = weightedvox(all_fstructs{1}, all_fstructs{3}, mask); %NB: DO NOT THRESHOLD AT ZERO because of interpolation in the reslicing.

disp('Image 1 is VBM MMN correlates, 0.001unc, Image 2 is Multiple Demand probablistic map 5% thresholded')
[total1,total2,total3] = weightedvox(all_fstructs{1}, all_fstructs{4}, mask);

% Again, control analysis with VBM
disp('Image 1 is VBM overall, FWE 0.05, Image 2 is Language probablistic map 5% thresholded')
[total1,total2,total3] = weightedvox(all_fstructs{2}, all_fstructs{3}, mask); %NB: DO NOT THRESHOLD AT ZERO because of interpolation in the reslicing.

disp('Image 1 is VBM overall, FWE 0.05, Image 2 is Multiple Demand probablistic map 5% thresholded')
[total1,total2,total3] = weightedvox(all_fstructs{2}, all_fstructs{4}, mask);


