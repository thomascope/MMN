Overlap_fnames = {'./VBM_bycond/Positive_Correlates_001unc.nii';
    './VBM/all.nii';
    './Dworetsky networks/HCP_Probabilistic_Network_Maps_t88_333/Language.nii'
    './Dworetsky networks/HCP_Probabilistic_Network_Maps_t88_333/Fronto-parietal.nii'
    './Dworetsky networks/HCP_Probabilistic_Network_Maps_t88_333/Cingulo-opercular.nii'
    './Dworetsky networks/HCP_Probabilistic_Network_Maps_t88_333/Salience.nii'
%     './Dworetsky networks/Probabilistic_Network_Maps_t88_333/Language.nii'
%     './Dworetsky networks/Probabilistic_Network_Maps_t88_333/Fronto-parietal.nii'
%     './Dworetsky networks/Probabilistic_Network_Maps_t88_333/Cingulo-opercular.nii'
%     './Dworetsky networks/Probabilistic_Network_Maps_t88_333/Salience.nii'
    };
spm_reslice(Overlap_fnames)
for i = 1:length(Overlap_fnames)
    [PATHSTR,NAME,EXT] = fileparts(Overlap_fnames{i});
    Overlap_fnames{i} = [PATHSTR filesep 'r' NAME EXT];
end
all_fstructs = spm_vol(Overlap_fnames);

%First calculate proportion of voxels that overlap with the probablistic maps thresholded at 5%
disp('Image 1 is VBM MMN correlates, 0.001unc, Image 2 is Language probablistic map 5% thresholded')
[totvox1,totvox2,totvox3] = propvox(all_fstructs{1}, 0.1, all_fstructs{3}, 5); %NB: DO NOT THRESHOLD AT ZERO because of interpolation in the reslicing.

disp('Image 1 is VBM MMN correlates, 0.001unc, Image 2 is Multiple Demand probablistic map 5% thresholded')
[totvox1,totvox2,totvox3] = propvox(all_fstructs{1}, 0.1, all_fstructs{4}, 5);

disp('Image 1 is VBM MMN correlates, 0.001unc, Image 2 is Cingulo-opercular probablistic map 5% thresholded')
[totvox1,totvox2,totvox3] = propvox(all_fstructs{1}, 0.1, all_fstructs{5}, 5);

disp('Image 1 is VBM MMN correlates, 0.001unc, Image 2 is Salience probablistic map 5% thresholded')
[totvox1,totvox2,totvox3] = propvox(all_fstructs{1}, 0.1, all_fstructs{6}, 5);

%Now check specificiity of this, by looking at overall atrophy map thresholded FWE 0.05, (t = 4.75)
disp('Image 1 is VBM overall, FWE 0.05, Image 2 is Language probablistic map 5% thresholded')
[totvox1,totvox2,totvox3] = propvox(all_fstructs{2}, 4.75, all_fstructs{3}, 5);

disp('Image 1 is VBM oveall, FWE 0.05, Image 2 is Multiple Demand probablistic map 5% thresholded')
[totvox1,totvox2,totvox3] = propvox(all_fstructs{2}, 4.75, all_fstructs{4}, 5);

disp('Image 1 is VBM overall, FWE 0.05, Image 2 is Cingulo-opercular probablistic map 5% thresholded')
[totvox1,totvox2,totvox3] = propvox(all_fstructs{2}, 4.75, all_fstructs{5}, 5);

disp('Image 1 is VBM oveall, FWE 0.05, Image 2 is Salience probablistic map 5% thresholded')
[totvox1,totvox2,totvox3] = propvox(all_fstructs{2}, 4.75, all_fstructs{6}, 5);

%Now do weighted multiplications to get overall association strength scores
mask = './control_majority_unsmoothed_mask_c1_thr0.05_cons0.8.img'; %5% grey matter mask

disp('Image 1 is VBM MMN correlates, 0.001unc, Image 2 is Language probablistic map')
[total1,total2,total3] = weightedvox(all_fstructs{1}, all_fstructs{3}, mask);

disp('Image 1 is VBM MMN correlates, 0.001unc, Image 2 is Multiple Demand probablistic map')
[total1,total2,total3] = weightedvox(all_fstructs{1}, all_fstructs{4}, mask);

disp('Image 1 is VBM MMN correlates, 0.001unc, Image 2 is Cingulo-opercular probablistic map')
[total1,total2,total3] = weightedvox(all_fstructs{1}, all_fstructs{5}, mask); 

disp('Image 1 is VBM MMN correlates, 0.001unc, Image 2 is Salience probablistic map')
[total1,total2,total3] = weightedvox(all_fstructs{1}, all_fstructs{6}, mask);

% Again, control analysis with VBM
disp('Image 1 is VBM overall, FWE 0.05, Image 2 is Language probablistic map')
[total1,total2,total3] = weightedvox(all_fstructs{2}, all_fstructs{3}, mask); 

disp('Image 1 is VBM overall, FWE 0.05, Image 2 is Multiple Demand probablistic map')
[total1,total2,total3] = weightedvox(all_fstructs{2}, all_fstructs{4}, mask);

disp('Image 1 is VBM overall, FWE 0.05, Image 2 is Cingulo-opercular probablistic map')
[total1,total2,total3] = weightedvox(all_fstructs{2}, all_fstructs{5}, mask);

disp('Image 1 is VBM overall, FWE 0.05, Image 2 is Salience probablistic map')
[total1,total2,total3] = weightedvox(all_fstructs{2}, all_fstructs{6}, mask);


