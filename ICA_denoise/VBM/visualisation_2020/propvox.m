function [totvox1,totvox2,totvox3] = propvox(im1, thresh1, im2, thresh2)
% calculate proportion of voxels that overlap with given thresholds
% Written by TEC Apr 2021

temp1 = spm_read_vols(spm_vol(im1));
temp2 = spm_read_vols(spm_vol(im2));

overthresh1 = temp1>thresh1;
overthresh2 = temp2>thresh2;

totvox1 = sum(sum(sum(overthresh1)));
totvox2 = sum(sum(sum(overthresh2)));

spm_imcalc([im1,im2],'temp.nii',['(i1>' num2str(thresh1) ').*(i2>' num2str(thresh2) ')']);
temp3 = spm_read_vols(spm_vol('temp.nii'));

overthresh3 = temp3>0;
totvox3 = sum(sum(sum(overthresh3)));

disp([num2str(totvox1) ' voxels above threshold in image 1, ' num2str(totvox2) ' voxels above threshold in image 2'])
disp([num2str(totvox3) ' voxels in thresholded overlap map'])
disp([num2str(100*totvox3/totvox1) ' percent of voxels from image 1 in overlap'])
disp([num2str(100*totvox3/totvox2) ' percent of voxels from image 2 in overlap'])
end
