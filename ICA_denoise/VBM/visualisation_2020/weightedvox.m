function [total1,total2,total3,num] = weightedvox(im1, im2, mask)
% calculate the normalised ratio of activation intensities within a mask,
% normalised to the activation extents
% Written by TEC Apr 2021

temp1 = spm_read_vols(spm_vol(im1));
temp2 = spm_read_vols(spm_vol(im2));

total1 = nansum(nansum(nansum(temp1)));
total2 = nansum(nansum(nansum(temp2)));

spm_imcalc([im1,im2],'temp.nii','i1.*i2');
temp3 = spm_read_vols(spm_vol('temp.nii'));

total3 = nansum(nansum(nansum(temp3)));

mask_vol = spm_read_vols(spm_vol(mask));
mask_size = sum(sum(sum(mask_vol)));

normalised_ratio = (total3*mask_size)/(total1*total2);

disp(['Normalised ratio of ' num2str(normalised_ratio) ' multiplying image 1 by image 2'])
end
