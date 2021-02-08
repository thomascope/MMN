function eximgs = bet_TEC(imgs, fthresh)
% wrapper for Steve Smith's Brain Extraction Tool
% FORMAT eximgs = bet(imgs, fthresh)
%
% imgs   - char or vol struct list of images to run BET on
% fthresh - fractional threshold to use (default 0.5)
% eximgs - vol struct array of extracted images
%
% Matthew Brett 28/10/00
% 
% Updated by Thomas Cope for new cluster 02082021

  
if nargin < 1
  imgs = spm_get(Inf, 'img', 'Images to skull strip');
end
if ischar(imgs)
  imgs = spm_vol(imgs);
end
if nargin < 2
  fthresh = spm_input('Fractional threshold', 1, 'r', 0.5, 1, [0 1]);
end



prefix = 'ss_';
eximgs = [];
for i=1:length(imgs)
  img = imgs(i);
  [pn fn e] = fileparts(img.fname);
  if isempty(pn)
      pn = '.';
  end
  unix(sprintf('FSLOUTPUTTYPE=NIFTI; bet %s%c%s %s%c%s%s -f %0.2f ', ...
	       pn,filesep,fn,pn,filesep, ...
	       prefix,fn, fthresh));
  outfil = fullfile(pn, [prefix fn e]);
  if ~exist(outfil)
          warning(['No output file created for ' img.fname]);
  else
    spm_get_space(outfil, img.mat);
    eximgs = [eximgs; spm_vol(outfil)];
  end
end






