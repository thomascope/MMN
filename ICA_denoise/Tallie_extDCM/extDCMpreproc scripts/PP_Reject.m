function D = PP_Reject(D)

%addpath(genpath('/home/as08/old_spm12/'));

% tmp = load(D);
D = spm_eeg_load(D);
S   = []; 
S.D = D.fullfile;

eog_thr = [1000];%[100e-6];
MEGMAGthr = 3000;
MEGPLANthr = 100;
EEGthr = 500; % the norm for this might be 150 - so this is a very high threshold, i.e. they have more freedom of movement
% fac = 10;
% i = find(arrayfun(@(x) strcmp('MEGMAG',x.type),tmp.D.channels));
% d = squeeze(tmp.D.data(i,:,:));
% MEGMAGthr = mean(d(:)) + (fac*std(d(:)));
% i = find(arrayfun(@(x) strcmp('MEGPLANAR',x.type),tmp.D.channels));
% d = squeeze(tmp.D.data(i,:,:));
% MEGPLANthr = mean(d(:)) + (fac*std(d(:)));
% i = find(arrayfun(@(x) strcmp('EEG',x.type),tmp.D.channels));
% d = squeeze(tmp.D.data(i,:,:));
% EEGthr = mean(d(:)) + (fac*std(d(:)));
% clear tmp d


S.methods(1).fun = 'flat';
S.methods(1).channels = 'MEG';
S.methods(1).settings.threshold = 0;
S.methods(1).settings.seqlength = 4;

S.methods(end+1).fun = 'flat';
S.methods(end).channels = 'EEG';
S.methods(end).settings.threshold = 0;
S.methods(end).settings.seqlength = 4;

%  S.methods(end+1).fun = 'peak2peak';
S.methods(end+1).fun = 'threshchan';
S.methods(end).channels = 'EOG';
S.methods(end).settings.threshold = eog_thr;

S.methods(end+1).fun = 'peak2peak';
S.methods(end).channels = 'MEGMAG';
S.methods(end).settings.threshold = MEGMAGthr;%20e-2;%20e-12

S.methods(end+1).fun = 'peak2peak';
S.methods(end).channels = 'MEGPLANAR';
S.methods(end).settings.threshold = MEGPLANthr;%200e-1;%200e-12

S.methods(end+1).fun = 'peak2peak';
S.methods(end).channels = 'EEG';
S.methods(end).settings.threshold = EEGthr;%100e-2;%100e-6


D = spm_eeg_artefact(S);

try nbadchan = length(D.badchannels);end
try nrejects = sum(D.reject);        end

try fprintf('setting %d bad channels\n',nbadchan);end
try fprintf('rejecting %d bad trials\n',nrejects);end



end