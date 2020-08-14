function D = PP_AllFilter(D)


D = spm_eeg_load(D);

S   = [];
S.D = D;
S.type  = 'butterworth';
S.order = 5;
S.band  = 'high';
S.dir = 'twopass';
S.freq  = 1;
S.prefix= '';

D = spm_eeg_filter(S);

S   = [];
S.D = D;
S.type  = 'butterworth';
S.order = 5;
S.band  = 'low';
S.dir = 'twopass';
S.freq  = 180;
S.prefix= 'F';

D = spm_eeg_filter(S);


end

