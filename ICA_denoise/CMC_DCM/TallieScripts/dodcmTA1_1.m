function dodcmTA1(fn,KK)


suffix = ['_Mod_' num2str(KK)]; % extra addition to outfile name if desired
cond = {'Dev'}; % really don't get this. afterall, I'm looking at two conditions, but this is what DCM produces left to it's own auto-generative devises..
model_number = KK; % when you beuild he model list with MODELSTA (below), this allows you to select the one you want, rather than having to rebuild every time you want a change

% [A,C,L,B] = MODELSTA;
% A = A{model_number};
% A{3} = L{model_number};
% B = B{model_number};
% C = C{model_number};
% clear L
[M,C,L,B] = MODELSTA;
A{1} = M{model_number};
A{1} = A{1}.*tril(ones(size(A{1},1)));
A{2} = M{model_number};
A{2} = A{2}.*tril(ones(size(A{1},1)));%-eye(size(A{1},1));
A{3} = L{model_number};
B = B(model_number);
C = C{model_number};


%% setup:

%clear DCM

options.trials      = [1 2 7];                                              % condition list indices to run
options.analysis    = 'ERP';                                                % analysis method - ERP will average trials by default!
options.model       = 'CMC';                                                % chosen neuronal model
options.spatial     = 'ECD';                                                % chosen spatial model
options.Tdcm        = [-100 300];                                           % time of interest (ms from stim trigger - asumed to be at 0ms)
options.Fdcm        = [4 40];                                               % frequencies to bother with..
options.Rft         = 5;                                                    % ???
options.onset       = 64;                                                   % prior assumed onset of activity - how much wiggle room there is, I don't know, I think it's in the docs
options.dur         = 16;                                                   % prior assumed onset of activity - how much wiggle room there is, I don't know, I think it's in the docs
options.Nmodes      = 8;                                                    % I think this is spectral modes .. er, see documentation
options.h           = 1;                                                    % style of detrending, I think this is linear (see documentation)
options.han         = 1;                                                    % if you want a hanning over the whole trial window
options.D           = 1;                                                    % ???
options.lock        = 0;                                                    % ???
options.multiC      = 0;                                                    % ???
options.location    = 0;                                                    % ???
options.symmetry    = 0;                                                    % ???
options.DoData      = 1;                                                    % custom addition - alex

xU.X                = [0 1 0; 0 0 1; 1 0 0];                                % a design matrix for trial comparisons
xU.name             = {'rp1' 'rp6' 'dev'};                                  % names for the rows of the design matrix

%fn = [nout(1,@fileparts,F{k}) filesep 'ICMmaiefff' nout(2,@fileparts,F{k}) '.mat']; %'G:\Roving_MMN_TGB\C1\TAtest\LFP6_2.mat'; % alter as appropriate
d = load(fn);

fi = find(arrayfun(@(x) strcmp('MEGMAG',x.type),d.D.channels));
lpos = d.D.sensors.meg.chanpos;
lpos = lpos(fi,:);
nam = arrayfun(@(x) x.label,d.D.channels,'Uni',0);
nam = nam(fi);

Fs = d.D.Fsample;
t = [1/Fs:1/Fs:1/Fs*d.D.Nsamples]+d.D.timeOnset;
sten = [findnearest(options.Tdcm(1),t*1000) findnearest(options.Tdcm(end),t*1000)];

xY.Dfile            = fn;
xY.modality         = 'MEG';                                                % assuming fn contains reconstructed sources of your ROIs
xY.name             = nam; %arrayfun(@(x) x.label,d.D.channels,'Uni',0);          % names of your ROIs (e.g. leftSTG, rightA1, etc)
xY.Ic               = size(lpos,1); %1:size(d.D.data,1);
xY.Time             = t;
xY.dt               = 1/Fs;
xY.coor2D           = [d.D.channels.X_plot2D ; d.D.channels.Y_plot2D];
xY.pst              = t(sten(1):sten(end));
xY.It               = sten(1):sten(end);
xY.nt               = length(find(strcmp(cond{1},arrayfun(@(x) x.label,d.D.trials,'Uni',0))));

DCM.options         = options;
DCM.xU              = xU;
DCM.xY              = xY;
DCM.Sname           = xY.name';
DCM.Lpos            = lpos'; %double.empty(3,0);                                    % I believe these are only needed if you are doing ECD and haven't inputted LFPs
DCM.A               = A;                                                    % extrinsic connectivity matrix (sources x sources)
DCM.B               = B;                                                    % things you allow to be modulated differently between the chosen conditions (sources x sources)
DCM.C               = C;                                                    % external input matrix (sources x 1)
DCM.name            = [fn(1:end-4) 'dcm' suffix '.mat'];                    % output file name


%% AS stuff method 1:
% 
% SYM_TF = 0;
% DCM.M.U            = sparse(diag(ones(size(D.data,2),1)));  ... ignore [modes]
% DCM.options.DoData = 1;             ... leave on [custom]
% DCM.options.Bdcm   = [-200 0];      ... baseline times [new!]
% DCM.CUSTOM.nograph = 0;
% DCM.symmetricTF    = SYM_TF;
% 
% 
%% AS stuff method 2:

Ns              = size(DCM.A{1},1);
DCM.CUSTOM.f    = 'spm_fx_cmcTA'; %KRISH13';          ... generative model  (.m)
DCM.CUSTOM.pE.G = zeros(Ns,13);                 ... [local coupling]
DCM.CUSTOM.pC.G = zeros(Ns,13);                 ... variance [off]
Self = find(diag(speye(Ns).*( DCM.A{1}+DCM.A{2} ))); % SP gain only if in model
DCM.CUSTOM.pC.G(Self,[1 4 7 10]) = 1/8;
DCM.CUSTOM.pE.T = zeros(Ns,4);                  ... population time const
DCM.CUSTOM.pC.T = zeros(Ns,4)+1/8;              ... variances
DCM.CUSTOM.gE.J = sparse([1 3 7],1,[.2 .8 .2],8,1)'; ... contributing states
DCM.CUSTOM.gC.J = sparse([1 3 7],1,1/8       ,8,1)'; ... variances
%DCM.symmetricTF = 0; % custom addition - TA mod of alex


%% run:

spm_dcm_erp(DCM);


end

