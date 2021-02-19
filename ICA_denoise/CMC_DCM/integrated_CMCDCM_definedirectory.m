function integrated_CMCDCM_definedirectory(megfilename,time_window,p)

% %% First load Tallie's version of SPM, with her additional scripts, if not already loaded
% old_path = path;
% cleanupObj = onCleanup(@()restore_env(old_path));
%
% spmpath = '/group/language/data/thomascope/MMN/ICA_denoise/Tallie_extDCM/spm12_latestTA/';
% thisspm = which('spm');
% if ~strcmp(thisspm(1:end-5), spmpath)
%     rmpath(genpath('/imaging/local/software/spm_cbu_svn/releases/spm12_latest/'));
%     rmpath(genpath('/group/language/data/thomascope/spm12_fil_r6906/'));
%     addpath(spmpath)
%     spm eeg
% end
%
% %addpath('/group/language/data/thomascope/MMN/ICA_denoise/Tallie_extDCM/mfiles_also_needed');
% addpath('/group/language/data/thomascope/MMN/ICA_denoise/CMC_DCM/TallieScripts/')

%%
DCM.outdir = p.CMC_DCM_outdir;


%% DCM

datafit = {};
for i = 2:length(p.conditions) %Assume STD is first
    datafit{i-1} = [p.conditions{1} '_' p.conditions{i}];
    datacons(i-1,:) = [1,i];
end

Twin = time_window;

%LF = 2; HF = 30;

%%
for dc = 1:length(datafit)
    for tw = 1:size(Twin,1)
        
        % Parameters and options used for setting up model.
        %-------------------------------------------------------
        DCMbase.options.analysis = 'ERP'; % analyze evoked responses
        DCMbase.options.model = 'CMC'; % ERP model
        DCMbase.options.spatial = 'LFP'; % spatial model OR ECD??? Manual recommends ECD
        DCMbase.options.trials  = datacons(dc,:);  % index of ERPs within ERP/ERF file
        DCMbase.options.Tdcm(1) = Twin(tw,1);      % start of peri-stimulus time to be modelled
        DCMbase.options.Tdcm(2) = Twin(tw,2);    % end of peri-stimulus time to be modelled
        DCMbase.options.Fdcm   = [4 48];
        DCMbase.options.Nmodes  = 8;      % nr of modes for data selection
        DCMbase.options.h       = 4;      % nr of DCT components
        DCMbase.options.onset   = 60;     % selection of onset (prior mean)
        DCMbase.options.D       = 1;      % downsampling
        DCMbase.options.han     = 1;      % Hanning
        DCMbase.options.DoData = 1;
        
        DCMbase.options.lock     = 0;      % (if want modulations to be in same direction for all connections)
        DCMbase.options.location = 0;      % (Not relevant; only for ECD)
        DCMbase.options.symmetry = 0;      % (Not relevant; only for ECD)
        
        DCMbase.options.nograph  = 0;
        
        % Alex custom options [added to spm_dcm_erp_data]
        %DCMbase.options.Bdcm     = [-100 0]; % baseline times
        DCMbase.options.Fltdcm   = [1 15];   % bandpass filter
        
        % location priors for dipoles
        %----------------------------------------------------------
        
        %DCMbase.Lpos = [[-42; -22; 7] [-61; -32; 8]  [-46; 20; 8] [-49; -38; 38] [46; -14; 8] [59; -25; 8] [46; 20; 8]  [57; -38; 42]];
        DCMbase.Sname = {'left AI', 'left STG', 'left IFG', 'left IPC', 'right A1',  'right STG', 'right IFG',  'right IPC'};
        Nareas  = 8;%size(DCMbase.Lpos,2);
        Ns = Nareas; % In AS script Ns = 6, which I think corresponds to number of nodes
        %----------------------------------------------------------
        % between trial effects
        %----------------------------------------------------------
        DCMbase.xU.X = [0 1]'; % [std dvt]'
        DCMbase.xU.name = {'rare'};
        
        %----------------------------------------------------------
        % specify connectivity model
        %----------------------------------------------------------
        %%
        model = [];
        
        % Bilateral inputs into A1
        m = 1;
        model(m).A{1} = zeros(8);             % forward connection
        model(m).A{2} = zeros(8);             % backward connection
        model(m).A{3} = zeros(8);             % lateral
        model(m).B{1} = zeros(8);             % Null model
        model(m).C    = [1 0 0 0 1 0 0 0]';         % Inputs into A1
        
        % add intrinsic connections
        model(m).B{1}(1,1)=1; %A1 modulation on itself
        model(m).B{1}(5,5)=1; %A1 modulation on itself
        
        % Add STG
        model(m).A{1}(2,1) = 1; % LA1 forward connection on LSTG
        model(m).A{1}(6,5) = 1; % RA1 forward connection on RSTG
        model(m).A{2}(1,2) = 1; % LA1 backward connection on LSTG
        model(m).A{2}(5,6) = 1; % RA1 backward connection on RSTG
        
        model(m).B{1}(2,1) = 1; % LA1 modulation on LSTG forward
        model(m).B{1}(6,5) = 1; % RA1 modulation on LSTG foward
        model(m).B{1}(1,2) = 1; % LA1 modulation on LSTG Backward
        model(m).B{1}(5,6) = 1; % RA1 modulation on RSTG backward
        
        % Add RIFG
        model(m).A{1}(7,6) = 1; % RSTG forward connection on RIFG
        model(m).A{2}(6,7) = 1; % RIFG backward connection on RSTG
        model(m).B{1}(7,6) = 1; %!
        model(m).B{1}(6,7) = 1; %!
        
        model(m).A{1}(1,1) = 1; % A1 self connections
        model(m).A{2}(1,1) = 1;
        model(m).A{1}(5,5) = 1;
        model(m).A{2}(5,5) = 1;
        
        % MODEL 2 lIFG and LSTG
        m=2;
        model(m) = model(1);
        model(m).A{1}(3,2) = 1; % LSTG forward connection on LIFG
        model(m).A{2}(2,3) = 1; % LIFG backward connection on LSTG
        model(m).B{1}(3,2) = 1; %!
        model(m).B{1}(2,3) = 1; %!
        
        % MODEL 3 add RP (8) on RSTG (6)
        m=3;
        model(m) = model(1);
        model(m).A{1}(8,6) = 1; % RSTG forward connection on RP
        model(m).A{2}(6,8) = 1; % RP backward connection on RSTG
        model(m).B{1}(8,6) = 1; %!
        model(m).B{1}(6,8) = 1; %!
        
        % MODEL 4 LR IFG and Ps (7 8) on STGs (3 4 resp)
        m=4;
        model(m) = model(3);
        model(m).A{1}(4,2) = 1; % LSTG forward connection on LP
        model(m).A{2}(2,4) = 1; % LP backward connection on LSTG
        model(m).B{1}(4,2) = 1; %!
        model(m).B{1}(2,4) = 1; %!
        
        model(m).A{1}(3,2) = 1; % LSTG forward connection on LIFG
        model(m).A{2}(2,3) = 1; % LIFG backward connection on LSTG
        model(m).B{1}(3,2) = 1; %!
        model(m).B{1}(2,3) = 1; %!
        
        % MODEL 5 RIFG and RP (8) on RIFG (6)
        m=5;
        model(m) = model(1);
        %     model(m).A{1}(8,7) = 1; % RIFG forward connection on RP
        %     model(m).A{2}(7,8) = 1; % RP backward connection on RIFG
        model(m).A{1}(7,8) = 1; % RIFG BACKWARD connection on RP
        model(m).A{2}(8,7) = 1; % RP FORWARD connection on RIFG
        model(m).B{1}(8,7) = 1; %!
        model(m).B{1}(7,8) = 1; %!
        
        % MODEL 6 LR IFG and Ps (7 8) on IFGs (5 6 resp)
        m=6;
        model(m) = model(5);
        %     model(m).A{1}(4,3) = 1; % LIFG forward connection on LP
        %     model(m).A{2}(3,4) = 1; % LP backward connection on LIFG
        model(m).A{1}(3,4) = 1; % LIFG BACKWARD connection on LP
        model(m).A{2}(4,3) = 1; % LP FORWARD connection on LIFG
        model(m).B{1}(4,3) = 1; %!
        model(m).B{1}(3,4) = 1; %!
        
        model(m).A{1}(3,2) = 1; % LSTG forward connection on LIFG
        model(m).A{2}(2,3) = 1; % LIFG backward connection on LSTG
        model(m).B{1}(3,2) = 1; %!
        model(m).B{1}(2,3) = 1; %!
        
        % MODEL 7 RIFG and RP (8) on RSTG (4) and RIFG (6)
        m=7;
        model(m) = model(3);
        %     model(m).A{1}(8,7) = 1; % RIFG forward connection on RP
        %     model(m).A{2}(7,8) = 1; % RP backward connection on RIFG
        model(m).A{1}(7,8) = 1; % RIFG BACKWARD connection on RP
        model(m).A{2}(8,7) = 1; % RP FORWARD connection on RIFG
        model(m).B{1}(8,7) = 1; %!
        model(m).B{1}(7,8) = 1; %!
        
        % MODEL 8 LR IFG and Ps (7 8) on STGs (3 4 resp) and IFGs (5 6 resp)
        m=8;
        model(m) = model(4); %#ok<*SAGROW>
        %     model(m).A{1}(8,7) = 1; % RIFG forward connection on RP
        %     model(m).A{2}(7,8) = 1; % RP backward connection on RIFG
        model(m).A{1}(7,8) = 1; % RIFG BACKWARD connection on RP
        model(m).A{2}(8,7) = 1; % RP FORWARD connection on RIFG
        model(m).B{1}(8,7) = 1; %!
        model(m).B{1}(7,8) = 1; %!
        %     model(m).A{1}(4,3) = 1; % LIFG forward connection on LP
        %     model(m).A{2}(3,4) = 1; % LP backward connection on LIFG
        model(m).A{1}(3,4) = 1; % LIFG BACKWARD connection on LP
        model(m).A{2}(4,3) = 1; % LP FORWARD connection on LIFG
        model(m).B{1}(4,3) = 1; %!
        model(m).B{1}(3,4) = 1; %!
        
        model(m).A{1}(3,2) = 1; % LSTG forward connection on LIFG
        model(m).A{2}(2,3) = 1; % LIFG backward connection on LSTG
        model(m).B{1}(3,2) = 1; %!
        model(m).B{1}(2,3) = 1; %!
        
        % Repeat for inputs into IFG
        m = 9;
        model(m) = model(1);
        model(m).C    = [1 0 0 0 1 0 0 1]';   % Inputs into A1 & right IFG
        
        % Repeat for inputs into IFG
        m = 10;
        model(m) = model(2);
        model(m).C    = [1 0 0 1 1 0 0 1]';   % Inputs into A1 & right IFG
        
        % Repeat for inputs into IFG
        m = 11;
        model(m) = model(3);
        model(m).C    = [1 0 0 0 1 0 0 1]';   % Inputs into A1 & right IFG
        
        % Repeat for inputs into IFG
        m = 12;
        model(m) = model(4);
        model(m).C    = [1 0 0 1 1 0 0 1]';   % Inputs into A1 & right IFG
        
        % Repeat for inputs into IFG
        m = 13;
        model(m) = model(5);
        model(m).C    = [1 0 0 0 1 0 0 1]';   % Inputs into A1 & right IFG
        
        % Repeat for inputs into IFG
        m = 14;
        model(m) = model(6);
        model(m).C    = [1 0 0 1 1 0 0 1]';   % Inputs into A1 & right IFG
        
        % Repeat for inputs into IFG
        m = 15;
        model(m) = model(7);
        model(m).C    = [1 0 0 0 1 0 0 1]';   % Inputs into A1 & right IFG
        
        % Repeat for inputs into IFG
        m = 16;
        model(m) = model(8);
        model(m).C    = [1 0 0 1 1 0 0 1]';   % Inputs into A1 & right IFG
        
        % Repeat first 16 models with lateral connections
        
        % all models have STG lateral connections
        for m = 17:32
            model(m) = model(m-16);
            model(m).A{1}(2,6) = 1;   % Lateral connections STG
            model(m).A{1}(6,2) = 1;
            model(m).A{2}(2,6) = 1;   % Lateral connections STG
            model(m).A{2}(6,2) = 1;
            model(m).B{1}(2,6) = 1;   % Lateral connections STG
            model(m).B{1}(6,2) = 1;
        end
        
        % Even models also have IPC and IFG lateral connections
        for m = 18:2:32
            model(m) = model(m);
            model(m).A{1}(3,7) = 1; % Lateral connections
            model(m).A{1}(7,3) = 1;
            model(m).A{1}(4,8) = 1;
            model(m).A{1}(8,4) = 1;
            model(m).A{2}(3,7) = 1; % Lateral connections
            model(m).A{2}(7,3) = 1;
            model(m).A{2}(4,8) = 1;
            model(m).A{2}(8,4) = 1;
            model(m).B{1}(3,7) = 1;   % Lateral connections
            model(m).B{1}(7,3) = 1;
            model(m).B{1}(4,8) = 1;
            model(m).B{1}(8,4) = 1;
        end
        
        
        LogEvd=[]; DCMname={};
        
        DCMsub = DCMbase;
        
        data_subj = sprintf([pwd filesep megfilename]);
        
        S=[]; S.D = data_subj;
        
        DCMsub.xY.Dfile = S.D;
        DCMsub.xY.modality = 'LFP';
        %DCMsub.M.dipfit.sens = 'meg';
        
        
        DCMsub.options.gety = 0;
        DCMsub.options.nograph  = 1;
        
        for n=1:numel(model)
            %n
            DCM      = DCMsub;
            DCM.name = sprintf([p.CMC_DCM_outdir 'mod_' num2str(n) '_' nout(2,@fileparts,[megfilename(1:end-4) '_dcm']) '_' datafit{dc} '.mat']);
            
            DCM.A = model(n).A;
            DCM.B = model(n).B;
            DCM.C = model(n).C;
            
            DCM = CustomPriors(DCM,Ns);
            
            DCM = spm_dcm_erp_dataTA(DCM,DCM.options.h);
            %DCM = spm_dcm_erp_dipfit(DCM, 0);
            
            if exist(DCM.name, 'file')
                disp(['file exists for model ' num2str(n) ', ' megfilename(1:end-4) ', ' datafit{dc} ', moving on.']);
            else
                disp(['Processing model ' num2str(n) ', ' megfilename(1:end-4) ', ' datafit{dc} '.']);
                
                DCM   = spm_dcm_erp_AS(DCM);
                
            end
            %LogEvd(n) = DCM.F;
            
        end % end of models
        
    end % End of time windows
        
end %End of contrasts


