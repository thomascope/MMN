%% DCM script for MEG MMN 

clear all
clc

cwd = '/imaging/mlr/users/tc02/pnfa_mmn/dcm/CMC_DCM_all_subjs_together_camcanHCs/2019_TC_LFPs_customPriors_32models/';
cd(cwd)
%%

addpath(genpath('/imaging/local/software/mne'));
addpath(genpath('/imaging/rowe/archive/users/hp02/spm12b'));
%addpath(genpath('/imaging/local/software/spm_cbu_svn/releases/spm12_fil_r6906'));

% Adding AS08 scripts for specialised adaptations for CMC modelling
% addpath(genpath('/imaging/local/software/spm_cbu_svn/releases/spm12_latest'));
% addpath('/imaging/as08/Roving/NEW_DCM_ERP_SPM12/');
% addpath('/imaging/as08/Roving/NEW_DCM_ERP_SPM12/cons1/');
addpath('/imaging/rowe/archive/users/hp02/pnfa_mmn/dcm/CMC_DCM_all_subjs_together_camcanHCs/TallieScripts');

%% Define your subject data and MRIs
prefix1 = 'fmraedfff';
path_file = '/imaging/rowe/archive/users/hp02/pnfa_mmn/dcm/CMC_DCM_all_subjs_together_camcanHCs/2019_TC_LFPs_customPriors_32models';
%group_folders = {'ftd', 'pca', 'vespa_patients', 'vespa_hcs'};

addpath(path_file)
%forward_subject_names;
load('/imaging/mlr/users/tc02/Holly_MMN/ICA_denoise_longwindow/LFPs/Participant.mat');
%cd(cwd)
% Note: you may have to copy your MRI first (e.g. from /imaging/local/structurals/cbu)

nr_subs = length(Participant);
fprintf(1, 'Going to process %d data sets\n', nr_subs);


%% Copy source loc files

%% Matlab pool stuff
% Start MATLAB pool
openparallel = 0;


 if openparallel && numel(gcp('nocreate')) == 0
     %if matlabpool('size')==0
         
         cbupool(80,'--mem-per-cpu=8G --time=167:00:00')

     %end
 end


%% DCM

datafit = {'std_dvt', 'std_loc','std_int','std_dur','std_gap','std_frq'};
datacons = [1,2; 1,3; 1,4; 1,5; 1,6; 1,7];

Twin = [0 300];


LF = 2; HF = 30;

%%
for dc = 5:length(datafit)
    for tw = 1:size(Twin,1)
        
        % Parameters and options used for setting up model.
        %-------------------------------------------------------
        DCMbase.options.analysis = 'ERP'; % analyze evoked responses
        DCMbase.options.model = 'CMC'; % ERP model
        DCMbase.options.spatial = 'LFP'; % spatial model OR ECD??? Manual recommends ECD
        DCMbase.options.trials  = [1, dc+1];  % index of ERPs within ERP/ERF file
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
        DCMbase.options.Bdcm     = [-100 0]; % baseline times
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
        
        %%
        %mods_to_run = 1:8;
        parfor ss = 1:nr_subs
            datafit = {'std_dvt', 'std_loc','std_int','std_dur','std_gap','std_frq'};
            
            ss
            DCMsub = DCMbase;
            
            data_subj = sprintf('/imaging/mlr/users/tc02/Holly_MMN/ICA_denoise_longwindow/LFPs/%s/8LFP_s_-100_500_LOR_fmbraedfffM%s.mat', Participant{ss}.diag,  Participant{ss}.name);
            
            
            S=[]; S.D = data_subj;
            
            DCMsub.xY.Dfile = S.D;
            DCMsub.xY.modality = 'LFP';
            %DCMsub.M.dipfit.sens = 'meg';
            
            
            DCMsub.options.gety = 0;
            DCMsub.options.nograph  = 1;
            
            for n=1:numel(model)
                %n
                DCM      = DCMsub;
                DCM.name = sprintf('mod%d_%s_%s',n,Participant{ss}.name,datafit{dc});
                
                DCM.A = model(n).A;
                DCM.B = model(n).B;
                DCM.C = model(n).C;
                
                DCM = CustomPriors(DCM,Ns);
                
                DCM = spm_dcm_erp_dataTA(DCM,DCM.options.h);
                %DCM = spm_dcm_erp_dipfit(DCM, 0);
                
                if exist([cwd '/' DCM.name '.mat'], 'file')
                    disp('file exists');
                else
                    n
                    DCM   = spm_dcm_erp_AS(DCM);
                    
                end
                %LogEvd(n) = DCM.F;
                
            end % end of models
            
            
            % put here because otherwise the parfor gets unhappy
            for n = 1:numel(model)
                DCM      = DCMsub;
                DCM.name = sprintf('mod%d_%s_%s',n,Participant{ss}.name,datafit{dc});
                DCMname{ss,m} = DCM.name;
                
            end
            
        end % End of subject
        
        
    end %end of the time windows
end



