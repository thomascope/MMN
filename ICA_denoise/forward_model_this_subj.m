function forward_model_this_subj(megpath, mripath, inv_meth, time_wind_path, windows)
%forward models a cellstring of megpaths against a cellstring of mripaths
addpath(genpath('/imaging/local/software/mne'));
addpath(genpath('/imaging/hp02/mmn_08/analysis_spm/new_spm_functions'));

%% Define loops

% % the inversion methods to be performed
inver_meth_path = { ''};%, 'MEGgrad_MSP/', 'MEGgrad_sLoreta/', 'MEGgrad_Beamf/' };%, 'MEGgrad_Beam/'};
inv_trls        = {'STD', 'DVT'};

%% Define steps to be done
forward_modeling = 1;
segment_only = 0;    % Flag whether only MRI segmentation (up to, not including, coregistration) shall be done (1: only segment; 0: do it all)
redoreg_flag = 0;    % Flag whether to redo MRI and registration model (eg if want to change fids)
redoformod_flag = 0; % Flag whether to redo forward modelling (eg if only one for_typ above)
display_flag = 0;    % display results on the fly
write3Dflag = 1;     % write SPM volumes of results
DoInv_contrast = 1;

%% Define subjects' preprocessed files

nr_subs = length(megpath);
fprintf(1, 'Going to process %d data sets\n', nr_subs);


%% Define Analysis Parameters/Processing options/ Inversion options that don't change:
% Names of conditions after averaging
% inv_trls = {'SPEC_ALL', 'FREE_ALL', 'NULL', 'ACTION'};    % Labels of conditions corresponding to trigger codes    % Labels of conditions corresponding to trigger codes

freq_start = 0.1;
freq_end = 100;

% Forward model options, for_typ{EEG/MEG}
for_typ{1} = 'Single Shell'; % MEG (mags+grads)
for_typ{2} = 'EEG BEM'; % EEG

% Which sensor types to use
inv_mods = {'MEGMAG';'MEGPLANAR'};   % Not EEG - sunbalanced acquisition by group

% Cortical surface option
mesh_size = 2;   % [1-3] for coarse-normal-fine

% Whether to use headshape for coregistration (in datareg)
use_headshape = 1;



%% Loop around each inversion method and each time window:

for inv_cnt = 1:length(inv_meth)
    
    % Define Analysis Parameters
    % Which inversion method to use
    inv_typ = inv_meth{inv_cnt};   % Minimum Norm Least-Squares %inv_typ = 'GS';   % (MSP)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Source localisation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    clear D
    if forward_modeling == 1
        for ss = 1:nr_subs
            
            % Make sure in right starting place:
            [sourceloc_path, megfname, ~] = fileparts(megpath{ss});
            cd(sourceloc_path)
            
            D = spm_eeg_load(megpath{ss});
            
            % Initialise... (if want to be safe!)
            val = inv_cnt;
            D.inv{val} = {struct('mesh', [], 'gainmat', [])};

            
            %% MRI processing and coregistration
            if redoreg_flag == 1 || ~isfield(D.inv{1},'mesh') || ~isfield(D.inv{1},'datareg')
                
                D.val = val;
                D.inv{val}.date    = strvcat(date,datestr(now,15));
                D.inv{val}.comment = {sprintf('%s_%s',inv_typ,char(inv_mods)')};    % remember inversion type and sensor configuration
                D.inv{val}.mesh    = [];
                D.inv{val}.datareg = [];
                D.inv{val}.forward = [];
                
                % Locations/names of MRIs
                assert(exist(mripath,'file')>0,['The MRI does not exist at ' mripath ])
                D.inv{val}.mesh.sMRI = mripath;
                
                %% Normalise sMRI (if not done already), and create inverse-normalised surfaces
                D.inv{val}.mesh = spm_eeg_inv_mesh(D.inv{val}.mesh.sMRI, mesh_size);
                if display_flag, spm_eeg_inv_checkmeshes(D); end
                D.save;
                
                % If coregistration not yet done...
                if segment_only,
                    continue;   % Stop processing for this subject here, continue with next one
                end;
                
                % If coregistration done manually and saved...
                newmrifid           = [];
                newmrifid.fid.pnt   = D.inv{val}.mesh.fid.fid.pnt(1:3,:);  % These are the fids after "Save"...
                newmrifid.fid.label = {'Nasion';'LPA';'RPA'};
                newmrifid.pnt       = D.inv{val}.mesh.fid.pnt;        % Scalp mesh points from MRI above
                
                meegfid = D.fiducials;
                %% Remove nose points (assume those for which y>0 and z<0)
                meegfid.pnt(find(meegfid.pnt(:,2)>0 & meegfid.pnt(:,3)<0),:) = [];
                
                fprintf(1, 'Coregistering\n');
                %D = spm_eeg_inv_datareg_noui_260416_pnfa(D, val, meegfid, newmrifid,use_headshape);
                D = spm_eeg_inv_datareg_ui(D, val, meegfid, newmrifid,use_headshape);
                fprintf(1, 'Done Coregistering\n');
                
                for ind = 1:length(D.inv{val}.datareg)
                    fprintf(1, '%d ', ind);
                    d = D.inv{val}.datareg(ind).fid_eeg.fid.pnt - D.inv{val}.datareg(ind).fid_mri.fid.pnt;
                    %err(ss,ind,val) = mean(sqrt(sum(d.^2,2)));
                end
                fprintf(1, '\n');
                
                redoreg_flag = 0;
                
            else
                D.inv{val}.mesh    = D.inv{1}.mesh;
                D.inv{val}.datareg = D.inv{1}.datareg;
            end % if redoreg_flag...
            
            %% Computing forward model/leadfield
            if redoformod_flag == 1 || ~isfield(D.inv{val},'forward') || ~isfield(D.inv{val}, 'gainmat')
                
                fprintf(1, 'Creating forward model\n');
                %% Create forward model (BEM) (could conditionalise this bit on modality inverted...)
                
                D.inv{val}.forward = struct([]);
                if length(D.inv{val}.datareg) == 1
                    clear for_typ
                    for_typ{1} = 'Single Shell'; % MEG (mags+grads)
                    
                elseif length(D.inv{val}.datareg) ==2
                    clear for_typ
                    for_typ{1} = 'EEG BEM'; % EEG
                    for_typ{2} = 'Single Shell'; % MEG (mags+grads)
                end
                
                for ind = 1:length(D.inv{val}.datareg)
                    D.inv{val}.forward(ind).voltype = for_typ{ind};
                end
                
                fprintf(1, 'Computing leadfield\n');
                D = spm_eeg_inv_forward(D,val);
                
                %             if display_flag
                %                 for ind = 1%:length(D.inv{val}.datareg) % !!!!!!!!!!just MEG, not EEG
                %                     spm_eeg_inv_checkforward_hp_130614(D, val, ind); %pause(3);
                %                 end
                %             end
                current_formod = val;
                redoformod_flag = 0;
                D.save;
            else
                D.inv{val}.forward = D.inv{current_formod}.forward;
                D.inv{val}.gainmat = D.inv{current_formod}.gainmat;
                D.save;
            end
        end
    end
end

for inv_cnt = 1:length(inv_meth)
    
    % Define Analysis Parameters
    % Which inversion method to use
    inv_typ = inv_meth{inv_cnt};   % Minimum Norm Least-Squares %inv_typ = 'GS';   % (MSP)
    val = inv_cnt;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Source localisation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    clear D
    if forward_modeling == 1
        for ss = 1:nr_subs
            %% Invert & Contrast
            if DoInv_contrast == 1
                
                for wind_cnt = 1:length(time_wind_path)
                    S={};
                    % Make a copy of the forward model file and rename for each
                    % time window:
                    [filepath,name,ext] = fileparts(megpath{ss});
                    if inv_cnt == 1
                        S.D = megpath{ss};
                    else
                        S.D = [filepath '/' 's_' time_wind_path{wind_cnt} '_' inv_meth{inv_cnt-1} '_' name];
                    end
                    S.outfile = [filepath '/' 's_' time_wind_path{wind_cnt} '_' inv_meth{inv_cnt} '_' name];
                    D = spm_eeg_copy(S);
                    
                    
                    D.inv{val}.inverse = [];   % Clear to be safe!
                    D.inv{val}.inverse.trials = inv_trls;
                    D.inv{val}.inverse.type   = inv_typ;
                    
                    D.inv{val}.inverse.woi    = [windows{wind_cnt}(1) windows{wind_cnt}(2)];
                    %D.inv{val}.inverse.lpf    = freq_start;
                    %D.inv{val}.inverse.hpf    = freq_end;
                    
                    D.inv{val}.inverse.modality = inv_mods;
                    
%                     try
                        D = spm_eeg_invert(D,val);
%                     catch
%                         D.inv{val}.inverse.modality = {'MEGPLANAR'} % Sometimes fails if there are insufficient spatial modes in the MEGMAGs after maxfilter
%                         D = spm_eeg_invert(D,val)
%                     end
                    
                    
                    D.inv{val}.contrast.woi  = [windows{wind_cnt}(1) windows{wind_cnt}(2)];
                    %D.inv{val}.contrast.fboi = [freq_start freq_end];
                    D.inv{val}.contrast.type = 'evoked';
                    
                    D = spm_eeg_inv_results(D);
                    
                    
                    %% Create images for stats for different smoothing:
                    
                    
                    % Write result to SPM volumes
                    if write3Dflag
                        D.inv{val}.contrast.smoothing = 8;%sm*4;
                        
                        D = spm_eeg_inv_Mesh2Voxels(D);
                        %SourceImgs{val}{ss} = strvcat(D.inv{val}.contrast.fname);
                    end
                    
                    
                    D.save;
                end
            end
        end
    end
end