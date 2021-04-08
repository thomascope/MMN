% PEB analysis
% 21/10/2019

clear all
clc
addpath(genpath('/imaging/local/software/spm_cbu_svn/releases/spm12_latest'));
addpath(genpath('/imaging/rowe/archive/users/hp02/pnfa_mmn/dcm/CMC_DCM_all_subjs_together_camcanHCs/customPriors/LFP/PEB/plottingTA'));

cd('/imaging/rowe/archive/users/hp02/pnfa_mmn/dcm/CMC_DCM_all_subjs_together_camcanHCs/2019_TC_LFPs_customPriors_32models/PEB')

% Load Participants
load('/imaging/rowe/archive/users/hp02/pnfa_mmn/dcm/CMC_DCM_all_subjs_together_camcanHCs/2019_TC_LFPs_customPriors_32models/Participant.mat');

% path to DCMs:
dcm_path = '/imaging/rowe/archive/users/hp02/pnfa_mmn/dcm/CMC_DCM_all_subjs_together_camcanHCs/2019_TC_LFPs_customPriors_32models/';
dvts = {'std_dvt', 'std_loc','std_int','std_dur','std_gap','std_frq'};


% Set up GCM - the fully connected DCM model - 32

% Comparisons to choose between:
comparison = {'All', 'HCvsAllPats', 'HCvsBV', 'HCvsPCA','HCvsPPA','HCvsMCI', 'ADvsFTLD','BVvsPCA','BVvsPPA','BVvsMCI', 'PCAvsPPA','PCAvsMCI','PPAvsMCI'};

%%
for c = 2:length(comparison);
    comparison{c}
    for d = 2:6
        dvts{d}
        clear X groups M
        % Groups and contrast model:
        [X, groups] = Mcontrast(comparison{c}, Participant);
        
        % Create a GCM for all groups together to begin with:
        GCM = {};
        for ss = 1:size(groups,1)
            clear DCM
            load([dcm_path, dvts{d}, '/mod32_',  groups{ss,1}, '_' dvts{d}, '.mat']);
            GCM{end+1,1} = DCM;
        end
        
        M = GCM{1}.M;
        M.X = X;
        
        
        
        [PEB, DCM2] = spm_dcm_peb(GCM,M, {'B'});
        BMA = spm_dcm_peb_bmc(PEB);
        
        %spm_dcm_peb_review(BMA,GCM)% produces matrices for connectivities of each A matrix
        
        % DCM average
        DCMav = spm_dcm_average(DCM2, 'DCMaverage');
        
        
        %     BMA.Ep = reshape(BMA.Ep,[],size(M.X,2));
        %     BMA.Pp = reshape(BMA.Pp,[],size(M.X,2));
        %
        %     n = [];
        %     for k=1:length(PEB.Pnames)
        %         str = strsplit(PEB.Pnames{k},'(');
        %         str = strsplit(str{2},')');
        %         str = strsplit(str{1},',');
        %         n(k,:) = cellfun(@(x) str2num(x),str);
        %     end
        
        %% Extract significant connections from the BMA output
        % Have to convert BMA arrays into matrix to extract
        
        
        HVAL = BMA.Pp >.95;
        clear loc
        for k = 1:length(BMA.Pnames)
            
            if strcmp('A',BMA.Pnames{k}(1))
                loc(k,:) = [str2num(BMA.Pnames{k}(3)) str2num(BMA.Pnames{k}(6)) str2num(BMA.Pnames{k}(8))];
            elseif strcmp('B',BMA.Pnames{k}(1))
                loc(k,:) = [(str2num(BMA.Pnames{k}(3))+4) str2num(BMA.Pnames{k}(6)) str2num(BMA.Pnames{k}(8))];
            elseif strcmp('C',BMA.Pnames{k}(1))
                loc(k,:) = [6 1 str2num(BMA.Pnames{k}(3))];
            
            end
            
            
        end
        hval = {zeros(8,8)};% zeros(8,8) zeros(8,8) zeros(8,8) zeros(8,8) zeros(1,8)};
        for j = 1:size(loc,1)
            hval{loc(j,1)-4}(loc(j,2),loc(j,3)) = HVAL(j+24);%+98 for ABC
        end
        
        savefolder = [dvts{d}, '/PEB_mod32_B_', comparison{c}, '.mat'];
        
        save(savefolder,'PEB','GCM','DCM2', 'DCMav','BMA', 'hval')

    end
end



function [X, groups] = Mcontrast(Comparison, P)

% M contrast matrices
% Comparing HCs to each patient group when all are included
 clear X
 ss=0;
 
switch Comparison
    
    case 'All' % Compare all possible combinations within all groups
       
        for i = 1:length(P)
            
            if ~strcmp(P{1,i}.diag, 'MCI-no')
                ss = ss+1;
               
                X(ss,1) = 1; % first column is ones
                
                groups{ss,1} = P{1,i}.name;
                % separate groups out
                
                switch P{1,i}.diag
                    case 'Control'
                        groups{ss,2} = 1;
                        
                        % Comparisons: CvsPat, CvsBV, CvsPCA,CvsPPA,CvsMCI, 
                        %ADvsFTLD,'BVvsPCA','BVvsPPA','BVvsMCI', 'PCAvsPPA','PCAvsMCI','PPAvsMCI' 
                        X(ss, 2:13) = [ 1  1  1  1  1  0  0  0  0  0  0  0];
                    
                    case 'bvFTD'
                        groups{ss,2} = 2;
                        X(ss, 2:13) = [-1 -1  0  0  0 -1  1  1  1  0  0  0];
                        
                    case 'pca'
                        groups{ss,2} = 3;
                        X(ss, 2:13) = [-1  0 -1  0  0  1 -1  0  0  1  1  0];
                        
                    case 'nfvppa'
                        groups{ss,2} = 4;
                        X(ss, 2:13) = [-1  0  0 -1  0 -1  0 -1  0 -1  0  1];
                        
                    case 'MCI-yes'
                        groups{ss,2} = 5;
                        X(ss, 2:13) = [-1  0  0  0 -1  1  0  0 -1  0 -1 -1];
                end
                
                
            end
        end
    case 'HCvsAllPats'

        for i = 1:length(P)
            if ~strcmp(P{1,i}.diag, 'MCI-no')
                ss = ss+1;
                X(ss,1) = 1; % first column is ones
                
                groups{ss,1} = P{1,i}.name;
                % separate groups out
                
                switch P{1,i}.diag
                    case 'Control'
                        groups{ss,2} = 1;
                        
                        % Comparisons: CvsPat 
                        X(ss, 2) = [ 1];
                    
                    otherwise
                        groups{ss,2} = 2;
                        X(ss, 2) = [-1];
                end
            end
        end
        
    case 'HCvsBV'

        for i = 1:length(P)
            switch P{1,i}.diag
                case 'Control'
                    ss = ss+1;
                    X(ss,1) = 1; % first column is ones
                    
                    % Comparisons: CvsBV
                    X(ss, 2) = [ 1];
                    
                    groups{ss,1} = P{1,i}.name;
                    groups{ss,2} = 1;

                case 'bvFTD'
                    ss = ss+1;
                    X(ss,1) = 1; % first column is ones
                    
                    % Comparisons: CvsBV
                    X(ss, 2) = [-1];
                    
                    groups{ss,1} = P{1,i}.name;
                    groups{ss,2} = 2;
            end
        end
        
        
    case 'HCvsPCA'
         for i = 1:length(P)
            switch P{1,i}.diag
                case 'Control'
                    ss = ss+1;
                    X(ss,1) = 1; % first column is ones
                    
                    % Comparisons: CvsPCA
                    X(ss, 2) = [ 1];
                    
                    groups{ss,1} = P{1,i}.name;
                    groups{ss,2} = 1;

                case 'pca'
                    ss = ss+1;
                    X(ss,1) = 1; % first column is ones
                    
                    % Comparisons: CvsPCA
                    X(ss, 2) = [-1];
                    
                    groups{ss,1} = P{1,i}.name;
                    groups{ss,2} = 2;
            end
         end
        
    case 'HCvsPPA'
        for i = 1:length(P)
            switch P{1,i}.diag
                case 'Control'
                    ss = ss+1;
                    X(ss,1) = 1; % first column is ones
                    
                    % Comparisons: CvsnfvPPa
                    X(ss, 2) = [ 1];
                    
                    groups{ss,1} = P{1,i}.name;
                    groups{ss,2} = 1;

                case 'nfvppa'
                    ss = ss+1;
                    X(ss,1) = 1; % first column is ones
                    
                    % Comparisons: CvsnfvPPA
                    X(ss, 2) = [-1];
                    
                    groups{ss,1} = P{1,i}.name;
                    groups{ss,2} = 2;
            end
        end
        
    case 'HCvsMCI'
        for i = 1:length(P)
            switch P{1,i}.diag
                case 'Control'
                    ss = ss+1;
                    X(ss,1) = 1; % first column is ones
                    
                    % Comparisons: CvsPCA
                    X(ss, 2) = [ 1];
                    
                    groups{ss,1} = P{1,i}.name;
                    groups{ss,2} = 1;

                case 'MCI-yes'
                    ss = ss+1;
                    X(ss,1) = 1; % first column is ones
                    
                    % Comparisons: CvsPCA
                    X(ss, 2) = [-1];
                    
                    groups{ss,1} = P{1,i}.name;
                    groups{ss,2} = 2;
            end
        end
    case 'ADvsFTLD' % AD pathologies vs FTLD pathologies
        for i = 1:length(P)
            switch P{1,i}.diag
                case 'bvFTD'
                    ss = ss+1;
                    X(ss,1) = 1; % first column is ones
                    
                    % Comparisons: AD vs FTLD
                    X(ss, 2) = [ -1];
                    
                    groups{ss,1} = P{1,i}.name;
                    groups{ss,2} = 2;
                    
                case 'pca'
                    ss = ss+1;
                    X(ss,1) = 1; % first column is ones
                    
                    % Comparisons:AD vs FTLD
                    X(ss, 2) = [1];
                    
                    groups{ss,1} = P{1,i}.name;
                    groups{ss,2} = 1;
                case 'nfvppa'
                    ss = ss+1;
                    X(ss,1) = 1; % first column is ones
                    
                    % Comparisons:AD vs FTLD
                    X(ss, 2) = [-1];
                    
                    groups{ss,1} = P{1,i}.name;
                    groups{ss,2} = 2;
                    
                case 'MCI-yes'
                    ss = ss+1;
                    X(ss,1) = 1; % first column is ones
                    
                    % Comparisons: AD vs FTLD
                    X(ss, 2) = [1];
                    
                    groups{ss,1} = P{1,i}.name;
                    groups{ss,2} = 1;
            end
            
        end
        
    case 'BVvsPCA'
        for i = 1:length(P)
            switch P{1,i}.diag
                case 'bvFTD'
                    ss = ss+1;
                    X(ss,1) = 1; % first column is ones
                    
                    % Comparisons: BVvsPCA
                    X(ss, 2) = [ 1];
                    
                    groups{ss,1} = P{1,i}.name;
                    groups{ss,2} = 1;

                case 'pca'
                    ss = ss+1;
                    X(ss,1) = 1; % first column is ones
                    
                    % Comparisons: BVvsPCA
                    X(ss, 2) = [-1];
                    
                    groups{ss,1} = P{1,i}.name;
                    groups{ss,2} = 2;
            end
        end
    case 'BVvsPPA'
        for i = 1:length(P)
            switch P{1,i}.diag
                case 'bvFTD'
                    ss = ss+1;
                    X(ss,1) = 1; % first column is ones
                    
                    % Comparisons: BVvsPPA
                    X(ss, 2) = [ 1];
                    
                    groups{ss,1} = P{1,i}.name;
                    groups{ss,2} = 1;

                case 'nfvppa'
                    ss = ss+1;
                    X(ss,1) = 1; % first column is ones
                    
                    % Comparisons: BVvsPPA
                    X(ss, 2) = [-1];
                    
                    groups{ss,1} = P{1,i}.name;
                    groups{ss,2} = 2;
            end
        end
    case 'BVvsMCI'
        for i = 1:length(P)
            switch P{1,i}.diag
                case 'bvFTD'
                    ss = ss+1;
                    X(ss,1) = 1; % first column is ones
                    
                    % Comparisons: BVvsMCI
                    X(ss, 2) = [ 1];
                    
                    groups{ss,1} = P{1,i}.name;
                    groups{ss,2} = 1;

                case 'MCI-yes'
                    ss = ss+1;
                    X(ss,1) = 1; % first column is ones
                    
                    % Comparisons: BVvsMCI
                    X(ss, 2) = [-1];
                    
                    groups{ss,1} = P{1,i}.name;
                    groups{ss,2} = 2;
            end
        end
    case 'PCAvsPPA'
        for i = 1:length(P)
            switch P{1,i}.diag
                case 'pca'
                    ss = ss+1;
                    X(ss,1) = 1; % first column is ones
                    
                    % Comparisons: PCAvsPPA
                    X(ss, 2) = [ 1];
                    
                    groups{ss,1} = P{1,i}.name;
                    groups{ss,2} = 1;

                case 'nfvppa'
                    ss = ss+1;
                    X(ss,1) = 1; % first column is ones
                    
                    % Comparisons: PCAvsPPA
                    X(ss, 2) = [-1];
                    
                    groups{ss,1} = P{1,i}.name;
                    groups{ss,2} = 2;
            end
        end
    case 'PCAvsMCI'
         for i = 1:length(P)
            switch P{1,i}.diag
                case 'pca'
                    ss = ss+1;
                    X(ss,1) = 1; % first column is ones
                    
                    % Comparisons: PCAvsMCI
                    X(ss, 2) = [ 1];
                    
                    groups{ss,1} = P{1,i}.name;
                    groups{ss,2} = 1;

                case 'MCI-yes'
                    ss = ss+1;
                    X(ss,1) = 1; % first column is ones
                    
                    % Comparisons: PCAvsMCI
                    X(ss, 2) = [-1];
                    
                    groups{ss,1} = P{1,i}.name;
                    groups{ss,2} = 2;
            end
        end
    case 'PPAvsMCI'
         for i = 1:length(P)
            switch P{1,i}.diag
                case 'nfvppa'
                    ss = ss+1;
                    X(ss,1) = 1; % first column is ones
                    
                    % Comparisons: PPAvsMCI
                    X(ss, 2) = [ 1];
                    
                    groups{ss,1} = P{1,i}.name;
                    groups{ss,2} = 1;

                case 'MCI-yes'
                    ss = ss+1;
                    X(ss,1) = 1; % first column is ones
                    
                    % Comparisons: PPAvsMCI
                    X(ss, 2) = [-1];
                    
                    groups{ss,1} = P{1,i}.name;
                    groups{ss,2} = 2;
            end
         end
    otherwise
        disp('Comparision does not exist');
end
P = P;