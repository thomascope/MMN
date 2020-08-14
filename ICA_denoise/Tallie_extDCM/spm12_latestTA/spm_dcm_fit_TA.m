function [P]   = spm_dcm_fit_nonx(P)
% Bayesian inversion of DCMs using Variational Laplace
% FORMAT [DCM] = spm_dcm_fit(P)
%
% P    - {N x M} DCM structure array (or filenames) from N subjects
%
% DCM  - Inverted (1st level) DCM structures with posterior densities
%__________________________________________________________________________
%
% This routine is just a wrapper that calls the appropriate dcm inversion
% routine for a set a pre-specifed DCMs.
%
% If called with a cell array, each column is assumed to contain 1st level
% DCMs inverted under the same model. Each row contains a different data
% set (or subject).
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_dcm_fit.m 6716 2016-02-08 18:21:37Z peter $


% get filenames and set up
%--------------------------------------------------------------------------
if ~nargin
    [P, sts] = spm_select([1 Inf],'^DCM.*\.mat$','Select DCM*.mat files');
    if ~sts, return; end
end
if ischar(P),   P = cellstr(P);  end
if isstruct(P), P = {P};         end

% Number of subjects (data) and models (of those data)
%--------------------------------------------------------------------------
[Ns,Nm] = size(P);

% Find model class and modality
%==========================================================================
try, load(P{1}); catch, DCM = P{1}; end

model = spm_dcm_identify(DCM);

if isempty(model)    
    warning('unknown inversion scheme');
    return   
end

% get data structure for each subject (column)
%------------------------------------------------------------------
for i = 1:Ns
    for j = 2:Nm
        switch model
            
            case{'DEM'}
                P{i, j}.xY = P{i, 1}.Y;
            otherwise
                P{i, j}.xY = P{i, 1}.xY;
        end
    end
end

%matlabpool open 7

% loop over subjects (columns)
%--------------------------------------------------------------------------
pp = cbupool(length(P)); % 49); % 
pp.ResourceTemplate = '-l nodes=^N^,mem=16GB,walltime=96:00:00';
parpool(pp,length(P)); %  49); % length(P)); %

parfor i = 1:numel(P)
    
    % loop over models (rows)
    %----------------------------------------------------------------------
    
    % Get model specification
    %------------------------------------------------------------------
    try, P{i} = load(P{i}); catch, P{i} = P{i}; end
    
    
    % invert and save
    %==================================================================
    switch model
        
                % conventional neural-mass and mean-field models
                %----------------------------------------------------------
            case{'ERP'}
                P{i} = spm_nlsi_NplusTA(P{i});
                
                % cross-spectral density model (complex)
                %----------------------------------------------------------
            case{'CSD'}
                P{i} = spm_dcm_csd_as(P{i});
                
                % phase coupling
                %----------------------------------------------------------
            case{'PHA'}
                P{i} = spm_dcm_phase(P{i});
                
                % generic nonlinear system identification
                %----------------------------------------------------------
            case{'NLSI'}
                [Ep,Cp,Eh,F] = spm_nlsi_GN(P{i}.M,P{i}.xU,P{i}.xY);
                P{i}.Ep       = Ep;
                P{i}.Eh       = Eh;
                P{i}.Cp       = Cp;
                P{i}.F        = F;
                
        otherwise
            try
                P{i} = feval(model, P{i});
            catch
                error('unknown DCM');
            end
    end
    
    % place inverted model in output array
    %------------------------------------------------------------------
    P{i} = P{i};
end


%matlabpool close