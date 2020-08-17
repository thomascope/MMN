function [L] = spm_lx_erp_cmmNMDA_symflexTA(P,dipfit)
%
% TA: ASSUMES J IS 16 IN THE SECOND DIM
%
% observer matrix for a neural mass model: y = G*x
% FORMAT [G] = spm_lx_erp(P,dipfit)
% FORMAT [G] = spm_lx_erp(P,M)
%
% M.dipfit - spatial model specification
%
% G        - where y = L*x; G = dy/dx
% x        - state vector
%__________________________________________________________________________
%
% David O, Friston KJ (2003) A neural mass model for MEG/EEG: coupling and
% neuronal dynamics. NeuroImage 20: 1743-1755
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_lx_erp.m 5939 2014-04-06 17:13:50Z karl $

% extract dipfit from model if necessary
%--------------------------------------------------------------------------
if  isfield(dipfit,'dipfit'), dipfit = dipfit.dipfit; end
if ~isfield(dipfit,'type'),   dipfit = 'LFP';         end

% parameterised lead field times source contribution to ECD
%--------------------------------------------------------------------------
L       = spm_erp_L(P,dipfit);               % lead field per source % the gain of each source electrode for LFP
% L       = kron(P.J,L);                       % lead-field per state

l = kron(P.J(1,:),L);

% Alex code: we want sep J [contributions] for each node but symmatrically
% confined, so 3 nodes (IFG, STG & A1). The order of nodes is:
% 1 = LIFG, 2 = LSTG, 3=LA1, 4 = RIFG, 5 = RSTG, 6 = RA1

try
    o = struct();
    for k = 1:size(P.J,1)
        o.(['L' k]) = kron(P.J(k,:),L);
    end

    L       = l*0; % get sparse size
    n = (size(L,1)/2);
    for k = 1:n
        L(k,k)  = o.(['L' k])(k,k);
        L(k+n,k+n) = o.(['L' k])(k+n,k+n);

        L(k,k+6) = o.(['L' k])(k,k+6);
        L(k+n,k+n+6) = o.(['L' k])(k+n,k+n+6);
        
        L(k,k+18) = o.(['L' k])(k,k+18);
        L(k+n,k+n+18) = o.(['L' k])(k+n,k+n+18);
    end
catch err
    if isnumeric(P.J)
        L = kron(P.J,L);                       % lead-field per state
    else

        % construct lead field for each source
        %----------------------------------------------------------------------
        for i = 1:numel(P.J)
            G{i} = L(:,i)*P.J{i};
        end
        L = spm_cat(G);
    end
end

