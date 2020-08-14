function [Q,M] = spm_gen_QTAint33(P,X,varargin)
% Helper routine for spm_gen routines
% FORMAT [Q] = spm_gen_Q(P,X)
%
% P - parameters
% X - vector of between trial effects
% c - trial in question
%
% Q - trial or condition-specific parameters
%
% This routine computes the parameters of a DCM for a given trial, where
% trial-specific effects are deployed according to a design vector X. The
% parameterisation follows a standard naming protocol where, for example,
% X(1)*P.B{1} + X(2)*P.B{2}... adjusts P.A for all (input) effects encoded
% in P.B.
% P.BN and P.AN operate at NMDA receptors along extrinsic connections
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_gen_Q.m 6725 2016-02-19 19:14:25Z karl $


% condition or trial specific parameters
%==========================================================================
if isfield(P,'B')
    Q = rmfield(P,'B');
else
    Q = P;
end


% trial-specific effects on C (first effect only)
%--------------------------------------------------------------------------
try
    Q.C = Q.C(:,:,1) + X(1)*P.C(:,:,2);
end

% trial-specific effects on A (connections)
%--------------------------------------------------------------------------
for i = 1:length(X)
    
    % extrinsic (driving) connections
    %----------------------------------------------------------------------
    for j = 1:length(Q.A)
        
        Q.A{j} = Q.A{j} + X(i)*P.B{i};
        
        % CMM-NMDA specific modulation on extrinsic NMDA connections
        %------------------------------------------------------------------
        if isfield(P,'AN')
            Q.AN{j} = Q.AN{j} + X(i)*P.BN{i};
        end
        
    end
    
    % modulatory connections
    %----------------------------------------------------------------------
    if isfield(P,'M')
        Q.M  = Q.M + X(i)*P.N{i};
    end
    
    % intrinsic connections
    %----------------------------------------------------------------------
    % TA:
    if isfield(Q,'G')
        if any(P.B{i}(:)) && X(i)~=0
            Q.G = Q.G + X(i);
            Q.G(end-(12-1):end,1) = Q.G(end-(12-1):end,1) - X(i);
        end
    end
    
    % intrinsic connections
    %----------------------------------------------------------------------
    % TA:
    if nargin==3
        M = varargin{1};
        if isfield(M,'int')
            for j = 1:numel(M.int)
                %if isfield(M.int{j},'B')
                    M.int{j}.G = X(i)*M.int{j}.B;
                %else
                %    M.int{j}.G(:,1) = M.int{j}.G(:,1) + X(i)*P.B{i}(j,j);
                %end
            end
        end
    end
    
end
