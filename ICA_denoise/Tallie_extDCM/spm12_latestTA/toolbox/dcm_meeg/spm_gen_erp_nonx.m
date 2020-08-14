function [y,pst] = spm_gen_erp(P,M,U)
% Generates a prediction of trial-specific source activity
% FORMAT [y,pst] = spm_gen_erp(P,M,U)
%
% P - parameters
% M - neural-mass model structure
% U - trial-effects
%   U.X  - between-trial effects (encodes the number of trials)
%   U.dt - time bins for within-trial effects
%
% y   - {[ns,nx];...} - predictions for nx states {trials}
%                     - for ns samples
% pst - peristimulus time (seconds)
%
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_gen_erp.m 6427 2015-05-05 15:42:35Z karl $

% default inputs - one trial (no between-trial effects)
%--------------------------------------------------------------------------
if nargin < 3, U.X = sparse(1,0); end

% check input u = f(t,P,M) and switch off full delay operator
%--------------------------------------------------------------------------
try, M.fu; catch, M.fu  = 'spm_erp_u'; end
try, M.ns; catch, M.ns  = 128;         end
try, M.N;  catch, M.N   = 0;           end
try, U.dt; catch, U.dt  = 0.004;       end

% peristimulus time
%--------------------------------------------------------------------------
if nargout > 1
    pst = (1:M.ns)*U.dt - M.ons/1000;
end

% within-trial (exogenous) inputs TA MOVED INSIDE C LOOP BELOW
%==========================================================================
% if ~isfield(U,'u')
%     
%     % peri-stimulus time inputs
%     %----------------------------------------------------------------------
%     if ~iscell(M.fu)
%         U.u = feval(M.fu,(1:M.ns)*U.dt,P,M);
%     else U.u = feval(M.fu{:});
%     end
%     
% end
% 
% if isfield(M,'u')
%     
%     % remove M.u to preclude endogenous input
%     %----------------------------------------------------------------------
%     M = rmfield(M,'u');
%     
% end

% between-trial (experimental) inputs
%==========================================================================
if isfield(U,'X')
    X = U.X;
else
    X = sparse(1,0);
end

if ~size(X,1)
    X = sparse(1,0);
end


% cycle over trials
%==========================================================================
y      = cell(size(X,1),1);
for  c = 1:size(X,1)
    
    % within-trial (exogenous) inputs
    %==========================================================================
    %if ~isfield(U,'u')
        % peri-stimulus time inputs
        if ~iscell(M.fu)
            U.u = feval(M.fu,(1:M.ns)*U.dt,P,M);
        else
            if nargin(M.fu{1})==1
                U.u = feval(M.fu{:},c);
            else U.u = feval(M.fu{1},c,(1:M.ns)*U.dt,P,M);
            end
        end

    %end
    if isfield(M,'u')
        % remove M.u to preclude endogenous input
        M = rmfield(M,'u');
    end
    
    % condition-specific parameters
    %----------------------------------------------------------------------
    Q   = feval(M.genQ,P,X(c,:));
    
    % solve for steady-state - for each condition
    %----------------------------------------------------------------------
    M.x  = spm_dcm_neural_x_ignored(Q,M);
    
    % integrate DCM - for this condition
    %----------------------------------------------------------------------
    y{c} = spm_int_L(Q,M,U);
    
end

