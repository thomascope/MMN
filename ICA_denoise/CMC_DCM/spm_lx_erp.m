function [L] = spm_lx_erp(P,dipfit)


%%HP HAS MODIFIED THIS ONE - GO TO SPM_LX_ERP_ORIG FOR ORIGINAL



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
L       = spm_erp_L(P,dipfit);               % lead field per source

% L       = kron(P.J,L);                       % lead-field per state

l = kron(P.J(1,:),L);

% Alex code: we want sep J [contributions] for each node but symmatrically
% confined, so 3 nodes (IFG, STG & A1). The order of nodes is:
% 1 = LA1, 2 = RA1, 3=LSTG, 4 = RSTG, 5 = LIFG, 6 = RIFG, 7 = LIPC, 8 =
% RIPC

L1 = kron(P.J(1,:),L); % A1's
L2 = kron(P.J(2,:),L); % STG's
L3 = kron(P.J(3,:),L); % IFG's
L4 = kron(P.J(4,:),L); % IPC's

L       = l*0; % get sparse size
L(1,1)  = L1(1,1);
L(2,2)  = L2(2,2);
L(3,3)  = L3(3,3);
L(4,4)  = L4(4,4);
L(5,5)  = L1(5,5);
L(6,6)  = L2(6,6);
L(7,7)  = L3(7,7);
L(8,8)  = L4(8,8);
L(1,17) = L1(1,17);
L(2,18) = L2(2,18);
L(3,19) = L3(3,19);
L(4,20) = L1(4,20);
L(5,21) = L2(5,21);
L(6,22) = L3(6,22);
L(7,23) = L2(7,23);
L(8,24) = L3(8,24);

L(1,49) = L1(1,49);
L(2,50) = L2(2,50);
L(3,51) = L3(3,51);
L(4,52) = L1(4,52);
L(5,53) = L2(5,53);
L(6,54) = L3(6,54);
L(7,53) = L2(7,53);
L(8,54) = L3(8,54);