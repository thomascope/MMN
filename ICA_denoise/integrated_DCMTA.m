function DCM = integrated_DCMTA(megfilename,time_window,condition,flipdipoles)
% An original script by Tallie Adams to apply an extended neuronal
% model DCM to mismatch negativity data, subtly modified by Thomas Cope to be
% more generalisable in application, but with the mechanics unchanged. It
% still explicitly specifies the connectivity model.
% Condition can be specified either with a number or with a code
% Additional modification by TEC- introduction of the flipdipoles argument, which
% specifies a time period for every source that should be positive, for
% example to specify that all M100s in A1, STG and IPC should be positive
% and all second deflections in IFG positive (M100 not always reliably seen in IFG)
% for my data I would specify a vector [30 80; 30 80; 90 165; 30 80]

%% First load Tallie's version of SPM, with her additional scripts, if not already loaded
old_path = path;
cleanupObj = onCleanup(@()restore_env(old_path));

spmpath = '/group/language/data/thomascope/MMN/ICA_denoise/Tallie_extDCM/spm12_latestTA/';
thisspm = which('spm');
if ~strcmp(thisspm(1:end-5), spmpath)
    rmpath(genpath('/imaging/local/software/spm_cbu_svn/releases/spm12_latest/'));
    rmpath(genpath('/group/language/data/thomascope/spm12_fil_r6906/'));
    addpath(spmpath)
    spm eeg
end

addpath('/group/language/data/thomascope/MMN/ICA_denoise/Tallie_extDCM/mfiles_also_needed');

%% CHOOSE NEURONAL MODEL:
%%

disp(' ')
disp(' -- loading your specified neuronal model')

    function [f,J,Q] = fx(x,u,P,M)
        if nargin~=0
            ns   = size(M.x,1);                      % number of sources
            np   = size(M.x,2);                      % number of populations per source
            nk   = size(M.x,3);                      % number of states per population
            x    = reshape(x,ns,np,nk);              % hidden states
            
            % extrinsic connection strengths
            % exponential transform to ensure positivity constraints
            aA{1}  = exp(P.A{1});                      % forward
            aA{2}  = exp(P.A{2});                      % backward
            aN{1} = exp(P.AN{1});                     % forward
            aN{2} = exp(P.AN{2});                     % backward
            Cc     = exp(P.C);                         % subcortical
            
            % detect and reduce the strength of reciprocal (lateral) connections
            for ii = 1:length(aA)
                Ll    = (aA{ii} > exp(-8)) & (aA{ii}' > exp(-8));
                aA{ii} = aA{ii}./(1 + 8*Ll);
            end
            
            % intrinsic connection strengths
            % condition specific effects: Inhibition of SPGP's % TA: this is condition specific because P originates from Q, from spm_gen_Q, which accounts for the modulation matrix B.
            G = full(P.H);
            if any(P.G(:))
                for kn = 1:size(G,3)
                    ind = (np*kn)-(np-1):(np*kn);
                    for k = 1:size(G,4)
                        pg = reshape(P.G(ind,k),1,size(G,1)); % 1 % size(G,3)
                        for ki = 1:size(G,1)
                            G(ki,ki,kn,k) = squeeze(G(ki,ki,kn,k)) + pg(ki); % 1:size(G,3), % P.G
                        end
                    end
                    ind1 = size(P.G,1)-(4*kn)+1 : size(P.G,1)-(4*kn)+2;
                    ind2 = size(P.G,1)-(4*kn)+3 : size(P.G,1)-(4*kn)+4;
                    % GABAA:
                    G([1 2],3,kn,3) = G([1 2],3,kn,3) + P.G(ind1,3); % super
                    G([4 6],5,kn,3) = G([4 6],5,kn,3) + P.G(ind2,3); % deep
                    % GABAB:
                    G([1 2],3,kn,4) = G([1 2],3,kn,4) + P.G(ind1,4); % super
                    G([4 6],5,kn,4) = G([4 6],5,kn,4) + P.G(ind2,4); % deep
                    % NMDA:
                    G([2 3],2,kn,2) = G([2 3],2,kn,2) + P.G(ind1,1); % super-super
                    G([4 6],2,kn,2) = G([4 6],2,kn,2) + P.G(ind2,1); % super-deep
                    G([1 2],6,kn,2) = G([1 2],6,kn,2) + P.G(ind1,2); % deep-super
                    G([4 5],6,kn,2) = G([4 5],6,kn,2) + P.G(ind2,2); % deep-deep
                end
            end
            G = exp(G);
        end
        % connectivity switches
        % extrinsic connections (F B) - from superficial and deep pyramidal cells
        SA   = [2   0   2
            0   1   0
            0   2   0
            0   0   0
            0   0   0
            0   0   0]/8;
        
        % extrinsic NMDA-mediated connections (F B) - from superficial and deep pyramidal cells
        SNMDA =[2   0   2
            0   1   0
            0   0   0
            0   0   0
            0   0   0
            0   0   0]/8;
        
        %         SA = zeros(6,3);
        %         SNMDA = zeros(6,3);
        
        % intrinsic connections (np x np)
        GE   = [ 1    0    0    0    0    0
            2    1    0    0    0    0 % TA made ss to sp 8 instead of 4
            1    2    0    0    0    0
            1    1    0    1    0    0
            0    0    0    2    0    2
            0    1    0    0    0    2];
        
        GEn  = [ 1    0    0    0    0    0
            2    1    0    0    0    0 % TA made ss to sp 8 instead of 4
            0    0    0    0    0    0
            1    1    0    1    0    0
            0    0    0    0    0    0
            0    1    0    0    0    2];
        
        GI   = [ 8    0    2    0    0    0
            0   16    2    0    0    0
            0    0    8    0    0    0 % 32
            0    0    0   64    2    0
            0    0    0    0   16    0
            0    0    0    0    2   64]; % TA changed from 128 to 64
        
        GIb  = [ 8    0    2    0    0    0
            0   16    2    0    0    0
            0    0    8    0    0    0 % 32
            0    0    0   64    2    0
            0    0    0    0   16    0
            0    0    0    0    2   64]/2; % TA changed from 128 to 64
        
        GIm  = [ 0    0    0    0    0    0
            0    0    0    0    0    0
            0    0    0    0    0    0 % 32
            0    0    0    0    0    0
            0    0    0    0    0    0
            0    0    0    0    0    4]; % TA changed from 128 to 64
        
        GIh  = [ 0    0    0    0    0    0
            0    0    0    0    0    0
            0    0    0    0    0    0 % 32
            0    0    0    0    0    0
            0    0    0    0    0    0
            0    0    0    0    0    4]; % TA changed from 128 to 64
        if nargin==0
            f = struct;
            f.SA = SA;
            f.SNMDA = SNMDA;
            f.GE = GE;
            f.GEn = GEn;
            f.GI = GI;
            f.GIb = GIb;
            f.GIm = GIm;
            f.GIh = GIh;
        end
        if nargin~=0
            % rate constants (ns x np) (excitatory 4ms, inhibitory 16ms)
            KE    = exp(-P.T(:,1))*1000/4;                       % excitatory rate constants (AMPA)
            KNMDA = exp(-(P.T(:,2) ))*1000/100;                  % excitatory rate constants (NMDA)
            KI    = exp(-P.T(:,3))*1000/16;                      % inhibitory rate constants (GABAa)
            KIb   = exp(-P.T(:,4))*1000/200;                     % inhibitory rate constants (GABAb)
            KM    = ones(ns,1)*1000/160;
            KH    = ones(ns,1)*1000/100;
            
            % Voltages
            VL   = -70;                               % reversal  potential leak (K)
            VE   =  60;                               % reversal  potential excite (Na)
            VI   = -90;                               % reversal  potential inhib (Cl)
            VIb  = -100;                               % reversal  potential inhib (Cl)
            VR   = -40; % [-55 -55 -50 -55 -50 -55];                               % threshold potential
            VN   =  10;                               % reversal Ca(NMDA)
            VM   = -70;
            VH   = -30;
            
            CV   = exp(P.CV).*[128 128 50 400 50 200]/1000;  % membrane capacitance % 128 128 256 32 % [200 150 50 400 50 200]
            GL   = 1;                                 % leak conductance
            
            % mean-field effects:
            % neural-mass approximation to covariance of states: trial specific
            %                    ss   sp   ii   dp   di   tp
            Vx   = exp(P.S)*32.*[ 1    1    1    1    1    1]; % variance of the presynaptic pop
            
            % mean population firing and afferent extrinsic input
            h       = 1-spm_Ncdf_jdw(x(:,:,1),-100,300);
            m       = spm_Ncdf_jdw(x(:,:,1),VR,Vx);     % mean firing rate  % TA: (nxp) cumulative distribution of normal distributions (only V of x are used, leaving (nxp)
            %m(:,6) = m(:,6)*10;
            a(:,1)  = aA{1}*m(:,2);                      % forward afference  AMPA % TA: (nx2) 2 = pops sp & dp, because II is intrinsic only, and SS is input to the local area as well)
            a(:,2)  = aA{2}*m(:,4);                      % backward afference AMPA % TA: (nx2) 2 = pops sp & dp, because II is intrinsic only, and SS is input to the local area as well)
            a(:,3)  = (aA{1}*m(:,6)) + (aA{2}*m(:,6)/10);
            an(:,1) = aN{1}*m(:,2);                     % forward afference  NMDA % TA: (nx2) 2 = pops sp & dp, because II is intrinsic only, and SS is input to the local area as well)
            an(:,2) = aN{2}*m(:,4);                     % backward afference NMDA % TA: (nx2) 2 = pops sp & dp, because II is intrinsic only, and SS is input to the local area as well)
            an(:,3) = (aN{1}*m(:,6)) + (aN{2}*m(:,6)/10);
            
            BE     = exp(P.E)*0.8; % Averge background activity and exogenous input
            
            % input
            if isfield(M,'u')
                U = u(:);
            else U = Cc*u(:);
            end
            
            f     = x;
            for ii = 1:ns % nodes
                % intrinsic coupling % TA: made up of intrinsics x their weights [or existences] x firing rates. (These are all used to make the three conductance states)
                if ii<=(ns/2)
                    E      = (G(:,:,ii,1)  .*GE )*m(ii,:)';
                    ENMDA  = (G(:,:,ii,2)  .*GEn)*m(ii,:)';
                    I      = (G(:,:,ii,3)  .*GI )*m(ii,:)';
                    Ib     = (G(:,:,ii,4)  .*GIb)*m(ii,:)';
                    Im     =              GIm *m(ii,:)';
                    Ih     =              GIh *h(ii,:)';
                else
                    E      = (G(:,:,ii-(ns/2),1)  .*GE )*m(ii,:)';
                    ENMDA  = (G(:,:,ii-(ns/2),2)  .*GEn)*m(ii,:)';
                    I      = (G(:,:,ii-(ns/2),3)  .*GI )*m(ii,:)';
                    Ib     = (G(:,:,ii-(ns/2),4)  .*GIb)*m(ii,:)';
                    Im     =              GIm *m(ii,:)';
                    Ih     =              GIh *h(ii,:)';
                end
                
                % extrinsic coupling (excitatory only) and background activity TA: made up of the intrinsics (above) + weigthed background activity + weighted, firing-rate-determined excitatory connections, masked by A matrix
                E     = (E +  BE + SA*a(ii,:)')*2;
                ENMDA = (ENMDA + BE  + SNMDA*an(ii,:)')*2;
                
                % and exogenous input(U) TA: made up of the intrinsics (abovex2) + extrinsics (above) + input
                try E(1) = E(1)  + U(ii); catch  E(1) = E(1)  + U; end
                try ENMDA(1) = ENMDA(1) + U(ii); catch ENMDA(1) = ENMDA(1) + U; end
                
                % Voltage
                f(ii,:,1) =    (GL*(VL - x(ii,:,1))+...
                    x(ii,:,2).*(VE - x(ii,:,1))+...
                    x(ii,:,3).*(VI - x(ii,:,1))+...
                    x(ii,:,5).*(VIb - x(ii,:,1))+...
                    x(ii,:,6).*(VM - x(ii,:,1))+...
                    x(ii,:,7).*(VH - x(ii,:,1))+...
                    x(ii,:,4).*(VN - x(ii,:,1)).*mg_switch(x(ii,:,1)))./CV;
                
                % Conductance
                f(ii,:,2) = (E' - x(ii,:,2)).*KE(ii,:);
                f(ii,:,3) = (I' - x(ii,:,3)).*KI(ii,:);
                f(ii,:,5) = (Ib' - x(ii,:,5)).*KIb(ii,:);
                f(ii,:,6) = (Im' - x(ii,:,6)).*KM(ii,:);
                f(ii,:,7) = (Ih' - x(ii,:,7)).*KH(ii,:);
                f(ii,:,4) = (ENMDA' - x(ii,:,4))*KNMDA(ii,:) ;
            end
            f = spm_vec(f);
            
            if nargout < 2, return, end
            
            J = spm_cat(spm_diff(M.f,x,u,P,M,1)); % self-referential spm_fx_cmm_NMDA
            
            if nargout < 3, return, end
            
            % Delays
            Dd  = [2 16 80];
            d  = -Dd.*exp(P.D)/1000;
            Sp = kron(ones(nk,nk),kron( eye(np,np),eye(ns,ns)));  % states: same pop.
            Ss = kron(ones(nk,nk),kron(ones(np,np),eye(ns,ns)));  % states: same source
            
            yt = zeros(ns,np,nk);
            yt(:,6,:) = 1;
            yt = ~~reshape(yt,[],1);
            yesthal = repmat(yt,1,length(yt));
            yesthal(:,yt) = 1;
            
            Dp = ~Ss;                            % states: different sources
            Ds = ~Sp & Ss;                       % states: same source different pop.
            Ds(yesthal) = 0;
            Dtp = ~Sp & Ss & yesthal;
            Dd  = d(2)*Dp + d(1)*Ds + d(3)*Dtp;
            Q  = spm_inv(speye(length(J)) - Dd.*J);
        end
    end


%% SELECT INPUT:
%%
disp(' -- chosen input is a simple gaussian')
    function uu = fu(c,varargin)
        t = varargin{1}{1};
        PP = varargin{1}{2};
        Mm = varargin{1}{3};
        try
            if length(Mm.dur) ~= length(Mm.ons)
                Mm.dur = Mm.dur(1) + Mm.ons - Mm.ons;
            end
        catch err
            Mm.dur = 32 + Mm.ons - Mm.ons;
        end
        % check sustained input (0,1)
        try
            if length(Mm.sus) ~= length(Mm.ons)
                Mm.sus = Mm.sus(1) + Mm.ons - Mm.ons;
            end
        catch er
            Mm.sus = 0 + Mm.ons - Mm.ons;
        end
        % stimulus � Gaussian (subcortical) impulse
        nu    = length(Mm.ons);
        uu     = sparse(length(t),nu);
        t     = t*1000;
        for ii = 1:nu
            % Gaussian bump function
            delay  = Mm.ons(ii) + 128*PP.R(ii,1);
            scale1  = Mm.dur(ii) * exp(PP.R(ii,2));
            U      = exp(-(t - delay).^2/(2*scale1^2));
            % sustained inputs
            try prop = Mm.sus(ii)*exp(P.R(ii,3));
            catch err
                prop = Mm.sus(ii);
            end
            U      = prop*cumsum(U)/sum(U) + U*(1 - prop);
            uu(:,ii) = 32*U;
        end
    end


%% SELECT OPTIONS AND PRIORS (yours to change as you will):
%%

disp(' -- building DCM structure to your specification')

D         = spm_eeg_load(megfilename);
these_data = D(:,:,:); % Don't mess around with the original data
if exist('flipdipoles','var') %Flip the dipoles if this argument is specified
    assert(size(flipdipoles,1)==length(D.chanlabels)/2,'The number of rows in the time window vector must be half the number of sources')
    
    for i = 1:size(flipdipoles,1)
        if abs(max(these_data(i,D.time>=flipdipoles(i,1)/1000&D.time<=flipdipoles(i,2)/1000,1))) < abs(min(these_data(i,D.time>=flipdipoles(i,1)/1000&D.time<=flipdipoles(i,2)/1000,1)))
            these_data(i,:,:) = -these_data(i,:,:);
        end
        if abs(max(these_data(size(flipdipoles,1)+i,D.time>=flipdipoles(i,1)/1000&D.time<=flipdipoles(i,2)/1000,1))) < abs(min(these_data(size(flipdipoles,1)+i,D.time>=flipdipoles(i,1)/1000&D.time<=flipdipoles(i,2)/1000,1)))
            these_data(size(flipdipoles,1)+i,:,:) = -these_data(size(flipdipoles,1)+i,:,:);
        end
    end
end

plot_flipped_data = 0; %For checking flipping
if plot_flipped_data
    figure
    for i = 1:8
        subplot(4,4,8+i)
        plot(D.time,these_data(i,:,1),'k')
        hold on
        plot(D.time,these_data(i,:,2),'r')
        title([D.chanlabels{i} ' flipped'])
    end
    D         = spm_eeg_load(megfilename);
    for i = 1:8
        subplot(4,4,i)
        plot(D.time,D(i,:,1),'k')
        hold on
        plot(D.time,D(i,:,2),'r')
        title([D.chanlabels{i} ' original'])
        if i == 4 
            legend({'STD','DVT'})
        end
    end
end

ti        = time(D,[],'ms');
cond      = D.condlist;
mode      = 'LFP';



if isstr(condition)
    DM_names                = condition;
    DM_trials               = find(strcmp(cond,condition));
else
    DM_trials               = condition;
    DM_names                = cond(DM_trials);
end
if ~isscalar(DM_trials)
    error('A single condition needs to be specified')
end

% Need to specify which model. On the basis of Holly's first pass classic
% DCM, a fully connected model is specified. Specifically:
% The order of the nodes is as follows: left A1, left STG, left IFG, left IPC, right A1, right STG, right IFG, right IPC
% 1) The model is symmetrical
% 2) A1 projects forwards to STG, and also has intrinsic connectivity and external input
% 3) STG projects forwards to IFG and IPC, backwards to A1, and cross-hemispherically to itself
% 4) IFG projects forwards to IPC, backwards to STG, cross-hemispherically to itself, and has external input
% 5) IPC projects backwards to IFG and STG, and cross-hemispherically to itself
% However, additional intrinsic connectivity is specified because of the
% extDCM (essentially eye(num_nodes))

% forward
A{1} = [     0     0     0     0     0     0     0     0
    1     0     0     0     0     1     0     0
    0     1     0     0     0     0     1     0
    0     1     1     0     0     0     0     1
    0     0     0     0     0     0     0     0
    0     1     0     0     1     0     0     0
    0     0     1     0     0     1     0     0
    0     0     0     1     0     1     1     0];


% backward
A{2} = [     0     1     0     0     0     0     0     0
    0     0     1     1     0     1     0     0
    0     0     0     1     0     0     1     0
    0     0     0     0     0     0     0     1
    0     0     0     0     0     1     0     0
    0     1     0     0     0     0     1     1
    0     0     1     0     0     0     0     1
    0     0     0     1     0     0     0     0];

% modulatory - INGORE
A{3} = [     0     0     0     0     0     0     0     0
    0     0     0     0     0     0     0     0
    0     0     0     0     0     0     0     0
    0     0     0     0     0     0     0     0
    0     0     0     0     0     0     0     0
    0     0     0     0     0     0     0     0
    0     0     0     0     0     0     0     0
    0     0     0     0     0     0     0     0];

% across trial &/ intrinsic
B{1} = [     1     0     0     0     0     0     0     0
    0     1     0     0     0     0     0     0
    0     0     1     0     0     0     0     0
    0     0     0     1     0     0     0     0
    0     0     0     0     1     0     0     0
    0     0     0     0     0     1     0     0
    0     0     0     0     0     0     1     0
    0     0     0     0     0     0     0     1];

C = [1 0 1 0 1 0 1 0]'; % input

DCM.A = A;
DCM.B = B;
DCM.C = C;

A{1}  = DCM.A{1} | DCM.A{3};
A{2}  = DCM.A{2} | DCM.A{3};
for i = 1:2
    DCM.pE.A{i} = A{i}*32 - 32;
    DCM.pC.A{i} = A{i}/8;
end
for i = 1:2
    DCM.pE.AN{i} = A{i}*32 - 32;
    DCM.pC.AN{i} = A{i}/8;
end
for i = 1:length(DCM.B)
    B{i} = ~~DCM.B{i};
    DCM.pE.B{i} = B{i} - B{i};
    DCM.pC.B{i} = B{i}/8;
end
for i = 1:length(DCM.B)
    B{i} = ~~DCM.B{i};
    DCM.pE.BN{i} = B{i} - B{i};
    DCM.pC.BN{i} = B{i}/8;
end
C         = ~~DCM.C;
DCM.pE.C  = C*32 - 32;
DCM.pC.C  = zeros(size(C)); % C/8;

DM_design               = [1];
suffix                  = ''; % extra addition to outfile name if desired

%time_window                     = [0 400];

DCM.options.trials      = DM_trials;
DCM.options.analysis    = 'ERP';
DCM.options.model       = 'CMM_NMDAint3ext2';
DCM.options.spatial     = 'LFP';
DCM.options.Tdcm        = time_window;                                         % time window if interest (specified above)
DCM.options.Fdcm        = [0 48];                                      % frequency band of interest
DCM.options.Rft         = 5;
DCM.options.onset       = 60;                                          % gaussian input onset time (ms)
DCM.options.dur         = 8;                                           % standard deviation of the gaussian kernel
DCM.options.Nmodes      = 8;                                           % number of modes in the system (default 8)
DCM.options.h           = 1;
DCM.options.han         = 1;                                           % if you have a hanning window
DCM.options.D           = 1;
DCM.options.lock        = 0;
DCM.options.multiC      = 0;                                           % used if your input matrix has more than one column
DCM.options.location    = 0;
DCM.options.symmetry    = 0;
DCM.options.Nmax        = 64;                                          % Max number of iterations - can be ANY number you want, if you're lucky it'll converge after 30.
disp(['       time window [' num2str(time_window(1)) ' ' num2str(time_window(2)) ']'])
disp(['       running a max of ' num2str(DCM.options.Nmax) ' iterations'])

Nns                     = 7;                                           % neuronal states (V, gE, gI, gN, ...)
Np                      = 6;                                           % populations
Nn                      = size(A{1},1);                                % nodes

DCM.M.hE                = 6;       %Changed from 6 because of two conditions where NaNs were creeping into the state space matrix
DCM.M.hC                = 1/128;
DCM.M.h                 = [];
DCM.M.Xp                = [];

DCM.M.x                 = ones(Nn,Np,Nns)/8;                           % starting states
DCM.M.x(:,:,1)          = -50; % -65; % -60; % -40; %
disp(['       starting potential of ' num2str(DCM.M.x(1,1,1))])


DCM.M.f                 = @(x,u,P,M)fx(x,u,P,M); % 'spm_fx_cmm_NMDAint3ext22';
DCM.M.G                 = 'spm_lx_erp_cmmNMDA_symflexTA'; %'spm_lx_erp_cmmNMDATA'
DCM.M.genQ              = 'spm_gen_QTAint33';
DCM.M.IS                = 'spm_gen_erp_nonx'; % spm_int_L
DCM.M.FS                = 'spm_fy_erp';
DCM.M.fu                = {@(c,varargin)fu(c,varargin)}; % {@(c)fu(c)}; %

% g:
DCM.gE.J                = sparse(1,Np*Nns);
DCM.gE.J(1:Np)          = [.2 .8 0 .2 0 .2]; %[.2 0 0 0 .8 0 0 0 0 0 0 0 .2 0 0 0];
DCM.gE.J                = repmat(DCM.gE.J,[Nn/2 1]);
DCM.gC.J                = sparse(1,Np*Nns);
DCM.gC.J(1:Np)          = [1/8 1/8 1/16 1/8 1/16 1/8]; % ones(1,Np)/8;
DCM.gC.J                = repmat(DCM.gC.J,[Nn/2 1]);

DCM.gE.Lpos             = sparse(3,0);
DCM.gC.Lpos             = sparse(3,0);

DCM.gE.L                = ones(1,Nn);
DCM.gC.L                = ones(1,Nn)/45;

% p:
DCM.pE.H                = zeros(Np,Np,Nn/2,4); % ,Nn

DCM.pC.H                = zeros(Np); % eye(Np)/8;
DCM.pC.H                = repmat(DCM.pC.H,1,1,Nn/2);
DCM.pC.H                = repmat(DCM.pC.H,1,1,1,4);

DCM.pC.H([2 3 4],1,:,1) = 1/8;
DCM.pC.H([2 4],1,:,2)   = 1/8;
DCM.pC.H([3 4 6],2,:,1) = 1/8;
DCM.pC.H([4 6],2,:,2)   = 1/8;
DCM.pC.H(5,[4 6],:,1)   = 1/8;

DCM.pC.H([1 2],3,:,[3 4]) = 1/8;
DCM.pC.H([4 6],5,:,[3 4]) = 1/8;
for kkk = 1:6, DCM.pC.H(kkk,kkk,:,[3 4]) = 1/8; end

DCM.pE.T                = zeros(Nn,4);
DCM.pC.T                = ones(Nn,4)/16; % [ones(Nn,2)/16 ones(Nn,2)/8];

DCM.pE.D                = [0 0 0];
DCM.pC.D                = [1 1 1]/45;

DCM.pE.R                = zeros(size(C,2),2);
DCM.pC.R                = ones(size(C,2),1)*[1/30 1/30];

DCM.pE.S                = 0;
DCM.pC.S                = 1/45;

DCM.pE.G                = zeros((Np*(Nn/2))+(4*(Nn/2)),4);
DCM.pC.G                = zeros((Np*(Nn/2))+(4*(Nn/2)),4); % ones((Np*(Nn/2))+(4*(Nn/2)),4)/8;
DCM.pC.G(1:Np*(Nn/2),1) = 0;

DCM.pE.CV               = zeros(1,Np);
DCM.pC.CV               = ones(1,Np)/30;

DCM.pE.E                = 0;
DCM.pC.E                = 1/45;

DCM.pE.Lpos             = sparse(3,0);
DCM.pC.Lpos             = sparse(3,0);


%% ALLOW CALCULATIONS OF LEFT OVER BITS (uses the above specification to fill in other necessary inormation):
%%

disp(' -- calculating other bits of DCM structure')
% stuff decided on above &/ calculated from the above before starting
% proper; i.e. stuff you can't really change:

Ic                      = D.indchantype(mode,'GOOD');
Nc                      = length(Ic);                               % number of channels
T1                      = DCM.options.Tdcm(1);
T2                      = DCM.options.Tdcm(2);
[~, T1]                 = min(abs(ti - T1));
[~, T2]                 = min(abs(ti - T2));
It                      = (T1:DCM.options.D:T2)';
Ns                      = length(It);                               % number of samples
Nr                      = size(C,1);                                % number of sources
Nu                      = size(C,2);                                % number of exogenous inputs
Nm                      = max(DCM.options.Nmodes,Nr);               % check the number of modes is greater or equal to the number of sources

DCM.M.U                 = sparse(diag(ones(Ns,1)));

if DCM.options.multiC
    DCM.pE.C(:,:,2) = DCM.pE.C(:,:,1);
    DCM.pC.C(:,:,2) = DCM.pC.C(:,:,1);
end

DCM.xY.Dfile            = megfilename;
DCM.xY.modality         = mode;
DCM.xY.name             = D.chanlabels(Ic);
DCM.xY.Ic               = Ic;
DCM.xY.Time             = ti;
DCM.xY.dt               = 1/D.fsample;
DCM.xY.coor2D           = D.coor2D(Ic);
DCM.xY.pst              = DCM.xY.Time(It);
DCM.xY.dt               = DCM.options.D/D.fsample;
DCM.xY.It               = It;
DCM.xY.code             = cond(DCM.options.trials);
DCM.xY.scale            = 1;
DCM.xY.Q                = {spm_Q(1/2,Ns,1)};

DCM.xU.X                = DM_design;
DCM.xU.name             = DM_names;
DCM.xU.dt               = DCM.xY.dt;

DCM.Sname               = DCM.xY.name';
DCM.Lpos                = double.empty(3,0);
DCM.name                = [megfilename(1:end-4) '_dcm' suffix '.mat'];
DCM.h                   = [];

DCM.M.dipfit.type       = DCM.options.spatial;
DCM.M.dipfit.location   = DCM.options.location;
DCM.M.dipfit.symmetry   = DCM.options.symmetry;
DCM.M.dipfit.model      = DCM.options.model;
DCM.M.dipfit.Ns         = length(DCM.Sname);
DCM.M.dipfit.Nc         = length(DCM.xY.Ic);
DCM.M.dipfit.modality   = 'LFP';

DCM.M.ons               = DCM.options.onset - DCM.xY.pst(1);
DCM.M.dur               = DCM.options.dur;

% disp(' -- spm_dcm_erp_data')

if DCM.options.h == 0
    X0 = sparse(Ns,1);
else X0 = spm_dctmtx(Ns,DCM.options.h);
end
R = speye(Ns) - X0*X0';

if DCM.options.han
    if Ns < 2048
        R = R*sparse(diag((tukeywin(Ns,.3)+.025)))*R; % R*sparse(diag(hanning(Ns)))*R; % TA altered this to be more inclusive of the epoch of interest
        disp('       Using TUKEY') % disp('       Using HANNING, not tukey') %
    else R = sparse(diag(Ns))*R;
    end
end

disp(['       Using -100ms -> 0ms, for baseline correction'])
tv0_i = [findnearest(-100,D.time) findnearest(0,D.time)];
for i = 1:length(DCM.options.trials)
    c = D.indtrial(cond(DCM.options.trials(i)), 'GOOD');
    Nt = length(c);
    DCM.xY.nt(i) = Nt;
    Y = zeros(Ns,Nc);
    for j = 1:Nt
        Y = Y + R*(these_data(Ic,It,c(j)) -  nanmean(these_data(Ic,tv0_i(1):tv0_i(2),c(j)),2))';
    end
    DCM.xY.y{i} = -Y/Nt;
end
disp(['       Y is -ve'])
DCM.xY.X0 = X0;
DCM.xY.R = R;
DCM.M.R = DCM.xY.R;

DCM.M.U        = spm_dcm_eeg_channelmodes(DCM.M.dipfit,Nm);

scale          = 30; % std(spm_vec(spm_fy_erp(DCM.xY.y,DCM.M)));       % TA altered this to stop the removal of amplitude differences across my groups
DCM.xY.y       = spm_unvec(spm_vec(DCM.xY.y)*scale,DCM.xY.y);
DCM.xY.scale   = DCM.xY.scale*scale;
disp('       Using my compromise scaling of 30') % disp('       Using original scaling') %

if strcmp('spm_gen_erp_nonx',DCM.M.IS)
    DCM.M.x = spm_dcm_neural_x_ignored(DCM.pE,DCM.M);
else DCM.M.x = spm_dcm_neural_x(DCM.pE,DCM.M);
end

DCM.M.pE       = DCM.pE;
DCM.M.pC       = DCM.pC;
DCM.M.gE       = DCM.gE;
DCM.M.gC       = DCM.gC;
DCM.M.m        = Nu;
DCM.M.n        = length(spm_vec(DCM.M.x));
DCM.M.l        = Nc;
DCM.M.ns       = Ns;
DCM.M.Nmax     = DCM.options.Nmax;


%%
disp(' -- RUNNING THE INVERSION')

DCM = spm_nlsi_N_TEC_modified(DCM);

function restore_env(old_path)
path(old_path);
end

end


