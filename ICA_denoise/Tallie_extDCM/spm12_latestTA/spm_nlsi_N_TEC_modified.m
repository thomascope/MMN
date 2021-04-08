function DCM = spm_nlsi_N_TEC_modified(DCM)
% Bayesian inversion of a linear-nonlinear model of the form F(p)*G(g)'
% FORMAT [Ep,Eg,Cp,Cg,S,F,L]= spm_nlsi_N(M,U,Y)
%
% Generative model
%__________________________________________________________________________
% 
% M.IS - IS(p,M,U) A prediction generating function name; usually an 
%        integration scheme for state-space models of the form
%
%        M.f  - f(x,u,p,M) - state equation:  dxdt = f(x,u)
%
%        that returns hidden states - x; however, it can be any nonlinear
%        function of the inputs u. I.e., x = IS(p,M,U)
%
% M.G  - G(g,M) - linear observer: y = (x - M.x')*G(g,M)'
%
% M.FS - function name f(y,M) - feature selection
%        This [optional] function performs feature selection assuming the
%        generalized model y = FS(y,M) = FS(x*G',M) + X0*P0 + e
%
% M.x  - The expansion point for the states (i.e., the fixed point)
%
% M.P  - starting estimates for model parameters [ states � optional]
% M.Q  - starting estimates for model parameters [ observer � optional]
%
% M.pE - prior expectation  - of model parameters - f(x,u,p,M)
% M.pC - prior covariance   - of model parameters - f(x,u,p,M)
%
% M.gE - prior expectation  - of model parameters - G(g,M)
% M.gC - prior covariance   - of model parameters - G(g,M)
%
% M.hE - prior expectation  - E{h}   of log-precision parameters
% M.hC - prior covariance   - Cov{h} of log-precision parameters
%
% U.u  - inputs
% U.dt - sampling interval
%
% Y.y  - {[ns,nx],...} - [ns] samples x [nx] channels x {trials}
% Y.X0 - Confounds or null space
% Y.dt - sampling interval for outputs
% Y.Q  - error precision components
%
%
% Parameter estimates
%--------------------------------------------------------------------------
% Ep  - (p x 1)         conditional expectation  E{p|y}
% Cp  - (p x p)         conditional covariance   Cov{p|y}
%
% Eg  - (p x 1)         conditional expectation  E{g|y}
% Cg  - (p x p)         conditional covariance   Cov{g|y}
%
% S   - (v x v)         [Re]ML estimate of error Cov{e(h)}
%
% log evidence
%--------------------------------------------------------------------------
% F   - [-ve] free energy F = log evidence = p(y|m)
% 
%     L(1) = - ey'*iS*ey/2;             accuracy of states
%     L(2) = - ep'*ipC*ep/2;            accuracy of parameters (f)
%     L(3) = - eg'*igC*eg/2;            accuracy of parameters (g)
%     L(4) = - eu'*iuC*eu/2;            accuracy of parameters (u)
%     L(5) = - eh'*ihC*eh/2;            accuracy of precisions (u)
%     L(6) = - ns*nr*log(8*atan(1))/2;  constant
%     L(7) = - nq*spm_logdet(S)/2;      precision
%     L(8) = spm_logdet(ibC*Cb)/2;      parameter complexity
%     L(9) = spm_logdet(ihC*Ch)/2;      precision complexity
%
%__________________________________________________________________________
% Returns the moments of the posterior p.d.f. of the parameters of a
% nonlinear model specified by IS(P,M,U) under Gaussian assumptions. Usually,
% IS would be an integrator of a dynamic MIMO input-state-output model 
%
%              dx/dt = f(x,u,p)
%              y     = G(g)*x  + X0*B + e
%
% The E-Step uses a Fisher-Scoring scheme and a Laplace
% approximation to estimate the conditional expectation and covariance of P
% If the free-energy starts to increase, a Levenberg-Marquardt scheme is
% invoked.  The M-Step estimates the precision components of e, in terms
% of [Re]ML point estimators of the log-precisions.
% An optional feature selection can be specified with parameters M.FS
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_nlsi_N.m 7143 2017-07-29 18:50:38Z karl $
 
% options
%--------------------------------------------------------------------------
% TA adapted inputs:
disp('    -- setup')
M = DCM.M;
U = DCM.xU;
Y = DCM.xY;

try, M.nograph; catch, M.nograph = 0;  end
try, M.Nmax;    catch, M.Nmax    = 64; end
try, M.Gmax;    catch, M.Gmax    = 8;  end
try, M.Hmax;    catch, M.Hmax    = 4;  end

% figure (unless disabled)
%--------------------------------------------------------------------------
% if ~M.nograph
%     %hf = fig;
%     Fsi = spm_figure('GetWin','SI');
% end
 
% check integrator
%--------------------------------------------------------------------------
try
    IS = M.IS;
catch
    IS = 'spm_int_U';
end
 
% check observer has not been accidentally specified
%--------------------------------------------------------------------------
try
    M = rmfield(M,'g');
end
 
% composition of feature selection and prediction (usually an integrator)
%--------------------------------------------------------------------------
if isfield(M,'FS')
 
    % FS(y,M)
    %----------------------------------------------------------------------
    try
        y  = feval(M.FS,Y.y,M);
        try
            FS = inline([M.FS '(y,M)'],'y','M');
        catch
            FS = M.FS;
        end
 
    % FS(y)
    %----------------------------------------------------------------------
    catch
        y  = feval(M.FS,Y.y);
        FS = inline([M.FS '(y)'],'y','M');
 
    end
else
 
    % y
    %----------------------------------------------------------------------
    y  = Y.y;
    FS = inline('y','y','M');
end
 
% size of data (usually samples x channels)
%--------------------------------------------------------------------------
if iscell(y)
    
    % concatenate samples over cell, ensuring the same for predictions
    %----------------------------------------------------------------------
    ns = size(y{1},1);
    y  = spm_cat(y(:));             
    IS = inline(['spm_cat(' IS '(P,M,U))'],'P','M','U');
    
else
    ns = size(y,1);
end
ny   = length(spm_vec(y));
nr   = ny/ns;                            % number of samples and responses
M.ns = ns;                               % store in M.ns for integrator
 
% initial states
%--------------------------------------------------------------------------
try
    M.x;
catch
    try
        M.n;
    catch
        M.n = 0;
    end
    M.x = sparse(M.n,1);
end
 
% input
%--------------------------------------------------------------------------
try
    U;
catch
    U = [];
end
 
% initial parameters
%--------------------------------------------------------------------------
try
    spm_vec(M.P) - spm_vec(M.pE);
    if ~M.nograph
        fprintf('\n(state) parameter initialisation successful\n')
    end
catch
    M.P = M.pE;
end
try
    spm_vec(M.Q) - spm_vec(M.gE);
    if ~M.nograph
        fprintf('\n(observer) parameter initialisation successful\n')
    end
catch
    M.Q = M.gE;
end
 
 
% time-step
%--------------------------------------------------------------------------
try
    Y.dt;
catch
    Y.dt = 1;
end
 
 
% precision components Q
%--------------------------------------------------------------------------
try
    Q = Y.Q;
catch
    Q = spm_Ce(ns*ones(1,nr));
end
nh    = length(Q);          % number of precision components
nt    = length(Q{1});       % number of time bins
nq    = nr*ns/nt;           % for compact Kronecker form of M-step

 
% confounds (if specified)
%--------------------------------------------------------------------------
try
    if isempty(Y.X0)
        Y.X0 = sparse(ns,0);
    end
    dgdu = kron(speye(nr,nr),Y.X0);
catch
    dgdu = sparse(ns*nr,0);
end
 
% hyperpriors - expectation (and initialize hyperparameters)
%--------------------------------------------------------------------------
try
    hE  = M.hE;
    if length(hE) ~= nh
        hE = hE + sparse(nh,1);
    end
catch
    hE  = sparse(nh,1) - log(var(spm_vec(y))) + 4;
end
h       = hE;

% hyperpriors - covariance
%--------------------------------------------------------------------------
try
    ihC = spm_inv(M.hC);
    if length(ihC) ~= nh
        ihC = ihC*speye(nh,nh);
    end
catch
    ihC = speye(nh,nh)*exp(4);
end


% unpack prior covariances
%--------------------------------------------------------------------------
if isstruct(M.pC); M.pC = spm_diag(spm_vec(M.pC)); end
if isstruct(M.gC); M.gC = spm_diag(spm_vec(M.gC)); end
if isvector(M.pC); M.pC = spm_diag(M.pC); end
if isvector(M.gC); M.gC = spm_diag(M.gC); end

% dimension reduction of parameter space
%--------------------------------------------------------------------------
Vp    = spm_svd(M.pC,0);
Vg    = spm_svd(M.gC,0);
np    = size(Vp,2);                   % number of parameters (f)
ng    = size(Vg,2);                   % number of parameters (g)
nu    = size(dgdu,2);                 % number of parameters (u)

 
% prior moments
%--------------------------------------------------------------------------
pE    = M.pE;
gE    = M.gE;
uE    = sparse(nu,1);
 
% second-order moments (in reduced space)
%--------------------------------------------------------------------------
sw    = warning('off','all');
pC    = Vp'*M.pC*Vp;
gC    = Vg'*M.gC*Vg;
uC    = speye(nu,nu)*exp(16);
ipC   = spm_inv(pC);                           % p - state parameters
igC   = spm_inv(gC);                           % g - observer parameters
iuC   = spm_inv(uC);                           % u - fixed parameters
ibC   = spm_cat(spm_diag({ipC,igC,iuC}));      % all parameters
bC    = speye(size(ibC))*exp(-16);

 
% initialize conditional density
%--------------------------------------------------------------------------
Ep    = M.P;
Eg    = M.Q;
Eu    = spm_pinv(dgdu)*spm_vec(y);

% expansion point
%--------------------------------------------------------------------------
if ~isempty(M.x)
    x0 = ones(size(y,1),1)*spm_vec(M.x)';
else
    x0 = 0;
end


% EM
%==========================================================================
warning(sw); sw = warning('off','all');
criterion       = [0 0 0 0];

C.F   = -Inf;                                   % free energy
v     = -4;                                     % log ascent rate
dgdp  = zeros(ny,np);
dgdg  = zeros(ny,ng);
dFdh  = zeros(nh,1);
dFdhh = zeros(nh,nh);

fprintf('.      Vp: [%s] \n.      Vg: [%s] \n',num2str(size(Vp)),num2str(size(Vg)))
 
% Optimize p: parameters of f(x,u,p)
%==========================================================================
strR = nout(2,@fileparts,DCM.name);
strRa = datestr(datetime);
strRa = strRa(1:end-3);
mkdir(['/imaging/mlr/users/tc02/Holly_MMN/DCM_REPORTS/REPORTS_' strRa])
disp(['    -- reporting incrementally saved in /imaging/mlr/users/tc02/Holly_MMN/DCM_REPORTS/REPORTS_' strRa])
%disp(['    -- EM optimize p ' strR])
EP     = [];
for ip = 1:M.Nmax
 
    % time
    %----------------------------------------------------------------------  
    Ti = tic;
    
    % predicted hidden states (x) and dxdp
    %----------------------------------------------------------------------
    [dxdp,x] = spm_diff(IS,Ep,M,U,1,{Vp}); % spm_gen_erp + spm_int_L + [spm_fx_cmcTA OR spm_fx_cmm_NMDA]
    % DCM.M.Xp = x(51,:); % TA removed - historical test from me
    % check for inital iterations and dissipative dynamics
    %----------------------------------------------------------------------
    if all(isfinite(spm_vec(x)))
        Gmax = M.Gmax;
        if ip < 8
            vg = -4;
        else
            vg = 2;
        end
    else
        Gmax = 0;
    end
    
       
    % Optimize g: parameters of G(g)
    %======================================================================
    %if ip==1 || 2, disp(['    -- EM optimize g ' datestr(datetime)]), end
    for ig = 1:Gmax
        
        % prediction yp = G(g)*x
        %------------------------------------------------------------------
        [dGdg,G] = spm_diff(M.G,Eg,M,1,{Vg}); % spm_erp_L -> spm_lx_erp (your symmetrical build of L based on your inputted J's)
        yp       = FS((x - x0)*G',M); % spm_fy_erp
        
        % prediction errors - states
        %==================================================================
        ey    = spm_vec(y)  - spm_vec(yp) - dgdu*Eu;
 
        % prediction errors - parameters
        %------------------------------------------------------------------
        ep    = Vp'*(spm_vec(Ep) - spm_vec(pE));
        eg    = Vg'*(spm_vec(Eg) - spm_vec(gE));
        eu    =      spm_vec(Eu) - spm_vec(uE);
 
        % gradients
        %------------------------------------------------------------------
        for i = 1:np
            dgdp(:,i) = spm_vec(FS(dxdp{i}*G',M));
        end
        try
            for i = 1:ng
                dgdg(:,i) = spm_vec(FS((x - x0)*dGdg{i}',M));
            end
        catch
            dgdg = FS((x - x0)*dGdg,M);
        end
  
        % Optimize F(h): parameters of iS(h)
        %==================================================================   
        %if ip==1 || 2, disp(['    -- EM optimize h ' datestr(datetime)]), end
        dgdb   = [dgdp dgdg dgdu];           
        for ih = 1:M.Hmax
 
            % precision
            %--------------------------------------------------------------
            iS    = speye(nt,nt)*exp(-32);
            for i = 1:nh
                iS = iS + Q{i}*exp(h(i));
            end
            S     = spm_inv(iS);
            iS    = kron(speye(nq),iS);
            dFdbb = dgdb'*iS*dgdb + ibC;
            Cb    = spm_inv(dFdbb) + bC;
            
            % precision operators for M-Step
            %--------------------------------------------------------------
            for i = 1:nh
                P{i}  = Q{i}*exp(h(i));
                PS{i} = P{i}*S;
                P{i}  = kron(speye(nq),P{i});
            end
 
            % derivatives: dLdh = dL/dh,...
            %--------------------------------------------------------------
            for i = 1:nh
                dFdh(i,1)      =   trace(PS{i})*nq/2 ...
                                 - real(ey'*P{i}*ey)/2 ...
                                 - spm_trace(Cb,dgdb'*P{i}*dgdb)/2;
                for j = i:nh
                    dFdhh(i,j) = - spm_trace(PS{i},PS{j})*nq/2;
                    dFdhh(j,i) =   dFdhh(i,j);
                end
            end
 
            % add hyperpriors
            %--------------------------------------------------------------
            eh    = h     - hE;
            dFdh  = dFdh  - ihC*eh;
            dFdhh = dFdhh - ihC;
            Ch    = spm_inv(-dFdhh);
            
            % M-Step: update ReML estimate of h
            %--------------------------------------------------------------
            dh    = spm_dx(dFdhh,dFdh,{4});
            h     = h + min(max(dh,-2),2);
 
            % convergence
            %--------------------------------------------------------------
            if dFdh'*dh < exp(-2), break, end
 
        end
 
        % E-step: optimise F(g,u)
        %==================================================================
        % update gradients and curvature - counfounds
        %------------------------------------------------------------------
        dFdu  =  dgdu'*iS*ey   - iuC*eu;
        dFduu = -dgdu'*iS*dgdu - iuC;

        % Conditional updates of confounds (u)
        %------------------------------------------------------------------
        du    = spm_dx(dFduu,dFdu,{4});
        Eu    = Eu + du;
        
        % update gradients and curvature - parameters
        %------------------------------------------------------------------
        dFdg  =  dgdg'*iS*ey   - igC*eg;
        dFdgg = -dgdg'*iS*dgdg - igC;
 
        % Conditional updates of parameters (g)
        %------------------------------------------------------------------
        dg    = spm_dx(dFdgg,dFdg,{vg});
        Eg    = spm_unvec(spm_vec(Eg) + Vg*dg,Eg);
         
        % convergence
        %------------------------------------------------------------------
        dG    = dFdg'*dg;
        if ig > 1 && dG < exp(-2), break, end
        
    end
    
    % optimise objective function: F(p) = log-evidence - divergence
    %======================================================================
    L(1) = - ey'*iS*ey/2;            % accuracy
    L(2) = - ep'*ipC*ep/2;           % complexity
    L(3) = - eg'*igC*eg/2;           % complexity
    L(4) = - eu'*iuC*eu/2;           % complexity
    L(5) = - eh'*ihC*eh/2;           % complexity
    L(6) = - ns*nr*log(8*atan(1))/2; % accuracy (2*pi)
    L(7) = - nq*spm_logdet(S)/2;     % accuracy
    L(8) = spm_logdet(ibC*Cb)/2;     % complexity
    L(9) = spm_logdet(ihC*Ch)/2;     % complexity
    F    = sum(L);
    
    % record increases and reference log-evidence for reporting
    %----------------------------------------------------------------------
    try
        F0;
        fprintf(' actual: %.3e (%.2f sec)\n',full(F - C.F),toc(Ti))
    catch
        F0 = F;
    end
     
    % if F has increased, update gradients and curvatures for E-Step
    %----------------------------------------------------------------------
    if F > C.F || ip < 4
        
        % update gradients and curvature
        %------------------------------------------------------------------
        dFdp  =  dgdp'*iS*ey   - ipC*ep;
        dFdpp = -dgdp'*iS*dgdp - ipC;
 
        % decrease regularization
        %------------------------------------------------------------------
        v     = min(v + 1/2,4);
        str   = 'EM(+)';
 
        % accept current estimates
        %------------------------------------------------------------------
        C.Cb  = Cb;                               % conditional covariance
        C.Ep  = Ep;                               % and expectations
        C.Eg  = Eg;
        C.Eu  = Eu;
        C.h   = h;
        C.F   = F;
        C.L   = L;   
        
    else
 
        % reset expansion point
        %------------------------------------------------------------------
        Cb    = C.Cb;                             % conditional covariance
        Ep    = C.Ep;                             % and expectations
        Eg    = C.Eg;
        Eu    = C.Eu;
        h     = C.h;
 
        % and increase regularization
        %------------------------------------------------------------------
        v     = min(v - 2,-4);
        str   = 'EM(-)';
 
    end
 
    % Optimize p: parameters of f(x,u,p)
    %======================================================================
    dp    = spm_dx(dFdpp,dFdp,{v});
    Ep    = spm_unvec(spm_vec(Ep) + Vp*dp,Ep);
    
    % diagnostic
    %----------------------------------------------------------------------
    % EP(:,end + 1) = spm_vec(Ep);
 
    
    % subplot times
    %----------------------------------------------------------------------
    try
        if length(Y.pst) == size(yp,1)
            yt = Y.pst;
        else
            yt = (1:size(yp,1))*Y.dt*1000;
        end
    catch
        yt = (1:size(yp,1))*Y.dt*1000;
    end
    
    
    % graphics
    %----------------------------------------------------------------------

        % convergence
        %----------------------------------------------------------------------
        dF  = dFdp'*dp;
        ig  = max([0 ig]);
        fprintf('%-6s: %-2i (%i,%i) %4s %-6.3e %6s %6.3e  %s ',str,ip,ig,ih,'F:',full(C.F - F0),'dF predicted:',full(dF),strR)
        criterion = [(dF < 1e-1) criterion(1:end - 1)];
        R.F(ip) = full(C.F); R.F0(ip) = full(F0); R.dF(ip) = full(dF);
        R.x(:,:,ip) = full(x); R.yp(:,:,ip) = full(yp); R.y(:,:,ip) = full(cat(1,DCM.xY.y{:}));
        R.plab = {'A1' 'A2' 'AN1' 'AN2' 'B' 'BN' 'H E' 'H IA' 'HIB' 'C, CV, D R S' 'G' 'T' 'J' 'L'};
        R.p{1}(:,:,ip) = Ep.A{1};
        R.p{2}(:,:,ip) = Ep.A{2};
        R.p{3}(:,:,ip) = Ep.AN{1};
        R.p{4}(:,:,ip) = Ep.AN{2};
        R.p{5}(:,:,ip) = Ep.B{1};
        R.p{6}(:,:,ip) = Ep.BN{1};
        R.p{7}(:,:,ip) = Ep.H(:);%(:,:,1);
        R.p{8}(:,:,ip) = Ep.H(:);%(:,:,2);
        R.p{9}(:,:,ip) = Ep.H(:);%(:,:,3);
        R.p{10}(:,:,ip) = zeros(1,1);%[Ep.C Ep.CV' [Ep.D';Ep.R';Ep.S]];
        R.p{11}(:,:,ip) = Ep.G(:);
        R.p{12}(:,:,ip) = Ep.T;
        R.p{13}(:,:,ip) = full(Eg.J);
        R.p{14}(:,:,ip) = Eg.L;
        
        R.p{10} = permute(R.p{10},[2 1 3]);
        R.p{11} = permute(R.p{11},[2 1 3]);
        
        save(['/imaging/mlr/users/tc02/Holly_MMN/DCM_REPORTS/REPORTS_' strRa '/' strR],'R')
        
        %NODES = {'LIFG' 'LSTG' 'LAud' 'RIFG' 'RSTG' 'RAud'};
        NODES = DCM.Sname';
        POPS = {'ss' 'sp' 'si' 'dp' 'di' 'tp'};
        lab =  {{'AMPA connectivity';'Forward'};...
                {'AMPA connectivity';'Backward'};...
                {'NMDA connectivity';'Forward'};...
                {'NMDA connectivity';'Backward'};...
                {'Trial-dependent';'AMPA difference'};...
                {'Trial-dependent';'NMDA difference'};...
                {'NMDA connectivity';'Intrinsic'};...
                {'GABAA connectivity';'Intrinsic'};...
                {'GABAB connectivity';'Intrinsic'};...
                '(left overs)';...
                {'Trial-dependent';'intrinsic differences'};...
                'Time Constants';...
                'Hyperprior J';...
                'Hyperprior L'};
        labx = {NODES NODES NODES NODES NODES NODES POPS POPS POPS {''} [strcat('self .',POPS) {'si.ss' 'si.sp' 'di.dp' 'di.tp'}] {'AMPA' 'NMDA' 'GABAA' 'GABAB'} '' ''};
        laby = {NODES NODES NODES NODES NODES NODES POPS POPS POPS {'C' 'CV' 'D1 D2 D3 R1 R2 S'} {'AMPA' 'NMDA' 'GABAA' 'GABAB'} NODES '' ''};
        these_times = DCM.xY.Time(DCM.xY.Time>=DCM.options.Tdcm(1)&DCM.xY.Time<=DCM.options.Tdcm(2));
        
        sp = {1 2 3 4 5 6 8 9 10 [] [11 12] 7};
        cm = colormap_tor([0 0 1],[1 0 0],[0 1 1],[1 1 1],[1 1 0]);
        cm((end/2)-2:(end/2)+2,:) = [];
        clr = [repmat([1 0 0],6,1) ; repmat([0 1 0],6,1) ; repmat([0 0 1],6,1) ; repmat([.7 .7 0],6,1) ; repmat([0 1 1],6,1) ; repmat([1 0 1],6,1)];
        node_colours = distinguishable_colors(length(NODES));
        fm; % set(gcf,'Color','w')
        k = ip;
        for k1 = [1:9 11:12]
            Z = R.p{k1}(:,:,k);
            Z(Z==-32) = 0;
            subplot(3,7,sp{k1}), surf([Z zeros(size(R.p{k1},1),1) ; zeros(1,size(R.p{k1},2)+1)],'EdgeColor','w','LineWidth',2), view(2)
            set(gca,'YDir','reverse')
            set(gca,'XTick',[1:length(labx{k1})]+.5,'XTickLabel',labx{k1},'XTickLabelRotation',90,'YTick',[1:length(laby{k1})]+.5,'YTickLabel',laby{k1})
            title(lab{k1}), caxis([-1 1]), colormap(cm), axis tight equal
            if k1==11, colorbar, end
        end

        subplot(3,7,[15 16]), hold on
        for k1 = [1:3:36], plot(these_times,squeeze(R.x(:,k1,k)),'Color',clr(k1,:)), end
        xlabel('time (ms)'), axis tight, ylim([-90 10]), grid on
        title('IFG pop voltages')
        subplot(3,7,[17 18]), hold on
        for k1 = [2:3:36], plot(these_times,squeeze(R.x(:,k1,k)),'Color',clr(k1,:)), end
        xlabel('time (ms)'), axis tight, ylim([-90 10]), grid on
        title('STG pop voltages')
        subplot(3,7,[19 20]), hold on
        for k1 = [3:3:36], plot(these_times,squeeze(R.x(:,k1,k)),'Color',clr(k1,:)), end
        xlabel('time (ms)'), axis tight, ylim([-90 10]), grid on
        title('Aud pop voltages')
        legend({'L ss' 'R ss' 'L sp' 'R sp' 'L si' 'R si' 'L dp' 'R dp' 'L di' 'R di' 'L tp' 'R tp'},'FontSize',6,'Position',[0.8 0.15 0.04 0.14])

        subplot(3,7,[13 14]), hold on
        for k1 = 1:size(R.yp,2)
            if DCM.options.han
                plot(these_times,(R.yp(:,k1,k)-nanmean(R.yp(:,k1,k),1)).*[(tukeywin(ceil(size(R.yp,1)/2),.3)+(1/3))/(4/3)';(tukeywin(floor(size(R.yp,1)/2),.3)+(1/3))/(4/3)'],'Color',node_colours(k1,:),'LineWidth',2.5)
            else
                plot(these_times,(R.yp(:,k1,k)-nanmean(R.yp(:,k1,k),1)),'Color',node_colours(k1,:),'LineWidth',2.5)
            end % [hanning(size(R.yp,1)/2)' hanning(size(R.yp,1)/2)']'
            plot(these_times,R.y(:,k1,k),'Color',node_colours(k1,:))
        end
        axis tight, xlabel('time (ms)'), title('data LFP & predicted LFP'), grid on
        this_legend = {};
        for this_node = 1:length(NODES)
            this_legend = [this_legend, ['yp ' NODES{this_node}], ['y ' NODES{this_node}]];
        end
        legend(this_legend,'FontSize',6,'Position',[0.92 0.44 0.05 0.14])

        subplot(9,14,112)
        ax = plotyy(1:k,R.F(1:k),1:k,R.dF(1:k));
        xlabel(ax(2),'iterations'), ylabel(ax(1),'F'), ylabel(ax(2),'dF'), set(ax,'YTick','')
        xlim(ax,[0 length(R.F)]), try ylim(ax(1),[min(R.F) max(R.F)]), ylim(ax(2),[0 max(R.dF)]), end

        drawnow, close(gcf)
        
        R.p{10} = permute(R.p{10},[2 1 3]);
        R.p{11} = permute(R.p{11},[2 1 3]);
 
    if all(criterion), fprintf(' convergence\n'), break, end
    
end
% if exist('Fsi', 'var')
%     spm_figure('Focus', Fsi)
% end

% outputs
%--------------------------------------------------------------------------
DCM.Ce     = S;
DCM.Ep     = C.Ep;
DCM.Eg     = C.Eg;
DCM.Cp     = Vp*C.Cb((1:np),     (1:np)     )*Vp';
DCM.Cg     = Vg*C.Cb((1:ng) + np,(1:ng) + np)*Vg';
DCM.F      = C.F;
DCM.L      = C.L;
warning(sw);

% diagnostic
%--------------------------------------------------------------------------
% save('spm_nlsi_N_Ep','EP')

% TA added this last section from the bottom of spm_dcm_erp:

Ns                      = size(DCM.xY.y{1},1);                      % number of samples
Nr                      = size(DCM.C,1);                                % number of sources

DCM.ID = spm_data_id(feval(M.FS,DCM.xY.y,M));

% Bayesian inference
sw  = warning('off','SPM:negativeVariance');
dp  = spm_vec(DCM.Ep) - spm_vec(DCM.pE);
DCM.Pp  = spm_unvec(1 - spm_Ncdf(0,abs(dp),diag(DCM.Cp)),DCM.Ep);
warning(sw);

% neuronal and sensor responses (x and y)
L       = feval(DCM.M.G,DCM.Eg,DCM.M);                 % get gain matrix
DCM.x   = feval(DCM.M.IS,DCM.Ep,DCM.M,DCM.xU);         % prediction (source space)

% trial-specific responses (in mode, channel and source space)
try j = find(kron(DCM.Eg.J,ones(1,Nr)));               % Indices of contributing states
catch err
    j = find(spm_cat(DCM.Eg.J));
end
x0  = ones(Ns,1)*spm_vec(DCM.M.x)';                    % expansion point for states
for i = 1:length(DCM.xY.y)
    DCM.K{i} = DCM.x{i} - x0;                          % centre on expansion point
    DCM.H{i} = DCM.M.R*DCM.K{i}*L'*DCM.M.U;            % prediction (sensor space)
    DCM.R{i} = DCM.M.R*DCM.xY.y{i}*DCM.M.U - DCM.H{i}; % residuals  (sensor space)
    DCM.K{i} = DCM.K{i}(:,j);                          % Depolarization in sources
end
DCM.M = rmfield(DCM.M,'dipfit');
if ~exist('/imaging/mlr/users/tc02/Holly_MMN/extDCMs/','dir')
    mkdir('/imaging/mlr/users/tc02/Holly_MMN/extDCMs/')
end
save(['/imaging/mlr/users/tc02/Holly_MMN/extDCMs/' nout(2,@fileparts,DCM.name) '_' DCM.xY.code{1} '.mat'],'DCM')

return
