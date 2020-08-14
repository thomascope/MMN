function DCM = spm_dcm_csd_as(DCM)
% Estimate parameters of a DCM of (complex) cross-spectral density
% FORMAT DCM = spm_dcm_csd(DCM)
%
% DCM
%    name: name string
%       xY: data   [1x1 struct]
%       xU: design [1x1 struct]
%
%   Sname: cell of source name strings
%       A: {[nr x nr double]  [nr x nr double]  [nr x nr double]}
%       B: {[nr x nr double], ...}   Connection constraints
%       C: [nr x 1 double]
%
%   options.Nmodes       - number of spatial modes
%   options.Tdcm         - [start end] time window in ms
%   options.Fdcm         - [start end] Frequency window in Hz
%   options.D            - time bin decimation       (usually 1 or 2)
%   options.spatial      - 'ECD', 'LFP' or 'IMG'     (see spm_erp_L)
%   options.model        - 'ERP', 'SEP', 'CMC', 'LFP', 'NMM' or 'MFM'
%
% Esimates:
%--------------------------------------------------------------------------
% DCM.dtf                   - directed transfer functions (source space)
% DCM.ccf                   - cross covariance functions (source space)
% DCM.coh                   - cross coherence functions (source space)
% DCM.fsd                   - specific delay functions (source space)
% DCM.pst                   - peristimulus time
% DCM.Hz                    - frequency
%
% DCM.Ep                    - conditional expectation
% DCM.Cp                    - conditional covariance
% DCM.Pp                    - conditional probability
% DCM.Hc                    - conditional responses (y), channel space
% DCM.Rc                    - conditional residuals (y), channel space
% DCM.Hs                    - conditional responses (y), source space
% DCM.Ce                    - eML error covariance
% DCM.F                     - Laplace log evidence
% DCM.ID                    -  data ID
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_dcm_csd.m 5975 2014-05-07 18:07:42Z karl $
 
 
% check options
%==========================================================================
drawnow
clear spm_erp_L
name = sprintf('DCM_%s',date);
DCM.options.analysis  = 'CSD';
 
% Filename and options
%--------------------------------------------------------------------------
try, DCM.name;                      catch, DCM.name = name;      end
try, model   = DCM.options.model;   catch, model    = 'NMM';     end
try, spatial = DCM.options.spatial; catch, spatial  = 'LFP';     end
try, Nm      = DCM.options.Nmodes;  catch, Nm       = 8;         end
try, DATA    = DCM.options.DATA;    catch, DATA     = 1;         end
 
% Spatial model
%==========================================================================
DCM.options.Nmodes = Nm;
DCM.M.dipfit.model = model;
DCM.M.dipfit.type  = spatial;

if ~DATA || ~DCM.options.DONE
    DCM  = spm_dcm_csd_data_as(DCM);                   % data
end
if DATA
    DCM  = spm_dcm_erp_dipfit(DCM, 1);              % spatial model
end
Ns   = length(DCM.A{1});                            % number of sources


% Design model and exogenous inputs
%==========================================================================
if ~isfield(DCM,'xU'),   DCM.xU.X = sparse(1 ,0); end
if ~isfield(DCM.xU,'X'), DCM.xU.X = sparse(1 ,0); end
if ~isfield(DCM,'C'),    DCM.C    = sparse(Ns,0); end
if isempty(DCM.xU.X),    DCM.xU.X = sparse(1 ,0); end
if isempty(DCM.xU.X),    DCM.C    = sparse(Ns,0); end

% Neural mass model
%==========================================================================
 
% prior moments on parameters
%--------------------------------------------------------------------------
[pE,pC]  = spm_dcm_neural_priors(DCM.A,DCM.B,DCM.C,model);
  
% check to see if neuronal priors have already been specified
%--------------------------------------------------------------------------
try
    if length(spm_vec(DCM.M.pE)) == length(spm_vec(pE));
        pE = DCM.M.pE;
        pC = DCM.M.pC;
        fprintf('Using existing priors\n')
    end
end
 
% augment with priors on spatial model
%--------------------------------------------------------------------------
[pE,pC]  = spm_L_priors(DCM.M.dipfit,pE,pC);
 
% augment with priors on endogenous inputs (neuronal) and noise
%--------------------------------------------------------------------------
[pE,pC]  = spm_ssr_priors(pE,pC);
pE.J = DCM.CUSTOM.gE.J;
pC.J = DCM.CUSTOM.gC.J;
 
% initial states and equations of motion
%--------------------------------------------------------------------------
[x,f]    = spm_dcm_x_neural(pE,model,DCM.CUSTOM.f);
 
% create DCM
%--------------------------------------------------------------------------
DCM.M.IS = 'spm_csd_mtf';
DCM.M.FS = 'spm_fs_csd';
DCM.M.g  = 'spm_gx_erp';
DCM.M.f  = f;
DCM.M.x  = x;
DCM.M.n  = length(spm_vec(x));
DCM.M.pE = pE;
DCM.M.pC = pC;
DCM.M.hE = 6;
DCM.M.hC = 1/64;
DCM.M.m  = Ns;
DCM.xY.Hz = 1:48;

% specify M.u - endogenous input (fluctuations) and intial states
%--------------------------------------------------------------------------
DCM.M.u  = sparse(Ns,1);

%-Feature selection using principal components (U) of lead-field
%==========================================================================
 
% Spatial modes
%--------------------------------------------------------------------------
try
    DCM.M.U  = spm_dcm_eeg_channelmodes(DCM.M.dipfit,Nm);
end
 

% Alex modify: Check for CUSTOM functions and priors and insert
%--------------------------------------------------------------------------
Str = 'Checking for custom parameter priors.....\n';
Str = [Str , repmat('-',[1 length(Str)-2]) '\n' ]; fprintf(Str);
try CUS   = DCM.CUSTOM;
    PARAM = fieldnames(CUS.pE);
    TARGT = fieldnames(pE);
    for i = 1:length(PARAM)
        if ismember(PARAM{i},TARGT)
            DCM.M.pE.(PARAM{i}) = CUS.pE.(PARAM{i});
            fprintf('Using custom prior parameter: %s\n',(PARAM{i}));
        end
    end
end
try CUS   = DCM.CUSTOM;
    PARAM = fieldnames(CUS.pC);
    TARGT = fieldnames(pC);
    for i = 1:length(PARAM)
        if ismember(PARAM{i},TARGT)
            DCM.M.pC.(PARAM{i}) = CUS.pC.(PARAM{i});
            fprintf('Using custom prior variance parameter: %s\n',(PARAM{i}));
        end
    end
end    
    
try DCM.M.f    = CUS.f;   
    fprintf('\nCustom state equations %s',CUS.f);
end

try DCM.M.IS    = CUS.IS;   
    fprintf('\nCustom transfer (noise) function %s\n\n',CUS.IS);
end


% re-initialise states given new priors
[x,f]   = spm_dcm_x_neural(DCM.M.pE,model,DCM.CUSTOM.f);
DCM.M.x = x;

% get data-features (in reduced eigenspace)
%==========================================================================
if ~DATA
   % DCM  = spm_dcm_csd_data(DCM);
   try fdata = DCM.options.data;
   catch
       fdata = 'spm_dcm_csd_data_as';
   end
   fdata
    %DCM  = spm_dcm_csd_data_as(DCM);
    DCM = feval(fdata,DCM);
end
 
% scale data features (to a variance of about 8)
%--------------------------------------------------------------------------
ccf      = spm_csd2ccf(DCM.xY.y,DCM.xY.Hz);
scale    = max(spm_vec(ccf));
DCM.xY.y = spm_unvec(8*spm_vec(DCM.xY.y)/scale,DCM.xY.y);



% complete model specification and invert
%==========================================================================
Nm       = size(DCM.M.U,2);                    % number of spatial modes
DCM.M.l  = Nm;
DCM.M.Hz = DCM.xY.Hz;
DCM.M.dt = DCM.xY.dt;
 
% precision of noise: AR(1/2)
%--------------------------------------------------------------------------
% y     = spm_fs_csd(DCM.xY.y,DCM.M);
% DCM.M.FS = 'spm_fs_csd';


% 
% Alex new features selection - just abs(csd)
%--------------------------------------------------------------------------
FS = inline('spm_unvec(real(spm_vec(x)),x)','x','y');
y  = feval(FS,DCM.xY.y,DCM.M);
DCM.M.FS = FS;
% 
% FS = inline('spm_unvec(log(spm_vec(x)),x)','x','y');
% y  = feval(FS,DCM.xY.y,DCM.M);
% DCM.M.FS = FS;

% y     = spm_fs_csd_as(DCM.xY.y,DCM.M);
% DCM.M.FS = 'spm_fs_csd_as';


for i = 1:length(y)
    n      = size(y{i},1);
    m      = size(y{i},2)*size(y{i},3);
    q      = spm_Q(1/2,n,1);  
    
    % alex - focus on higher frequencies: Q*diag(Hz)*Q
    if  n == length(DCM.M.Hz)
        H = DCM.M.Hz;
    else
        H = linspace(DCM.M.Hz(1),DCM.M.Hz(2),n);
    end
    
    q = spm_Q(1/2,n,1)*diag(H)*spm_Q(1/2,n,1); % as per ssr gamma
    Q{i,i} = kron(speye(m,m),q);
end
DCM.xY.Q  = spm_cat(Q);
DCM.xY.X0 = sparse(size(Q,1),0);



try NGP = DCM.CUSTOM.nograph;
    if NGP
        DCM.M.nograph = 1;
    end
end
try DCM.M.ShowNodes = CUS.ShowNodes; end

% Variational Laplace: model inversion
%==========================================================================
[Qp,Cp,Eh,F] = spm_nlsi_GN(DCM.M,DCM.xU,DCM.xY);


% Two-step fit
%--------------------------------------------------------------------------
try DCM.options.TwoStep; catch; DCM.options.TwoStep = 0; end

if  DCM.options.TwoStep
    close;
    fprintf('Implementing two-step fitting procedure...\n');
    DCM.M.pE = Qp;
    DCM.M.pC = spm_unvec( diag(Cp)*0, DCM.M.pC);

    DCM.M.pE.A{1} = ones(size(DCM.M.pE.A{1}));
    DCM.M.pE.A{2} = ones(size(DCM.M.pE.A{2}));

    DCM.M.pC.A{1} = DCM.M.pC.A{1} + 2;
    DCM.M.pC.A{2} = DCM.M.pC.A{2} + 2;

    DCM.M.pC.H = DCM.M.pC.H + 2;
    DCM.M.pC.a = DCM.M.pC.a + 1/8;
    DCM.M.pC.b = DCM.M.pC.b + 1/8;
    %DCM.M.pC.T = DCM.M.pC.T + 2;

    [Qp,Cp,Eh,F] = spm_nlsi_GN(DCM.M,DCM.xU,DCM.xY);
    close;
end


% Draw [spatial] prediction and error [Alex add]
%--------------------------------------------------------------------------
% Hc = feval(DCM.M.IS,Qp,DCM.M,DCM.xU);
% PCSD(DCM.xY.y{:},Hc{:},'abs',DCM.xY.Hz);drawnow
% PCSD(abs(DCM.xY.y{:})-abs(Hc{:}),[],'abs',DCM.xY.Hz,'col','--k');drawnow


% Data ID
%--------------------------------------------------------------------------
try
    try
        ID = spm_data_id(feval(DCM.M.FS,DCM.xY.y,DCM.M));
    catch
        ID = spm_data_id(feval(DCM.M.FS,DCM.xY.y));
    end
catch
    ID = spm_data_id(DCM.xY.y);
end
 
pE = DCM.M.pE;
pC = DCM.M.pC;

% Bayesian inference {threshold = prior} NB Prior on A,B and C = exp(0) = 1
%==========================================================================
warning('off','SPM:negativeVariance');
dp  = spm_vec(Qp) - spm_vec(pE);
Pp  = spm_unvec(1 - spm_Ncdf(0,abs(dp),diag(Cp)),Qp);
warning('on', 'SPM:negativeVariance');
 
 
% predictions (csd) and error (sensor space)
%--------------------------------------------------------------------------
Hc  = spm_csd_mtf_as(Qp,DCM.M,DCM.xU);                   % prediction
Ec  = spm_unvec(spm_vec(DCM.xY.y) - spm_vec(Hc),Hc);     % prediction error
 
 
% predictions (source space - cf, a LFP from virtual electrode)
%--------------------------------------------------------------------------
M             = rmfield(DCM.M,'U'); 
M.dipfit.type = 'LFP';

M.U         = 1; 
M.l         = Ns;
qp          = Qp;
qp.L        = ones(1,Ns);             % set virtual electrode gain to unity
qp.b        = qp.b - 32;              % and suppress non-specific and
qp.c        = qp.c - 32;              % specific channel noise

[Hs Hz dtf] = spm_csd_mtf_as(qp,M,DCM.xU);
[ccf pst]   = spm_csd2ccf(Hs,DCM.M.Hz);
[coh fsd]   = spm_csd2coh(Hs,DCM.M.Hz);
DCM.dtf     = dtf;
DCM.ccf     = ccf;
DCM.coh     = coh;
DCM.fsd     = fsd;
DCM.pst     = pst;
DCM.Hz      = Hz;

 
% store estimates in DCM
%--------------------------------------------------------------------------
DCM.Ep = Qp;                   % conditional expectation
DCM.Cp = Cp;                   % conditional covariance
DCM.Pp = Pp;                   % conditional probability
DCM.Hc = Hc;                   % conditional responses (y), channel space
DCM.Rc = Ec;                   % conditional residuals (y), channel space
DCM.Hs = Hs;                   % conditional responses (y), source space
DCM.Ce = exp(-Eh);             % ReML error covariance
DCM.F  = F;                    % Laplace log evidence
DCM.ID = ID;                   % data ID
 
% and save
%--------------------------------------------------------------------------
DCM.options.Nmodes = Nm;
 
[fp fn fe] = fileparts(DCM.name); % ensure local saving
if isempty(fp); 
    namename = [fn fe];
else
    namename = [fp '/' fn fe];
end
DCM.name = namename

save(DCM.name, 'DCM', spm_get_defaults('mat.format'));
