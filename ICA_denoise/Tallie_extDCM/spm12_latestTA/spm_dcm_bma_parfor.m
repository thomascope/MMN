function bma = spm_dcm_bma_parfor(post,post_indx,subj,Nsamp,oddsr)
% Model-independent samples from DCM posterior
% FORMAT BMA = spm_dcm_bma(DCM)
% FORMAT bma = spm_dcm_bma(post,post_indx,subj,Nsamp,oddsr)
%
% DCM   - {subjects x models} cell array of DCMs over which to average
% ---------------------------------------------------------------------
%     DCM{i,j}.Ep - posterior expectation
%     DCM{i,j}.Cp - posterior covariances
%     DCM{i,j}.F  - free energy
%
% BMA   - Baysian model average structure
% ---------------------------------------------------------------------
%     BMA.Ep      - BMA posterior mean
%     BMA.Cp      - BMA posterior VARIANCE
%     BMA.F       - Accumulated free energy over subjects;
%     BMA.P       - Posterior model probability over subjects;
%
%     BMA.SUB.Ep  - subject specific BMA posterior mean
%     BMA.SUB.Sp  - subject specific BMA posterior variance
%     BMA.nsamp   - Number of samples
%     BMA.Nocc    - number of models in Occam's window
%     BMA.Mocc    - index of models in Occam's window
%
% If DCM is an array, Bayesian model averaging will be applied over 
% subjects (i.e., over columns) using FFX Baysian parameter averaging
%
%--------------------------------------------------------------------------
% OR
%--------------------------------------------------------------------------
%
% post      [Ni x M] vector of posterior model probabilities
%           If Ni > 1 then inference is based on subject-specific RFX posterior
% post_indx models to use in BMA (position of models in subj structure)
% subj      subj(n).sess(s).model(m).fname: DCM filename
% Nsamp     Number of samples (default = 1e3)
% oddsr     posterior odds ratio for defining Occam's window (default=0, ie
%           all models used in average)
%
% bma       Returned data structure contains
%
%           .nsamp  Number of samples
%           .oddsr  odds ratio
%           .Nocc   number of models in Occam's window
%           .Mocc   index of models in Occam's window
%           .indx   subject specific indices of models in Occam's window
%
%           For `Subject Parameter Averaging (SPA)':
%
%           .mEp    posterior mean 
%           .sEp    posterior SD           
%           .mEps   subject specific posterior mean 
%           .sEps   subject specific posterior SD
%
%           use the above values in t-tests, ANOVAs to look for significant
%           effects in the group
%
%           For `Group Parameter Averaging (GPA)':
%
%           The following structures contain samples of the DCM A,B,C and D
%           matrices from the group posterior density. See pages 6 and 7 of [1]
%
%           .a [dima x Nsamp] 
%           .b [dima x Nsamp] 
%           .c [dima x Nsamp] 
%           .d [dima x Nsamp] 
%                       
%           Use these to make inferences using the group posterior density approach. 
%           Essentially, for each parameter, GPA gets a sample which is the average 
%           over subjects. The collection of samples then constitutes a distribution of
%           the group mean from which inferences can be made directly. This is to
%           be contrasted with SPA where, for each subject, we average over
%           samples to get a mean for that subject. Group level inferences
%           are then made using classical inference. SPA is the standard
%           approach.
%
%
%           For RFX BMA, different subject can have different models in
%           Occam's window (and different numbers of models in Occam's
%           window)
%
% This routine implements Bayesian averaging over models and subjects
%
% See [1] W Penny, K Stephan, J. Daunizeau, M. Rosa, K. Friston, T. Schofield 
% and A Leff. Comparing Families of Dynamic Causal Models. 
% PLoS Computational Biology, Mar 2010, 6(3), e1000709.
%__________________________________________________________________________
% Copyright (C) 2009 Wellcome Trust Centre for Neuroimaging

% Will Penny
% $Id: spm_dcm_bma.m 7081 2017-05-27 19:36:09Z karl $

% defaults
%--------------------------------------------------------------------------
if nargin < 4 || isempty(Nsamp)
    Nsamp = 1e3; % TA removed
    %Nsamp = length(post)*100; % TA replaced above with this
end
if nargin < 5 || isempty(oddsr)
    oddsr = 0;
end

% inputs are DCMs � assemble input arguments
%--------------------------------------------------------------------------
% try
%     pp = cbupool(80);
%     pp.ResourceTemplate = '-l nodes=^N^,mem=16GB,walltime=96:00:00';
%     parpool(pp,80);
% end

if nargin == 1
    
    if ~iscell(post), post = {post}; end
    
    DCM   = post;
    [n,m] = size(DCM); disp('PARFOR 3')
    for i = 1:n
        parfor j = 1:m
            if ~isfield(DCM{i,j}, 'Ep')
                error(['Could not average: subject %d model %d ' ...
                       'not estimated'], i, j);
            end            
            subj{i,1,j} = DCM{i,j}.Ep; % subj(i).sess(1).model(j).Ep
            subjCp{i,1,j} = DCM{i,j}.Cp;
            F(i,j) = DCM{i,j}.F;
        end
        disp('PARFOR 3 END')
    end
    
    % (FFX) posterior over models
    %----------------------------------------------------------------------
    F    = sum(F,1);
    F    = F - max(F);
    P    = exp(F);
    post = P/sum(P);
    indx = 1:m;
    
    % BMA (and BPA)
    %----------------------------------------------------------------------
    bma        = spm_dcm_bma_parfor(post,indx,{subj subjCp},Nsamp);
    BMA.Ep     = bma.mEp;
    BMA.Cp     = spm_unvec(spm_vec(bma.sEp).^2,bma.sEp);
    
    BMA.nsamp  = bma.nsamp;
    BMA.Nocc   = bma.Nocc;
    BMA.Mocc   = bma.Mocc;
    BMA.F      = F;
    BMA.P      = P;
    
    for i = 1:n
        BMA.SUB(i).Ep = bma.mEps{i};
        BMA.SUB(i).Cp = spm_unvec(spm_vec(bma.sEps{i}).^2,bma.sEps{i});
    end
    bma        = BMA;
    return
end

subjCp = subj{2};
subj = subj{1};
Nsub = size(subj,1);
Nses = size(subj,2);

% Number of regions
%--------------------------------------------------------------------------
% TA removed try catch as not using
dcm_fmri = 0;

firstsub  = 1;
firstmod  = 1;

Ep  = [];
[Ni,M] = size(post);

rfx = 1;

if rfx
    
    for i = 1:Ni

        mp          = max(post(i,:));
        post_ind{i} = find(post(i,:)>mp*oddsr);
        Nocc(i)     = length(post_ind{i});
        disp(' ');
        disp(sprintf('Subject %d has %d models in Occams window',i,Nocc(i)));

        if Nocc(i) == 0
            return;
        end

%         for occ = 1:Nocc(i),
%             m = post_ind{i}(occ);
%             disp(sprintf('Model %d, <p(m|Y>=%1.2f',m,post(i,m)));
%         end

        % Renormalise post prob to Occam group
        %------------------------------------------------------------------
        renorm{i} = post(i,post_ind{i});
        sp             = sum(renorm{i},2);
        renorm{i} = renorm{i}./(sp*ones(1,Nocc(i)));

        % Load DCM posteriors for models in Occam's window
        %------------------------------------------------------------------
        params = repmat({struct},Ni,Nocc);
        hwb = waitbar(0,[num2str(i) ' ...']);
        for kk = 1:Nocc(i)

            sel     = post_indx(post_ind{i}(kk));

            params{i,kk}.Ep  = subj{i,1,sel};
            params{i,kk}.vEp = spm_vec(params{i,kk}.Ep);
            params{i,kk}.Cp  = full(subjCp{i,1,sel});

%             if dcm_fmri
%                 dimDtmp = size(params{i,kk}.Ep.D,3);
%                 if dimDtmp ~= 0, dimD = dimDtmp; firstsub = i; firstmod = kk;end
%             end
            
            % Average sessions
            %--------------------------------------------------------------
            % TA REMOVED as unused
            %--------------------------------------------------------------

            [evec, eval] = eig(params{i,kk}.Cp);
            deig         = diag(eval);

            params{i,kk}.dCp = deig;
            params{i,kk}.vCp = evec;
            waitbar(kk/Nocc(i))
        end
        close(hwb)
    end
    
else % Use an FFX
    % TA removed as not using
end

% Pre-allocate sample arrays
%--------------------------------------------------------------------------
Np = length(params{firstsub,firstmod}.vEp);

% get dimensions of a b c d parameters
%--------------------------------------------------------------------------
if dcm_fmri

    Nr      = nreg*nreg;
    nmods   = size(DCM.Ep.B,3);
    
    Etmp.A  = zeros(nreg,nreg,Nsamp);
    Etmp.B  = zeros(nreg,nreg,nmods,Nsamp);
    Etmp.C  = zeros(nreg,min,Nsamp);
    Etmp.D  = zeros(nreg,nreg,dimD,Nsamp);

    dima    = Nr;
    dimb    = Nr+Nr*nmods;
    dimc    = Nr+Nr*nmods+nreg*min;

end

clear Ep
disp('')
disp('Averaging models in Occams window...')

Ep_all = zeros(Np,Nsub);
Ep_sbj = repmat({zeros(Np,Nsub,Nsamp)},Nsamp,1);
Ep     = zeros(Np,Nsamp);
m = zeros(Nsamp,1);
dsig = cell(Nsamp,1);
vsig = cell(Nsamp,1);
Ep_all = cell(Nsamp,1);
Ep_subj = cell(Nsamp,1);
hwb = waitbar(0,'distr build');
for i=1:Nsamp % iterates for a reasonable number of samples to create a decent distribution, where models with higher evidences are included more often and therefore have higher weight
    
    % Pick a model
    %----------------------------------------------------------------------
%     if ~rfx % TA removed as not using
%         m(i) = spm_multrnd(post,1);
%     end
    % Pick parameters from model for each subject
    %----------------------------------------------------------------------
    for n=1:Nsub

        %clear mu dsig vsig
        
        %if rfx
            m(i) = spm_multrnd(renorm{n},1);
        %end

        mu                = params{n,m(i)}.vEp;
        nmu               = length(mu);
        dsig{i}              = params{n,m(i)}.dCp(1:nmu,1);
        vsig{i}(:,:)         = params{n,m(i)}.vCp(1:nmu,1:nmu);

        tmp               = spm_normrnd(mu,{dsig{i},vsig{i}},1);
        
        Ep_all{i}(1:nmu,n)   = tmp(:);
        Ep_sbj{i}(1:nmu,n) = Ep_all{i}(1:nmu,n); 

    end
    
    % Average over subjects
    %----------------------------------------------------------------------
    Ep(:,i) = mean(Ep_all{i},2);
    waitbar(i/Nsamp)
end
close(hwb)
%disp('done parfor 3/3')

% save mean parameters
%--------------------------------------------------------------------------
Ep_avg     = mean(Ep,2);
Ep_std     = std(Ep,0,2);
Ep_avg     = spm_unvec(Ep_avg,params{1,1}.Ep);
Ep_std     = spm_unvec(Ep_std,params{1,1}.Ep);
bma.mEp    = Ep_avg;
bma.sEp    = Ep_std;

Ep_avgsbj  = mean(cat(3,Ep_sbj{:}),3);
Ep_stdsbj  = std(cat(3,Ep_sbj{:}),0,3);

for is=1:Nsub
    bma.mEps{is}=spm_unvec(Ep_avgsbj(:,is),params{1,1}.Ep);
    bma.sEps{is}=spm_unvec(Ep_stdsbj(:,is),params{1,1}.Ep);
end

if dcm_fmri
    bma.a  = spm_unvec(Ep(1:dima,:),Etmp.A);
    bma.b  = spm_unvec(Ep(dima+1:dimb,:),Etmp.B);
    bma.c  = spm_unvec(Ep(dimb+1:dimc,:),Etmp.C);
    if dimD ~=0
        bma.d  = spm_unvec(Ep(dimc+1:dimc+Nr*dimD,:),Etmp.D);
    else
        bma.d  = Etmp.D;
    end
end

% storing parameters
% -------------------------------------------------------------------------
bma.nsamp = Nsamp;
bma.oddsr = oddsr;
bma.Nocc  = Nocc;
bma.Mocc  = post_ind;
