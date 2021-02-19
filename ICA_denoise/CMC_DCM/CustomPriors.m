function DCM = CustomPriors(DCM,Ns)

%% Custom options for CMCs
        
        % This will put a bunch of custom (non-DCM/SPM) functions or priors into a
        % new structure in the DCM. These will be copied over in the call to
        % spm_dcm_csd. As such this will require a custom version of spm_dcm_csd
        % which includes a line to detect whether this struct has been provided and
        % to copy over it's contents prior to inversion.
        
        DCM.CUSTOM    = [];
        DCM.CUSTOM.f  = 'spm_fx_cmcTA';
        
        DCM.CUSTOM.pE   = [];
        DCM.CUSTOM.pE.G = zeros(Ns,13);                 ... intrinsic connectivity
        DCM.CUSTOM.pC.G = zeros(Ns,13); % off
        
        Self = find(diag(speye(Ns).*( DCM.A{1}+DCM.A{2} ))); % SP gain only if in model
        DCM.CUSTOM.pC.G(Self,7) = 1/8;
        
        
        DCM.CUSTOM.pE.T = zeros(Ns,4);                  ... population time const
        DCM.CUSTOM.pC.T = zeros(Ns,4)+1/8;
        
        % Below are the priors for the forward projection - I have hacked
        % spm_lx_erp.m to permit 3 deifferent sets of contirbuting-states: 1 for
        % frontal sources, 1 for parietal and 1 for auditory, so that we can assess
        % whether reduction in the contibution of frontal layers 2/3 to signal
        % without it being blurred out by averaging over all nodes.
        
        DCM.CUSTOM.gE.J = sparse([1 3 7],1,[.2 .8 .2],8,1)'; ... spatial model
        DCM.CUSTOM.gC.J = sparse([1 3 7],1,1/8       ,8,1)';
        
        DCM.CUSTOM.gC.J = repmat(DCM.CUSTOM.gC.J,[4 1]);
        %DCM.CUSTOM.gC.J = repmat(DCM.CUSTOM.gC.J,[3 1]);
        
        DCM.CUSTOM.gE.J = repmat(DCM.CUSTOM.gE.J,[4 1]);
        
end