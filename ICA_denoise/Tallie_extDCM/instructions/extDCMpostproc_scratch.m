
%
% extDCMpostproc_notes
%

%% PEB 2 - the second level design matrix

% load all the reduced DCMs created in the above step into a single
% structure:
GCM = {};
X1 = [];
for k = 1:length(subdirname_DCM)
    GCMtmp = load([dirname_DCM subdirname_DCM{k} nam]);
    GCMtmp = GCMtmp.DCM;
    GCM(end+1:end+length(GCMtmp),1) = GCMtmp;
    X1(end+1:end+length(GCMtmp),end+1) = ones(length(GCMtmp),1);
end

% Make a model structure, M, with a design matrix specified in .X. The
% matrix needs a column of ones first, which represents the group mean. In
% this example the rest of the columns represent the condition differences,
% (from the differen DCM_all folders) with an example interaction in the
% final column:
X = [X1(:,1)-X1(:,2) X1(:,3)-X1(:,4)]; % condition differences
X(:,4) = X(:,2).*X(:,3); % interaction
M = GCM{1}.M;
M.X = [ones(size(GCM,1),1) X];
% N.B. the above is just an illustrative example. In my own data,
% conditions had to be further separated within these groups. Be aware you
% will need to build your own design matrix, X, for your individual needs.

% Now you have your DCMs and a design matrix, you need to decide which
% parameters you a interested in. This is basically a list of connections, 
% either extrinsic or intrinsic, that you are interested in. e.g. It could 
% just be all the extrinsics: ... 
params = {'A'};
% ... or it could be a specific list of intrinsics - in the below example 
% all GABA-A synapses:
params = {'H(1,1,1,3)'
          'H(2,2,1,3)'
          'H(3,3,1,3)'
          'H(4,4,1,3)'
          'H(5,5,1,3)'
          'H(6,6,1,3)'
          'H(2,3,1,3)'
          'H(1,3,1,3)'
          'H(4,5,1,3)'
          'H(6,5,1,3)'};

% Run the second level PEB and save the results. The .Ep field in the PEB
% has a column for each column of the design matrix:
[PEB_2,DCM_2] = spm_dcm_peb(GCM,M,params);
save('PEB_2','PEB_2','DCM_2')

% The above results can then be used to find the posterior p-vals for the
% parameters in this PEB (found in .pP). The relevant fields need reshaping
% to the size of the design matrix:
BMA_2 = spm_dcm_peb_bmc(PEB_2);
BMA_2.Ep = reshape(BMA_2.Ep,[],size(M.X,2));
BMA_2.Cp = reshape(BMA_2.Cp,[],size(M.X,2));
BMA_2.Pp = reshape(BMA_2.Pp,[],size(M.X,2));
save('PEB_2','BMA_2','-append')

% You are now ready to view some results. I have created a plotting script 
% based on my 6-node, 6-population model. It shows the design matrix and 
% then a row of graph plots, one plot for each column of the design matrix 
% (not including the first column, as this is how the PEB has adjusted the 
% parameters for the group mean - not the constrasts you are interested in).
% To use it for a single BMA to view both the parameter estimates and their
% associated P-value:
clear B
B.BMA.bma = BMA;
plotPEBBMA(B,'Pp','cax',.5)

% To use it for multiple BMAs in one figure (only allowing you to look at 
% the parameter estimates, not their P-values):
clear B
B.BMA1.bma = BMA1;
B.BMA2.bma = BMA2;
plotPEBBMA(B,'Ep','cax',.5)


