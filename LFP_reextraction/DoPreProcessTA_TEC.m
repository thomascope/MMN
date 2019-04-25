function DoPreProcessTA(varargin)


%% files:

[F,S] = FilesTA;
if ~isempty(varargin)
    F = F(varargin{1});
    S = S(varargin{1});
end


%% jobs:

job = {...
    %'DoNotchTA'        ;... (remove mains noise)
    %'DoHiFilterTA'     ;... (apply 1 Hz high pass)
    %'DoLoFilterTA'     ;... (apply 48 Hz low pass)
    %'DoEpochRovTA'     ;... (epoch all 10 conditions)                           % %'FastICATA';...
    %'DoReject2TA'      ;... (trial rejection by AS recipe)
    %'DoCoregTA'        ;... (coregister individuals MRI & MEG; STILL YET TO MAKE: use template MRI if no MRI available)
    %'DoInvTA'          ;... (create IID inversion state)
    %'DoInvTA'          ;... (create MSP inversion state)
    %'DoAverage2TA'     ;... (robust averaging of the trials)
    %'DoExtractmTA'     ;... (extract the 6 sources of interest using IID)
    %'DoExtractmTA'     ;... (extract the 6 sources of interest using MSP)
    'CombPlanTA'       ;... (combine gradiometers)
    'DoInvTA'          ;... (create IID inversion state on combined gradiometers)
    'DoInvTA'          ;... (create MSP inversion state on combined gradiometers)
    'DoExtractmTA'     ;... (extract the 6 sources of interest from just the combine gradiometers using IID)
    'DoExtractmTA'     ;... (extract the 6 sources of interest from just the combine gradiometers using MSP)
    %'DoTFTA2'          ;... (time-frequency analysis on the mean trial using hilbert)
    %'DoTFTA2'          ;... (time-frequency analysis on all trials using hilbert)
    %'DoInvResultsTA'   ;... (ERP IID source-space stats, across subjects and conditions)
    %'DoInvResultsTA'   ;... (theta IID source-space meshes & images across subjects and conditions)
    %'DoInvResultsTA'   ;... (alpha IID source-space meshes & images across subjects and conditions)
    %'DoInvResultsTA'   ;... (beta1 IID source-space meshes & images across subjects and conditions)
    %'DoInvResultsTA'   ;... (beta2 IID source-space meshes & images across subjects and conditions)
    %'DoInvResultsTA'   ;... (beta IID source-space meshes & images across subjects and conditions)
    %'DoInvResultsTA'   ;... (gamma1 IID source-space meshes & images across subjects and conditions)
    %'DoInvResultsTA'   ;... (ERP MSP source-space stats, across subjects and conditions)
    %'DoInvResultsTA'   ;... (theta MSP source-space meshes & images across subjects and conditions)
    %'DoInvResultsTA'   ;... (alpha MSP source-space meshes & images across subjects and conditions)
    %'DoInvResultsTA'   ;... (beta1 MSP source-space meshes & images across subjects and conditions)
    %'DoInvResultsTA'   ;... (beta2 MSP source-space meshes & images across subjects and conditions)
    %'DoInvResultsTA'   ;... (beta MSP source-space meshes & images across subjects and conditions)
    %'DoInvResultsTA'   ;... (gamma1 MSP source-space meshes & images across subjects and conditions)
    %'dodcmTA2'         ;... (DCM ERP custom CMC on the 6 source LFPs created using IID for the full model (21))
    %'dodcmTA2'         ;... (DCM ERP custom CMC on the 6 source LFPs created using MSP for the full model (21))
    };


%% inputs:

% inputs for DoInvResultsTA:
InvVal = cell(length(F),1); [InvVal{:}] = deal('InvVal');
FB = cell(length(F),1); [FB{:}] = deal('FB');
FBt = cell(length(F),1); [FBt{:}] = deal([4 8]);
FBa = cell(length(F),1); [FBa{:}] = deal([8 12]);
FBb1 = cell(length(F),1); [FBb1{:}] = deal([13 19]);
FBb2 = cell(length(F),1); [FBb2{:}] = deal([20 29]);
FBb = cell(length(F),1); [FBb{:}] = deal([13 29]);
FBg1 = cell(length(F),1); [FBg1{:}] = deal([30 48]);

% inputs for dodcmTA1:
od = cd('/imaging/na01/LFPs');

nam1 = strcat('/imaging/na01/LFPs/',gfTA('inv1'));
nam1(cellfun(@(x) isempty(regexp(x,'.mat','once')),nam1)) = [];
dcm1 = cellfun(@(y) nout(2,@fileparts,y),nam1(cellfun(@(x) ~isempty(regexp(x,'dcm','once')),nam1)),'Uni',0);
%dcm1 = cellfun(@(x) x(1:end-11),dcm1,'Uni',0); % use to ignore ones you've already done
for k = 1:length(dcm1)
    nam1(cellfun(@(x) ~isempty(regexp(x,dcm1{k},'once')),nam1)) = [];
end

nam2 = strcat('/imaging/na01/LFPs/',gfTA('inv2'));
nam2(cellfun(@(x) isempty(regexp(x,'.mat','once')),nam2)) = [];
dcm2 = cellfun(@(y) nout(2,@fileparts,y),nam2(cellfun(@(x) ~isempty(regexp(x,'dcm','once')),nam2)),'Uni',0);
%dcm2 = cellfun(@(x) x(1:end-11),dcm2,'Uni',0); % use to ignore ones you've already done
for k = 1:length(dcm2)
    nam2(cellfun(@(x) ~isempty(regexp(x,dcm2{k},'once')),nam2)) = [];
end
cd(od)

% all collated inputs:
inputs = {...
    %{F};...
    %{strcat(cellfun(@(x) fileparts(x),F,'Uni',0),filesep,'f',cellfun(@(x) nout(2,@fileparts,x),F,'Uni',0),'.mat')};...
    %{strcat(cellfun(@(x) fileparts(x),F,'Uni',0),filesep,'ff',cellfun(@(x) nout(2,@fileparts,x),F,'Uni',0),'.mat')};...        % %{strcat(cellfun(@(x) fileparts(x),F,'Uni',0),filesep,'ff',cellfun(@(x) nout(2,@fileparts,x),F,'Uni',0),'.mat')};...
    %{strcat(cellfun(@(x) fileparts(x),F,'Uni',0),filesep,'fff',cellfun(@(x) nout(2,@fileparts,x),F,'Uni',0),'.mat')};...        % %{strcat(cellfun(@(x) fileparts(x),F,'Uni',0),filesep,'efff',cellfun(@(x) nout(2,@fileparts,x),F,'Uni',0),'.mat')};...
    %{strcat(cellfun(@(x) fileparts(x),F,'Uni',0),filesep,'efff',cellfun(@(x) nout(2,@fileparts,x),F,'Uni',0),'.mat')};...
    %{strcat(cellfun(@(x) fileparts(x),F,'Uni',0),filesep,'aefff',cellfun(@(x) nout(2,@fileparts,x),F,'Uni',0),'.mat'),  S};...
    %{strcat(cellfun(@(x) fileparts(x),F,'Uni',0),filesep,'aefff',cellfun(@(x) nout(2,@fileparts,x),F,'Uni',0),'.mat')};...
    %{strcat(cellfun(@(x) fileparts(x),F,'Uni',0),filesep,'aefff',cellfun(@(x) nout(2,@fileparts,x),F,'Uni',0),'.mat'), num2cell(2*ones(length(F),1))};...
    %{strcat(cellfun(@(x) fileparts(x),F,'Uni',0),filesep,'aefff',cellfun(@(x) nout(2,@fileparts,x),F,'Uni',0),'.mat')};...
    %{strcat(cellfun(@(x) fileparts(x),F,'Uni',0),filesep,'maefff',cellfun(@(x) nout(2,@fileparts,x),F,'Uni',0),'.mat')};...
    %{strcat(cellfun(@(x) fileparts(x),F,'Uni',0),filesep,'maefff',cellfun(@(x) nout(2,@fileparts,x),F,'Uni',0),'.mat'), num2cell(2*ones(length(F),1))};...        % %{strcat(cellfun(@(x) fileparts(x),F,'Uni',0),filesep,'aefff',cellfun(@(x) nout(2,@fileparts,x),F,'Uni',0),'.mat')};...        % %{strcat(cellfun(@(x) fileparts(x),F,'Uni',0),filesep,'maefff',cellfun(@(x) nout(2,@fileparts,x),F,'Uni',0),'.mat')};...
    {strcat(cellfun(@(x) fileparts(x),F,'Uni',0),filesep,'maefff',cellfun(@(x) nout(2,@fileparts,x),F,'Uni',0),'.mat')};...
    {strcat(cellfun(@(x) fileparts(x),F,'Uni',0),filesep,'Pmaefff',cellfun(@(x) nout(2,@fileparts,x),F,'Uni',0),'.mat')};...
    {strcat(cellfun(@(x) fileparts(x),F,'Uni',0),filesep,'Pmaefff',cellfun(@(x) nout(2,@fileparts,x),F,'Uni',0),'.mat'), num2cell(2*ones(length(F),1))};...
    {strcat(cellfun(@(x) fileparts(x),F,'Uni',0),filesep,'Pmaefff',cellfun(@(x) nout(2,@fileparts,x),F,'Uni',0),'.mat')};...
    {strcat(cellfun(@(x) fileparts(x),F,'Uni',0),filesep,'Pmaefff',cellfun(@(x) nout(2,@fileparts,x),F,'Uni',0),'.mat'), num2cell(2*ones(length(F),1))};...        % %{strcat(cellfun(@(x) fileparts(x),F,'Uni',0),filesep,'aefff',cellfun(@(x) nout(2,@fileparts,x),F,'Uni',0),'.mat')};...        % %{strcat(cellfun(@(x) fileparts(x),F,'Uni',0),filesep,'maefff',cellfun(@(x) nout(2,@fileparts,x),F,'Uni',0),'.mat')};...
    %{strcat(cellfun(@(x) fileparts(x),F,'Uni',0),filesep,'aefff',cellfun(@(x) nout(2,@fileparts,x),F,'Uni',0),'.mat')};...
    %{strcat(cellfun(@(x) fileparts(x),F,'Uni',0),filesep,'maefff',cellfun(@(x) nout(2,@fileparts,x),F,'Uni',0),'.mat')};...
    %{strcat(cellfun(@(x) fileparts(x),F,'Uni',0),filesep,'maefff',cellfun(@(x) nout(2,@fileparts,x),F,'Uni',0),'.mat')};...
    %{strcat(cellfun(@(x) fileparts(x),F,'Uni',0),filesep,'maefff',cellfun(@(x) nout(2,@fileparts,x),F,'Uni',0),'.mat'), FB, FBt};...
    %{strcat(cellfun(@(x) fileparts(x),F,'Uni',0),filesep,'maefff',cellfun(@(x) nout(2,@fileparts,x),F,'Uni',0),'.mat'), FB, FBa};...
    %{strcat(cellfun(@(x) fileparts(x),F,'Uni',0),filesep,'maefff',cellfun(@(x) nout(2,@fileparts,x),F,'Uni',0),'.mat'), FB, FBb1};...
    %{strcat(cellfun(@(x) fileparts(x),F,'Uni',0),filesep,'maefff',cellfun(@(x) nout(2,@fileparts,x),F,'Uni',0),'.mat'), FB, FBb2};...
    %{strcat(cellfun(@(x) fileparts(x),F,'Uni',0),filesep,'maefff',cellfun(@(x) nout(2,@fileparts,x),F,'Uni',0),'.mat'), FB, FBb};...
    %{strcat(cellfun(@(x) fileparts(x),F,'Uni',0),filesep,'maefff',cellfun(@(x) nout(2,@fileparts,x),F,'Uni',0),'.mat'), FB, FBg1};...
    %{strcat(cellfun(@(x) fileparts(x),F,'Uni',0),filesep,'maefff',cellfun(@(x) nout(2,@fileparts,x),F,'Uni',0),'.mat'), InvVal, num2cell(2*ones(length(F),1))};...
    %{strcat(cellfun(@(x) fileparts(x),F,'Uni',0),filesep,'maefff',cellfun(@(x) nout(2,@fileparts,x),F,'Uni',0),'.mat'), FB, FBt, InvVal, num2cell(2*ones(length(F),1))};...
    %{strcat(cellfun(@(x) fileparts(x),F,'Uni',0),filesep,'maefff',cellfun(@(x) nout(2,@fileparts,x),F,'Uni',0),'.mat'), FB, FBa, InvVal, num2cell(2*ones(length(F),1))};...
    %{strcat(cellfun(@(x) fileparts(x),F,'Uni',0),filesep,'maefff',cellfun(@(x) nout(2,@fileparts,x),F,'Uni',0),'.mat'), FB, FBb1, InvVal, num2cell(2*ones(length(F),1))};...
    %{strcat(cellfun(@(x) fileparts(x),F,'Uni',0),filesep,'maefff',cellfun(@(x) nout(2,@fileparts,x),F,'Uni',0),'.mat'), FB, FBb2, InvVal, num2cell(2*ones(length(F),1))};...
    %{strcat(cellfun(@(x) fileparts(x),F,'Uni',0),filesep,'maefff',cellfun(@(x) nout(2,@fileparts,x),F,'Uni',0),'.mat'), FB, FBb, InvVal, num2cell(2*ones(length(F),1))};...
    %{strcat(cellfun(@(x) fileparts(x),F,'Uni',0),filesep,'maefff',cellfun(@(x) nout(2,@fileparts,x),F,'Uni',0),'.mat'), FB, FBg1, InvVal, num2cell(2*ones(length(F),1))};...
    %{nam1, num2cell(21*ones(length(nam1),1))};...
    %{nam2, num2cell(21*ones(length(nam2),1))};...
    };


%% the 'doing' bit:

if strcmp('8.5.0.197613 (R2015a)',version)
    fprintf('\n\nto run:')
    disp(table(job))
    fprintf('\n\non files:')
    table([cellfun(@(x) nout(2,@fileparts,x),F,'Uni',0),S]')
    fprintf('\n\n')
end

% open pool:
try matlabpool(cbupool(length(F)));
catch err
    if strcmp('parallel:lang:matlabpool:OpenDeprecation',err.identifier)
        P = parpool(cbupool(length(F)));
    end
end

% run:
try sr = strrep(strrep(datestr(datetime),':',''),' ','_'); %alphanumTA(5);
catch err, sr = date;
end
mkdir(['/imaging/na01/N_' sr])
for k = 1:length(job)
    j = job{k};
    i = inputs{k};
    if strcmp('dodcmTA1',j)
        N = length(nam1);
    else N = length(F);
    end
    parfor k2 = 1:N
        i2 = cellfun(@(x) x{k2},i,'Uni',0);
        try
            feval(j,i2{:});
            disp(['job ' num2str(k) ' ; meg ' num2str(k2)])
            n = 1;
            save_parTA(['/imaging/na01/N_' sr '/N' num2str(k) '_' num2str(k2) '.mat'],'n')
        catch err
        end
    end
end

% close pool:
if exist('P','var'), delete(P), end % parpool close force CBU_Cluster % matlabpool close force CBU_Cluster

% delete unneeded files:
% for k = 1:length(F)
%     od = cd(fileparts(F{k}));
%     [~,f1,f2] = fileparts(F{k});
%     delete(['ff' f1 f2],['f' f1 f2]) % keep montaged-averaged & pre-averaged; delete all the others % ['iefff*' f1 f2],['efff' f1 f2],['fff' f1 f2],
%     delete(['ff' f1 '.dat'],['f' f1 '.dat']) % ['iefff*' f1 '.dat'],['efff' f1 '.dat'],['fff' f1 '.dat'],
%     cd(od)
% end


end

