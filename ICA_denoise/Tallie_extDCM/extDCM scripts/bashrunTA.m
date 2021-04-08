function bashrunTA(n,RouteDir)


%% setup (edit for your imaging space requirements):

% create directory for saving the results:
RouteDir = char(RouteDir);
if ~exist(RouteDir,'dir')
    mkdir(RouteDir)
end
if ~exist([RouteDir 'tDCMTA'],'dir')
    mkdir([RouteDir 'tDCMTA'])
end


%% load the input files (edit for your needs):

% % This might be the contents of a folder you specify:
% cd('/imaging/your/path/to/extracted/LFPs')
% F = gfTA;

% % Or it might be a list you've already complied in an mfile:
% F = FileList;

% If you've used my preproccessed data it might be:
meg = FilesTA;
meg_file = cellfun(@(x) [nout(2,@fileparts,x) nout(3,@fileparts,x)],meg,'Uni',0);
F = cellfun(@(x,y) ['/imaging/rowe/users/na01/LFPs/aeF/LFP6inv2_aeF' x],meg_file,meg,'Uni',0);


%% run DCM for specified file:

% run it:
cd(RouteDir)
disp('running DCMTA')
DCM = DCMTA(F{n});

% append the exact script used to the outputted data structure:
r = fileread('DCMTA.m');
DCM.options.MFILE = r;

% save it:
disp('saving DCM structure...')
save([RouteDir 'tDCMTA/' nout(2,@fileparts,F{n}) '_dcm.mat'],'DCM')


end

