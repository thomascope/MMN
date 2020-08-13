function bashrunTA(n)


% load files:
RouteDir = '/imaging/na01/';
if ~exist([RouteDir 'DCMTA'],'dir')
    mkdir([RouteDir 'DCMTA'])
end

cd([RouteDir 'DCM_LFPs'])
F = gfTA('.mat');
cd(RouteDir)

% run DCM for specified file:
disp('running DCMTA')
DCM = DCMTA(F{n});
r = fileread('DCMTA.m');
DCM.options.MFILE = r;

% save DCM:
disp('saving DCM structure...')
save([RouteDir 'DCMTA/' nout(2,@fileparts,F{n}) '_dcm.mat'],'DCM')


end

