function O = MaxFiltSPM_batchinfo_XLSXinput(xlsx,rawpath,mfpath)
%
% This function takes your spreadsheet of MEG raw-file locations and
% desired maxfiltered-file locations, along with the bad channels for each
% file, and preparing it for use in MaxFiltSPM_x.m.
%
% INPUTS:
% xlsx    = full path and filename of an xlsx with 4 columns:
%           col1 = paths of files to maxfilter inside 'rawpath'
%           col2 = paths of file locations after maxfilter inside 'mfpath'
%           col3 = bad MEG channels (comma-separately)
%           col4 = bad EEG channels (comma-separately)
% rawpath = outer path to the raw datafiles (e.g. '/megdata/cbu/ftdrug/')
% mfpath  = outer path to the maxfiltered datafiles (e.g. '/imaging/group/rowelab/MEM/resting_closed/maxfilt/')
%

% import data:
d = importdata(xlsx);
df = [cellfun(@(x) strsplit(strrep(fileparts(x),rawpath,''),'/'),d.textdata(:,1),'Uni',0) ...
      cellfun(@(x) strsplit(strrep(fileparts(x),mfpath,''),'/'),d.textdata(:,2),'Uni',0)];
dn = cellfun(@(x) nout(2,@fileparts,x),d.textdata(:,1:2),'Uni',0);

spm('defaults','EEG');

% organize meta structures:
subjects   = cellfun(@(x) x{1},df(:,1),'Uni',0); % {'' ''}; % subject folder on the CBUs /megdata
dates      = cellfun(@(x) x{2},df(:,1),'Uni',0); % {'' ''}; % date folder inside the subject folder
subj_new   = cellfun(@(x) x{1},df(:,2),'Uni',0); % TA added this. Much better!
blocksin   = cellfun(@(x) {strrep(x,'_raw','')},dn(:,1),'Uni',0); % {{''} {''}}; % orig file name (without the raw suffix)
blocksout  = cellfun(@(x) {x},dn(:,2),'Uni',0); % {{''} {''}}; % what you want it to be called

badchannels = cellfun(@(y) cellfun(@(x) ['MEG' x],y,'Uni',0),cellfun(@(x) strsplit(x,','),d.textdata(:,3),'Uni',0),'Uni',0);
badeeg = cellfun(@(x) strsplit(x,','),d.textdata(:,4),'Uni',0);

for i = 1:length(subjects) % cnt
    if strcmp(class(subjects{i}),'char')
        subjects1{i}.newname = subj_new{i};
        subjects1{i}.oldname = subjects{i};
    else subjects1{i}.newname = subjects{i}.newname;
        subjects1{i}.oldname = subj_new{i}.oldname;
    end
end
subjects = subjects1;

% create output structure for maxfilter script
O.subjects = subjects;
O.dates = dates;
O.blocksin = blocksin;
O.blocksout = blocksout;
O.badchannels = badchannels;
O.badeeg = badeeg;

end

