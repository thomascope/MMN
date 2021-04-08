function all_outfiles = maxfilter_this_participant(rawdatapath,subjfolder,outfilename,repeatmaxfilter)
%rawdatapath is a cell array of strings, each containing a raw data path


%% Add paths

%addpath(genpath('/imaging/rowe/archive/users/hp02/spm12b'));
addpath(genpath('/imaging/local/software/mne'));
addpath('/imaging/rowe/archive/users/hp02/finger_tapping08/analysis_spm/new_functions');
addpath('/imaging/rowe/archive/users/hp02/pnfa_mmn/');

addpath(genpath('/imaging/local/meg_misc'))
addpath(genpath('/neuro/meg_pd_1.2/'))
addpath(genpath('/imaging/rowe/archive/users/hp02/mmn_08/analysis_spm/new_spm_functions'))

% Define destination
datapath = [subjfolder outfilename];

if ~exist(datapath,'dir')
    mkdir(datapath)
end

%% Now run Maxfilter with Rik's parameters
mvcomp_fail = ones(1,max(size(rawdatapath)),1,1);  % group, subject, experiment, run, Turn off all mvcomp, since seems to fail randomly!

movfile = 'trans_move.txt'; % This file will record translations between runs

basestr = ' -ctc /neuro/databases/ctc/ct_sparse.fif -cal /neuro/databases/sss/sss_cal.dat';
basestr = [basestr ' -linefreq 50 -hpisubt amp'];
basestr = [basestr ' -force'];
%maxfstr = '!/neuro/bin/util/x86_64-pc-linux-gnu/maxfilter-2.2 '
maxfstr = '!/neuro/bin/util/maxfilter-2.2';

for this_file = 1:max(size(rawdatapath))
    raw_file = rawdatapath{this_file};
    [~,raw_stem,~] = fileparts(raw_file);
    
    if ~exist(fullfile(datapath,sprintf('%s_trans1st.fif',raw_stem)),'file') || repeatmaxfilter
        
        
        % Fit sphere (since better than MaxFilter does)
        incEEG = 0;
        if exist(fullfile(datapath,'fittmp.txt')); delete(fullfile(datapath,'fittmp.txt')); end
        if exist(fullfile(datapath,sprintf('run_%02d_hpi.txt',this_file))); delete(fullfile(datapath,sprintf('run_%02d_hpi.txt',this_file)));  end
        [orig,rad,fit] = meg_fit_sphere(raw_file,datapath,sprintf('%s_hpi.txt',raw_stem),incEEG);
        delete(fullfile(datapath,'fittmp.txt'));
        
        origstr = sprintf(' -origin %d %d %d -frame head',orig(1),orig(2),orig(3))
        badstr  = sprintf(' -autobad %d -badlimit %d',900,7); % 900s is 15mins - ie enough for whole recording!
        
        outfile = fullfile(datapath,sprintf('%s_trans1st',raw_stem))
        
        % 1. Bad channel detection (this email says important if doing tSSS later https://www.jiscmail.ac.uk/cgi-bin/webadmin?A2=NEUROMEG;d3f363f3.1205)
        
        outfile = fullfile(datapath,sprintf('%s_bad',raw_stem));
        filestr = sprintf(' -f %s -o %s.fif',raw_file,outfile);
        
        % Write out movements too...
        posfile = fullfile(datapath,sprintf('%s_headpos.txt',raw_stem));
        compstr = sprintf(' -headpos -hpistep 10 -hp %s',posfile);
        
        finstr = [maxfstr filestr origstr basestr badstr compstr sprintf(' -v | tee %s.log',outfile)]
        
        rik_eval(finstr);
        
        
        delete(sprintf('%s.fif',outfile));
        
        % Pull out bad channels from logfile:
        badfile = sprintf('%s.txt',outfile); delete(badfile);
        rik_eval(sprintf('!cat %s.log | sed -n -e ''/Detected/p'' -e ''/Static/p'' | cut -f 5- -d '' '' > %s',outfile,badfile));
        %if s ==1
        %    tmp=dlmread(badfile); Nbuf = size(tmp,1);
        %else
        tmp=dlmread(badfile,' '); Nbuf = size(tmp,1);
        %end
        tmp=reshape(tmp,1,prod(size(tmp)));
        tmp=tmp(tmp>0); % Omit zeros (padded by dlmread):
        
        % Get frequencies (number of buffers in which chan was bad):
        [frq,allbad] = hist(tmp,unique(tmp));
        
        % Mark bad based on threshold (currently ~5% of buffers (assuming 500 buffers)):
        badchans = allbad(frq>0.05*Nbuf);
        if isempty(badchans)
            badstr = '';
        else
            badstr = sprintf(' -bad %s',num2str(badchans))
        end
        
        
        %% 2. tSSS and trans to first file (ie, align within subject if multiple runs)
        
        tSSSstr = ' -st 10 -corr 0.98'; %'tSSSstr = '';
        
        if mvcomp_fail(1,this_file,1,1) == 1
            compstr = '';
        else
            compstr = sprintf(' -movecomp inter');
        end
        
        outfile = fullfile(datapath,sprintf('%s_trans1st',raw_stem))
        transfstfile = [outfile '.fif'];
        transtr = '';
        
        dsstr = ' -ds 4';   % downsample to 250Hz
        %dsstr = ' -ds 1';   % don't downsample
        
        filestr = sprintf(' -f %s -o %s.fif',raw_file,outfile);
        finstr = [maxfstr filestr basestr badstr tSSSstr compstr origstr transtr dsstr sprintf(' -v | tee %s.log',outfile)]
        rik_eval(finstr);
        
        rik_eval(sprintf('!echo ''Trans 1st...'' >> %s',movfile));
        rik_eval(sprintf('!cat %s.log | sed -n ''/Position change/p'' | cut -f 7- -d '' '' >> %s',outfile,movfile));
        
        
        %% 3. trans to default helmet space (align across subjects)
        
        transfile = fullfile(datapath,sprintf('%s_trans1stdef',raw_stem))
        transtr = sprintf(' -trans default -origin %d %d %d -frame head -force',orig+[0 -13 6])
        filestr = sprintf(' -f %s.fif -o %s.fif',outfile,transfile);
        finstr = [maxfstr filestr transtr sprintf(' -v | tee %s.log',outfile)]
        rik_eval(finstr);
        
        rik_eval(sprintf('!echo ''Trans def...'' >> %s',movfile));
        rik_eval(sprintf('!cat %s.log | sed -n ''/Position change/p'' | cut -f 7- -d '' '' >> %s',outfile,movfile));
        disp(['Maxfiltering complete for file ' num2str(this_file)])
        
    else
        disp(['Maxfiltering skipped for file ' num2str(this_file)])
        cd(datapath)
    end
end


delete tempdata_*

%%        Convert and merge if necessary:
if max(size(rawdatapath))==1
    raw_file = rawdatapath{this_file};
    [~,raw_stem,~] = fileparts(raw_file);
    S = [];
    S.dataset = fullfile(datapath,sprintf('%s_trans1st.fif',raw_stem));
    %S.channels = 'all';
    S.outfile  = outfilename;
    all_outfiles = S.outfile;
    if ~exist(sprintf('f%s.mat',subjfolder), 'file')
        D = spm_eeg_convert(S);
    else
        D = spm_eeg_load(sprintf('%s.mat',subjfolder));
    end
else
    premerged_files = {};
    for this_file = 1:max(size(rawdatapath))
        raw_file = rawdatapath{this_file};
        [~,raw_stem,~] = fileparts(raw_file);
        S = [];
        S.dataset = fullfile(datapath,sprintf('%s_trans1st.fif',raw_stem));
        %S.channels = 'all';
        S.outfile  = [outfilename '_' num2str(this_file)];
        all_outfiles{this_file} = S.outfile;
        if ~exist(sprintf('f%s.mat',subjfolder), 'file')
            D = spm_eeg_convert(S);
        else
            D = spm_eeg_load(sprintf('%s.mat',subjfolder));
        end
    end
end

disp(['Files converted'])

