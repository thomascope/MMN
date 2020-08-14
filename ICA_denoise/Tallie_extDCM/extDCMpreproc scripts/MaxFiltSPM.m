%Maxfilter final edit by TEC, after edited by ES, based on script by Jason Taylor

% In parfor loop run like this
% parfor i = 1:n
% es_batch_/neuro/bin/util/maxfilter_saber_followup(i)
% end

function MaxFiltSPM(do_subjects)

MaxFiltSPM_batchinfo

% es edit- Assemble info needed to locate files for each subject
for s=1:length(subjects)
    subject{s} = { subjects{s}.newname dates{s} };
end

nr_sbj = length(subject);

try do_subjects;    % if do_subjects not defined, do all subjects
catch
    do_subjects = [1:nr_sbj];
end;

% Check file names and paths
checkflag = 0;
for ss = do_subjects,
    nr_bls = length( blocksin{ss} );
    if length(blocksin{ss}) ~= length(blocksout{ss}),
        checkflag = 1;
        fprintf(1, 'Different number of input and output names for subject %d (%s, %s)\n', ss, subject{ss}{1}, subject{ss}{2});
    end;
    for bb = 1:nr_bls,
        try
            rawpath = fullfile( rawpathstem, subjects{ss}.oldname, subject{ss}{2} );
        catch
            error(rawpath)
            %rawpath = fullfile( rawpathstem, subject{ss}{1}, subject{ss}{2} );
        end
        rawfname = fullfile( rawpath, [blocksin{ss}{bb} '_raw.fif'] );
        outpath = fullfile( pathstem, subject{ss}{1} );
        if ~exist( outpath, 'dir' ),
            success = mkdir( outpath );
            if ~success,
                checkflag = 1;
                fprintf(1, 'Could not create directory %s\n', outpath);
            end;
        end;
        if ~exist( rawfname, 'file' ),
            checkflag = 1;
            fprintf(1, '%s does not exist\n', rawfname);
        end;
    end;
end;
if checkflag,
    fprintf(1, 'You''ve got some explaining to do.\n');
    return;
end;

try
    for ss = do_subjects,
        nr_bls = length( blocksin{ss} );

        for bb = 1:nr_bls,

            try
                rawpath = fullfile( rawpathstem, subjects{ss}.oldname, subject{ss}{2} );
            catch
                rawpath = fullfile( rawpathstem, subject{ss}{1}, subject{ss}{2} );
            end
            rawfname = fullfile( rawpath, [blocksin{ss}{bb} '_raw.fif'] );

            outpath = fullfile( pathstem, subject{ss}{1} );

            outfname1 = fullfile( outpath, [blocksout{ss}{bb} '_tmp.fif'] );    % files after bad channel check
            logfname1 = fullfile( outpath, [blocksout{ss}{bb} '_tmp.log'] );

            outfname2 = fullfile( outpath, [blocksout{ss}{bb} '_sss.fif'] );    % files after SSS+ST
            logfname2 = fullfile( outpath, [blocksout{ss}{bb} '_sss.log'] );

            outfname3 = fullfile( outpath, [blocksout{ss}{bb} '_ssst.fif'] );   % files after interpolation to first specified session
            logfname3 = fullfile( outpath, [blocksout{ss}{bb} '_ssst.log'] );

            posfname = fullfile( outpath, [blocksout{ss}{bb} '_hpi.pos'] );     % HPI info

            badfname = fullfile( outpath, [blocksout{ss}{bb} '_bad.txt'] );     % bad channel info

            markbadfname = fullfile( outpath, [blocksout{ss}{bb} '_markbad.fif'] );

            if bb == 1
                b1file = outfname2 %#ok<NOPRT>
            end

            if length(badchannels)<ss, badchannels{ss} = []; end
            fprintf(1, '\n Now processing %s with %d pre-specified bad channels.\n', rawfname, length( badchannels{ss, bb} ) );

            if exist(outfname3,'file')
                disp(['File ' outfname3 ' already exists, moving on...'])
            else
                %% (2) Convert data

                skipint = '0 20';
                %         mfcmd2=[
                %         '/neuro/bin/util//neuro/bin/util/maxfilter-2.2 -f ' [rawfname] ' -o ' [outfname1],...
                %         ' -autobad 20 -skip ' [skipint] ' -v | tee ' [logfname1]
                %         ];
                %       Skipping removed, because errors, because file does not start at
                %       zero for some reason.
                mfcmd2=[
                    '/neuro/bin/util//neuro/bin/util/maxfilter -f ' [rawfname] ' -o ' [outfname1],...
                    ' -autobad 20 -v | tee ' [logfname1]
                    ];

                fprintf(1, '\n\n%s\n\n', mfcmd2);

                if exist(outfname2,'file')
                    disp(['File ' outfname2 ' already exists, moving on...'])
                else
                eval([' ! ' mfcmd2])
                delete( outfname1 );

                %% Get bad channels
                % Get bad channels from log file, store in file:
                badcmd=[
                    'cat ' [logfname1] ' | sed -n ''/Static/p'' | cut -f 5- -d '' '' > ' [badfname]
                    ];
                fprintf(1, 'Looking for bad channels\n');
                fprintf(1, '\n%s\n', badcmd);
                eval([' ! ' badcmd]);


                % Read bad channels in to matlab variable:
                fprintf(1, '\nReading bad channel information\n');
                try
                    x=dlmread([badfname],' ');
                catch %For case where no bad channels
                    x  = 0;
                end
                x=reshape(x,1,prod(size(x))); %#ok<PSIZE>
                x=x(x>0); % Omit zeros (padded by dlmread):


                % Get frequencies (number of buffers in which chan was bad):
                [frq,allbad] = hist(x,unique(x));


                % Mark bad based on threshold (currently 5 buffers):
                bads=allbad(frq>5);
                badstxt = sprintf('%s%s%s',num2str(bads)) %#ok<NOPRT>
                if sum(badstxt)>0
                    dlmwrite([markbadfname],badstxt,'delimiter',' ');
                else
                    eval(['! touch ' [markbadfname] ])
                end

                % If extra bad channels defined, append them here
                if ~isempty( badchannels{ss,bb} ),
                    for i=1:length(badchannels{ss,bb}),
                        badstxt = [badstxt ' ' badchannels{ss,bb}{i}];
                    end;
                end;
                fprintf(1, '\nThe following channels are marked as bad: %s\n\n', badstxt);

                %% (3) Maxfilter incl. ST and Movecomp
                % -- MAXFILTER ARGUMENTS --:

                % ORIGIN and FRAME:
                orgcmd=sprintf('  -frame head -origin 0 0 45');


                % BAD CHANNELS:
                if ~isempty(badstxt)
                    badcmd=['  -bad ', badstxt];
                else
                    badcmd='';
                end


                % HPI ESTIMATION/MOVEMENT COMPENSATION:
                hpistep=200;hpisubt='amp';
                hpicmd=sprintf('  -hpistep %d -hpisubt %s -movecomp -hp %s',hpistep,hpisubt,posfname);

                % SSS with ST:
                stwin=4;
                stcorr=0.980;
                stcmd=sprintf('  -st %d -corr %g',stwin,stcorr);

                % Downsampling (es edit- downsampling turned off)
                %         dsval = 4;
                %         dscmd=sprintf('  -ds %d', dsval');


                % -- MAXFILTER COMMAND --

                if exist(outfname2),
                    fprintf(1, 'Deleting %s\n', outfname2);
                    delete( outfname2 );
                end;

    % TA:            
    %             mfcmd3=[
    %                 ' /neuro/bin/util//neuro/bin/util/maxfilter -f ' [rawfname] ' -o ' [outfname2],...
    %                 '  -ctc /neuro/databases/ctc/ct_sparse.fif' ' ',...
    %                 '  -cal /neuro/databases/sss/sss_cal.dat' ' ',...
    %                 '  -autobad off ',...
    %                 stcmd,...       % temporal SSS
    %                 badcmd,...      % bad channels
    %                 orgcmd,...      % head frame and origin
    %                 hpicmd,...      % movement compensation
    %                 '  -format short ',...
    %                 '  -v | tee ' [logfname2]
    %                 ];

    % TA:
                mfcmd3 = ['/neuro/bin/util/maxfilter -f ',[rawfname],...           input % mf_str & f %TA
                    ' -o ',[outfname2],    ...   output
                    ' -ctc /neuro/databases/ctc/ct_sparse.fif',...
                    ' -cal /neuro/databases/sss/sss_cal.dat',...
                    ' -st false',...
                    ' -frame head',...
                    stcmd,...       % temporal SSS %TA
                    badcmd,...      % bad channels %TA:                 ' -autobad on',...                    autobad chans
                    orgcmd,...      % head frame and origin %TA
                    hpicmd,...      % movement compensation %TA:                 ' -movecomp ',...                     movement comp %TA:                 ' -v | tee ',[logfname2],...
                    ' -format short',...
                    ' -force',...
                    ' -corr 0.9'];




                % removed: dscmd,...       % downsampling
                % removed: '  -skip 0 20 ',...  Errors because doesn't start at zero for
                % some reason

                fprintf(1, '\nMaxfiltering... (SSS+ST)\n');
                fprintf(1, '\n\n%s\n\n', mfcmd3);
                eval([' ! ' mfcmd3 ]);
                end

                % (4) %%%%%%%%%%%%%%%%%%%%%%%%%


                % TRANSFORMATION (all but first file, block 1):
                if bb>1

                    trcmd=sprintf(' -trans %s -frame head -origin 0 0 45',b1file);

                    mfcmd4=[
                        '/neuro/bin/util//neuro/bin/util/maxfilter -f ' [outfname2] ' -o ' [outfname3],...
                        '  -autobad off ', trcmd, ' -force -v | tee ' logfname3
                        ];
                    fprintf(1, '\nMaxfiltering... -trans\n');
                    fprintf(1, '%s\n', mfcmd4);
                    eval([' ! ' mfcmd4 ])

                else

                    b1file = outfname2;     % file used for future "trans"
                    copyfile( outfname2, outfname3 );

                end;    % if bb>1
            end

            [f1,f2] = fileparts(outfname2);
            try
                spm_eeg_convert(struct('dataset',outfname2,'mode','continuous','outfile',[f1 filesep 'S_' f2],'channels','All'))
            catch err
                warning('if this person has no EEG, comment out ln 519-537 in spm_eeg_convert')
                rethrow(err)
            end

        end;    % blocks

    end;    % subjects
catch err
    return
end