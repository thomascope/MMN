%Latest stats requires Matlab2015a or greater with statistics toolbox

thisdir = pwd;


datapathstem = '/imaging/tc02/Holly_MMN/Coherence_Connectivity_ICA_LOR/crosshem/'; %With minimum norm source reconstruction
addpath(['/group/language/data/thomascope/vespa/SPM12version/Standalone preprocessing pipeline/tc_source_stats/ojwoodford-export_fig-216b30e']);
addpath(['/group/language/data/thomascope/MMN/ICA_denoise/stdshade']);

%groupstodo = {'Control' 'pca' 'bvFTD' 'nfvppa' 'MCI'};
for i = 1:length(Participant)
all_diags{i} = Participant{i}.diag;
end
groupstodo = unique(all_diags,'stable');
cmap = colormap(parula(length(groupstodo)));

these_corr_sig_pairs = {1:2; 2:3; [2,4]; 3:4; 5:6; 6:7; [6,8]; 7:8; [1,5]; [2,6]; [3,7]; [4,8]};
sources = {'left A1';
    'left STG';
    'left IFG';
    'left IPC';
    'right A1';
    'right STG';
    'right IFG';
    'right IPC'};

averagesubtracted = 1;
highfreq = 1;
timewins = [0 500];
%timewins = [32 296; 300 564; 636 900];
topfreqband = 49;

closeafter = 1;

%analysis_type = 'Granger';
analysis_type = 'icoh';

switch(analysis_type)
    case 'Granger'
        datapathstem = [datapathstem 'granger/'];
    case 'icoh'
        datapathstem = [datapathstem 'coh/'];
end

addpath('/imaging/tc02/toolboxes/rsatoolbox/Engines/')

cd(datapathstem)
filenames = {};
all_subjs = [];

runningtotal = 1;
for i = 1:length(groupstodo)
    cd(groupstodo{i});
    thesefiles = dir('*overall*mat');
    all_subjs = [all_subjs; [repmat(i,size(thesefiles)), (runningtotal:(runningtotal+size(thesefiles,1)-1))']];
    runningtotal = runningtotal+size(thesefiles,1)';
    filenames = [filenames, thesefiles.name];
    cd(datapathstem)
end

cd(thisdir)

group = all_subjs(:,1)';

for t = 1:size(timewins)
    
    start_times = timewins(t,1);
    end_times = timewins(t,2);
    
    all_granger_data = zeros(8,8,1,topfreqband,2,length(group));
    all_mismatch_contrasts = zeros(8,8,1,topfreqband,length(group));
    all_absolute_mismatch_contrasts = zeros(8,8,1,topfreqband,length(group));
    all_clarity_contrasts = zeros(8,8,1,topfreqband,length(group));
    all_absolute_clarity_contrasts = zeros(8,8,1,topfreqband,length(group));
    all_random_granger_data = zeros(8,8,1,topfreqband,2,100,length(group));
    all_random_mismatch_contrasts = zeros(8,8,1,topfreqband,100,length(group));
    all_random_clarity_contrasts = zeros(8,8,1,topfreqband,100,length(group));
    all_absolute_random_mismatch_contrasts = zeros(8,8,1,topfreqband,100,length(group));
    all_absolute_random_clarity_contrasts = zeros(8,8,1,topfreqband,100,length(group));
    all_interaction_contrasts = zeros(8,8,1,topfreqband,length(group));
    all_absolute_interaction_contrasts = zeros(8,8,1,topfreqband,length(group));
    all_random_interaction_contrasts = zeros(8,8,1,topfreqband,100,length(group));
    all_absolute_random_interaction_contrasts = zeros(8,8,1,topfreqband,100,length(group));
    all_controls_granger_data = [];
    all_patients_granger_data = [];
    all_controls_mismatch_contrasts = [];
    all_patients_mismatch_contrasts = [];
    all_controls_clarity_contrasts = [];
    all_patients_clarity_contrasts = [];
    all_controls_interaction_contrasts = [];
    all_patients_interaction_contrasts = [];
    
    demeaned_all_granger_data = zeros(8,8,1,topfreqband,2,length(group));
    demeaned_all_mismatch_contrasts = zeros(8,8,1,topfreqband,length(group));
    demeaned_all_clarity_contrasts = zeros(8,8,1,topfreqband,length(group));
    demeaned_all_controls_granger_data = [];
    demeaned_all_patients_granger_data = [];
    demeaned_all_controls_mismatch_contrasts = [];
    demeaned_all_patients_mismatch_contrasts = [];
    demeaned_all_controls_clarity_contrasts = [];
    demeaned_all_patients_clarity_contrasts = [];
    all_random_controls_granger_data = [];
    all_random_controls_mismatch_contrasts = [];
    all_random_controls_clarity_contrasts = [];
    all_random_patients_granger_data = [];
    all_random_patients_mismatch_contrasts = [];
    all_random_patients_clarity_contrasts = [];
    all_random_controls_interaction_contrasts = [];
    all_random_patients_interaction_contrasts = [];
    all_absolute_controls_mismatch_contrasts = [];
    all_absolute_patients_mismatch_contrasts = [];
    all_absolute_random_controls_mismatch_contrasts = [];
    all_absolute_random_patients_mismatch_contrasts = [];
    
    for s = 1:length(group)
        %load([datapathstem 's' num2str(s) '_evoked_grangerdata_' num2str(start_times) '_' num2str(end_times) '_overall']); %For evoked data
        if averagesubtracted == 1
            if highfreq == 1
                load([datapathstem groupstodo{group(s)} '/s' num2str(all_subjs(s,2)) '_grangerdata_highfreq_averagesubtracted_100_' num2str(start_times) '_' num2str(end_times) '_overall']);
            else
                error('I havent done this analysis')
            end
        else
            if highfreq == 1
                error('I havent done this analysis')
            else
                error('I havent done this analysis')
            end
        end
        %     if rejecteeg{s} == 1 %Biases analysis
        %         continue
        %     end
        
        %Trim final frequency (or all higher frequencies) - seems artefactually high ?filtering
        foi = foi(1:topfreqband);
        granger_data = granger_data(:,:,:,1:topfreqband,:,:);
        
        %     foi = foi(1:end-1);
        %     granger_data = granger_data(:,:,:,1:end-1,:,:);
        switch(analysis_type)
            case 'icoh'
                granger_data = abs(granger_data);
        end
        
        all_granger_data(:,:,:,:,:,s) = granger_data(:,:,:,:,:,201);
        all_random_granger_data(:,:,:,:,:,:,s) = granger_data(:,:,:,:,:,101:200); %Random permutations by trial
        demeaned_all_granger_data(:,:,:,:,:,s) = granger_data(:,:,:,:,:,201) / mean(mean(mean(nanmean(granger_data(:,:,:,:,:,201)))));
        all_mismatch_contrasts(:,:,:,:,s) = mean((granger_data(:,:,:,:,1,201)-granger_data(:,:,:,:,2,201)),5);
        all_absolute_mismatch_contrasts(:,:,:,:,s) = mean((abs(granger_data(:,:,:,:,1,201))-abs(granger_data(:,:,:,:,2,201))),5);
        all_random_mismatch_contrasts(:,:,:,:,:,s) = mean((granger_data(:,:,:,:,1,1:100)-granger_data(:,:,:,:,2,1:100)),5);
        all_absolute_random_mismatch_contrasts(:,:,:,:,:,s) = mean((abs(granger_data(:,:,:,:,1,1:100))-abs(granger_data(:,:,:,:,2,1:100))),5);
        demeaned_all_mismatch_contrasts(:,:,:,:,s) = mean((granger_data(:,:,:,:,1,201)-granger_data(:,:,:,:,2,201)),5) / mean(mean(mean(nanmean(granger_data(:,:,:,:,:,201)))));
%         all_clarity_contrasts(:,:,:,:,s) = mean((granger_data(:,:,:,:,[5,6],201)-granger_data(:,:,:,:,[1,2],201))/2,5);
%         all_absolute_clarity_contrasts(:,:,:,:,s) = mean((abs(granger_data(:,:,:,:,[5,6],201))-abs(granger_data(:,:,:,:,[1,2],201)))/2,5);
%         demeaned_all_clarity_contrasts(:,:,:,:,s) = mean((granger_data(:,:,:,:,[5,6],201)-granger_data(:,:,:,:,[1,2],201))/2,5) / mean(mean(mean(nanmean(granger_data(:,:,:,:,:,201)))));
%         all_random_clarity_contrasts(:,:,:,:,:,s) = mean((granger_data(:,:,:,:,[5,6],1:100)-granger_data(:,:,:,:,[1,2],1:100)),5);
%         all_absolute_random_clarity_contrasts(:,:,:,:,:,s) = mean((abs(granger_data(:,:,:,:,1,1:100))-abs(granger_data(:,:,:,:,2,1:100))),5);
%         all_interaction_contrasts(:,:,:,:,s) = mean((granger_data(:,:,:,:,[2,5],201)-granger_data(:,:,:,:,[1,6],201)),5); %Match4-MisMatch4+MisMatch16-Match16
%         all_absolute_interaction_contrasts(:,:,:,:,s) = mean((abs(granger_data(:,:,:,:,[2,5],201))-abs(granger_data(:,:,:,:,[1,6],201))),5);
%         all_random_interaction_contrasts(:,:,:,:,:,s) = mean((granger_data(:,:,:,:,[2,5],1:100)-granger_data(:,:,:,:,[1,6],1:100)),5);
%         all_absolute_random_interaction_contrasts(:,:,:,:,:,s) = mean((abs(granger_data(:,:,:,:,[2,5],1:100))-abs(granger_data(:,:,:,:,[1,6],1:100))),5);
        if group(s) == 1
            all_controls_granger_data(:,:,:,:,:,end+1) = granger_data(:,:,:,:,:,201);
            all_random_controls_granger_data(:,:,:,:,:,:,end+1) = granger_data(:,:,:,:,:,101:200); %Random permutations by trial
            demeaned_all_controls_granger_data(:,:,:,:,:,end+1) = granger_data(:,:,:,:,:,201)/ mean(mean(mean(nanmean(granger_data(:,:,:,:,:,201)))));
            all_controls_mismatch_contrasts(:,:,:,:,end+1) = mean((granger_data(:,:,:,:,1,201)-granger_data(:,:,:,:,2,201)),5);
            all_absolute_controls_mismatch_contrasts(:,:,:,:,end+1) = mean((abs(granger_data(:,:,:,:,1,201))-abs(granger_data(:,:,:,:,2,201))),5);
            all_absolute_random_controls_mismatch_contrasts(:,:,:,:,:,end+1) = mean((abs(granger_data(:,:,:,:,1,1:100))-abs(granger_data(:,:,:,:,2,1:100))),5);
            all_random_controls_mismatch_contrasts(:,:,:,:,:,end+1) = mean((granger_data(:,:,:,:,1,1:100)-granger_data(:,:,:,:,2,1:100)),5);
            demeaned_all_controls_mismatch_contrasts(:,:,:,:,end+1) = mean((granger_data(:,:,:,:,1,201)-granger_data(:,:,:,:,2,201)),5)/ mean(mean(mean(nanmean(granger_data(:,:,:,:,:,201)))));
%             all_controls_clarity_contrasts(:,:,:,:,end+1) = mean((granger_data(:,:,:,:,[5,6],201)-granger_data(:,:,:,:,[1,2],201))/2,5);
%             all_random_controls_clarity_contrasts(:,:,:,:,:,end+1) = mean((granger_data(:,:,:,:,[5,6],1:100)-granger_data(:,:,:,:,[1,2],1:100)),5);
%             demeaned_all_controls_clarity_contrasts(:,:,:,:,end+1) = mean((granger_data(:,:,:,:,[5,6],201)-granger_data(:,:,:,:,[1,2],201))/2,5)/ mean(mean(mean(nanmean(granger_data(:,:,:,:,:,201)))));
%             all_controls_interaction_contrasts(:,:,:,:,end+1) = mean((granger_data(:,:,:,:,[2,5],201)-granger_data(:,:,:,:,[1,6],201)),5); %Match4-MisMatch4+MisMatch16-Match16
%             all_random_controls_interaction_contrasts(:,:,:,:,:,end+1) = mean((granger_data(:,:,:,:,[2,5],1:100)-granger_data(:,:,:,:,[1,6],1:100)),5);
        elseif group(s) >= 2
            all_patients_granger_data(:,:,:,:,:,end+1) = granger_data(:,:,:,:,:,201);
            all_random_patients_granger_data(:,:,:,:,:,:,end+1) = granger_data(:,:,:,:,:,101:200); %Random permutations by trial
            demeaned_all_patients_granger_data(:,:,:,:,:,end+1) = granger_data(:,:,:,:,:,201)/ mean(mean(mean(nanmean(granger_data(:,:,:,:,:,201)))));
            all_patients_mismatch_contrasts(:,:,:,:,end+1) = mean((granger_data(:,:,:,:,1,201)-granger_data(:,:,:,:,2,201)),5);
            all_absolute_patients_mismatch_contrasts(:,:,:,:,end+1) = mean((abs(granger_data(:,:,:,:,1,201))-abs(granger_data(:,:,:,:,2,201))),5);
            all_absolute_random_patients_mismatch_contrasts(:,:,:,:,:,end+1) = mean((abs(granger_data(:,:,:,:,1,1:100))-abs(granger_data(:,:,:,:,2,1:100))),5);
            all_random_patients_mismatch_contrasts(:,:,:,:,:,end+1) = mean((granger_data(:,:,:,:,1,1:100)-granger_data(:,:,:,:,2,1:100)),5);
            demeaned_all_patients_mismatch_contrasts(:,:,:,:,end+1) = mean((granger_data(:,:,:,:,1,201)-granger_data(:,:,:,:,2,201)),5)/ mean(mean(mean(nanmean(granger_data(:,:,:,:,:,201)))));
%             all_patients_clarity_contrasts(:,:,:,:,end+1) = mean((granger_data(:,:,:,:,[5,6],201)-granger_data(:,:,:,:,[1,2],201))/2,5);
%             all_random_patients_clarity_contrasts(:,:,:,:,:,end+1) = mean((granger_data(:,:,:,:,[5,6],1:100)-granger_data(:,:,:,:,[1,2],1:100)),5);
%             demeaned_all_patients_clarity_contrasts(:,:,:,:,end+1) = mean((granger_data(:,:,:,:,[5,6],201)-granger_data(:,:,:,:,[1,2],201))/2,5)/ mean(mean(mean(nanmean(granger_data(:,:,:,:,:,201)))));
%             all_patients_interaction_contrasts(:,:,:,:,end+1) = mean((granger_data(:,:,:,:,[2,5],201)-granger_data(:,:,:,:,[1,6],201)),5); %Match4-MisMatch4+MisMatch16-Match16
%             all_random_patients_interaction_contrasts(:,:,:,:,:,end+1) = mean((granger_data(:,:,:,:,[2,5],1:100)-granger_data(:,:,:,:,[1,6],1:100)),5);
        end
        
    end
    
     %Account for stupid zero matrices in first line of end+1
    these_lines_to_strip = {
        'all_controls_granger_data(:,:,:,:,:,end+1) = granger_data(:,:,:,:,:,201);'
        'all_random_controls_granger_data(:,:,:,:,:,:,end+1) = granger_data(:,:,:,:,:,101:200);'
        'demeaned_all_controls_granger_data(:,:,:,:,:,end+1) = granger_data(:,:,:,:,:,201)/ mean(mean(mean(nanmean(granger_data(:,:,:,:,:,201)))));'
        'all_controls_mismatch_contrasts(:,:,:,:,end+1) = mean((granger_data(:,:,:,:,1,201)-granger_data(:,:,:,:,2,201)),5);'
        'all_absolute_controls_mismatch_contrasts(:,:,:,:,end+1) = mean((abs(granger_data(:,:,:,:,1,201))-abs(granger_data(:,:,:,:,2,201))),5);'
        'all_absolute_random_controls_mismatch_contrasts(:,:,:,:,:,end+1) = mean((abs(granger_data(:,:,:,:,1,1:100))-abs(granger_data(:,:,:,:,2,1:100))),5);'
        'all_random_controls_mismatch_contrasts(:,:,:,:,:,end+1) = mean((granger_data(:,:,:,:,1,1:100)-granger_data(:,:,:,:,2,1:100)),5);'
        'demeaned_all_controls_mismatch_contrasts(:,:,:,:,end+1) = mean((granger_data(:,:,:,:,1,201)-granger_data(:,:,:,:,2,201)),5)/ mean(mean(mean(nanmean(granger_data(:,:,:,:,:,201)))));'
        'all_controls_clarity_contrasts(:,:,:,:,end+1) = mean((granger_data(:,:,:,:,[5,6],201)-granger_data(:,:,:,:,[1,2],201))/2,5);'
        'all_random_controls_clarity_contrasts(:,:,:,:,:,end+1) = mean((granger_data(:,:,:,:,[5,6],1:100)-granger_data(:,:,:,:,[1,2],1:100)),5);'
        'demeaned_all_controls_clarity_contrasts(:,:,:,:,end+1) = mean((granger_data(:,:,:,:,[5,6],201)-granger_data(:,:,:,:,[1,2],201))/2,5)/ mean(mean(mean(nanmean(granger_data(:,:,:,:,:,201)))));'
        'all_patients_granger_data(:,:,:,:,:,end+1) = granger_data(:,:,:,:,:,201);'
        'all_random_patients_granger_data(:,:,:,:,:,:,end+1) = granger_data(:,:,:,:,:,101:200);'
        'demeaned_all_patients_granger_data(:,:,:,:,:,end+1) = granger_data(:,:,:,:,:,201)/ mean(mean(mean(nanmean(granger_data(:,:,:,:,:,201)))));'
        'all_patients_mismatch_contrasts(:,:,:,:,end+1) = mean((granger_data(:,:,:,:,1,201)-granger_data(:,:,:,:,2,201)),5);'
        'all_absolute_patients_mismatch_contrasts(:,:,:,:,end+1) = mean((abs(granger_data(:,:,:,:,1,201))-abs(granger_data(:,:,:,:,2,201))),5);'
        'all_absolute_random_patients_mismatch_contrasts(:,:,:,:,:,end+1) = mean((abs(granger_data(:,:,:,:,1,1:100))-abs(granger_data(:,:,:,:,2,1:100))),5);'
        'all_random_patients_mismatch_contrasts(:,:,:,:,:,end+1) = mean((granger_data(:,:,:,:,1,1:100)-granger_data(:,:,:,:,2,1:100)),5);'
        'demeaned_all_patients_mismatch_contrasts(:,:,:,:,end+1) = mean((granger_data(:,:,:,:,1,201)-granger_data(:,:,:,:,2,201)),5)/ mean(mean(mean(nanmean(granger_data(:,:,:,:,:,201)))));'
        'all_patients_clarity_contrasts(:,:,:,:,end+1) = mean((granger_data(:,:,:,:,[5,6],201)-granger_data(:,:,:,:,[1,2],201))/2,5);'
        'all_random_patients_clarity_contrasts(:,:,:,:,:,end+1) = mean((granger_data(:,:,:,:,[5,6],1:100)-granger_data(:,:,:,:,[1,2],1:100)),5);'
        'demeaned_all_patients_clarity_contrasts(:,:,:,:,end+1) = mean((granger_data(:,:,:,:,[5,6],201)-granger_data(:,:,:,:,[1,2],201))/2,5)/ mean(mean(mean(nanmean(granger_data(:,:,:,:,:,201)))));'
        'all_controls_interaction_contrasts(:,:,:,:,end+1) = mean((granger_data(:,:,:,:,[2,5],201)-granger_data(:,:,:,:,[1,6],201)),5);'
        'all_patients_interaction_contrasts(:,:,:,:,end+1) = mean((granger_data(:,:,:,:,[2,5],201)-granger_data(:,:,:,:,[1,6],201)),5);'
        'all_random_controls_interaction_contrasts(:,:,:,:,:,end+1) = mean((granger_data(:,:,:,:,[2,5],1:100)-granger_data(:,:,:,:,[1,6],1:100)),5);'
        'all_random_patients_interaction_contrasts(:,:,:,:,:,end+1) = mean((granger_data(:,:,:,:,[2,5],1:100)-granger_data(:,:,:,:,[1,6],1:100)),5);'
        };
    
    for this_line = 1:size(these_lines_to_strip,1)
        target_string = '(:,:,';
        variable_blocks = strsplit(these_lines_to_strip{this_line},target_string);
        variable_indices = strsplit(variable_blocks{2},'=');
        variable_indices{1}(ismember(variable_indices{1},' ')) = [];
        eval([variable_blocks{1} '=' variable_blocks{1} target_string variable_indices{1}(1:end-6) '2:end);'])
    end
    
    for this_pair = 1:length(these_corr_sig_pairs)
        from = these_corr_sig_pairs{this_pair}(1);
        to = these_corr_sig_pairs{this_pair}(2);
        from_name = sources{from};
        to_name = sources{to};
        
        if closeafter == 1; close all; end; figure
        linehandle(1) = stdshade_TEC(squeeze(mean(all_granger_data(from,to,:,:,:,:),5))',0.2,'g',foi,1,1);
        hold on
        linehandle(2) = stdshade_TEC(squeeze(mean(all_granger_data(to,from,:,:,:,:),5))',0.2,'b',foi,1,1);
        plot(foi,squeeze(mean(mean(all_random_granger_data(from,to,:,:,:,:,:),5),7)),'k')
        plot(foi,squeeze(mean(mean(all_random_granger_data(to,from,:,:,:,:,:),5),7)),'k')
        stdshade_TEC(squeeze(mean(all_granger_data(from,to,:,:,:,:),5))',0.2,'g',foi,1,1)
        stdshade_TEC(squeeze(mean(all_granger_data(to,from,:,:,:,:),5))',0.2,'b',foi,1,1)
        title(['All ' analysis_type ' time ' num2str(start_times) '-' num2str(end_times)])
        legend(linehandle(1:2), {[from_name '-' to_name],[to_name '-' from_name]})
        for i = 1:topfreqband %Compare across group permutation
            p_tf(i) = 1-relRankIn_includeValue_lowerBound(squeeze(mean(mean(all_random_granger_data(from,to,:,i,:,:,:),5),7)),squeeze(mean(mean(all_granger_data(from,to,:,i,:,:),5),6)))
            p_ft(i) = 1-relRankIn_includeValue_lowerBound(squeeze(mean(mean(all_random_granger_data(to,from,:,i,:,:,:),5),7)),squeeze(mean(mean(all_granger_data(to,from,:,i,:,:),5),6)))
            if p_tf(i)<=0.05 || p_tf(i)>=0.95
                plot(foi(i),0,'g*')
            end
            if p_ft(i)<=0.05 || p_ft(i)>=0.95
                plot(foi(i),0,'bx')
            end
        end
        savestring = ['./figures/' from_name '_' to_name '_All_' analysis_type '_time_' num2str(start_times) '-' num2str(end_times) '.pdf'];
        savestring = strrep(savestring,' ','_');
        eval(['export_fig ' savestring ' -transparent']); eval(['export_fig ' savestring(1:end-3) 'png -transparent']);
        
        
        switch(analysis_type)
            case 'Granger'
                
                
                if closeafter == 1; close all; end; figure
                linehandle(1) = stdshade_TEC(squeeze(mean(all_granger_data(from,to,:,:,:,:),5))',0.2,'g',foi,1,1);
                hold on
                linehandle(2) = stdshade_TEC(squeeze(mean(all_granger_data(to,from,:,:,:,:),5))',0.2,'b',foi,1,1);
                title(['By direction ' analysis_type ' time ' num2str(start_times) '-' num2str(end_times)])
                legend(linehandle(1:2),{[from_name '-' to_name],[to_name '-' from_name]})
                for i = 1:topfreqband
                    thisarray = [repmat(sort(repmat([1:(runningtotal-1)],1,2))',2,1),repmat([1:2]',(runningtotal-1)*2,1),[ones((runningtotal-1)*2,1);2*ones((runningtotal-1)*2,1)],[reshape(all_granger_data(from,to,:,i,:,:),[2*(runningtotal-1),1]);reshape(all_granger_data(to,from,:,i,:,:),[2*(runningtotal-1),1])]];
                    t1=array2table(thisarray,'VariableNames',{'Subject','Condition','Direction','Response'});
                    thistest = fitglm(t1,'Response ~ Direction + Condition + Subject');
                    p_tf(i) = thistest.Coefficients{4,4};
                    if p_tf(i) < 0.05
                        if mean(mean(all_granger_data(from,to,:,i,:,:),5),6)>mean(mean(all_granger_data(to,from,:,i,:,:),5),6)
                            plot(foi(i),0,'g*')
                        else
                            plot(foi(i),0,'bx')
                        end
                    end
                end
                savestring = ['./figures/' from_name '_' to_name '_All_mean_' analysis_type '_time_' num2str(start_times) '-' num2str(end_times) '.pdf'];
                savestring = strrep(savestring,' ','_');
                eval(['export_fig ' savestring ' -transparent']); eval(['export_fig ' savestring(1:end-3) 'png -transparent']);
                
                %             if closeafter == 1; close all; end; figure
                %             stdshade_TEC(squeeze(median(all_granger_data(from,to,:,:,:,:),5))',0.2,'g',foi,1,1)
                %             hold on
                %             stdshade_TEC(squeeze(median(all_granger_data(to,from,:,:,:,:),5))',0.2,'b',foi,1,1)
                %             plot(foi,squeeze(median(mean(all_random_granger_data(from,to,:,:,:,:,:),5),7)),'k')
                %             plot(foi,squeeze(median(mean(all_random_granger_data(to,from,:,:,:,:,:),5),7)),'k')
                %             stdshade_TEC(squeeze(median(all_granger_data(from,to,:,:,:,:),5))',0.2,'g',foi,1,1)
                %             stdshade_TEC(squeeze(median(all_granger_data(to,from,:,:,:,:),5))',0.2,'b',foi,1,1)
                %             title(['All median ' analysis_type ' time ' num2str(start_times) '-' num2str(end_times)])
                %             legend({[from_name '-' to_name],[to_name '-' from_name]})
                %             for i = 1:topfreqband %Compare across group permutation
                %                 p_tf(i) = 1-relRankIn_includeValue_lowerBound(squeeze(median(mean(all_random_granger_data(from,to,:,i,:,:,:),5),7)),squeeze(median(mean(all_granger_data(from,to,:,i,:,:),5),6)))
                %                 p_ft(i) = 1-relRankIn_includeValue_lowerBound(squeeze(median(mean(all_random_granger_data(to,from,:,i,:,:,:),5),7)),squeeze(median(mean(all_granger_data(to,from,:,i,:,:),5),6)))
                %                 if p_tf(i)<=0.05 || p_tf(i)>=0.95
                %                     plot(foi(i),0,'g*')
                %                 end
                %                 if p_ft(i)<=0.05 || p_ft(i)>=0.95
                %                     plot(foi(i),0,'bx')
                %                 end
                %             end
                %             savestring = ['./figures/' from_name '_' to_name '_All_median_' analysis_type '_time_' num2str(start_times) '-' num2str(end_times) '.pdf'];
                %             savestring = strrep(savestring,' ','_');
                %             eval(['export_fig ' savestring ' -transparent']); eval(['export_fig ' savestring(1:end-3) 'png -transparent']);
                
                if closeafter == 1; close all; end; figure
                linehandle(1) = stdshade_TEC(squeeze(mean(demeaned_all_granger_data(from,to,:,:,:,:),5))',0.2,'g',foi,1,1);
                hold on
                linehandle(2) = stdshade_TEC(squeeze(mean(demeaned_all_granger_data(to,from,:,:,:,:),5))',0.2,'b',foi,1,1);
                title(['All demeaned ' analysis_type ' time ' num2str(start_times) '-' num2str(end_times)])
                legend(linehandle(1:2),{[from_name '-' to_name],[to_name '-' from_name]}) %Not yet clear how to do stats on this - XXX
                %             for i = 1:topfreqband %Rough initial parametric stats - better to use permutation
                %                 %[h_tf(i) p_tf(i)] = ttest(mean(demeaned_all_granger_data(from,to,:,i,:,:),5)-mean(demeaned_all_granger_data(to,from,:,i,:,:),5));
                %                 [h_tf(i) p_tf(i)] = ttest(reshape(demeaned_all_granger_data(from,to,:,i,:,:),[126,1]),reshape(demeaned_all_granger_data(to,from,:,i,:,:),[126,1]))
                %                 if h_tf(i) == 1
                %                     if mean(mean(demeaned_all_granger_data(from,to,:,i,:,:),5),6)>mean(mean(demeaned_all_granger_data(to,from,:,i,:,:),5),6)
                %                         plot(foi(i),0,'g*')
                %                     else
                %                         plot(foi(i),0,'bx')
                %                     end
                %                 end
                %             end
                % Fit GLMs with subject and condition. NB: I know that this
                % isn't quite right as subject isn't treated as a random
                % effect, but I can't figure out how to do this and the
                % Bonferroni correction is stringent and should more than
                % compensate
                for i = 1:topfreqband
                    thisarray = [repmat(sort(repmat([1:(runningtotal-1)],1,2))',2,1),repmat([1:2]',(runningtotal-1)*2,1),[ones((runningtotal-1)*2,1);2*ones((runningtotal-1)*2,1)],[reshape(demeaned_all_granger_data(from,to,:,i,:,:),[(runningtotal-1)*2,1]);reshape(demeaned_all_granger_data(to,from,:,i,:,:),[(runningtotal-1)*2,1])]];
                    t1=array2table(thisarray,'VariableNames',{'Subject','Condition','Direction','Response'});
                    thistest = fitglm(t1,'Response ~ Direction + Condition + Subject');
                    p_tf(i) = thistest.Coefficients{4,4};
                    SEM(i) = thistest.Coefficients{4,2};
                    if p_tf(i) < 0.05
                        if mean(mean(demeaned_all_granger_data(from,to,:,i,:,:),5),6)>mean(mean(demeaned_all_granger_data(to,from,:,i,:,:),5),6)
                            plot(foi(i),0,'g*')
                        else
                            plot(foi(i),0,'bx')
                        end
                    end
                end
                savestring = ['./figures/' from_name '_' to_name '_All_demeaned_' analysis_type '_time_' num2str(start_times) '-' num2str(end_times) '.pdf'];
                savestring = strrep(savestring,' ','_');
                eval(['export_fig ' savestring ' -transparent']); eval(['export_fig ' savestring(1:end-3) 'png -transparent']);
                
                if closeafter == 1; close all; end; figure
                difference_demeaned = demeaned_all_granger_data(from,to,:,:,:,:)-demeaned_all_granger_data(to,from,:,:,:,:);
                hold on
                %fill([foi,fliplr(foi)],[squeeze(mean(mean(difference_demeaned(1,1,:,:,:,:),5),6))-SEM';flipud(squeeze(mean(mean(difference_demeaned(1,1,:,:,:,:),5),6))+SEM')]','g');
                stdshade_TEC(squeeze(mean(difference_demeaned(1,1,:,:,:,:),5))',0.2,'b',foi,1,1);
                %ylim([-1 1]);
                plot([0 40],[0 0],'k--','LineWidth',1);
                savestring = ['./figures/' from_name '_' to_name '_Difference_demeaned_' analysis_type '_time_' num2str(start_times) '-' num2str(end_times) '.pdf'];
                title(['Direction difference ' analysis_type ' ' from_name ' to ' to_name])
                savestring = strrep(savestring,' ','_');
                eval(['export_fig ' savestring ' -transparent']); eval(['export_fig ' savestring(1:end-3) 'png -transparent']);
                
        end
    
    %             if closeafter == 1; close all; end; figure
    %             stdshade_TEC(squeeze(mean(demeaned_all_controls_granger_data(from,to,:,:,:,:),5))',0.2,'g:',foi,1,1)
    %             hold on
    %             stdshade_TEC(squeeze(mean(demeaned_all_patients_granger_data(from,to,:,:,:,:),5))',0.2,'g--',foi,1,1)
    %             stdshade_TEC(squeeze(mean(demeaned_all_controls_granger_data(to,from,:,:,:,:),5))',0.2,'b:',foi,1,1)
    %             stdshade_TEC(squeeze(mean(demeaned_all_patients_granger_data(to,from,:,:,:,:),5))',0.2,'b--',foi,1,1)
    %             legend({[from_name '-' to_name ' control'],[from_name '-' to_name ' patient'],[to_name '-' from_name ' control'],[to_name '-' from_name ' patient']})
    %             title(['By group demeaned ' analysis_type ' time ' num2str(start_times) '-' num2str(end_times)])
    %             legend({[from_name '-' to_name],[to_name '-' from_name]}) %Not yet clear how to do stats on this - XXX
    %             for i = 1:topfreqband %Rough initial parametric stats - better to use permutation
    %                 [h_tf(i) p_tf(i)] = ttest2(mean(demeaned_all_controls_granger_data(from,to,:,i,:,:),5),mean(demeaned_all_patients_granger_data(from,to,:,i,:,:),5));
    %                 [h_ft(i) p_ft(i)] = ttest2(mean(demeaned_all_controls_granger_data(to,from,:,i,:,:),5),mean(demeaned_all_patients_granger_data(to,from,:,i,:,:),5));
    %                 if h_tf(i) == 1
    %                     plot(foi(i),0,'g*')
    %                 end
    %                 if h_ft(i) == 1
    %                     plot(foi(i),0,'bx')
    %                 end
    %             end
    
    
    if closeafter == 1; close all; end; figure
    hold on
    legend_text = cell(1,2*length(groupstodo));
    for grp = 1:2:2*length(groupstodo)
        linehandle(grp) = stdshade_TEC_cmap(squeeze(mean(all_granger_data(from,to,:,:,:,group==ceil(grp/2)),5))',0.2,cmap(ceil(grp/2),:),foi,1,1,':');
        linehandle(grp+1) = stdshade_TEC_cmap(squeeze(mean(all_granger_data(to,from,:,:,:,group==ceil(grp/2)),5))',0.2,cmap(ceil(grp/2),:),foi,1,1,'--');
        legend_text{grp} = [from_name '-' to_name ' ' groupstodo{ceil(grp/2)}]; 
        legend_text{grp+1} = [to_name '-' from_name ' ' groupstodo{ceil(grp/2)}];
    end
    legend(linehandle,legend_text)
    title(['By group ' analysis_type ' time ' num2str(start_times) '-' num2str(end_times)])
    
    for i = 1:topfreqband %Rough initial parametric stats - better to use permutation
        % all_granger_data = zeros(8,8,1,topfreqband,2,length(group));
        %Key: regions, regions, timepoints, frequencies, conditions, group
        p_tf(i) = anova1(squeeze(mean(all_granger_data(from,to,:,i,:,:),5)),group,'off');
        p_ft(i) = anova1(squeeze(mean(all_granger_data(to,from,:,i,:,:),5)),group,'off');
        if p_tf(i) <= 0.05
            plot(foi(i),0,'g*')
        end
        if p_ft(i) <= 0.05
            plot(foi(i),0,'bx')
        end
    end
    
%     for i = 1:topfreqband %Rough initial parametric stats - better to use permutation
%         [h_tf(i) p_tf(i)] = ttest2(mean(all_controls_granger_data(from,to,:,i,:,:),5),mean(all_patients_granger_data(from,to,:,i,:,:),5));
%         [h_ft(i) p_ft(i)] = ttest2(mean(all_controls_granger_data(to,from,:,i,:,:),5),mean(all_patients_granger_data(to,from,:,i,:,:),5));
%         if h_tf(i) == 1
%             plot(foi(i),0,'g*')
%         end
%         if h_ft(i) == 1
%             plot(foi(i),0,'bx')
%         end
%     end
    
    
    %             control_test = [];
    %             patient_test = [];
    %             for i = 2:34 % A little bit of 'frequency smoothing'
    %                 control_test(i-1,:) = squeeze(mean(mean(all_controls_granger_data(from,to,:,i-1:i+1,:,:),5),4));
    %                 patient_test(i-1,:) = squeeze(mean(mean(all_patients_granger_data(from,to,:,i-1:i+1,:,:),5),4));
    %                 [h_tf(i) p_tf(i)] = ttest2(control_test(i-1,:),patient_test(i-1,:));
    %                 control_test(i-1,:) = squeeze(mean(mean(all_controls_granger_data(to,from,:,i-1:i+1,:,:),5),4));
    %                 patient_test(i-1,:) = squeeze(mean(mean(all_patients_granger_data(to,from,:,i-1:i+1,:,:),5),4));
    %                 [h_ft(i) p_ft(i)] = ttest2(control_test(i-1,:),patient_test(i-1,:));
    %
    %
    %                 if h_tf(i) == 1
    %                     plot(foi(i),0,'g*')
    %                 end
    %                 if h_ft(i) == 1
    %                     plot(foi(i),0,'bx')
    %                 end
    %             end
    savestring = ['./figures/' from_name '_' to_name '_By_group_mean_' analysis_type '_time_' num2str(start_times) '-' num2str(end_times) '.pdf'];
    savestring = strrep(savestring,' ','_');
    eval(['export_fig ' savestring ' -transparent']); eval(['export_fig ' savestring(1:end-3) 'png -transparent']);
    
    
%     if closeafter == 1; close all; end; figure
%     hold on
%     patient_SEMs = std(squeeze(mean(all_patients_granger_data(from,to,:,:,:,2:end),5))')/sqrt(10);
%     control_SEMs = std(squeeze(mean(all_controls_granger_data(from,to,:,:,:,2:end),5))')/sqrt(11);
%     fill([foi,fliplr(foi)],[squeeze(mean(mean(all_controls_granger_data(to,from,:,:,:,2:end),5),6))-control_SEMs';flipud(squeeze(mean(mean(all_controls_granger_data(to,from,:,:,:,2:end),5),6))+control_SEMs')]','c');
%     fill([foi,fliplr(foi)],[squeeze(mean(mean(all_patients_granger_data(to,from,:,:,:,2:end),5),6))-patient_SEMs';flipud(squeeze(mean(mean(all_patients_granger_data(to,from,:,:,:,2:end),5),6))+patient_SEMs')]','m');
%     stdshade_TEC(squeeze(mean(all_controls_granger_data(to,from,:,:,:,2:end),5))',0.2,'k',foi,1,1)
%     stdshade_TEC(squeeze(mean(all_patients_granger_data(to,from,:,:,:,2:end),5))',0.2,'y',foi,1,1)
%     legend({[from_name '-' to_name ' control'],[from_name '-' to_name ' patient'],[to_name '-' from_name ' control'],[to_name '-' from_name ' patient']})
%     title(['By group ' analysis_type ' time ' num2str(start_times) '-' num2str(end_times)])
%     for i = 1:topfreqband %Rough initial parametric stats - better to use permutation
%         [h_tf(i) p_tf(i)] = ttest2(mean(all_controls_granger_data(from,to,:,i,:,:),5),mean(all_patients_granger_data(from,to,:,i,:,:),5));
%         [h_ft(i) p_ft(i)] = ttest2(mean(all_controls_granger_data(to,from,:,i,:,:),5),mean(all_patients_granger_data(to,from,:,i,:,:),5));
%         if h_tf(i) == 1
%             plot(foi(i),0,'g*')
%         end
%         if h_ft(i) == 1
%             plot(foi(i),0,'bx')
%         end
%     end
%     savestring = ['./figures/' from_name '_' to_name '_By_group_mean_withSEM_' analysis_type '_time_' num2str(start_times) '-' num2str(end_times) '.pdf'];
%     savestring = strrep(savestring,' ','_');
%     eval(['export_fig ' savestring ' -transparent']); eval(['export_fig ' savestring(1:end-3) 'png -transparent']);
%     
    
%     if closeafter == 1; close all; end; figure
%     hold on
%     stdshade_TEC(squeeze(median(all_controls_granger_data(from,to,:,:,:,:),5))',0.2,'g:',foi,1,1)
%     stdshade_TEC(squeeze(median(all_patients_granger_data(from,to,:,:,:,:),5))',0.2,'g--',foi,1,1)
%     stdshade_TEC(squeeze(median(all_controls_granger_data(to,from,:,:,:,:),5))',0.2,'b:',foi,1,1)
%     stdshade_TEC(squeeze(median(all_patients_granger_data(to,from,:,:,:,:),5))',0.2,'b--',foi,1,1)
%     legend({[from_name '-' to_name ' control'],[from_name '-' to_name ' patient'],[to_name '-' from_name ' control'],[to_name '-' from_name ' patient']})
%     title(['By group median ' analysis_type ' time ' num2str(start_times) '-' num2str(end_times)])
%     for i = 1:topfreqband %Rough initial parametric stats - better to use permutation
%         [p_tf(i) h_tf(i)] = ranksum(squeeze(mean(all_controls_granger_data(from,to,:,i,:,:),5)),squeeze(mean(all_patients_granger_data(from,to,:,i,:,:),5)));
%         [p_ft(i) h_ft(i)] = ranksum(squeeze(mean(all_controls_granger_data(to,from,:,i,:,:),5)),squeeze(mean(all_patients_granger_data(to,from,:,i,:,:),5)));
%         if h_tf(i) == 1
%             plot(foi(i),0,'g*')
%         end
%         if h_ft(i) == 1
%             plot(foi(i),0,'bx')
%         end
%     end
%     
%     savestring = ['./figures/' from_name '_' to_name '_By_group_median_' analysis_type '_time_' num2str(start_times) '-' num2str(end_times) '.pdf'];
%     savestring = strrep(savestring,' ','_');
%     eval(['export_fig ' savestring ' -transparent']); eval(['export_fig ' savestring(1:end-3) 'png -transparent']);
    
    
%     if closeafter == 1; close all; end; figure
%     hold on
%     stdshade_TEC(squeeze(mean(demeaned_all_controls_granger_data(from,to,:,:,:,2:end),5))',0.2,'g:',foi,1,1)
%     stdshade_TEC(squeeze(mean(demeaned_all_patients_granger_data(from,to,:,:,:,2:end),5))',0.2,'g--',foi,1,1)
%     stdshade_TEC(squeeze(mean(demeaned_all_controls_granger_data(to,from,:,:,:,2:end),5))',0.2,'b:',foi,1,1)
%     stdshade_TEC(squeeze(mean(demeaned_all_patients_granger_data(to,from,:,:,:,2:end),5))',0.2,'b--',foi,1,1)
%     
%     title(['By group demeaned ' analysis_type ' time ' num2str(start_times) '-' num2str(end_times)])
%     legend({[from_name '-' to_name ' control'],[from_name '-' to_name ' patient'],[to_name '-' from_name ' control'],[to_name '-' from_name ' patient']})
%     
    %             con_tf_lowfreq = squeeze(mean(median(mean(all_controls_granger_data(from,to,:,4:5,:,:),5),6),4));
    %             pat_tf_lowfreq = squeeze(mean(median(mean(all_patients_granger_data(from,to,:,4:5,:,:),5),6),4));
    %             con_ft_lowfreq = squeeze(mean(median(mean(all_controls_granger_data(to,from,:,4:5,:,:),5),6),4));
    %             pat_ft_lowfreq = squeeze(mean(median(mean(all_patients_granger_data(to,from,:,4:5,:,:),5),6),4));
    %
    %             con_tf_highfreq = squeeze(mean(median(mean(all_controls_granger_data(from,to,:,7:10,:,:),5),6),4));
    %             pat_tf_highfreq = squeeze(mean(median(mean(all_patients_granger_data(from,to,:,7:10,:,:),5),6),4));
    %             con_ft_highfreq = squeeze(mean(median(mean(all_controls_granger_data(to,from,:,7:10,:,:),5),6),4));
    %             pat_ft_highfreq = squeeze(mean(median(mean(all_patients_granger_data(to,from,:,7:10,:,:),5),6),4));
    
%     continue
%     
%     config = figure;
%     patfig = figure;
%     con_contrastfig = figure;
%     pat_contrastfig = figure;
%     con_interactionfig = figure;
%     pat_interactionfig = figure;
%     this_con = 0;
%     this_pat = 0;
%     
%     for s = 1:length(group)
%         if group(s) == 1
%             this_con = this_con+1;
%             figure(config);
%             h = subplot(4,3,this_con);
%             stdshade_TEC(squeeze(mean(demeaned_all_patients_granger_data(to,from,:,:,:,2:end),5))',0.2,'b--',foi,1,1)
%             hold on
%             stdshade_TEC(squeeze(mean(demeaned_all_patients_granger_data(to,from,:,:,:,2:end),5))',0.2,'b--',foi,1,1)
%             for i = 1:topfreqband %Compare within subj permutation
%                 p_tf(i) = 1-relRankIn_includeValue_lowerBound(squeeze(mean(all_random_granger_data(from,to,:,i,:,:,s),5)),squeeze(mean(all_granger_data(from,to,:,i,:,s),5)));
%                 p_ft(i) = 1-relRankIn_includeValue_lowerBound(squeeze(mean(all_random_granger_data(to,from,:,i,:,:,s),5)),squeeze(mean(all_granger_data(to,from,:,i,:,s),5)));
%                 if p_tf(i)<=0.05 || p_tf(i)>=0.95
%                     plot(foi(i),0,'g*')
%                 end
%                 if p_ft(i)<=0.05 || p_ft(i)>=0.95
%                     plot(foi(i),0,'bx')
%                 end
%             end
%             
%             figure(con_contrastfig);
%             h = subplot(4,3,this_con);
%             hold on
%             stdshade_TEC(squeeze(all_mismatch_contrasts(from,to,:,:,s))',0.2,'g',foi,1,1)
%             
%             stdshade_TEC(squeeze(all_mismatch_contrasts(to,from,:,:,s))',0.2,'b',foi,1,1)
%             plot(foi,squeeze(mean(all_random_mismatch_contrasts(from,to,:,:,:,s),6)),'k')
%             plot(foi,squeeze(mean(all_random_mismatch_contrasts(to,from,:,:,:,s),6)),'k')
%             stdshade_TEC(squeeze(all_mismatch_contrasts(from,to,:,:,s))',0.2,'g',foi,1,1)
%             
%             stdshade_TEC(squeeze(all_mismatch_contrasts(to,from,:,:,s))',0.2,'b',foi,1,1)
%             
%             %legend({[from_name '-' to_name],[to_name '-' from_name]})
%             for i = 1:topfreqband %Compare across group permutation
%                 p_tf(i) = 1-relRankIn_includeValue_lowerBound(squeeze(mean(all_random_mismatch_contrasts(from,to,:,i,:,s),6)),squeeze(mean(all_mismatch_contrasts(from,to,:,i,s),5)));
%                 p_ft(i) = 1-relRankIn_includeValue_lowerBound(squeeze(mean(all_random_mismatch_contrasts(to,from,:,i,:,s),6)),squeeze(mean(all_mismatch_contrasts(to,from,:,i,s),5)));
%                 if p_tf(i)<=0.05 || p_tf(i)>=0.95
%                     plot(foi(i),0,'g*')
%                 end
%                 if p_ft(i)<=0.05 || p_ft(i)>=0.95
%                     plot(foi(i),0,'bx')
%                 end
%             end
%             
%             figure(con_interactionfig);
%             h = subplot(4,3,this_con);
%             hold on
%             stdshade_TEC(squeeze(all_interaction_contrasts(from,to,:,:,s))',0.2,'g',foi,1,1)
%             
%             stdshade_TEC(squeeze(all_interaction_contrasts(to,from,:,:,s))',0.2,'b',foi,1,1)
%             plot(foi,squeeze(mean(all_random_interaction_contrasts(from,to,:,:,:,s),6)),'k')
%             plot(foi,squeeze(mean(all_random_interaction_contrasts(to,from,:,:,:,s),6)),'k')
%             stdshade_TEC(squeeze(all_interaction_contrasts(from,to,:,:,s))',0.2,'g',foi,1,1)
%             
%             stdshade_TEC(squeeze(all_interaction_contrasts(to,from,:,:,s))',0.2,'b',foi,1,1)
%             
%             %legend({[from_name '-' to_name],[to_name '-' from_name]})
%             for i = 1:topfreqband %Compare across group permutation
%                 p_tf(i) = 1-relRankIn_includeValue_lowerBound(squeeze(mean(all_random_interaction_contrasts(from,to,:,i,:,s),6)),squeeze(mean(all_interaction_contrasts(from,to,:,i,s),5)));
%                 p_ft(i) = 1-relRankIn_includeValue_lowerBound(squeeze(mean(all_random_interaction_contrasts(to,from,:,i,:,s),6)),squeeze(mean(all_interaction_contrasts(to,from,:,i,s),5)));
%                 if p_tf(i)<=0.05 || p_tf(i)>=0.95
%                     plot(foi(i),0,'g*')
%                 end
%                 if p_ft(i)<=0.05 || p_ft(i)>=0.95
%                     plot(foi(i),0,'bx')
%                 end
%             end
%             
%             
%             
%         elseif group(s) == 2
%             this_pat = this_pat+1;
%             figure(patfig);
%             h = subplot(4,3,this_pat);
%             stdshade_TEC(squeeze(all_interaction_contrasts(to,from,:,:,s))',0.2,'b',foi,1,1)
%             hold on
%             stdshade_TEC(squeeze(all_interaction_contrasts(to,from,:,:,s))',0.2,'b',foi,1,1)
%             for i = 1:topfreqband %Compare within subj permutation
%                 p_tf(i) = 1-relRankIn_includeValue_lowerBound(squeeze(mean(all_random_granger_data(from,to,:,i,:,:,s),5)),squeeze(mean(all_granger_data(from,to,:,i,:,s),5)));
%                 p_ft(i) = 1-relRankIn_includeValue_lowerBound(squeeze(mean(all_random_granger_data(to,from,:,i,:,:,s),5)),squeeze(mean(all_granger_data(to,from,:,i,:,s),5)));
%                 if p_tf(i)<=0.05 || p_tf(i)>=0.95
%                     plot(foi(i),0,'g*')
%                 end
%                 if p_ft(i)<=0.05 || p_ft(i)>=0.95
%                     plot(foi(i),0,'bx')
%                 end
%             end
%             
%             figure(pat_contrastfig);
%             h = subplot(4,3,this_pat);
%             hold on
%             stdshade_TEC(squeeze(all_mismatch_contrasts(from,to,:,:,s))',0.2,'g',foi,1,1)
%             
%             stdshade_TEC(squeeze(all_mismatch_contrasts(to,from,:,:,s))',0.2,'b',foi,1,1)
%             plot(foi,squeeze(mean(all_random_mismatch_contrasts(from,to,:,:,:,s),6)),'k')
%             plot(foi,squeeze(mean(all_random_mismatch_contrasts(to,from,:,:,:,s),6)),'k')
%             stdshade_TEC(squeeze(all_mismatch_contrasts(from,to,:,:,s))',0.2,'g',foi,1,1)
%             
%             stdshade_TEC(squeeze(all_mismatch_contrasts(to,from,:,:,s))',0.2,'b',foi,1,1)
%             
%             %legend({[from_name '-' to_name],[to_name '-' from_name]})
%             for i = 1:topfreqband %Compare across group permutation
%                 p_tf(i) = 1-relRankIn_includeValue_lowerBound(squeeze(mean(all_random_mismatch_contrasts(from,to,:,i,:,s),6)),squeeze(mean(all_mismatch_contrasts(from,to,:,i,s),5)));
%                 p_ft(i) = 1-relRankIn_includeValue_lowerBound(squeeze(mean(all_random_mismatch_contrasts(to,from,:,i,:,s),6)),squeeze(mean(all_mismatch_contrasts(to,from,:,i,s),5)));
%                 if p_tf(i)<=0.05 || p_tf(i)>=0.95
%                     plot(foi(i),0,'g*')
%                 end
%                 if p_ft(i)<=0.05 || p_ft(i)>=0.95
%                     plot(foi(i),0,'bx')
%                 end
%             end
%             
%             figure(pat_interactionfig);
%             h = subplot(4,3,this_pat);
%             hold on
%             stdshade_TEC(squeeze(all_interaction_contrasts(from,to,:,:,s))',0.2,'g',foi,1,1)
%             
%             stdshade_TEC(squeeze(all_interaction_contrasts(to,from,:,:,s))',0.2,'b',foi,1,1)
%             plot(foi,squeeze(mean(all_random_interaction_contrasts(from,to,:,:,:,s),6)),'k')
%             plot(foi,squeeze(mean(all_random_interaction_contrasts(to,from,:,:,:,s),6)),'k')
%             stdshade_TEC(squeeze(all_interaction_contrasts(from,to,:,:,s))',0.2,'g',foi,1,1)
%             
%             stdshade_TEC(squeeze(all_interaction_contrasts(to,from,:,:,s))',0.2,'b',foi,1,1)
%             
%             %legend({[from_name '-' to_name],[to_name '-' from_name]})
%             for i = 1:topfreqband %Compare across group permutation
%                 p_tf(i) = 1-relRankIn_includeValue_lowerBound(squeeze(mean(all_random_interaction_contrasts(from,to,:,i,:,s),6)),squeeze(mean(all_interaction_contrasts(from,to,:,i,s),5)));
%                 p_ft(i) = 1-relRankIn_includeValue_lowerBound(squeeze(mean(all_random_interaction_contrasts(to,from,:,i,:,s),6)),squeeze(mean(all_interaction_contrasts(to,from,:,i,s),5)));
%                 if p_tf(i)<=0.05 || p_tf(i)>=0.95
%                     plot(foi(i),0,'g*')
%                 end
%                 if p_ft(i)<=0.05 || p_ft(i)>=0.95
%                     plot(foi(i),0,'bx')
%                 end
%             end
%             
%         end
%     end
%     figure(config)
%     title(['Controls ' analysis_type ' time ' num2str(start_times) '-' num2str(end_times)])
%     figure(patfig)
%     title(['Patients ' analysis_type ' time ' num2str(start_times) '-' num2str(end_times)])
%     figure(con_contrastfig);
%     title(['Controls Match-Mismatch ' analysis_type])
%     figure(pat_contrastfig);
%     title(['Patients Match-Mismatch ' analysis_type])
%     figure(con_interactionfig);
%     title(['Controls Interaction ' analysis_type])
%     figure(pat_interactionfig);
%     title(['Patients Interaction ' analysis_type])
    
    
    if closeafter == 1; close all; end; figure
    hold on
    linehandle(1) = stdshade_TEC(squeeze(all_mismatch_contrasts(from,to,:,:,:))',0.2,'g',foi,1,1);
    
    linehandle(2) = stdshade_TEC(squeeze(all_mismatch_contrasts(to,from,:,:,:))',0.2,'b',foi,1,1);
    plot(foi,squeeze(mean(all_random_mismatch_contrasts(from,to,:,:,:,:),6)),'k')
    plot(foi,squeeze(mean(all_random_mismatch_contrasts(to,from,:,:,:,:),6)),'k')
%     stdshade_TEC(squeeze(all_mismatch_contrasts(from,to,:,:,:))',0.2,'g',foi,1,1)
%     
%     stdshade_TEC(squeeze(all_mismatch_contrasts(to,from,:,:,:))',0.2,'b',foi,1,1)
    title(['All Standard-Deviant ' analysis_type ' time ' num2str(start_times) '-' num2str(end_times)])
    legend(linehandle(1:2),{[from_name '-' to_name],[to_name '-' from_name]})
    for i = 1:topfreqband %Compare across group permutation
        p_tf(i) = 1-relRankIn_includeValue_lowerBound(squeeze(mean(all_random_mismatch_contrasts(from,to,:,i,:,:),6)),squeeze(mean(all_mismatch_contrasts(from,to,:,i,:),5)));
        p_ft(i) = 1-relRankIn_includeValue_lowerBound(squeeze(mean(all_random_mismatch_contrasts(to,from,:,i,:,:),6)),squeeze(mean(all_mismatch_contrasts(to,from,:,i,:),5)));
        if p_tf(i)<=0.05 || p_tf(i)>=0.95
            plot(foi(i),0,'g*')
        end
        if p_ft(i)<=0.05 || p_ft(i)>=0.95
            plot(foi(i),0,'bx')
        end
    end
    savestring = ['./figures/' from_name '_' to_name '_All_mean_std-dvt_' analysis_type '_time_' num2str(start_times) '-' num2str(end_times) '.pdf'];
    savestring = strrep(savestring,' ','_');
    eval(['export_fig ' savestring ' -transparent']); eval(['export_fig ' savestring(1:end-3) 'png -transparent']);
    
    if closeafter == 1; close all; end; figure
    hold on
    linehandle(1) = stdshade_TEC_median(squeeze(all_mismatch_contrasts(from,to,:,:,:))',0.2,'g',foi,1,1);
    
    linehandle(2) = stdshade_TEC_median(squeeze(all_mismatch_contrasts(to,from,:,:,:))',0.2,'b',foi,1,1);
    plot(foi,squeeze(median(all_random_mismatch_contrasts(from,to,:,:,:,:),6)),'k')
    plot(foi,squeeze(median(all_random_mismatch_contrasts(to,from,:,:,:,:),6)),'k')
%     stdshade_TEC(squeeze(all_mismatch_contrasts(from,to,:,:,:))',0.2,'g',foi,1,1)
%     
%     stdshade_TEC(squeeze(all_mismatch_contrasts(to,from,:,:,:))',0.2,'b',foi,1,1)
    title(['All median Standard-Deviant ' analysis_type ' time ' num2str(start_times) '-' num2str(end_times)])
    legend(linehandle(1:2), {[from_name '-' to_name],[to_name '-' from_name]})
    for i = 1:topfreqband %Compare across group permutation
        p_tf(i) = 1-relRankIn_includeValue_lowerBound(squeeze(median(all_random_mismatch_contrasts(from,to,:,i,:,:),6)),squeeze(median(all_mismatch_contrasts(from,to,:,i,:),5)));
        p_ft(i) = 1-relRankIn_includeValue_lowerBound(squeeze(median(all_random_mismatch_contrasts(to,from,:,i,:,:),6)),squeeze(median(all_mismatch_contrasts(to,from,:,i,:),5)));
        if p_tf(i)<=0.05 || p_tf(i)>=0.95
            plot(foi(i),0,'g*')
        end
        if p_ft(i)<=0.05 || p_ft(i)>=0.95
            plot(foi(i),0,'bx')
        end
    end
    savestring = ['./figures/' from_name '_' to_name '_All_median_std-dvt_' analysis_type '_time_' num2str(start_times) '-' num2str(end_times) '.pdf'];
    savestring = strrep(savestring,' ','_');
    eval(['export_fig ' savestring ' -transparent']); eval(['export_fig ' savestring(1:end-3) 'png -transparent']);
    
%     if closeafter == 1; close all; end; figure
%     hold on
%     stdshade_TEC(squeeze(all_clarity_contrasts(from,to,:,:,:))',0.2,'g',foi,1,1)
%     
%     stdshade_TEC(squeeze(all_clarity_contrasts(to,from,:,:,:))',0.2,'b',foi,1,1)
%     plot(foi,squeeze(mean(all_random_clarity_contrasts(from,to,:,:,:,:),6)),'k')
%     plot(foi,squeeze(mean(all_random_clarity_contrasts(to,from,:,:,:,:),6)),'k')
%     stdshade_TEC(squeeze(all_clarity_contrasts(from,to,:,:,:))',0.2,'g',foi,1,1)
%     
%     stdshade_TEC(squeeze(all_clarity_contrasts(to,from,:,:,:))',0.2,'b',foi,1,1)
%     title(['All Clear-Unclear ' analysis_type ' time ' num2str(start_times) '-' num2str(end_times)])
%     legend({[from_name '-' to_name],[to_name '-' from_name]})
%     for i = 1:topfreqband %Compare across group permutation
%         p_tf(i) = 1-relRankIn_includeValue_lowerBound(squeeze(mean(all_random_clarity_contrasts(from,to,:,i,:,:),6)),squeeze(mean(all_clarity_contrasts(from,to,:,i,:),5)));
%         p_ft(i) = 1-relRankIn_includeValue_lowerBound(squeeze(mean(all_random_clarity_contrasts(to,from,:,i,:,:),6)),squeeze(mean(all_clarity_contrasts(to,from,:,i,:),5)));
%         if p_tf(i)<=0.05 || p_tf(i)>=0.95
%             plot(foi(i),0,'g*')
%         end
%         if p_ft(i)<=0.05 || p_ft(i)>=0.95
%             plot(foi(i),0,'bx')
%         end
%     end
%     
%     if closeafter == 1; close all; end; figure
%     hold on
%     stdshade_TEC(squeeze(all_clarity_contrasts(from,to,:,:,:))',0.2,'g',foi,1,1)
%     
%     stdshade_TEC(squeeze(all_clarity_contrasts(to,from,:,:,:))',0.2,'b',foi,1,1)
%     plot(foi,squeeze(median(all_random_clarity_contrasts(from,to,:,:,:,:),6)),'k')
%     plot(foi,squeeze(median(all_random_clarity_contrasts(to,from,:,:,:,:),6)),'k')
%     stdshade_TEC(squeeze(all_clarity_contrasts(from,to,:,:,:))',0.2,'g',foi,1,1)
%     
%     stdshade_TEC(squeeze(all_clarity_contrasts(to,from,:,:,:))',0.2,'b',foi,1,1)
%     title(['All median Clear-Unclear ' analysis_type ' time ' num2str(start_times) '-' num2str(end_times)])
%     legend({[from_name '-' to_name],[to_name '-' from_name]})
%     for i = 1:topfreqband %Compare across group permutation
%         p_tf(i) = 1-relRankIn_includeValue_lowerBound(squeeze(median(all_random_clarity_contrasts(from,to,:,i,:,:),6)),squeeze(median(all_clarity_contrasts(from,to,:,i,:),5)));
%         p_ft(i) = 1-relRankIn_includeValue_lowerBound(squeeze(median(all_random_clarity_contrasts(to,from,:,i,:,:),6)),squeeze(median(all_clarity_contrasts(to,from,:,i,:),5)));
%         if p_tf(i)<=0.05 || p_tf(i)>=0.95
%             plot(foi(i),0,'g*')
%         end
%         if p_ft(i)<=0.05 || p_ft(i)>=0.95
%             plot(foi(i),0,'bx')
%         end
%     end
%     
    if closeafter == 1; close all; end; figure
    hold on
    legend_text = cell(1,2*length(groupstodo));
    for grp = 1:2:2*length(groupstodo)
        linehandle(grp) = stdshade_TEC_cmap(squeeze(all_mismatch_contrasts(from,to,:,:,group==ceil(grp/2)))',0.2,cmap(ceil(grp/2),:),foi,1,1,':');
        linehandle(grp+1) = stdshade_TEC_cmap(squeeze(all_mismatch_contrasts(to,from,:,:,group==ceil(grp/2)))',0.2,cmap(ceil(grp/2),:),foi,1,1,'--');
        legend_text{grp} = [from_name '-' to_name ' ' groupstodo{ceil(grp/2)}]; 
        legend_text{grp+1} = [to_name '-' from_name ' ' groupstodo{ceil(grp/2)}];
    end
    plot(foi,zeros(1,length(foi)),'k--');
    legend(linehandle,legend_text)
    title(['By group Mean Standard-Deviant ' analysis_type ])
    savestring = ['./figures/' from_name '_' to_name '_Group_mean_std-dvt_' analysis_type '_time_' num2str(start_times) '-' num2str(end_times) '.pdf'];
    savestring = strrep(savestring,' ','_');
    eval(['export_fig ' savestring ' -transparent']); eval(['export_fig ' savestring(1:end-3) 'png -transparent']);

%     for i = 1:topfreqband %Rough initial parametric stats - better to use permutation
%         [h_tf(i) p_tf(i)] = ttest2(squeeze(all_patients_mismatch_contrasts(from,to,:,i,:)),squeeze(all_controls_mismatch_contrasts(from,to,:,i,:)));
%         [p_tf(i) h_tf(i)] = ttest2(squeeze(all_absolute_patients_mismatch_contrasts(from,to,:,i,:)),squeeze(all_absolute_controls_mismatch_contrasts(from,to,:,i,:)));
%         [h_ft(i) p_ft(i)] = ttest2(squeeze(all_patients_mismatch_contrasts(to,from,:,i,:)),squeeze(all_controls_mismatch_contrasts(to,from,:,i,:)));
%         if h_tf(i) == 1
%             plot(foi(i),0,'g*')
%         end
%         if h_ft(i) == 1
%             plot(foi(i),0,'bx')
%         end
%     end
%     
%     if closeafter == 1; close all; end; figure
%     hold on
%     stdshade_TEC(squeeze(all_controls_mismatch_contrasts(from,to,:,:,:))',0.2,'g:',foi,1,1)
%     stdshade_TEC(squeeze(all_patients_mismatch_contrasts(from,to,:,:,:))',0.2,'g--',foi,1,1)
%     stdshade_TEC(squeeze(all_controls_mismatch_contrasts(to,from,:,:,:))',0.2,'b:',foi,1,1)
%     stdshade_TEC(squeeze(all_patients_mismatch_contrasts(to,from,:,:,:))',0.2,'b--',foi,1,1)
%     title(['By group Median Standard-Deviant ' analysis_type ])
%     legend({[from_name '-' to_name ' control'],[from_name '-' to_name ' patient'],[to_name '-' from_name ' control'],[to_name '-' from_name ' patient']})
%     for i = 1:topfreqband %Rough initial parametric stats - better to use permutation
%         [p_tf(i) h_tf(i)] = ranksum(squeeze(all_patients_mismatch_contrasts(from,to,:,i,:)),squeeze(all_controls_mismatch_contrasts(from,to,:,i,:)));
%         [p_ft(i) h_ft(i)] = ranksum(squeeze(all_patients_mismatch_contrasts(to,from,:,i,:)),squeeze(all_controls_mismatch_contrasts(to,from,:,i,:)));
%         if h_tf(i) == 1
%             plot(foi(i),0,'g*')
%         end
%         if h_ft(i) == 1
%             plot(foi(i),0,'bx')
%         end
%     end
%     
%     if closeafter == 1; close all; end; figure
%     stdshade_TEC(squeeze(all_interaction_contrasts(from,to,:,:,:))',0.2,'g',foi,1,1)
%     hold on
%     stdshade_TEC(squeeze(all_interaction_contrasts(to,from,:,:,:))',0.2,'b',foi,1,1)
%     plot(foi,squeeze(mean(all_random_interaction_contrasts(from,to,:,:,:,:),6)),'k')
%     plot(foi,squeeze(mean(all_random_interaction_contrasts(to,from,:,:,:,:),6)),'k')
%     stdshade_TEC(squeeze(all_interaction_contrasts(from,to,:,:,:))',0.2,'g',foi,1,1)
%     stdshade_TEC(squeeze(all_interaction_contrasts(to,from,:,:,:))',0.2,'b',foi,1,1)
%     title(['All Interaction ' analysis_type ' time ' num2str(start_times) '-' num2str(end_times)])
%     legend({[from_name '-' to_name],[to_name '-' from_name]})
%     for i = 1:topfreqband %Compare across group permutation
%         p_tf(i) = 1-relRankIn_includeValue_lowerBound(squeeze(mean(all_random_interaction_contrasts(from,to,:,i,:,:),6)),squeeze(mean(all_interaction_contrasts(from,to,:,i,:),5)));
%         p_ft(i) = 1-relRankIn_includeValue_lowerBound(squeeze(mean(all_random_interaction_contrasts(to,from,:,i,:,:),6)),squeeze(mean(all_interaction_contrasts(to,from,:,i,:),5)));
%         if p_tf(i)<=0.05 || p_tf(i)>=0.95
%             plot(foi(i),0,'g*')
%         end
%         if p_ft(i)<=0.05 || p_ft(i)>=0.95
%             plot(foi(i),0,'bx')
%         end
%     end
%     
%     if closeafter == 1; close all; end; figure
%     stdshade_TEC(squeeze(all_interaction_contrasts(from,to,:,:,:))',0.2,'g',foi,1,1)
%     hold on
%     stdshade_TEC(squeeze(all_interaction_contrasts(to,from,:,:,:))',0.2,'b',foi,1,1)
%     plot(foi,squeeze(median(all_random_interaction_contrasts(from,to,:,:,:,:),6)),'k')
%     plot(foi,squeeze(median(all_random_interaction_contrasts(to,from,:,:,:,:),6)),'k')
%     stdshade_TEC(squeeze(all_interaction_contrasts(from,to,:,:,:))',0.2,'g',foi,1,1)
%     
%     stdshade_TEC(squeeze(all_interaction_contrasts(to,from,:,:,:))',0.2,'b',foi,1,1)
%     title(['All median Interaction ' analysis_type ' time ' num2str(start_times) '-' num2str(end_times)])
%     legend({[from_name '-' to_name],[to_name '-' from_name]})
%     for i = 1:topfreqband %Compare across group permutation
%         p_tf(i) = 1-relRankIn_includeValue_lowerBound(squeeze(median(all_random_interaction_contrasts(from,to,:,i,:,:),6)),squeeze(median(all_interaction_contrasts(from,to,:,i,:),5)));
%         p_ft(i) = 1-relRankIn_includeValue_lowerBound(squeeze(median(all_random_interaction_contrasts(to,from,:,i,:,:),6)),squeeze(median(all_interaction_contrasts(to,from,:,i,:),5)));
%         if p_tf(i)<=0.05 || p_tf(i)>=0.95
%             plot(foi(i),0,'g*')
%         end
%         if p_ft(i)<=0.05 || p_ft(i)>=0.95
%             plot(foi(i),0,'bx')
%         end
%     end
%     
%     if closeafter == 1; close all; end; figure
%     hold on
%     stdshade_TEC(squeeze(all_interaction_contrasts(to,from,:,:,:))',0.2,'b',foi,1,1)
%     stdshade_TEC(squeeze(all_interaction_contrasts(to,from,:,:,:))',0.2,'b',foi,1,1)
%     stdshade_TEC(squeeze(all_interaction_contrasts(to,from,:,:,:))',0.2,'b',foi,1,1)
%     stdshade_TEC(squeeze(all_interaction_contrasts(to,from,:,:,:))',0.2,'b',foi,1,1)
%     plot(foi,squeeze(mean(all_random_controls_mismatch_contrasts(from,to,:,:,:,:),6)),'k:')
%     plot(foi,squeeze(mean(all_random_patients_mismatch_contrasts(from,to,:,:,:,:),6)),'k--')
%     plot(foi,squeeze(mean(all_random_controls_mismatch_contrasts(to,from,:,:,:,:),6)),'k:')
%     plot(foi,squeeze(mean(all_random_patients_mismatch_contrasts(to,from,:,:,:,:),6)),'k--')
%     stdshade_TEC(squeeze(all_interaction_contrasts(to,from,:,:,:))',0.2,'b',foi,1,1)
%     stdshade_TEC(squeeze(all_interaction_contrasts(to,from,:,:,:))',0.2,'b',foi,1,1)
%     stdshade_TEC(squeeze(all_interaction_contrasts(to,from,:,:,:))',0.2,'b',foi,1,1)
%     stdshade_TEC(squeeze(all_interaction_contrasts(to,from,:,:,:))',0.2,'b',foi,1,1)
%     title(['By group sig Mean Standard-Deviant ' analysis_type ])
%     legend({[from_name '-' to_name ' control'],[from_name '-' to_name ' patient'],[to_name '-' from_name ' control'],[to_name '-' from_name ' patient']})
%     for i = 1:topfreqband %Rough initial parametric stats - better to use permutation
%         p_tf_c(i) = 1-relRankIn_includeValue_lowerBound(squeeze(mean(all_random_controls_mismatch_contrasts(from,to,i,:,:,:),6)),squeeze(mean(all_controls_mismatch_contrasts(from,to,:,i,:),5)));
%         p_ft_c(i) = 1-relRankIn_includeValue_lowerBound(squeeze(mean(all_random_controls_mismatch_contrasts(to,from,i,:,:,:),6)),squeeze(mean(all_controls_mismatch_contrasts(to,from,:,i,:),5)));
%         if p_tf_c(i)<=0.05 || p_ft_c(i)>=0.95
%             plot(foi(i),0,'go')
%         end
%         if p_ft_c(i)<=0.05 || p_ft_c(i)>=0.95
%             plot(foi(i),0,'b.')
%         end
%         p_tf_p(i) = 1-relRankIn_includeValue_lowerBound(squeeze(mean(all_random_patients_mismatch_contrasts(from,to,i,:,:,:),6)),squeeze(mean(all_patients_mismatch_contrasts(from,to,:,i,:),5)));
%         p_ft_p(i) = 1-relRankIn_includeValue_lowerBound(squeeze(mean(all_random_patients_mismatch_contrasts(to,from,i,:,:,:),6)),squeeze(mean(all_patients_mismatch_contrasts(to,from,:,i,:),5)));
%         if p_tf_p(i)<=0.05 || p_ft_p(i)>=0.95
%             plot(foi(i),0,'gp')
%         end
%         if p_ft_p(i)<=0.05 || p_ft_p(i)>=0.95
%             plot(foi(i),0,'bs')
%         end
%     end
%     
%     if closeafter == 1; close all; end; figure
%     hold on
%     stdshade_TEC(squeeze(all_interaction_contrasts(to,from,:,:,:))',0.2,'b',foi,1,1)
%     stdshade_TEC(squeeze(all_interaction_contrasts(to,from,:,:,:))',0.2,'b',foi,1,1)
%     stdshade_TEC(squeeze(all_interaction_contrasts(to,from,:,:,:))',0.2,'b',foi,1,1)
%     stdshade_TEC(squeeze(all_interaction_contrasts(to,from,:,:,:))',0.2,'b',foi,1,1)
%     plot(foi,squeeze(median(all_random_controls_mismatch_contrasts(from,to,:,:,:,:),6)),'k:')
%     plot(foi,squeeze(median(all_random_patients_mismatch_contrasts(from,to,:,:,:,:),6)),'k--')
%     plot(foi,squeeze(median(all_random_controls_mismatch_contrasts(to,from,:,:,:,:),6)),'k:')
%     plot(foi,squeeze(median(all_random_patients_mismatch_contrasts(to,from,:,:,:,:),6)),'k--')
%     stdshade_TEC(squeeze(all_interaction_contrasts(to,from,:,:,:))',0.2,'b',foi,1,1)
%     stdshade_TEC(squeeze(all_interaction_contrasts(to,from,:,:,:))',0.2,'b',foi,1,1)
%     stdshade_TEC(squeeze(all_interaction_contrasts(to,from,:,:,:))',0.2,'b',foi,1,1)
%     stdshade_TEC(squeeze(all_interaction_contrasts(to,from,:,:,:))',0.2,'b',foi,1,1)
%     title(['By group sig Median Standard-Deviant ' analysis_type ])
%     legend({[from_name '-' to_name ' control'],[from_name '-' to_name ' patient'],[to_name '-' from_name ' control'],[to_name '-' from_name ' patient']})
%     for i = 1:topfreqband %Rough initial parametric stats - better to use permutation
%         p_tf_c(i) = 1-relRankIn_includeValue_lowerBound(squeeze(median(all_random_controls_mismatch_contrasts(from,to,i,:,:,:),6)),squeeze(median(all_controls_mismatch_contrasts(from,to,:,i,:),5)));
%         p_ft_c(i) = 1-relRankIn_includeValue_lowerBound(squeeze(median(all_random_controls_mismatch_contrasts(to,from,i,:,:,:),6)),squeeze(median(all_controls_mismatch_contrasts(to,from,:,i,:),5)));
%         if p_tf_c(i)<=0.05 || p_ft_c(i)>=0.95
%             plot(foi(i),0,'go')
%         end
%         if p_ft_c(i)<=0.05 || p_ft_c(i)>=0.95
%             plot(foi(i),0,'b.')
%         end
%         p_tf_p(i) = 1-relRankIn_includeValue_lowerBound(squeeze(median(all_random_patients_mismatch_contrasts(from,to,i,:,:,:),6)),squeeze(median(all_patients_mismatch_contrasts(from,to,:,i,:),5)));
%         p_ft_p(i) = 1-relRankIn_includeValue_lowerBound(squeeze(median(all_random_patients_mismatch_contrasts(to,from,i,:,:,:),6)),squeeze(median(all_patients_mismatch_contrasts(to,from,:,i,:),5)));
%         if p_tf_p(i)<=0.05 || p_ft_p(i)>=0.95
%             plot(foi(i),0,'gp')
%         end
%         if p_ft_p(i)<=0.05 || p_ft_p(i)>=0.95
%             plot(foi(i),0,'bs')
%         end
%     end
%     
%     if closeafter == 1; close all; end; figure
%     hold on
%     stdshade_TEC(squeeze(all_interaction_contrasts(to,from,:,:,:))',0.2,'b',foi,1,1)
%     stdshade_TEC(squeeze(all_interaction_contrasts(to,from,:,:,:))',0.2,'b',foi,1,1)
%     stdshade_TEC(squeeze(all_interaction_contrasts(to,from,:,:,:))',0.2,'b',foi,1,1)
%     stdshade_TEC(squeeze(all_interaction_contrasts(to,from,:,:,:))',0.2,'b',foi,1,1)
%     plot(foi,squeeze(mean(all_random_controls_interaction_contrasts(from,to,:,:,:,:),6)),'k:')
%     plot(foi,squeeze(mean(all_random_patients_interaction_contrasts(from,to,:,:,:,:),6)),'k--')
%     plot(foi,squeeze(mean(all_random_controls_interaction_contrasts(to,from,:,:,:,:),6)),'k:')
%     plot(foi,squeeze(mean(all_random_patients_interaction_contrasts(to,from,:,:,:,:),6)),'k--')
%     stdshade_TEC(squeeze(all_interaction_contrasts(to,from,:,:,:))',0.2,'b',foi,1,1)
%     stdshade_TEC(squeeze(all_interaction_contrasts(to,from,:,:,:))',0.2,'b',foi,1,1)
%     stdshade_TEC(squeeze(all_interaction_contrasts(to,from,:,:,:))',0.2,'b',foi,1,1)
%     stdshade_TEC(squeeze(all_interaction_contrasts(to,from,:,:,:))',0.2,'b',foi,1,1)
%     title(['By group Interaction ' analysis_type ])
%     legend({[from_name '-' to_name ' control'],[from_name '-' to_name ' patient'],[to_name '-' from_name ' control'],[to_name '-' from_name ' patient']})
%     for i = 1:topfreqband %Rough initial parametric stats - better to use permutation
%         p_tf_c(i) = 1-relRankIn_includeValue_lowerBound(squeeze(mean(all_random_controls_interaction_contrasts(from,to,i,:,:,:),6)),squeeze(mean(all_controls_interaction_contrasts(from,to,:,i,:),5)));
%         p_ft_c(i) = 1-relRankIn_includeValue_lowerBound(squeeze(mean(all_random_controls_interaction_contrasts(to,from,i,:,:,:),6)),squeeze(mean(all_controls_interaction_contrasts(to,from,:,i,:),5)));
%         if p_tf_c(i)<=0.05 || p_tf_c(i)>=0.95
%             plot(foi(i),0,'go')
%         end
%         if p_ft_c(i)<=0.05 || p_ft_c(i)>=0.95
%             plot(foi(i),0,'b.')
%         end
%         p_tf_p(i) = 1-relRankIn_includeValue_lowerBound(squeeze(mean(all_random_patients_interaction_contrasts(from,to,i,:,:,:),6)),squeeze(mean(all_patients_interaction_contrasts(from,to,:,i,:),5)));
%         p_ft_p(i) = 1-relRankIn_includeValue_lowerBound(squeeze(mean(all_random_patients_interaction_contrasts(to,from,i,:,:,:),6)),squeeze(mean(all_patients_interaction_contrasts(to,from,:,i,:),5)));
%         if p_tf_p(i)<=0.05 || p_ft_p(i)>=0.95
%             plot(foi(i),0,'gp')
%         end
%         if p_ft_p(i)<=0.05 || p_ft_p(i)>=0.95
%             plot(foi(i),0,'bs')
%         end
%     end
%     
%     if closeafter == 1; close all; end; figure
%     hold on
%     stdshade_TEC(squeeze(all_interaction_contrasts(to,from,:,:,:))',0.2,'b',foi,1,1)
%     stdshade_TEC(squeeze(all_interaction_contrasts(to,from,:,:,:))',0.2,'b',foi,1,1)
%     stdshade_TEC(squeeze(all_interaction_contrasts(to,from,:,:,:))',0.2,'b',foi,1,1)
%     stdshade_TEC(squeeze(all_interaction_contrasts(to,from,:,:,:))',0.2,'b',foi,1,1)
%     plot(foi,squeeze(mean(all_random_controls_clarity_contrasts(from,to,:,:,:,:),6)),'k:')
%     plot(foi,squeeze(mean(all_random_patients_clarity_contrasts(from,to,:,:,:,:),6)),'k--')
%     plot(foi,squeeze(mean(all_random_controls_clarity_contrasts(to,from,:,:,:,:),6)),'k:')
%     plot(foi,squeeze(mean(all_random_patients_clarity_contrasts(to,from,:,:,:,:),6)),'k--')
%     stdshade_TEC(squeeze(all_interaction_contrasts(to,from,:,:,:))',0.2,'b',foi,1,1)
%     stdshade_TEC(squeeze(all_interaction_contrasts(to,from,:,:,:))',0.2,'b',foi,1,1)
%     stdshade_TEC(squeeze(all_interaction_contrasts(to,from,:,:,:))',0.2,'b',foi,1,1)
%     stdshade_TEC(squeeze(all_interaction_contrasts(to,from,:,:,:))',0.2,'b',foi,1,1)
%     title(['By group Clarity ' analysis_type ])
%     legend({[from_name '-' to_name ' control'],[from_name '-' to_name ' patient'],[to_name '-' from_name ' control'],[to_name '-' from_name ' patient']})
%     for i = 1:topfreqband %Rough initial parametric stats - better to use permutation
%         p_tf_c(i) = 1-relRankIn_includeValue_lowerBound(squeeze(mean(all_random_controls_clarity_contrasts(from,to,i,:,:,:),6)),squeeze(mean(all_controls_clarity_contrasts(from,to,:,i,:),5)));
%         p_ft_c(i) = 1-relRankIn_includeValue_lowerBound(squeeze(mean(all_random_controls_clarity_contrasts(to,from,i,:,:,:),6)),squeeze(mean(all_controls_clarity_contrasts(to,from,:,i,:),5)));
%         if p_tf_c(i)<=0.05 || p_tf_c(i)>=0.95
%             plot(foi(i),0,'go')
%         end
%         if p_ft_c(i)<=0.05 || p_ft_c(i)>=0.95
%             plot(foi(i),0,'b.')
%         end
%         p_tf_p(i) = 1-relRankIn_includeValue_lowerBound(squeeze(mean(all_random_patients_clarity_contrasts(from,to,i,:,:,:),6)),squeeze(mean(all_patients_clarity_contrasts(from,to,:,i,:),5)));
%         p_ft_p(i) = 1-relRankIn_includeValue_lowerBound(squeeze(mean(all_random_patients_clarity_contrasts(to,from,i,:,:,:),6)),squeeze(mean(all_patients_clarity_contrasts(to,from,:,i,:),5)));
%         if p_tf_p(i)<=0.05 || p_ft_p(i)>=0.95
%             plot(foi(i),0,'gp')
%         end
%         if p_ft_p(i)<=0.05 || p_ft_p(i)>=0.95
%             plot(foi(i),0,'bs')
%         end
%     end
%     %
%     
%     %         case 'icoh'
%     %
%     %             if closeafter == 1; close all; end; figure
%     %             stdshade_TEC(squeeze(mean(abs(all_granger_data(from,to,:,:,:,:)),5))',0.2,'g',foi,1,1)
%     %             hold on
%     %             stdshade_TEC(squeeze(mean(abs(all_granger_data(to,from,:,:,:,:)),5))',0.2,'b',foi,1,1)
%     %             plot(foi,squeeze(mean(mean(abs(all_random_granger_data(from,to,:,:,:,:,:)),5),7)),'k')
%     %             plot(foi,squeeze(mean(mean(abs(all_random_granger_data(to,from,:,:,:,:,:)),5),7)),'k')
%     %             stdshade_TEC(squeeze(mean(abs(all_granger_data(from,to,:,:,:,:)),5))',0.2,'g',foi,1,1)
%     %             stdshade_TEC(squeeze(mean(abs(all_granger_data(to,from,:,:,:,:)),5))',0.2,'b',foi,1,1)
%     %
%     %             title(['All abs ' analysis_type ' time ' num2str(start_times) '-' num2str(end_times)])
%     %             legend({[from_name '-' to_name],[to_name '-' from_name]})
%     %             for i = 1:topfreqband %Compare across group permutation
%     %                 p_tf(i) = 1-relRankIn_includeValue_lowerBound(squeeze(mean(mean(abs(all_random_granger_data(from,to,:,i,:,:,:)),5),7)),squeeze(mean(mean(abs(all_granger_data(from,to,:,i,:,:)),5),6)))
%     %                 p_ft(i) = 1-relRankIn_includeValue_lowerBound(squeeze(mean(mean(abs(all_random_granger_data(to,from,:,i,:,:,:)),5),7)),squeeze(mean(mean(abs(all_granger_data(to,from,:,i,:,:)),5),6)))
%     %                 if p_tf(i)<=0.05 || p_tf(i)>=0.95
%     %                     plot(foi(i),0,'g*')
%     %                 end
%     %                 if p_ft(i)<=0.05 || p_ft(i)>=0.95
%     %                     plot(foi(i),0,'bx')
%     %                 end
%     %             end
%     %
%     %             if closeafter == 1; close all; end; figure
%     %             stdshade_TEC(squeeze(median(abs(all_granger_data(from,to,:,:,:,:)),5))',0.2,'g',foi,1,1)
%     %             hold on
%     %             stdshade_TEC(squeeze(median(abs(all_granger_data(to,from,:,:,:,:)),5))',0.2,'b',foi,1,1)
%     %             plot(foi,squeeze(median(mean(abs(all_random_granger_data(from,to,:,:,:,:,:)),5),7)),'k')
%     %             plot(foi,squeeze(median(mean(abs(all_random_granger_data(to,from,:,:,:,:,:)),5),7)),'k')
%     %             stdshade_TEC(squeeze(median(abs(all_granger_data(from,to,:,:,:,:)),5))',0.2,'g',foi,1,1)
%     %             stdshade_TEC(squeeze(median(abs(all_granger_data(to,from,:,:,:,:)),5))',0.2,'b',foi,1,1)
%     %
%     %             title(['All median ' analysis_type ' time ' num2str(start_times) '-' num2str(end_times)])
%     %             legend({[from_name '-' to_name],[to_name '-' from_name]})
%     %             for i = 1:topfreqband %Compare across group permutation
%     %                 p_tf(i) = 1-relRankIn_includeValue_lowerBound(squeeze(median(mean(abs(all_random_granger_data(from,to,:,i,:,:,:)),5),7)),squeeze(median(mean(abs(all_granger_data(from,to,:,i,:,:)),5),6)))
%     %                 p_ft(i) = 1-relRankIn_includeValue_lowerBound(squeeze(median(mean(abs(all_random_granger_data(to,from,:,i,:,:,:)),5),7)),squeeze(median(mean(abs(all_granger_data(to,from,:,i,:,:)),5),6)))
%     %                 if p_tf(i)<=0.05 || p_tf(i)>=0.95
%     %                     plot(foi(i),0,'g*')
%     %                 end
%     %                 if p_ft(i)<=0.05 || p_ft(i)>=0.95
%     %                     plot(foi(i),0,'bx')
%     %                 end
%     %             end
%     %
%     %             if closeafter == 1; close all; end; figure
%     %             hold on
%     %             stdshade_TEC(squeeze(mean(abs(all_controls_granger_data(from,to,:,:,:,:)),5))',0.2,'g:',foi,1,1)
%     %             stdshade_TEC(squeeze(mean(abs(all_patients_granger_data(from,to,:,:,:,:)),5))',0.2,'g--',foi,1,1)
%     %             stdshade_TEC(squeeze(mean(abs(all_controls_granger_data(to,from,:,:,:,:)),5))',0.2,'b:',foi,1,1)
%     %             stdshade_TEC(squeeze(mean(abs(all_patients_granger_data(to,from,:,:,:,:)),5))',0.2,'b--',foi,1,1)
%     %             legend({[from_name '-' to_name ' control'],[from_name '-' to_name ' patient'],[to_name '-' from_name ' control'],[to_name '-' from_name ' patient']})
%     %             title(['By group ' analysis_type ' time ' num2str(start_times) '-' num2str(end_times)])
%     %             for i = 1:topfreqband %Rough initial parametric stats - better to use permutation
%     %                 [h_tf(i) p_tf(i)] = ttest2(mean(abs(all_controls_granger_data(from,to,:,i,:,:)),5),mean(abs(all_patients_granger_data(from,to,:,i,:,:)),5));
%     %                 [h_ft(i) p_ft(i)] = ttest2(mean(abs(all_controls_granger_data(to,from,:,i,:,:)),5),mean(abs(all_patients_granger_data(to,from,:,i,:,:)),5));
%     %                 if h_tf(i) == 1
%     %                     plot(foi(i),0,'g*')
%     %                 end
%     %                 if h_ft(i) == 1
%     %                     plot(foi(i),0,'bx')
%     %                 end
%     %             end
%     %
%     %             if closeafter == 1; close all; end; figure
%     %             hold on
%     %             stdshade_TEC(squeeze(median(abs(all_controls_granger_data(from,to,:,:,:,:)),5))',0.2,'g:',foi,1,1)
%     %             stdshade_TEC(squeeze(median(abs(all_patients_granger_data(from,to,:,:,:,:)),5))',0.2,'g--',foi,1,1)
%     %             stdshade_TEC(squeeze(median(abs(all_controls_granger_data(to,from,:,:,:,:)),5))',0.2,'b:',foi,1,1)
%     %             stdshade_TEC(squeeze(median(abs(all_patients_granger_data(to,from,:,:,:,:)),5))',0.2,'b--',foi,1,1)
%     %             legend({[from_name '-' to_name ' control'],[from_name '-' to_name ' patient'],[to_name '-' from_name ' control'],[to_name '-' from_name ' patient']})
%     %             title(['By group median ' analysis_type ' time ' num2str(start_times) '-' num2str(end_times)])
%     %             for i = 1:topfreqband %Rough initial parametric stats - better to use permutation
%     %                 [p_tf(i) h_tf(i)] = ranksum(squeeze(mean(abs(all_controls_granger_data(from,to,:,i,:,:)),5)),squeeze(mean(abs(all_patients_granger_data(from,to,:,i,:,:)),5)));
%     %                 [p_ft(i) h_ft(i)] = ranksum(squeeze(mean(abs(all_controls_granger_data(to,from,:,i,:,:)),5)),squeeze(mean(abs(all_patients_granger_data(to,from,:,i,:,:)),5)));
%     %                 if h_tf(i) == 1
%     %                     plot(foi(i),0,'g*')
%     %                 end
%     %                 if h_ft(i) == 1
%     %                     plot(foi(i),0,'bx')
%     %                 end
%     %             end
%     %
%     %             config = figure;
%     %             patfig = figure;
%     %             this_con = 0;
%     %             this_pat = 0;
%     %
%     %             for s = 1:length(group)
%     %                 if group(s) == 1
%     %                     this_con = this_con+1;
%     %                     figure(config);
%     %                     h = subplot(4,3,this_con);
%     %                     stdshade_TEC(squeeze(median(abs(all_patients_granger_data(to,from,:,:,:,:)),5))',0.2,'b--',foi,1,1)
%     %                     hold on
%     %                     stdshade_TEC(squeeze(median(abs(all_patients_granger_data(to,from,:,:,:,:)),5))',0.2,'b--',foi,1,1)
%     %                     for i = 1:topfreqband %Compare within subj permutation
%     %                         p_tf(i) = 1-relRankIn_includeValue_lowerBound(squeeze(mean(abs(all_random_granger_data(from,to,:,i,:,:,s)),5)),squeeze(mean(abs(all_granger_data(from,to,:,i,:,s)),5)));
%     %                         p_ft(i) = 1-relRankIn_includeValue_lowerBound(squeeze(mean(abs(all_random_granger_data(to,from,:,i,:,:,s)),5)),squeeze(mean(abs(all_granger_data(to,from,:,i,:,s)),5)));
%     %                         if p_tf(i)<=0.05 || p_tf(i)>=0.95
%     %                             plot(foi(i),0,'g*')
%     %                         end
%     %                         if p_ft(i)<=0.05 || p_ft(i)>=0.95
%     %                             plot(foi(i),0,'bx')
%     %                         end
%     %                     end
%     %
%     %                 elseif group(s) == 2
%     %                     this_pat = this_pat+1;
%     %                     figure(patfig);
%     %                     h = subplot(4,3,this_pat);
%     %                     stdshade_TEC(squeeze(median(abs(all_patients_granger_data(to,from,:,:,:,:)),5))',0.2,'b--',foi,1,1)
%     %                     hold on
%     %                     stdshade_TEC(squeeze(median(abs(all_patients_granger_data(to,from,:,:,:,:)),5))',0.2,'b--',foi,1,1)
%     %                     for i = 1:topfreqband %Compare within subj permutation
%     %                         p_tf(i) = 1-relRankIn_includeValue_lowerBound(squeeze(mean(abs(all_random_granger_data(from,to,:,i,:,:,s)),5)),squeeze(mean(abs(all_granger_data(from,to,:,i,:,s)),5)));
%     %                         p_ft(i) = 1-relRankIn_includeValue_lowerBound(squeeze(mean(abs(all_random_granger_data(to,from,:,i,:,:,s)),5)),squeeze(mean(abs(all_granger_data(to,from,:,i,:,s)),5)));
%     %                         if p_tf(i)<=0.05 || p_tf(i)>=0.95
%     %                             plot(foi(i),0,'g*')
%     %                         end
%     %                         if p_ft(i)<=0.05 || p_ft(i)>=0.95
%     %                             plot(foi(i),0,'bx')
%     %                         end
%     %                     end
%     %                 end
%     %             end
%     %             figure(config)
%     %             title(['Controls ' analysis_type ' time ' num2str(start_times) '-' num2str(end_times)])
%     %             figure(patfig)
%     %             title(['Patients ' analysis_type ' time ' num2str(start_times) '-' num2str(end_times)])
%     %
%     %
%     %             if closeafter == 1; close all; end; figure
%     %             stdshade_TEC(squeeze(all_absolute_mismatch_contrasts(from,to,:,:,:))',0.2,'g',foi,1,1)
%     %             hold on
%     %             stdshade_TEC(squeeze(all_absolute_mismatch_contrasts(to,from,:,:,:))',0.2,'b',foi,1,1)
%     %             plot(foi,squeeze(mean(all_absolute_random_mismatch_contrasts(from,to,:,:,:,:),6)),'k')
%     %             plot(foi,squeeze(mean(all_absolute_random_mismatch_contrasts(to,from,:,:,:,:),6)),'k')
%     %             stdshade_TEC(squeeze(all_absolute_mismatch_contrasts(from,to,:,:,:))',0.2,'g',foi,1,1)
%     %
%     %             stdshade_TEC(squeeze(all_absolute_mismatch_contrasts(to,from,:,:,:))',0.2,'b',foi,1,1)
%     %
%     %             title(['All Standard-Deviant ' analysis_type ' time ' num2str(start_times) '-' num2str(end_times)])
%     %             legend({[from_name '-' to_name],[to_name '-' from_name]})
%     %             for i = 1:topfreqband %Compare across group permutation
%     %                 p_tf(i) = 1-relRankIn_includeValue_lowerBound(squeeze(mean(all_absolute_random_mismatch_contrasts(from,to,:,i,:,:),6)),squeeze(mean(all_absolute_mismatch_contrasts(from,to,:,i,:),5)));
%     %                 p_ft(i) = 1-relRankIn_includeValue_lowerBound(squeeze(mean(all_absolute_random_mismatch_contrasts(to,from,:,i,:,:),6)),squeeze(mean(all_absolute_mismatch_contrasts(to,from,:,i,:),5)));
%     %                 if p_tf(i)<=0.05 || p_tf(i)>=0.95
%     %                     plot(foi(i),0,'g*')
%     %                 end
%     %                 if p_ft(i)<=0.05 || p_ft(i)>=0.95
%     %                     plot(foi(i),0,'bx')
%     %                 end
%     %             end
%     %
%     %             if closeafter == 1; close all; end; figure
%     %             stdshade_TEC(squeeze(all_absolute_mismatch_contrasts(from,to,:,:,:))',0.2,'g',foi,1,1)
%     %             hold on
%     %             stdshade_TEC(squeeze(all_absolute_mismatch_contrasts(to,from,:,:,:))',0.2,'b',foi,1,1)
%     %             plot(foi,squeeze(median(all_absolute_random_mismatch_contrasts(from,to,:,:,:,:),6)),'k')
%     %             plot(foi,squeeze(median(all_absolute_random_mismatch_contrasts(to,from,:,:,:,:),6)),'k')
%     %             stdshade_TEC(squeeze(all_absolute_mismatch_contrasts(from,to,:,:,:))',0.2,'g',foi,1,1)
%     %
%     %             stdshade_TEC(squeeze(all_absolute_mismatch_contrasts(to,from,:,:,:))',0.2,'b',foi,1,1)
%     %
%     %             title(['All median Standard-Deviant ' analysis_type ' time ' num2str(start_times) '-' num2str(end_times)])
%     %             legend({[from_name '-' to_name],[to_name '-' from_name]})
%     %             for i = 1:topfreqband %Compare across group permutation
%     %                 p_tf(i) = 1-relRankIn_includeValue_lowerBound(squeeze(median(all_absolute_random_mismatch_contrasts(from,to,:,i,:,:),6)),squeeze(median(all_absolute_mismatch_contrasts(from,to,:,i,:),5)));
%     %                 p_ft(i) = 1-relRankIn_includeValue_lowerBound(squeeze(median(all_absolute_random_mismatch_contrasts(to,from,:,i,:,:),6)),squeeze(median(all_absolute_mismatch_contrasts(to,from,:,i,:),5)));
%     %                 if p_tf(i)<=0.05 || p_tf(i)>=0.95
%     %                     plot(foi(i),0,'g*')
%     %                 end
%     %                 if p_ft(i)<=0.05 || p_ft(i)>=0.95
%     %                     plot(foi(i),0,'bx')
%     %                 end
%     %             end
%     %
%     %             if closeafter == 1; close all; end; figure
%     %             stdshade_TEC(squeeze(all_absolute_clarity_contrasts(from,to,:,:,:))',0.2,'g',foi,1,1)
%     %             hold on
%     %             stdshade_TEC(squeeze(all_absolute_clarity_contrasts(to,from,:,:,:))',0.2,'b',foi,1,1)
%     %             plot(foi,squeeze(mean(all_absolute_random_clarity_contrasts(from,to,:,:,:,:),6)),'k')
%     %             plot(foi,squeeze(mean(all_absolute_random_clarity_contrasts(to,from,:,:,:,:),6)),'k')
%     %             stdshade_TEC(squeeze(all_absolute_clarity_contrasts(from,to,:,:,:))',0.2,'g',foi,1,1)
%     %
%     %             stdshade_TEC(squeeze(all_absolute_clarity_contrasts(to,from,:,:,:))',0.2,'b',foi,1,1)
%     %
%     %             title(['All Clear-Unclear ' analysis_type ' time ' num2str(start_times) '-' num2str(end_times)])
%     %             legend({[from_name '-' to_name],[to_name '-' from_name]})
%     %             for i = 1:topfreqband %Compare across group permutation
%     %                 p_tf(i) = 1-relRankIn_includeValue_lowerBound(squeeze(mean(all_absolute_random_clarity_contrasts(from,to,:,i,:,:),6)),squeeze(mean(all_absolute_clarity_contrasts(from,to,:,i,:),5)));
%     %                 p_ft(i) = 1-relRankIn_includeValue_lowerBound(squeeze(mean(all_absolute_random_clarity_contrasts(to,from,:,i,:,:),6)),squeeze(mean(all_absolute_clarity_contrasts(to,from,:,i,:),5)));
%     %                 if p_tf(i)<=0.05 || p_tf(i)>=0.95
%     %                     plot(foi(i),0,'g*')
%     %                 end
%     %                 if p_ft(i)<=0.05 || p_ft(i)>=0.95
%     %                     plot(foi(i),0,'bx')
%     %                 end
%     %             end
%     %
%     %             if closeafter == 1; close all; end; figure
%     %             stdshade_TEC(squeeze(all_absolute_clarity_contrasts(from,to,:,:,:))',0.2,'g',foi,1,1)
%     %             hold on
%     %             stdshade_TEC(squeeze(all_absolute_clarity_contrasts(to,from,:,:,:))',0.2,'b',foi,1,1)
%     %             plot(foi,squeeze(median(all_absolute_random_clarity_contrasts(from,to,:,:,:,:),6)),'k')
%     %             plot(foi,squeeze(median(all_absolute_random_clarity_contrasts(to,from,:,:,:,:),6)),'k')
%     %             stdshade_TEC(squeeze(all_absolute_clarity_contrasts(from,to,:,:,:))',0.2,'g',foi,1,1)
%     %
%     %             stdshade_TEC(squeeze(all_absolute_clarity_contrasts(to,from,:,:,:))',0.2,'b',foi,1,1)
%     %
%     %             title(['All median Clear-Unclear ' analysis_type ' time ' num2str(start_times) '-' num2str(end_times)])
%     %             legend({[from_name '-' to_name],[to_name '-' from_name]})
%     %             for i = 1:topfreqband %Compare across group permutation
%     %                 p_tf(i) = 1-relRankIn_includeValue_lowerBound(squeeze(median(all_absolute_random_clarity_contrasts(from,to,:,i,:,:),6)),squeeze(median(all_absolute_clarity_contrasts(from,to,:,i,:),5)));
%     %                 p_ft(i) = 1-relRankIn_includeValue_lowerBound(squeeze(median(all_absolute_random_clarity_contrasts(to,from,:,i,:,:),6)),squeeze(median(all_absolute_clarity_contrasts(to,from,:,i,:),5)));
%     %                 if p_tf(i)<=0.05 || p_tf(i)>=0.95
%     %                     plot(foi(i),0,'g*')
%     %                 end
%     %                 if p_ft(i)<=0.05 || p_ft(i)>=0.95
%     %                     plot(foi(i),0,'bx')
%     %                 end
%     %             end
%     %
%     %             if closeafter == 1; close all; end; figure
%     %             stdshade_TEC(squeeze(all_absolute_interaction_contrasts(from,to,:,:,:))',0.2,'g',foi,1,1)
%     %             hold on
%     %             stdshade_TEC(squeeze(all_absolute_interaction_contrasts(to,from,:,:,:))',0.2,'b',foi,1,1)
%     %             title(['All Interaction ' analysis_type ' time ' num2str(start_times) '-' num2str(end_times)])
%     %             legend({[from_name '-' to_name],[to_name '-' from_name]})
%     %             for i = 1:topfreqband %Compare across group permutation
%     %                 p_tf(i) = 1-relRankIn_includeValue_lowerBound(squeeze(mean(all_absolute_random_interaction_contrasts(from,to,:,i,:,:),6)),squeeze(mean(all_absolute_interaction_contrasts(from,to,:,i,:),5)));
%     %                 p_ft(i) = 1-relRankIn_includeValue_lowerBound(squeeze(mean(all_absolute_random_interaction_contrasts(to,from,:,i,:,:),6)),squeeze(mean(all_absolute_interaction_contrasts(to,from,:,i,:),5)));
%     %                 if p_tf(i)<=0.05 || p_tf(i)>=0.95
%     %                     plot(foi(i),0,'g*')
%     %                 end
%     %                 if p_ft(i)<=0.05 || p_ft(i)>=0.95
%     %                     plot(foi(i),0,'bx')
%     %                 end
%     %             end
%     %
%     %             if closeafter == 1; close all; end; figure
%     %             stdshade_TEC(squeeze(all_absolute_interaction_contrasts(from,to,:,:,:))',0.2,'g',foi,1,1)
%     %             hold on
%     %             stdshade_TEC(squeeze(all_absolute_interaction_contrasts(to,from,:,:,:))',0.2,'b',foi,1,1)
%     %             title(['All median Interaction ' analysis_type ' time ' num2str(start_times) '-' num2str(end_times)])
%     %             legend({[from_name '-' to_name],[to_name '-' from_name]})
%     %             for i = 1:topfreqband %Compare across group permutation
%     %                 p_tf(i) = 1-relRankIn_includeValue_lowerBound(squeeze(median(all_absolute_random_interaction_contrasts(from,to,:,i,:,:),6)),squeeze(median(all_absolute_interaction_contrasts(from,to,:,i,:),5)));
%     %                 p_ft(i) = 1-relRankIn_includeValue_lowerBound(squeeze(median(all_absolute_random_interaction_contrasts(to,from,:,i,:,:),6)),squeeze(median(all_absolute_interaction_contrasts(to,from,:,i,:),5)));
%     %                 if p_tf(i)<=0.05 || p_tf(i)>=0.95
%     %                     plot(foi(i),0,'g*')
%     %                 end
%     %                 if p_ft(i)<=0.05 || p_ft(i)>=0.95
%     %                     plot(foi(i),0,'bx')
%     %                 end
%     %             end
%     %
%     %             if closeafter == 1; close all; end; figure
%     %             hold on
%     %             stdshade_TEC(squeeze(all_absolute_interaction_contrasts(to,from,:,:,:))',0.2,'b',foi,1,1)
%     %             stdshade_TEC(squeeze(all_absolute_interaction_contrasts(to,from,:,:,:))',0.2,'b',foi,1,1)
%     %             stdshade_TEC(squeeze(all_absolute_interaction_contrasts(to,from,:,:,:))',0.2,'b',foi,1,1)
%     %             stdshade_TEC(squeeze(all_absolute_interaction_contrasts(to,from,:,:,:))',0.2,'b',foi,1,1)
%     %
%     %             title(['Group by Mean Standard-Deviant ' analysis_type ])
%     %             legend({[from_name '-' to_name ' control'],[from_name '-' to_name ' patient'],[to_name '-' from_name ' control'],[to_name '-' from_name ' patient']})
%     %             for i = 1:topfreqband %Rough initial parametric stats - better to use permutation
%     %                 [p_tf(i) h_tf(i)] = ttest2(squeeze(all_absolute_patients_mismatch_contrasts(from,to,:,i,:)),squeeze(all_absolute_controls_mismatch_contrasts(from,to,:,i,:)));
%     %                 [p_ft(i) h_ft(i)] = ttest2(squeeze(all_absolute_patients_mismatch_contrasts(to,from,:,i,:)),squeeze(all_absolute_controls_mismatch_contrasts(to,from,:,i,:)));
%     %                 if h_tf(i) == 1
%     %                     plot(foi(i),0,'g*')
%     %                 end
%     %                 if h_ft(i) == 1
%     %                     plot(foi(i),0,'bx')
%     %                 end
%     %             end
%     %
%     %             if closeafter == 1; close all; end; figure
%     %             hold on
%     %             stdshade_TEC(squeeze(all_absolute_interaction_contrasts(to,from,:,:,:))',0.2,'b',foi,1,1)
%     %             stdshade_TEC(squeeze(all_absolute_interaction_contrasts(to,from,:,:,:))',0.2,'b',foi,1,1)
%     %             stdshade_TEC(squeeze(all_absolute_interaction_contrasts(to,from,:,:,:))',0.2,'b',foi,1,1)
%     %             stdshade_TEC(squeeze(all_absolute_interaction_contrasts(to,from,:,:,:))',0.2,'b',foi,1,1)
%     %
%     %             title(['Group by Median Standard-Deviant ' analysis_type ])
%     %             legend({[from_name '-' to_name ' control'],[from_name '-' to_name ' patient'],[to_name '-' from_name ' control'],[to_name '-' from_name ' patient']})
%     %             for i = 1:topfreqband %Rough initial parametric stats - better to use permutation
%     %                 [p_tf(i) h_tf(i)] = ranksum(squeeze(all_absolute_patients_mismatch_contrasts(from,to,:,i,:)),squeeze(all_absolute_controls_mismatch_contrasts(from,to,:,i,:)));
%     %                 [p_ft(i) h_ft(i)] = ranksum(squeeze(all_absolute_patients_mismatch_contrasts(to,from,:,i,:)),squeeze(all_absolute_controls_mismatch_contrasts(to,from,:,i,:)));
%     %                 if h_tf(i) == 1
%     %                     plot(foi(i),0,'g*')
%     %                 end
%     %                 if h_ft(i) == 1
%     %                     plot(foi(i),0,'bx')
%     %                 end
%     %             end
%     %
%     %             if closeafter == 1; close all; end; figure
%     %             hold on
%     %             stdshade_TEC(squeeze(all_absolute_interaction_contrasts(to,from,:,:,:))',0.2,'b',foi,1,1)
%     %             stdshade_TEC(squeeze(all_absolute_interaction_contrasts(to,from,:,:,:))',0.2,'b',foi,1,1)
%     %             stdshade_TEC(squeeze(all_absolute_interaction_contrasts(to,from,:,:,:))',0.2,'b',foi,1,1)
%     %             stdshade_TEC(squeeze(all_absolute_interaction_contrasts(to,from,:,:,:))',0.2,'b',foi,1,1)
%     %             plot(foi,squeeze(mean(all_absolute_random_controls_mismatch_contrasts(from,to,:,:,:,:),6)),'k:')
%     %             plot(foi,squeeze(mean(all_absolute_random_patients_mismatch_contrasts(from,to,:,:,:,:),6)),'k--')
%     %             plot(foi,squeeze(mean(all_absolute_random_controls_mismatch_contrasts(to,from,:,:,:,:),6)),'k:')
%     %             plot(foi,squeeze(mean(all_absolute_random_patients_mismatch_contrasts(to,from,:,:,:,:),6)),'k--')
%     %             stdshade_TEC(squeeze(all_absolute_interaction_contrasts(to,from,:,:,:))',0.2,'b',foi,1,1)
%     %             stdshade_TEC(squeeze(all_absolute_interaction_contrasts(to,from,:,:,:))',0.2,'b',foi,1,1)
%     %             stdshade_TEC(squeeze(all_absolute_interaction_contrasts(to,from,:,:,:))',0.2,'b',foi,1,1)
%     %             stdshade_TEC(squeeze(all_absolute_interaction_contrasts(to,from,:,:,:))',0.2,'b',foi,1,1)
%     %
%     %             title(['By group sig Mean Standard-Deviant ' analysis_type ])
%     %             legend({[from_name '-' to_name ' control'],[from_name '-' to_name ' patient'],[to_name '-' from_name ' control'],[to_name '-' from_name ' patient']})
%     %             for i = 1:topfreqband %Rough initial parametric stats - better to use permutation
%     %                 p_tf_c(i) = 1-relRankIn_includeValue_lowerBound(squeeze(mean(all_absolute_random_controls_mismatch_contrasts(from,to,i,:,:,:),6)),squeeze(mean(all_absolute_controls_mismatch_contrasts(from,to,:,i,:),5)));
%     %                 p_ft_c(i) = 1-relRankIn_includeValue_lowerBound(squeeze(mean(all_absolute_random_controls_mismatch_contrasts(to,from,i,:,:,:),6)),squeeze(mean(all_absolute_controls_mismatch_contrasts(to,from,:,i,:),5)));
%     %                 if p_tf_c(i)<=0.05 || p_ft_p(i)>=0.95
%     %                     plot(foi(i),0,'go')
%     %                 end
%     %                 if p_ft_c(i)<=0.05 || p_ft_c(i)>=0.95
%     %                     plot(foi(i),0,'b.')
%     %                 end
%     %                 p_tf_p(i) = 1-relRankIn_includeValue_lowerBound(squeeze(mean(all_absolute_random_patients_mismatch_contrasts(from,to,i,:,:,:),6)),squeeze(mean(all_absolute_patients_mismatch_contrasts(from,to,:,i,:),5)));
%     %                 p_ft_p(i) = 1-relRankIn_includeValue_lowerBound(squeeze(mean(all_absolute_random_patients_mismatch_contrasts(to,from,i,:,:,:),6)),squeeze(mean(all_absolute_patients_mismatch_contrasts(to,from,:,i,:),5)));
%     %                 if p_tf_p(i)<=0.05 || p_ft_p(i)>=0.95
%     %                     plot(foi(i),0,'gp')
%     %                 end
%     %                 if p_ft_p(i)<=0.05 || p_ft_p(i)>=0.95
%     %                     plot(foi(i),0,'bs')
%     %                 end
%     %             end
%     %
%     %             if closeafter == 1; close all; end; figure
%     %             hold on
%     %             stdshade_TEC(squeeze(all_absolute_interaction_contrasts(to,from,:,:,:))',0.2,'b',foi,1,1)
%     %             stdshade_TEC(squeeze(all_absolute_interaction_contrasts(to,from,:,:,:))',0.2,'b',foi,1,1)
%     %             stdshade_TEC(squeeze(all_absolute_interaction_contrasts(to,from,:,:,:))',0.2,'b',foi,1,1)
%     %             stdshade_TEC(squeeze(all_absolute_interaction_contrasts(to,from,:,:,:))',0.2,'b',foi,1,1)
%     %             plot(foi,squeeze(median(all_absolute_random_controls_mismatch_contrasts(from,to,:,:,:,:),6)),'k:')
%     %             plot(foi,squeeze(median(all_absolute_random_patients_mismatch_contrasts(from,to,:,:,:,:),6)),'k--')
%     %             plot(foi,squeeze(median(all_absolute_random_controls_mismatch_contrasts(to,from,:,:,:,:),6)),'k:')
%     %             plot(foi,squeeze(median(all_absolute_random_patients_mismatch_contrasts(to,from,:,:,:,:),6)),'k--')
%     %             stdshade_TEC(squeeze(all_absolute_interaction_contrasts(to,from,:,:,:))',0.2,'b',foi,1,1)
%     %             stdshade_TEC(squeeze(all_absolute_interaction_contrasts(to,from,:,:,:))',0.2,'b',foi,1,1)
%     %             stdshade_TEC(squeeze(all_absolute_interaction_contrasts(to,from,:,:,:))',0.2,'b',foi,1,1)
%     %             stdshade_TEC(squeeze(all_absolute_interaction_contrasts(to,from,:,:,:))',0.2,'b',foi,1,1)
%     %
%     %             title(['By group sig Median Standard-Deviant ' analysis_type ])
%     %             legend({[from_name '-' to_name ' control'],[from_name '-' to_name ' patient'],[to_name '-' from_name ' control'],[to_name '-' from_name ' patient']})
%     %             for i = 1:topfreqband %Rough initial parametric stats - better to use permutation
%     %                 p_tf_c(i) = 1-relRankIn_includeValue_lowerBound(squeeze(median(all_absolute_random_controls_mismatch_contrasts(from,to,i,:,:,:),6)),squeeze(median(all_absolute_controls_mismatch_contrasts(from,to,:,i,:),5)));
%     %                 p_ft_c(i) = 1-relRankIn_includeValue_lowerBound(squeeze(median(all_absolute_random_controls_mismatch_contrasts(to,from,i,:,:,:),6)),squeeze(median(all_absolute_controls_mismatch_contrasts(to,from,:,i,:),5)));
%     %                 if p_tf_c(i)<=0.05 || p_ft_p(i)>=0.95
%     %                     plot(foi(i),0,'go')
%     %                 end
%     %                 if p_ft_c(i)<=0.05 || p_ft_c(i)>=0.95
%     %                     plot(foi(i),0,'b.')
%     %                 end
%     %                 p_tf_p(i) = 1-relRankIn_includeValue_lowerBound(squeeze(median(all_absolute_random_patients_mismatch_contrasts(from,to,i,:,:,:),6)),squeeze(median(all_absolute_patients_mismatch_contrasts(from,to,:,i,:),5)));
%     %                 p_ft_p(i) = 1-relRankIn_includeValue_lowerBound(squeeze(median(all_absolute_random_patients_mismatch_contrasts(to,from,i,:,:,:),6)),squeeze(median(all_absolute_patients_mismatch_contrasts(to,from,:,i,:),5)));
%     %                 if p_tf_p(i)<=0.05 || p_ft_p(i)>=0.95
%     %                     plot(foi(i),0,'gp')
%     %                 end
%     %                 if p_ft_p(i)<=0.05 || p_ft_p(i)>=0.95
%     %                     plot(foi(i),0,'bs')
%     %                 end
%     %             end
%     %
%     %     end
    end
end
    