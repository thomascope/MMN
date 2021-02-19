addpath('/group/language/data/thomascope/MMN/ICA_denoise/Tallie_extDCM/mfiles_also_needed')

filelist = dir([p.extDCM_outdir '/*DVT*']);
%filelist = dir('/imaging/tc02/Holly_MMN/extDCMs_first_attempt/*DVT*');
%Inverted_Conditions = {'STD','DVT','location','intensity','duration','gap','frequency'};
Inverted_Conditions = {'STD','DVT'};

all_modelled = zeros(101,8,length(filelist),length(Inverted_Conditions));
all_data = zeros(101,8,length(filelist),length(Inverted_Conditions));

clear D
figure('Position',[10 10 1010 610])
for j = 1:length(filelist)
    for i = 1:length(Inverted_Conditions)
        try
            D{i} = load([filelist(j).folder filesep strrep(filelist(j).name,'DVT',Inverted_Conditions{i})]);
        catch
            disp(['Missing file ' strrep(filelist(j).name,'DVT',Inverted_Conditions{i})])
        end
        %         figure('Position',[10 10 1010 610])
        %         for i = 1:8
        %             subplot(4,4,i)
        %             plot(0:4:400, D.DCM.H{1}(:,i))
        %             hold on
        %             plot(0:4:400, D.DCM.xY.y{1}(:,i))
        %             title(D.DCM.xY.name{i})
        %         end
        %         for i = 1:8
        %             subplot(4,4,i+8)
        %             plot(0:4:400, S.DCM.H{1}(:,i))
        %             hold on
        %             plot(0:4:400, S.DCM.xY.y{1}(:,i))
        %             title(S.DCM.xY.name{i})
        %         end
        if any(any(isnan(D{i}.DCM.H{1}))) %Debugging of NaNs
            disp(['NaNs found in file ' strrep(filelist(j).name,'DVT',Inverted_Conditions{i})])
            %pause
        end
        
        all_modelled(:,:,j,i) = D{i}.DCM.H{1};
        all_data(:,:,j,i) = D{i}.DCM.xY.y{1};
        
        for k = 1:length(Participant)
            if contains(filelist(j).name, Participant{k}.namepostmerge)
                all_dcm_groups{j} = Participant{k}.diag;
                continue
            end
        end
                
        
%         if D{1}.DCM.xY.y{1}(findnearest(D{1}.DCM.xY.pst,100),2)<0
%             all_modelled(:,:,j,i) = D{i}.DCM.H{1};
%             all_data(:,:,j,i) = D{i}.DCM.xY.y{1};
%             for k = 1:8
%                 subplot(2,4,k)
%                 plot(D{i}.DCM.xY.pst, squeeze(squeeze(all_data(:,k,j,i))))
%                 hold on
%                 plot(D{i}.DCM.xY.pst, squeeze(squeeze(all_modelled(:,k,j,i))))
%                 title(D{1}.DCM.xY.name{k})
%                 sgtitle([Inverted_Conditions{i} ' subject ' num2str(j)])
%             end
%             drawnow
%             pause
%         else
%             if i == 1
%                 disp(['M100 positive for subject ' num2str(j) ', flipping data'])
%             end
%             all_modelled(:,:,j,i) = -D{i}.DCM.H{1};
%             all_data(:,:,j,i) = -D{i}.DCM.xY.y{1};
%         end
        
        %pause
    end
end

if any(isempty(all_dcm_groups))
    error('Some group memberships not found')
end

%For each condition individually if required
% for i = 1:length(Inverted_Conditions)
%     figure('Position',[10 10 1010 610])
%     for j = 1:8
%         subplot(2,4,j)
%         plot(D{i}.DCM.xY.pst, squeeze(squeeze(nanmean(all_data(:,j,:,i),3))))
%         hold on
%         plot(D{i}.DCM.xY.pst, squeeze(squeeze(nanmean(all_modelled(:,j,:,i),3))))
%         title(D{1}.DCM.xY.name{j})
%         sgtitle(Inverted_Conditions{i})
%     end
% end

%First plot STD and DVT
figure('Position',[10 10 1010 610])
for i = 1:2
    for j = 1:8
        subplot(6,4,j+((i-1)*8))
        plot(D{i}.DCM.xY.pst, squeeze(squeeze(nanmean(all_data(:,j,:,i),3))),'r-')
        hold on
        plot(D{i}.DCM.xY.pst, squeeze(squeeze(nanmean(all_modelled(:,j,:,i),3))),'b-')
        title([D{1}.DCM.xY.name{j} ' ' Inverted_Conditions{i}])
    end
end
%Now plot MMN (looks weird if not mean centred)
for j = 1:8
    subplot(6,4,j+(i*8))
    plot(D{i}.DCM.xY.pst, squeeze(squeeze(nanmean(all_data(:,j,:,i-1)-all_data(:,j,:,i),3))),'r-')
    hold on
    plot(D{i}.DCM.xY.pst, squeeze(squeeze(nanmean(all_modelled(:,j,:,i-1)-all_modelled(:,j,:,i),3))),'b-')
    title([D{1}.DCM.xY.name{j} ' Mismatch'])
end

%% Now plot STD and DVT with errorbars
figure('Position',[10 10 1010 610])
clear linehandle
for i = 1:2
    for j = 1:8
        subplot(4,4,j+((i-1)*8))
        hold on
        linehandle(1) = stdshade_TEC(squeeze(all_data(:,j,:,i))',0.2,'r',D{1}.DCM.xY.pst,1,1);
        linehandle(2) = stdshade_TEC(squeeze(all_modelled(:,j,:,i))',0.2,'b',D{1}.DCM.xY.pst,1,1);
        
        if j==4 && i==1
            hleg = legend(linehandle,{'Data','Model'});
        end
        title([D{1}.DCM.xY.name{j} ' ' Inverted_Conditions{i}])
    end
end
saveas(gcf,'./figures/Overall_DCM_fits.png')
saveas(gcf,'./figures/Overall_DCM_fits.pdf')
set(hleg,'visible','off')
saveas(gcf,'./figures/Overall_DCM_fits_nolegend.png')
saveas(gcf,'./figures/Overall_DCM_fits_nolegend.pdf')

%% Now plot STD and DVT with errorbars by group
for ss = 1:length(Participant)
    diagnosis{ss} = Participant{ss}.diag;
end
[dcm_group_names,~, group_inds] = unique(diagnosis,'stable');
cmap = colormap(parula(length(dcm_group_names)));

for k = 1:length(unique(all_dcm_groups))
    figure('Position',[10 10 1010 610])
    clear linehandle
    for i = 1:2
        for j = 1:8
            subplot(4,4,j+((i-1)*8))
            hold on
            linehandle(1) = stdshade_TEC(squeeze(all_data(:,j,strcmp(all_dcm_groups,dcm_group_names{k}),i))',0.2,'r',D{1}.DCM.xY.pst,1,1);
            linehandle(2) = stdshade_TEC(squeeze(all_modelled(:,j,strcmp(all_dcm_groups,dcm_group_names{k}),i))',0.2,'b',D{1}.DCM.xY.pst,1,1);
            
            if j==4
                legend(linehandle,{'Data','Model'})
            end
            title([D{1}.DCM.xY.name{j} ' ' Inverted_Conditions{i}])
        end
    end
    sgtitle(dcm_group_names{k})
    saveas(gcf,['./figures/' dcm_group_names{k} '_DCM_fits.png'])
    saveas(gcf,['./figures/' dcm_group_names{k} '_DCM_fits.pdf'])
    set(hleg,'visible','off')
    saveas(gcf,['./figures/' dcm_group_names{k} '_DCM_fits_nolegend.png'])
    saveas(gcf,['./figures/' dcm_group_names{k} '_DCM_fits_nolegend.pdf'])
end

%% Now break down model residual by group.
for ss = 1:length(Participant)
diagnosis{ss} = Participant{ss}.diag;
end
[dcm_group_names,~, group_inds] = unique(diagnosis,'stable');
cmap = colormap(parula(length(dcm_group_names)));
figure('Position',[10 10 1010 610])
clear linehandle

for i = 1:2
    for j = 1:8
        subplot(4,4,j+((i-1)*8))
        hold on
        for k = 1:length(unique(all_dcm_groups))
        linehandle(k) = stdshade_TEC(squeeze(all_data(:,j,strcmp(all_dcm_groups,dcm_group_names{k}),i))'-squeeze(all_modelled(:,j,strcmp(all_dcm_groups,dcm_group_names{k}),i))',0.2,cmap(k,:),D{1}.DCM.xY.pst,1,1);
        end
        if j==4 && i==1
            hleg = legend(linehandle,dcm_group_names);
        end
        title([D{1}.DCM.xY.name{j} ' ' Inverted_Conditions{i}])
    end
end
sgtitle('Model Residuals')

saveas(gcf,'./figures/Model_Residuals.png')
saveas(gcf,'./figures/Model_Residuals.pdf')
set(hleg,'visible','off')
saveas(gcf,'./figures/Model_Residuals_nolegend.png')
saveas(gcf,'./figures/Model_Residuals_nolegend.pdf')

%% Now break STD and DVT down by group.
for ss = 1:length(Participant)
diagnosis{ss} = Participant{ss}.diag;
end
[dcm_group_names,~, group_inds] = unique(diagnosis,'stable');
cmap = colormap(parula(length(dcm_group_names)));
figure('Position',[10 10 1010 610])
clear linehandle
for j = 1:8
    subplot(4,4,j)
    hold on
    for i = 1:length(unique(all_dcm_groups))
        linehandle(i) = stdshade_TEC(squeeze(all_modelled(:,j,strcmp(all_dcm_groups,dcm_group_names{i}),1)-all_modelled(:,j,strcmp(all_dcm_groups,dcm_group_names{i}),2))',0.2,cmap(i,:),D{1}.DCM.xY.pst,1,1);
    end
    if j==4
        legend(linehandle,dcm_group_names)
    end
    title([D{1}.DCM.xY.name{j} ' MMN Model'])
end

for j = 1:8
    subplot(4,4,j+8)
    hold on
    for i = 1:length(unique(all_dcm_groups))
        linehandle(i) = stdshade_TEC(squeeze(all_data(:,j,strcmp(all_dcm_groups,dcm_group_names{i}),1)-all_data(:,j,strcmp(all_dcm_groups,dcm_group_names{i}),2))',0.2,cmap(i,:),D{1}.DCM.xY.pst,1,1);
    end
    if j==4
        legend(linehandle,dcm_group_names)
    end
    title([D{1}.DCM.xY.name{j} ' MMN data'])
end


%% Now break mismatch down by group.
for ss = 1:length(Participant)
diagnosis{ss} = Participant{ss}.diag;
end
[dcm_group_names,~, group_inds] = unique(diagnosis,'stable');
cmap = colormap(parula(length(dcm_group_names)));
figure('Position',[10 10 1010 610])
clear linehandle
for j = 1:8
    subplot(4,4,j)
    hold on
    for i = 1:length(unique(all_dcm_groups))
        linehandle(i) = stdshade_TEC(squeeze(all_modelled(:,j,strcmp(all_dcm_groups,dcm_group_names{i}),1)-all_modelled(:,j,strcmp(all_dcm_groups,dcm_group_names{i}),2))',0.2,cmap(i,:),D{1}.DCM.xY.pst,1,1);
    end
    if j==4
        legend(linehandle,dcm_group_names)
    end
    title([D{1}.DCM.xY.name{j} ' MMN Model'])
end

for j = 1:8
    subplot(4,4,j+8)
    hold on
    for i = 1:length(unique(all_dcm_groups))
        linehandle(i) = stdshade_TEC(squeeze(all_data(:,j,strcmp(all_dcm_groups,dcm_group_names{i}),1)-all_data(:,j,strcmp(all_dcm_groups,dcm_group_names{i}),2))',0.2,cmap(i,:),D{1}.DCM.xY.pst,1,1);
    end
    if j==4
        legend(linehandle,dcm_group_names)
    end
    title([D{1}.DCM.xY.name{j} ' MMN data'])
end

%% Do it again but trying to account for mean centring - I don't really understand why this doesn't work for the model.

Ns = size(all_modelled,1);
X0 = spm_dctmtx(Ns,1);
R = speye(Ns) - X0*X0';
R = R*sparse(diag((tukeywin(Ns,.3)+.025)))*R; % R*sparse(diag(hanning(Ns)))*R; % TA altered this to be more inclusive of the epoch of interest

for i = 1:size(all_modelled,3)
    for j = 1:size(all_modelled,2)
        for cond = 1:length(Inverted_Conditions)
            all_aligned_modelled(:,j,i,cond) = all_modelled(:,j,i,cond)'/R;
            all_aligned_data(:,j,i,cond) = all_data(:,j,i,cond)'/R;
        end
    end
end

figure('Position',[10 10 1010 610])
clear linehandle
for j = 1:8
    subplot(4,4,j)
    hold on
    for i = 1:length(unique(all_dcm_groups))
        linehandle(i) = stdshade_TEC(squeeze(all_aligned_modelled(:,j,strcmp(all_dcm_groups,dcm_group_names{i}),1)-all_aligned_modelled(:,j,strcmp(all_dcm_groups,dcm_group_names{i}),2))',0.2,cmap(i,:),D{1}.DCM.xY.pst,1,1);
    end
    if j==4
        legend(linehandle,dcm_group_names)
    end
    title([D{1}.DCM.xY.name{j} ' MMN Model'])
end

for j = 1:8
    subplot(4,4,j+8)
    hold on
    for i = 1:length(unique(all_dcm_groups))
        linehandle(i) = stdshade_TEC(squeeze(all_aligned_data(:,j,strcmp(all_dcm_groups,dcm_group_names{i}),1)-all_aligned_data(:,j,strcmp(all_dcm_groups,dcm_group_names{i}),2))',0.2,cmap(i,:),D{1}.DCM.xY.pst,1,1);
    end
    if j==4
        legend(linehandle,dcm_group_names)
    end
    title([D{1}.DCM.xY.name{j} ' MMN data'])
end

%% Now quantify MMNs
%first work out the mean difference between 125 and 175ms (this seems to be
%where the majority of the difference in the raw data is)
modelled_MMN = mean(squeeze(all_modelled(33:44,1,:,1)-all_modelled(33:44,1,:,2)));
data_MMN = mean(squeeze(all_data(33:44,1,:,1)-all_data(33:44,1,:,2)));

%modelled_MMN = max(abs(squeeze(all_modelled(33:44,1,:,1)-all_modelled(33:44,1,:,2))));
%data_MMN = max(abs(squeeze(all_data(33:44,1,:,1)-all_data(33:44,1,:,2))));

figure('Position',[10 10 610 1010])
clear linehandle
%sgtitle('MMN amplitude')
subplot(2,1,1)
hold on
for i = 1:length(unique(all_dcm_groups))
    scatter(repmat(i+0.1,1,sum(strcmp(all_dcm_groups,dcm_group_names{i}))),modelled_MMN(strcmp(all_dcm_groups,dcm_group_names{i})));
    
    
    [~,p] = ttest2(modelled_MMN(strcmp(all_dcm_groups,'Control')),modelled_MMN(strcmp(all_dcm_groups,dcm_group_names{i})));
    if p<0.05 
        errorbar(i-0.1,mean(modelled_MMN(strcmp(all_dcm_groups,dcm_group_names{i}))),std(modelled_MMN(strcmp(all_dcm_groups,dcm_group_names{i})))/sqrt(sum(strcmp(all_dcm_groups,dcm_group_names{i}))),'-s','MarkerSize',10,'MarkerEdgeColor','red','MarkerFaceColor','red','LineWidth',1,'color','black')
    else
        errorbar(i-0.1,mean(modelled_MMN(strcmp(all_dcm_groups,dcm_group_names{i}))),std(modelled_MMN(strcmp(all_dcm_groups,dcm_group_names{i})))/sqrt(sum(strcmp(all_dcm_groups,dcm_group_names{i}))),'-s','MarkerSize',10,'MarkerEdgeColor','black','MarkerFaceColor','black','LineWidth',1,'color','black')
    end
end
xlim([0 length(unique(all_dcm_groups))+1])
xticks([1:length(unique(all_dcm_groups))])
xticklabels(dcm_group_names)
xtickangle(30)
title('Modelled MMN','FontSize',34)
ylabel('Modelled MMN Amplitude (AU)')

subplot(2,1,2)
hold on
for i = 1:length(unique(all_dcm_groups))
    scatter(repmat(i+0.1,1,sum(strcmp(all_dcm_groups,dcm_group_names{i}))),data_MMN(strcmp(all_dcm_groups,dcm_group_names{i})));
    
    
    [~,p] = ttest2(data_MMN(strcmp(all_dcm_groups,'Control')),data_MMN(strcmp(all_dcm_groups,dcm_group_names{i})));
    if p<0.05 
        errorbar(i-0.1,mean(data_MMN(strcmp(all_dcm_groups,dcm_group_names{i}))),std(data_MMN(strcmp(all_dcm_groups,dcm_group_names{i})))/sqrt(sum(strcmp(all_dcm_groups,dcm_group_names{i}))),'-s','MarkerSize',10,'MarkerEdgeColor','red','MarkerFaceColor','red','LineWidth',1,'color','black')
    else
        errorbar(i-0.1,mean(data_MMN(strcmp(all_dcm_groups,dcm_group_names{i}))),std(data_MMN(strcmp(all_dcm_groups,dcm_group_names{i})))/sqrt(sum(strcmp(all_dcm_groups,dcm_group_names{i}))),'-s','MarkerSize',10,'MarkerEdgeColor','black','MarkerFaceColor','black','LineWidth',1,'color','black')
    end
end
xlim([0 length(unique(all_dcm_groups))+1])
xticks([1:length(unique(all_dcm_groups))])
xticklabels(dcm_group_names)
xtickangle(30)
title('Data MMN','FontSize',34)
ylabel('Data MMN Amplitude (AU)')

saveas(gcf,'./figures/DCM_MMN_Amplitudes_LA1.png')
saveas(gcf,'./figures/DCM_MMN_Amplitudes_LA1.pdf')