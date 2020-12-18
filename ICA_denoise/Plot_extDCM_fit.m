addpath('/group/language/data/thomascope/MMN/ICA_denoise/Tallie_extDCM/mfiles_also_needed')

filelist = dir('/imaging/tc02/Holly_MMN/extDCMs/*DVT*');
%filelist = dir('/imaging/tc02/Holly_MMN/extDCMs_first_attempt/*DVT*');
Inverted_Conditions = {'STD','DVT','location','intensity','duration','gap','frequency'};

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

for i = 1:length(Inverted_Conditions)
    figure('Position',[10 10 1010 610])
    for j = 1:8
        subplot(2,4,j)
        plot(D{i}.DCM.xY.pst, squeeze(squeeze(nanmean(all_data(:,j,:,i),3))))
        hold on
        plot(D{i}.DCM.xY.pst, squeeze(squeeze(nanmean(all_modelled(:,j,:,i),3))))
        title(D{1}.DCM.xY.name{j})
        sgtitle(Inverted_Conditions{i})
    end
end

figure('Position',[10 10 1010 610])
for i = 1:2
    for j = 1:8
        subplot(4,4,j+((i-1)*8))
        plot(D{i}.DCM.xY.pst, squeeze(squeeze(nanmean(all_data(:,j,:,i),3))),'r-')
        hold on
        plot(D{i}.DCM.xY.pst, squeeze(squeeze(nanmean(all_modelled(:,j,:,i),3))),'b-')
        title(D{1}.DCM.xY.name{j})
        sgtitle(Inverted_Conditions{i})
    end
end