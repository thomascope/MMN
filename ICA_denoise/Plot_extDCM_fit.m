filelist = dir('/imaging/tc02/Holly_MMN/extDCMs/*DVT*');
all_STDs_modelled = zeros(101,8,length(filelist));
all_STDs_data = zeros(101,8,length(filelist));
all_DVTs_modelled = zeros(101,8,length(filelist));
all_DVTs_data = zeros(101,8,length(filelist));

for j = 1:length(filelist)
    D = load([filelist(j).folder filesep strrep(filelist(j).name,'DVT','STD')]);
    S = load([filelist(j).folder filesep filelist(j).name]);
    %     figure('Position',[10 10 1010 610])
    %     for i = 1:8
    %         subplot(4,4,i)
    %         plot(0:4:400, D.DCM.H{1}(:,i))
    %         hold on
    %         plot(0:4:400, D.DCM.xY.y{1}(:,i))
    %         title(D.DCM.xY.name{i})
    %     end
    %     for i = 1:8
    %         subplot(4,4,i+8)
    %         plot(0:4:400, S.DCM.H{1}(:,i))
    %         hold on
    %         plot(0:4:400, S.DCM.xY.y{1}(:,i))
    %         title(S.DCM.xY.name{i})
    %     end
    if j == 59 || j == 67 %Debugging of NaNs
        pause
    end
    all_STDs_modelled(:,:,j) = D.DCM.H{1};
    all_STDs_data(:,:,j) = D.DCM.xY.y{1};
    all_DVTs_modelled(:,:,j) = S.DCM.H{1};
    all_DVTs_data(:,:,j) = S.DCM.xY.y{1};
    
    %pause
end

figure('Position',[10 10 1010 610])
for i = 1:8
    subplot(4,4,i)
    plot(0:4:400, squeeze(nanmean(all_STDs_modelled(:,i,:),3)))
    hold on
    plot(0:4:400, squeeze(mean(all_STDs_data(:,i,:),3)))
    title(D.DCM.xY.name{i})
end
for i = 1:8
    subplot(4,4,i+8)
    plot(0:4:400, squeeze(mean(all_DVTs_modelled(:,i,:),3)))
    hold on
    plot(0:4:400, squeeze(mean(all_DVTs_data(:,i,:),3)))
    title(S.DCM.xY.name{i})
end