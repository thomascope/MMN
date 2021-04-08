function Plot_DCM_Report(start_date,end_date)
reports_directory = '/imaging/mlr/users/tc02/Holly_MMN/DCM_REPORTS'
start_date = datetime(start_date);
end_date = datetime(end_date);
all_folders = dir(reports_directory);
for i = 3:size(all_folders,1) %Ignore . and ..
    if all_folders(i).date >= start_date && all_folders(i).date <= end_date
        these_files = dir([reports_directory filesep all_folders(i).name]);
        if exist('all_files','var')
            all_files = [all_files; these_files(3:end)];
        else
            all_files = these_files(3:end);
        end
    end
end
disp([num2str(size(all_files,1)) ' files found, proceeding'])
nodelist = {'left A1', 'left STG', 'left IFG', 'left IPC', 'right A1', 'right STG', 'right IFG', 'right IPC'};
figure
for i = size(all_files,1):-1:1
    clf
    disp(['plotting ' all_files(i).folder filesep all_files(i).name])
    all_files(i).name
    load([all_files(i).folder filesep all_files(i).name]);
    for j=1:8
        
        subplot(4,2,j)
        hold on
        plot(0:4:400, R.y(:,j,end))
        plot(0:4:400, R.yp(:,j,end), 'r')
        title(nodelist{j})
    end
    pause
end
