function [subjcondpair] = find_failed_extDCM_integrations(extDCM_directory,Inverted_Conditions,Participant)

addpath('/group/language/data/thomascope/MMN/ICA_denoise/Tallie_extDCM/mfiles_also_needed')

filelist = dir(extDCM_directory);
%filelist = dir('/imaging/tc02/Holly_MMN/extDCMs/*DVT*');
%filelist = dir('/imaging/tc02/Holly_MMN/extDCMs_first_attempt/*DVT*');
%Inverted_Conditions = {'STD','DVT','location','intensity','duration','gap','frequency'};

for j = 1:length(filelist)
    if contains(filelist(j).name,Participant{1}.namepostmerge) && contains(filelist(j).name,Inverted_Conditions{1})
        first_subj_file = j;
        break
    end
end
subjcondpair = {};
for j = 1:length(Participant)
    for i = 1:length(Inverted_Conditions)
        try
            D = load([filelist(j).folder filesep strrep(strrep(filelist(first_subj_file).name,Participant{1}.namepostmerge,Participant{j}.namepostmerge),Inverted_Conditions{1},Inverted_Conditions{i})]);
        catch
            disp(['Missing file ' strrep(strrep(filelist(first_subj_file).name,Participant{1}.namepostmerge,Participant{j}.namepostmerge),Inverted_Conditions{1},Inverted_Conditions{i})])
            subjcondpair(end+1,:) = {Participant{j}.namepostmerge,Inverted_Conditions{i}};
        end
        if any(any(isnan(D.DCM.H{1}))) %Debugging of NaNs
            disp(['NaNs found in file ' strrep(strrep(filelist(first_subj_file).name,Participant{1}.namepostmerge,Participant{j}.namepostmerge),Inverted_Conditions{1},Inverted_Conditions{i})])
            subjcondpair(end+1,:) = {Participant{j}.namepostmerge,Inverted_Conditions{i}};
        end
    end
end

mkdir('./extDCM_failures')

try
    T = cell2table(subjcondpair,'VariableNames',{'Subject','Condition'});
    todaysdate = datevec(date); 
writetable(T,['./extDCM_failures/failed_integrations' join([num2str(todaysdate(1)),num2str(todaysdate(2)),num2str(todaysdate(3))],'') '.csv']);
catch
    if isempty(subjcondpair)
        disp('Yay, no NaNs or missing files')
    else
        error('Something went wrong with writing the table but it was not that there were no dodgy files')
    end
end
