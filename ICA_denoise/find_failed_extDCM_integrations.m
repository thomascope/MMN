function [subjcondpair] = find_failed_extDCM_integrations(extDCM_directory,Inverted_Conditions,Participant)

addpath('/group/language/data/thomascope/MMN/ICA_denoise/Tallie_extDCM/mfiles_also_needed')

filelist = dir(extDCM_directory);
%filelist = dir('/imaging/tc02/Holly_MMN/extDCMs/*DVT*');
%filelist = dir('/imaging/tc02/Holly_MMN/extDCMs_first_attempt/*DVT*');
%Inverted_Conditions = {'STD','DVT','location','intensity','duration','gap','frequency'};

for j = 1:length(filelist)
    if contains(filelist(j).name,Participant{1}.name) && contains(filelist(j).name,Inverted_Conditions{1})
        first_subj_file = j;
        break
    end
end
subjcondpair = {};
for j = 1:length(Participant)
    for i = 1:length(Inverted_Conditions)
        try
            D = load([filelist(j).folder filesep strrep(strrep(filelist(first_subj_file).name,Participant{1}.name,Participant{j}.name),Inverted_Conditions{1},Inverted_Conditions{i})]);
        catch
            disp(['Missing file ' strrep(strrep(filelist(first_subj_file).name,Participant{1}.name,Participant{j}.name),Inverted_Conditions{1},Inverted_Conditions{i})])
            subjcondpair(end+1,:) = {Participant{j}.name,Inverted_Conditions{i}};
        end
        if any(any(isnan(D.DCM.H{1}))) %Debugging of NaNs
            disp(['NaNs found in file ' strrep(filelist(j).name,'DVT',Inverted_Conditions{i})])
            subjcondpair(end+1,:) = {Participant{j}.name,Inverted_Conditions{i}};
        end
    end
end