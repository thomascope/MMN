function Inter_region(dirname_DCM,diagnosis_list,source_names,thresh)
%A script for plotting the results of extDCM across all diagnoses by
%inter-regional connection
addpath('/group/language/data/thomascope/MMN/ICA_denoise/Helperfiles')
thisdir = pwd;
mkdir([thisdir '/circuit_diagrams'])
cd([dirname_DCM 'PEB_secondlevel'])

load('PEB_A_Overall.mat')
template_PEB = PEB_Overall;

assert(all(template_PEB.M.X(:,1)==1),'The first column of the PEB of PEBs contrast should be all ones, check please.')

for this_contrast = 2:size(template_PEB.M.X,2)
    these_differences = find(BMA_Overall.Pp(:,this_contrast)>thresh);
    
    for this_difference = 1:length(these_differences)
        
        this_connection = BMA_Overall.Pnames{these_differences(this_difference)};
        Condition_Split = strsplit(this_connection,'Covariate ');
        condition = str2num(Condition_Split{2}(1));
        
        Connection_Split = strsplit(this_connection,'A{');
        direction = Connection_Split{2}(1);
        if strcmp(direction,'1')
            direction = 'forwards';
        elseif strcmp(direction,'2')
            direction = 'backwards';
        end
        to = str2num(Connection_Split{2}(4));
        from = str2num(Connection_Split{2}(6));
        
        
        if condition == 1
            if BMA_Overall.Ep(these_differences(this_difference),this_contrast)>0
                disp([diagnosis_list{find(template_PEB.M.X(:,this_contrast)==1)} ' stronger than ' diagnosis_list{find(template_PEB.M.X(:,this_contrast)==-1)} ' for ' direction ' ' source_names{from} ' to ' source_names{to}])
            else
                disp([diagnosis_list{find(template_PEB.M.X(:,this_contrast)==-1)} ' stronger than ' diagnosis_list{find(template_PEB.M.X(:,this_contrast)==1)} ' for ' direction ' ' source_names{from} ' to ' source_names{to}])
            end
        else
            disp(['Interaction between ' diagnosis_list{find(template_PEB.M.X(:,this_contrast)==-1)} ' and ' diagnosis_list{find(template_PEB.M.X(:,this_contrast)==1)} ' for ' direction ' ' source_names{from} ' to ' source_names{to}])
        end
        
    end
end
cd(thisdir)
