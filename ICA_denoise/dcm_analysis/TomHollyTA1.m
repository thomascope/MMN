% A script by Tallie to combine the A+B matrices for DCM connection
% strengths

%% LABELS ACCORDING TO THE SUBJECT DATA MATRIX (collated with the help of ConStrength_A_README.m and ConStrength_B_README.m:

Alab_for1 = {'LA1-LSTG SP->SS' 'RA1-RSTG SP->SS' 'LSTG-LIFG SP->SS' 'RSTG-RIFG SP->SS' 'LSTG-LIPC SP->SS' 'RSTG-RIPC SP->SS' 'LIPC-LIFG SP->SS' 'RIPC-RIFG SP->SS'};
Alab_for2 = {'LA1-LSTG SP->DP' 'RA1-RSTG SP->DP' 'LSTG-LIFG SP->DP' 'RSTG-RIFG SP->DP' 'LSTG-LIPC SP->DP' 'RSTG-RIPC SP->DP' 'LIPC-LIFG SP->DP' 'RIPC-RIFG SP->DP'};
Alab_bac1 = {'LSTG-LA1 DP->SP' 'RSTG-RA1 DP->SP' 'LIFG-LSTG DP->SP' 'RIFG-RSTG DP->SP' 'LIPC-LSTG DP->SP' 'RIPC-RSTG DP->SP' 'LIFG-LIPC DP->SP' 'RIFG-RIPC DP->SP'};
Alab_bac2 = {'LSTG-LA1 DP->II' 'RSTG-RA1 DP->II' 'LIFG-LSTG DP->II' 'RIFG-RSTG DP->II' 'LIPC-LSTG DP->II' 'RIPC-RSTG DP->II' 'LIFG-LIPC DP->II' 'RIFG-RIPC DP->II'};
Alab_lat1 = {'RSTG-LSTG SP->SS' 'LSTG-RSTG SP->SS' 'RIFG-LIFG SP->SS' 'LIFG-RIFG SP->SS' 'RIPC-LIPC SP->SS' 'LIPC-RIPC SP->SS'};
Alab_lat2 = {'RSTG-LSTG SP->DP' 'LSTG-RSTG SP->DP' 'RIFG-LIFG SP->DP' 'LIFG-RIFG SP->DP' 'RIPC-LIPC SP->DP' 'LIPC-RIPC SP->DP'};
Alab_lat3 = {'RSTG-LSTG DP->SP' 'LSTG-RSTG DP->SP' 'RIFG-LIFG DP->SP' 'LIFG-RIFG DP->SP' 'RIPC-LIPC DP->SP' 'LIPC-RIPC DP->SP'};
Alab_lat4 = {'RSTG-LSTG DP->II' 'LSTG-RSTG DP->II' 'RIFG-LIFG DP->II' 'LIFG-RIFG DP->II' 'RIPC-LIPC DP->II' 'LIPC-RIPC DP->II'};

Blab = {'LA1-LSTG' 'RA1-RSTG' 'LSTG-LIFG' 'RSTG-RIFG' 'LSTG-LIPC' 'RSTG-RIPC' 'LIPC-LIFG' 'RIPC-RIFG',...
    'LSTG-LA1' 'RSTG-RA1' 'LIFG-LSTG' 'RIFG-RSTG' 'LIPC-LSTG' 'RIPC-RSTG' 'LIFG-LIPC' 'RIFG-RIPC',...
    'RSTG-LSTG' 'LSTG-RSTG' 'RIFG-LIFG' 'LIFG-RIFG' 'RIPC-LIPC' 'LIPC-RIPC'};


%% THE CONNECTION MATRIX:

%          1      2      3      4      5      6      7      8
order = {'LA1' 'LSTG' 'LIFG' 'LIPC' 'RA1' 'RSTG' 'RIFG' 'RIPC'};
% 1 LA1   [0      0      0      0      0      0      0      0
% 2 LSTG   0      0      0      0      0      0      0      0
% 3 LIFG   0      0      0      0      0      0      0      0
% 4 LIPC   0      0      0      0      0      0      0      0
% 5 RA1    0      0      0      0      0      0      0      0
% 6 RSTG   0      0      0      0      0      0      0      0
% 7 RIFG   0      0      0      0      0      0      0      0
% 8 RIPC   0      0      0      0      0      0      0      0];


%% load the DATA MATRIX:

% first cd to the relevant DCM_params folder, e.g.:
path_dcm = '/imaging/hp02/pnfa_mmn/dcm/CMC_DCM_all_subjs_together_camcanHCs/2019_TC_LFPs_customPriors_32models/';
datafit = {'std_dvt', 'std_loc','std_int','std_dur','std_gap','std_frq'};

for dvt = 1:length(datafit)
    
    % then load the relevant A and B data matrices for all subjects, e.g.:
    %load('con_strength_std_frq_A.mat')
    %load('con_strength_std_frq_B.mat')
    load(fullfile(path_dcm, datafit{dvt}, 'DCM_params', sprintf('con_strength_%s_A.mat', datafit{dvt})));
    load(fullfile(path_dcm, datafit{dvt}, 'DCM_params', sprintf('con_strength_%s_B.mat', datafit{dvt})));
    
    %% find where each column of the DATA MATRIX should go in the CONNECTION MATRIX :
    
    % first find the index from the label list above:
    clear Bi Ai_for1 Ai_for2 Ai_bac1 Ai_bac2 Ai_lat1 Ai_lat2 Ai_lat3 Ai_lat4
    for i = 1:length(order)
        for j = 1:length(order)
            
            Bi(i,j) = max([find(cellfun(@(x) ~isempty(regexp(x,[order{i} '-' order{j}],'once')),Blab)) 0]);
            
            Ai_for1(i,j) = max([find(cellfun(@(x) ~isempty(regexp(x,[order{i} '-' order{j}],'once')),Alab_for1)) 0]);
            Ai_for2(i,j) = max([find(cellfun(@(x) ~isempty(regexp(x,[order{i} '-' order{j}],'once')),Alab_for2)) 0]);
            Ai_bac1(i,j) = max([find(cellfun(@(x) ~isempty(regexp(x,[order{i} '-' order{j}],'once')),Alab_bac1)) 0]);
            Ai_bac2(i,j) = max([find(cellfun(@(x) ~isempty(regexp(x,[order{i} '-' order{j}],'once')),Alab_bac2)) 0]);
            
            Ai_lat1(i,j) = max([find(cellfun(@(x) ~isempty(regexp(x,[order{i} '-' order{j}],'once')),Alab_lat1)) 0]);
            Ai_lat2(i,j) = max([find(cellfun(@(x) ~isempty(regexp(x,[order{i} '-' order{j}],'once')),Alab_lat2)) 0]);
            Ai_lat3(i,j) = max([find(cellfun(@(x) ~isempty(regexp(x,[order{i} '-' order{j}],'once')),Alab_lat3)) 0]);
            Ai_lat4(i,j) = max([find(cellfun(@(x) ~isempty(regexp(x,[order{i} '-' order{j}],'once')),Alab_lat4)) 0]);
            
        end
    end
    
    % then add the starting index from the column ordering given in ConStrength_A_README.m and ConStrength_B_README.m:
    Bi(Bi>0) = Bi(Bi>0) + 2;
    Ai_for1(Ai_for1>0) = Ai_for1(Ai_for1>0) + 2;
    Ai_for2(Ai_for2>0) = Ai_for2(Ai_for2>0) + 10;
    Ai_bac1(Ai_bac1>0) = Ai_bac1(Ai_bac1>0) + 18;
    Ai_bac2(Ai_bac2>0) = Ai_bac2(Ai_bac2>0) + 26;
    Ai_lat1(Ai_lat1>0) = Ai_lat1(Ai_lat1>0) + 34;
    Ai_lat2(Ai_lat2>0) = Ai_lat2(Ai_lat2>0) + 40;
    Ai_lat3(Ai_lat3>0) = Ai_lat3(Ai_lat3>0) + 46;
    Ai_lat4(Ai_lat4>0) = Ai_lat4(Ai_lat4>0) + 52;
    
    
    %% put each subjects std-trial parameter value into a connection matrix (region x region x subject)
    
    clear A_for1 A_for2 A_bac1 A_bac2 A_lat1 A_lat2 A_lat3 A_lat4
    for k = 1:size(connectivity_EpA,1)
        
        tmp = Ai_for1;
        tmp(tmp>0) = connectivity_EpA(k,tmp(tmp>0));
        A_for1(:,:,k) = tmp;
        
        tmp = Ai_for2;
        tmp(tmp>0) = connectivity_EpA(k,tmp(tmp>0));
        A_for2(:,:,k) = tmp;
        
        tmp = Ai_bac1;
        tmp(tmp>0) = connectivity_EpA(k,tmp(tmp>0));
        A_bac1(:,:,k) = tmp;
        
        tmp = Ai_bac2;
        tmp(tmp>0) = connectivity_EpA(k,tmp(tmp>0));
        A_bac2(:,:,k) = tmp;
        
        tmp = Ai_lat1;
        tmp(tmp>0) = connectivity_EpA(k,tmp(tmp>0));
        A_lat1(:,:,k) = tmp;
        
        tmp = Ai_lat2;
        tmp(tmp>0) = connectivity_EpA(k,tmp(tmp>0));
        A_lat2(:,:,k) = tmp;
        
        tmp = Ai_lat3;
        tmp(tmp>0) = connectivity_EpA(k,tmp(tmp>0));
        A_lat3(:,:,k) = tmp;
        
        tmp = Ai_lat4;
        tmp(tmp>0) = connectivity_EpA(k,tmp(tmp>0));
        A_lat4(:,:,k) = tmp;
        
    end
    
    
    %% put each subjects difference-trial parameter value into a connection matrix (region x region x subject)
    
    clear B_all
    for k = 1:size(connectivity_EpB,1)
        tmp = Bi;
        tmp(tmp>0) = connectivity_EpB(k,tmp(tmp>0));
        B_all(:,:,k) = tmp;
    end
    
    %% combine the results of the above two sections to get the dvt-trial  parameter value into a connection matrix (region x region x subject)
    
    AB_for1 = zeros(size(A_for1));
    AB_for1(A_for1>0) = A_for1(A_for1>0) .* B_all(A_for1>0);
    
    AB_for2 = zeros(size(A_for2));
    AB_for2(A_for2>0) = A_for2(A_for2>0) .* B_all(A_for2>0);
    
    AB_bac1 = zeros(size(A_bac1));
    AB_bac1(A_bac1>0) = A_bac1(A_bac1>0) .* B_all(A_bac1>0);
    
    AB_bac2 = zeros(size(A_bac2));
    AB_bac2(A_bac2>0) = A_bac2(A_bac2>0) .* B_all(A_bac2>0);
    
    AB_lat1 = zeros(size(A_lat1));
    AB_lat1(A_lat1>0) = A_lat1(A_lat1>0) .* B_all(A_lat1>0);
    
    AB_lat2 = zeros(size(A_lat2));
    AB_lat2(A_lat2>0) = A_lat2(A_lat2>0) .* B_all(A_lat2>0);
    
    AB_lat3 = zeros(size(A_lat3));
    AB_lat3(A_lat3>0) = A_lat3(A_lat3>0) .* B_all(A_lat3>0);
    
    AB_lat4 = zeros(size(A_lat4));
    AB_lat4(A_lat4>0) = A_lat4(A_lat4>0) .* B_all(A_lat4>0);
    
    
    %% example check plot:
    
    figure; imagesc(mean(AB_for1,3)), colorbar
    title('mean dvt-trial extrinsic forward connectivity (SP->SS)')
    
    
    %% Added section for statistical comparison
    all_matrices = {'A_for1','A_for2','A_bac1','A_bac2','A_lat1','A_lat2','A_lat3','A_lat4','B_all','AB_for1','AB_for2','AB_bac1','AB_bac2','AB_lat1','AB_lat2','AB_lat3','AB_lat4'};
    all_matrices_labels = {'Alab_for1','Alab_for2','Alab_bac1','Alab_bac2','Alab_lat1','Alab_lat2','Alab_lat3','Alab_lat4','Blab','Alab_for1','Alab_for2','Alab_bac1','Alab_bac2','Alab_lat1','Alab_lat2','Alab_lat3','Alab_lat4'};
    
    
    % disp(['***' num2str(dvt) '***'])
    p = zeros(length(datafit),length(all_matrices),8,8);
    for this_matrix = 1:length(all_matrices)
        for from = 1:8
            for to = 1:8
                
                
                [p(dvt,this_matrix,from,to),~,~] = eval(['anova1(log(squeeze(' all_matrices{this_matrix} '(' num2str(from) ',' num2str(to) ',:))),connectivity_EpA(:,2),''off'');']);
                if p(dvt,this_matrix,from,to)<0.01
                    disp([datafit{dvt} ' from ' order{from} ' to ' order{to} ' ' all_matrices{this_matrix} ' p = ' num2str(p(dvt,this_matrix,from,to))])
                    %disp(eval([all_matrices_labels{this_matrix} '(' num2str(sum(sum(~isnan(p(dvt,this_matrix,:,:)))) - sum(sum(p(dvt,this_matrix,:,:)==0)) + 1) ')']));
                    [p(dvt,this_matrix,from,to),~,~] = eval(['anova1(log(squeeze(' all_matrices{this_matrix} '(' num2str(from) ',' num2str(to) ',:))),connectivity_EpA(:,2),''on'');']);
                    figure
                    eval(['scatter(connectivity_EpA(:,2),log(squeeze(' all_matrices{this_matrix} '(' num2str(from) ',' num2str(to) ',:))))']);
                    xlim([0 6])
                    hold on
                    plot([0,6],[0,0],'k--')
                    group_mean = zeros(1,5);
                    group_ste = zeros(1,5);
                    for grp = 1:max(connectivity_EpA(:,2))
                        eval(['scatter(grp,mean(log(squeeze(' all_matrices{this_matrix} '(' num2str(from) ',' num2str(to) ',connectivity_EpA(:,2)==grp)))),48,''r'',''filled'')']);
                        eval(['group_mean(grp) = mean(log(squeeze(' all_matrices{this_matrix} '(' num2str(from) ',' num2str(to) ',connectivity_EpA(:,2)==grp))));']);
                        eval(['group_ste(grp) = std(log(squeeze(' all_matrices{this_matrix} '(' num2str(from) ',' num2str(to) ',connectivity_EpA(:,2)==grp))))/sqrt(sum(connectivity_EpA(:,2)==grp));']);
                    end
                    figure
                    errorbar([1:5],group_mean,group_ste,'kx')
                    set(gca,'XTick',[1:5])
                    set(gca,'XTickLabel',{'Control','bvFTD','PCA','nfvPPA','MCI'})
                    xlim([0 6])
                    hold on
                    plot([0,6],[0,0],'k--')
                    pause
                end
                close all
                
            end
        end
    end
end
% end

