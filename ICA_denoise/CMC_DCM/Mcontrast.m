function [X, groups] = Mcontrast(Comparison, P)

% M contrast matrices
% Comparing HCs to each patient group when all are included
 clear X
 ss=0;
 
switch Comparison
    
    case 'All' % Compare all possible combinations within all groups
       
        for i = 1:length(P)
            
            if ~strcmp(P{1,i}.diag, 'MCI-no')
                ss = ss+1;
               
                X(ss,1) = 1; % first column is ones
                
                groups{ss,1} = P{1,i}.name;
                % separate groups out
                
                switch P{1,i}.diag
                    case 'Control'
                        groups{ss,2} = 1;
                        
                        % Comparisons: CvsPat, CvsBV, CvsPCA,CvsPPA,CvsMCI, 
                        %ADvsFTLD,'BVvsPCA','BVvsPPA','BVvsMCI', 'PCAvsPPA','PCAvsMCI','PPAvsMCI' 
                        X(ss, 2:13) = [ 1  1  1  1  1  0  0  0  0  0  0  0];
                    
                    case 'bvFTD'
                        groups{ss,2} = 2;
                        X(ss, 2:13) = [-1 -1  0  0  0 -1  1  1  1  0  0  0];
                        
                    case 'pca'
                        groups{ss,2} = 3;
                        X(ss, 2:13) = [-1  0 -1  0  0  1 -1  0  0  1  1  0];
                        
                    case 'nfvppa'
                        groups{ss,2} = 4;
                        X(ss, 2:13) = [-1  0  0 -1  0 -1  0 -1  0 -1  0  1];
                        
                    case 'MCI-yes'
                        groups{ss,2} = 5;
                        X(ss, 2:13) = [-1  0  0  0 -1  1  0  0 -1  0 -1 -1];
                end
                
                
            end
        end
    case 'HCvsAllPats'

        for i = 1:length(P)
            if ~strcmp(P{1,i}.diag, 'MCI-no')
                ss = ss+1;
                X(ss,1) = 1; % first column is ones
                
                groups{ss,1} = P{1,i}.name;
                % separate groups out
                
                switch P{1,i}.diag
                    case 'Control'
                        groups{ss,2} = 1;
                        
                        % Comparisons: CvsPat 
                        X(ss, 2) = [ 1];
                    
                    otherwise
                        groups{ss,2} = 2;
                        X(ss, 2) = [-1];
                end
            end
        end
        
    case 'HCvsBV'

        for i = 1:length(P)
            switch P{1,i}.diag
                case 'Control'
                    ss = ss+1;
                    X(ss,1) = 1; % first column is ones
                    
                    % Comparisons: CvsBV
                    X(ss, 2) = [ 1];
                    
                    groups{ss,1} = P{1,i}.name;
                    groups{ss,2} = 1;

                case 'bvFTD'
                    ss = ss+1;
                    X(ss,1) = 1; % first column is ones
                    
                    % Comparisons: CvsBV
                    X(ss, 2) = [-1];
                    
                    groups{ss,1} = P{1,i}.name;
                    groups{ss,2} = 2;
            end
        end
        
        
    case 'HCvsPCA'
         for i = 1:length(P)
            switch P{1,i}.diag
                case 'Control'
                    ss = ss+1;
                    X(ss,1) = 1; % first column is ones
                    
                    % Comparisons: CvsPCA
                    X(ss, 2) = [ 1];
                    
                    groups{ss,1} = P{1,i}.name;
                    groups{ss,2} = 1;

                case 'pca'
                    ss = ss+1;
                    X(ss,1) = 1; % first column is ones
                    
                    % Comparisons: CvsPCA
                    X(ss, 2) = [-1];
                    
                    groups{ss,1} = P{1,i}.name;
                    groups{ss,2} = 2;
            end
         end
        
    case 'HCvsPPA'
        for i = 1:length(P)
            switch P{1,i}.diag
                case 'Control'
                    ss = ss+1;
                    X(ss,1) = 1; % first column is ones
                    
                    % Comparisons: CvsnfvPPa
                    X(ss, 2) = [ 1];
                    
                    groups{ss,1} = P{1,i}.name;
                    groups{ss,2} = 1;

                case 'nfvppa'
                    ss = ss+1;
                    X(ss,1) = 1; % first column is ones
                    
                    % Comparisons: CvsnfvPPA
                    X(ss, 2) = [-1];
                    
                    groups{ss,1} = P{1,i}.name;
                    groups{ss,2} = 2;
            end
        end
        
    case 'HCvsMCI'
        for i = 1:length(P)
            switch P{1,i}.diag
                case 'Control'
                    ss = ss+1;
                    X(ss,1) = 1; % first column is ones
                    
                    % Comparisons: CvsPCA
                    X(ss, 2) = [ 1];
                    
                    groups{ss,1} = P{1,i}.name;
                    groups{ss,2} = 1;

                case 'MCI-yes'
                    ss = ss+1;
                    X(ss,1) = 1; % first column is ones
                    
                    % Comparisons: CvsPCA
                    X(ss, 2) = [-1];
                    
                    groups{ss,1} = P{1,i}.name;
                    groups{ss,2} = 2;
            end
        end
    case 'ADvsFTLD' % AD pathologies vs FTLD pathologies
        for i = 1:length(P)
            switch P{1,i}.diag
                case 'bvFTD'
                    ss = ss+1;
                    X(ss,1) = 1; % first column is ones
                    
                    % Comparisons: AD vs FTLD
                    X(ss, 2) = [ -1];
                    
                    groups{ss,1} = P{1,i}.name;
                    groups{ss,2} = 2;
                    
                case 'pca'
                    ss = ss+1;
                    X(ss,1) = 1; % first column is ones
                    
                    % Comparisons:AD vs FTLD
                    X(ss, 2) = [1];
                    
                    groups{ss,1} = P{1,i}.name;
                    groups{ss,2} = 1;
                case 'nfvppa'
                    ss = ss+1;
                    X(ss,1) = 1; % first column is ones
                    
                    % Comparisons:AD vs FTLD
                    X(ss, 2) = [-1];
                    
                    groups{ss,1} = P{1,i}.name;
                    groups{ss,2} = 2;
                    
                case 'MCI-yes'
                    ss = ss+1;
                    X(ss,1) = 1; % first column is ones
                    
                    % Comparisons: AD vs FTLD
                    X(ss, 2) = [1];
                    
                    groups{ss,1} = P{1,i}.name;
                    groups{ss,2} = 1;
            end
            
        end
        
    case 'BVvsPCA'
        for i = 1:length(P)
            switch P{1,i}.diag
                case 'bvFTD'
                    ss = ss+1;
                    X(ss,1) = 1; % first column is ones
                    
                    % Comparisons: BVvsPCA
                    X(ss, 2) = [ 1];
                    
                    groups{ss,1} = P{1,i}.name;
                    groups{ss,2} = 1;

                case 'pca'
                    ss = ss+1;
                    X(ss,1) = 1; % first column is ones
                    
                    % Comparisons: BVvsPCA
                    X(ss, 2) = [-1];
                    
                    groups{ss,1} = P{1,i}.name;
                    groups{ss,2} = 2;
            end
        end
    case 'BVvsPPA'
        for i = 1:length(P)
            switch P{1,i}.diag
                case 'bvFTD'
                    ss = ss+1;
                    X(ss,1) = 1; % first column is ones
                    
                    % Comparisons: BVvsPPA
                    X(ss, 2) = [ 1];
                    
                    groups{ss,1} = P{1,i}.name;
                    groups{ss,2} = 1;

                case 'nfvppa'
                    ss = ss+1;
                    X(ss,1) = 1; % first column is ones
                    
                    % Comparisons: BVvsPPA
                    X(ss, 2) = [-1];
                    
                    groups{ss,1} = P{1,i}.name;
                    groups{ss,2} = 2;
            end
        end
    case 'BVvsMCI'
        for i = 1:length(P)
            switch P{1,i}.diag
                case 'bvFTD'
                    ss = ss+1;
                    X(ss,1) = 1; % first column is ones
                    
                    % Comparisons: BVvsMCI
                    X(ss, 2) = [ 1];
                    
                    groups{ss,1} = P{1,i}.name;
                    groups{ss,2} = 1;

                case 'MCI-yes'
                    ss = ss+1;
                    X(ss,1) = 1; % first column is ones
                    
                    % Comparisons: BVvsMCI
                    X(ss, 2) = [-1];
                    
                    groups{ss,1} = P{1,i}.name;
                    groups{ss,2} = 2;
            end
        end
    case 'PCAvsPPA'
        for i = 1:length(P)
            switch P{1,i}.diag
                case 'pca'
                    ss = ss+1;
                    X(ss,1) = 1; % first column is ones
                    
                    % Comparisons: PCAvsPPA
                    X(ss, 2) = [ 1];
                    
                    groups{ss,1} = P{1,i}.name;
                    groups{ss,2} = 1;

                case 'nfvppa'
                    ss = ss+1;
                    X(ss,1) = 1; % first column is ones
                    
                    % Comparisons: PCAvsPPA
                    X(ss, 2) = [-1];
                    
                    groups{ss,1} = P{1,i}.name;
                    groups{ss,2} = 2;
            end
        end
    case 'PCAvsMCI'
         for i = 1:length(P)
            switch P{1,i}.diag
                case 'pca'
                    ss = ss+1;
                    X(ss,1) = 1; % first column is ones
                    
                    % Comparisons: PCAvsMCI
                    X(ss, 2) = [ 1];
                    
                    groups{ss,1} = P{1,i}.name;
                    groups{ss,2} = 1;

                case 'MCI-yes'
                    ss = ss+1;
                    X(ss,1) = 1; % first column is ones
                    
                    % Comparisons: PCAvsMCI
                    X(ss, 2) = [-1];
                    
                    groups{ss,1} = P{1,i}.name;
                    groups{ss,2} = 2;
            end
        end
    case 'PPAvsMCI'
         for i = 1:length(P)
            switch P{1,i}.diag
                case 'nfvppa'
                    ss = ss+1;
                    X(ss,1) = 1; % first column is ones
                    
                    % Comparisons: PPAvsMCI
                    X(ss, 2) = [ 1];
                    
                    groups{ss,1} = P{1,i}.name;
                    groups{ss,2} = 1;

                case 'MCI-yes'
                    ss = ss+1;
                    X(ss,1) = 1; % first column is ones
                    
                    % Comparisons: PPAvsMCI
                    X(ss, 2) = [-1];
                    
                    groups{ss,1} = P{1,i}.name;
                    groups{ss,2} = 2;
            end
         end
    otherwise
        disp('Comparision does not exist');
end
P = P;