function [psc] = se_PSC_area(xSPM,targets)

global st;
global MAP;
global SPM;
global CLUSTER;
global displayType;
global group

% define number of sessions etc.
switch st.SPM;
    case 'SPM99'
        sess = size(xSPM.Sess,2);		%number of sessions
        
%        sps = size(xSPM.Sess{1}.row,2);		%scans per session
        sps = []; for xi=1:numel(SPM.Sess); sps(xi) = size(xSPM.Sess{xi}.row,2); end; 

        regressors = size(xSPM.Sess{1}.col,2);	%number of regressors per session (incl. movement regressors)
        allreg = (regressors);
        condSeq = [1:regressors];
    otherwise
        sess = size(xSPM.Sess,2);		%number of sessions
        
%        sps = size(xSPM.Sess(1).row,2);		%scans per session
        sps = []; for xi=1:numel(xSPM.Sess); sps(xi) = size(xSPM.Sess(xi).row,2); end;
        
        regressors = size(group(1).xSPM.Sess(1).Fc,2);	%number of regressors per session (incl. movement regressors)
        allreg = max(xSPM.Sess(1).col);
        condSeq =[];
        for i=1:regressors
            condSeq = [condSeq xSPM.Sess(1).Fc(i).i(1)];
        end
end

PSC = [];
    xyz = inv(xSPM.xVOL.M)*[targets;ones(1,size(targets,2))];

    betaAll  = [];
    for i = 1:length(xSPM.Vbeta)
    	betaAll = [betaAll spm_sample_vol(xSPM.Vbeta(i),xyz(1,:),xyz(2,:),xyz(3,:),0)];
    end
    betaAll = reshape(betaAll,size(targets,2),length(xSPM.Vbeta));
    dimBeta = size(betaAll,2);
for tar = 1:size(targets,2)
    beta = betaAll(tar,:);
    for c = 1:size(condSeq,2)
        k = condSeq(c);
    	seq = []; for i = 1:sess;  seq = [seq (i-1)*allreg+k]; 	end;  s = 0;
        
        cnt = 0;
	    for i = seq,
		    cnt = cnt+1; beta_xyz(c,cnt) = beta(i);
       		signChg(1:sps(cnt)) = beta_xyz(c,cnt)*xSPM.xX.xKXs.X(xSPM.Sess(cnt).row,i);
    		if beta_xyz(c,cnt)>0
			    delta(c,cnt) = max(signChg);	% positive signal change
    		else
			    delta(c,cnt) = min(signChg);	% negative signal change
    		end;
	    	sMSI(cnt) = beta(dimBeta-(sess-cnt));
	    	PSCtmp(c,cnt) = (delta(c,cnt)/sMSI(cnt))*100;     % k = konditions, s = sessions
	    end;
    end
    
    PSC = [PSC PSCtmp];
end


PSC = PSC(:,~isnan(sum(PSC)));

psc = struct('beta',beta_xyz,'delta',delta,'sMSI',sMSI,'PSC',PSC,'number',size(targets,2),'Sessions',size(condSeq,2));
       
        %-Reset title
%-----------------------------------------------------------------------
spm('Pointer','Arrow')

