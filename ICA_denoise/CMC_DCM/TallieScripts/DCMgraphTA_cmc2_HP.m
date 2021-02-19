function O = DCMgraphTA_cmc2_HP(BMA,varargin)


% setup:
V = varpairs2structTA(varargin{:});
if ~isfield(V,'plot_field'), V.plot_field = 'B'; end % 'A' for baseline trial input, 'A+B' for comparison trial input, 'B' for difference
if ~isfield(V,'DMn'), V.DMn = 1; end % design matrix column number
if ~isfield(V,'nodes'), V.nodes = {'left A1', 'left STG', 'left IPC', 'left IFG', 'right A1',  'right STG',  'right IPC', 'right IFG'}; end % {'L1' 'L2' 'L3' 'L4' 'R1' 'R2' 'R3' 'R4'};
if ~isfield(V,'pops'), V.pops = {'ss' 'sp' 'dp' 'ii'}; end
if ~isfield(V,'diffTF'), V.diffTF = 1; end
if ~isfield(V,'coldiffTF'), V.coldiffTF = 1; end
if ~isfield(V,'clim'), V.clim = []; end
if ~isfield(V,'conn_for'), V.conn_for = [2 1; 2 3]; end % i.e. sp->ss ; sp->dp
if ~isfield(V,'conn_bac'), V.conn_bac = [3 2; 3 4]; end % i.e. dp->sp ; dp->ii
if ~isfield(V,'conn_DA'), V.conn_DA = [1 1; 2 2; 3 3; 4 4]; end % i.e. dp->sp ; dp->ii
if ~isfield(V,'conn_for_type'), V.conn_for_type = [1 1]; end % i.e. the two conn_for are positive (use zero for modulatory)
if ~isfield(V,'conn_bac_type'), V.conn_bac_type = [-1 -1]; end % i.e. the two conn_bac are negative (use zero for modulatory)
if ~isfield(V,'conn_DA_type'), V.conn_DA_type = [0 0 0 0]; end % i.e. the two conn_bac are negative (use zero for modulatory)
if ~isfield(V,'manual_lats'), V.manual_lats = [2 6; 6 2; 3 7; 7 3; 4 8; 8 4]; end % leave empty if not using this
if ~isfield(V,'pop_order'), V.pop_order = [2 1 4 3]; end % i.e. sp ss dp ii (used after conn_for, conn_bac & conn_type)
if ~isfield(V,'node_order'), V.node_order = [1 2 4 3 5 6 8 7]; end
if ~isfield(V,'signifTF'), V.signifTF = 1; end
if ~isfield(V,'plotTF'), V.plotTF = 1; end
if ~isfield(V,'cmap'), V.cmap = []; end

nn = length(V.nodes);
np = length(V.pops);


% create node+pop labels and start ordering them:
order = [];
for k = 1:nn
    order = [order k:nn:(nn*np)];
end
lab = {};
for k = 1:np
    nodes_tmp = strcat(V.nodes,[' - ' V.pops{k}]);
    lab = [lab nodes_tmp{:}];
end
Lab = lab(order);

if ~V.coldiffTF, BMA.Ep(:,1) = 0; end


% locate values to plot
for kb = 1:length(V.DMn)
    s = strsplit(V.plot_field,'+');
    
    % AMPA:
    clear loca locan locm locc
    if V.signifTF
        HVAL = full(BMA.Pp(:,V.DMn(kb))) >=.05;
        hvala = {zeros(nn,nn) zeros(nn,nn) zeros(nn,nn) zeros(nn,nn)};
        hvalan = {zeros(nn,nn) zeros(nn,nn) zeros(nn,nn) zeros(nn,nn)};
        hvalm = zeros(nn,nn);
        hvalc = zeros(nn,1);
    end
    aval = {zeros(nn,nn) zeros(nn,nn) zeros(nn,nn) zeros(nn,nn)};
    anval = {zeros(nn,nn) zeros(nn,nn) zeros(nn,nn) zeros(nn,nn)};
    mval = zeros(nn,nn);
    cval = zeros(nn,1);
    for k = 1:length(BMA.Pnames)
        if strcmp(s{1},BMA.Pnames{k}(1)) && ~strcmp([s{1} 'N'],BMA.Pnames{k}(1:2)) % AMPA
            loca(k,:) = [str2num(BMA.Pnames{k}(3)) str2num(BMA.Pnames{k}(6)) str2num(BMA.Pnames{k}(8))];
            if V.signifTF, hvala{loca(k,1)}(loca(k,2),loca(k,3)) = HVAL(k); end
            aval{loca(k,1)}(loca(k,2),loca(k,3)) = BMA.Ep(k,1) + BMA.Ep(k,V.DMn(kb)); % adding the group of interest to the first column (i.e. group mean column)
        end
        if strcmp([s{1} 'N'],BMA.Pnames{k}(1:2)) % NMDA:
            locan(k,:) = [str2num(BMA.Pnames{k}(4)) str2num(BMA.Pnames{k}(7)) str2num(BMA.Pnames{k}(9))];
            if V.signifTF, hvalan{locan(k,1)}(locan(k,2),locan(k,3)) = HVAL(k); end
            anval{locan(k,1)}(locan(k,2),locan(k,3)) = BMA.Ep(k,1) + BMA.Ep(k,V.DMn(kb));
        end
        if strcmp('M',BMA.Pnames{k}(1)) % top-down DA
            locm(k,:) = [str2num(BMA.Pnames{k}(3)) str2num(BMA.Pnames{k}(5))];
            if V.signifTF, hvalm(locm(k,1),locm(k,2)) = HVAL(k); end
            mval(locm(k,1),locm(k,2)) = BMA.Ep(k,1) + BMA.Ep(k,V.DMn(kb));
        end
        try
            if strcmp('C',BMA.Pnames{k}(1)) % input
                locc(k) = str2num(BMA.Pnames{k}(3));
                if V.signifTF, hvalc(locc(k)) = HVAL(k); end
                cval(locc(k)) = BMA.Ep(k,1) + BMA.Ep(k,V.DMn);
            end
        catch
            warning('C matrix not used')
        end
        if length(s)==2
            for j = 1:4
                if strcmp(s{2},BMA.Pnames{k}(1)) && ~strcmp([s{2} 'N'],BMA.Pnames{k}(1:2))
                    loca(k,:) = [str2num(BMA.Pnames{k}(3)) str2num(BMA.Pnames{k}(6)) str2num(BMA.Pnames{k}(8))];
                    if V.signifTF, hvala{loca(k,1)}(loca(k,2),loca(k,3)) = hvala{loca(k,1)}(loca(k,2),loca(k,3)) || HVAL(k); end
                    aval{loca(k,1)+j-1}(loca(k,2),loca(k,3)) = aval{loca(k,1)+j-1}(loca(k,2),loca(k,3)) + BMA.Ep(k,V.DMn(kb));
                end
                if strcmp([s{2} 'N'],BMA.Pnames{k}(1:2))
                    locan(k,:) = [str2num(BMA.Pnames{k}(4)) str2num(BMA.Pnames{k}(7)) str2num(BMA.Pnames{k}(9))];
                    if V.signifTF, hvalan{locan(k,1)}(locan(k,2),locan(k,3)) = hvalan{locan(k,1)}(locan(k,2),locan(k,3)) || HVAL(k); end
                    anval{locan(k,1)}(locan(k,2),locan(k,3)) = anval{locan(k,1)}(locan(k,2),locan(k,3)) + BMA.Ep(k,V.DMn(kb));
                end
                if strcmp('M',BMA.Pnames{k}(1))
                    locm(k,:) = [str2num(BMA.Pnames{k}(3)) str2num(BMA.Pnames{k}(5))];
                    if V.signifTF, hvalm(locm(k,1),locm(k,2)) = hvalm(locm(k,1),locm(k,2)) || HVAL(k); end
                    mval(locm(k,1),locm(k,2)) = mval(locm(k,1),locm(k,2)) + BMA.Ep(k,V.DMn(kb));
                end
            end
        end
    end
    
    aval = cellfun(@(x) x(V.node_order,V.node_order),aval,'Uni',0);
    anval = cellfun(@(x) x(V.node_order,V.node_order),anval,'Uni',0);
    aval2 = aval;
    anval2 = anval;

    
    
    if strcmp('B',V.plot_field)
        aval{3} = triu(aval2{1});
        aval{4} = triu(aval2{1});
        aval{1} = tril(aval2{1});
        aval{2} = tril(aval2{1});
    elseif strcmp('A+B',V.plot_field)
        aval{3} = triu(aval2{3});
        aval{4} = triu(aval2{4});
        aval{1} = tril(aval2{1});
        aval{2} = tril(aval2{2});
    end

    if V.signifTF
        hvala = cellfun(@(x) x(V.node_order,V.node_order),hvala,'Uni',0);
         if strcmp('B',V.plot_field)
            hvala{3} = triu(hvala{1});
            hvala{4} = triu(hvala{1});
            hvala{1} = tril(hvala{1});
            hvala{2} = tril(hvala{1});
         elseif strcmp('A+B',V.plot_field)
            hvala{3} = triu(hvala{3});
            hvala{4} = triu(hvala{4});
            hvala{1} = tril(hvala{1});
            hvala{2} = tril(hvala{2});
         end
    end
    
    if ~isempty(V.manual_lats)
        for kk = 1:length(aval)
            for k = 1:size(V.manual_lats,1)
                if strcmp('B',V.plot_field)
                    aval{kk}(V.manual_lats(k,1),V.manual_lats(k,2)) = aval2{1}(V.manual_lats(k,1),V.manual_lats(k,2));
                else aval{kk}(V.manual_lats(k,1),V.manual_lats(k,2)) = aval2{kk}(V.manual_lats(k,1),V.manual_lats(k,2));
                end
            end
        end
        la = zeros(size(aval{1}));
        for k = 1:size(V.manual_lats,1)
            la(V.manual_lats(k,1),V.manual_lats(k,2)) = 1;
        end
        lan = zeros(size(anval{1}));
        for k = 1:size(V.manual_lats,1)
            lan(V.manual_lats(k,1),V.manual_lats(k,2)) = 1;
        end
    else
        la = cellfun(@logical,aval2,'Uni',0);
        la = la{1} & la{1}';
        lan = cellfun(@logical,anval2,'Uni',0);
        lan = lan{1} & lan{1}';
    end
    la2 = repelem(la,6,6);
    lan2 = repelem(lan,6,6);
    
    % fill matrix with the above located values:
    CG1 = repmat({nan(nn*np,nn*np)},1,5);
%     if V.signifTF; C(:) = cval .* hvalc; else C(:) = cval; end % inputs
    for kc = 1:size(V.conn_for,1) % forward AMPA connections
        c = V.conn_for(kc,:)-1;
        for ki = 1:nn
            for kj = 1:nn
                if V.signifTF
                    CG1{1}((ki*np)-np+1+c(2),(kj*np)-np+1+c(1)) = aval{1+kc-1}(ki,kj) .* hvala{1+kc-1}(ki,kj);
                else CG1{1}((ki*np)-np+1+c(2),(kj*np)-np+1+c(1)) = aval{1+kc-1}(ki,kj);
                end
            end
        end
    end
    for kc = 1:size(V.conn_bac,1) % backward AMPA connections
        c = V.conn_bac(kc,:)-1;
        for ki = 1:nn
            for kj = 1:nn
                if V.signifTF
                    CG1{2}((ki*np)-np+1+c(2),(kj*np)-np+1+c(1)) = aval{3+kc-1}(ki,kj) .* hvala{3+kc-1}(ki,kj);
                else CG1{2}((ki*np)-np+1+c(2),(kj*np)-np+1+c(1)) = aval{3+kc-1}(ki,kj);
                end
            end
        end
    end
    for kc = 1:size(V.conn_DA,1) % top-down-DA connections
        c = V.conn_DA(kc,:)-1;
        for ki = 1:nn
            for kj = 1:nn
                if V.signifTF
                    CG1{5}((ki*np)-np+1+c(2),(kj*np)-np+1+c(1)) = mval(ki,kj) .* hvalm(ki,kj);
                else CG1{5}((ki*np)-np+1+c(2),(kj*np)-np+1+c(1)) = mval(ki,kj);
                end
            end
        end
    end
    
    %CG1{1} = nanmax(cat(3,CG1{1},CG1{2}),[],3);
    %CG1{3} = nanmax(cat(3,CG1{3},CG1{4}),[],3);

    try
        C(C==0) = NaN;
        C = exp(C);
        Ccat{kb} = C;
        Ccat{kb}(isnan(Ccat{kb})) = 0;
    end

    ct = zeros(1,np);
    ctf1 = ct;
    ctf1(V.conn_for(:,1)) = V.conn_for_type;
    ctb1 = ct;
    ctb1(V.conn_bac(:,1)) = V.conn_bac_type;
    ctl = ct;
    ctl(V.conn_DA(:,1)) = V.conn_DA_type;
    hct = {[ones(nn*np,np).*ctf1] [ones(nn*np,np).*ctb1]};% ones(nn*np,np).*ctl};
    KK = {V.conn_for(:,1) V.conn_bac(:,1)};% V.conn_DA(:,1)};
    tit = {'Forward' 'Backward'};% 'top-down-DA'};
    for K = 1:length(KK) % for backward on one plot and forward on the next, etc
        if ~isempty(KK{K})
            % re-order for display. Also: sp ss dp ii i.e. swap the first two populations for display purposes - makes more intuative sense
            % in second dimension:
            CG2 = [];
            head = [];
            Lab2 = {};
            CG2tmp = [];
            for ki = 1:nn
                tmp = CG1{K}(:,V.pop_order+(ki*np)-np);
                i = 1:size(tmp,2); 
                i(V.pop_order(KK{K})) = [];
                tmp(:,i) = 0;
                CG2tmp = [CG2tmp tmp zeros(nn*np,2)];
                Lab2 = [Lab2 Lab(V.pop_order+(ki*np)-np) {'' ''}];
                head = [head hct{K}(:,V.pop_order) zeros(nn*np,2)];
            end
            CG2 = cat(3,CG2,CG2tmp);
            % in first dimension:
            CG3 = [];
            head2 = [];
            CG3tmp = [];
            for ki = 1:nn
                tmp = CG2(V.pop_order+(ki*np)-np,:);
                CG3tmp = [CG3tmp ; tmp ; zeros(2,((nn*np)+(nn*2)))];
                head2 = [head2 ; head(V.pop_order+(ki*np)-np,[1:((nn*np)+(nn*2))]) ; zeros(2,((nn*np)+(nn*2)))];
            end
            if min(min(head2))==-1
                head2 = head2 + (head2 & la2);
            else head2 = head2 - (head2 & la2);
            end
            CG3 = cat(3,CG3,CG3tmp);

            CG = CG3;
            CG(CG3==-32) = NaN;
            CG(CG3==0) = NaN;
            CG = nanmean(CG,3);
            CG = exp(CG);
            CG(isnan(CG)) = 0;

            NodeClrs = zeros(((nn*np)+(nn*2)),3);
            NodeClrs([(np+1):(np+2):end np+2:np+2:end],:) = 1;

            % flip & rotate:
            order2 = [1:(((nn*np)+(nn*2))/2) fliplr(((((nn*np)+(nn*2))/2)+1):((nn*np)+(nn*2)))];
            CG = CG(order2,order2);
            head2 = head2(order2,order2);
            NodeClrs = NodeClrs(order2,:);
            Lab2 = Lab2(order2);
            if isint((((nn*np)+(nn*2))/4))
                order3 = [(((nn*np)+(nn*2))/4):((nn*np)+(nn*2)) ones(1,np)*(np+1) 1:((((nn*np)+(nn*2))/4)-1)];
            else order3 = ceil([(((nn*np)+(nn*2))/4):((nn*np)+(nn*2)) ones(1,np)*(np+1) 1:((((nn*np)+(nn*2))/4))]);
            end
            CG = CG(order3,order3);
            head2 = head2(order3,order3);
            NodeClrs = NodeClrs(order3,:);
            Lab2 = Lab2(order3);

            H{K} = head2;
            CGcat{kb}{K} = CG;
        else CGcat{kb}{K} = nan(size(CGcat{kb}{1}));
            H{K} = nan(size(H{1}));
        end
    end

end


% finally, plot:
if isempty(V.clim)
    cg_tmp = [];
    for K = 1:length(CGcat{1})
        if length(CGcat)==2
            cg_tmp = [cg_tmp(:) ; CGcat{1}{K}(:)-CGcat{2}{K}(:)];
        else cg_tmp = [cg_tmp(:) ; CGcat{1}{K}(:)];
        end
    end
%     if length(CGcat)==2
%         cg_tmp = [cg_tmp(:) ; Ccat{1}(:)-Ccat{2}(:)];
%     else cg_tmp = [cg_tmp(:) ; Ccat{1}(:)];
%     end
    V.clim = [nanmin(cg_tmp) nanmax(cg_tmp)];
end
if V.plotTF, figure, set(gcf,'Units','norm','Position',[0 0 1 1]), end
for K = 1:length(KK)
    if length(CGcat)==2
        z = CGcat{1}{K} - CGcat{2}{K};
    else z = CGcat{1}{K};
    end
    if any(z(:))
        z(np+1,np+1) = V.clim(1); % arbitrarily set two values that you know will not be seen as the max and min so that it scales to these
        z(np+2,np+2) = V.clim(2);

        if V.plotTF
            subplot(2,2,K)
            title(tit{K})
            if length(CGcat)==2
                for kz = 1:2
                    ztmp = z;
                    scaler = 100;
                    if kz==1
                        ztmp(ztmp<0) = 0;
                        cm = flipud(cbrewer('seq','YlOrRd',max([ceil(max(ztmp(:))*scaler)+1 5])));
                    else, ztmp(ztmp>0) = 0;
                        ztmp = -ztmp;
                        cm = flipud(cbrewer('seq','Blues',max([floor(max(ztmp(:))*scaler)+1 5])));
                    end
                    if ~isempty(V.cmap), cm = V.cmap; end
                    hc = circularGraphArr(ztmp,'Label',Lab2,'Colormap',NodeClrs,'Head',H{K},'ScaledClrTF',1,'ConnOrNodeClr','Conn','ConnClrs',cm,'scaler',scaler);
                end
            else
                if V.clim(2) < 100
                    scaler = 100;
                else, scaler = 10;
                end
                if V.diffTF
                    cm = [flipud(winter(100)) ; [.7 .7 .7] ; autumn(ceil(V.clim(2)*scaler)-100)];
                else, cm = flipud(autumn(ceil(V.clim(2)*scaler)+1)); cm(1,:) = [.7 .7 .7];
                end
                if ~isempty(V.cmap), cm = V.cmap; end
                hc = circularGraphArr(z,'Label',Lab2,'Colormap',NodeClrs,'Head',H{K},'ScaledClrTF',1,'ConnOrNodeClr','Conn','ConnClrs',cm,'scaler',scaler);
            end
            try % inputs
                x = linspace(-pi,pi,360);
                deg = (360-(360/nn/4))/nn;
                for k = 1:nn
                    st = round((deg*k)-(deg/3)-(deg/4));
                    en = round((deg*k)+(deg/3)-(deg/4));
                    xx{k} = x(st:en);
                end
                for k = 1:nn
                    if length(Ccat)==2
                        cval = ceil((Ccat{1}(k)-Ccat{2}(k))*100);
                        if cval~=0
                            for kz = 1:2
                                if kz==1
                                    cm = colormap_tor([1 1 1],[0 1 1],[0 0 .5],[0 0 1],[0 .5 1]);
                                    cval(cval<0) = 0;
                                    if cval~=0, plot(sin(xx{k})*1.1,cos(xx{k})*1.1,'-k','LineWidth',3,'Color',cm(cval,:)), end
                                else, cm = colormap_tor([1 1 1],[1 1 0],[0 0 .5],[1 0 0],[1 .5 0]);
                                    cval(cval>0) = 0;
                                    cval = -cval;
                                    if cval~=0, plot(sin(xx{k})*1.1,cos(xx{k})*1.1,'-k','LineWidth',3,'Color',cm(cval,:)), end
                                end
                            end
                        end
                    else cval = ceil(Ccat{1}(k)*100);
                        if cval~=0, plot(sin(xx{k})*1.1,cos(xx{k})*1.1,'-k','LineWidth',3,'Color',cm(cval,:)), end
                    end
                end
            end
            drawnow
        end
    else
        disp(['none in ' tit{K}])
    end
end

O.CGcat = CGcat;
%O.Ccat = Ccat;
O.Lab2 = Lab2;
O.H = H;
if exist('cm','var')
    O.cmap = cm;
else O.cmap = V.cmap;
end


end

