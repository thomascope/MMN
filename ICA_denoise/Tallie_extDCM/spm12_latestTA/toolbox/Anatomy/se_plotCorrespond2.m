function [PQ, Z, nP, index] = se_plotCorrespond(CLUSTER)

global st;
global MAP;
global SPM;
global xSPM;

FS   = spm('FontSizes');                 PF   = spm_platform('fonts');
side = sign(mean(CLUSTER.XYZmm(1,:))); 
if side == 0, side = sign(median(CLUSTER.XYZmm(1,:))); end
if side == 0, side = 1; end
if side == -1; seite = 'Left'; else, seite = 'Right';end;


try 
    PQ      = CLUSTER.PQ;
    Z       = CLUSTER.xZ;
    index   = CLUSTER.index;
catch
    PQ    = [];
    Z     = [];
    index = [];
end

if isempty(PQ) || isempty(Z) || isempty(index)
    for i=1:size(MAP,2)
        Z(i,:)  = (spm_get_data(MAP(i).PMap,CLUSTER.gtmp));
        if any(Z(i,:)>0)
            PQ(i)   = mean(Z(i,Z(i,:)>0))/MAP(i).AM(side+2);
            PV(i)   = sum(Z(i,:)>0);
            PP(i)   = max(Z(i,Z(i,:)>0));
        else
            PV(i)   = 0;
            PQ(i)   = 0;
            PP(i)   = 0;
        end
    end
    Pv = PV/size(Z,2)*100;
    Areas = {}; PQ(~isreal(PQ)) = 0;
    [buffer index] = sort(PQ); index = fliplr(index);
    cut = min(numel(find(PQ>.25 & Pv > 10)),10);
    index = index(1:cut); Z = Z(index,:); PQ = PQ(index);
end




if any(PQ)
    
    try; delete(st.FX);             end
    if 1
        FS   = spm('FontSizes');                 PF   = spm_platform('fonts');
        screen = get(0,'ScreenSize');
        st.FX = figure(...
        'Tag','SATB',                             'Position',[10 screen(4)-590 screen(3)*.35 500],...
        'Resize','on',                            'MenuBar','figure',...
        'Color','w',                              'ColorMap',gray(64),...
        'DefaultTextColor','k',                   'DefaultTextInterpreter','tex',...
        'DefaultTextFontName',PF.helvetica,       'DefaultTextFontSize',FS(12),...
        'DefaultAxesColor','w',                   'DefaultAxesXColor','k',...
        'DefaultAxesYColor','k',                  'DefaultAxesZColor','k',...
        'DefaultAxesFontName',PF.helvetica,       'DefaultPatchFaceColor','k',...
        'DefaultPatchEdgeColor','k',              'DefaultSurfaceEdgeColor','k',...
        'DefaultLineColor','k',                   'DefaultUicontrolFontName',PF.helvetica,...
        'DefaultUicontrolFontSize',FS(12),        'DefaultUicontrolInterruptible','on',...
        'PaperType','A4',                         'PaperUnits','normalized',...
        'PaperPosition',[.0726 .0644 .854 .870],  'InvertHardcopy','off',...
        'Renderer','zbuffer',                     'Visible','on');

        set(st.FX,'Name',[spm('ver') ' Anatomy Toolbox'],'NumberTitle','off')
    end
    [B,INDEX] = sort(CLUSTER.maxZ); INDEX = fliplr(INDEX);
    if numel(INDEX)>8, INDEX = INDEX(1:8); end
    xyz = inv(MAP(1).MaxMap.mat) * xSPM.xVOL.M * [CLUSTER.maxM(1,INDEX); CLUSTER.maxM(2,INDEX); CLUSTER.maxM(3,INDEX); ones(1,size(INDEX,2))];		%-Reorientation matrix
    T   = CLUSTER.maxZ(1,INDEX);
    Area = {MAP(index).name};


    try
        nP  = CLUSTER.nP;
        nFx = CLUSTER.nFx;
    catch
        nP  = [];
        nFx = [];
    end
    wicht = [.1:.1:1];
    if isempty(nP), do = 1; else; do = 0; end
    for i=1:numel(index)
        nA = MAP(index(i)).nA{side+2};
        if do
            [nF,xout] = histc(Z(i,Z(i,:)>0),[0.01:.1:1.01]); nF=nF(1:end)/sum(nF)*100; nF = nF(1:end-1);
            nP(i,:)   = zeros(1,10)*NaN; nP(i,nA>0) = (nF(nA>0)-nA(nA>0))./nA(nA>0);
            nP2(i,:)   = zeros(1,10)*NaN; nP2(i,nA>0) = (nF(nA>0)./nA(nA>0));
        end
        nA(nA == 0) = NaN; nAp(i,:)  = nA;
        nG(i) = sum(nP2(i,~isnan(nP2(i,:))) .* wicht(~isnan(nP2(i,:)))  )/(sum(nP2(i,~isnan(nP2(i,:))))) * 10;
    end

    nPplot = nP;  scal = [min(min(nP(~isnan(nP)))) max(max(nP(~isnan(nP))))];    
    nPplot(nP<0) = round((nP(nP<0)-scal(1)/(0-scal(1)))*100)+1;;
    nPplot(nP>0) = round((nP(nP>0)/(scal(2)))*100)+100;;
    nPplot(isnan(nP)) = 301;
    
    

    
    ax=axes('Parent',st.FX,'units','normalized', 'Position',[0.16, 0.075, .63, .8], 'Visible','off');
    image(nPplot), colormap(MAP(1).cmap)
    set(gca,'YTick',1:numel(index),'XTick',[1:10]+.5,'YTickLabel',Area,'XTickLabel',[10:10:100],'FontSize',FS(9),'FontWeight','bold','TickLength',[0 0])

    ax=axes('Parent',st.FX,'units','normalized', 'Position',[0.16, .9, .63, .05], 'Visible','off');
    image(repmat(linspace(1,200,1000),20,1)), colormap(MAP(1).cmap)
    set(gca,'XTick',[1 500 1000],'XTickLabel',{[int2str(scal(1)*100) '%'];'observed vs. expected probability represent.';['+' int2str(scal(2)*100) '%']},'YTick',[],'XAxisLocation','top')
    
    ax=axes('Parent',st.FX,'units','normalized', 'Position',[0.83, .075, .14, .8], 'Visible','off','XLim',[-1 1],'YLim',[.5 numel(index)+.5],'YDir','reverse');
    for i=1:numel(index)
        
        AN = MAP(index(i)).nA{side+2}.*[.1:.1:1];
        text(-.6, i,[num2str(mean(Z(i,Z(i,:)>0))*100,'%3.1f') '%'],'Color','k','FontSize',FS(11),'HorizontalAlignment','center','VerticalAlignment','bottom')
        text(-.6, i, [num2str(sum(AN(~isnan(AN))),'%3.1f') '%'],'Color','k','FontSize',FS(11),'HorizontalAlignment','center','VerticalAlignment','top')
        line([-1 -.3],[i i],'Color','k')
        
        text(.45, i, ['= ' num2str(mean(Z(i,Z(i,:)>0))*100/sum(AN(~isnan(AN))),'%3.2f')],'Color','k','FontSize',FS(11),'HorizontalAlignment','center','VerticalAlignment','middle','FontWeight','bold')
    end


%     maxima = zeros(numel(T),10,numel(index));
%     exceed = ones(numel(T),10,numel(index))*NaN;
% 	for i=1:numel(index)
%         p = round(spm_get_data(MAP(index(i)).PMap,xyz)/25); cnt = 0;
%         for ii=1:min(numel(T),5)
%             if p(ii) > 0
%                 maxima(ii,p(ii),i) = 1;
%                 exceed(ii,p(ii),i) = (sum(MAP(index(i)).allZ(MAP(index(i)).allLR == side)>p(ii)))  /sum(MAP(index(i)).allLR == side)*100;
%            end
%         end
%     end
%     
%     
%     indx2 = find(squeeze(sum(sum(maxima,1),2))>0);
% 
%     if numel(indx2)>0
%     maxima = maxima(:,:,indx2); exceed = exceed(:,:,indx2);
%     indxx = index(indx2);       Area2 = {MAP(indxx).name};
%     
%     nAp   = nAp(indx2,:);       
%     scal2 = [0 max(max(nAp))];      
%     nAp   = round((nAp-scal2(1))/(scal2(2)-scal2(1))*99)+201; 
%     nAp(isnan(nAp)) = 301;
%     
%     ax=axes('Parent',st.FX,'units','normalized', 'Position',[0.16, .09, .68, .37], 'Visible','off');
%     image(nAp), colormap(MAP(1).cmap)
%     set(gca,'YTick',1:numel(indxx),'XTick',[1:10],'YTickLabel',Area2,'XTickLabel',[10:10:100],'FontSize',FS(9),'FontWeight','bold','TickLength',[0 0])
%     set(get(gca,'XLabel'),'String','cytoarchitectonic probability in %','FontSize',12,'FontWeight','bold')
%  
% 	for i=1:numel(indx2)
%         for ii=1:10
%             mxStr = {};
%             mxStr2 = {};            
%             welche = find(maxima(:,ii,i)>0);
%             for iii=1:numel(welche)
%                 
%                 if exceed(welche(iii),ii,i) == min(min(squeeze(exceed(welche(iii),:,:))))
%                     mxStr{end+1,1}  = [int2str(welche(iii)) ') ' num2str(exceed(welche(iii),ii,i),'%2.0f')];
% %                    mxStr{end+1,1}  = ['(' int2str(welche(iii)) ') ' num2str(exceed(welche(iii),ii,i),'%2.0f')];
%                     mxStr2{end+1,1} = [''];
%                 else
%                     mxStr{end+1,1}  = [''];
%                     mxStr2{end+1,1} = [int2str(welche(iii)) ') ' num2str(exceed(welche(iii),ii,i),'%2.0f')];
%                 end
%             end
%             text(ii,i,mxStr,'Color','k','FontSize',FS(10),'HorizontalAlignment','center','VerticalAlignment','middle','FontWeight','bold')
%             text(ii,i,mxStr2,'Color','k','FontSize',FS(9),'HorizontalAlignment','center','VerticalAlignment','middle')
%         end
%     end
%     
%     ax=axes('Parent',st.FX,'units','normalized', 'Position',[0.88, .09, .05, .37], 'Visible','off');
%     image(flipud(round(repmat(linspace(201,300)',1,20)))), colormap(MAP(1).cmap)
%     set(gca,'YTick',[1 100],'YTickLabel',{[int2str(scal2(2)) '%'];'0%'},'XTick',[],'YAxisLocation','right','TickLength',[0 0])
%     end
end

