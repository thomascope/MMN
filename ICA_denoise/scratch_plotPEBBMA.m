function scratch_plotPEBBMA(B,EpPp,diagnosis_list,source_names)
%
% This script is designed for looking at PEB results for my original extDCM
% model. If you change the number of intrinsic or extrinsic connections, 
% this script will need altering accordingly! - Edited by TEC for
% multi-group MMN Oct 2020
% 
% B must be a structure with your collection of BMAs in (see example in 
% extDCMpostproc_notes.m'). - Edited by TEC to allow a cell array of
% structs
% 
% EpPp is either 'Ep' or 'Pp' depending on whether you just want to plot 
% the Ep fields (use if plotting multiple BMAs together) or whether you 
% want to plot Ep with it's corresponding Pp in the row below it (only use
% if B contains a single BMA).
%

% Additional notes specific to TEC analysis:
% First half of first level contrast fields are STD+DVT together, second half of
% fields are contrast of STD-DVT. Therefore PEB of PEB results in the first
% half are main effects of group, and in the second half are group by
% condition interactions
% diagnosis_list: {'Control'  'bvFTD'  'pca'  'nfvppa'  'ADMCI'}
% source_names: {'left A1', 'left STG', 'left IFG', 'left IPC', 'right A1', 'right STG', 'right IFG', 'right IPC'}

if nargin==1, EpPp = 'Ep'; end

if iscell(B)
    flds = fieldnames(B{1});
    disp('Multiple BMAs in cell array - assuming all have the same fields')
else
    flds = fieldnames(B);
    B = {B};
end

if ~isempty(regexp(EpPp,'Pp','once')), p = .75; else p = .95; end

%V.cax = .25;
V.cax = .5;

for this_BMA = 1:length(B)
% sort out which parameters your've looking at:
% which are H1:
h1 = cellfun(@(x) ~isempty(regexp(x,'H','once')) && ~isempty(regexp(x,',1)','once')),B{this_BMA}.(flds{1}).bma.Pnames);
s1hi = find(h1);
s1 = cellfun(@(x) strsplit(x,'('),B{this_BMA}.(flds{1}).bma.Pnames(h1),'Uni',0);
s1 = cellfun(@(x) strsplit(x{2},','),s1,'Uni',0);
s1 = cellfun(@(x) [x{1} strsplit(x{2},')')],s1,'Uni',0);
s1 = cellfun(@(y) cell2mat(cellfun(@(x) str2num(x),y,'Uni',0)),s1,'Uni',0);
s1h = cat(1,s1{:});
% which are H3:
h3 = cellfun(@(x) ~isempty(regexp(x,'H','once')) && ~isempty(regexp(x,',7)','once')),B{this_BMA}.(flds{1}).bma.Pnames);
s3hi = find(h3);
s1 = cellfun(@(x) strsplit(x,'('),B{this_BMA}.(flds{1}).bma.Pnames(h3),'Uni',0);
s1 = cellfun(@(x) strsplit(x{2},','),s1,'Uni',0);
s1 = cellfun(@(x) [x{1} strsplit(x{2},')')],s1,'Uni',0);
s1 = cellfun(@(y) cell2mat(cellfun(@(x) str2num(x),y,'Uni',0)),s1,'Uni',0);
s3h = cat(1,s1{:});
% which are A:
ha = cellfun(@(x) ~isempty(regexp(x,'H','once')),B{this_BMA}.(flds{1}).bma.Pnames);
sahi = find(~ha);
s1 = cellfun(@(x) strsplit(x,'('),B{this_BMA}.(flds{1}).bma.Pnames(~ha),'Uni',0);
s1 = cellfun(@(x) strsplit(x{2},','),s1,'Uni',0);
s1 = cellfun(@(x) [x{1} strsplit(x{2},')')],s1,'Uni',0);
s1 = cellfun(@(y) cell2mat(cellfun(@(x) str2num(x),y,'Uni',0)),s1,'Uni',0);
sa = cat(1,s1{:});
end
fmx = max(cellfun(@(x) size(B{this_BMA}.(x).bma.Ep,2),flds));
np = fmx; % 8;
doneEp = 0;
donePp = 0;
ln = length(flds);
nf = ceil(ln/np);

cm = flipud(bone(np));
cm = [ones(24,3).*cm(2,:) ; ones(20,3).*cm(5,:) ; ones(4,3).*cm(7,:) ; ones(1,3).*cm(8,:)]; % THIS IS FOR P OF .5 UPWARDS
% cm = [ones(20,3).*cm(3,:) ; ones(4,3).*cm(5,:) ; ones(1,3).*cm(8,:)]; % THIS IS FOR P .75 UPWARDS

c0 = 0;
c1 = 1:min([np ln]);
for kf = 1:nf
    fm; set(gcf,'Color','w')
    for k = c1
        subplot(4,np,1+(np*(k-c0))-np), imagesc(B.(flds{k}).bma.M.X), axis tight square, colormap(gca,bone(256))
        title([labelTA(flds{k});{'Design Matrix'}])

        wmp = 2:size(B.(flds{k}).bma.Ep,2); % size(B.(flds{k}).bma.M.X,2);
        for kw = 1:length(wmp) % for each covariate
            H1 = zeros(length(source_names)); H1a = zeros(length(source_names));
            for k1 = 1:size(s1h,1)
                if B.(flds{k}).bma.Pp(s1hi(k1),wmp(kw))>=p % find the significantly correlated parameters
                    H1(s1h(k1,2),s1h(k1,1)) = B.(flds{k}).bma.Ep(s1hi(k1),wmp(kw)); % fill a pop*pop matrix with their posterior estimates
                    H1a(s1h(k1,2),s1h(k1,1)) = B.(flds{k}).bma.Pp(s1hi(k1),wmp(kw));
                end
            end
            H3 = zeros(length(source_names)); H3a = zeros(length(source_names));
            for k1 = 1:size(s3h,1)
                if B.(flds{k}).bma.Pp(s3hi(k1),wmp(kw))>=p % find the significantly correlated parameters
                    H3(s3h(k1,2),s3h(k1,1)) = B.(flds{k}).bma.Ep(s3hi(k1),wmp(kw)); % fill a pop*pop matrix with their posterior estimates
                    H3a(s3h(k1,2),s3h(k1,1)) = B.(flds{k}).bma.Pp(s3hi(k1),wmp(kw));
                end
            end
            A = zeros(length(source_names)); Aa = zeros(length(source_names));
            for k1 = 1:size(sa,1)
                if B.(flds{k}).bma.Pp(sahi(k1),wmp(kw))>=p % find the significantly correlated parameters
                    A(sa(k1,2),sa(k1,1)) = A(sa(k1,1),sa(k1,2)) + B.(flds{k}).bma.Ep(sahi(k1),wmp(kw)); % fill a pop*pop matrix with their posterior estimates
                    Aa(sa(k1,2),sa(k1,1)) = Aa(sa(k1,1),sa(k1,2)) + B.(flds{k}).bma.Pp(sahi(k1),wmp(kw));
                end
            end
            
            H1H3A = zeros(3*length(source_names));
            H1H3A(1:length(source_names),1:length(source_names)) = H1;
            H1H3A((length(source_names)+1):2*length(source_names),(length(source_names)+1):2*length(source_names)) = H3;
            H1H3A((2*length(source_names)+1):3*length(source_names),(2*length(source_names)+1):3*length(source_names)) = A;
                        
            H1H3Aa = zeros(3*length(source_names));
            H1H3Aa(1:length(source_names),1:length(source_names)) = H1a;
            H1H3Aa((length(source_names)+1):2*length(source_names),(length(source_names)+1):2*length(source_names)) = H3a;
            H1H3Aa((2*length(source_names)+1):3*length(source_names),(2*length(source_names)+1):3*length(source_names)) = Aa;
            
            g = digraph(H1H3A);
            H1H3Ac = double.empty([0 1]);
            H1H3Ac_a = double.empty([0 1]);
            for i = 1:size(g.Edges{:,1},1)
                H1H3Ac(i) = H1H3A(g.Edges{:,1}(i,1),g.Edges{:,1}(i,2));
                H1H3Ac_a(i) = H1H3Aa(g.Edges{:,1}(i,1),g.Edges{:,1}(i,2));
            end
            
            subplot(4,np,2+(kw-1)+(np*(k-c0))-np)
            plot(g,'XData',[1 1.35 .6 1 .6 1.35   3 3.35 2.6 3 2.6 3.35   1.5 1.2 1.5 2.5 2.8 2.5],...
                'YData',[1.5 2.5 2.35 .35 .15 -.15   1.5 2.5 2.35 .35 .15 -.15   -3 -4 -5 -3 -4 -5],...
                'MarkerSize',5,'ShowArrows','on','NodeColor',[.75 .75 .75],'NodeLabel',{},'EdgeCData',H1H3Ac,'EdgeAlpha',1,'LineWidth',2)
            axis equal tight, xlim([0 4]), set(gca,'Visible','off'), colormap(gca,cbrewer('div','RdBu',256)), caxis([-V.cax V.cax])
            if any(ha), mx_y = 3; else mx_y = -.5; end
            if any(~ha), mn_y = -6; else mn_y = -2.5; end
            ylim([mn_y mx_y])
            if any(h1), mn_x = 0; else mn_x = 2; end
            if any(h3), mx_x = 4; else mx_x = 2; end
            xlim([mn_x mx_x])
            if doneEp==0, colorbar('Location','southoutside','Position',[0.5628 0.0853 0.0526 0.0160]), doneEp = 1; end
            
            if ~isempty(regexp(EpPp,'Pp','once'))
                subplot(4,np,np+2+(kw-1)+(np*(k-c0))-np)
                plot(g,'XData',[1 1.35 .6 1 .6 1.35   3 3.35 2.6 3 2.6 3.35   1.5 1.2 1.5 2.5 2.8 2.5],...
                    'YData',[1.5 2.5 2.35 .35 .15 -.15   1.5 2.5 2.35 .35 .15 -.15   -3 -4 -5 -3 -4 -5],...
                    'MarkerSize',5,'ShowArrows','on','NodeColor',[.75 .75 .75],'NodeLabel',{},'EdgeCData',H1H3Ac_a,'EdgeAlpha',1,'LineWidth',2)
                axis equal tight, xlim([0 4]), set(gca,'Visible','off'), colormap(gca,cm), caxis([p 1])
                if any(ha), mx_y = 3; else mx_y = -.5; end
                if any(~ha), mn_y = -6; else mn_y = -2.5; end
                ylim([mn_y mx_y])
                if any(h1), mn_x = 0; else mn_x = 2; end
                if any(h3), mx_x = 4; else mx_x = 2; end
                xlim([mn_x mx_x])
                if donePp==0, colorbar('Location','southoutside','Position',[0.6628 0.0853 0.0526 0.0160]); doneEp = 1; end
            end
        end
        drawnow
    end
    c0 = c1(end);
    c1 = c1(end)+1:min([length(flds) c1(end)+np]);
end

