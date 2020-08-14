function se_render(cl)

global st
global CLUSTER
global MAP


try
    dat = CLUSTER(cl).dat;
catch
    dat = [];
end
if isempty(dat)
    dat = struct('XYZ',CLUSTER(cl).gtmp,'t',CLUSTER(cl).gZ,'mat',MAP(1).MaxMap.mat,'dim',MAP(1).MaxMap.dim(1:3));
    dat.mat(2,4) = dat.mat(2,4) + 4;
    dat.mat(3,4) = dat.mat(3,4) - 5;
    
    CLUSTER(cl).dat = dat;
end


try
    ren = dat.ren;
catch
    ren = {};
end
if isempty(ren)
    
    
    
    rend = MAP(1).rend;
    
    mx = zeros(length(rend),1)+eps;
    mn = zeros(length(rend),1);
    
    for j=1:length(dat),
        XYZ = dat(j).XYZ;
        t   = dat(j).t;
        dim = dat(j).dim;
        mat = dat(j).mat;
        depView = zeros(1,6);
        
        for i=1:length(rend),
            M1  = rend{i}.M*dat(j).mat;
            zm  = sum(M1(1:2,1:3).^2,2).^(-1/2);
            M2  = diag([zm' 1 1]);
            M  = M2*M1;
            cor = [1 1 1 ; dim(1) 1 1 ; 1 dim(2) 1; dim(1) dim(2) 1 ;
                1 1 dim(3) ; dim(1) 1 dim(3) ; 1 dim(2) dim(3); dim(1) dim(2) dim(3)]';
            tcor= M(1:3,1:3)*cor + M(1:3,4)*ones(1,8);
            off = min(tcor(1:2,:)');
            M2  = spm_matrix(-off+1)*M2;
            M  = M2*M1;
            xyz = (M(1:3,1:3)*XYZ + M(1:3,4)*ones(1,size(XYZ,2)));
            d2  = ceil(max(xyz(1:2,:)'));
            dep = spm_slice_vol(rend{i}.dep,spm_matrix([0 0 1])*inv(M2),d2,1);
            z1  = dep(round(xyz(1,:))+round(xyz(2,:)-1)*size(dep,1));
            msk = find(xyz(3,:) < (z1+50) & xyz(3,:) > (z1-2));
            if ~isempty(msk),
                depView(i) = mean(xyz(3,msk)-(z1(msk)+35));
                xyz = xyz(:,msk); t0  = t(msk);
                X0  = full(sparse(round(xyz(1,:)), round(xyz(2,:)), t0, d2(1), d2(2)));
                hld = 0;
                X   = spm_slice_vol(X0,spm_matrix([0 0 1])*M2,size(rend{i}.dep),hld);
                msk = find(X<0);
                X(msk) = 0;
            else,
                X = zeros(size(rend{i}.dep));
            end;
            mx(j) = max([mx(j) max(max(X))]);
            mn(j) = min([mn(j) min(min(X))]);
            rend{i}.data{j} = X;
        end;
        
    end;
    
    
    mxmx = max(mx); mnmn = min(mn);
    
    [B I] = sort(depView); I = I(1:2);
    
    startFigure
    for xi=1:2,
        ren = rend{I(xi)}.ren;
        X   = (rend{I(xi)}.data{1}-mnmn)/(mxmx-mnmn);
        msk = find(X);
        ren(msk) = X(msk)+(1+1.51/64);
        if xi == 1
            ax=axes('Parent',st.Fgraph,'units','normalized', 'Position',[0, 0, 1, .5], 'Visible','off');
        else
            ax=axes('Parent',st.Fgraph,'units','normalized', 'Position',[0, .5, 1, .5], 'Visible','off');
        end
        image(ren*64,'Parent',ax);
        CLUSTER(cl).dat.ren{xi} = ren*64;
        set(ax,'DataAspectRatio',[1 1 1], 'PlotBoxAspectRatioMode','auto', 'YTick',[],'XTick',[],'XDir','normal','YDir','normal');
    end;
    
    
else
    startFigure
    for xi=1:2,
        if xi == 1
            ax=axes('Parent',st.Fgraph,'units','normalized', 'Position',[0, 0, 1, .5], 'Visible','off');
        else
            ax=axes('Parent',st.Fgraph,'units','normalized', 'Position',[0, .5, 1, .5], 'Visible','off');
        end
        image(CLUSTER(cl).dat.ren{xi},'Parent',ax);
        set(ax,'DataAspectRatio',[1 1 1], 'PlotBoxAspectRatioMode','auto', 'YTick',[],'XTick',[],'XDir','normal','YDir','normal');
    end;
end








spm('Pointer')
return;









function startFigure
global st

st.Fgraph = findobj(get(0,'Children'),'Flat','Tag','xSATB');
if ~isempty(st.Fgraph)
    try
        set(0,'CurrentFigure',st.Fgraph);
        clf
    catch
        st.Fgraph = [];
    end
end
if isempty(st.Fgraph)
    FS   = spm('FontSizes');                 PF   = spm_platform('fonts');
    
    screen = get(0,'ScreenSize');
    st.Fgraph = figure(...
        'Tag','xSATB',                             'Position',[screen(3)-255 screen(4)-590 250 500],...
        'Resize','off',                           'MenuBar','figure',...
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
    set(st.Fgraph,'Name',[spm('ver') ' Anatomy Toolbox'],'NumberTitle','off')
    
else
    figure(st.Fgraph)
end
ax=axes('Parent',st.Fgraph,'units','normalized','Position',[0, 0, 1, 1],'Visible','off'); image(0,'Parent',ax); set(ax,'YTick',[],'XTick',[]);
load Split;
colormap(split);
