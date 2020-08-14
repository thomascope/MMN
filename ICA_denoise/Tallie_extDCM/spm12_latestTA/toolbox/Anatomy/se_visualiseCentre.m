function se_visualiseCentre(op)

global MAP
global st

try
    if isempty(op)
        op = '';
    end
catch
        op = '';
end

if strcmp(op,'shopos')

    
    
    
elseif strcmp(op,'exit')
    fg = spm_figure('GetWin','Graphics'); spm_figure('Clear','Graphics'); clear
    Anatomy;

else

    MapName = spm_select(1,'mat',['Select Map'],[],spm('Dir','se_anatomy'),'MPM.',1);
    se_getMap('anat',MapName);

    P      = deblank([MapName(1:end-8) '.img']);

    fg = se_figure('GetWin','Graphics');
    if isempty(fg), error('Can''t create graphics window'); end
    se_figure('Clear','Graphics');

    se_orthviews('Reset');
    switch spm('FnBanner')
        case 'SPM99'
            se_orthviews('Image', P, [0.0 0.22 1 .8]);     % Image appears
        otherwise
            se_orthviews('Image', P, [0.0 0.2 1 .8]);     % Image appears
    end
    if isempty(st.vols{1}), return; end;

    se_orthviews('MaxBB');
    st.callback = 'se_visualiseCentre(''shopos'');';

    st.B = [0 0 0  0 0 0  1 1 1  0 0 0];


    LRnamen = {'left';'';'right'};
    raus = -1;

    while raus < 1

        if raus > -1
            raus = spm_input_ui('Continue?','!+1','b','Yes|No, Exit',[0 1],1);
        else
            raus = 0;
        end

        if raus < 1
            spm_figure('Clear','Interactive');
            [B,index] = sortrows(char(MAP.name));
            locNr = index(spm_input('Which area ?', '!+1', 'm',str2mat(B)));
            side = spm_input_ui('Side','!+1','b','Left|Right',[-1 1],1,'Select side');
            mode = spm_input_ui('Mode','!+1','m','All assigned|Highest probability',[0 1],1,'Select ROI definition mode');


            if mode == 1

                total = sum(MAP(locNr).allLR == side); perc = zeros(1,10);
                for i=1:10, perc(i) = sum(MAP(locNr).allLR == side & MAP(locNr).allZ >= i) / total*100; end
                perc2 = abs(perc-10); thres = find(perc2 == min(perc2));


                QQ = (MAP(locNr).allZ >= thres) & (MAP(locNr).allLR == side);

                se_orthviews('rmblobs',1); try, delete(tx), end
                se_orthviews('addcolouredblobs',1,MAP(locNr).allXYZ(:,QQ),ones(1,sum(QQ)),MAP(1).MaxMap.mat,[1 0 0]);
                center = mean(MAP(1).MaxMap.mat * [MAP(locNr).allXYZ(:,QQ); ones(1,size(MAP(locNr).allXYZ(:,QQ),2))],2);
                se_orthviews('Reposition',center(1:3));
                se_orthviews('Xhairs','on')
                fg = spm_figure('GetWin','Graphics'); figure(fg)
                tx = text(100,-50,{[LRnamen{side+2} '  ' MAP(locNr).name]; [int2str(sum(QQ)) ' voxel'];['Threshold: prob. >= ' int2str(thres*10) '%']},'HorizontalAlignment','center','FontSize',20);

            else

                QQ = MAP(locNr).LR == side;

                se_orthviews('rmblobs',1); try, delete(tx), end
                se_orthviews('addcolouredblobs',1,MAP(locNr).XYZ(:,QQ),ones(1,sum(QQ)),MAP(1).MaxMap.mat,[1 0 0]);
                center = mean(MAP(1).MaxMap.mat * [MAP(locNr).XYZ(:,QQ); ones(1,size(MAP(locNr).XYZ(:,QQ),2))],2);
                se_orthviews('Reposition',center(1:3));
                se_orthviews('Xhairs','on')
                fg = spm_figure('GetWin','Graphics'); figure(fg)
                tx = text(100,-50,{[LRnamen{side+2} '  ' MAP(locNr).name]; [int2str(sum(QQ)) ' voxel'];['Threshold: MPM']},'HorizontalAlignment','center','FontSize',20);
            end
        else
            se_visualiseCentre('exit')
        end
    end
end
