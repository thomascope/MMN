function circuit_diagram(dirname_DCM,diagnosis_list,source_names,conductances,thresh)
%A script for plotting the results of extDCM across all diagnoses and
%sources. PEB of PEBs on H matrix.
addpath(['/group/language/data/thomascope/MMN/ICA_denoise/extDCM_visualisation/ojwoodford-export_fig-216b30e'])
addpath('/group/language/data/thomascope/MMN/ICA_denoise/Helperfiles')
thisdir = pwd;
mkdir([thisdir '/circuit_diagrams'])
cd([dirname_DCM 'PEB_secondlevel'])

all_folders = strsplit(dirname_DCM,'/');
all_folders = all_folders(~cellfun(@isempty,all_folders)); %For save location

cleanupObj = onCleanup(@()cd(thisdir));

%Load any input DCM and PEB of PEBs as template
this_DCM = load(['PEB_H_' source_names{1} '_' conductances{1} '_' diagnosis_list{1} '.mat']);
this_DCM = this_DCM.DCM{1};

template_PEB = load(['PEB_H_' source_names{1} '_' conductances{1} '_Overall.mat']);
template_PEB = template_PEB.PEB_Overall;

assert(all(template_PEB.M.X(:,1)==1),'The first column of the PEB of PEBs contrast should be all ones, check please.')

xw = 400; %X_width
yw = 600; %Y_width

for this_contrast = 0:size(template_PEB.M.X,2)
    for this_source = 1:length(source_names)
        
        %Create figure
        circuit_diagram = figure('Position',[1, 1, xw, yw+20]);
        xlim([-40 xw+40])
        ylim([0 yw+25])
        set(gca,'visible','off')
        hold on
        
        coords_top = [];
        coords_bot = [];
        coords_cen = [];
        
        %Draw Spiny Stellate - cellpop 1
        % Define parameters
        numberOfPoints = 5;
        rotationAngle = pi/10;
        xCenter = (xw/2+xw/8);
        yCenter = yw/2;
        % Determine the angles that the arm tips are at
        theta = (0 : (numberOfPoints-1)/numberOfPoints*pi : (numberOfPoints-1)*pi) + rotationAngle;
        % Define distance from the arm tip to the center of star.
        amplitude = xw/10;
        % Get x and y coordinates of the arm tips.
        x_star = amplitude .* cos(theta) + xCenter;
        y_star = amplitude .* sin(theta) + yCenter;
        % Draw it.
        spin_stel = fill(x_star,y_star,[.7 .7 .7]);
        coords_top = [coords_top; (xw/2+xw/8),yw/2+xw/10];
        coords_cen = [coords_cen; (xw/2+xw/8),yw/2];
        
        % Calculate bottom crossing of star
        left_point = [x_star(find(x_star==min(x_star))), y_star(find(x_star==min(x_star)))];
        right_point = [x_star(find(x_star==max(x_star))), y_star(find(x_star==max(x_star)))];
        right_point = right_point([1,3]);
        bottom_points = [x_star(find(y_star==min(y_star))), y_star(find(y_star==min(y_star)))];
        
        [x,y]=lineintersect([left_point,bottom_points([1,3])],[right_point,bottom_points([2,4])]);
        coords_bot = [coords_bot;x,y];
        
        %Draw superficial pyramidal - cellpop 2
        sup_py = fill([xw/2,xw/2-xw/10,xw/2+xw/10],[yw,yw-yw/10,yw-yw/10],'b');
        coords_top = [coords_top; xw/2,yw];
        coords_bot = [coords_bot; xw/2,yw-yw/10];
        coords_cen = [coords_cen; xw/2,yw-yw/20];
        
        %Draw Superficial interneuron - cellpop 3
        circle_radius = xw/15;
        deep_int = rectangle('Position',[xw/8 - circle_radius, yw-yw/4  - circle_radius, circle_radius*2, circle_radius*2],...
            'Curvature',[1,1],...
            'FaceColor','r');
        coords_top = [coords_top; xw/8, yw-yw/4  + circle_radius];
        coords_bot = [coords_bot; xw/8, yw-yw/4  - circle_radius];
        coords_cen = [coords_cen; xw/8,yw-yw/4];
        
        %Draw deep pyramidal - cellpop 4
        deep_py = fill([(xw/2-xw/6),(xw/2-xw/6)-xw/10,(xw/2-xw/6)+xw/10],[yw/6,yw/6-yw/10,yw/6-yw/10],'g');
        coords_top = [coords_top; xw/2-xw/6,yw/6];
        coords_bot = [coords_bot; xw/2-xw/6,yw/6-yw/10];
        coords_cen = [coords_cen; xw/2-xw/6,yw/6-yw/20];
        
        %Draw Deep interneuron - cellpop 5
        circle_radius = xw/15;
        deep_int = rectangle('Position',[xw/8 - circle_radius, yw/3  - circle_radius, circle_radius*2, circle_radius*2],...
            'Curvature',[1,1],...
            'FaceColor','r');
        coords_top = [coords_top; xw/8, yw/3  + circle_radius];
        coords_bot = [coords_bot; xw/8, yw/3  - circle_radius];
        coords_cen = [coords_cen; xw/8,yw/3];
        
        %Draw Thalamic projection - cellpop 6
        thal_proj = fill([(xw-xw/8),(xw-xw/8)-xw/10,(xw-xw/8)+xw/10],[yw/6,yw/6-yw/10,yw/6-yw/10],'k');
        coords_top = [coords_top; (xw-xw/8),yw/6];
        coords_bot = [coords_bot; (xw-xw/8),yw/6-yw/10];
        coords_cen = [coords_cen; (xw-xw/8),yw/6-yw/20];
        
        % %Test locations (comment out in final code)
        % plot(coords_top(:,1),coords_top(:,2),'*')
        % plot(coords_bot(:,1),coords_bot(:,2),'s')
        % plot(coords_cen(:,1),coords_cen(:,2),'o')
        
        %Now add links between populations - requires arrow3 from file exchange
        % Population Order: AMPA, NMDA, GABAa, GABAb
        population_colors=['kbrg'];
        
        if this_contrast == 0
            for i = 1:size(this_DCM.Ep.H,4)
                line_code=[population_colors(i) '--'];
                for from = 1:6
                    for to = 1:6
                        if this_DCM.Ep.H(to,from,this_source,i)~=0&&sum(this_DCM.Ep.H(to,from,this_source,:)~=0)==1
                            draw_arrow(from, to, coords_top, coords_bot, coords_cen, xw, yw, line_code)
                        elseif this_DCM.Ep.H(to,from,this_source,i)~=0
                            overlapping_lines_before = sum(this_DCM.Ep.H(to,from,this_source,1:i-1)~=0);
                            this_coords_bot = coords_bot-(overlapping_lines_before*(xw/200));
                            this_coords_top = coords_top-(overlapping_lines_before*(xw/200));
                            draw_arrow(from, to, this_coords_top, this_coords_bot, coords_cen, xw, yw, line_code)
                        end
                    end
                end
            end
            daspect([1 1 1])
            
            this_title = 'Template circuit';
            title(this_title)
            set(findall(gca, 'type', 'text'), 'visible', 'on')
            saveas(circuit_diagram,[thisdir '/circuit_diagrams/' this_title '.jpg'])
            saveas(circuit_diagram,[thisdir '/circuit_diagrams/' this_title '.pdf'])
              eval(['export_fig ''' thisdir '/circuit_diagrams/' this_title '.png'' -transparent -nocrop'])
            break
        else
            for condition = 1:2
                all_to_froms = [];
                for this_conductance = 1:length(conductances)
                    load(['PEB_H_' source_names{this_source} '_' conductances{this_conductance} '_Overall.mat'])
                    these_differences = find(BMA_Overall.Pp(:,this_contrast)>thresh);
                    
                    for this_difference = 1:length(these_differences)
                        if condition == 1
                            if BMA_Overall.Ep(these_differences(this_difference),this_contrast)>0 %Controls greater, dashed line
                                line_code=[population_colors(this_conductance) '--'];
                            else
                                line_code=[population_colors(this_conductance) '-']; %Patients greater, solid line
                            end
                        else
                            line_code=[population_colors(this_conductance) ':']; %Interaction, dotted line
                        end
                        
                        this_connection = BMA_Overall.Pnames{these_differences(this_difference)};
                        Condition_Split = strsplit(this_connection,'Covariate ');
                        if str2num(Condition_Split{2}(1))~=condition %Not the condition we want
                            continue
                        else
                            Connection_Split = strsplit(this_connection,'H(');
                            to = str2num(Connection_Split{2}(1));
                            from = str2num(Connection_Split{2}(3));
                            if isempty(all_to_froms)||~any(ismember(all_to_froms,[from,to],'rows'))
                                draw_arrow(from, to, coords_top, coords_bot, coords_cen, xw, yw, line_code)
                            else %Make slight offset to avoid overlapping lines
                                overlapping_lines_before = sum(ismember(all_to_froms,[from,to],'rows'));
                                this_coords_bot = coords_bot-(overlapping_lines_before*(xw/200));
                                this_coords_top = coords_top-(overlapping_lines_before*(xw/200));
                                draw_arrow(from, to, this_coords_top, this_coords_bot, coords_cen, xw, yw, line_code)
                            end
                            all_to_froms = [all_to_froms; from,to];
                        end
                    end
                end
                daspect([1 1 1])
                
                
            end
            if this_contrast == 1
                this_title = ['All positives ' source_names{this_source}];
            else
            this_title = [diagnosis_list{find(template_PEB.M.X(:,this_contrast)==1)} ' minus ' diagnosis_list{find(template_PEB.M.X(:,this_contrast)==-1)} ' ' source_names{this_source}];
            end
            title(this_title)
            set(findall(gca, 'type', 'text'), 'visible', 'on')
            
            %eval(['export_fig ' thisdir '/circuit_diagrams/' this_title '.pdf -transparent'])
            %             saveas(circuit_diagram,[thisdir '/circuit_diagrams/' this_title '.png'])
            saveas(circuit_diagram,[thisdir '/circuit_diagrams/' this_title '.pdf'])
            eval(['export_fig ''' thisdir '/circuit_diagrams/' this_title '.png'' -transparent -nocrop'])
            close all
        end
    end
end
cd(thisdir)

function draw_arrow(from, to, coords_top, coords_bot, coords_cen, xw, yw, line_code) %assumes a total of 6 populations for offsetting
if from == to %Special case, self inhibition
    % Draw an arc between two points
    a=coords_bot(from,:);
    b=coords_top(to,:);
    
    % Define centre
    Center=(a+b)./2;
    xc=Center(1);
    yc=Center(2);
    % Radius
    r=norm((a-b))/2;
    % Subtended angle
    theta = pi/2;
    arc_distance = [coords_top(to,2)-coords_bot(from,2)];
    r = arc_distance/sqrt(2-2*cos(theta));
    horiz_offset=sqrt(r^2-(arc_distance/2)^2);  %By Pythagoras
    
    all_thetas = linspace(pi+(theta/2),3*pi-(theta/2),200); %all angles
    
    circle_xs = [xc + horiz_offset + r*cos(all_thetas)];
    circle_ys = [yc + r*sin(all_thetas)];
    
    plot(circle_xs,circle_ys,line_code) % Plot line
    fill([circle_xs(100),circle_xs(100)-r/4,circle_xs(100),circle_xs(100)+r/4],[circle_ys(100),circle_ys(100)-r/2,circle_ys(100)-r/3,circle_ys(100)-r/2],line_code(1))
    
else
    
    %Make a little offset to prevent arrows overlapping
    if from<to
        offset_amount = (to-4)*xw/100;
    elseif from>to
        offset_amount = (to-3)*xw/100;
    end
    
    %First draw vertical component then draw arrow
    if coords_bot(from,2) > coords_top(to,2) % Is above, go straight down
        plot([coords_bot(from,1),coords_bot(from,1)]+offset_amount,[coords_bot(from,2), coords_top(to,2)],line_code)
        arrow3([coords_bot(from,1)+offset_amount,coords_top(to,2)],coords_top(to,:),line_code)
        %quiver([coords_bot(from,1)+offset_amount],[coords_top(to,2)],coords_top(to,1)-coords_bot(from,1)-offset_amount,0,0,line_code,'filled')
    elseif coords_bot(from,2) < coords_top(to,2) % Is below, go down then around
        plot([coords_bot(from,1),coords_bot(from,1)]+offset_amount,[coords_bot(from,2), coords_bot(from,2) - xw/15 - offset_amount],line_code) %Down
        plot([coords_bot(from,1),coords_bot(from,1)-xw/6]+offset_amount,[coords_bot(from,2) - xw/15 - offset_amount,coords_bot(from,2) - xw/15 - offset_amount],line_code) %Left
        plot([coords_bot(from,1)-xw/6,coords_bot(from,1)-xw/6]+offset_amount,[coords_bot(from,2) - xw/15 - offset_amount,coords_top(to,2)],line_code) %Up
        arrow3([coords_bot(from,1)-xw/6+offset_amount,coords_top(to,2)],coords_top(to,:),line_code) % note arrow3 syntax takes point pairs rather than x then y values
        %quiver([coords_bot(from,1)+offset_amount],[coords_top(to,2)],coords_top(to,1)-coords_bot(from,1)-offset_amount,0,0,line_code,'filled')
    end
    
end
