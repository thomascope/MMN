function circuit_diagram(this_DCM,this_BMA,diagnosis_list,source_names)
%A script for plotting the results of extDCM

xw = 400; %X_width
yw = 600; %Y_width

for j = 1:length(source_names)
    
    %Create figure
    circuit_diagram = figure('Position',[1, 1, xw, yw+20]);
    xlim([0 xw+40])
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
    
    for i = 1:size(this_DCM.Ep.H,4)
        line_code=[population_colors(i) '--'];
        for from = 1:6
            for to = 1:6
                if this_DCM.Ep.H(from,to,j,i)~=0
                    draw_arrow(from, to, coords_top, coords_bot, coords_cen, xw, yw, line_code)
                end
            end
        end
    end
    daspect([1 1 1])
    
    title(source_names{j})
    set(findall(gca, 'type', 'text'), 'visible', 'on')
end

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
    %Special case of stellate - origin centre
%     if from == 1
%         coords_bot = coords_cen;
%     end
    
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
        plot([coords_bot(from,1),coords_bot(from,1)-xw/10]+offset_amount,[coords_bot(from,2) - xw/15 - offset_amount,coords_bot(from,2) - xw/15 - offset_amount],line_code) %Left
        plot([coords_bot(from,1)-xw/10,coords_bot(from,1)-xw/10]+offset_amount,[coords_bot(from,2) - xw/15 - offset_amount,coords_top(to,2)],line_code) %Up
        arrow3([coords_bot(from,1)-xw/10+offset_amount,coords_top(to,2)],coords_top(to,:),line_code) % note arrow3 syntax takes point pairs rather than x then y values
        %quiver([coords_bot(from,1)+offset_amount],[coords_top(to,2)],coords_top(to,1)-coords_bot(from,1)-offset_amount,0,0,line_code,'filled')
    end
    
end

