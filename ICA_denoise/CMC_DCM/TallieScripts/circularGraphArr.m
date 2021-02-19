classdef circularGraphArr < handle
% CIRCULARGRAPH Plot an interactive circular graph to illustrate connections in a network.
%
%% Syntax
% circularGraphArr(X)
% circularGraphArr(X,'PropertyName',propertyvalue,...)
% h = circularGraphArr(...)
%
%% Description
% A 'circular graph' is a visualization of a network of nodes and their
% connections. The nodes are laid out along a circle, and the connections
% are drawn within the circle. Click on a node to make the connections that
% emanate from it more visible or less visible. Click on the 'Show All'
% button to make all nodes and their connections visible. Click on the
% 'Hide All' button to make all nodes and their connections less visible.
%
% Required input arguments.
% X : A symmetric matrix of numeric or logical values.
%
% Optional properties.
% Colormap : A N by 3 matrix of [r g b] triples, where N is the 
%            length(adjacenyMatrix).
% Label    : A cell array of N strings.
%%
% Copyright 2016 The MathWorks, Inc.
  properties
    Node = node(0,0); % Array of nodes
    ColorMap;         % Colormap
    ConnClrs;
    ConnOrNodeClr;
    ScaledClrTF;
    scaler
    Head;             % Matrix of -1 and 1 for arrow or dot (excit / inhib)
    Label;            % Cell array of strings
    ShowButton;       % Turn all nodes on
    HideButton;       % Turn all nodes off
  end
  
  methods
    function this = circularGraphArr(adjacencyMatrix,varargin)
      % Constructor
      p = inputParser;
      
      defaultConnOrNodeClr = 'Node';
      defaultColorMap = parula(length(adjacencyMatrix));
      defaultConnClrs = parula(length(find(adjacencyMatrix(:)))+1);
      defaultScaledClrTF = 0;
      defaultScaler = 100;
      defaultLabel = cell(length(adjacencyMatrix));
      for i = 1:length(defaultLabel)
        defaultLabel{i} = num2str(i);
      end
      defaultHead = zeros(size(adjacencyMatrix));
      
      addRequired(p,'adjacencyMatrix',@(x)(isnumeric(x) || islogical(x)));
      addParameter(p,'ColorMap',defaultColorMap,@(colormap)length(colormap) == length(adjacencyMatrix));
      addParameter(p,'Label'   ,defaultLabel   ,@iscell);
      addParameter(p,'Head'    ,defaultHead     ,@(x)(isnumeric(x)));
      addParameter(p,'ConnClrs',defaultConnClrs ,@(x)(isnumeric(x)));
      addParameter(p,'ConnOrNodeClr',defaultConnOrNodeClr);
      addParameter(p,'ScaledClrTF',defaultScaledClrTF);
      addParameter(p,'scaler',defaultScaler);
      
      parse(p,adjacencyMatrix,varargin{:});
      this.ColorMap = p.Results.ColorMap;
      this.Label    = p.Results.Label;
      this.Head     = p.Results.Head;
      this.ConnClrs     = p.Results.ConnClrs;
      this.ConnOrNodeClr     = p.Results.ConnOrNodeClr;
      this.ScaledClrTF     = p.Results.ScaledClrTF;
      this.scaler     = p.Results.scaler;
      
      this.ShowButton = uicontrol(...
        'Style','pushbutton',...
        'Position',[0 40 80 40],...
        'String','Show All',...
        'Callback',@circularGraphArr.showNodes,...
        'UserData',this);
      
      this.HideButton = uicontrol(...
        'Style','pushbutton',...
        'Position',[0 0 80 40],...
        'String','Hide All',...
        'Callback',@circularGraphArr.hideNodes,...
        'UserData',this);
      
      fig = gcf;
      set(fig,...
        'UserData',this,...
        'CloseRequestFcn',@circularGraphArr.CloseRequestFcn);
      
      % Draw the nodes
      delete(this.Node);
      t = linspace(-pi,pi,length(adjacencyMatrix) + 1).'; % theta for each node
      extent = zeros(length(adjacencyMatrix),1);
      for i = 1:length(adjacencyMatrix)
        this.Node(i) = node(cos(t(i)),sin(t(i)));
        this.Node(i).Color = this.ColorMap(i,:);
        this.Node(i).Label = this.Label{i};
      end
      
      % Find non-zero values of s and their indices
      [row,col,vi] = find(adjacencyMatrix);
      
      % Calculate line widths based on values of s (stored in v).
      if this.ScaledClrTF
          minLineWidth  = 2;
      else, minLineWidth  = .5;
      end
      lineWidthCoef = 5;
      L = vi;%./nanmax(vi);
      if this.ScaledClrTF
          lineWidth = ones(size(vi));
      else lineWidth = L;
      end
      if sum(lineWidth) == numel(lineWidth) % all lines are the same width.
        lineWidth = repmat(minLineWidth,numel(lineWidth),1);
      else % lines of variable width.
        lineWidth = lineWidthCoef*lineWidth + minLineWidth;
      end
      
      for i = 1:length(vi)
       if ~isnan(vi(i))
          if strcmp('Node',this.ConnOrNodeClr)
              clr = this.ColorMap(row(i),:);
          else
              if this.ScaledClrTF
                  try clr = this.ConnClrs(ceil(L(i)*this.scaler),:);
                  catch
                      disp('um.. had to use black for a connection - not sure what''s going on with the clim')
                      clr = [0 0 0];
                  end
              else, clr = this.ConnClrs(i,:);
              end
          end
          HL = 10;
          HW = 10;
        if this.Head(row(i),col(i))<0
            HeadStyle = 'ellipse';
        elseif this.Head(row(i),col(i))>0
            HeadStyle = 'plain';
        else
            HeadStyle = 'none';
        end
        if row(i) ~= col(i)
          if abs(row(i) - col(i)) - length(adjacencyMatrix)/2 == 0 
            % points are diametric, so draw a straight line
            u = [cos(t(row(i)));sin(t(row(i)))]; % start node
            v = [cos(t(col(i)));sin(t(col(i)))]; % end node
            x0 = -(u(2)-v(2))/(u(1)*v(2)-u(2)*v(1));
            y0 =  (u(1)-v(1))/(u(1)*v(2)-u(2)*v(1));
            lx = linspace(u(1),v(1),100);
            ly = linspace(u(2),v(2),100);
            this.Node(row(i)).Connection(end+1) = line(...
              lx,...[u(1);v(1)]
              ly,...[u(2);v(2)]
              'LineWidth', lineWidth(i),...
              'LineStyle','-',...
              'Color', clr,...
              'PickableParts','none');
          
            ip = 10;
            hold on
            fx = findnearest(this.Node(i).Position(1),lx); fx = fx(1);
            fy = findnearest(this.Node(i).Position(2),ly); fy = fy(1);
            if any([fx fy] > 95)
                xyuv = [lx(end-ip),...
                  ly(end-ip),...
                  diff([lx(end-ip) lx(end)]),...
                  diff([ly(end-ip) ly(end)])];
            elseif any([fx fy] < 5)
                xyuv = [lx(ip),...
                  ly(ip),...
                  diff([lx(ip) lx(1)]),...
                  diff([ly(ip) ly(1)])];
            else fx, fy
                error('I fear something is wrong')
            end
            ah = annotation('arrow','HeadStyle',HeadStyle,'HeadLength',HL,'HeadWidth',HW,'Color',clr);
            set(ah,'Parent',gca);
            set(ah,'Position',xyuv);
          
          else % points are not diametric, so draw an arc
            u  = [cos(t(row(i)));sin(t(row(i)))]; % start node
            v  = [cos(t(col(i)));sin(t(col(i)))]; % end node
            x0 = -(u(2)-v(2))/(u(1)*v(2)-u(2)*v(1));
            y0 =  (u(1)-v(1))/(u(1)*v(2)-u(2)*v(1));
            r  = sqrt(x0^2 + y0^2 - 1);
            thetaLim(1) = atan2(u(2)-y0,u(1)-x0);
            thetaLim(2) = atan2(v(2)-y0,v(1)-x0);
            
            if u(1) >= 0 && v(1) >= 0 
              % ensure the arc is within the unit disk
              theta = [linspace(max(thetaLim),pi,50),...
                       linspace(-pi,min(thetaLim),50)].';
            else
              theta = linspace(thetaLim(1),thetaLim(2)).';
            end
            
            this.Node(row(i)).Connection(end+1) = line(...
              r*cos(theta)+x0,...
              r*sin(theta)+y0,...
              'LineWidth', lineWidth(i),...
              'LineStyle','-',...
              'Color', clr,...
              'PickableParts','none');
          
            ip = 10;
            hold on
            fx = findnearest(this.Node(row(i)).Position(1),r*cos(theta)+x0); fx = fx(1);
            fy = findnearest(this.Node(row(i)).Position(2),r*sin(theta)+y0); fy = fy(1);
            if any([fx fy] > 95)
                xyuv = [r*cos(theta(end-ip))+x0,...
                  r*sin(theta(end-ip))+y0,...
                  diff([r*cos(theta(end-ip))+x0 r*cos(theta(end))+x0]),...
                  diff([r*sin(theta(end-ip))+y0 r*sin(theta(end))+y0])];
            elseif any([fx fy] < 5)
                xyuv = [r*cos(theta(ip))+x0,...
                  r*sin(theta(ip))+y0,...
                  diff([r*cos(theta(ip))+x0 r*cos(theta(1))+x0]),...
                  diff([r*sin(theta(ip))+y0 r*sin(theta(1))+y0])];
            else, warning(['I fear something is wrong with i=' num2str(i) ' fx=' num2str(fx) ' fy=' num2str(fy)])
            end
            ah = annotation('arrow','HeadStyle',HeadStyle,'HeadLength',HL,'HeadWidth',HW,'Color',clr);
            set(ah,'Parent',gca);
            set(ah,'Position',xyuv);
            
          end
        end
       end
      end
      
      axis image;
      ax = gca;
      for i = 1:length(adjacencyMatrix)
        extent(i) = this.Node(i).Extent;
      end
      extent = max(extent(:));
      ax.XLim = ax.XLim + extent*[-1 1];
      fudgeFactor = 1.75; % Not sure why this is necessary. Eyeballed it.
      ax.YLim = ax.YLim + fudgeFactor*extent*[-1 1];
      ax.Visible = 'off';
      set(findall(ax,'Type','text'),'Visible','on')
      ax.SortMethod = 'depth';
      
      fig = gcf;
      fig.Color = [1 1 1];
    end
    
  end
  
  methods (Static = true)
    function showNodes(this,~)
      % Callback for 'Show All' button
      n = this.UserData.Node;
      for i = 1:length(n)
        n(i).Visible = true;
      end
    end
    
    function hideNodes(this,~)
      % Callback for 'Hide All' button
      n = this.UserData.Node;
      for i = 1:length(n)
        n(i).Visible = false;
      end
    end
    
    function CloseRequestFcn(this,~)
      % Callback for figure CloseRequestFcn
      c = this.UserData;
      for i = 1:length(c.Node)
        delete(c.Node(i));
      end
      delete(gcf);
    end
    
  end
  
end