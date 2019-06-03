function linehandle = stdshade_TEC_cmap(amatrix,alpha,acolor,F,smth,steflag,line_style)
% usage: stdshading(amatrix,alpha,acolor,F,smth,steflag)
% plot mean and sem/std coming from a matrix of data, at which each row is an
% observation. sem/std is shown as shading.
% - acolor defines the used color (default is red) 
% - F assignes the used x axis (default is steps of 1).
% - alpha defines transparency of the shading (default is no shading and black mean line)
% - smth defines the smoothing factor (default is no smooth)
% smusall 2010/4/23
% TEC_added steflag to plot standard error instead of stdev

if exist('acolor','var')==0 || isempty(acolor)
    acolor='r'; 
end

if exist('F','var')==0 || isempty(F); 
    F=1:size(amatrix,2);
end

if exist('smth','var'); if isempty(smth); smth=1; end
else smth=1;
end  

if exist('steflag','var'); if isempty(steflag); smth=0; end
else steflag=0;
end  


if ne(size(F,1),1)
    F=F';
end

amean=smooth(nanmean(amatrix),smth)';
astd=nanstd(amatrix); % to get std shading
% astd=nanstd(amatrix)/sqrt(size(amatrix,1)); % to get sem shading
if steflag==1
    astd = astd/sqrt(size(amatrix,1));
end
    

if exist('alpha','var')==0 || isempty(alpha) 
    fill([F fliplr(F)],[amean+astd fliplr(amean-astd)],acolor,'linestyle','none','HandleVisibility','off');
    acolor='k';
else fill([F fliplr(F)],[amean+astd fliplr(amean-astd)],acolor, 'FaceAlpha', alpha,'linestyle','none','HandleVisibility','off');    
end

if ishold==0
    check=true; else check=false;
end

hold on;
linehandle = plot(F,amean,'Color', acolor,'LineStyle',line_style,'linewidth',1.5); %% change color or linewidth to adjust mean line

if check
    hold off;
end

end



