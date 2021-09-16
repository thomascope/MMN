%Load prevoiusly extracted whitened and filtered ys
load('MMN_amplitudes_forgraph','MMN_amplitudes_forgraph','R_frontal_volume','R_frontal_volume_mixed','R_parietal_volume','R_parietal_volume_mixed','L_frontal_volume','L_frontal_volume_mixed','L_parietal_volume','L_parietal_volume_mixed','group_forgraph','all_diagnoses')
Graph_Regions = {'R_frontal_volume','R_frontal_volume_mixed','R_parietal_volume','R_parietal_volume_mixed','L_frontal_volume','L_frontal_volume_mixed','L_parietal_volume','L_parietal_volume_mixed'};

figure
set(gcf,'Position',[100 100 800 800]);
set(gcf,'color','w')
these_colors = jet(length(unique(group_forgraph)));
for i = 1:8
    subplot(4,2,i)
    %eval(['scatter(' Graph_Regions{i} ',MMN_amplitudes_forgraph,12,group_forgraph,''filled'')'])
    eval(['gscatter(' Graph_Regions{i} ',MMN_amplitudes_forgraph,group_forgraph,these_colors,''o'',6,''off'')'])
    hold on
    title(Graph_Regions{i},'interpreter','none')
    ylabel('MMN Amplitude')
    xlabel('Normalised volume')
    % Plot single line for overall, separate line for mised models
    if mod(i,2)==1
        coefficients = eval(['polyfit(' Graph_Regions{i} ', MMN_amplitudes_forgraph'', 1)']);
        % Create a new x axis with exactly 1000 points (or whatever you want).
        xFit = eval(['linspace(min(' Graph_Regions{i} '), max(' Graph_Regions{i} '), 1000)']);
        % Get the estimated yFit value for each of those 1000 new x locations.
        yFit = polyval(coefficients , xFit);
        % Plot everything.
        plot(xFit, yFit, 'k-', 'LineWidth', 2); % Plot fitted line.
    else
        for j = 1:length(unique(group_forgraph))
            coefficients = eval(['polyfit(' Graph_Regions{i} '(group_forgraph==j), MMN_amplitudes_forgraph(group_forgraph==j)'', 1)']);
            % Create a new x axis with exactly 1000 points (or whatever you want).
            xFit = eval(['linspace(min(' Graph_Regions{i} '(group_forgraph==j)), max(' Graph_Regions{i} '(group_forgraph==j)), 1000)']);
            % Get the estimated yFit value for each of those 1000 new x locations.
            yFit = polyval(coefficients , xFit);
            % Plot everything.
            plot(xFit, yFit, 'Color',these_colors(j,:), 'LineWidth', 2); % Plot fitted line.
        end
    end
    
end
legend(unique(all_diagnoses,'stable'))
drawnow
saveas(gcf,'Regional Correlations.png')
set(gcf, 'PaperPositionMode', 'auto');
saveas(gcf,'Regional Correlations.pdf')
