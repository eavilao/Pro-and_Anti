function plot_ProAnti_eye(units, plotType, cellNum, recArea)

switch plotType
    case 'eye_move'
        % LOAD FILE SENT BY NICO
        
        t = units(cellNum).pro.timepointsEye/1000; % time
        
        % extract eye trace around saccade onset
        % pro
        c = 1;
        for i = 1:length(units(cellNum).pro.behav.trial)
            sacc_onset = units(cellNum).pro.behav.trial(i).saccadeOnset;
            align_sacc = t - sacc_onset-0.012; t_sacc = align_sacc(align_sacc >= -0.150 & align_sacc <= 0.150);
            x = units(cellNum).pro.behav.trial(i).eyePositionX; y = units(cellNum).pro.behav.trial(i).eyePositionY;
            if ~isnan(x)
                heye(c,:) = x(align_sacc >= -0.150 & align_sacc <= 0.150);
                veye(c,:) = y(align_sacc >= -0.150 & align_sacc <= 0.150);
                c = c+1;
            end
        end
        
        figure; hold on;
        plot(t_sacc,heye','Color',[0.18 0.83 0.159], 'LineWidth', 2);
        plot(t_sacc,veye','Color',[0.117 0.44  1], 'LineWidth', 2);
        set(gca, 'yTick',[], 'TickDir', 'out', 'FontSize',18);
        title(['Eye trace pro cell ' num2str(cellNum)])
        
        % anti
        c = 1;
        for i = 1:length(units(cellNum).anti.behav.trial)
            sacc_onset = units(cellNum).anti.behav.trial(i).saccadeOnset;
            align_sacc = t - sacc_onset-0.012; t_sacc = align_sacc(align_sacc >= -0.150 & align_sacc <= 0.150);
            x = units(cellNum).anti.behav.trial(i).eyePositionX; y = units(cellNum).anti.behav.trial(i).eyePositionY;
            if ~isnan(x)
                heye(c,:) = x(align_sacc >= -0.150 & align_sacc <= 0.150);
                veye(c,:) = y(align_sacc >= -0.150 & align_sacc <= 0.150);
                c = c+1;
            end
        end
        
        figure; hold on;
        plot(t_sacc,heye','Color',[0.18 0.83 0.159], 'LineWidth', 2);
        plot(t_sacc,veye','Color',[0.117 0.44  1], 'LineWidth', 2);
        set(gca, 'yTick',[], 'TickDir', 'out', 'FontSize',18);
        title(['Eye trace anti cell ' num2str(cellNum)])
        
        % print('eye_trace','-depsc2', '-painters', '-cmyk')
        
end