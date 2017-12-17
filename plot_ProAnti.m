function plot_ProAnti(trialData, i, plottype)

% Needed: trialData.mat

switch plottype
    
    case 'raster'

%% Raster

% aligned to ??
figure; hold on; 
for i=1:length(trialData(1).trial)
    if ~isempty(trialData(1).trial(i).tspk_SS)
        plot(trialData(1).trial(i).saccadeOnset,i,'.r');
        plot(trialData(1).trial(i).tspk_SS(1:5:end),i,'.k');
    end
end
xlim([-1 1]);

% saccade aligned
figure; hold on; 
for j=1:length(trialData(i).trial)
    if ~isempty(trialData(i).trial(j).tspk_SS)
        % plot(trialData(1).trial(i).tspk_SS(1:5:end) - trialData(1).trial(i).saccadeOnset,i,'.k');
        plot(trialData(i).trial(j).tspk_SS - trialData(i).trial(j).saccadeOnset,j,'.k');
    end
end
xlim([-0.5 0.5]);
vline(0) 


%% PSTH 
    case 'psth'


%gather
t= trialData(1).pro.ts; 
r_pro= trialData(i).pro.nspk - ;
r_anti = trialData(i).anti.nspk; 

% smooth 
r_pro_smooth = smooth(r_pro,50);
r_anti_smooth = smooth(r_anti,50);
% Plot single units aligned to saccade
figure; 
plot(t, r_pro_smooth, 'LineWidth', 2); hold on
plot(t, r_anti_smooth, 'LineWidth', 2); 
xlim([-1 1]); 
  

% Try smoothing smoothFiringRate = ksdensity(newTrialStruct(trialNum).alignedSpikes{1},testTimeVector,'width',binSize);

%% colormaps


end
