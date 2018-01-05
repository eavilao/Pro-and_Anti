function [hRast] = createCombinedRaster(alignedRaster,timeWindow,varargin)
%function that takes an aligned raster (cell (numTrials,1)) with each row
%contianig spike times aligned to some timepoint which is now time 0(ouput
%of createRaster).
%Will plot the resulting raster and a smoothed sdf (if flagged)
p = inputParser;
p.addParamValue('spikeMarkSize',3,@(x) isnumeric(x));
p.addParamValue('showKsDensity',1,@(x) isnumeric(x));
p.parse(varargin{:});
spikeMarkSize = p.Results.spikeMarkSize;
showKsDensity = p.Results.showKsDensity;

numTriggers = size(alignedRaster,1);
%create figure
hF = figure('name','Raster Plot'); 
hRast = axes('Parent',hF);
xlabel(hRast,'Time (s)')
ylabel(hRast,'Trial Number')
set(hRast,'Xlim',[-timeWindow(1) timeWindow(2)] ,'Ylim',[1 max(numTriggers,2)+0.5])

%line to mark 0
line([0 0],[1 numTriggers],'Color','g','LineWidth',2,'parent',hRast)

for triggerNum = 1:numTriggers
    trialSpikeTimes = alignedRaster{triggerNum};
    line(trialSpikeTimes,triggerNum*ones(1,length(trialSpikeTimes)),...
        'color','k','linestyle','none','marker','s','markerfacecolor','k',...
        'markersize',spikeMarkSize)
end
    
if showKsDensity 
    hApsth = axes('XAxisLocation','bottom','YAxisLocation','right','color','none',...
    'Xcolor',[1 1 1],'Ycolor',[1 0 0],'XTick',[],'parent',hF,'position',get(hRast,'position')); %,'XLim',[-15 15],'YLim',[-15 15]
    allSpikes = vertcat(alignedRaster{:});
    timeVector = [-timeWindow(1):0.001:timeWindow(2)];
binSize = 0.01;
smoothFiringRate = ksdensity(allSpikes,timeVector,'width',binSize)*(length(allSpikes)/numTriggers);
line(timeVector,smoothFiringRate,'parent',hApsth,'color','r','linewidth',2)
linkaxes([hRast hApsth],'x')
end