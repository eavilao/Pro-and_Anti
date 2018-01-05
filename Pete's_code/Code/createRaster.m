function [hRast alignedRaster] = createRaster(spikeTimes,triggerTimes,windowSize,varargin)
%fucntion taken from old createPsth.m but now just produces the raster and
%returns the matrix and a handle to the axis

p = inputParser;
p.addParamValue('eventTimes',[],@(x) isnumeric(x));
p.addParamValue('complexSpikes',[],@(x) isnumeric(x));
p.addParamValue('spikeMarkSize',3,@(x) isnumeric(x));
p.addParamValue('showKsDensity',0,@(x) isnumeric(x));
p.addParamValue('axHandle',[],@(x) ishandle(x));
p.parse(varargin{:});
eventTimes = p.Results.eventTimes; %[start(i) end(i);start(ii) end(ii)] matrix of toher events to mark as patches on raster
spikeMarkSize = p.Results.spikeMarkSize;
showKsDensity = p.Results.showKsDensity;
axHandle = p.Results.axHandle;
complexSpikes = p.Results.complexSpikes;
%ensure column vectors
spikeTimes = spikeTimes(:);
triggerTimes = triggerTimes(:);

plotTitle = 'test';


calcWindowSize = windowSize*1.2; %for calculation to stop end point errors
numTriggers = length(triggerTimes); 

rasterCell = cell(numTriggers,1);

for triggerNum = 1:numTriggers
     rasterCell{triggerNum} = spikeTimes(spikeTimes>(triggerTimes(triggerNum)-calcWindowSize(1)) & spikeTimes<(triggerTimes(triggerNum)+calcWindowSize(2)));
%         if ~isempty(eventTimes)
%             eventCell{triggerNum} = eventTimes(eventTimes>(triggerTimes(triggerNum)-calcWindowSize(1)) & eventTimes<(triggerTimes(triggerNum)+calcWindowSize(2)));
%         end    
end    
alignedRaster = cellfun(@minus,rasterCell,num2cell(triggerTimes),'UniformOutput',0);

if ~isempty(complexSpikes)
    complexSpikes =complexSpikes(:);
    complexSpikeRaster = cell(numTriggers,1);
    for triggerNum = 1:numTriggers
        complexSpikeRaster{triggerNum} = complexSpikes(complexSpikes>(triggerTimes(triggerNum)-calcWindowSize(1)) & complexSpikes<(triggerTimes(triggerNum)+calcWindowSize(2)));
        %         if ~isempty(eventTimes)
        %             eventCell{triggerNum} = eventTimes(eventTimes>(triggerTimes(triggerNum)-calcWindowSize(1)) & eventTimes<(triggerTimes(triggerNum)+calcWindowSize(2)));
        %         end
    end
   complexSpikeRaster = cellfun(@minus,complexSpikeRaster,num2cell(triggerTimes),'UniformOutput',0);
end

if ~isempty(eventTimes)
    %eventTimes = eventTimes(:);
    eventCell = cell(numTriggers,1);
    %each cel of eventCell should have a copy of all events aligned to this
    %trigger
    %TODO in future trim events to save memory
    for triggerNum = 1:numTriggers
            %eventCell{triggerNum} = eventTimes(eventTimes>(triggerTimes(triggerNum)-calcWindowSize(1)) & eventTimes<(triggerTimes(triggerNum)+calcWindowSize(2)));
%eventTimes(:,1)>(triggerTimes(triggerNum)-calcWindowSize(1)) | eventTimes(:,2)<(triggerTimes(triggerNum)+calcWindowSize(2))
        eventCell{triggerNum} = eventTimes;
    end 
    alignedEvents = cellfun(@minus,eventCell,num2cell(triggerTimes),'UniformOutput',0);
end

if isempty(axHandle)
hF = figure('name',['Raster Plot ' plotTitle]); 
hRast = axes('Parent',hF);
else
hRast = axHandle;
end
xlabel(hRast,'Time (s)')
ylabel(hRast,'Trial Number')
set(hRast,'Xlim',[-windowSize(1) windowSize(2)] ,'Ylim',[1 max(numTriggers,2)+0.5])

%line to mark 0
line([0 0],[1 numTriggers],'Color','g','LineWidth',2,'parent',hRast)
for triggerNum = 1:numTriggers
    trialSpikeTimes = alignedRaster{triggerNum}; %rasterCell{triggerNum} - (triggerTimes(triggerNum));
    %numTrialSpikes = length(trialSpikeTimes);
    %     for trialSpikesNum = 1:numTrialSpikes
    %         rectangle('Position',[trialSpikeTimes(trialSpikesNum) triggerNum 0.001 1],...
    %             'Parent',hRast,'FaceColor','k','EdgeColor','k');
    %
    %
    %     end
    line(trialSpikeTimes,triggerNum*ones(1,length(trialSpikeTimes)),...
        'color','k','linestyle','none','marker','s','markerfacecolor','k',...
        'markersize',spikeMarkSize)
    %old method for events specified with single number
    %     if ~isempty(eventTimes)
    %         trialEventTimes = alignedEvents{triggerNum};
    %         if isempty(trialEventTimes)
    %             for trialEventNum = 1:length(trialEventTimes)
    %                 rectangle('Position',[trialEventTimes(trialEventNum) triggerNum-0.5 0.001 1],...
    %                     'Parent',hRast,'FaceColor','r','EdgeColor','r');
    %             end
    %         end
    %     end
    
    if ~isempty(eventTimes)
        theseEvents = alignedEvents{triggerNum};
        numTheseEvents = size(alignedEvents{triggerNum},1);
        %now draw a patch i light red showing the start and end of each event
        for eventNum = 1:numTheseEvents
            
            patch([theseEvents(eventNum,2) theseEvents(eventNum,2) theseEvents(eventNum,1) theseEvents(eventNum,1)],...
                [triggerNum+0.5 triggerNum-0.5 triggerNum-0.5 triggerNum+0.5],...
                [1 0 0],'parent',hRast,'edgecolor','none')
        end
    end
    
    
    
    if ~isempty(complexSpikes)
        trialComplexSpikeTimes = complexSpikeRaster{triggerNum}; 
        line(trialComplexSpikeTimes,triggerNum*ones(1,length(trialComplexSpikeTimes)),...
        'color','r','linestyle','none','marker','o','markerfacecolor','r',...
        'markersize',spikeMarkSize*2)
        
    end
end

%if flagged then calculate the smooth firing rate and plot over the top on
%second axis on right hand side
if showKsDensity 
    hApsth = axes('XAxisLocation','bottom','YAxisLocation','right','color','none',...
    'Xcolor',[1 1 1],'Ycolor',[1 0 0],'XTick',[],'parent',hF,'position',get(hRast,'position')); %,'XLim',[-15 15],'YLim',[-15 15]
    allSpikes = vertcat(alignedRaster{:});
    timeVector = [-windowSize(1):0.001:windowSize(2)];
binSize = 0.01;
smoothFiringRate = ksdensity(allSpikes,timeVector,'width',binSize)*(length(allSpikes)/numTriggers);
line(timeVector,smoothFiringRate,'parent',hApsth,'color',[1 0 0],'linewidth',2)
linkaxes([hRast hApsth],'x')
end