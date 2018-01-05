function [hAx1 hAx2 hA2] = rasterFromTrials(trialStructure,trialsToPlot,timeWindow,origTimeWindow,varargin)
%function to take eric's trial structure and plot all the rasters from each
%trial on top of each other as well as the eye movement
%origTimeWindow is needed to know the oringial limits of teh eyeMovmeent
%epochs
%TODO should also be able to realign to diffrent time point
ccSpikes = [0 0 0;1 0 0]; %spike colourmap
spikeMarkerSizes = [2 6]; %different size markers fo the different spikes
spikeMarkers = {'s','o'};
ccEvents = [1 0 1;0 1 0]; %events colourmap

eventNum = 1; %TODO for multiple event markers


%use dynamic field names to switch event types
% sortBy = 'saccadeTime';
% alignTo = 'goCueTime';
% markEvent = 'saccadeTime';

p = inputParser;
p.addParamValue('sortBy','goCueTime',@(x) ischar(x));
p.addParamValue('alignTo','saccadeTime',@(x) ischar(x));
p.addParamValue('markEvent','goCueTime',@(x) ischar(x));
p.addParamValue('burstMarkingFlag',[],@(x) isnumeric(x));
p.addParamValue('pauseMarkingFlag',[],@(x) isnumeric(x));%turns on marking of pauses within trials, separate colours for cs and non cs
p.addParamValue('spikeMarkingFlag',1,@(x) isnumeric(x));
p.parse(varargin{:});
sortBy = p.Results.sortBy;
alignTo = p.Results.alignTo;
markEvent = p.Results.markEvent;
burstMarkingFlag = p.Results.burstMarkingFlag;
pauseMarkingFlag = p.Results.pauseMarkingFlag;
spikeMarkingFlag = p.Results.spikeMarkingFlag;
% sortBy = 'goCueTime';
% alignTo = 'saccadeTime';
% markEvent = 'goCueTime';
%trialStructure().bit6time = 0;
[trialStructure.bit6time] = deal(0);

%sorting by event time must happen after realignment
alignmentTimes = {trialStructure(trialsToPlot).(alignTo)};
% indCounter=1;
% tempCell = {trialStructure(indCounter).alignedSpikes};
% while cellfun('isempty',tempCell)
% tempCell = {trialStructure(indCounter).alignedSpikes}; %need to find the numbers of the acftiveUnits by findng non empty cells
% indCounter =indCounter+1;
% end
% activeUnits = find(~cellfun('isempty',tempCell{:}));
%method above didnt work if the first trials had no complex spikes, but
%other did.
tempCell = vertcat(trialStructure.alignedSpikes);
activeUnits = find(sum(~cellfun('isempty',tempCell))>1);

numTrials = length(trialsToPlot);
stInit = cell(numTrials,1); %TODO should be numUnits
alignedTrials = struct('alignedSpikes',stInit,'trialInFile',stInit,...
    'goCueTime',stInit,'saccadeTime',stInit,'eyePositionX',stInit,...
    'eyePositionY',stInit,'bit6time',stInit,'burstTime',stInit)
for trialNum = 1:numTrials %TODO replace with cellfun?
    
    alignedTrials(trialNum).alignedSpikes = cellfun(@(x,y) x-y,trialStructure(trialsToPlot(trialNum)).alignedSpikes,repmat(alignmentTimes(trialNum),1,size(trialStructure(trialsToPlot(trialNum)).alignedSpikes,2)),'uniformoutput',false)
    alignedTrials(trialNum).trialInFile = trialsToPlot(trialNum);
    %alignedTrials(trialNum).bit2time = trialStructure(trialsToPlot(trialNum)).bit2time - alignmentTimes{trialNum}; %0.1 is to change go cue to relate to 
    %alignedTrials(trialNum).bit2time = trialStructure(trialsToPlot(trialNum)).bit3time - alignmentTimes{trialNum}; %0.1 is to change go cue to relate to 
    alignedTrials(trialNum).goCueTime = trialStructure(trialsToPlot(trialNum)).goCueTime - alignmentTimes{trialNum}; %0.1 is to change go cue to relate to 
    alignedTrials(trialNum).saccadeTime = trialStructure(trialsToPlot(trialNum)).saccadeTime - alignmentTimes{trialNum};
    alignedTrials(trialNum).timeShift = -alignmentTimes{trialNum}; %TODO enable realignemnt later and for eye traces
    alignedTrials(trialNum).bit6time = -alignmentTimes{trialNum};
    
    alignedTrials(trialNum).eyePositionX = trialStructure(trialsToPlot(trialNum)).eyePositionX;
    alignedTrials(trialNum).eyePositionY = trialStructure(trialsToPlot(trialNum)).eyePositionY;
    
%     %also do for bursts
%     if ~isempty(burstMarkingFlag)
%         
%         %alignedTrials(trialNum).burstToMark =
%         %trialStructure(trialsToPlot(trialNum)).incBurstStruct
%         
%         if ~isnan(trialStructure(trialsToPlot(trialNum)).overlapFr)
%             burstToUse = 1;
%             if size(trialStructure(trialsToPlot(trialNum)).incBurstStruct,2)>1
%                 %TODO currently we only mark the biggest burst as indicated
%                 %bythe highest surprise index, in future it should take all
%                 %bursts
%                 [val burstToUse] = max([trialStructure(trialsToPlot(trialNum)).incBurstStruct.SI]);
%                 
%             end
%                 
%             alignedTrials(trialNum).burstTime =  [trialStructure(trialsToPlot(trialNum)).incBurstStruct(burstToUse).startTime, trialStructure(trialsToPlot(trialNum)).incBurstStruct(burstToUse).endTime] - (trialStructure(trialsToPlot(trialNum)).trialStart + alignmentTimes{trialNum});
%        
%            
%         end
%     end
if ~isempty(pauseMarkingFlag)

    
         %trialStructure(trialsToPlot(trialNum)).incBurstStruct
    
    if ~isempty(trialStructure(trialsToPlot(trialNum)).relativePauseInfo)


        alignedTrials(trialNum).relativePauseInfo = trialStructure(trialsToPlot(trialNum)).relativePauseInfo;
        
                numPauses = numel(trialStructure(trialsToPlot(trialNum)).relativePauseInfo);
        for pauseNum = 1:numPauses
        alignedTrials(trialNum).relativePauseInfo(pauseNum).pauseStart = alignedTrials(trialNum).relativePauseInfo(pauseNum).pauseStart-alignmentTimes{trialNum};
        alignedTrials(trialNum).relativePauseInfo(pauseNum).pauseEnd = alignedTrials(trialNum).relativePauseInfo(pauseNum).pauseEnd-alignmentTimes{trialNum};
       
        end
        %[alignedTrials(trialNum).relativePauseInfo.pauseStart] = [alignedTrials(trialNum).relativePauseInfo.pauseStart]-alignmentTimes{trialNum}
    else
        
    alignedTrials(trialNum).relativePauseInfo = trialStructure(trialsToPlot(trialNum)).relativePauseInfo;
    end
end
end


%sort trials by chosen event
%[b byEventOrder] = sort([trialStructure(trialsToPlot).(sortBy)]); 
%trialPlotOrder = trialsToPlot(byEventOrder);
if ~isempty(sortBy)
[b, byEventOrder] = sort([alignedTrials.(sortBy)]);
else
    byEventOrder = 1:numTrials;
end
trialPlotOrder = byEventOrder;



hF = figure;
hAx1 = subplot(2,1,1,'parent',hF); %for eye traces
hAx2 = subplot(2,1,2,'parent',hF);%for raster
set([hAx1 hAx2],'xlim',[-timeWindow(1) timeWindow(2)])


set(hAx2,'ylim',[0 numTrials])

hPatches = cell(numTrials,1);
%loop over trials
for trialNum = 1:numTrials
    
    
    %thisTrial = trialStructure(trialPlotOrder(trialNum));
    thisTrial = alignedTrials(trialPlotOrder(trialNum));
    if ~isnan(alignmentTimes{trialPlotOrder(trialNum)})
    
    %plot both eye traces %TODO also plot target position
    eyeT = linspace(-origTimeWindow(1),origTimeWindow(2),length(thisTrial.eyePositionX));
    eyeT = eyeT+thisTrial.timeShift;
    line(eyeT,thisTrial.eyePositionX,'color','r','linewidth',2,'parent',hAx1)
    line(eyeT,thisTrial.eyePositionY,'color','b','linewidth',2,'parent',hAx1)
    
    if spikeMarkingFlag
    for unitNum = 1:length(activeUnits)
    thisUnitNum = activeUnits(unitNum);
        %plot spikes in correct row
        theseSpikes = thisTrial.alignedSpikes{thisUnitNum};
        
        line(theseSpikes,trialNum*ones(1,length(theseSpikes)),'parent',hAx2,...
            'color',ccSpikes(unitNum,:),'linestyle','none','marker',spikeMarkers{unitNum},'MarkerFaceColor',ccSpikes(unitNum,:),...
            'markersize',spikeMarkerSizes(unitNum),'tag',['tMark' num2str(unitNum)])
    end
    end
    %also plot selected events
    
    thisEvent = thisTrial.(markEvent);
    line(thisEvent,trialNum*ones(1,length(thisEvent)),'parent',hAx2,...
        'color',ccEvents(eventNum,:),'linestyle','none','marker','s','MarkerFaceColor',ccEvents(eventNum,:),...
        'markersize',2)
      
    
    %also show bursts as light green patches
    if ~isempty(pauseMarkingFlag) && ~isempty(thisTrial.relativePauseInfo)
    numPauses = numel(thisTrial.relativePauseInfo);
    thesePatches = nan(numPauses,1);
    for pauseNum = 1:numPauses
%         hPatches(trialNum) =  patch([thisTrial.burstTime(2) thisTrial.burstTime(2) thisTrial.burstTime(1) thisTrial.burstTime(1)],...
%             [trialNum+0.25 trialNum-0.25 trialNum-0.25 trialNum+0.25],...
%             [0 1 0],'parent',hAx2,'edgecolor','none','facealpha',0.5)
if thisTrial.relativePauseInfo(pauseNum).csFlag
        thesePatches(pauseNum) =  patch([thisTrial.relativePauseInfo(pauseNum).pauseEnd thisTrial.relativePauseInfo(pauseNum).pauseEnd thisTrial.relativePauseInfo(pauseNum).pauseStart thisTrial.relativePauseInfo(pauseNum).pauseStart],...
            [trialNum+0.25 trialNum-0.25 trialNum-0.25 trialNum+0.25],...
            [0 0 1],'parent',hAx2,'edgecolor','none','facealpha',1);
else
     thesePatches(pauseNum) =  patch([thisTrial.relativePauseInfo(pauseNum).pauseEnd thisTrial.relativePauseInfo(pauseNum).pauseEnd thisTrial.relativePauseInfo(pauseNum).pauseStart thisTrial.relativePauseInfo(pauseNum).pauseStart],...
            [trialNum+0.25 trialNum-0.25 trialNum-0.25 trialNum+0.25],...
            [0 1 0],'parent',hAx2,'edgecolor','none','facealpha',1);
end
    end
    hPatches{trialNum} = thesePatches;
    end
    end
end
hAllPatches = vertcat(hPatches{:});
uistack(hAllPatches(~isnan(hAllPatches)),'bottom')
%add smoothed psth for unit 1
%place the complex spike marks on top
uistack(findobj(hAx2,'tag','tMark2'),'top')

hA2 = axes('XAxisLocation','bottom','YAxisLocation','right','color','none',...
    'Xcolor',[1 1 1],'Ycolor',[1 0 0],'XTick',[],'parent',hF,'position',get(hAx2,'position')); %,'XLim',[-15 15],'YLim',[-15 15]
%set(get(hA2,'YLabel'),'string','Firing Rate (Hz)')

thisUnitNum = activeUnits(1); %for simple spikes
test = {alignedTrials.alignedSpikes};
alignedRaster = cellfun(@(x) x(thisUnitNum),test);
allSpikes = vertcat(alignedRaster{:});

timeVector = [-timeWindow(1):0.001:timeWindow(2)];
binSize = 0.01;
smoothFiringRate = ksdensity(allSpikes,timeVector,'width',binSize)*(length(allSpikes)/numTrials);

line(timeVector,smoothFiringRate,'parent',hA2,'color',[0 0.7 0],'linewidth',2)
linkaxes([hAx1 hAx2 hA2],'x')
%set([hAx1 hAx2 hA2],'tag','keep')