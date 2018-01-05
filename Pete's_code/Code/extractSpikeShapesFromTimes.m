function [spikeShapes timeVector] = extractSpikeShapesFromTimes(data,spikeTimes,eventWindow,fs,varargin)
%extract the spike shapes for each of the given times and returna s celll
%array
%
%first varargin is 
p = inputParser;
p.addParamValue('cc',lines(12),@(x) isnumeric(x));
p.addParamValue('doPlotFlag',0,@(x) isnumeric(x));
p.addParamValue('exSpikeShapeNum',[],@(x) isnumeric(x));
p.addParamValue('hAx',[],@(x) isnumeric(x));
p.parse(varargin{:});
cc = p.Results.cc; 
doPlotFlag = p.Results.doPlotFlag;
exSpikeShapeNum = p.Results.exSpikeShapeNum;
hAx = p.Results.hAx;


data = data(:); %ensure column
if iscell(spikeTimes);
    activeUnits  = find(~cellfun('isempty',spikeTimes));
else
    activeUnits = 1;
    spikeTimes = {spikeTimes};
end


eventWindowSamp = [floor(eventWindow(1)*fs/1000) ceil(eventWindow(2)*fs/1000)]; %convert window to samples
eventLength = sum(eventWindowSamp)+1;  

spikeShapes = cell(1,max(activeUnits));
%loop over units and get shapes for each
for unitNum = activeUnits
    theseSpikeTimes = spikeTimes{unitNum};

     eventTimesSamp = round(theseSpikeTimes*fs); %convert triggerTimes to samples
    
    %remove those too close to end or beginning
    
    eventTimesSamp(eventTimesSamp>length(data)-eventWindowSamp(2)-1) = [];
    eventTimesSamp(eventTimesSamp<eventWindowSamp(1)+1) =[];
    
    %loop over and extract shapes
    numEvents = length(eventTimesSamp);
    
    theseSpikeShapes = zeros(eventLength,numEvents);
    for eventNum = 1:numEvents
        theseSpikeShapes(:,eventNum) = data(eventTimesSamp(eventNum)-eventWindowSamp(1):eventTimesSamp(eventNum)+eventWindowSamp(2));
    end
    
    spikeShapes{unitNum} = theseSpikeShapes;
end
timeVector = linspace(-eventWindowSamp(1)/fs,eventWindowSamp(2)/fs,eventLength);
%plot if asked
if doPlotFlag
    %now do plot
    if isempty(hAx)
        hF = figure('color','w');
        hAx = axes('parent',hF);
        
    end
    for unitNum = fliplr(activeUnits)
        
        if isempty(exSpikeShapeNum)
            theseSpikes = spikeShapes{unitNum};
        else
            
            theseSpikes = spikeShapes{unitNum};
            if exSpikeShapeNum<size(theseSpikes,2)
                spikesToShow = randperm(size(theseSpikes,2),exSpikeShapeNum);
            else
                spikesToShow = 1:size(theseSpikes,2);
            end
            theseSpikes = theseSpikes(:,spikesToShow);
        end
        
        line(timeVector,theseSpikes,'color',cc(unitNum,:),'parent',hAx)
        
        
    end
    allUnitLabels = {'unit1','unit2','unit3','unit4','unit5','unit6'};
    hL = createLegend(cc(activeUnits,:),allUnitLabels(activeUnits),repmat({'-'},1,length(activeUnits)),hAx,[0.7 0.7 0.2 0.2]);
end
end