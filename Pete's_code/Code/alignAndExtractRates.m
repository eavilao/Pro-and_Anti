function [alignedTrials muRate] = alignAndExtractRates(trialStructure,selectedTrials,alignmentPoint,varargin)
%function to get a raster from eric' data realign each trial and then
%calculate a spike density functino for each row and in total.

p = inputParser;
p.addParamValue('unitNum',1,@(x) isnumeric(x));
p.addParamValue('timeWindow',[4 4],@(x) isnumeric(x));
p.addParamValue('doPlot',0,@(x) isnumeric(x));

p.addParamValue('binSize',0.01,@(x) isnumeric(x));
p.parse(varargin{:});
unitNum = p.Results.unitNum;
timeWindow = p.Results.timeWindow;
doPlot = p.Results.doPlot;
binSize = p.Results.binSize;


theseStableTrials = trialStructure(selectedTrials);
alignedTrials = theseStableTrials;
numTrials = size(alignedTrials,2);
if strcmpi(alignmentPoint,'fixation')
    %if fixation is chosen then no realignemnt required
    alignmentTimes = cell(1,numTrials);
    alignmentTimes(:) = {0};
else
    alignmentTimes = {theseStableTrials.(alignmentPoint)};
end


for trialNum = 1:numTrials
    alignedTrials(trialNum).alignedSpikes = cellfun(@(x,y) x-y,theseStableTrials(trialNum).alignedSpikes,repmat(alignmentTimes(trialNum),1,size(theseStableTrials(trialNum).alignedSpikes,2)),'uniformoutput',false);
    %alignedTrials(trialNum).alignedSpikes = cellfun(@(x,y) x-y,theseStableTrials(trialNum).alignedSpikes,repmat(alignmentTimes(trialNum),1,2),'uniformoutput',false);
end



test = {alignedTrials.alignedSpikes};
alignedRaster = cellfun(@(x) x(unitNum),test);
allSpikes = vertcat(alignedRaster{:});



timeVector = [-timeWindow(1):0.001:timeWindow(2)];

if ~isempty(allSpikes)
    smoothFiringRate = ksdensity(allSpikes,timeVector,'width',binSize)*(length(allSpikes)/numTrials);
else
    smoothFiringRate = zeros(1,length(timeVector));
end
%ksdensity has error with emptys so only run on non empty the fill rows
%with zeros of right size

rateCell = repmat({zeros(1,length(timeVector))},1,numTrials);
fSdf = @(x) ksdensity(x,timeVector,'width',binSize)*length(x);
trialsWithSpikes =  find(~cellfun('isempty',alignedRaster));
rateCell(trialsWithSpikes) = cellfun(@(x) fSdf(x),alignedRaster(trialsWithSpikes),'uniformoutput',false);


rateMat = cell2mat(rateCell')';
muRate = mean(rateMat,2);
for trialNum = 1:numTrials
    alignedTrials(trialNum).sdf = rateCell{trialNum};
    alignedTrials(trialNum).sdfT = timeVector;
    
    
end
if doPlot
    figure;
    hA2 = axes;
    plot(timeVector,rateMat,'k')
    line(timeVector,smoothFiringRate,'parent',hA2,'color','r','linewidth',2)
    
    
end
