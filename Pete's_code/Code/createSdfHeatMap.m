function [hA] = createSdfHeatMap(trialStructure,sortByString,hA,varargin)
%function to create a plot which represents the amplitude of the sdf in
%each row with a bar whos coulpur at each points represents the value of
%teh sdf a tthat time point
%TODO enable passing fixed scale as input arg
p = inputParser;
%p.addParamValue('displayFlag',0,@(x) isnumeric(x));
p.addParamValue('numBins',0,@(x) isnumeric(x));
p.addParamValue('frRange',[-inf inf],@(x) isnumeric(x));

p.parse(varargin{:});
%displayFlag = p.Results.displayFlag;
numBins = p.Results.numBins;
frRange = p.Results.frRange;


cc = jet(100);
if isempty(hA)
    hF = figure;
    hA = axes('parent',hF,'units','normalized','position',[0.05 0.1 0.7 0.8]);
    hA2 = axes('parent',hF,'units','normalized',...
        'position',[0.75 0.1 0.2 0.8],'yaxislocation','right','ydir','reverse');
end



%should also allow sorting by some feature (field of input struct)and
%binning in theis field as well

rateMat = cell2mat({trialStructure.sdf}');
[val srtOrder] = sort([trialStructure.(sortByString)]);

rateMat = rateMat(srtOrder,:);
trialLength = size(rateMat,2); 
if numBins==0
    %if no binning required
    
    hI = imagesc(rateMat,'parent',hA);
    set(hI,'alphadata',~isnan(rateMat));
    
    line(val,1:length(srtOrder),'parent',hA2)
    set(get(hA2,'xlabel'),'string',sortByString)
    set(hA2,'ylim',[0 max(length(srtOrder),1)])
    
else
    
    %if binning needed then divide sortBy measure into correct numBins
    
    binEdges = linspace(min(val),max(val),numBins);
    
    %now for each bin go through and get mean of all trials within that bin
    binnedResults = cell(numBins,1);
    binnnedVals = cell(numBins,1);
    binnedCounts = zeros(numBins,1);
    for binNum = 1:numBins-1
        
        theseTrials = find(val>=binEdges(binNum) & val<binEdges(binNum+1));
        if ~isempty(theseTrials)
            
            if length(theseTrials)>1
                binnedResults{binNum} = mean(rateMat(theseTrials,:));
                binnnedVals{binNum} = mean(val(theseTrials));
                 binnedCounts(binNum) = length(theseTrials);
            else
                binnedResults{binNum} = rateMat(theseTrials,:);
                binnnedVals{binNum} = val(theseTrials);
                binnedCounts(binNum) = 1;
            end
            
            
        else
            binnedResults{binNum} = nan(1,trialLength);
            binnnedVals{binNum} = NaN;
        end
    end
   
    
    
    hI = imagesc(cell2mat(binnedResults),'parent',hA);
    set(hI,'alphadata',~isnan(cell2mat(binnedResults)));
    
    
    %line(cell2mat(binnnedVals),1:numBins-1,'parent',hA2)
    %replave drawing a straight line of bin vlues with the cumsum of bins
    line(cumsum(binnedCounts),1:length(cumsum(binnedCounts)),'parent',hA2)
    
    set(get(hA,'ylabel'),'string',sortByString)
    set(get(hA2,'xlabel'),'string','Num Trials')
    %set(hA2,'ylim',[0 length(srtOrder)])
    axTickMark = [min(binEdges) mean(binEdges) max(binEdges)];
   
  %set(hA,'ylim',[axTickMark(1) axTickMark(3)])
   
    set(hA2,'ytick',[0 numBins/2 numBins],'yticklabel',cellstr(num2str(axTickMark')))
end

%show scale

set(hA,'xlim',[1750 2250])
set(hA,'xtick',[1750:250:2250],'xticklabel',{'-0.5','0','0.5'})

%set the scale if necessary
if ~isinf(frRange(2))
    set(hA,'clim',frRange)
    
end
colorbar('peer',hA)

