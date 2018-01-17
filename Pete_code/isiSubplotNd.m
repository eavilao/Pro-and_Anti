function hF = isiSubplotNd(spikeTimes,varargin)
% fucntion for creatng nxn grid of subplots of isi histograms
% on diagonal and with raster triggered unit 1v unit2 in corners
%INPUTS:
%   spikeTimes - cell (1,numUnits) with times of each spike in s
%       eg. {[1 2 3 4 5],[0.1 0.5 3.5]} is a set of two units with 5 spikes
%       in the first at the given times and 3 in the second unit at given
%       times
%   Param/Value pairs
%           'rastWinWidth' - 1,2 double - time before and after trigger to
%           include in raster 
%               eg [0.2 0.2] 200ms before and after (default)
%OUTPUTS:
%   hF - handle to the figure


p = inputParser;
p.addParamValue('rastWinWidth',[0.2 0.2],@(x) isnumeric(x));
p.parse(varargin{:});


rastWinWidth = p.Results.rastWinWidth; 
%%
hF = figure;

numUnits = length(spikeTimes);
numPlots = numUnits^2;
%TODO switch to subaxis for more control
%create matrix with subplot style indexing
temp = zeros(numUnits,numUnits);
temp(:) = 1:numUnits^2;
temp  = temp'; %to match subplot indexing

%loop over all plots
for plotNum =1:numPlots
    
    hTemp = subplot(numUnits,numUnits,plotNum,'parent',hF);
    %find which triggger and response units are required for this plot
    [triggerUnit respUnit] = ind2sub([numUnits numUnits],temp(plotNum));
    
    if respUnit == triggerUnit %ismember(plotNum,identLine)
        %do isi histogram if  unit x v unit x
        isiVec = diff(spikeTimes{respUnit});
        [isiHist bins] = hist(isiVec,500);
        isiHist = 100*isiHist/length(isiVec);
        bar(bins,isiHist,'parent',hTemp);
        xlabel('time (s)')
        ylabel('precentage spikes')
        title(['Unit ' num2str(respUnit)])
    else
        %do raster if unit x v unit y
        unitSpikes = spikeTimes{respUnit};
        triggerSpikes = spikeTimes{triggerUnit};
        [hRast alignedRaster] = createRaster(unitSpikes,triggerSpikes,rastWinWidth);
        test = copyobj(hRast,hF);  
        delete(get(hRast,'parent'))
        set(test,'position',get(hTemp,'position'))
        delete(hTemp)
        hTemp = test;
    end   
end

%add some labels
hA = axes('units','normalized','position',[0.02 0.4 0.05 0.2],'visible','off');
text('Parent',hA,'string','Trigger Unit','rotation',90)
hA = axes('units','normalized','position',[0.4 0.95 0.2 0.05],'visible','off');
text('Parent',hA,'string','Response Unit','rotation',0)

end