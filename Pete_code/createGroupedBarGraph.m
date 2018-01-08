function [hPatches] = createGroupedBarGraph(xIn,erIn,hAx,cc,barWidth)
%function to create grouped bar graph, need 2 groups of 3 with individual
%colours 
%this is a modified version of createTafGroupedBarGraph but should now take
% different sizes of pairs(should really read groups) and bars per group
%the size of xIn and erIn should be the same and should be (numPairs,numBars) in size 

%xIn = reshape(nanmean(graphRelGains),2,[])';
%haX - axis handle


if isempty(hAx)
hF = figure; 
hAx = axes('parent',hF);
end
if isempty(cc)
cc = lines(12); 
end
%xIn is arranged in 2 columns of 4, each row should be plotted as a pair of
%bars 
%erIn is the same but error bars
numPairs = size(xIn,1);
numBars = size(xIn,2);
if isempty(barWidth)
barWidth = 0.15;
end
hatSize = 0.05;
%numPairs = 4;
hPatches = nan(numPairs,numBars);
barCount = 1;
for pairNum =1:numPairs
    
    %draw two patches for teh bars
    for barNum = 1:numBars
        %xIn(pairNum,barNum)
        startX = pairNum+barWidth*(barNum-2);
        endX = pairNum+barWidth*(barNum-1);
        hPatches(pairNum,barNum) = patch([endX endX startX startX],[xIn(pairNum,barNum) 0 0 xIn(pairNum,barNum)],ones(1,4),...
            'parent',hAx,'facecolor',cc(barCount,:));
        
        %draw vertical error bar on top
        %mean([startX endX])
        line([mean([startX endX]) mean([startX endX])],[xIn(pairNum,barNum)-erIn(pairNum,barNum) xIn(pairNum,barNum)+erIn(pairNum,barNum)],...
            'color','k','parent',hAx)
        %draw 'hat line on top off error bars
        line(repmat([mean([startX endX])-hatSize mean([startX endX])+hatSize],2,1)',repmat([xIn(pairNum,barNum)-erIn(pairNum,barNum); xIn(pairNum,barNum)+erIn(pairNum,barNum)],1,2)',...
            'color','k','parent',hAx)
        barCount = barCount+1;
    end
    
end
marg = 0.4;
set(hAx,'xlim',[marg (pairNum+1)-marg])