function [MI confLims sigLevel] = miAndConfLims(labels,data,montCarlIterations)
%function to calculate the mutual information provided about labels by data
%(sould be column vectors of teh same size. Then performs montecarlo style
%resampling to get confidence limits on distribution of MIs given by
%chance. 


MI = mutualinfo(labels,data);

numLabels = length(labels);
montCarlRes = nan(montCarlIterations,1);
%tic
for MontCarlIterNum = 1:montCarlIterations
    shuffledLabels = labels(randperm(numLabels));
    montCarlRes(MontCarlIterNum) = mutualinfo(shuffledLabels,data);
end
%toc
confLims = prctile(montCarlRes,[2.5 97.5]);



[val ord] = sort([montCarlRes; MI]);
posInList = find(ord==montCarlIterations+1);
sigLevel = posInList/(montCarlIterations+1);

sigLevel = 0.5-abs(0.5-sigLevel);
% hF = figure;
% hAx = axes('parent',hF);
% [counts,centers] = hist(montCarlRes,100);
% bar(centers,counts,'parent',hAx)
% 
% line([confLims(1) confLims(1)],[0 100],'yliminclude','off','parent',hAx,'color','r')
% line([confLims(2) confLims(2)],[0 100],'yliminclude','off','parent',hAx,'color','r')
% line([MI MI],[0 100],'yliminclude','off','parent',hAx,'color','g')
% set(get(hAx,'title'),'string',[num2str(sigLevel)])
