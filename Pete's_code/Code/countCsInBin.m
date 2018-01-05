function [csCount] = countCsInBin(allTrials,trialsToAnalyse,unitNum,countPeriod)
%very simple function to take Eric's trials structure and count the number
%of cs in the bin provided. Should return a vector same length as
%trialsToAnalyse with CS count per trial.


%cat all trials we want
allAlignedSpikes = {allTrials(trialsToAnalyse).alignedSpikes};
csAligned = cellfun(@(x) x(unitNum),allAlignedSpikes);
%now 
saccadePeriodCs = cellfun(@(x) x(x>countPeriod(1) & x<countPeriod (2)),csAligned,'uniformOutput',false);
csCount = cellfun(@(x) numel(x),saccadePeriodCs,'uniformOutput',false);
csCount = cell2mat(csCount);