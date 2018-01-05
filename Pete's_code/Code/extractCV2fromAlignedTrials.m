function [cv2mat hFC hAxH trialStructure] = extractCV2fromAlignedTrials(trialStructure,timeVector,unitNum,sortByString)
%function that takes a (xnumTrials) structure containing the field
%alignedSpikes which has the times of spikes in seconds aligned to (time 0)
%some point
%TODO also add a flag for figure generation
hFC = figure;
hAxC = subplot(2,1,1,'parent',hFC);
hAxC2 = subplot(2,1,2,'parent',hFC);

theseTrials = trialStructure;
numTrials = size(theseTrials,2);

cv2mat = nan(length(timeVector),numTrials);
for trialNum = 1:numTrials
    disp(num2str(trialNum))
    theseSpikes = theseTrials(trialNum).alignedSpikes(unitNum);
    theseSpikes = [theseSpikes{:}];
    theseSpikes = unique(theseSpikes); %this line was addded as there wa an error with one file with 2 spikes at the same time
    thisCV2 = calculateCV2(theseSpikes);
    
    line(theseSpikes(1:end-1),thisCV2,'parent',hAxC)
    
    %cv2 is calcualted per spike, later on when we want
    %to compare across neurons we need to interp to a 1ms
    %grid or soemthing
    %for now I just want to plot the CV2 over time for
    %every trial
    resampledCV2 = interp1(theseSpikes(2:end-1),thisCV2(2:end),timeVector);
    
    line(timeVector,resampledCV2,'parent',hAxC2)
    cv2mat(:,trialNum) = resampledCV2;
    
    trialStructure(trialNum).cv2 = resampledCV2;
    trialStructure(trialNum).cv2T = timeVector;
end

%now reorder the trials in sorted order
[val srtOrder] = sort([trialStructure.(sortByString)]);
cv2mat = cv2mat(:,srtOrder);

%now also draw the mean of the CV2matrix
line(timeVector,nanmean(cv2mat,2),'color','k','linewidth',2,'parent',hAxC2)
set([hAxC hAxC],'xlim',[-1 1])
cv2mat = cv2mat';
hFh = figure;
hAxH = axes('parent',hFh,'ydir','reverse')
%TODO now try drawing a heat map style plot
ccHeat = jet(100);
hI = imagesc(cv2mat,'parent',hAxH);
set(hI,'alphadata',~isnan(cv2mat));
set(hAxH,'xlim',[3000 5000])