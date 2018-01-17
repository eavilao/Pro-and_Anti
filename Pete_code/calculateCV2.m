function [CV2vector] = calculateCV2(spikeTimes)
%function to calculate the CV2 of a spike train, usually the median of this
%vector will be taken as the final result
fCV2 = @(x,y) 2*abs(x-y)/(x+y);
if size(spikeTimes,1)>=3
    theseISIs = diff(spikeTimes);
    
    
    
    isiMat = [circshift(theseISIs,1) theseISIs];
    
    
    CV2vector = arrayfun(@(x,y) fCV2(x,y),isiMat(:,1),isiMat(:,2));
    CV2vector(1) = nan; %first spike cv2 is not real
    %cv2Med = nanmedian(cv2Vec);
    
    %thisNeuronResults(slidingStepNum,trialNum) = cv2Med;
else
    CV2vector = [];
    
end