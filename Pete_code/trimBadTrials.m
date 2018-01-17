function [trialsToKeep rejectionCodes percentTrialsRemoved] = trimBadTrials(trialsToTrim,trimCriteria)
%removes bad trials from erics data based on the fields of trimCriteria
%where the field name is the name of the field in the trialStrucutre to
%test and the value is [max min] allowed for that parameter.

%find what we want to test
rejectionCriteria = fieldnames(trimCriteria);

numCriteria = size(rejectionCriteria,1);
numTrials = size(trialsToTrim,2);

rejectionCodes = zeros(numTrials,numCriteria);
%loop over each 
criteriumCounter = 1;
for thisCriterium = rejectionCriteria'
    
    thisCriteriumLow = trimCriteria.(char(thisCriterium))(1);
    thisCriteriumHigh = trimCriteria.(char(thisCriterium))(2);
    %get value of thsi for all trials
    trialValues = [trialsToTrim.(char(thisCriterium))];
    %find those outside our limits
   rejectionCodes(:,criteriumCounter) = trialValues<thisCriteriumLow | trialValues>thisCriteriumHigh;
    criteriumCounter = criteriumCounter+1;
end

%now find all trials which violated at least one
trialsToRemove = find(sum(rejectionCodes,2));

numTrialsRemoved = length(trialsToRemove);

trialsToKeep = trialsToTrim;
trialsToKeep(trialsToRemove) = [];

percentTrialsRemoved = 100*(numTrialsRemoved/numTrials);