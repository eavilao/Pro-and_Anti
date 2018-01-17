function normData = normToBaseInds(data,baseInd)
%data should be trials x samples, inds specify wherre to take the mean and
%sd from, performs a z-score style normalisation trial by trail of input
%matrix
%is set up to ignore nans
  baseVals = nanmean(data(:,baseInd(1):baseInd(2)),2);
  baseSds = nanstd(data(:,baseInd(1):baseInd(2)),0,2);
  meanCorrected = bsxfun(@minus,data,baseVals);
  normData = bsxfun(@rdivide,meanCorrected,baseSds);
  