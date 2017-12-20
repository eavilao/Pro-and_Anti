function [nspk,timepoints] = Spiketimes2RateTrial(thisTrial_spk,timepoints,binwidth,analyse_sacc_win)

ntrls = length(thisTrial_spk);
timepoints = [timepoints(1)-binwidth timepoints timepoints(end)+binwidth];

if analyse_sacc_win
    [nspk,~] = hist(cell2mat({thisTrial_spk}'),timepoints);
else
    [nspk,~] = hist(cell2mat({thisTrial_spk.tspk_SS}'),timepoints);
end
% throw away histogram edges
nspk = nspk(2:end-1);
timepoints = timepoints(2:end-1);


