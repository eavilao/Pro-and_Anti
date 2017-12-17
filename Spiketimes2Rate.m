function [nspk,timepoints] = Spiketimes2Rate(thisTrial_spk,timepoints,binwidth)

ntrls = length(thisTrial_spk);
timepoints = [timepoints(1)-binwidth timepoints timepoints(end)+binwidth];

[nspk,~] = hist(cell2mat({thisTrial_spk.tspk_SS}'),timepoints);

% throw away histogram edges
nspk = nspk(2:end-1);
timepoints = timepoints(2:end-1);
% trial-average firing rates in units of spikes/s
nspk = nspk/(ntrls*binwidth);
