function [nspk,timepoints] = Spiketimes2Rate(thisTrial_spk,timepoints,binwidth,analyse_sacc_win,id)

ntrls = length(thisTrial_spk);
timepoints = [timepoints(1)-binwidth timepoints timepoints(end)+binwidth];

if strcmp(id,'SS') % either SS or CS
    
    if analyse_sacc_win
        [nspk,~] = hist(cell2mat({thisTrial_spk.tspk_SS_align_sacc}'),timepoints);
    else
        [nspk,~] = hist(cell2mat({thisTrial_spk.tspk_SS}'),timepoints);
    end
    
else
    if analyse_sacc_win
        [nspk,~] = hist(cell2mat({thisTrial_spk.tspk_CS_align_sacc}'),timepoints);
    else
        [nspk,~] = hist(cell2mat({thisTrial_spk.tspk_CS}'),timepoints);
    end
    
end


% throw away histogram edges
nspk = nspk(2:end-1);
timepoints = timepoints(2:end-1);

% trial-average firing rates in units of spikes/s
nspk = nspk/(ntrls*binwidth);


