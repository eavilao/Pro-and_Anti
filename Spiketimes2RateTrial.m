function [nspk,timepoints] = Spiketimes2RateTrial(thisTrial_spk,timepoints,binwidth,analyse_sacc_win,id)

ntrls = length(thisTrial_spk);
timepoints = [timepoints(1)-binwidth timepoints timepoints(end)+binwidth];

if strcmp(id,'SS') % either SS or CS
    
    if analyse_sacc_win
        [nspk,~] = hist(cell2mat({thisTrial_spk}'),timepoints);
    else
        [nspk,~] = hist(cell2mat({thisTrial_spk.tspk_SS}'),timepoints);
         nspk(nspk>1)=1;
    end
else
    
    if analyse_sacc_win
        [nspk,~] = hist(cell2mat({thisTrial_spk}'),timepoints);
    else
        [nspk,~] = hist(cell2mat({thisTrial_spk.tspk_CS}'),timepoints);
        nspk(nspk>1)=1;
    end

end

% throw away histogram edges
nspk = nspk(2:end-1);
timepoints = timepoints(2:end-1);



