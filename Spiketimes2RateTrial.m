function [nspk,timepoints] = Spiketimes2RateTrial(thisTrial_spk,timepoints,binwidth,analyse_sacc_align,id)

ntrls = length(thisTrial_spk);
timepoints = [timepoints(1)-binwidth timepoints timepoints(end)+binwidth];

if strcmp(id,'SS') % either SS or CS
    
    if analyse_sacc_align
        [nspk,~] = hist(cell2mat({thisTrial_spk}'),timepoints);
       
    else
        [nspk,~] = hist(cell2mat({thisTrial_spk}'),timepoints);
        
    end
else
    
    if analyse_sacc_align
        [nspk,~] = hist(cell2mat({thisTrial_spk}'),timepoints);
        
    else
        [nspk,~] = hist(cell2mat({thisTrial_spk}'),timepoints);
      
    end

end

% throw away histogram edges
nspk = nspk(2:end-1);
timepoints = timepoints(2:end-1);

% rate 
nspk = sum(nspk)/1.1; 
