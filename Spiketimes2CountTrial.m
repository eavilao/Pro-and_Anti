function [nspk,timepoints] = Spiketimes2CountTrial(thisTrial_spk,timepoints,binwidth,analyse_sacc_align,analyse_instrDir_align,id)

ntrls = length(thisTrial_spk);
timepoints = [timepoints(1)-binwidth timepoints timepoints(end)+binwidth];

if strcmp(id,'SS') % either SS or CS
    
    if analyse_sacc_align
        [nspk,~] = hist(cell2mat({thisTrial_spk.tspk_SS_align_sacc}'),timepoints);
        nspk(nspk>1)=1;
    elseif analyse_instrDir_align
        [nspk,~] = hist(cell2mat({thisTrial_spk.tspk_SS_align_instrDir}'),timepoints);
        nspk(nspk>1)=1;
    else
        [nspk,~] = hist(cell2mat({thisTrial_spk.tspk_SS}'),timepoints);
        nspk(nspk>1)=1;
    end
else
    if analyse_sacc_align
        [nspk,~] = hist(cell2mat({thisTrial_spk.tspk_CS_align_sacc}'),timepoints);
        nspk(nspk>1)=1;
    else
        [nspk,~] = hist(cell2mat({thisTrial_spk.tspk_CS}'),timepoints);
        nspk(nspk>1)=1;
    end

end

% throw away histogram edges
nspk = nspk(2:end-1);
timepoints = timepoints(2:end-1);



