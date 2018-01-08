%% find_crossings: function description
function [ups, downs, pauses] = find_crossings(aev, threshold, min_pause)
    ups = find((aev(1:end-1)<threshold) & (aev(2:end)>=threshold))-1;
    downs = find((aev(1:end-1)>threshold) & (aev(2:end)<=threshold))+1;
    

keyboard


    if numel(ups)~=numel(downs)
        warning('Breaking ups and downs do not match.');
    end
    
    if ups(1) >= downs(1)
        warning('The first breaking down happens befor the first breaking up');
    end

    
    pauses = ups(2:end)-downs(1:end-1);
    N = numel(pauses);
    marks = find(pauses<=min_pause);
    ups(marks+1) = [];
    downs(marks) = [];
    pauses = ups(2:end)-downs(1:end-1);




end