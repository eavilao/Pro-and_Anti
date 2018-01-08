function [out] = find_blinks(pupilsize)
% finds blinks based on the assumption that lighting is constant, and thus only minor pupil changes are expected (median + 1*std)

mp = median(pupilsize(~isnan(pupilsize)));
sp = std(pupilsize(~isnan(pupilsize)));

t = 1:length(pupilsize);
duringblinks = find(pupilsize < (mp-1*sp));
nanblinks = find(isnan(pupilsize));
duringblinks = sort([duringblinks;  nanblinks]);



D = duringblinks;

if D(1) == 1 % if trial starts with the eyes closed

	eyeopen =  [ D(find(diff([ D ; 0 ])>2))];
	eyeclose = [ D(find(diff([ 0 ; D])>2)) ];
else

	eyeopen =  [0 ; D(find(diff([ D ; 0 ])>2))];
	eyeclose = [ D(find(diff([ 0 ; D])>2)) ];

end


	out.open_intervals = [eyeopen eyeclose ];

plotblinks = 0;
if plotblinks
	figure
	hold on
	line(out.open_intervals', ones(size(out.open_intervals))'* (mp-1*sp)/100,'linewidth', 5)
	psz = pupilsize; psz = (psz/100); 
	plot(psz,'r')
	
end