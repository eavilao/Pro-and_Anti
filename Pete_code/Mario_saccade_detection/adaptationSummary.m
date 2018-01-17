function out = adaptationSummary(varargin)


p = inputParser;
p.addRequired('sacs')  % a matrix with two columns or a cell array with two cells;
p.addRequired('data')  % a matrix with two columns or a cell array with two cells;
p.addParamValue('trialintervals', 0)
p.addParamValue('allinone', 1) % stdandard deviation criterion for offset threshold

p.parse(varargin{:});

sacs = p.Results.sacs;
data = p.Results.data;
allinone = p.Results.allinone;
trialintervals = p.Results.trialintervals;

% trialintervals = 
lm = sacs.saccades.landmarks.v; 
gazeerror = data.target - data.eyexy;
xy = data.eyexy;

vprofiles = sacs.saccades.vProfiles.v;

notail = find(lm(:,6)==0); 
withtail = find(lm(:,6)~=0);


sacend_t = lm(:,6);
sacend_t(notail) = lm(notail,4);

% onset gaze error
gazeerror_on(:,1) = gazeerror(int32(lm(:,1)),1);
gazeerror_on(:,2) = gazeerror(int32(lm(:,1)),2);

% offset gaze error
gazeerror_off(:,1) = gazeerror(int32(sacend_t),1);
gazeerror_off(:,2) = gazeerror(int32(sacend_t),2);

gazeerror_off(:,1) = gazeerror(int32(sacend_t),1);
gazeerror_off(:,2) = gazeerror(int32(sacend_t),2);

% offset gaze error at velocity threshold crossing (similar to maximum deceleration)
gazeerror_off_early(:,1) = gazeerror(int32(lm(:,4)),1);
gazeerror_off_early(:,2) = gazeerror(int32(lm(:,4)),2);

% error from the target for saccades with or without tails
gazeerror_sac_end = hypot(gazeerror_off(:,1), gazeerror_off(:,2));
gazeerror_sac_early_end = hypot(gazeerror_off_early(:,1), gazeerror_off_early(:,2));

angle_on = sacs.saccades.angle.onset;
	

% gaze position (two column matrix with x and y of onset eye positions and offset positions)
pos_on(:,1) =  xy(int32(lm(:,1)),1);
pos_on(:,2) =  xy(int32(lm(:,1)),2);
pos_off(:,1) = xy(int32(lm(:,4)),1);
pos_off(:,2) = xy(int32(lm(:,4)),2);
	
pos_late_off = pos_off;
pos_late_off(withtail,1) = xy(int32(lm(withtail,6)),1);
pos_late_off(withtail,2) = xy(int32(lm(withtail,6)),2);


trialintervals = data.intervals;
sacs_in_trial = find_sacs_in_trial(lm,trialintervals);

primary_sacs = find(sacs_in_trial(:,2)==1);
secondary_sacs = find(sacs_in_trial(:,2)==2);


% amplitude of saccade at onset and offset, assuming cue is at 0,0 deg origin

amplitude_primary = sacs.saccades.amplitude(secondary_sacs - 1, 1);
amplitude_secondary = sacs.saccades.amplitude(secondary_sacs  , 1);	


% saccade indices where the position onset is according to the following coordinate intervals
oncue  = find(hypot(pos_on(:,1), pos_on(:,2)) <=3);
offcue = find(hypot(pos_on(:,1), pos_on(:,2)) >=3);

trials = 1:length(gazeerror_sac_early_end)';

trials_oncue_primary = trials(intersect(primary_sacs,oncue));
trials_oncue_secondary = trials(secondary_sacs);
	
	
amp_on = hypot(pos_on(:,1),pos_on(:,2));
amp_off = hypot(pos_off(:,1),pos_off(:,2));
amp_late_off = hypot(pos_late_off(:,1),pos_late_off(:,2));

X = gazeerror_sac_end;

fwdwin = 20;
framerange = [1:length(gazeerror_sac_end)-fwdwin];

qlevels = [.1 .9]
wk = [10];

for frame = framerange

			accum_pkg = amp_on(frame:frame+fwdwin);
			accum_pkg = accum_pkg(:);

			q(frame,:) = quantile(accum_pkg,qlevels);
			g(frame,:) = ksdensity(accum_pkg,[0:.1:20],'width',wk);


			
end




% trials = data.intervals(:,1);





% 	   _   __    _____                           
% 	  (_)_/_/   / __(_)___ ___  __________  _____
% 	   _/_/    / /_/ / __ `/ / / / ___/ _ \/ ___/
% 	 _/_/_    / __/ / /_/ / /_/ / /  /  __(__  ) 
% 	/_/ (_)  /_/ /_/\__, /\__,_/_/   \___/____/  
% 	               /____/                        

if allinone
%	figure
%	axespositions = [
%         0.13          0.8         0.84         0.13;
%         0.13         0.67         0.84         0.13;
%         0.13         0.21         0.19         0.14;
%         0.13         0.02         0.19         0.14;
%         0.35         0.03         0.16         0.16;
%         0.35         0.21         0.16         0.13;
%         0.13         0.52         0.83         0.14;
%  		 0.13         0.38         0.83         0.13
%         ];
    

%	a(1) = axes('position',axespositions(1,:));

	figure;
	ax(1) = subplot(211);
	line([trials(trials_oncue_primary) ; trials(trials_oncue_primary)], [amp_on(trials_oncue_primary)  amp_off(trials_oncue_primary)]','linestyle','-', 'color', [.5 0 0])
		line([trials(trials_oncue_primary) ; trials(trials_oncue_primary)], [amp_off(trials_oncue_primary)  amp_late_off(trials_oncue_primary)]','linestyle','-', 'color', [ 0 .5 0],'linewidth',2)
		line(trials(trials_oncue_primary), amp_late_off(trials_oncue_primary),'linestyle','none', 'color', [.5 0 0],'markersize',10,'marker','o')

	
	title('Primary saccade and glissade amplitude')
	ylabel('amplitude of saccade in degrees')
	xlabel('trial')

	ax(2) = subplot(212);
	line([trials(trials_oncue_secondary) ; trials(trials_oncue_secondary)], [amp_on(trials_oncue_secondary)  amp_off(trials_oncue_secondary)]','linestyle','-', 'color', [.5 0 0])
	line([trials(trials_oncue_secondary) ; trials(trials_oncue_secondary)], [amp_off(trials_oncue_secondary)  amp_late_off(trials_oncue_secondary)]','linestyle','-', 'color', [ 0 .5 0],'linewidth',2)
	line(trials(trials_oncue_secondary), amp_late_off(trials_oncue_secondary),'linestyle','none', 'color', [.5 0 0],'markersize',10,'marker','o')

	title('Secondary saccade and glissade amplitude')
	ylabel('amplitude')

	linkaxes(ax,'x')


%	a(3) = axes('position',axespositions(3,:));
	
	% line(trials(oncue), amp_late_off(oncue),'linestyle','none', 'color', [0 0 0],'markersize',10,'marker','^')
	% line(intersect(offcue, withtail),amp_off(intersect(offcue, withtail)),'marker', '.' ,'linestyle','none')

figure;	
	line([trials(trials_oncue_primary) ; trials(trials_oncue_primary)], [zeros(size(trials_oncue_primary));  (amp_off(trials_oncue_primary)- amp_on(trials_oncue_primary))'] )
	

	% title('vector amplitude for primary')

	% line([pos_on(:,1) ; pos_off(:,1)]', [pos_on(:,2) ; pos_off(:,2)]' )
	
	% line(pos_off(oncue,1), pos_off(oncue,2),'linestyle','none','markersize', 15,'marker' ,'.' )
	
	% axis tight




%	a(4)  = axes('position',axespositions(4,:));
figure;
	rose(angle_on(oncue))
	title('saccade counts per direction')
	
%	a(5) = axes('position',axespositions(5,:));
figure;
	    muV = mean(vprofiles)';
        stdV = std(vprofiles)';

        
        plot(vprofiles','color', [.7 .7 .7])
        hold on
        plot(muV,'color','r','linewidth',2)
        plot(muV+stdV,'color', 'r')
        plot(muV-stdV,'color', 'r')


        title('5')



%	a(6) = axes('position',axespositions(6,:))	

figure;
	plot([pos_on(oncue,1) pos_off(oncue,1)]', [pos_on(oncue,2) pos_off(oncue,2) ]' )

	% subplot(223)
	% line(gazeerror_off_early(notail,1), gazeerror_off(notail,1),'linestyle','none')
	% subplot(224)
	% line(gazeerror_off_early(withtail,1), gazeerror_off(notail,1),'linestyle','none')

	corr([pos_on(oncue,1) pos_off(oncue,1) pos_on(oncue,2) pos_off(oncue,2) ] )
	title('6')


%	a(7) = axes('position',axespositions(7,:));
figure;
	X=  [trials(oncue) ; trials(oncue)];
	Y = [zeros(size(trials(oncue))) ; (amp_late_off(oncue)-amp_off(oncue))'];
	line(X, Y,'linestyle','-', 'color', [0 .5 0],'linewidth',2)
	title('glisssade amplitudes')


%	a(8) = axes('position',axespositions(8,:));
figure;	
	line([trials(oncue) ; trials(oncue)], [zeros(size(amp_off(oncue))) amp_on(oncue)  ]','linestyle','-', 'color', [ 0 0 .5])
	title('cue to onset position')


	% line([trials ; trials], [pos_on pos_off]','linestyle','-', 'color', [.5 0 0])
	% line(trials, pos_off,'linestyle','none', 'color', [.5 0 0],'markersize',10,'marker','^')
	% hold on
	% plot(trials, amp_off,'.y')
	% plot(notail,amp_off(notail),'.b')

	% plot([pos_on(:,1) ; pos_off(:,1)]', [pos_off(:,1) ; pos_off(:,2)]' )
	% there seems to be a correlation between onset position and offset trajectory?

figure;
	
line(trials(trials_oncue_primary),(amp_off(trials_oncue_primary))', 'linestyle','none', 'color', [0 0 1],'markersize',10,'marker','o')



else

	figcount = 1;

	fig(figcount) = figure; figcount = figcount +1;


	line(amplitude_primary, amplitude_secondary,'linestyle','none','marker','o' )

	fig(figcount) = figure; figcount = figcount +1;


	line(amplitude_primary, amplitude_secondary,'linestyle','none','marker','o' )

% 	   _   __                _                          __      __ 
	% 	  (_)_/_/   ____ ___  __(_)   _____  _____   ____  / /___  / /_
	% 	   _/_/    / __ `/ / / / / | / / _ \/ ___/  / __ \/ / __ \/ __/
	% 	 _/_/_    / /_/ / /_/ / /| |/ /  __/ /     / /_/ / / /_/ / /_  
	% 	/_/ (_)   \__, /\__,_/_/ |___/\___/_/     / .___/_/\____/\__/  
	% 	            /_/                          /_/                   	% 	fig(figcount) = figure; figcount = figcount +1;
	% 	% -1 trick in the index to select primary saccades that had a corrective saccade
	% 	quiver(x, y, u, v);

	% 	x = [1:length(secondary_sacs-1)];
	% 	y = zeros(1,length(secondary_sacs-1));
	% 	u = (pos_on(secondary_sacs,2) - pos_on(secondary_sacs,1) )' ;
	% 	v = (pos_off(secondary_sacs,2) - pos_off(secondary_sacs,1))' ;
	% quiver(x, y, u, v);

% 	   _   __               __           _ __                                     __      __  _                __         __                              _____           __                     __                                  __
	% 	  (_)_/_/   _   _____  / /___  _____(_) /___  __   _________  _____________  / /___ _/ /_(_)___  ____     / /_  ___  / /__      _____  ___  ____     / __(_)_________/ /_   ____ _____  ____/ /  ________  _________  ____  ____/ /
	% 	   _/_/    | | / / _ \/ / __ \/ ___/ / __/ / / /  / ___/ __ \/ ___/ ___/ _ \/ / __ `/ __/ / __ \/ __ \   / __ \/ _ \/ __/ | /| / / _ \/ _ \/ __ \   / /_/ / ___/ ___/ __/  / __ `/ __ \/ __  /  / ___/ _ \/ ___/ __ \/ __ \/ __  / 
	% 	 _/_/_     | |/ /  __/ / /_/ / /__/ / /_/ /_/ /  / /__/ /_/ / /  / /  /  __/ / /_/ / /_/ / /_/ / / / /  / /_/ /  __/ /_ | |/ |/ /  __/  __/ / / /  / __/ / /  (__  ) /_   / /_/ / / / / /_/ /  (__  )  __/ /__/ /_/ / / / / /_/ /  
	% 	/_/ (_)    |___/\___/_/\____/\___/_/\__/\__, /   \___/\____/_/  /_/   \___/_/\__,_/\__/_/\____/_/ /_/  /_.___/\___/\__/ |__/|__/\___/\___/_/ /_/  /_/ /_/_/  /____/\__/   \__,_/_/ /_/\__,_/  /____/\___/\___/\____/_/ /_/\__,_/   
	% 	                                       /____/                                                                                                                                                                                      
	% fig(figcount) = figure; figcount = figcount +1;

	% V = sacs.saccades.peakV.v;
	% line(V(secondary_sacs-1),V(secondary_sacs),'linestyle','none','marker','o' )


		fig(figcount) = figure; figcount = figcount +1;	
		% line([1:length(q)],q,'linewidth',5)

		line([trials(trials_oncue_primary) ; trials(trials_oncue_primary)], [amp_on(trials_oncue_primary)  amp_off(trials_oncue_primary)]','linestyle','-', 'color', [.5 0 0])
		line([trials(trials_oncue_primary) ; trials(trials_oncue_primary)], [amp_off(trials_oncue_primary)  amp_late_off(trials_oncue_primary)]','linestyle','-', 'color', [ 0 .5 0],'linewidth',2)
		line(trials(trials_oncue_primary), amp_late_off(trials_oncue_primary),'linestyle','none', 'color', [.5 0 0],'markersize',10,'marker','o')

		
		figure
		line([trials(trials_oncue_secondary) ; trials(trials_oncue_secondary)], [amp_on(trials_oncue_secondary)  amp_off(trials_oncue_secondary)]','linestyle','-', 'color', [.5 0 0])
		line([trials(trials_oncue_secondary) ; trials(trials_oncue_secondary)], [amp_off(trials_oncue_secondary)  amp_late_off(trials_oncue_secondary)]','linestyle','-', 'color', [ 0 .5 0],'linewidth',2)
		line(trials(trials_oncue_secondary), amp_late_off(trials_oncue_secondary),'linestyle','none', 'color', [.5 0 0],'markersize',10,'marker','o')

		

		axis tight
		title('saccade and glissade amplitude')
		ylabel('amplitude of saccade in degrees')
		xlabel('trial')





	% line([1:length(q)],q,'linewidth',5)
	% line([trials ; trials], [gazeerror_sac_early_end  gazeerror_sac_end]','linestyle','-', 'color', [.5 0 0])
	% line(trials, gazeerror_sac_early_end,'linestyle','none', 'color', [.5 0 0],'markersize',10,'marker','^')
	% hold on
	% plot(trials, gazeerror_sac_early_end,'.y')
	% plot(notail,gazeerror_sac_end(notail),'.b')


	% fig(figcount) = figure; figcount = figcount +1;
	% line([trials(oncue) ; trials(oncue)], [gazeerror_sac_early_end(oncue) gazeerror_sac_end(oncue)]','linestyle','-', 'color', [.5 0 0])
	% line(trials(oncue), gazeerror_sac_early_end(oncue),'linestyle','none', 'color', [.5 0 0],'markersize',10,'marker','^')
	% hold on
	% plot(trials(oncue), gazeerror_sac_early_end(oncue),'.y')
	% plot(intersect(notail, oncue),gazeerror_sac_end(intersect(notail,oncue)),'.b')

	% fig(figcount) = figure; figcount = figcount +1;
	% line([trials(offcue) ; trials(offcue)], [gazeerror_sac_early_end(offcue) gazeerror_sac_end(offcue)]','linestyle','-', 'color', [.5 0 0])
	% line(trials(offcue), gazeerror_sac_early_end(offcue),'linestyle','none', 'color', [.5 0 0],'markersize',10,'marker','^')
	% hold on
	% plot(trials(offcue), gazeerror_sac_early_end(offcue),'.y')
	% plot(intersect(notail, offcue),gazeerror_sac_end(intersect(notail,offcue)),'.b')


	fig(figcount) = figure; figcount = figcount +1;

	line([1:length(q)],q,'linewidth',5)
	line([trials ; trials], [amp_on  amp_off]','linestyle','-', 'color', [.5 0 0])
	line([trials ; trials], [amp_off  amp_late_off]','linestyle','-', 'color', [ 0 .5 0],'linewidth',2)
	line(trials, amp_late_off,'linestyle','none', 'color', [.5 0 0],'markersize',10,'marker','^')
	hold on
	plot(trials, amp_late_off,'.y')
	plot(notail,amp_off(notail),'.b')
	plot(offcue, amp_late_off(offcue),'marker','x','color','g','markersize',10,'linestyle','none')
	axis tight
	title('saccade and glissade amplitude')
	ylabel('amplitude of saccade in degrees')
	xlabel('trial')



	fig(figcount) = figure; figcount = figcount +1;

	line([1:length(q)],q,'linewidth',2)
	line([trials(oncue) ; trials(oncue)], [amp_on(oncue)  amp_off(oncue)]','linestyle','-', 'color', [ 0 0 .5])
	line([trials(oncue) ; trials(oncue)], [amp_off(oncue)  amp_late_off(oncue)]','linestyle','-', 'color', [0 .5 0],'linewidth',2)

	line(trials(oncue), amp_late_off(oncue),'linestyle','none', 'color', [.5 0 0],'markersize',10,'marker','^')
	plot(trials(oncue), amp_late_off(oncue),'.y')
	plot(intersect(oncue, notail),amp_off(intersect(oncue,notail)),'.b')
	axis tight


	hold on
	plot(trials(oncue), amp_late_off(oncue),'.y')



	fig(figcount) = figure; figcount = figcount +1;
	subplot(221)
	line([pos_on(:,1) ; pos_off(:,1)]', [pos_on(:,2) ; pos_off(:,2)]' )
	line( pos_off(:,1), pos_off(:,2),'linestyle','none','markersize', 15,'marker' ,'.' )

	subplot(222)
	line([pos_on(:,1) ; pos_off(:,1)]', [pos_on(:,2) ; pos_off(:,2)]' )
	line( pos_off(:,1), pos_off(:,2),'linestyle','none','markersize', 15,'marker' ,'.' )

	subplot(223)
	plot([pos_on(oncue,1) pos_off(oncue,1)]', [pos_on(oncue,2) pos_off(oncue,2) ]' )

	subplot(224)
	plot([pos_on(offcue,1) pos_off(offcue,1)]', [pos_on(offcue,2) pos_off(offcue,2) ]' )

	% subplot(223)
	% line(gazeerror_off_early(notail,1), gazeerror_off(notail,1),'linestyle','none')
	% subplot(224)
	% line(gazeerror_off_early(withtail,1), gazeerror_off(notail,1),'linestyle','none')

	corr([pos_on(oncue,1) pos_off(oncue,1) pos_on(oncue,2) pos_off(oncue,2) ] )


	fig(figcount) = figure; figcount = figcount +1;
	subplot(211)
	X=  [trials(oncue) ; trials(oncue)];
	Y = [zeros(size(trials(oncue))) ; (amp_late_off(oncue)-amp_off(oncue))'];
	line(X, Y,'linestyle','-', 'color', [0 .5 0],'linewidth',2)
	title('glisssade amplitudes')

	subplot(212)
	line([trials(oncue) ; trials(oncue)], [zeros(size(amp_off(oncue))) amp_on(oncue)  ]','linestyle','-', 'color', [ 0 0 .5])
	title('cue to onset position')


	% line([trials ; trials], [pos_on pos_off]','linestyle','-', 'color', [.5 0 0])
	% line(trials, pos_off,'linestyle','none', 'color', [.5 0 0],'markersize',10,'marker','^')
	% hold on
	% plot(trials, amp_off,'.y')
	% plot(notail,amp_off(notail),'.b')

	% plot([pos_on(:,1) ; pos_off(:,1)]', [pos_off(:,1) ; pos_off(:,2)]' )
	% there seems to be a correlation between onset position and offset trajectory?


end


	out.offcue = offcue;
	out.oncue = oncue;
	out.position_saccade_onset = pos_on;
	out.position_saccade_early_offset = pos_off;
	out.position_saccade_late_offset = pos_late_off;
	out.gain = 'disp not yet computed'
	out.sac_ordered_in_trial = sacs_in_trial;
	out.number_of_corrective_saccades = length(secondary_sacs);
	



function out = find_sacs_in_trial(lm, intervals )

	sac_in_trial = @ (c)intervals( find(lm(c,1) > intervals(:,1) & lm(c,4) <= intervals(:,5),1,'first'),1)	;
	trial_of_sac = arrayfun(sac_in_trial, [1:size(lm,1)])';

	sac_order_in_trial(1) = 1;
	for t = 2:length(trial_of_sac)
		if trial_of_sac(t-1) == trial_of_sac(t);
			sac_order_in_trial(t) = sac_order_in_trial(t-1)+1;
		else
			sac_order_in_trial(t) = 1;
		end
	end


out = [trial_of_sac sac_order_in_trial'];













