function [out selected] = filterSaccadeStruct(varargin)



p = inputParser;
p.addRequired('sacs')  % the saccade structure as output from saccadeDetection
p.addParamValue('selection', 'good') % stdandard deviation criterion for offset threshold
p.addParamValue('problemstoremove', [1:17]) %[0 0 0 0 0 0 1 1 1   1 1 1 1 0]
p.addParamValue('validintervals', [0 Inf]) % 
p.addParamValue('minvel', 0)
p.addParamValue('minamp', 0)
p.addParamValue('resdidualthreshold', 0.02)

p.parse(varargin{:});

in = p.Results.sacs;
selection = p.Results.selection;
probtorem = p.Results.problemstoremove;
minamp = p.Results.minamp;
minvel = p.Results.minvel;
resthr = p.Results.resdidualthreshold;
valint = p.Results.validintervals;

filterproblematic = 1;

% transform problem indices in binary yes/no selection vector for problems
% to remove
ptr = zeros(17,1);

if ~isempty(probtorem)
    ptr(probtorem)  = 1;
end

% problems to remove, marked with 1

		% 1 riseonset = 0
		% 2 minimum velocity not reached
		% 3 velocity larger than maximum accepted
		% 4 falloff very slow
		% 5 mean preceding velocity too large
		% 6 nan's in window
		% 7 too many peaks (>100deg/s) within detection window
		% 8 too small risetime
		% 9 saccade with negative duration?
		
		% 10 pospk of velocity - prewin <= 0 
		% 11 problematic onset: acceleration inflexion (to positive) with velocity larger than threshold onset 
		% 12 large motion before onset, using only low velocity samples
		% 13 no clear end of the main sequence: time to maximum deceleration larger than rise time
		% 14 no clear end of the saccade: velocity always larger than offset (raised threshold)



toosmallvel			= 	find(in.saccades.problematic(:,1)*ptr(1)	)';
landmarkmissing		= 	find(in.saccades.problematic(:,2)*ptr(2)	)';
toohighvel  		= 	find(in.saccades.problematic(:,3)*ptr(3)	)';
falloffslow	  		= 	find(in.saccades.problematic(:,4)*ptr(4)	)';
highvelpre  		= 	find(in.saccades.problematic(:,5)*ptr(5)	)';
nansinwin			= 	find(in.saccades.problematic(:,6)*ptr(6)	)';
problem8			= 	find(in.saccades.problematic(:,8)*ptr(8)	)';

brokenprepk 		= 	find(in.saccades.problematic(:,10)*ptr(10)	)'; %== 12
brokenonset 		= 	find(in.saccades.problematic(:,11)*ptr(11)	)'; %== 12
mainseqbroken 		= 	find(in.saccades.problematic(:,12)*ptr(12)	)'; %== 12
noend 				= 	find(in.saccades.problematic(:,14)*ptr(14)	)'; %== 12

toosmallamp 		= 	find(in.saccades.amplitude(:,1)< minamp )'*ptr(15)    ; % smaller than D degrees
hugeresiduals 		= 	find(in.saccades.fits.residual >= resthr );





accepted = [];
for s = 1:length(in.saccades.landmarks.v(:,1))

    S = in.saccades.landmarks.v(s,1);
    for i = 1:size(valint,1)
        if S >= valint(i,1) & S <= valint(i,2)
            accepted = [accepted s];
            continue
        end
    end
end

            
bad = unique([toosmallvel toohighvel falloffslow landmarkmissing highvelpre brokenprepk mainseqbroken  brokenonset toosmallamp problem8 nansinwin noend hugeresiduals]);
good = setdiff([1: length(in.saccades.amplitude)], bad);




good = intersect( good, accepted);

switch selection
	case 'good'
		selected = good;
	case 'bad'
		selected = bad;
	case 'nores'
		selected = setdiff([1: length(in.saccades.amplitude)], hugeresiduals);
end

disp([selection ' saccades selected: ' num2str(length(selected)) ' from ' num2str(length(in.saccades.amplitude(:,1)))])

% Hack for bug in saccadeDetection
selected(end) = [];

out.description = {[selection ' saccades selected: ' num2str(length(selected)) ' from ' num2str(length(in.saccades.amplitude(:,1)))]};

out.saccades.landmarks.v 			= in.saccades.landmarks.v 				(selected,:);
out.saccades.amplitude	    		= in.saccades.amplitude	    			(selected,:);
out.saccades.peakV.v	    		= in.saccades.peakV.v	    			(selected,:);
out.saccades.peakV.x	    		= in.saccades.peakV.x	    			(selected,:);
out.saccades.peakV.y	    		= in.saccades.peakV.y	    			(selected,:);
out.saccades.threshold.offset 		= in.saccades.threshold.offset 			(selected)	;
out.saccades.threshold.onset  		= in.saccades.threshold.onset						;

out.saccades.origin = in.saccades.origin(selected,:);
out.saccades.angle.overall = in.saccades.angle.overall(selected);


out.saccades.landmarkCoordinates.x =in.saccades.landmarkCoordinates.x(selected,:);
out.saccades.landmarkCoordinates.y =in.saccades.landmarkCoordinates.y(selected,:);


out.saccades.duration       		= in.saccades.duration       			(selected)	;
out.saccades.angle.onset  			= in.saccades.angle.onset  				(selected)	;
out.saccades.angle.offset  			= in.saccades.angle.offset  			(selected)	;
out.saccades.tail.duration 			= in.saccades.tail.duration 			(selected)	;
out.saccades.tail.amplitude  		= in.saccades.tail.amplitude  			(selected)	;
out.saccades.tail.SNR 				= in.saccades.tail.SNR 								;
out.saccades.problematic 			= in.saccades.problematic 				(selected,:);

out.saccades.fits.residual 			= in.saccades.fits.residual 			(selected)	;
out.saccades.fits.parameters 		= in.saccades.fits.parameters;
out.summary.interSaccadeInterval    = diff(in.saccades.landmarks.v(selected,1));
out.summary.SNR 			  		= in.summary.SNR;
out.summary.trialDuration	 		= in.summary.trialDuration;
out.summary.Fs				 		= in.summary.Fs;

if ~isempty(valint)
	out.summary.analyzedIntervals 		= valint;
else
	try
	out.summary.analyzedIntervals 		= in.summary.analyzedIntervals;
	catch
		disp('summary.analyzedIntervals field missing')
	end	
end


out.traces.position.raw 			= in.traces.position.raw;
out.traces.velocity.v 				= in.traces.velocity.v;
out.traces.velocity.xy 				= in.traces.velocity.xy;

out.summary.filteredSaccades 		= length(selected);
out.summary.unfilteredSaccades		= length(in.saccades.amplitude(:,1));
out.summary.yield = [' saccades selected: ' num2str(length(selected)) ' from ' num2str(length(in.saccades.amplitude(:,1)))];



out.saccades.peakV.description 		= in.saccades.peakV.description 		;
out.traces.remark 					= in.traces.remark 					;
out.saccades.threshold.description	= in.saccades.threshold.description;
out.saccades.landmarks.description  = in.saccades.landmarks.description ;

% out.summary.statistics = statisticsFromSaccadeStruct(out)



try
	in.saccades.vProfiles.v;
	out.saccades.vProfiles.v = in.saccades.vProfiles.v(selected,:);
	out.saccades.vProfiles.x = in.saccades.vProfiles.x(selected,:);
	out.saccades.vProfiles.y = in.saccades.vProfiles.y(selected,:);
catch
	warning('warning: no v profiles in saccade struct.')
end


    
