function out = saccadeDetection(varargin)
% ======= Saccadic Detection v0.93 ======
% 
% DETECTS SACCADES AND CALCULATES SACCADIC PROPERTIES
% 
% detects saccades that meet certain criteria:
% 1 have radial velocity peaks that cross an adaptive threshold 
% 2 have a sharp onset, a fast fall (main sequence) and a gradual falloff 
% 3 the velocity peak is the maximum peak within a window
% 
% Also, fits the saccades' main sequence (normalized) with gamma functions, to verify 
% whether they are well formed. Empirically, a fit is acceptable when it's residual is smaller than .3 degrees (cummulative over fitted saccade)
% 	* the output structure contains the fit results as well as fit parameters
% 
% 
%  ALGORITHM:
% 	 0. Determine adaptive threshold
% 	 1. Find Velocity Peaks
%    2. Find onset shoulder (left threshold crossing)
%	 3. Determine offset threshold as a function of eye velocity preceding onset 
% 	 4. Find main sequence end (maximum deceleration after velocity peak)
% 	 5. Find Standard end (threshold crossing) 
% 	 6. Find tail end, by following the humps in velocity until there are no bumps that cross offset threshold (calculated for every saccade)
% 	 7. Fit a gamma function to the main sequence, obtain the residual for the whole saccade
% 	 8. If the residual is above .35 mark the saccade as ill-formed
% 
%  REQUIRES:
%	 summaryFromSaccadeStruct.m
% 
%  INPUTS:
% 	eye position :
% 		- accepts either : 
% 			a cell array with two cells {X_position ; Y_position})
% 			a N x 2 matrix
%
%  PARAMETERS:
% 	parameters are passed as  parameter-value pairs (functionarguments,....,  'param1' , 'value1',  'param2' , 'value2')
% 
% 			'PARAMETER' , DEFAULT
%			'fixedthreshold', 20 (deg/s)
% 					fixed threshold to compare with standard literature
% 			'threshold', 4 
% 					threshold is a parameter for the adaptive threshold and for the calculation of the offset threshold
% 			'sgolayfilter',true
% 					whether to use a savitzky-golay filter (matlab help: help msgolay)
% 			'sgolaywindow', 10
% 			'maxsacfreq',10
% 					max saccadic frequency accepted - drop if more
% 			'detectionwindow',[90 190];
% 					window around the onset where to check for saccade parameters
% 			'method', 'veloc'
% 			'artifactamplitude', 800
% 			'Fs', 1000
% 					sampling rate
% 		**'debug', false
% 					-> show every single saccade and fit!
% 			'offsetthreshfactors', [.5 .5] %default .3 .7
% 			'output_traces', 1
% 			'removeartifacts', true
% 			'wellformed', true
% 			'plotresults', true
%				  plots all detected saccades, ordered by selected parameter
% 			'artifactmask', []
% 
%  OUTPUT:
% 
%   A structure with a number of self-explanatory >:D fields. Importantly:
% 
% 	a matrix where rows are saccades and columns are saccadic landmarks (interesting points)
% 	saccadic landmarks are given in trial_time in the following form:
% 		landmarks = [ONSET MAXVELOCITY NEXT_ACCELERATION_MINIMUM FIXED_THRESHOLD_CROSSING TAIL_END] 
%		size = 6 X no_saccades]
% 
% 
% 	DETECTION PROBLEMS:
% 	 	out.saccade.problems  - columns
% 
		% 1 
		% 2 
		% 3 
		% 4 
		% 5 
		
		% 7 
		% 8 
		% 9 
		
		% 10 pospk of velocity - prewin <= 0 
		% 11 problematic onset: acceleration inflexion (to positive) with velocity larger than threshold onset 
		% 12 large motion before onset, using only low velocity samples
		% 13 no clear end of the main sequence: time to maximum deceleration larger than rise time
		% 14 no clear end of the saccade: velocity always larger than offset (raised threshold)
% 
%
% 
% TODO: make interactive
% TODO: plot summary
% TODO: add the possibility to analyze selected intervals
% TODO: more docu
% 
% v.0.9.0 - May/2013  - Mario Negrello - mnegrello@gmail.com
% v.0.9.1 - May/2013  - Mario Negrello - mnegrello@gmail.com
% v.0.9.2 - Sept/2013  - Mario Negrello - mnegrello@gmail.com
% v.0.9.3 - Sept/2013  - Mario Negrello - mnegrello@gmail.com
% 
% -change computation of theta_offset when there are large fluctuations before eye motion
% -tag saccades where we had problems recognizing particular aspects in field out.saccades.problematic


p = inputParser;
p.addRequired('EyeChannels')  % a matrix with two columns or a cell array with two cells;
p.addParamValue('threshold', 4) % stdandard deviation criterion for offset threshold
p.addParamValue('fixedthreshold', 20) 
p.addParamValue('sgolayfilter',true)
p.addParamValue('sgolaywindow', 10)
p.addParamValue('sgolayall', false) % for very noisy data, this may enhance recognition (at the expense of some temporal imprecisions)
p.addParamValue('maxsacfreq',15)
p.addParamValue('detectionwindow',[100 150]);
p.addParamValue('artifactamplitude', 1300)
p.addParamValue('Fs', 1000)
p.addParamValue('debug', false)
p.addParamValue('offsetthreshfactors', [.3 .5]) %default .3 * mean .7 *std
p.addParamValue('output_traces', 1)
p.addParamValue('removeartifacts', true)
p.addParamValue('wellformed', false) %only return well formed saccades?
p.addParamValue('plotresults', true)
p.addParamValue('artifactmask', [])
p.addParamValue('verbose',false)
p.addParamValue('accept_problematic',true)
p.addParamValue('onsetthresholdtype','fancy');
p.addParamValue('detectionintervals',[]);
p.addParamValue('colorscheme','light');

adaptivethreshold = 1;

%%========================================================     Method Parameters
p.parse(varargin{:});

	EyeChannels = p.Results.EyeChannels;
	threshold   = p.Results.threshold;
	fixedthresh = p.Results.fixedthreshold;
	sgolaywindow = p.Results.sgolaywindow;
	sgolayfilter = p.Results.sgolayfilter;
	offsetthreshfactors = p.Results.offsetthreshfactors;
	detectionwindow = p.Results.detectionwindow;
	plotresults = p.Results.plotresults;
	artifactamplitude = p.Results.artifactamplitude;
	sgolayall = p.Results.sgolayall;
	wellformed = p.Results.wellformed;
	displaywarnings = p.Results.verbose;
	accept_problematic = p.Results.accept_problematic;
	onsetthresholdtype = p.Results.onsetthresholdtype;
	detectionintervals = p.Results.detectionintervals;
	colorscheme = p.Results.colorscheme;

	% method = p.Results.method;
	Fs = p.Results.Fs;

		testopts.max_acceleration_accepted  = 40000;				
		testopts.max_velocity_accepted 	    = 1100;				
		testopts.min_velocity_accepted      = 90;				
		testopts.accept_problematic = accept_problematic;				
		
		testopts.detectionwindow = detectionwindow;				

	removeartifacts = p.Results.removeartifacts;

	debugmode = p.Results.debug;
	output_traces = p.Results.output_traces;

	maxsacfreq = p.Results.maxsacfreq;

	fitemphasis = 'mainseq';
	% fitemphasis = 'tail';

if debugmode
	plot_problems = true;
	displaywarnings = true;
else
	plot_problems = false;
	displaywarnings = false;
end






%%============================================================     preliminaries
	if ~displaywarnings
		warning off
	else
		warning on
	end

	sgolaywindow = sgolaywindow*Fs/1000;
	winspanl = detectionwindow(1)*Fs/1000;
	winspanr = detectionwindow(2)*Fs/1000;

	if iscell(EyeChannels)
		% if not row, make it so
			if size(EyeChannels{1},2) == 1
					EyeChannels{1} = EyeChannels{1}; EyeChannels{2} = EyeChannels{2}; 
					xy(:,1) = EyeChannels{1};
					xy(:,2) = EyeChannels{2};
					else
					EyeChannels{1} = EyeChannels{1}'; EyeChannels{2} = EyeChannels{2}';
					xy(:,1) = EyeChannels{1};
					xy(:,2) = EyeChannels{2};
			end

		else
			if size(EyeChannels,2)==2
					xy = EyeChannels;
				elseif size(EyeChannels,1)==2
					xy = EyeChannels';
				else
					disp('Critical: eye data format not recognized. check input formats.')
					return
			end
	end

%%============================================     Savitzky-Golay Filters (position, velocity and acceleration)

	LoT = length(xy)/Fs;

	if sgolayfilter & not(sgolaywindow==0)

			s = [xy];
			t = [1:length(s)]';
			nans = isnan(prod(s')');
			
			% fix nans to speed msgolay - we reverse that afterwards
			s(nans,1) = interp1(t(~nans) ,s(~nans,1),t(nans));
			s(nans,2) = interp1(t(~nans) ,s(~nans,2),t(nans));

			xy = mssgolay(t,s,'span', sgolaywindow);
			
			dxy = [0 0 ; diff(xy) ];	
			if sgolayall
				dxy = mssgolay(t,dxy,'span', sgolaywindow);
			end

			linV =  sqrt(sum(dxy.^2,2));	

			ddxy = [diff(dxy) ; 0 0 ];
			linA = [diff(linV) ; 0];


			if sgolayall
				linA = mssgolay(t,linA,'span', sgolaywindow);
				ddxy = mssgolay(t,ddxy,'span', sgolaywindow);
			end

			% to compare with original signal:
			% plot([VxY [0 0 ; diff(xy)] ])

			noise = mean(std(xy - s));
		
			else

			s = [xy];
			t = [1:length(s)]';
			
			dxy = [0 0 ; diff(xy) ];		
			ddxy = [ diff(dxy) ; 0 0 ];

			linV =  [0; sqrt(sum(dxy.^2,2))];
			linA = [diff(linV) ; 0];	
			
			noise = mean(std(xy-s));
			W = [];	

			s(nans,:) = NaN;

	end	

	linV = linV*Fs;
	dxy = dxy*Fs; % in deg/s
	% linA = linA*Fs;


%%=======================================================     ADAPTIVE THRESHOLD
	% for all samples with velocities lower than AdapThr calculate the meand and std. Then update threshold to Median+6*sigma. 
	% -- converges to low enough to detect as many saccades as the noise level allows while reducing false positives.
	% See Bernquist and Nystrom 2010.

	if adaptivethreshold
		theta(1) = adaptiveThreshold(linV,Fs,threshold);
		theta(2) = adaptiveThreshold(abs(dxy(:,1)),Fs,threshold);
		theta(3) = adaptiveThreshold(abs(dxy(:,2)),Fs,threshold);
		else 
		theta(1,2,3) = fixedthreshold;
	end


% 
% ___   __                   _______           __   ______                ___     __      __           
% _/_(_)_/_/_________       / ____(_)___  ____/ /  / ____/___ _____  ____/ (_)___/ /___ _/ /____  _____
% _/  _/_//____/____/      / /_  / / __ \/ __  /  / /   / __ `/ __ \/ __  / / __  / __ `/ __/ _ \/ ___/
% _ _/_/_/____/____/      / __/ / / / / / /_/ /  / /___/ /_/ / / / / /_/ / / /_/ / /_/ / /_/  __(__  ) 
% _)_/ (_)               /_/   /_/_/ /_/\__,_/   \____/\__,_/_/ /_/\__,_/_/\__,_/\__,_/\__/\___/____/  
% 

	% mpd = round((1/maxsacfreq)*Fs); 
	% mph = theta(1);
	mpd = 150;
	mph = 100;
	[val_ candidates] = findpeaks(linV, 'minpeakheight', mph, 'minpeakdistance', mpd);


% 	% find valid trials
% 	accepted = [];
% 	for s = candidates'
% 	
% 		for i = 1:size(detectionintervals,1)
% 			% pause(.1)
% 	        if s >= detectionintervals(i,1) & s <= detectionintervals(i,2)
% 	            accepted = [accepted s];
% 
% 	            continue
% 	        end
% 	    end
% 	end
% 	
% 
% 	candidates = unique(accepted);
	val_ = linV(candidates);





	%%========================================================     artifact handling

		if removeartifacts
			peak_vel_accepted = 1100;
			candidates = candidates(find(val_<peak_vel_accepted));
			val_ = val_(find(val_<peak_vel_accepted));
		end


		disp(['number of tentative peaks:' num2str(length(candidates))])
		disp(['average of tentative peaks:' num2str(mean(val_))])
		
		meanvalofpeaks = mean(val_);
		if meanvalofpeaks>600
			warning('Either Monkey is The Flash or we have conversion problems: mean(peakvelocity)>1000. Check conversion to degrees.')
			
		end

		% prune

		todiscardL = find(candidates-winspanl<=0)';
		todiscardR = find(candidates+winspanr > length(xy))';

		
		todiscard = [todiscardL todiscardR];
		warning(['discarding ' num2str(length(todiscard)) ' candidate(s) (out of bounds).'])
		candidates(todiscard) = [];

	%%=======================================================     Candidate Saccade Boundaries
	indxSelect = @(n) round([candidates(n)-winspanl : candidates(n)+winspanr]);



%% =====================================================================   
	%     __  ___      _          __                    
	%    /  |/  /___ _(_)___     / /   ____  ____  ____ 
	%   / /|_/ / __ `/ / __ \   / /   / __ \/ __ \/ __ \
	%  / /  / / /_/ / / / / /  / /___/ /_/ / /_/ / /_/ /
	% /_/  /_/\__,_/_/_/ /_/  /_____/\____/\____/ .___/ 
	%                                          /_/      
%% =====================================================================   



c = 0;SO=1; lastcandidate = candidates(1);
while SO<=length(candidates)  %%	iterate through candidates
	if ~plot_problems
			clc
			progress = ['progress: ' num2str(SO) ' of ' num2str(length(candidates))];
			disp(progress);
		else
			disp(['[===============================================]'])
			disp(['[===== analyzing:' num2str(SO)])
			disp(['[===== candiate:' num2str(candidates(SO))])
	end


		current_candidate = candidates(SO);
		V_i = linV(indxSelect(SO));
		A_i = linA(indxSelect(SO));
		Vx_i = dxy(indxSelect(SO),1);
		Vy_i = dxy(indxSelect(SO),2);
		XY = xy(indxSelect(SO),:);

		if length(find(V_i==NaN))>1
			SO = SO+1;
			continue
		end

		if max(V_i)<testopts.min_velocity_accepted
			SO = SO+1;
			continue
		end

		if max(V_i)>testopts.max_velocity_accepted;
			SO = SO+1;
			continue
		end



		
		%%%==================================================     find saccadic landmarks
		
		params.fs = Fs;
		params.theta_off_high = fixedthresh;
		params.theta_on = theta(1);
		params.detectionwindow = [winspanl winspanr];
		params.stdfactor = 2;
		params.onsetthresholdtype = onsetthresholdtype;


	

		try
			[landmarks thetaoffv problematic_v] = findLandmarks( V_i,  A_i,  params);
		catch
			disp('problem with findLandmarks')
			SO = SO+1;
			keyboard
			continue
		end
		
		warning(num2str(landmarks))
		disp(['theta onset:' num2str(theta(1))])

		try
			% find maximum of x and y components
			[max_v_x max_v_x_i] = max(abs(Vx_i (landmarks(1):landmarks(4)))) ; 
			max_v_x_i = max_v_x_i + landmarks(1);
			[max_v_y max_v_y_i] = max(abs(Vy_i (landmarks(1):landmarks(4))));
			max_v_y_i = max_v_y_i + landmarks(1);

		catch
			disp('some problem when finding max of V and H components.')
			keyboard
		end

		
		%    __            __      
		%   / /____  _____/ /______
		%  / __/ _ \/ ___/ __/ ___/
		% / /_/  __(__  ) /_(__  ) 
		% \__/\___/____/\__/____/  
		
		testopts.theta_on = theta(1);
		testopts.fs = Fs;
		[problem_tags] = testSaccade(V_i, A_i, landmarks, testopts);
			% 1 riseonset = 0
			% 2 minimum velocity not reached
			% 3 velocity larger than maximum accepted
			% 4 falloff very slow
			% 5 mean preceding velocity too large
			% 6 large velocities after offset
			% 7 too many peaks within detection window
			% 8 doubled landmarks
			% 9 saccade with negative duration?

				problems = [problem_tags problematic_v]
				problemcount = length(find(problem_tags));

			if ~problemcount
				% if passed all tests, count.
				c = c+1; % count saccades that passed all tests
				SO = SO+1;
				elseif problemcount & ~accept_problematic
				SO = SO+1;
				elseif problemcount & accept_problematic 
				c = c+1; % count saccades that passed all tests
				SO = SO+1;
			end

			problematic_saccades(c,:) = problems;

	        if problems(1)==1
	        	if plot_problems
	        		figure(1000), plot(V_i)
	        	end
	        	% issue  = 'argh'
	            landmarks = repairSaccade(V_i, A_i, landmarks, issue);
	        end


	        % fall faster than rise
	        if problems(3)==1
	        	if plot_problems
	        		figure(1000), plot(V_i)
	        	end
	        	issue  = 'fall faster than rise';
	        	if displaywarnings
	        		warning(issue)
	        	end
	        	
	            % landmarks = repairSaccade(V_i, A_i, landmarks, issue);
	        end


			if problems(5)==1

				issue = 'prec_vel_high';
				if displaywarnings
					warning(issue)
				end
	            landmarks = repairSaccade(V_i, A_i, landmarks, issue);
	        end


			if problems(6)==1

				issue = 'nans in samples';
				if displaywarnings
					warning(issue)
				end
	            % landmarks = repairSaccade(V_i, A_i, landmarks, issue);
	        end




	   %     _____________________ __
 	   %    / ____/  _/_  __/ ___// /
 	   %   / /_   / /  / /  \__ \/ / 
 	   %  / __/ _/ /  / /  ___/ /_/  
 	   % /_/   /___/ /_/  /____(_)   
         	                               
        % opts.fitemphasis = 'mainseq';
        opts.fitemphasis = 'peak-thresh';

        if ~problems([6 10 12 11]) & landmarks(3)-landmarks(2)>2

        		fits = fit2gamma(V_i, landmarks, opts);

		        residual(c) = fits.residual;
		        pset 	    = fits.pset;
		        fitspecs    = fits.fitspecs;
		        predicted_sac = fits.predicted_fit;
				observed_sac  = fits.observed_sac; 


        elseif ~isempty(find(V_i)==NaN)
           
				disp('cannot fit gamma - NaN in velocity profile.')
				residual(c) = NaN;
		        pset 	    = NaN;
		        fitspecs    = NaN;
		        predicted_sac = NaN;
				observed_sac  = NaN; 

			else

				fits = fit2gamma(V_i, [1 80 101 150 0 0] , opts)

				residual(c) = fits.residual;
		        pset 	    = fits.pset;
		        fitspecs    = fits.fitspecs;
		        predicted_sac = fits.predicted_fit;
				observed_sac  = fits.observed_sac; 
			
		end

		


			

		%%===================================================     Saccade properties

			landmarks_v(c,:) = landmarks + current_candidate - winspanl; landmarks_v(c,:) = landmarks_v(c,:) .* (landmarks>0);
			warning(num2str(landmarks_v(c,:)))

			peakVel_v(c) = V_i(landmarks(2));
			peakVel_x(c,:) = [max_v_x max_v_x_i];
			peakVel_y(c,:) = [max_v_y max_v_y_i];

			peakVelocity(c,1) = V_i(landmarks(2));

			vProfiles_v(c,:) = V_i;
			vProfiles_x(c,:) = Vx_i;
			vProfiles_y(c,:) = Vy_i;

			% sacXY{c} = xy(indxSelect(SO),:);

			%saccadeLandmarks expressed in time units
			saccadeLandmarks_v(c,:) = double(landmarks_v(c,:) )/(Fs)*1000;
			
			SL(c,:) = saccadeLandmarks_v(c,:);
			
			SL_ind(c,:) = landmarks_v(c,:);
			
			saccadeAmplitude(c,1) = sqrt((xy(landmarks_v(c,4),1) - xy(landmarks_v(c,1),1)).^2	+ (xy(landmarks_v(c,4),2) - xy(landmarks_v(c,1),2)).^2);
			saccadeAmplitude(c,2) = xy(landmarks_v(c,4),1) - xy(landmarks_v(c,1),1);
			saccadeAmplitude(c,3) = xy(landmarks_v(c,4),2) - xy(landmarks_v(c,1),2);

			saccadeOnsetXY (c,:)= [xy(landmarks_v(c,1),1)  xy(landmarks_v(c,1),2)];

			if landmarks_v(c,6)~=0
				landmarkCoordinates.x(c,:) = [ xy(landmarks_v(c,1),1) xy(landmarks_v(c,2),1) xy(landmarks_v(c,3),1) xy(landmarks_v(c,4),1) xy(landmarks_v(c,5),1) xy(landmarks_v(c,6),1) ];
				landmarkCoordinates.y(c,:) = [ xy(landmarks_v(c,1),2) xy(landmarks_v(c,2),2) xy(landmarks_v(c,3),2) xy(landmarks_v(c,4),2) xy(landmarks_v(c,5),2) xy(landmarks_v(c,6),2) ];
			else
				landmarkCoordinates.x(c,:) = [ xy(landmarks_v(c,1),1) xy(landmarks_v(c,2),1) xy(landmarks_v(c,3),1) xy(landmarks_v(c,4),1) 0 0 ];
				landmarkCoordinates.y(c,:) = [ xy(landmarks_v(c,1),2) xy(landmarks_v(c,2),2) xy(landmarks_v(c,3),2) xy(landmarks_v(c,4),2) 0 0 ];
			end
			
			
			if landmarks_v(c,5)~=0 & landmarks_v(c,6)~=0
				tailAmplitude(c) = sqrt((xy(landmarks_v(c,6),1) - xy(landmarks_v(c,4),1))^2 +(xy(landmarks_v(c,6),2) - xy(landmarks_v(c,4),2)).^2);
				tailDuration(c)  = saccadeLandmarks_v(c,6) - saccadeLandmarks_v(c,4);
			else
				tailAmplitude(c) = 0;
				tailDuration(c)  = 0;
			end



			angleAtOnset(c,:)  =  atan2( xy(landmarks_v(c,2),2) - xy(landmarks_v(c,1),2),  ...
									xy(landmarks_v(c,2),1) - xy(landmarks_v(c,1),1));

			angleAtOffset(c,:)  =  atan2( xy(landmarks_v(c,4),2) - xy(landmarks_v(c,2),2),  ...
			       				  	 xy(landmarks_v(c,4),1) - xy(landmarks_v(c,2),1));

			
			angleOnToOff(c,:) =   atan2( xy(landmarks_v(c,4),2) - xy(landmarks_v(c,1),2),  ...
			       				  	 xy(landmarks_v(c,4),1) - xy(landmarks_v(c,1),1));



			sacdur_v(c) = (landmarks_v(c,4) - landmarks_v(c,1))*Fs/1000;
			

			thetaoffset(c) = thetaoffv;
	
	
	%  DEBUGGING PLOTS
	if debugmode & ~sum(SL(c,:))==0;
		w = 100;
		blue =  [.6 .6 1];
		red  = [1 .3 .2];
		redder = [1 0 0];
		green = [.6 1 .8];
		yellow = 'y';
		cyan   = 'c';
		indc = (current_candidate-winspanl:current_candidate+winspanr);
		figure(100), clf
		set(0,'defaultaxescolor', [.4 .4 .4]);

		pv = peakVelocity(c,1);		
		xa = indc/Fs*1000;
		% xa = (winspanl:winspanr)/Fs;
		subplot(2,4,1:2); 
			line([SL(c,1) SL(c,1)], [-pv  pv],'color', green)
			line([SL(c,2) SL(c,2)], [-pv  pv],'color', yellow);
			line([SL(c,3) SL(c,3)], [-pv  pv],'color', red);
			line([SL(c,4) SL(c,4)], [-pv  pv],'color', redder);
			if SL(c,5)~=0
			line([SL(c,5) SL(c,5)], [-pv  pv],'color', blue);
			line([SL(c,6) SL(c,6)], [-pv  pv],'color', redder);
			end


			line(xa, linV(indc),'color', green, 'linewidth', 3);
			line(xa, dxy(indc,1), 'linewidth', 1,'color',red);
			line(xa, dxy(indc,2), 'linewidth', 1,'color',blue);
			line(xa, log(linV(indc)*1e3)/max(log(linV(indc)*1e3))*max(linV(indc)),'color', yellow)
			line(xa, linA(indc)/max(linA(indc))*max(linV(indc))+max(linV(indc)),'color', redder,'linewidth', 1); 
			line([xa(1) xa(end)]',[theta(1) theta(1)]','color', [.8 .8 .8]);
			title('velocity (deg/samples ) x Fs')
			axis tight


		subplot(2,4,3);
			exn = xy(indc,1);
			eyn = xy(indc,2);
			line([SL(c,1) SL(c,1)], [xy(SL_ind(c,1),1)-exn(1) xy(SL_ind(c,1),2)-eyn(1)],'color', green ,'markersize', 16,'marker', 'o')
			line([SL(c,2) SL(c,2)], [xy(SL_ind(c,2),1)-exn(1) xy(SL_ind(c,2),2)-eyn(1)],'color', yellow ,'markersize', 16,'marker', 'o')
			line([SL(c,3) SL(c,3)], [xy(SL_ind(c,3),1)-exn(1) xy(SL_ind(c,3),2)-eyn(1)],'color', red ,'markersize', 16,'marker', 'o')
			line([SL(c,4) SL(c,4)], [xy(SL_ind(c,4),1)-exn(1) xy(SL_ind(c,4),2)-eyn(1)],'color', redder ,'markersize', 16,'marker', 'o')

			line([SL(c,1)], xy(SL_ind(c,1),2)-eyn(1),'color', green	   ,'markersize', 6 ,'marker', 'o')
			line([SL(c,2)], xy(SL_ind(c,2),2)-eyn(1),'color', yellow ,'markersize', 6 ,'marker', 'o')
			line([SL(c,3)], xy(SL_ind(c,3),2)-eyn(1),'color', red  ,'markersize', 6 ,'marker', 'o')
			line([SL(c,4)], xy(SL_ind(c,4),2)-eyn(1),'color', redder  ,'markersize', 6 ,'marker', 'o')
			
			line(xa, exn-exn(1), 'linewidth', 2,'color',red); 
			line(xa, eyn-eyn(1), 'linewidth', 2,'color',blue); 
			
			line([SL(c,4) SL(c,4)], [xy(SL_ind(c,4),1)-exn(1) xy(SL_ind(c,4),2)-eyn(1)],'color', red  ,'markersize', 16,'marker', 'o')
			line([SL(c,4)], 		xy(SL_ind(c,4),2)-eyn(1),							'color', red  ,'markersize', 6 ,'marker', 'o')
			if SL(c,5)~=0 & SL(c,6)~=0 
				line([SL(c,5) SL(c,5)], [xy(SL_ind(c,5),1)-exn(1) xy(SL_ind(c,5),2)-eyn(1)],'color', cyan  ,'markersize', 16,'marker', 'o')
				line([SL(c,5)], 		xy(SL_ind(c,5),2)-eyn(1),							'color', redder  ,'markersize', 6 ,'marker', 'o')
			end
			title('position (deg)')
			axis tight



		subplot(2,4,4); 
			hold on
			plot(xy(indc,1), xy(indc,2))
			plot(xy(SL_ind(c,1),1), xy(SL_ind(c,1),2), 'o','markersize', 6,'color', green)
			plot(xy(SL_ind(c,2),1), xy(SL_ind(c,2),2), 'o','markersize', 6,'color', yellow)
			plot(xy(SL_ind(c,3),1), xy(SL_ind(c,3),2), 'o','markersize', 6,'color', red)
			
			if SL(c,4)~=0
			plot(xy(SL_ind(c,4),1), xy(SL_ind(c,4),2), 'o','markersize', 6,'color',redder)
			end

			if SL(c,5)~=0 & SL(c,6)~=0 
				plot(xy(SL_ind(c,5),1), xy(SL_ind(c,5),2), 'o','markersize', 6,'color', blue)
				plot(xy(SL_ind(c,6),1), xy(SL_ind(c,6),2), 'o','markersize', 6,'color', redder)
			end

			% axis([min(xy(:,1)) max(xy(:,1)) min(xy(:,2)) max(xy(:,2))])
			axis([min(xy(:)) max(xy(:)) min(xy(:)) max(xy(:))])
			axis([-20 20 -20 20 ])
			hold off
			title('trace')
			xlabel('deg')

		subplot(2,4,5)
			blax = landmarks(1):landmarks(2)+100;

			
				hold on
				% plot(V_i(landmarks(1):end)/max(V_i(landmarks(1):end)),'r')
				plot(observed_sac,'b')
				plot(predicted_sac,'g');
				legend({'observed' ; 'predicted'},'color',[1 1 1])
				hold off



			title(['residual:' num2str(residual(c)) '.'])
		

		subplot(2,4,6:8)
			if ~(current_candidate-winspanl*15<1) & ~(current_candidate+winspanr*15>length(linV))
				
				line([SL(c,1) SL(c,1)], [0 100], 'color', 'r')
				line([SL(c,2) SL(c,2)], [0 100], 'color', 'r')
				line([SL(c,3) SL(c,3)], [0 100], 'color', 'r')
				line([current_candidate-winspanl*15 current_candidate+winspanr*15]/Fs*1000,[theta(1) theta(1)]','color', [.8 .8 .8]);
				line([current_candidate-winspanl*15:current_candidate+winspanr*15]/Fs*1000, linV(current_candidate-winspanl*15:current_candidate+winspanr*15))
				
				ylim([0 1200])
				title(['event number: ' num2str(c)])
			end

			automated = 0;
			if ~automated
				text(0,0,'click anywhere to go on')
				waitforbuttonpress
			end

	end
	
end
	
%     __  ______    _____   __   __    ____  ____  ____  ____     _______   ______ 
%    /  |/  /   |  /  _/ | / /  / /   / __ \/ __ \/ __ \/ __ \   / ____/ | / / __ \
%   / /|_/ / /| |  / //  |/ /  / /   / / / / / / / / / / /_/ /  / __/ /  |/ / / / /
%  / /  / / ___ |_/ // /|  /  / /___/ /_/ / /_/ / /_/ / ____/  / /___/ /|  / /_/ / 
% /_/  /_/_/  |_/___/_/ |_/  /_____/\____/\____/\____/_/      /_____/_/ |_/_____/  
                                                                                 

	if ~exist('saccadeLandmarks_v')
		disp('no saccades found! ')
		return
	end

%          ___           
%   ____ _/ (_)___ _____ 
%  / __ `/ / / __ `/ __ \
% / /_/ / / / /_/ / / / /
% \__,_/_/_/\__, /_/ /_/ 
%          /____/        

	alignment = 'onset';
	switch alignment
		case 'onset'
				% peak index - onset index, in trial time
				d = round(saccadeLandmarks_v(:,2) - saccadeLandmarks_v(:,1));
				saccadeLandmarks_v = saccadeLandmarks_v - repmat(d,1,6);% + repmat(saccadeLandmarks_v(:,1),1,6);
				
				saccadeLandmarks_v(find(saccadeLandmarks_v<0)) = 0;
				
				shifts_row_v = @(row) circshift(vProfiles_v(row,:)', d(row));
				shifts_row_x = @(row) circshift(vProfiles_x(row,:)', d(row));
				shifts_row_y = @(row) circshift(vProfiles_y(row,:)', d(row));

				vProfiles_v = cell2mat(arrayfun(shifts_row_v, [1:length(d)],'uniformoutput', false))';
				vProfiles_x = cell2mat(arrayfun(shifts_row_x, [1:length(d)],'uniformoutput', false))';
				vProfiles_y = cell2mat(arrayfun(shifts_row_y, [1:length(d)],'uniformoutput', false))';
	end


%    ____  __  ____________  __  ________   _______________  __  ________________  ______  ______
%   / __ \/ / / /_  __/ __ \/ / / /_  __/  / ___/_  __/ __ \/ / / / ____/_  __/ / / / __ \/ ____/
%  / / / / / / / / / / /_/ / / / / / /     \__ \ / / / /_/ / / / / /     / / / / / / /_/ / __/   
% / /_/ / /_/ / / / / ____/ /_/ / / /     ___/ // / / _, _/ /_/ / /___  / / / /_/ / _, _/ /___   
% \____/\____/ /_/ /_/    \____/ /_/     /____//_/ /_/ |_|\____/\____/ /_/  \____/_/ |_/_____/   
% \____/\____/ /_/ /_/    \____/ /_/     /____//_/ /_/ |_|\____/\____/ /_/  \____/_/ |_/_____/





	out.saccades.amplitude	    = saccadeAmplitude;  % in a straightline from onset to offset - has three compoents: total distance covered, x and y.

	out.saccades.landmarks.description = {'onset' 'maxV' 'maX_neg_accel_post_pk' 'highthresh' 'tailpk' 'tailend'};
	out.saccades.landmarks.v = saccadeLandmarks_v;
	
	out.saccades.origin = saccadeOnsetXY;
	out.saccades.landmarkCoordinates.x = landmarkCoordinates.x;
	out.saccades.landmarkCoordinates.y = landmarkCoordinates.y;

	out.saccades.peakV.v	    = peakVelocity; % in a straightline from onset to offset
	out.saccades.peakV.x	    = peakVel_x; 
	out.saccades.peakV.y	    = peakVel_y; 
	out.saccades.peakV.description ={'max value' 'index'}

	out.saccades.threshold.offset = thetaoffset;
	out.saccades.threshold.onset  = theta(1);
	out.saccades.threshold.description = 'adaptive';

	out.saccades.angle.onset   = angleAtOnset;
	out.saccades.angle.offset  = angleAtOffset
	out.saccades.angle.overall = angleOnToOff;

	out.saccades.tail.duration = tailDuration;
	out.saccades.tail.amplitude  = tailAmplitude;
	out.saccades.tail.SNR = mean(tailAmplitude)/noise;

	out.saccades.problematic = problematic_saccades;

	out.saccades.duration       = sacdur_v;
	
	out.saccades.fits.residual = residual;
	out.saccades.fits.parameters = pset;
	out.saccades.fits.fitspecs = fitspecs;

	out.summary.analyzedIntervals = detectionintervals;
	out.summary.interSaccadeInterval    = diff(saccadeLandmarks_v(:,1));
	out.summary.SNR 			  		= mean(saccadeAmplitude)/noise; %TODO: correct s
	out.summary.trialDuration	 = length(xy)/Fs; % in seconds
	out.summary.Fs				 = Fs;
	out.summary.numberOfSaccades = c;

	out.traces.position.raw = xy;
	out.traces.velocity.v = linV;
	out.traces.velocity.xy = dxy;
	out.traces.remark = 'sgolay filtered';


	if output_traces
		out.saccades.vProfiles.v = vProfiles_v;
		out.saccades.vProfiles.x = vProfiles_x;
		out.saccades.vProfiles.y = vProfiles_y;
	end


%                                   _                                __         __      __  __  _            
%  _      ___________ _____  ____  (_)___  ____ _   ____ _____  ____/ /  ____  / /___  / /_/ /_(_)___  ____ _
% | | /| / / ___/ __ `/ __ \/ __ \/ / __ \/ __ `/  / __ `/ __ \/ __  /  / __ \/ / __ \/ __/ __/ / __ \/ __ `/
% | |/ |/ / /  / /_/ / /_/ / /_/ / / / / / /_/ /  / /_/ / / / / /_/ /  / /_/ / / /_/ / /_/ /_/ / / / / /_/ / 
% |__/|__/_/   \__,_/ .___/ .___/_/_/ /_/\__, /   \__,_/_/ /_/\__,_/  / .___/_/\____/\__/\__/_/_/ /_/\__, /  
%                  /_/   /_/            /____/                       /_/                            /____/   

	if wellformed
		out = filterSaccadeStruct(out);
	end

	if plotresults
		out = summaryFromSaccadeStruct(out,'userselection','all');
	end	


%%=======================================================     adaptive threshold

function [landmarks  theta_off problematic] = findLandmarks(dx,ddx,params)

	% maximum number of peaks
	% m_i = minimum(V)
	% go to first maximum
	% is it within window?
	% - find next minimum
	num_peaks = 0;
	theta_off = 0;
	problematic = false(1,5);
	stdfactor = params.stdfactor;

	theta_on = params.theta_on;
	theta_off_high = params.theta_off_high;
	theta_on_fixed = 30;
	fs = params.fs;
	detectionwindow = params.detectionwindow;
	onsetthresholdtype = params.onsetthresholdtype;

	maxwin = round(20*fs/1000); % for bump detection
	minwin = round(25*fs/1000); % for bump detection
	prewin = round(40*fs/1000); 
	pospkwin = round(45*fs/1000);

	% init
	landmarks = [0 0 0 0 0 0];

	% make v always positive
	dx = abs(dx);
	ddx = diff(dx);

	pv = dx(detectionwindow(1)+1);
	pospk_v = detectionwindow(1)+1;
	
	[pa pospk_a] = max(ddx(1:pospk_v-1));


	% LANDMARK 1: ONSET
	
	switch onsetthresholdtype
		case 'acceleration'
    		riseon =  find( ddx(pospk_v-prewin:pospk_v) > theta_on , 1,'first');	
    	case 'velocity'
    		riseon =  find( dx(pospk_v-prewin:pospk_v) > theta_on , 1,'first');	
    	case 'fancy'
    		riseon =  find( medfilt1(dx(pospk_v-prewin:pospk_v),7) > theta_on , 1,'first');	

    		%number of samples between onset and maxv
    		


    end


	if isempty(riseon)
			warning('onset is empty: problematic')
			problematic(1) = 1;

			riseon = find(dx(pospk_v-prewin:pospk_v) > 3*median(dx(pospk_v-prewin:pospk_v)), 1, 'first');

			landmarks(1) = riseon + pospk_v - prewin - 2;
		else
			landmarks(1) = riseon + pospk_v - prewin - 2;
	end

	% LANDMARK 2: PEAK V
	landmarks(2) = pospk_v;

	upstroke_duration = landmarks(2) - landmarks(1);


	% COMPUTE OFFSET THRESHOLD
	% are there large eye motion before the rise onset? 
	
	large_fluctuations_before_riseon = find(dx(landmarks(1)-prewin:landmarks(1)) > theta_on);

	if length(large_fluctuations_before_riseon) > 20*fs/1e3
		warning('large eye motion before onset: using fixed threshold');
		theta_off = theta_on_fixed;
		% theta_off = median(low_vel_samples) + stdfactor*std(low_vel_samples);
		problematic(2) = 1;
	else
		theta_off = median(dx(1:riseon-prewin)) + stdfactor*std(dx(1:riseon-prewin));
	end
	

	% LANDMARK 3: MAIN SEQUENCE END

	% MAIN SEQUENCE END - MAXIMUM DECELERATION
	mainseq_offset = find(ddx(landmarks(2):landmarks(2)+upstroke_duration) == min(ddx(landmarks(2):landmarks(2)+upstroke_duration)) ,1,'first');
	
	if mainseq_offset > 1.5 * upstroke_duration
		warning('mainsequence offset > 1.5 * upstroke duration: problematic')
		problematic(3) =1;
		mainseq_offset = ceil(.5 * upstroke_duration);
		landmarks(3) = landmarks(2)+mainseq_offset;
	elseif isempty(mainseq_offset)
		problematic(3) =2;
		mainseq_offset = ceil(.5 * upstroke_duration);
		landmarks(3) = landmarks(2)+mainseq_offset;
	else
		landmarks(3) = landmarks(2)+mainseq_offset;
	end


	% FIND OFFSET THRESHOLD CROSSING
	mn = find(dx(landmarks(3):end) < theta_off_high,1, 'first');

	if isempty(mn)
			warning('could not find saccade end, estimating falloff: problematic')
			
			landmarks(4) = landmarks(3) + ceil(1.5*upstroke_duration); 
			if landmarks(4)> length(dx)
				landmarks(4) == length(dx);
			end

			problematic(4) = 1;
		else
			landmarks(4) = mn + landmarks(3) - 1 ;
	end





	n = 4;
	ref = landmarks(n);
	sacw = dx(ref:end);
	lasth = false;
	mxv = theta_off;
	
	% call to find humps
	while ~lasth
		n = n +1;
	
		if maxwin >= length(sacw)
			lasth = true;
			break
		end
		
		mx = find_next_max_above_threshold(sacw, maxwin, mxv);

		if mx==0
			lasth = true;
			break
		else
			mxv = sacw(mx);
			sacw = sacw(mx+1:end);
			landmarks(n) = mx+ref-1;
			ref = landmarks(n);
			n = n +1;
		end
	
		if minwin >=length(sacw)
			lasth = true;
			break
		end

		mn = find_next_min_below_threshold(sacw,minwin,theta_off);
	
		if isempty(mn)
			lasth = true;
			landmarks(n) = length(sacw);
			break
		else
			ref = ref+mn;
			landmarks(n) = ref;
		end
	
		sacw = dx(ref+1:end);

	end

	
	if length(find(landmarks)) > 6
		landmarks = [landmarks(1:5) landmarks(end)];
	end




	%     _____           __   __                              
	%    / __(_)___  ____/ /  / /_  __  ______ ___  ____  _____
	%   / /_/ / __ \/ __  /  / __ \/ / / / __ `__ \/ __ \/ ___/
	%  / __/ / / / / /_/ /  / / / / /_/ / / / / / / /_/ (__  ) 
	% /_/ /_/_/ /_/\__,_/  /_/ /_/\__,_/_/ /_/ /_/ .___/____/  
	%                                           /_/            

function [theta] = adaptiveThreshold(V,Fs,param1)

	% cocr = 1/Fs; % convergence criterion 1deg/s
	% T0 = 50/Fs; % initial threshold 50deg/s

	cocr = 1; % convergence criterion 1deg/s
	T0 = 50; % initial threshold 50deg/s



	stdveye = std(V);	


	fixatsamples = find(V<T0);
	Told = T0;
	Tnew = median(V(fixatsamples)) + param1 * std(V(fixatsamples));

	while abs(Tnew-Told)>cocr 
		Told = Tnew;
		fixatsamples = find(V<Tnew);
		Tnew = median(V(fixatsamples)) + param1 * std(V(fixatsamples));
	end


	theta = Tnew;
                                                                                   
function nm = find_next_min_below_threshold(s, win,thr)

 	pointer = 0;
	increment = 2;
	x = s;

	if length(s)<win
		nm = 0;
		return
	end

	thismin = find(s(1:win)==min(s(1:win)),1,'last');
	
	while (thismin==win) & (length(s)-win-increment) > 0 
	
	
		thismin = find(s(1:win)==min(s(1:win)),1,'last');
		if s(thismin)>thr
			break
		end

		pointer = pointer+increment;	
		s(1:increment) = [];

		% line(pointer, x(pointer),'marker', 'o')
		% line(pointer+thismin, x(pointer+thismin),'marker', 'o','color', 'r')
		% drawnow
		% pause(.2)

	end

	nm = pointer + thismin;

function nmx = find_next_max_above_threshold(s, win,thr)

	pointer = 0;
	increment = 1;
	x = s;

	if length(s)<win
		nmx = 0;
		return
	end

	thismax = find(s(1:win)==max(s(1:win)),1,'first');
		
	while (thismax==win) & ((length(s)-win-increment) > 0 ) 

		pointer = pointer+increment;
		s(1:increment) = [];

		thismax = find(s(1:win)==max(s(1:win)),1,'first');
	
		% line(pointer, x(pointer),'marker', 'o')
		% line(pointer+thismax, x(pointer+thismax),'marker', 'o','color', 'r')
		% drawnow
		% waitforbuttonpress

	end

	if s(thismax)<thr
		nmx = 0;
	else
		nmx = pointer + thismax - 1;
	end

function [problem_tag] = testSaccade(V_i, A_i, landmarks, opts)
	%TODO: build test for Rise/Fall < 1

		max_acceleration_accepted = opts.max_acceleration_accepted ;% = 40000;
		max_velocity_accepted 	  = opts.max_velocity_accepted 	   ;% = 1400;
		min_velocity_accepted     = opts.min_velocity_accepted     ;% = 200;
		theta_on = opts.theta_on;
		fs = opts.fs;
		detwin = opts.detectionwindow;

		% fall_div_rise = 1;
		
		accept_problematic = opts.accept_problematic;

		problem_tag = false(1,9);
		% 1 riseonset = 0
		% 2 minimum velocity not reached
		% 3 falloff faster than rise
		% 4 falloff very slow
		% 5 mean preceding velocity too large
		% 6 large velocities after offset
		% 7 too many peaks (>100deg/s) within detection window
		% 8 too small risetime
		% 9 saccade with negative duration?


		if landmarks(1) == 0
			problem_tag(1) = 1;
			return
		end

		if max(V_i) < min_velocity_accepted;
			problem_tag(2) = 1;
		end

		if (landmarks(4) - landmarks(2))/(landmarks(2)-landmarks(1)) < .5
			problem_tag(3) = 1;
		end

		if landmarks(2)-landmarks(1) > 4* (landmarks(4) - landmarks(2))
			problem_tag(4) = 1;
		end

		% preceding velocity too high
		if mean(V_i(1:landmarks(1)))>3*theta_on
			problem_tag(5) = 1;

		end

		if length(find(isnan(V_i)))>0
			length(find(isnan(V_i)))
			problem_tag(6) = 1;

		end
			
		if  length(findpeaks(V_i, 'minpeakheight', 100*fs/1e3, 'minpeakdistance', 50*fs/1e3))>3
			problem_tag(7) = 1;
		end

		if landmarks(2)-landmarks(1)<3
			problem_tag(8) = 1;
		end
		

		if landmarks(3)-landmarks(1)<0
			problem_tag(9) = 1;
		end
		


function out = repairSaccade(V_i, A_i, landmarks, issue) % TODO: might want to repair other problems

    	switch issue
    		case 'prec_vel_high'
    			warning('REPAIR: Trying to repair onset.')
		        lw = 100; % max upstroke duration
		        ht = 35; % high threshold , deg/s
		        lm = landmarks;
		        lm(5:6) = 0;
        
		        if lm(2)-lw < 0
		        	lw = lm(2)-1;
		        end

		     	lm_1 = [];
				while isempty(lm_1)
			 		lm_1 = lm(2) - find(V_i(lm(2):-1:lm(2)-lw)<ht ,1,'first')+1;
			 		warning('raising threshold')
			 		ht = ht+1;
				end
				lm(1) = lm_1;
			case 'argh'

		end
       

        out = lm;
        
function out = fit2gamma(V_i, landmarks, opts)

        fitemphasis = opts.fitemphasis;

		% samples between relevant landmarks (fitting part)
		
		switch fitemphasis
			case 'mainseq'
				safp = (landmarks(3)-landmarks(2));
			case 'peak-thresh'
				safp = (landmarks(4)-landmarks(2));
			case 'tail'
				if landmarks(6)~=0
					safp = (landmarks(6)-landmarks(2))-1;
				else
					safp = (landmarks(4)-landmarks(2))-1;
				end
		end

		% verify well formed saccades according to gamma fit residuals

		ydata = V_i(landmarks(1):landmarks(2)+safp)' - V_i(landmarks(1));
		 normfact = sum(ydata); 
		 ydatanorm=ydata/normfact;
		xdata = [1:landmarks(2)-landmarks(1)+safp+1];
		% fitfunction = @(p,x) p(1)*gampdf(x,p(2),p(3));
		fitfunction = @(p,x) gampdf(x,p(1),p(2));
		
		% fitspecs.x0 = [1  3  10];
		% fitspecs.lb = [.7  1  5];
		% fitspecs.ub = [2.5 10 15];
		fitspecs.x0 = [3  10];
		fitspecs.lb = [1  5];
		fitspecs.ub = [10 15];

		
		% sigmd = @(p,x) p(1) + p(2)./(1+exp(-(x-p(3))./p(4)));
		
		fitopts = optimset('lsqcurvefit');
		fitopts.Display = 'off';
		fitspecs.objective = fitfunction;
		fitspecs.xdata = xdata;
		fitspecs.ydata = ydatanorm;
		fitspecs.solver = 'lsqcurvefit';
		fitspecs.options = fitopts;

	
		[gfit.pset gfit.residualnorm gfit.residual gfit.flag gfit.output gfit.lambda] = lsqcurvefit(fitspecs);
		pset = gfit.pset;


		% 1/ detect saccades
		% 2/ find the residual for each
		% 3/ use a threshold based on the standard deviation of the residual of the population
		
		l = length(landmarks(1):landmarks(4));
		fittedsac_pred = fitfunction(gfit.pset,[1:l])';
		fittedsac_obs  = (V_i(landmarks(1):landmarks(4))- V_i(landmarks(1))) / normfact; 
		
						
		% compare the motion residual at the tail

	
		% [gof.h gof.pval gof.stats] = chi2gof(fittedsac_obs,'cdf', @(z)gamcdf(z, gfit.pset(1) ,gfit.pset(2)),'nparams',2 );
		% [gof.h gof.pval gof.stats] = chi2gof(fittedsac_obs,'expected', fittedsac_pred,'nparams',2 );


		% compare the motion residual at the tail
		
		residual = sum(sqrt((fittedsac_pred - fittedsac_obs).^2)/length(fittedsac_obs));
		
	
		out.predicted_fit = fitfunction(gfit.pset,[1:l])';
		out.observed_sac  = (V_i(landmarks(1):landmarks(4))- V_i(landmarks(1))) / normfact; 
		out.residual = sum(residual);
		out.fitspecs = fitspecs;
		out.pset  = pset;
		% out.gof = gof
		% gof.stats





