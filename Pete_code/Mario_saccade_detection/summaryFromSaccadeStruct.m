function sacs = summaryFromSaccadeStruct(varargin)

	p = inputParser;
	p.addRequired('sacs')
	p.addParamValue('selection', [])
	p.addParamValue('order_profiles', 'none')
	p.addParamValue('userselection', 'user' )
	p.addParamValue('ND_scatter_handle', [])
	p.addParamValue('to_plot', [1 1 1 1 1 1 1 1 1]) % 1/0 toggles plots
	p.addParamValue('summary_type','all')
	p.addParamValue('groups',[])
	p.addParamValue('colorscheme','light')
	

% function sacs = summaryFromSaccadeStruct(varargin)
% 
%(c) Mario Negrello - 12/2013
% mnegrello@gmail.com
% 
% IN:
% 'selection', []) - selected saccades for plotting.
% 'order_profiles', 'none') - how to order the velocity profiles in the waterfall plot. See parameters below .
% 'userselection', 'user' ) - 
% 'ND_scatter_handle', [])
% 'to_plot', [1 1 1 1 1 1 1 1 1]) % 1/0 
% 'summary_type','all')
% 'groups',[])
	
% paramter:
% 	to_plot accepts a vector with these plots(which scatter dimensions to plot)
	% 1 'onset-pk (ms)' 		
	% 2 'onset-offset (ms)'	 (duration)
	% 3 'amplitude (dg)'		
	% 4 'Rise / Fall (ms/ms)'
	% 5 'peak velocity (rad)'
	% 6 'x-y pk (ms)'			
	% 7 'residual (rad)'		
	% 8  'onset angle (rad)'	
	% 9 'tail duration (ms)'	
% 
%  order_profiles
% 				case 'vel'
% 					order by velocity
% 				case 'dur'
% 					by saccade duration (without tail)
% 				case 'dir'
% 					by onset direction
% 				case 'amp'
% 					by amplitude
% 				case 'res'
% 					by residual
% 				case 'rise'
% 					by rise time
% 				otherwise

% % 
% summary_type
% 		case 'all'
% 			sumplot = [1 1 1];
% 		case 'prof'
% 			sumplot = [0 1 0];
% 		case 'scat'
% 			sumplot = [1 0 0];
% 		case 'traces'
% 			sumplot = [0 0 1];
% 		otherwise
% 			sumplot = [0 0 0];

% 	end


	% p.addParamValue('categorize', 'tails')
	
	p.parse(varargin{:});

	sacs = p.Results.sacs;
	selection = p.Results.selection;
	order_profiles = p.Results.order_profiles;
	userselection = p.Results.userselection;
	to_plot = p.Results.to_plot;
	LM_v = sacs.saccades.landmarks.v;
	summary_type = p.Results.summary_type;
	groups = p.Results.groups;
	colorscheme = p.Results.colorscheme;



% colorscheme

global prefs
	prefs = [];

	prefs.colorscheme = colorscheme;


	switch prefs.colorscheme
		case 'light'
			prefs.figurecolor = [1 1 1];
			prefs.axescolor = [ 1 1 1];
			prefs.axiscolor = [0 0 0 ];
			prefs.markercolor = 'k';
		case 'dark'
			prefs.figurecolor = [0 0 0];
			prefs.axescolor = [0 0 0];
			prefs.axiscolor = [1 1 1];
			prefs.markercolor = 'w';
		otherwise
			keyboard
			% prefs.figurecolor = [1 1 1];
			% prefs.axescolor = [ 1 1 1];
			% prefs.axiscolor = [0 0 0 ];
			% prefs.markercolor = 'k';
	end




	%               __                        _          
	%   _________ _/ /____  ____ _____  _____(_)___  ___ 
	%  / ___/ __ `/ __/ _ \/ __ `/ __ \/ ___/ /_  / / _ \
	% / /__/ /_/ / /_/  __/ /_/ / /_/ / /  / / / /_/  __/
	% \___/\__,_/\__/\___/\__, /\____/_/  /_/ /___/\___/ 
	%                    /____/                          



	% select saccades
	switch userselection
		case 'user'
			if isempty(selection)
				disp('selected saccades missing, using all.')
				selected = [1:size(LM_v,1)];
			else	
				selected = selection;
			end
			marker = 'b.';
		case 'all'
			selected = [1:size(LM_v,1)];

		otherwise
			if isempty(selection)
				selected = [1:size(LM_v,1)];
				marker = 'k.';
			else
				selected = selection;
			end

	end




	%                    _                                  __      
	%   ____ ___________(_)___ _____  ____ ___  ___  ____  / /______
	%  / __ `/ ___/ ___/ / __ `/ __ \/ __ `__ \/ _ \/ __ \/ __/ ___/
	% / /_/ (__  |__  ) / /_/ / / / / / / / / /  __/ / / / /_(__  ) 
	% \__,_/____/____/_/\__, /_/ /_/_/ /_/ /_/\___/_/ /_/\__/____/  
	%                  /____/                                       

	LM_v = sacs.saccades.landmarks.v(selected,:);

	LM_v = LM_v(:,[2 3 4 5 6]) - repmat(LM_v(:,1),1,5);
	LM_v(find(LM_v<0)) = 0;

	% Duration
	ONtoPK  = LM_v(:,1);
	PKtoMAX = LM_v(:,2) - LM_v(:,1); % from peak to maximum deceleration
	PKtoOFF = LM_v(:,3) - LM_v(:,1);
	ONtoOFF = LM_v(:,3) - 0;
	
	TD = (LM_v(:,5) - LM_v(:,3)).*(LM_v(:,5)>0) ;


	% Rise and Fall
	RiseXFall    = ONtoPK./PKtoOFF;
	RiseXFallOff = ONtoPK./PKtoMAX; % should be above 1, the rise is always faster than the fall

	% V peaks
	PKV_v = sacs.saccades.peakV.v(selected,:);
	PKV_x = sacs.saccades.peakV.x(selected,:);
	PKV_y = sacs.saccades.peakV.y(selected,:);

	% if X Y peak disparity > 2*median OntoPk, disregard
	PK_xydif = PKV_x(:,2) - PKV_y(:,2);
	PK_xydif(find(PKV_x(:,1) < 2*median(ONtoPK))) = NaN;

	% Tail Amplitude
	TA = sacs.saccades.tail.amplitude(selected);

	% Angles
	ON_ang = sacs.saccades.angle.onset(selected);
	OFF_ang = sacs.saccades.angle.offset(selected);

	Origin_XY = sacs.saccades.origin;

	% For a histogram
	offsetthresh = sacs.saccades.threshold.offset(selected);
	fitresidual = sacs.saccades.fits.residual(selected);
			
	% For a histogram		
	SA = sacs.saccades.amplitude(selected,:);

	% Residuals
	RES = sacs.saccades.fits.residual(:,selected);

	% saccade velocity profiles
	vprofs   = sacs.saccades.vProfiles.v(selected,:);
	vprofs_x = sacs.saccades.vProfiles.x(selected,:);
	vprofs_y = sacs.saccades.vProfiles.y(selected,:);

	sacproblems = sacs.saccades.problematic;


	% WHATEVER CATEGORIES ONE MAY DESIRE COME HERE

	catfun = @(where, criterion) ( find(where >= criterion));


	if ~isempty(groups)
		
		categories = zeros(length(LM_v),1);
		for g = 1:length(groups)
		categories(groups{g}) = deal(g);
		end

		tagged = zeros(length(RES),3);
		tagged(find(LM_v(:,4) == 0),1)   = 1;
		tagged(find(LM_v(:,5) > 0),2)	 = 1;
		% tagged(find(RES < .35),3)		 = 1;
		tags = {'no tail' ;'tail' ;'well formed'};

	else

		categories(find(LM_v(:,5) == 0))	= 1;
		categories(find(LM_v(:,5) >  0))	= 2;
		
		tagged = zeros(length(RES),3);
		tagged(find(LM_v(:,4) == 0),1)   = 1;
		tagged(find(LM_v(:,5) > 0),2)	 = 1;
		% tagged(find(RES < .35),3)		 = 1;
		tags = {'no tail' ;'tail' ;'well formed'};


	end


	if ~exist('linspecer','file')
		prefs.groupcolors = hsv(length(unique(categories)));
	else
		prefs.groupcolors = linspecer(length(unique(categories)));
	end





	%     ____  ___  _________       __________     _____ _________  ______________________
	%     ____  ___  _________       __________     _____ _________  ______________________ 
	%    / __ \/   |/_  __/   |     /_  __/ __ \   / ___// ____/   |/_  __/_  __/ ____/ __ \
	%   / / / / /| | / / / /| |      / / / / / /   \__ \/ /   / /| | / /   / / / __/ / /_/ /
	%  / /_/ / ___ |/ / / ___ |     / / / /_/ /   ___/ / /___/ ___ |/ /   / / / /___/ _, _/ 
	% /_____/_/  |_/_/ /_/  |_|    /_/  \____/   /____/\____/_/  |_/_/   /_/ /_____/_/ |_|  
	% /_____/_/  |_/_/ /_/  |_|    /_/  \____/   /____/\____/_/  |_/_/   /_/ /_____/_/ |_| 

	
	data{1,1} = ONtoPK		;
	data{2,1} = ONtoOFF		;
	data{3,1} = SA(:,1)		;
	data{4,1} = RiseXFall 	;
	data{5,1} = PKV_v 		;
	data{6,1} = PK_xydif 	;
	data{7,1} = RES 		;
	data{8,1} = ON_ang 		;
	data{9,1} = TD 			;
	data{10,1} = TA 		;
	
	data{1,2} = { 'onset-pk (ms)' 		};
	data{2,2} = { 'onset-offset (ms)'	};
	data{3,2} = { 'amplitude (dg)'		};
	data{4,2} = { 'Rise / Fall (ms/ms)'	};
	data{5,2} = { 'peak velocity (rad)'	};
	data{6,2} = { 'x-y pk (ms)'			};
	data{7,2} = { 'residual (rad)'		};
	data{8,2} = { 'onset angle (rad)'	};
	data{9,2} = { 'tail duration (ms)'	};
	data{10,2} = { 'tail amplitude (deg)'};


	dd = 0;
	for ds = 1:9
		if to_plot(ds)
			dd = dd+1;
			datatoplot{dd,1}  = data{ds,1};
			datatoplot{dd,2}  = data{ds,2};
		end
	end

	switch summary_type
		case 'all'
			sumplot = [1 1 1];
		case 'prof'
			sumplot = [0 1 0];
		case 'scat'
			sumplot = [1 0 0];
		case 'traces'
			sumplot = [0 0 1];
		otherwise
			sumplot = [0 0 0];
	end



	
	% annotate = {['trial duration: ' num2str(sacs.trialDuration)];
	 			% ['saccade count: ' num2str(sacs.numberOfSaccades)];
	 			% ['well formed:' num2str(length(find(RES<.35)))];
	 			% ['SN/R:' num2str(sacs.summary.SNR)]};

	
	% set(hands(3,8), 'xlim', [0 1], 'ylim', [0 1])
	% text(hands(1,8), 0, .5, annotate)
	% legend(hands(2,1), {'no tail' ; 'tail' ; 'high residual'})

	%            __      __                      _____ __         
	%     ____  / /___  / /_   ____  _________  / __(_) /__  _____
	%    / __ \/ / __ \/ __/  / __ \/ ___/ __ \/ /_/ / / _ \/ ___/
	%   / /_/ / / /_/ / /_   / /_/ / /  / /_/ / __/ / /  __(__  ) 
	%  / .___/_/\____/\__/  / .___/_/   \____/_/ /_/_/\___/____/  
	% /_/                  /_/                                    


	switch order_profiles
				case 'vel'
				[v o] = sort(PKV_v);
				case 'dur'
				[v o] = sort(ONtoOFF);
				case 'dir'
				[v o] = sort(ON_ang);
				case 'amp'
				[v o] = sort(SA(:,1));
				case 'res'
				[v o] = sort(residual_motion);
				case 'rise'
				[v o] = sort(RiseXFall);
				otherwise
				o = selected;
	end



	vprofs_sorted = vprofs(o, :);
	
	fvd = zeros(size(RES(o)));
	fvd(find(RES(o)<=.2)) 				= 0;
	fvd(find(RES(o)>.2 & RES(o)<.5)) 	= .2;
	fvd(find(RES(o)>.5 & RES(o)<.8)) 	= .8;
	fvd(find(RES(o)>.8)) 				= 1;
	fvd= fvd';

	
	residual_cmap = [fvd ones(length(fvd),1)*.3 ones(length(fvd),1)*.3];

	sacproblems = sacproblems(o,:);

	% scatter
	if sumplot(1)
		fighand(1) = figure('color', prefs.figurecolor);
		hands = ND_scat_subplot(datatoplot , categories);
	end


	% profiles
	if sumplot(2)
		fighand(2) = figure('color', prefs.figurecolor);
	
		landmarks = sacs.saccades.landmarks.v(o,:);
		plot_profiles( vprofs_sorted, residual_cmap , fighand(2) );
		plot_problems( fighand(2), sacproblems);


	end

	% traces
	if sumplot(3)
		fighand(3) = figure('color', prefs.figurecolor);

		f3a(1) = subplot(2,2,1);
			
			rose(ON_ang)
			title('saccade counts per direction')
			set(f3a(1),'xcolor', prefs.axiscolor,'ycolor', prefs.axiscolor,'color',prefs.axescolor)

		f3a(2) = subplot(2,2,2);
			
			for g = length(unique(categories))
				line(SA(find(categories==1),2)', SA(find(categories==1),3)','color', prefs.groupcolors(g,:),'linestyle', 'none')
			end
			title('saccade counts per direction')
			set(f3a(2),'xcolor', prefs.axiscolor,'ycolor', prefs.axiscolor,'color',prefs.axescolor)


		f3a(3) = subplot(2,2,3);
			
		% 	hold on
		% 	compass(f3a2, SA(find(categories==1),2)', SA(find(categories==1),3)','g')
		% 	compass(f3a2, SA(find(categories==2),2)', SA(find(categories==2),3)','y')
		% 	compass(f3a2, SA(find(categories==3),2)', SA(find(categories==3),3)','r')
		% 	hold off
		% 	title('saccade directions')
			fs = 1000;

				xx = cumsum(vprofs_x,2)/fs ; xx = xx(:);
				yy = cumsum(vprofs_y,2)/fs ; yy = yy(:);
				vv = vprofs   ; vv = vv(:)/fs;
				
				
				scatter(xx,yy,log(vv+10),vv,'filled') 
				colorbar
				title('saccades and instantaneous velocity')
				xlabel('H deg')
				ylabel('V deg')
				axis([-18 18 -18 18])
				axis equal
				set(f3a(3),'xcolor', prefs.axiscolor,'ycolor', prefs.axiscolor,'color',prefs.axescolor)

		f3a(4) = subplot(2,2,[4]);
			
			
			
				prof_x = vprofs_x/fs; prof_x(:,1) = prof_x(:,1)+ Origin_XY(:,1);
				prof_y = vprofs_y/fs; prof_y(:,1) = prof_y(:,1)+ Origin_XY(:,2);
				xx = cumsum(prof_x,2)  ; xx = xx(:);
				yy = cumsum(prof_y,2)  ; yy = yy(:);
				vv = vprofs   ; vv = vv(:);
				

				scatter(xx,yy,log(vv+10),vv,'filled') 
				colorbar
				title('saccades and instantaneous velocity')
				xlabel('H deg')
				ylabel('V deg')
				axis([-18 18 -18 18])
				axis equal
				set(f3a(4),'xcolor', prefs.axiscolor,'ycolor', prefs.axiscolor,'color',prefs.axescolor)

				% line(, vprofs_y(find(categories==1),:)','color','g');
				% line(vprofs_x(find(categories==2),:)', vprofs_y(find(categories==2),:)','color','y');
				% line(vprofs_x(find(categories==3),:)', vprofs_y(find(categories==3),:)','color','r');
				% title('velocity profiles')
		
		
		% f3a4 = subplot(2,2,[2]);
		% 		line(cumsum(vprofs_x(find(categories==1),:),2)', cumsum(vprofs_y(find(categories==1),:),2)','color','g');
		% 		line(cumsum(vprofs_x(find(categories==2),:),2)', cumsum(vprofs_y(find(categories==2),:),2)','color','y');
		% 		line(cumsum(vprofs_x(find(categories==3),:),2)', cumsum(vprofs_y(find(categories==3),:),2)','color','r');
		% 		title('saccade traces')


	end

	% if sumplot(4)
	% end





		%    _________ ___________   ____  __  __/ /_____  __  __/ /_
		%   / ___/ __ `/ ___/ ___/  / __ \/ / / / __/ __ \/ / / / __/
		%  (__  ) /_/ / /__(__  )  / /_/ / /_/ / /_/ /_/ / /_/ / /_  
		% /____/\__,_/\___/____/   \____/\__,_/\__/ .___/\__,_/\__/  
		%                                        /_/                 


% data{1,1} = ONtoPK	
% data{2,1} = ONtoOFF	
% data{3,1} = SA 		
% data{4,1} = RiseXFall
% data{5,1} = PKV_v 	
% data{6,1} = PK_xydif 
% data{7,1} = RES 	
% data{8,1} = ON_ang 	
% data{9,1} = TD 		



	sacs.saccades.categories = tagged;
	sacs.saccades.category_names = tags;

	


end

	%     __  __     __                   ______                 __  _                 
	%    / / / /__  / /___  ___  _____   / ____/_  ______  _____/ /_(_)___  ____  _____
	%    / / / /__  / /___  ___  _____   / ____/_  ______  _____/ /_(_)___  ____  _____
	%   / /_/ / _ \/ / __ \/ _ \/ ___/  / /_  / / / / __ \/ ___/ __/ / __ \/ __ \/ ___/
	%  / __  /  __/ / /_/ /  __/ /     / __/ / /_/ / / / / /__/ /_/ / /_/ / / / (__  ) 
	% /_/ /_/\___/_/ .___/\___/_/     /_/    \__,_/_/ /_/\___/\__/_/\____/_/ /_/____/  
	% /_/ /_/\___/_/ .___/\___/_/     /_/    \__,_/_/ /_/\___/\__/_/\____/_/ /_/____/  
	%             /_/                                                                  

function hand = ND_scat_subplot(data, groups)

	global prefs
	ngroups = length(unique(groups));
	groupcolors = prefs.groupcolors;


	N = length(data);

	l1 = linspace(0, .9, N+1); l1 = l1(1:end -1)+.05; 
	l2 = linspace(0, .9, N+1); l2 = l2(1:end -1)+.05; l2 = fliplr(l2);
	w  = .9/N;
	h  = .9/N;

	r = 0; c = 0;
	for rr = l1
		r = r+1;
		for cc = l2
			c = c+1;
			
			hand(r, c) = axes('position', [cc rr w h]);
			if r > c
				hold on
				for g = ngroups:-1:1
					line(data{c,1}(groups==g), data{r,1}(groups==g), ...
						 'markersize', 5, 'color', groupcolors(g,:), 'marker','.','linestyle','none','markersize',15)
				end
				hold off
			elseif r == c
				[hhh xxx] = hist(data{c,1},30);
				bar(xxx, hhh, 'facecolor',prefs.markercolor)
			end

			if c > r
				line(data{c,1}, data{r,1} , 'color', prefs.markercolor, 'marker','.','linestyle','none','markersize', 5)
			end


			if r == 1
				title(data{c,2})
			end
			if c == N
				ylabel(data{r,2})
			end
			if r ~= 1
				set(gca, 'xtick' , [])
			end
			if c ~= 1
				set(gca, 'ytick' , [])
			end
		
			set(gca,'color',prefs.axescolor,'xcolor', prefs.axiscolor, 'ycolor', prefs.axiscolor);
		end

		c = 0;
	end

	colormap(jet(5))

	pl_i = [1:length(l1)];
	for r = pl_i

		linkaxes(hand(r,pl_i(find(pl_i~=r))),'y')

	end



end

function f = plot_profiles(vprofs, colors, figh)

		global prefs

		f = figure(figh)

		h = waterfall(vprofs(:,:));
		set(gca,'color',prefs.axescolor,'xcolor',prefs.axiscolor,'ycolor',prefs.axiscolor)
		set(h, 'FaceColor', 'flat');
		set(h, 'FaceVertexCData', colors)
		set(h, 'FaceAlpha', 0.7);
		set(h, 'EdgeColor', prefs.markercolor);
		grid off
		view(39,57)
		line([0 0], [0 0], [0 100],'color', prefs.markercolor, 'linewidth', 10)


end


function plot_problems(figh, problemlist, vprofiles)

	ns = length(problemlist);
	cl = jet(15);

	for p = 1:size(problemlist,2)
		y = find(problemlist(:,p));
		ny = length(y);
		if ny>1

			line ( ones(ny,1) * (250+ (2*p) ),  y, zeros(ny,1) ,'linestyle', 'none', 'marker', '.' ,'color' ,cl(p,:)  )
		end

	end


end


% function hand = ND_regressions(N)

% titles = {'1';'1';'1';'1';'1'};

% l1 = linspace(0, .9, N+1); l1 = l1(1:end -1)+.05;
% l2 = linspace(0, .9, N+1); l2 = l2(1:end -1)+.05;
% w  = .9/N;
% h  = .9/N;

% r = 0; c = 0;
% for rr = l1
% 	r = r+1;
% 	for cc = l2
% 		c = c+1;
		
% 		hand(r, c) = axes('position', [cc rr w h]);

% 		scatter(data{r},data{c})

% 		if r == N
% 			title(titles{r})
% 		end
% 		if c == 1
% 			ylabel(titles{r})
% 		end
% 		if r ~= 1
% 			set(gca, 'xtick' , [])
% 		end
% 		if c ~= 1
% 			set(gca, 'ytick' , [])
% 		end
% 	end
% 	c = 0;
% end




