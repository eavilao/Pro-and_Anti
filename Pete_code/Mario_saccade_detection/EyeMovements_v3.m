%%==============================================================================
%%                                                   Eye Movement Reconstruction
%%==============================================================================

function [M] = EyeMovements_v3(varargin)


p = inputParser;
p.addRequired('xy')  % a matrix with two columns or a cell array with two cells;
p.addOptional('lfp',[])  % a matrix with two columns or a cell array with two cells;
p.addParamValue('Fs', 1000) % stdandard deviation criterion for offset threshold

p.parse(varargin{:});

xy = p.Results.xy;
Fs = p.Results.Fs;
lfp = p.Results.lfp;




%comment this out if no display_figuredefaults.m present
% display_figuredefaults.m; 

[r c] = size(xy);
if r==2 & c>=2
	xy = xy';
elseif r>=2 & c==2
	xy = xy;
else
	disp('unrecognized input format.')
end


colors = linspecer(5);

%%===============================================================     Parameters

snakesize = 200; %in ms if Fs = 1000


initialtime = 4000;
lengthofaxis = 1.5;

steps = 100; %in ms after resampling (42 = ~24fps)
L = 1 % proportion of the file to reconstruct: LengthOfTrial * L ( 0 < L < 1)
 
genMovie = false; % whether to output a movie variable
saveMovie = false; % whether to generate an avi


steps = 7 ; %in ms after resampling (42 = ~24fps)
self.steps = 0;
self.update = false;

L = 1; % proportion of the file to reconstruct: LengthOfTrial * L ( 0 < L < 1)

%%==================================================     Define Snake Body Parts

tailLength = snakesize*.4;
lowerBodyLength = snakesize*.3 + tailLength;
upperBodyLength = snakesize*.2 + lowerBodyLength;
headLength = snakesize*.1 + upperBodyLength;

magnifier = [1:snakesize];

%%==========================================================     Gather Eye Data

lengthOfTrial = length(xy(:,1));
% define period of reconstruction
lengthOfTrial = round(lengthOfTrial*L);


maxx  = max(xy(:,1));
maxy  = max(xy(:,2));
minx  = min(xy(:,1));
miny  = min(xy(:,2));


%%=====================================================     Compute Eye Velocity


eyeDisplacement = diff(xy);
eyeVelocity = hypot(eyeDisplacement(:,1), eyeDisplacement(:,2));



eyeVelocity(1:5) = 0;
minv = min(eyeVelocity);
maxv = max(eyeVelocity);
x = [1:length(eyeVelocity)]/Fs;





if ~exist('fig')
	figure1 = figure;
	set(figure1, 'position', [ 500 500 500 500])
	disp('Adjust position and press some key.')
	pause;
else
	clf
	disp('Adjust position and press some key.')
	pause;
	
end
	clc
	disp('==========  Keyboard Controls  ==========')
	disp('[+] increases the size of the step')
	disp('[0] keeps the size of the current step.')
	disp('[-] decreases the size of the time step')
	disp('[p] pauses.')
	disp('[q] quits.')
	disp('==========  Keyboard Controls  ==========')

%%====================================================     Callback for Keypress
set(gcf,'KeyPressFcn',@keypress);
set(gcf,'userdata',self);
% set(gcf,'Interruptible','on');


%%==============================================================     create Axes

% Create SnakeAxes axes
p0 =  [0.02         0.05         0.45       0.9123];
% SnakeAxes = axes('Parent',figure1,'YTick',zeros(1,0),'XTick',zeros(1,0),...
    % 'Position',p0);
SnakeAxes = axes('Parent',figure1,...
    'Position',p0);


box(SnakeAxes,'on');

% % Create Autocorr Axes
% ACAxes = axes('Parent',figure1,'YAxisLocation','right',...
%     'Position',[0.72 0.11 0.25 0.39]);
% box(ACAxes,'on');

% Create LinVelAxes

p1 = [0.47673      0.71974      0.45346      0.17974];
EyeVelAxes = axes('Parent',figure1,'YTick',1,'YAxisLocation','right',...
    'XTick',zeros(1,0),...
    'TickDir','in',...
    'TickLength',[0.001 0.025],...
    'Position',p1);
box(EyeVelAxes,'on');


p2 = [0.47673      0.53576      0.45227      0.17524];
% Create EyePositionAxes
EyePositionAxes = axes('Parent',figure1,'YTick',zeros(1,0),...
    'TickLength',[0.001 0.0025],...
    'Position',p2,...
    'FontSize',11);
box(EyePositionAxes,'on');


if ~isempty(lfp)

	p3 = [ 0.47673         0.11      0.45346         0.39];
	% Create lfpAxes
	lfpAxes = axes('Parent',figure1,'YTick',1,'YAxisLocation','right',...
	    'TickDir','in',...
	    'TickLength',[0.001 0.025],...
	    'Position',p3);
	box(lfpAxes,'on');
	xlabel('Time (s)')

end

%%==============================================================================
%%                                                         First Iteration Plots
%%==============================================================================


%%============================================================     Plot Autocorr
% axes(ACAxes)
% grid
% set(gca,'lineWidth',1.5)
% acax = line(ACT, ACR(1,:), 'color', [.35 .57 .671], 'linewidth', 1.2);
% axis([min(ACT) max(ACT) 100 1000])



%%============================================================     Subplot Snake

axes(SnakeAxes)
set(gca,'lineWidth',3)
% s2 = subplot('position', [0.55 0.60 0.34 0.33])
% subplot(3,1,1)
	axis([minx maxx miny maxy])
	% set(gca, 'xtick', [],'ytick', [], 'box' , 'on')
	set(gca, 'box' , 'on')
		t = 1;
		
		tail = line(xy(t:t+tailLength,1), xy(t:t+tailLength,2), ...
			'color', colors(2,:), 'linewidth', 1);
	
		lowbody = line(xy(t+tailLength+1:t+lowerBodyLength,1), xy(t+tailLength+1:t+lowerBodyLength,2), ...
			'color', colors(3,:) , 'linewidth', 1);
	
		upbody = line(xy(t+lowerBodyLength+1:t+upperBodyLength,1), xy(t+lowerBodyLength+1:t+upperBodyLength,2), ...
			'color', colors(4,:), 'linewidth', 3);
	
		head = line(xy(t+upperBodyLength+1:t+headLength,1), xy(t+upperBodyLength+1:t+headLength,2), ...
			'color', colors(5,:), 'linewidth', 4);

			xlabel('Eye Axis X')
			ylabel('Eye Axis Y')
			axis square
		titl = title([ num2str(t)  ' ms']);


%%=========================================================     Subplot Velocity

axes(EyeVelAxes)
set(gca,'lineWidth',3)
% s3 = subplot('position', [0.05 0.32 0.9 0.25])
% s3 = subplot(3,1,2)


	set(gca, 'xtick', [],'ytick', [], 'box' , 'on')
	bg  = line(x(1:lengthofaxis*headLength),eyeVelocity(1:lengthofaxis*headLength),'color', [.659 .533  .353]);
	hair = line([headLength headLength]/Fs,[-.1  maxv],'color', [.6 .6 .6]);

	
	vline  = line(x(1:lengthofaxis*headLength+1),ones(1,lengthofaxis*headLength+1)*.02,'color', [.659 .533  .353]);

	axis tight
	% set(s2, 'position', [ 0.1 0.4 0.8 0.3])	
	
%%=========================================================     Subplot Position

axes(EyePositionAxes)
set(gca,'lineWidth',3)

% s4 = subplot('position', [0.05 0.05 0.9 0.25])
% subplot(3,1,3)
		% set(fig,'CurrentAxes',et)
	set(gca, 'box' , 'on','ytick', [])
	% set(gca, 'position', [ 0.1 0.1 0.8 0.3])


	ex  = line(x(1:lengthofaxis*headLength),xy(1:lengthofaxis*headLength,1),'color', [.655 .45  .67],'lineWidth',2);
	ey  = line(x(1:lengthofaxis*headLength),xy(1:lengthofaxis*headLength,2),'color', [.55 .59 .70],'lineWidth',2);
	
	hair2 = line([headLength headLength]/Fs,[minx maxx],'color', [.4 .4 .4]);
	% xlabel('time (s)')
	axis tight

%%==============================================================     Subplot lfp

if ~isempty(lfp)
	axes(lfpAxes)
	set(gca,'lineWidth',3)
		lfpline = line(magnifier/Fs, lfp(magnifier),'color', [1 .81 .73]);
		lfpMinYaxis = min(lfp);
		lfpMaxYaxis = max(lfp);
		set(gca, 'ylim', [lfpMinYaxis lfpMaxYaxis])
end


%%===============================================     Annotation of Plot Details 
% details = annotation(figure1,'textbox',...
%     [0.029 0.52 0.27 0.07],...
%     'String',{['step: ' ],['cumsteps: ' ],['Index of AC: ' ]},...
%     'FitBoxToText','off');


		
%%==============================================================================
%%                                                               Main Loopy Loop
%%==============================================================================
cumsteps = 0;
nextAcUpdate = 2;
f = 0;
mf = 0;
idx = 1;

t=2;

try
	while t < lengthOfTrial-snakesize
	
		%%====================================================     Act upon the callback

		self = get(gcf,'userdata');
		
		if self.update == true
			cumsteps = cumsteps + self.steps;
			self.update = false;
			set(gcf,'userdata', self)
		elseif self.update =='pause'
			pause
		elseif self.update =='quit'
			break
		end
		t = t+ cumsteps;
		
		%  set(details, 'string',{['size of step: : ' num2str(cumsteps) 'ms'],['update: ' num2str(self.update)]} )
		% 	
		% %%=====================================     Autocorrelation Figure Update
		% 
		% 
		% if abs(t-ACIntervals(idx)) > res/2			
		% 	idx = find(ACIntervals > t , 1, 'first');
		% 	if isempty(idx)
		% 		idx = 1;
		% 	elseif length(idx)>1
		% 		idx = idx(1);
		% 	end
		% 	axes(ACAxes)
		% 	set(acax, 'ydata', ACR(idx,:))
		% 
		% end

	
		%%===========================================================      Snake
		
		f = f+1;
		axes(SnakeAxes)
		
		tailx = xy(t:t+tailLength,1);
				taily = xy(t:t+tailLength,2);
			
				lowbodyx =xy(t+tailLength+1:t+lowerBodyLength,1); 
				lowbodyy =xy(t+tailLength+1:t+lowerBodyLength,2);
			
				upbodyx = xy(t+lowerBodyLength+1:t+upperBodyLength,1);
				upbodyy = xy(t+lowerBodyLength+1:t+upperBodyLength,2);
			
				headx =   xy(t+upperBodyLength+1:t+headLength,1);
				heady =   xy(t+upperBodyLength+1:t+headLength,2);
			
				
				set(tail, 'Xdata',tailx, 'Ydata', taily);
				set(lowbody, 'Xdata',lowbodyx, 'Ydata', lowbodyy);
				set(upbody,'Xdata',upbodyx, 'Ydata', upbodyy);
				set(head,'Xdata',headx, 'Ydata', heady);
			
				set(titl, 'string', [ num2str(t+headLength)  ' ms'])
	
	%%====================================================== Update:     Eye Velocity
		axes(EyeVelAxes)

		xdatabg  = x(t:t+lengthofaxis*headLength);
		ydatabg  = eyeVelocity(t:t+lengthofaxis*headLength);
		% hair = line([t+headLength t+headLength]/Fs,[0 (meanv + 3*stdv)]);


		% xlabel('time (s)')

		set(bg, 'Xdata',xdatabg, 'Ydata', ydatabg);
		set(hair, 'Xdata',[t+headLength t+headLength]/Fs);
		set(vline,'Xdata',xdatabg)
		axis tight
		
	%%====================================================== Update:     Eye Position
		axes(EyePositionAxes)

		exxdata  = x(t:t+lengthofaxis*headLength);
		exydata  = xy(t:t+lengthofaxis*headLength,1);

		eyxdata  = x(t:t+lengthofaxis*headLength);
		eyydata  = xy(t:t+lengthofaxis*headLength,2);

		set(hair2, 'Xdata',[t+headLength t+headLength]/Fs);
		set(ex, 'Xdata',exxdata, 'Ydata', exydata);
		set(ey, 'Xdata',eyxdata, 'Ydata', eyydata);
		axis tight		
	
	%%==========================================================     Update lfp Plot 
		if ~isempty(lfp)
			axes(lfpAxes)
			set(lfpline, 'Xdata',(magnifier+t)/Fs ,'Ydata', lfp(magnifier+t));
			% axis tight
			% keyboard
			drawnow

			
			set(gca, 'xlim', [min(magnifier+t)/Fs max(magnifier+t)/Fs])
		end
	

		makemovie = 1;
		if makemovie
			if self.steps >0
			mf = mf +1;
			M(mf) = getframe;
			end
		end



	end



	% M(t) = getframe;

catch E
	keyboard

	% Try catch block to exit cleanly if user aborts execution
    if ~strcmp(E.identifier, 'MATLAB:class:InvalidHandle')||~strcmp(E.identifier, 'MATLAB:capture:FigureDestroyedDuringGetFrame')
        
        % rethrow(E);
	else
        rethrow(E);
		fig = gcf;
		if exist('fig')
			close(fig)
		end
    end
end

fig = gcf;
if exist('fig')
	close(fig)
end






%%==============================================================================
%%                                                         On Key Press Function
%%==============================================================================


function keypress(src,evnt)


increment = 1;

self = get(gcf,'userdata');

if isempty(evnt.Character)
	self.steps = 0;
	self.update = false;
    % return
end


if evnt.Character == '+' | evnt.Character == '='
    self.steps = increment;
	self.update = true;
    recognized = 1;
    set(gcf,'userdata',self);

elseif evnt.Character == '-' | evnt.Character == '_'
    self.steps = -increment;
	self.update = true;
    recognized = 1;
    set(gcf,'userdata',self);

elseif evnt.Character == '0'
    self.steps = 0;
	self.update = true;
    recognized = 1;
    set(gcf,'userdata',self);

elseif strcmp(evnt.Character, 'p')
	
	self.update = 'pause';
    set(gcf,'userdata',self);

elseif evnt.Character == 'q'

	self.update = 'quit';

    % close(gcf); 
   
end















% 
% %%==================================================     Define Snake Body Parts
% 
% 
% 
% 
% tailLength = snakesize*.4;
% lowerBodyLength = snakesize*.3 + tailLength;
% upperBodyLength = snakesize*.2 + lowerBodyLength;
% headLength = snakesize*.1 + upperBodyLength;
% 
% %%=====================================================     Compute Eye Velocity
% 
% 
% 
% 
% 	eyeDisplacement = sqrt(xy(:,1).^2+xy(:,2).^2);
% 	eyeDisplacement = conv(eyeDisplacement, gausswin(50))/sum(gausswin(50));
% 	eyeVelocity = abs(diff(eyeDisplacement));
% 	eyeVelocity(1:5) = 0;
% 	stdv = std(eyeVelocity);
% 	meanv = mean(eyeVelocity);
% 	x = [1:length(eyeVelocity)]/Fs;
% 	
% %%=============================================================     Prepare Plot
% 
% t=1;
% 
% subplot(1,3,1)
% 	axis([minx maxx miny maxy])
% 	set(gca, 'xtick', [],'ytick', [], 'box' , 'on')
% 	
% 
% 
% subplot(1,3,2)
% 	bg  = line(x(t:t+10*headLength),eyeVelocity(t:t+10*headLength),'color', [1 0 0]);
% 	hair = line([t+headLength t+headLength]/Fs,[0 (meanv + 3*stdv)]);
% 	axis tight
% 	xlabel('time (s)')	
% 	
% subplot(1,3,3)
% 		% set(fig,'CurrentAxes',et)
% 
% 	ex  = line(x(t:t+10*headLength),xy(:,1)(t:t+10*headLength),'color', [1 0 0]);
% 	ey  = line(x(t:t+10*headLength),xy(:,2)(t:t+10*headLength),'color', [1 1 0]);
% 	
% 	hair2 = line([t+headLength t+headLength]/Fs,[minx maxx]);
% 	axis tight
% 	xlabel('time (s)')
% 
% 	
% 		
% %%==============================================================================
% %%                                                               Main Loopy Loop
% %%==============================================================================
% 
% f = 0;
% try
% 	for t = [initialtime: steps: lengthOfTrial-snakesize]
% 	
% 	if gcf ~= fig
% 		break
% 	end
% 		
% 	
% 	
% 	%%====================================================================     Snake
% 		
% 	subplot(3,1,1)
% 		if ~mod(f,100)
% 			cla
% 		end
% 
% 		f = f+1;
% 
% 
% 		% 
% 		% line(xy(:,1)(t:t+tailLength), xy(:,2)(t:t+tailLength), ...
% 		% 	'color', colors(2,:), 'linewidth', 1);
% 		% 	
% 		% line(xy(:,1)(t+tailLength+1:t+lowerBodyLength), xy(:,2)(t+tailLength+1:t+lowerBodyLength), ...
% 		% 	'color', colors(3,:) , 'linewidth', 1);
% 		% 	
% 		% line(xy(:,1)(t+lowerBodyLength+1:t+upperBodyLength), xy(:,2)(t+lowerBodyLength+1:t+upperBodyLength), ...
% 		% 	'color', colors(4,:), 'linewidth', 3);
% 		% 	
% 		% line(xy(:,1)(t+upperBodyLength+1:t+headLength), xy(:,2)(t+upperBodyLength+1:t+headLength), ...
% 		% 	'color', colors(5,:), 'linewidth', 4);
% 		% line(xy(:,1)(t+headLength:t+headLength+1), xy(:,2)(t+headLength:t+headLength+1), ...
% 		% 		'color', [1 .1 .1], 'linewidth', 4.5);
% 		% 
% 		% drawnow
% 			
% 			
% 			
% 			tail = line(xy(:,1)(t:t+tailLength), xy(:,2)(t:t+tailLength), ...
% 				'color', colors(2,:), 'linewidth', 1);
% 				
% 			lowbody = line(xy(:,1)(t+tailLength+1:t+lowerBodyLength), xy(:,2)(t+tailLength+1:t+lowerBodyLength), ...
% 				'color', colors(3,:) , 'linewidth', 1);
% 				
% 			upbody = line(xy(:,1)(t+lowerBodyLength+1:t+upperBodyLength), xy(:,2)(t+lowerBodyLength+1:t+upperBodyLength), ...
% 				'color', colors(4,:), 'linewidth', 3);
% 				
% 			head = line(xy(:,1)(t+upperBodyLength+1:t+headLength), xy(:,2)(t+upperBodyLength+1:t+headLength), ...
% 				'color', colors(5,:), 'linewidth', 4);
% 
% 			xlabel('Eye Axis X')
% 			ylabel('Eye Axis Y')
% 			axis square
% 		title([ num2str(t)  ' ms'])
% 
% 
% 	%%====================================================== Update:     Eye Velocity
% 	
% 	subplot(3,1,2)
%  		% set(fig,'CurrentAxes',et)
% 
% 		xdatabg  = x(t:t+10*headLength);
% 		ydatabg  = eyeVelocity(t:t+10*headLength);
% 		% hair = line([t+headLength t+headLength]/Fs,[0 (meanv + 3*stdv)]);
% 
% 		% xlabel('time (s)')
% 
% 		set(bg, 'Xdata',xdatabg, 'Ydata', ydatabg);
% 		set(hair, 'Xdata',[t+headLength t+headLength]/Fs);
% 		axis tight
% 		
% 	%%====================================================== Update:     Eye Position
% 	
% 	subplot(3,1,3)
%  		% set(fig,'CurrentAxes',et)
% 
% 		% ex  = line(x(t:t+10*headLength),xy(:,1)(t:t+10*headLength),'color', [1 0 0]);
% 		% ey  = line(x(t:t+10*headLength),xy(:,2)(t:t+10*headLength),'color', [1 1 0]);
% 		exxdata  = x(t:t+10*headLength);
% 		exydata  = xy(:,1)(t:t+10*headLength);
% 
% 		eyxdata  = x(t:t+10*headLength);
% 		eyydata  = xy(:,2)(t:t+10*headLength);
% 
% 		
% 		% hair2 = line([t+headLength t+headLength]/Fs,[minx maxx]);
% 
% 		% xlabel('time (s)')
% 
% 		set(hair2, 'Xdata',[t+headLength t+headLength]/Fs);
% 		set(ex, 'Xdata',exxdata, 'Ydata', exydata);
% 		set(ey, 'Xdata',eyxdata, 'Ydata', eyydata);
% 		axis tight		
% 		
% 
% 		drawnow
% 
% 		if genMovie | saveMovie
% 		M(f) = getFrame;
% 		end
% 
% 	end
% 
% catch E
% 	% Try catch block to exit cleanly if user aborts execution
%     if ~strcmp(E.identifier, 'MATLAB:class:InvalidHandle')||~strcmp(E.identifier, 'MATLAB:capture:FigureDestroyedDuringGetFrame')
%         rethrow(E);
% 	else
% 		fig = gcf;
% 		if exist('fig')
% 			close(fig)
% 		end
%     end
% end
% 
% fig = gcf;
% if exist('fig')
% 	close(fig)
% end
% 
% 
% 
% 
