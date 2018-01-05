function sacs = statisticsFromSaccadeStruct(varargin)
% in the plan: pass categories and make distributions for each #TODO


	p = inputParser;
	p.addRequired('sacs')
	
	p.parse(varargin{:});

	sacs = p.Results.sacs;
	
	%   ____ ___________(_)___ _____  ____ ___  ___  ____  / /______
	%  / __ `/ ___/ ___/ / __ `/ __ \/ __ `__ \/ _ \/ __ \/ __/ ___/
	% / /_/ (__  |__  ) / /_/ / / / / / / / / /  __/ / / / /_(__  ) 
	% \__,_/____/____/_/\__, /_/ /_/_/ /_/ /_/\___/_/ /_/\__/____/  
	%                  /____/                                       

	LM_v = sacs.saccades.landmarks.v;
	LM_v = LM_v(:,[2 3 4 5 6]) - repmat(LM_v(:,1),1,5);
	LM_v(find(LM_v<0)) = 0;


	% stats.tail.ratio = 
	% stats.tail.duration = 
	% stats.tail.amplitude = 
	

% stats.boot.tail_duration = bootci(100, @(x) [mean(x) std(x)], S.saccades.tail.duration(find(S.saccades.tail.duration)));
% stats.boot.tail_duration = bootci(100, @(x) [mean(x) std(x)], S.saccades.tail.duration(find(S.saccades.tail.duration)));
