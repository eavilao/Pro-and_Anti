function [out] = trial_parser(txtfile)
%  V.0001 - parses the txt file with the exported variables from eye link data. 
%  Makes VERY rigid assumptions about the order of columns in C.
% 
%  MSG to GTS: Esta é uma mensagem autoreferencial!
% 
% Aqui que vamos explicar como a funcao funciona!
% 
% to calculate positions we assume a screensize of 1024 x 768.
% 
% 
% returns timestamps of successful trials
% (resets first trial as t = 0)
% 
%  TODO: improve mapping between pixels and degrees - currently, systematic distortions are expected!
% We've observed that the resolution distortion is within the realm of the reasonable. 12deg/px -> 16deg/px.
 % We are using the mean. We incur on a mild overestimation of the length of amplitude vectors when close to the center.


%1 RECORDING_SESSION_LABEL	
%2 TRIAL_INDEX	
%3 TIMESTAMP
%4 RIGHT_GAZE_X
%5 RIGHT_GAZE_Y
%6 RIGHT_VELOCITY_X
%7 RIGHT_VELOCITY_Y
%8 RESOLUTION_X
%9 RESOLUTION_Y
%10 RIGHT_PUPIL_SIZE
%11 RIGHT_IN_SACCADE
%12 RIGHT_IN_BLINK
%13 SAMPLE_MESSAGE

fid = fopen(txtfile);

if fid == -1
	disp(['could not open file :' txtfile])
	disp(['usual suspect: Path.'])
	return
end

% sk = '%s%u32%u32%u32%f%f%f%f%f%f%d8%d8%s';
sk = '%s%u32%u32%u32%f%f%f%f%f%f%d8%d8%s%f%f';
C = textscan(fid, sk, 'delimiter', '\t','HeaderLines',1,'TreatAsEmpty','.','EmptyValue',NaN);

fclose(fid);


% MSGS = C{13};
% TS = C{3};
% uTS = unique(TS);

% remove offset of time stamps


% M(uTS,1) = [C{2} C{3} C{4} C{5} C{8} C{9}];

% timestamps
TS = unique(C{3}); % TS - timestamp
TS = TS - TS(1)+1;

X = NaN(size(TS));
Y = NaN(size(TS));
X(TS) = double(C{4});
Y(TS) = double(C{5});


RX = NaN(size(TS));
RY = NaN(size(TS));
RX(TS) = double(C{8});
RY(TS) = double(C{9});

pupil = NaN(size(TS));
pupil(TS) = double(C{10});

% resolution (only valid for points close to center of screen!)
resolution(1) = mean(RY(~isnan(RY)));
resolution(2) = mean(RY(~isnan(RY)));

% center of screen
center(1) = 800/2 ;
center(2) = 600/2 ;

target = NaN(length(TS),2);
target(TS,1) = double(C{14});
target(TS,2) = double(C{15});

target(:,1) = (target(:,1) - center(1))/resolution(1);
target(:,2) = (target(:,2) - center(2))/resolution(2);

X = (X - center(1))/resolution(1);
Y = (Y - center(2))/resolution(2);

Y = Y*-1;
target(:,2) = target(:,2)*-1;

find_rwstr = @(s) (strcmp(s, 'SEND_REWARD_TTL'));
find_targetjump = @(s) (strcmp(s, 'DISP_SHIFTED_TARGET'));
find_targstr= @(s) (strcmp(s, 'DISP_TARG'));

successtrials = unique(C{2}(find(cellfun(find_rwstr, C{13}))));


rewardttl_msg = find(cellfun(find_rwstr, C{13}));
rewardttl = zeros(length(C{2}),1);
rewardttl(rewardttl_msg) = 1;

disp_target_msg = find(cellfun(find_targstr, C{13}));
disp_target = zeros(length(C{2}),1);
disp_target(disp_target_msg) = 1;

shift_target_msg = find(cellfun(find_targetjump, C{13}));
shift_target = zeros(length(C{2}),1);
shift_target(shift_target_msg) = 1;



D = [C{2} TS disp_target shift_target rewardttl];


% beginning of valid trials
bvt = @(tr)TS(find(C{2}==tr,1,'first'));
% end of valid trials
evt = @(tr)TS(find(C{2}==tr,1,'last'));
% target on
target_on    = @(tr) D(D(:,1)==tr  & D(:,3)==1,2);
% target jump
target_jump  = @(tr) D(D(:,1)==tr  & D(:,4)==1,2);
% reward ttl
reward_ttl = @(tr) D(D(:,1)==tr  & D(:,5)==1,2);


c= 0;
for t = successtrials'
	
	tb = bvt(t);
	te = evt(t);
	to = target_on(t);
	rw = reward_ttl(t);

	if ~isempty(tb) && ~isempty(te) && ~isempty(to) && ~isempty(rw)
		c = c+1;
		intervals(c,1) = t;
		intervals(c,2) = tb;
		intervals(c,3) = te;
		intervals(c,4) = to(end);
		intervals(c,5) = rw;
	end

	% find trials that start with DISP_TARGET and end with SEND_REWARD_TTL
end




disp(['using ' num2str(c) ' trials of '  num2str(length(unique(C{2}))) ])


out.intervals_fields = {'trial_number' , 'trial_begin' , 'trial_end', 'target_on', 'reward'};
out.intervals = intervals;
out.eyexy = [X Y];
out.successtrials = successtrials;
out.pupilsize = pupil;
out.blinks = find_blinks(pupil);
out.resolution = [RX RY];
out.timestamps = TS;
out.target = target;
out.yield = c / length(unique(C{2}));