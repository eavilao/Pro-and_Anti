% [onsets, rpt, fplot, fanimate] = well_formed_saccade(aev) finds the saccade
% onsets from the eye speed aev. rpt reports 
% '# of detected saccades \t # of well-formed saccades \t acceptance rate (in %)'.
% fplot and fanimate are for cheking the results.
% Written by Sungho Hong, CNS unit, OIST, 2013
function [onsets, rpt, fplot, fanimate] = well_formed_saccade(aev, landmarks)




T = length(aev);

Tpre = 75;
Tpost = 100;
msec = 1e-3;  % 1 ms

%%% lognormal model for eye speed.
lv = log10(aev); 
lv = lv(~isinf(lv));
lv = lv(~isnan(lv));


p005 = 10^(mean(lv) + std(lv)*1.96);
p001 = 10^(mean(lv) + std(lv)*2.576);  %% This seems to work best.
p0001 = 10^(mean(lv) + std(lv)*3.3);


    % M> Using my own landmarks to avoid unpaired pairs
    % 
    % [ups, downs, ~] = find_crossings(aev, p001, Tpre);



    ups = landmarks(:,1);
    downs = landmarks(:,6);

    ups = ups(ups<(T-Tpost))+1; % Shifting makes the alignment easier.
    downs = downs(1:length(ups))+1;
    ib = (ups>Tpre);
    ups = ups(ib);
    downs = downs(ib);






faev = sacfilt(aev);
cc = fsts(aev, ups, Tpre,Tpost);
ccx = fsts(faev, ups, Tpre,Tpost);

%% Find minima in the filtered eye speed within 15 ms, which become onsets.
wsize = 15;wbegin = Tpre-wsize;
[~, ix] = min(ccx(:,wbegin:Tpre)');
onsets0 = ups+ix'+wbegin-Tpre-2;
ccmx = fsts(faev,onsets0,Tpre,Tpost);
ccm = fsts(aev,onsets0,Tpre,Tpost);

%% The quied period = -22ms:-2ms from the onset.
t = max(ccm(:,(Tpre-22):(Tpre-2)),[],2);

% DEPRECATED
% lt = log10(t);
% t_thres = 10^(median(lt)+std(lt)*1.65);
% ic = (t<t_thres);
% DEPRECATED

% DEPRECATED
% ic = cluster(linkage(pdist(lt'),'average'),'maxclust',2);
% if mean(t(ic==1))<mean(t(ic==2))
%     cl = 1;
% else
%     cl = 2;
% end
% ic = (ic==cl);
% DEPRECATED

%%% Find the outliers that the eye moves significantly during the quiet period.
dm = median(t)-min(t);
t_thres = median(t)+dm*2;
ic = (t<t_thres);
onsets = onsets0(ic);

ccmi = fsts(aev,onsets,Tpre,Tpost);

%% Report generation.
rpt = sprintf('%d\t%d\t%.2f', length(ic), sum(ic), 100*sum(ic)/length(ic));

%% Help function to plot all.
    function r = plot_all()
        r = figure;
        
        f1 = subplot(231);
        plot((1:T)*msec,aev, 'k', ups*msec, aev(ups),'or', downs*msec, aev(downs),'og')
        axis tight
        xlabel('Time (s)')
        title(f1, 'Detect saccades via a threshold')
        
        f1 = subplot(232);
        plot((1:T)*msec, aev,'k', (1:T)*msec, faev, 'b', ...
             onsets*msec, faev(onsets),'or', ...
             onsets*msec, aev(onsets),'or') %, downs*msec, faev(downs),'og')
        axis tight
        xlabel('Time (s)')
        title(f1, 'Find transition points')
        
        f1 = subplot(233);
        plot(1:numel(t),sort(t),'o',1:numel(t),ones(size(t))*t_thres,'r')
        axis tight
        xlabel('Saccades')
        ylabel('Max eye speed during the quied period')
        title(f1, 'Prune out the not-well-formed saccades')
        
        f1 = subplot(234);
        plot(-Tpre:Tpost, cc'), axis tight
        xlabel('Time (ms)')
        title(f1, 'Eye speeds for detected saccades')
        
        f1 = subplot(235);
        plot(-Tpre:Tpost, ccm'), axis tight
        title(f1, 'Eye speeds aligned with the detected onsets')
        
        f1 = subplot(236);
        plot(-Tpre:Tpost, ccmi'), axis tight
        title(f1, 'Final result')

    end

fplot = @plot_all;

%% Help function to generate an animation.
    function animate(is_onset, is_ev, ev)
        switch is_onset
            case 'onset'
                q = fsts(ev, onsets, Tpre, Tpost);
            case 'thres'
                q = fsts(ev, ups(ic),Tpre,Tpost);
        end
        
        switch is_ev
            case 'ev'
                ep = q';
            case 'ep'
                ep = (cumsum(q'));
        end
        
        figure;
        plot(ep), axis tight
        a = axis();
        for i=1:(size(ep,1)-1)
            hold off
            plot(ep(i:(i+1),:),'linewidth',5)
            hold on
            plot(ep(1:(i+1),:),':')
            axis(a)
            title(num2str(i))
            pause(0.05)
        end
    end

fanimate = @animate;

end
%%% End of the main function


%% sacfilt cuts out low freq from eye speed
function fsig = sacfilt(sig, varargin)
p = inputParser;
p.addRequired('sig', @isvector);
p.addParamValue('Fs', 1000, @isscalar);
p.parse(sig,varargin{:});
inputs = p.Results;

% Cut out >250 Hz for smoothing. The result is robust with this freq.
cutoff_frequency = 250; 

Fs = inputs.Fs;
if Fs>1000
    order = 8;
else
    order = 16;
end

fsig1 = filt_common(sig, Fs, cutoff_frequency, cutoff_frequency+50, 1, 10); 
fsig = filt_common(fsig1, Fs, 15, 12, 1, 10); % Cut out <15 Hz.

end
