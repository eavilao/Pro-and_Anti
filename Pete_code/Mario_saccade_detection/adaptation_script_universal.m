ope[bla, cn] = system('hostname');

switch strtrim(cn)
    case 'nin297'
        try
            pathtodata = 'C:\Users\flierman\Dropbox\NHP_Data_and_analysis\DATA\Mickey';
            pathtores = 'C:\Users\flierman\Dropbox\NHP_Data_and_analysis\CODE_2_B_deleted';
        catch
            disp('path to data failed')
        end
        
%     case 'GTS-NIN-DESKTOP'
%         pathtodata = ['E:' filesep 'Dropbox' filesep 'NHP_Data_and_analysis' filesep 'DATA' filesep 'Mickey'];
%         pathtores = ['E:' filesep 'Dropbox' filesep 'NHP_Data_and_analysis' filesep 'CODE' filesep 'resources'];
% 
%     case 'gtss-mbp.home'
%             pathtodata = '/Users/gts/Dropbox/NHP_Data_and_analysis/DATA/Mickey';
%             pathtores = '/Users/gts/Dropbox/NHP_Data_and_analysis/CODE/resources';
%     case 'gtss-MacBook-Pro.local'
%             pathtodata = '/Users/gts/Dropbox/NHP_Data_and_analysis/DATA/Mickey';
%             pathtores = '/Users/gts/Dropbox/NHP_Data_and_analysis/CODE/resources';
    otherwise
        cn
end
addpath(genpath(pathtodata));
addpath(genpath(pathtores));

experiments_to_process = [44];

parse_data  = 0;
saccade_detection_and_filter = 0;
plot_summaries = 1;
    save_summary_figures = 0;  % plot_summaries must be run first
plot_experiment_and_trials = 0; %  saccade_detection_and_filter must be run first
plot_adaptation_summary = 0;
plot_selected_trial = 0;
plot_groups = 0;
plot_movement_during_tail = 0;
plot_averaged_profiles= 0;
plot_concatenated_profiles = 0;

txtfile{1} =  ['140409_Out_adap_LE' filesep '140409_Out_adap_LE.txt'];
txtfile{2} =  ['140409_Out_adap2_LE' filesep '140409_Out_adap2_LE.txt'];
txtfile{3} =  ['140410_Out_adap_Ctr_LE' filesep '140410_Out_adap_Ctr_LE.txt'];
txtfile{4} =  ['140410_Out_adap_Ctr2_LE' filesep '140410_Out_adap_Ctr2_LE.txt'];
txtfile{5} =  ['140415_Out_adap_LE' filesep '140415_Out_adap_LE.txt';];
txtfile{6} =  ['140415_Out_adap2_LE' filesep '140415_Out_adap2_LE.txt'];
txtfile{7} =  ['140416_Out_adap_LE' filesep '140416_Out_adap_LE.txt'];
txtfile{8} =  ['140416_Out_adap2_LE' filesep '140416_Out_adap2_LE.txt'];
txtfile{9} =  ['140417_OutRot_adap_LE' filesep '140417_OutRot_adap_LE.txt'];
txtfile{10} = ['140417_OutRot_adap2_LE' filesep '140417_OutRot_adap2_LE.txt'];
txtfile{11} = ['140506_Out_adap_LE' filesep '140506_Out_adap_LE.txt'];
txtfile{12} = ['140506_Out_adap2_LE' filesep '140506_Out_adap2_LE.txt'];
txtfile{13} = ['140507_Out_adap_LE' filesep '140507_Out_adap_LE.txt'];
txtfile{14} = ['140508_Out_adap_LE' filesep '140508_Out_adap_LE.txt'];
txtfile{15} = ['140508_Out_adap2_LE' filesep '140508_Out_adap2_LE.txt'];
txtfile{16} = ['140514_Out_adap_LE' filesep '140514_Out_adap_LE.txt'];
txtfile{17} = ['140515_Out_adap_LE' filesep '140515_Out_adap_LE.txt']; % ephys
txtfile{18} = ['140515_Out_adap2_LE' filesep '140515_Out_adap2_LE.txt']; % ephys
txtfile{19} = ['140528_Out_adap_LE' filesep '140528_Out_adap_LE.txt'];
txtfile{20} = ['140604_Out_adap_RE' filesep '140604_Out_adap_RE.txt'];
txtfile{21} = ['140605_Out_adap_RE' filesep '140605_Out_adap_RE.txt'];
txtfile{22} = ['140610_Out_adap_LE' filesep '140610_Out_adap_LE.txt'];
txtfile{23} = ['140618_switch_adap_LE' filesep '140618_switch_adap_LE.txt'];
txtfile{24} = ['140618_switch_adap2_LE' filesep '140618_switch_adap2_LE.txt'];
txtfile{25} = ['140619_switch_adap_LE' filesep '140619_switch_adap_LE.txt'];
txtfile{26} = ['140619_switch_adap2_LE' filesep '140619_switch_adap2_LE.txt'];
txtfile{27} = ['140624_switch_ctrl_LE' filesep '140624_switch_ctrl_LE.txt'];
txtfile{28} = ['140624_switch_ctrl2_LE' filesep '140624_switch_ctrl2_LE.txt'];
txtfile{29} = ['140625_switch_ctrl_LE' filesep '140625_switch_ctrl_LE.txt'];
txtfile{30} = ['140625_switch_ctrl2_LE' filesep '140625_switch_ctrl2_LE.txt'];
txtfile{31} = ['140625_switch_ctrl3_LE' filesep '140625_switch_ctrl3_LE.txt'];
txtfile{32} = ['140626_switch_ctrl_LE' filesep '140626_switch_ctrl_LE.txt'];
txtfile{33} = ['140626_switch_ctrl2_LE' filesep '140626_switch_ctrl2_LE.txt'];
txtfile{34} = ['140709_switch_LE' filesep '140709_switch_LE.txt'];
txtfile{35} = ['140709_switch2_LE' filesep '140709_switch2_LE.txt'];
txtfile{36} = ['140710_switch_LE' filesep '140710_switch_LE.txt'];
txtfile{37} = ['140710_switch2_LE' filesep '140710_switch2_LE.txt'];
txtfile{38} = ['140717_switch_ctrl_LE' filesep '140717_switch_ctrl_LE.txt'];
txtfile{39} = ['140717_switch_ctrl2_LE' filesep '140717_switch_ctrl2_LE.txt'];
txtfile{40} = ['140724_switch_ctrl_RE' filesep '140724_switch_ctrl_RE.txt'];
txtfile{41} = ['140724_switch_ctrl2_RE' filesep '140724_switch_ctrl2_RE.txt'];
txtfile{42} = ['140731_switch_ctrl_RE' filesep '140731_switch_ctrl_RE.txt'];
txtfile{43} = ['140731_switch_ctrl2_RE' filesep '140731_switch_ctrl2_RE.txt'];
txtfile{44} = ['140812_switch_ctrl_RE' filesep '140812_switch_ctrl_RE.txt'];
txtfile{45} = ['140812_switch_ctrl2_RE' filesep '140812_switch_ctrl2_RE.txt'];
txtfile{46} = ['140814_switch_ctrl_RE' filesep '140814_switch_ctrl_RE.txt'];
txtfile{47} = ['140814_switch_ctrl2_RE' filesep '140814_switch_ctrl2_RE.txt'];


%%

for f = experiments_to_process

    if parse_data
        
        [pth name{f} ext] = fileparts(txtfile{f});
        data{f} = trial_parser(txtfile{f});
        
    end    
        
    
    if saccade_detection_and_filter 

            %  intervals(:,[4 5]) are the timestamps of the target display and the delivery of the reward
       unfilteredsacs{f} = saccadeDetection(data{f}.eyexy,'debug', 0,'detectionintervals', data{f}.intervals(:,[4 5]),'plotresults',plot_summaries);
       filteredsaccades{f} = filterSaccadeStruct(unfilteredsacs{f});
        
    end
    %%
    % intervals(c,1) = t; trial
    % 		intervals(c,2) = tb; begin
    % 		intervals(c,3) = te; end
    % 		intervals(c,4) = to; targed_disp
    % 		intervals(c,5) = rw; rw_ttl

    
    if plot_summaries
       summaryFromSaccadeStruct(filteredsaccades{f},'colorscheme','dark'); 
       if save_summary_figures
           saveas(figure(1),[name{f} '_1']);saveas(figure(2),[name{f} '_2']);saveas(figure(3),[name{f} '_3']); close all;
       end
    end

    %%
    
    if plot_experiment_and_trials

        try
        lc = linspecer(7);
        catch
        lc = jet(7)
        end
         lm = filteredsaccades{f}.saccades.landmarks.v; 
        %lm = unfilteredsaccs{f}.saccades.landmarks.v; 
        % vtraces = unfilteredsacs{f}.traces.velocity.v;
         vtraces = filteredsaccades{f}.traces.velocity.v;
        intervals = data{f}.intervals;


        clf
        ax(1) = subplot(2,1,1);

        hold on
        % plot(filteredsaccades{f}.traces.velocity.v) 
        plot(vtraces) 
        line(lm(:,2), ones(size(lm(:,1)))*-10, 'color', lc(1,:),'linestyle', 'none','markersize', 2,'marker', '.','markersize', 10)
        line(lm(:,3), ones(size(lm(:,1)))*-20, 'color', lc(2,:),'linestyle', 'none','markersize', 2,'marker', '.','markersize', 10)
        line(lm(:,4), ones(size(lm(:,1)))*-30, 'color', lc(3,:),'linestyle', 'none','markersize', 2,'marker', '.','markersize', 10)
        line(lm(:,5), ones(size(lm(:,1)))*-40, 'color', lc(4,:),'linestyle', 'none','markersize', 2,'marker', '.','markersize', 10)
        line(lm(:,6), ones(size(lm(:,1)))*-50, 'color', lc(5,:),'linestyle', 'none','markersize', 2,'marker', '.','markersize', 10)

        % X = intervals(:,2);
        % Y = ones(length(intervals),1)*-59;

        text(double(intervals(:,2)), ones(length(intervals),1)*-70, num2str(intervals(:,1)));
        

        % TRIAL MARKERS
        % intervals 2 - trial on - green
        % intervals 3 - trial off - red
        % intervals 4 - disp target - cyan
        % intervals 5 - shift target - yellow
        % intervals 6 - reward ttl - blue
        
        line(intervals(:,2), ones(length(intervals))*-59, 'color', 'g','linestyle', 'none','markersize', 2,'marker', '+','markersize', 10 )
        line(intervals(:,3), ones(length(intervals))*-59, 'color', 'r','linestyle', 'none','markersize', 2,'marker', '+','markersize', 10 )
        line(intervals(:,4), ones(length(intervals))*-59, 'color', 'c','linestyle', 'none','markersize', 2,'marker', '+','markersize', 10 )
        line(intervals(:,5), ones(length(intervals))*-59, 'color', 'y','linestyle', 'none','markersize', 2,'marker', '+','markersize', 10 )
        % line(intervals(:,6), ones(length(intervals))*-59, 'color', 'b','linestyle', 'none','markersize', 2,'marker', '+','markersize', 10 )

        line(data{f}.blinks.open_intervals', ones(size(data{f}.blinks.open_intervals))'*60,'linewidth', 5,'color','k','linestyle','-')

        psz = data{f}.pupilsize; psz = (psz/100); 
        plot(psz,'r')


        colors = jet(8);

        ax(2) = subplot(2,1,2);
        line([1:length(data{f}.target)], data{f}.target(:,1),'color', colors(1,:),'linewidth',4)
        line([1:length(data{f}.target)], data{f}.target(:,2),'color', colors(7,:),'linewidth',4)

        line([1:length(data{f}.eyexy)], data{f}.eyexy(:,1),'color', colors(2,:),'linewidth',1)
        line([1:length(data{f}.eyexy)], data{f}.eyexy(:,2),'color', colors(8,:),'linewidth',1)
        


        linkaxes(ax, 'x')
    end

    if plot_adaptation_summary
        adaptationSummary(unfilteredsacs{f},data{f})
    end


    if plot_movement_during_tail
        
        lmarks_x = filteredsaccades{f}.saccades.landmarkCoordinates.x;
        lmarks_y = filteredsaccades{f}.saccades.landmarkCoordinates.y;

        % no tails:
        nt = find(lmarks_x(:,5)==0);
        lmarks_x(nt,5) = lmarks_x(nt,4);
        lmarks_x(nt,6) = lmarks_x(nt,4);
        lmarks_y(nt,5) = lmarks_y(nt,4);
        lmarks_y(nt,6) = lmarks_y(nt,4);


        fig{f} = figure;
        subplot(121)

         set(gca,'colororder',linspecer(6))
         line(lmarks_x,lmarks_y,'marker', '.','linestyle', 'none')
         axis equal
         
        subplot(122)
         line(lmarks_x(:,[4 6])',lmarks_y(:,[4 6])','linestyle','-')
         line(lmarks_x(:,[6])',lmarks_y(:,[6])','marker', '.','markersize',15,'color','r','linestyle','none')
            axis equal

         title(txtfile{f},'interpreter','none')

         % set(fig{f},'position',figpos ,'paperunits','centimeters','paperposition',figsize{ttp}, 'visible','off','color','none')

    end



    if plot_averaged_profiles;

        muV = mean(filteredsaccades{f}.saccades.vProfiles.v )';
        stdV = std(filteredsaccades{f}.saccades.vProfiles.v )';

        figure
        clf
        hold on
        plot(filteredsaccades{f}.saccades.vProfiles.v','color', [.7 .7 .7])
        plot(muV,'color','r','linewidth',2)
        plot(muV+stdV,'color', 'r')
        plot(muV-stdV,'color', 'r')

        %plot concatenated sacacdes in a sequence

    end




    if plot_concatenated_profiles
            figure
            plot(reshape(filteredsaccades{f}.saccades.vProfiles.v',1,[]))

    end




    
end


%% SELECT TRIAL AND PLOT SACCADIC TRACE

plot_selected_trial = 0
    if plot_selected_trial
        experiment = 6;
        trial = 52;

        % trial of interest
        toi = find(data{experiment}.intervals(:,1)==52);

        %time stamps of interest
        tsoi = data{experiment}.intervals(toi,:);

        pos_xy = data{experiment}.eyexy;
        % 2 - trial on
        % 3 - trial off
        % 4 - display target
        % 5 - reward
        X = pos_xy(tsoi(:,2):tsoi(:,3),1);
        Y = pos_xy(tsoi(:,2):tsoi(:,3),2);


        figure
        scatter(X , Y,'filled')

    end



if plot_groups

    allsacs = 1:length(filteredsaccades{f}.saccades.landmarks.v);

    groups{1} = 58;
    groups{2} = setdiff(allsacs, groups{1});

    summaryFromSaccadeStruct(filteredsaccades{f},'to_plot', [1 1 1 1 1 1 0 0 0 ], 'groups',groups); 
end


make_movie_of_trial = 0;
if make_movie_of_trial
    % find saccade after display target: index 3


    t78 = data{5}.intervals(find(data{5}.intervals(:,1)==77),[4 5]);
    t80 = data{5}.intervals(find(data{5}.intervals(:,1)==80),[2 3]);
    xy_t = data{5}.eyexy(t78(1):t80(2),:);
    py_t = data{5}.pupilsize(t78(1):t80(2));

end



% Questions:
% - see if glissades are functional
% plot glissade work as a function of trial number







