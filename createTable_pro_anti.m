function createTable_pro_anti(units)
% Create stats table
total_cells = length(units);

%% vermis 
for cellNum = 1:length(units)
    indx_area(cellNum) = strcmp(units(cellNum).area, 'vermis');
end
indx_vermis = find(indx_area);
total_vermis = length(indx_vermis); %%%%%%%%%%%%% to table

% how many significant diff for pro vs anti and base vs windows

for i=1:length(indx_vermis)
indx_sign_sacc(i)= units(indx_vermis(i)).stats.sacc.flags.proVsAnti_sacc_ks_nspk; 
indx_sign_instr(i)= units(indx_vermis(i)).stats.instr.flags.proVsAnti_instr_ks_nspk; 
indx_sign_sacc_base_pro(i)= units(indx_vermis(i)).stats.pro.flags.saccVSbase_nspk; 
indx_sign_instr_base_pro(i)= units(indx_vermis(i)).stats.pro.flags.instrVSbase_nspk;
indx_sign_sacc_base_anti(i)= units(indx_vermis(i)).stats.anti.flags.saccVSbase_nspk; 
indx_sign_instr_base_anti(i)= units(indx_vermis(i)).stats.anti.flags.instrVSbase_nspk;

num_corr_pro(i) = length(units(i).pro.indx_correctProTrials);  % calculate total number of trials for pro and anti
num_corr_anti(i) = length(units(i).anti.indx_correctAntiTrials);  


end 
vermis_diff_pro_anti_sacc = sum(indx_sign_sacc); %%%%%%%%%%%%% to table
proportion_vermis_sacc = vermis_diff_pro_anti_sacc/length(indx_vermis); %%%%%%%%%%%%% to table
vermis_diff_pro_anti_instr = sum(indx_sign_instr); %%%%%%%%%%%%% to table
proportion_vermis_instr = vermis_diff_pro_anti_instr/length(indx_vermis); %%%%%%%%%%%%% to table
vermis_diff_sacc_base_pro = sum(indx_sign_sacc_base_pro); %%%%%%%%%%%%% to table
vermis_diff_instr_base_pro = sum(indx_sign_instr_base_pro); %%%%%%%%%%%%% to table
vermis_diff_sacc_base_anti = sum(indx_sign_sacc_base_anti); %%%%%%%%%%%%% to table
vermis_diff_instr_base_anti = sum(indx_sign_instr_base_anti); %%%%%%%%%%%%% to table




%% lateral
clear indx_area
for cellNum = 1:length(units)
    indx_area(cellNum) = strcmp(units(cellNum).area, 'lateral');
end
indx_lateral = find(indx_area);
total_lateral = length(indx_lateral); %%%%%%%%%%%%% to table

for i=1:length(indx_lateral)
indx_sign_sacc_lat(i)= units(indx_lateral(i)).stats.sacc.flags.proVsAnti_sacc_ks_nspk; 
indx_sign_instr_lat(i)= units(indx_lateral(i)).stats.instr.flags.proVsAnti_instr_ks_nspk; 
indx_sign_sacc_base_pro_lat(i)= units(indx_lateral(i)).stats.pro.flags.saccVSbase_nspk; 
indx_sign_instr_base_pro_lat(i)= units(indx_lateral(i)).stats.pro.flags.instrVSbase_nspk;
indx_sign_sacc_base_anti_lat(i)= units(indx_lateral(i)).stats.anti.flags.saccVSbase_nspk; 
indx_sign_instr_base_anti_lat(i)= units(indx_lateral(i)).stats.anti.flags.instrVSbase_nspk;

num_corr_pro_lat(i) = length(units(i).pro.indx_correctProTrials);  % calculate total number of trials for pro and anti
num_corr_anti_lat(i) = length(units(i).anti.indx_correctAntiTrials);  

end 

lateral_diff_pro_anti_sacc = sum(indx_sign_sacc_lat); %%%%%%%%%%%%% to table
proportion_lateral_sacc = lateral_diff_pro_anti_sacc/length(indx_lateral); %%%%%%%%%%%%% to table
lateral_diff_pro_anti_instr = sum(indx_sign_instr_lat); %%%%%%%%%%%%% to table
proportion_lateral_instr = lateral_diff_pro_anti_instr/length(indx_lateral); %%%%%%%%%%%%% to table
lateral_diff_sacc_base_pro = sum(indx_sign_sacc_base_pro_lat); %%%%%%%%%%%%% to table
lateral_diff_instr_base_pro = sum(indx_sign_instr_base_pro_lat); %%%%%%%%%%%%% to table
lateral_diff_sacc_base_anti = sum(indx_sign_sacc_base_anti_lat); %%%%%%%%%%%%% to table
lateral_diff_instr_base_anti = sum(indx_sign_instr_base_anti_lat); %%%%%%%%%%%%% to table

% print table
disp('>>> VERMIS <<<')
disp(['Total num of cells = ' num2str(total_cells)]); 
disp(['Total num vermis = ' num2str(total_vermis)]); 
disp(['Sign pro vs anti instr = ' num2str(vermis_diff_pro_anti_instr)]);
disp(['Proportion sign pro vs anti instr = ' num2str(proportion_vermis_instr)]);
disp(['Sign pro vs anti sacc = ' num2str(vermis_diff_pro_anti_sacc)]);
disp(['Proportion sign pro vs anti sacc = ' num2str(proportion_vermis_sacc)]);
disp(['Sign diff base vs instr pro = ' num2str(vermis_diff_instr_base_pro)]);
disp(['Sign diff base vs sacc pro = ' num2str(vermis_diff_sacc_base_pro)]);
disp(['Sign diff base vs instr anti = ' num2str(vermis_diff_instr_base_anti)]);
disp(['Sign diff base vs sacc anti = ' num2str(vermis_diff_sacc_base_anti)]);
disp('                             ')
disp('>>> LATERAL <<<')
disp(['Total num of cells = ' num2str(total_cells)]); 
disp(['Total num lateral = ' num2str(total_lateral)]); 
disp(['Sign pro vs anti instr = ' num2str(lateral_diff_pro_anti_instr)]);
disp(['Proportion sign pro vs anti instr = ' num2str(proportion_lateral_instr)]);
disp(['Sign pro vs anti sacc = ' num2str(lateral_diff_pro_anti_sacc)]);
disp(['Proportion sign pro vs anti sacc = ' num2str(proportion_lateral_sacc)]);
disp(['Sign diff base vs instr pro = ' num2str(lateral_diff_instr_base_pro)]);
disp(['Sign diff base vs sacc pro = ' num2str(lateral_diff_sacc_base_pro)]);
disp(['Sign diff base vs instr anti = ' num2str(lateral_diff_instr_base_anti)]);
disp(['Sign diff base vs sacc anti = ' num2str(lateral_diff_sacc_base_anti)]);

%% Firing rate, Amplitude, RT, Peak vel.  -- Nico add CS

%% SACCADE PERIOD
%% Vermis
recArea = 'vermis';
for cellNum = 1:length(units)
    indx_area(cellNum) = strcmp(units(cellNum).area, recArea);
end
indx_area = find(indx_area);

for i = 1:length(indx_area)
    r_pro(i,:) = units(indx_area(i)).pro.neural.sacc.rate_pst;
    sem_pro(i,:) = std(units(indx_area(i)).pro.neural.sacc.rate_pst)/sqrt(length(units(indx_area(i)).pro.neural.trial));
    std_pro(i,:) = std(units(indx_area(i)).pro.neural.sacc.rate_pst);
end

%% Exc and sup
cnt_exc=1; cnt_sup=1;
for cellNum = 1:length(units)
    if strcmp(units(cellNum).area, recArea) && units(cellNum).pro.neural.exc==1
        indx_exc(cnt_exc) = cellNum; cnt_exc=cnt_exc+1;
    elseif strcmp(units(cellNum).area, recArea) && units(cellNum).pro.neural.sup==1
        indx_sup(cnt_sup) = cellNum; cnt_sup=cnt_sup+1;
    end
end

% get exc
t = units(1).pro.neural.sacc.ts_pst_win;
for i = 1:length(indx_exc)
    r_exc_pro(i,:) = units(indx_exc(i)).pro.neural.sacc.rate_pst_win;
    std_exc_pro(i) = std(units(indx_exc(i)).pro.neural.sacc.rate_pst_win);
    sem_exc_pro(i,:)= std(units(indx_exc(i)).pro.neural.sacc.rate_pst_win)/sqrt(length(indx_exc));
    r_exc_anti(i,:) = units(indx_exc(i)).anti.neural.sacc.rate_pst_win;
    std_exc_anti(i,:) = std(units(indx_exc(i)).anti.neural.sacc.rate_pst_win);
    sem_exc_anti(i,:)= std(units(indx_exc(i)).anti.neural.sacc.rate_pst_win)/sqrt(length(indx_sup));
end

% get supp
for i = 1:length(indx_sup)
    r_sup_pro(i,:) = units(indx_sup(i)).pro.neural.sacc.rate_pst_win;
    std_sup_pro(i,:) = std(units(indx_sup(i)).pro.neural.sacc.rate_pst_win);
    sem_sup_pro(i,:)= std(units(indx_sup(i)).pro.neural.sacc.rate_pst_win)/sqrt(length(indx_sup));
    r_sup_anti(i,:) = units(indx_sup(i)).anti.neural.sacc.rate_pst_win;
    std_sup_anti(i,:) = std(units(indx_sup(i)).anti.neural.sacc.rate_pst_win);
    sem_sup_anti(i,:)= std(units(indx_sup(i)).anti.neural.sacc.rate_pst_win)/sqrt(length(indx_sup));
end
%% Lateral 
recArea = 'lateral';


%% Table with FRs and kinematics
%% vermis 
omv = units(indx_vermis);
clear r_pro r_anti
%get fr for all neurons
for i = 1:length(omv)
    r_pro(i) = omv(i).pro.neural.sacc.rate_mu; 
    r_anti(i) = omv(i).anti.neural.sacc.rate_mu;
end 
pro_mu = nanmean(r_pro); pro_std = nanstd(r_pro);
anti_mu = nanmean(r_anti); anti_std = nanstd(r_anti);

%% lateral 
lat = units(indx_lateral);
clear r_pro r_anti
%get fr for all neurons
for i = 1:length(lat)
    r_pro(i) = lat(i).pro.neural.sacc.rate_mu; 
    r_anti(i) = lat(i).anti.neural.sacc.rate_mu;
end 
pro_mu = nanmean(r_pro); pro_std = nanstd(r_pro);
anti_mu = nanmean(r_anti); anti_std = nanstd(r_anti);

%% Plot

total_rec = [198 114]; total_rec_labels = {'OMV', 'Lateral'};
included = [117 46 81 68]; included_labels = {'OMV ns', 'Lat ns', 'OMV sign', 'Lat sign'}; 
type_activity = [47 34 41 27]; activity_labels = {'omv facil', 'omv supp', 'lat facil', 'lat supp'}; 

% plot
p1 = pie(total_rec); p1(2).String = [total_rec_labels{1} p1(2).String]; p1(4).String = [total_rec_labels{2} p1(4).String];
p2 = pie(included); p2(2).String = [included_labels{1} p2(2).String]; p2(4).String = [included_labels{2} p2(4).String]; p2(6).String = [included_labels{3} p2(6).String]; p2(8).String = [included_labels{4} p2(8).String]; 
p3 = pie(type_activity); p3(2).String = [activity_labels{1} p3(2).String]; p3(4).String = [activity_labels{2} p3(4).String]; p3(6).String = [activity_labels{3} p3(6).String]; p3(8).String = [activity_labels{4} p3(8).String];




