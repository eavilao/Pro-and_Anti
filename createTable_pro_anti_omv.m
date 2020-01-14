function selected = createTable_pro_anti_omv(units)
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
indx_sign_instrDir_instr_pro(i) = units(indx_vermis(i)).stats.pro.flags.instrDirVSinstr_nspk;
indx_sign_instrDir_instr_anti(i) = units(indx_vermis(i)).stats.anti.flags.instrDirVSinstr_nspk;
indx_sign_goCue_instrDir_pro(i) = units(indx_vermis(i)).stats.pro.flags.goCueVSinstrDir_nspk;
indx_sign_goCue_instrDir_anti(i) = units(indx_vermis(i)).stats.anti.flags.goCueVSinstrDir_nspk;
indx_sign_saccVSinstrDir_nspk_pro(i) = units(indx_vermis(i)).stats.pro.flags.saccVSinstrDir_nspk;
indx_sign_saccVSinstrDir_nspk_anti(i) = units(indx_vermis(i)).stats.anti.flags.saccVSinstrDir_nspk;
indx_exc_sacc_pro(i) = units(indx_vermis(i)).pro.neural.exc; 
indx_sup_sacc_pro(i) = units(indx_vermis(i)).pro.neural.sup; 
indx_exc_sacc_anti(i) = units(indx_vermis(i)).anti.neural.exc; 
indx_sup_sacc_anti(i) = units(indx_vermis(i)).anti.neural.sup;
indx_exc_instr_pro(i) = units(indx_vermis(i)).pro.neural.exc_instr;
indx_exc_instr_anti(i) = units(indx_vermis(i)).anti.neural.exc_instr;
indx_sup_instr_pro(i) = units(indx_vermis(i)).pro.neural.sup_instr;
indx_sup_instr_anti(i) = units(indx_vermis(i)).anti.neural.sup_instr;
indx_sacc_pro_exc_sup_sign(i) = units(indx_vermis(i)).pro.neural.categ.flag;
indx_sacc_anti_exc_sup_sign(i) = units(indx_vermis(i)).anti.neural.categ.flag;
pval_sacc_pro_exc_sup_sign(i) = units(indx_vermis(i)).pro.neural.categ.pval;
pval_sacc_anti_exc_sup_sign(i) = units(indx_vermis(i)).anti.neural.categ.pval;
indx_instr_pro(i) = units(indx_vermis(i)).stats.pro.flags.instrVSbase_nspk; 
indx_instr_anti(i) = units(indx_vermis(i)).stats.anti.flags.instrVSbase_nspk; 
indx_instr_back_pro(i) = units(indx_vermis(i)).stats.pro.flags.instr_backVSbase_nspk; 
indx_instr_back_anti(i) = units(indx_vermis(i)).stats.anti.flags.instr_backVSbase_nspk; 

num_corr_pro(i) = length(units(i).pro.indx_correctProTrials);  % calculate total number of trials for pro and anti
num_corr_anti(i) = length(units(i).anti.indx_correctAntiTrials);  
end 

pro_only_instr_omv = sum(indx_sign_instr_base_pro & ~indx_sign_instr_base_anti);
anti_only_instr_omv = sum(indx_sign_instr_base_anti & ~indx_sign_instr_base_pro);
pro_anti_omv = sum(indx_sign_instr_base_anti & indx_sign_instr_base_pro);
pro_only_sacc_omv = sum(indx_sign_sacc_base_pro & ~indx_sign_sacc_base_anti);
% pro_only_sacc_exc_omv = sum(indx_exc_sacc_pro & ~indx_sign_sacc_base_anti & indx_sacc_pro_exc_sup_sign);  indx_pro_only_sacc_exc_omv = indx_exc_sacc_pro & ~indx_sign_sacc_base_anti & indx_sacc_pro_exc_sup_sign;
% pro_only_sacc_sup_omv = sum(indx_sup_sacc_pro & ~indx_sign_sacc_base_anti  & indx_sacc_pro_exc_sup_sign);  indx_only_sacc_sup_omv = indx_sup_sacc_pro & ~indx_sign_sacc_base_anti  & indx_sacc_pro_exc_sup_sign;
anti_only_sacc_omv = sum(indx_sign_sacc_base_anti & ~indx_sign_sacc_base_pro);
pro_anti_sacc_omv = sum(indx_sign_sacc_base_anti & indx_sign_sacc_base_pro);
pro_only_instrDir_omv = sum(indx_sign_instrDir_instr_pro & ~indx_sign_instrDir_instr_anti);
anti_only_instrDir_omv = sum(indx_sign_instrDir_instr_anti & ~indx_sign_instrDir_instr_pro);
pro_anti_instrDir_omv = sum(indx_sign_instrDir_instr_pro & indx_sign_instrDir_instr_anti);
pro_only_goCue_omv = sum(indx_sign_goCue_instrDir_pro & ~indx_sign_goCue_instrDir_anti);
anti_only_goCue_omv = sum(indx_sign_goCue_instrDir_anti & ~indx_sign_goCue_instrDir_pro);
pro_anti_goCue_omv = sum(indx_sign_goCue_instrDir_pro & indx_sign_goCue_instrDir_anti);

p_val = 0.05; 
% pro_only_sacc_exc_omv = sum(indx_exc_sacc_pro & ~indx_sign_sacc_base_anti & indx_sacc_pro_exc_sup_sign);  indx_pro_only_sacc_exc_omv = indx_exc_sacc_pro & ~indx_sign_sacc_base_anti & indx_sacc_pro_exc_sup_sign;
% pro_only_sacc_exc_omv =  sum((pval_sacc_pro_exc_sup_sign<p_val) & ~(pval_sacc_anti_exc_sup_sign<p_val) & ~indx_sup_sacc_pro);  indx_pro_only_sacc_exc_omv = (pval_sacc_pro_exc_sup_sign<p_val) & ~(pval_sacc_anti_exc_sup_sign<p_val) & ~indx_sup_sacc_pro;
% pro_only_sacc_sup_omv = sum((pval_sacc_pro_exc_sup_sign<p_val) & ~(pval_sacc_anti_exc_sup_sign<p_val) & ~indx_exc_sacc_pro);  indx_pro_only_sacc_sup_omv = (pval_sacc_pro_exc_sup_sign<p_val) & ~(pval_sacc_anti_exc_sup_sign<p_val) & ~indx_exc_sacc_pro;
% anti_only_sacc_exc_omv = sum(~(pval_sacc_pro_exc_sup_sign<p_val) & (pval_sacc_anti_exc_sup_sign<p_val) & ~indx_sup_sacc_anti);  indx_anti_only_sacc_exc_omv = ~(pval_sacc_pro_exc_sup_sign<p_val) & (pval_sacc_anti_exc_sup_sign<p_val) & ~indx_sup_sacc_anti;
% anti_only_sacc_sup_omv = sum(~(pval_sacc_pro_exc_sup_sign<p_val) & (pval_sacc_anti_exc_sup_sign<p_val) & ~indx_exc_sacc_anti);  indx_anti_only_sacc_sup_omv = ~(pval_sacc_pro_exc_sup_sign<p_val) & (pval_sacc_anti_exc_sup_sign<p_val) & ~indx_exc_sacc_anti;
% pro_anti_sacc_exc_omv = sum(pval_sacc_pro_exc_sup_sign<p_val & (pval_sacc_anti_exc_sup_sign<p_val) & indx_exc_sacc_pro & ~indx_sup_sacc_anti);  indx_pro_anti_sacc_exc_omv = pval_sacc_pro_exc_sup_sign<p_val & (pval_sacc_anti_exc_sup_sign<p_val) & indx_exc_sacc_pro & ~indx_sup_sacc_anti; 
% pro_anti_sacc_sup_omv = sum(pval_sacc_pro_exc_sup_sign<p_val & (pval_sacc_anti_exc_sup_sign<p_val) & ~indx_exc_sacc_pro & indx_sup_sacc_anti);  indx_pro_anti_sacc_sup_omv = pval_sacc_pro_exc_sup_sign<p_val & (pval_sacc_anti_exc_sup_sign<p_val) & ~indx_exc_sacc_pro & indx_sup_sacc_anti;   

pro_only_sacc_exc_omv =  sum((pval_sacc_pro_exc_sup_sign<p_val) & ~indx_sacc_anti_exc_sup_sign & ~indx_sup_sacc_pro);  indx_pro_only_sacc_exc_omv = (pval_sacc_pro_exc_sup_sign<p_val) & ~indx_sacc_anti_exc_sup_sign & ~indx_sup_sacc_pro;
pro_only_sacc_sup_omv = sum((pval_sacc_pro_exc_sup_sign<p_val) & ~indx_sacc_anti_exc_sup_sign & ~indx_exc_sacc_pro);  indx_pro_only_sacc_sup_omv = (pval_sacc_pro_exc_sup_sign<p_val) & ~indx_sacc_anti_exc_sup_sign & ~indx_exc_sacc_pro;
anti_only_sacc_exc_omv = sum(~indx_sacc_pro_exc_sup_sign & (pval_sacc_anti_exc_sup_sign<p_val) & ~indx_sup_sacc_anti);  indx_anti_only_sacc_exc_omv = ~indx_sacc_pro_exc_sup_sign & (pval_sacc_anti_exc_sup_sign<p_val) & ~indx_sup_sacc_anti;
anti_only_sacc_sup_omv = sum(~indx_sacc_pro_exc_sup_sign & (pval_sacc_anti_exc_sup_sign<p_val) & ~indx_exc_sacc_anti);  indx_anti_only_sacc_sup_omv = ~indx_sacc_pro_exc_sup_sign & (pval_sacc_anti_exc_sup_sign<p_val) & ~indx_exc_sacc_anti;
pro_anti_sacc_exc_omv = sum(pval_sacc_pro_exc_sup_sign<p_val & (pval_sacc_anti_exc_sup_sign<p_val) & indx_exc_sacc_pro & ~indx_sup_sacc_anti);  indx_pro_anti_sacc_exc_omv = pval_sacc_pro_exc_sup_sign<p_val & (pval_sacc_anti_exc_sup_sign<p_val) & indx_exc_sacc_pro & ~indx_sup_sacc_anti; 
pro_anti_sacc_sup_omv = sum(pval_sacc_pro_exc_sup_sign<p_val & (pval_sacc_anti_exc_sup_sign<p_val) & ~indx_exc_sacc_pro & indx_sup_sacc_anti);  indx_pro_anti_sacc_sup_omv = pval_sacc_pro_exc_sup_sign<p_val & (pval_sacc_anti_exc_sup_sign<p_val) & ~indx_exc_sacc_pro & indx_sup_sacc_anti;   
both_anti_exc_sup_omv = sum(pval_sacc_pro_exc_sup_sign<p_val & (pval_sacc_anti_exc_sup_sign<p_val) & indx_exc_sacc_pro & indx_sup_sacc_anti);

selected.sacc.pro.exc = indx_vermis(indx_pro_only_sacc_exc_omv); 
selected.sacc.pro.sup = indx_vermis(indx_pro_only_sacc_sup_omv); 
selected.sacc.anti.exc = indx_vermis(indx_anti_only_sacc_exc_omv);
selected.sacc.anti.sup = indx_vermis(indx_anti_only_sacc_sup_omv);
selected.sacc.both.exc = indx_vermis(indx_pro_anti_sacc_exc_omv);
selected.sacc.both.sup = indx_vermis(indx_pro_anti_sacc_sup_omv);
selected.sacc.all.pro.exc = indx_vermis(pval_sacc_pro_exc_sup_sign<p_val & ~indx_sup_sacc_pro); 
selected.sacc.all.anti.exc = indx_vermis(pval_sacc_anti_exc_sup_sign<p_val & ~indx_sup_sacc_anti); 
selected.sacc.all.pro.sup = indx_vermis(pval_sacc_pro_exc_sup_sign<p_val & ~indx_exc_sacc_pro); 
selected.sacc.all.anti.sup = indx_vermis(pval_sacc_anti_exc_sup_sign<p_val & ~indx_exc_sacc_anti);


%% instruction too
instr_pro_only_sacc_exc_omv = sum(indx_pro_only_sacc_exc_omv & indx_instr_back_pro);
instr_pro_only_sacc_sup_omv = sum (indx_pro_only_sacc_sup_omv & indx_instr_back_pro);
instr_anti_only_sacc_exc_omv = sum(indx_anti_only_sacc_exc_omv & indx_instr_back_anti);
instr_anti_only_sacc_sup_omv = sum(indx_anti_only_sacc_sup_omv & indx_instr_back_anti);
instr_pro_anti_sacc_exc_omv = sum(indx_pro_anti_sacc_exc_omv & indx_instr_back_pro & indx_instr_back_anti);
instr_pro_anti_sacc_sup_omv = sum(indx_pro_anti_sacc_sup_omv & indx_instr_back_pro & indx_instr_back_anti);

selected.instr_back.pro.exc = indx_vermis(indx_pro_only_sacc_exc_omv & indx_instr_back_pro); 
selected.instr_back.pro.sup = indx_vermis(indx_pro_only_sacc_sup_omv & indx_instr_back_pro); 
selected.instr_back.anti.exc = indx_vermis(indx_anti_only_sacc_exc_omv & indx_instr_back_anti);
selected.instr_back.anti.sup = indx_vermis(indx_anti_only_sacc_sup_omv & indx_instr_back_anti);
selected.instr_back.both.exc = indx_vermis(indx_pro_anti_sacc_exc_omv & indx_instr_back_pro & indx_instr_back_anti);
selected.instr_back.both.sup = indx_vermis(indx_pro_anti_sacc_sup_omv & indx_instr_back_pro & indx_instr_back_anti);
selected.instr_back.all.pro.exc = indx_vermis(pval_sacc_pro_exc_sup_sign<p_val & ~indx_sup_sacc_pro & indx_instr_back_pro);
selected.instr_back.all.anti.exc = indx_vermis(pval_sacc_anti_exc_sup_sign<p_val & ~indx_sup_sacc_anti & indx_instr_back_anti);
selected.instr_back.all.pro.sup = indx_vermis(pval_sacc_pro_exc_sup_sign<p_val & ~indx_exc_sacc_pro & indx_instr_back_pro); 
selected.instr_back.all.anti.sup = indx_vermis(pval_sacc_anti_exc_sup_sign<p_val & ~indx_exc_sacc_anti & indx_instr_back_anti);


selected.instr.pro.exc = indx_vermis(indx_pro_only_sacc_exc_omv & indx_instr_pro); 
selected.instr.pro.sup = indx_vermis(indx_pro_only_sacc_sup_omv & indx_instr_pro); 
selected.instr.anti.exc = indx_vermis(indx_anti_only_sacc_exc_omv & indx_instr_anti);
selected.instr.anti.sup = indx_vermis(indx_anti_only_sacc_sup_omv & indx_instr_anti);
selected.instr.both.exc = indx_vermis(indx_pro_anti_sacc_exc_omv & indx_instr_pro & indx_instr_anti);
selected.instr.both.sup = indx_vermis(indx_pro_anti_sacc_sup_omv & indx_instr_pro & indx_instr_anti);
selected.instr.all.pro.exc = indx_vermis(pval_sacc_pro_exc_sup_sign<p_val & ~indx_sup_sacc_pro & indx_instr_pro);
selected.instr.all.anti.exc = indx_vermis(pval_sacc_anti_exc_sup_sign<p_val & ~indx_sup_sacc_anti & indx_instr_anti);
selected.instr.all.pro.sup = indx_vermis(pval_sacc_pro_exc_sup_sign<p_val & ~indx_exc_sacc_pro & indx_instr_pro); 
selected.instr.all.anti.sup = indx_vermis(pval_sacc_anti_exc_sup_sign<p_val & ~indx_exc_sacc_anti & indx_instr_anti);


indx_vermis(~indx_pro_only_sacc_exc_omv & ~indx_pro_only_sacc_sup_omv & ~indx_anti_only_sacc_exc_omv & ~indx_anti_only_sacc_sup_omv & ~indx_pro_anti_sacc_exc_omv & ~indx_pro_anti_sacc_sup_omv & indx_instr_pro & indx_instr_anti);
%% 
% pro_sacc_instrDir_omv = sum(indx_sign_saccVSinstrDir_nspk_pro & ~indx_sign_saccVSinstrDir_nspk_anti)
% anti_sacc_instrDir_omv = sum(indx_sign_saccVSinstrDir_nspk_anti & ~indx_sign_saccVSinstrDir_nspk_pro)
% pro_anti_sacc_instrDir_omv = sum(indx_sign_saccVSinstrDir_nspk_pro & indx_sign_saccVSinstrDir_nspk_anti)

vermis_diff_pro_anti_sacc = sum(indx_sign_sacc); %%%%%%%%%%%%% to table
proportion_vermis_sacc = vermis_diff_pro_anti_sacc/length(indx_vermis); %%%%%%%%%%%%% to table
vermis_diff_pro_anti_instr = sum(indx_sign_instr); %%%%%%%%%%%%% to table
proportion_vermis_instr = vermis_diff_pro_anti_instr/length(indx_vermis); %%%%%%%%%%%%% to table
vermis_diff_sacc_base_pro = sum(indx_sign_sacc_base_pro); %%%%%%%%%%%%% to table
vermis_diff_instr_base_pro = sum(indx_sign_instr_base_pro); %%%%%%%%%%%%% to table
vermis_diff_sacc_base_anti = sum(indx_sign_sacc_base_anti); %%%%%%%%%%%%% to table
vermis_diff_instr_base_anti = sum(indx_sign_instr_base_anti); %%%%%%%%%%%%% to table
vermis_diff_instrDir_base_pro = sum(indx_sign_instrDir_instr_pro); %%%%%%%%%%%%% to table
vermis_diff_instrDir_base_anti = sum(indx_sign_instrDir_instr_anti); %%%%%%%%%%%%% to table
vermis_diff_goCue_instrDir_pro = sum(indx_sign_goCue_instrDir_pro); %%%%%%%%%%%%% to table
vermis_diff_goCue_instrDir_anti = sum(indx_sign_goCue_instrDir_anti); %%%%%%%%%%%%% to table
vermis_diff_saccVSinstrDir_nspk_pro = sum(indx_sign_saccVSinstrDir_nspk_pro); %%%%%%%%%%%%% to table
vermis_diff_saccVSinstrDir_nspk_anti = sum(indx_sign_saccVSinstrDir_nspk_anti); %%%%%%%%%%%%% to table

%% instr
indx_vermis_sign = (indx_pro_only_sacc_exc_omv | indx_pro_only_sacc_sup_omv | indx_anti_only_sacc_exc_omv | indx_anti_only_sacc_sup_omv | indx_pro_anti_sacc_exc_omv | indx_pro_anti_sacc_sup_omv); 
instr_pro_pick = sum(indx_vermis_sign & indx_sign_instr_base_pro & ~indx_sign_instr_base_anti); 
instr_pro_exc_pick = sum(indx_vermis_sign & indx_sign_instr_base_pro & ~indx_sign_instr_base_anti & indx_exc_instr_pro);
instr_pro_sup_pick = sum(indx_vermis_sign & indx_sign_instr_base_pro & ~indx_sign_instr_base_anti & indx_sup_instr_pro);
instr_anti_pick = sum(indx_vermis_sign & ~indx_sign_instr_base_pro & indx_sign_instr_base_anti); 
instr_anti_exc_pick = sum(indx_vermis_sign & ~indx_sign_instr_base_pro & indx_sign_instr_base_anti & indx_exc_instr_anti); 
instr_anti_supp_pick = sum(indx_vermis_sign & ~indx_sign_instr_base_pro & indx_sign_instr_base_anti & indx_sup_instr_anti);
instr_pro_anti_pick = sum(indx_vermis_sign & indx_sign_instr_base_pro & indx_sign_instr_base_anti); 
instr_pro_anti_exc_pick = sum(indx_vermis_sign & indx_sign_instr_base_pro & indx_sign_instr_base_anti & indx_exc_instr_pro);
instr_pro_anti_sup_pick = sum(indx_vermis_sign & indx_sign_instr_base_pro & indx_sign_instr_base_anti & indx_sup_instr_pro);
instr_sign_pick = sum(indx_vermis_sign & indx_sign_instr );


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
% lat = units(indx_lateral);
% clear r_pro r_anti
% %get fr for all neurons
% for i = 1:length(lat)
%     r_pro(i) = lat(i).pro.neural.sacc.rate_mu; 
%     r_anti(i) = lat(i).anti.neural.sacc.rate_mu;
% end 
% pro_mu = nanmean(r_pro); pro_std = nanstd(r_pro);
% anti_mu = nanmean(r_anti); anti_std = nanstd(r_anti);

%% Plot
% 
% total_rec = [198 114]; total_rec_labels = {'OMV', 'Lateral'};
% included = [117 46 81 68]; included_labels = {'OMV ns', 'Lat ns', 'OMV sign', 'Lat sign'}; 
% type_activity = [47 34 41 27]; activity_labels = {'omv facil', 'omv supp', 'lat facil', 'lat supp'}; 
% 
% % plot
% p1 = pie(total_rec); p1(2).String = [total_rec_labels{1} p1(2).String]; p1(4).String = [total_rec_labels{2} p1(4).String];
% p2 = pie(included); p2(2).String = [included_labels{1} p2(2).String]; p2(4).String = [included_labels{2} p2(4).String]; p2(6).String = [included_labels{3} p2(6).String]; p2(8).String = [included_labels{4} p2(8).String]; 
% p3 = pie(type_activity); p3(2).String = [activity_labels{1} p3(2).String]; p3(4).String = [activity_labels{2} p3(4).String]; p3(6).String = [activity_labels{3} p3(6).String]; p3(8).String = [activity_labels{4} p3(8).String];
% 
% 


