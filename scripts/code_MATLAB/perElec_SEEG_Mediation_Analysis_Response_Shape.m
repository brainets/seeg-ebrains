%% Global settings
% Analysis settings:
ROI_responsiveness_p_value_threshold = 1; % Threshold to select ROI based on their responsivness (put to 1 to disable)

% Visualisation settings:
NORMALIZU_VISUALISATION= false; % Should the data be centered on their peak point and their starting baseline aligned to zero?
% Plot colors
color_easy       = [0,100,250]/255;
color_hard       = [250,65,0]/255;
color_easyRTlow  = [0,225,250]/255;
color_easyRThigh = [25,0,250]/255;
color_hardRTlow  = [250,130,0]/255;
color_hardRThigh = [250,0,20]/255;

RELOAD_DATA = true;

%% Loading & Preparing variables
if RELOAD_DATA
    load('paramshapeAllElec.mat')
    load('ROI_selection.mat')
    load('subject_XY_&_rmvd_trials.mat')
    load('formated_neural_data_best_elec.mat')
    load('formated_neural_data.mat')
    load('trial_time.mat')
end

[~, n_subj, n_roi] = size(subj_roi_P);

param_name = {'c ','A ','w1','p1','b1','w2','p2','b2'};
min_percentage_before_stim = min(percentage_time);
min_time_before_locking    = min(trial_time);

%% ROI selection (only to restrain even more than in the SEEG_Preproc_and_Shape_Fit.m script)
p_diff_FDR = p_diff*NaN;
ROI_selected = cell(3,1);
for ii_data=1:3
    % FDR correction
    [~,ps_idx]      = sort(p_diff(:,ii_data),'ascend');
    porder          = [];
    porder(ps_idx)  = 1:numel(ps_idx);
    n_non_nan_value = sum(~isnan(p_diff(:,ii_data)));
    p_diff_FDR(:,ii_data) = p_diff(:,ii_data) .* n_non_nan_value ./ porder';
    
    ROI_selected{ii_data} = find(p_diff_FDR(:,ii_data)<ROI_responsiveness_p_value_threshold);
end

ROI_selected = cellfun(@(x) {setdiff(x,83)},ROI_selected); % let's remove the "not in marsatlas ROI"

disp(['Selected ROI (prct time):  ',num2str(numel(all_rois_names(ROI_selected{1}))),'--', strjoin(all_rois_names(ROI_selected{1}),' ; ')])
disp(['Selected ROI (S-lkd time): ',num2str(numel(all_rois_names(ROI_selected{2}))),'--', strjoin(all_rois_names(ROI_selected{2}),' ; ')])
disp(['Selected ROI (R_lkd time): ',num2str(numel(all_rois_names(ROI_selected{3}))),'--', strjoin(all_rois_names(ROI_selected{3}),' ; ')])

%% Mediation analyses of response shape
[Shape_a,Shape_b] = deal(cell(3,n_subj,n_roi));
[pa_roiXparam,pb_roiXparam] = deal(NaN(3,n_roi,8));
for ii_data = 1:3
    %% med analysis estimation
    % For each subjects, for each selected ROI, for each shape parameter,
    % we will assess the effect of the trial type on the trial-per-trial
    % variation of the shape parameter (coefficient a) and then assess how
    % much this parameter is related (coefficient b) to the RT variation
    % (above and beyond any trial-type effect)
    for ii_subj=1:n_subj
        [X,Y,bad_idx]=deal(XYs_and_bad_index{ii_subj,:});
        %FYI: X = [cst, zscore([trial_number, trial_type])]; Y = zscore(RT); bad_idx = RT>2.5
        
        get_Xbeta_mat = pinv(X'*X)*X';
        for ii_roi =1:n_roi
            if ~isempty(subj_roi_P{ii_data,ii_subj,ii_roi})
                [Shape_a{ii_data,ii_subj,ii_roi},Shape_b{ii_data,ii_subj,ii_roi}] = deal(NaN(...
                    size(subj_roi_P{ii_data,ii_subj,ii_roi},1),size(subj_roi_P{ii_data,ii_subj,ii_roi}{1},1)));
                for ii_elec = 1:size(subj_roi_P{ii_data,ii_subj,ii_roi},1)
                    Ps= cat(2,subj_roi_P{ii_data,ii_subj,ii_roi}{ii_elec,:}); %subj_roi_P{ii_data,ii_subj,ii_roi}{ii_elec,ii_trial}
                    
                    for ii_param = 1:size(Ps,1)
                        % Extract the ii_param-th parameter for each trial
                        M = Ps(ii_param,:)';
                        
                        M(bad_idx)=[];
                        
                        X_beta   = get_Xbeta_mat*M; % FYI: X_beta(1): cst, X_beta(2) : effect of trial number, X_beta(3) effect of trial type
                        
                        XM = [X,zscore(M-X*X_beta)]; % removing all effect of X on sEEG data (mean, trial number & trial type)
                        XM_betas =  pinv(XM'*XM)*XM'*Y;
                        
                        Shape_a{ii_data,ii_subj,ii_roi}(ii_elec,ii_param)= X_beta(end);
                        Shape_b{ii_data,ii_subj,ii_roi}(ii_elec,ii_param)= XM_betas(end);
                    end
                end
                Shape_a{ii_data,ii_subj,ii_roi} = mean(Shape_a{ii_data,ii_subj,ii_roi},1)';
                Shape_b{ii_data,ii_subj,ii_roi} = mean(Shape_b{ii_data,ii_subj,ii_roi},1)';
            end
        end
    end
    
    %% med analysis uncorrected stat evaluation
    % For each ROI and each parameter we gather all the 'a' and 'b'
    % coefficients across the participant that where recorded in this area
    % and run a simple unsigned T-test on them.
    for ii_roi =1:n_roi
        roi_as = cat(2,Shape_a{ii_data,:,ii_roi});
        roi_bs = cat(2,Shape_b{ii_data,:,ii_roi});
        for ii_param = 1:size(roi_as,1)
            [~,pa_roiXparam(ii_data,ii_roi,ii_param)] = ttest(roi_as(ii_param,:),0,'Tail','Both');
            [~,pb_roiXparam(ii_data,ii_roi,ii_param)] = ttest(roi_bs(ii_param,:),0,'Tail','Both');
        end
    end
end

%% Multiple comparison correction
% For a mediation to be detected both 'a' and 'b' coefficient should be
% significant. But the more ROIs are tested the stronger the multiple
% comparison correction needs to be. Therefore we here use the test on one
% coefficient to constrain the region on which to test for the second.
% Because the less likely effect is an sEEG effect on the RT (above and
% beyond the task-effect) we will first use the 'b' coefficient to select a
% first set of statistically significant ROI, and then test them for the
% 'a' coefficient.
[ab_selected_rois] = deal(cell(3,1));
[selected_idx_a, selected_idx_b,stat_report] = deal(cell(3,n_roi));
[pa_roiXparam_FDR,pb_roiXparam_FDR,pa_roiXparam_FDRbFiltered]=deal(NaN(size(pa_roiXparam)));
for ii_data=1:3
    [all_pbs_FDR,all_pas_FDR,all_pas_bFDR_filtered_FDR] = deal(nan(n_roi,8));
    for ii_param=1:8
        all_pbs         = squeeze(pb_roiXparam(ii_data,:,ii_param));
        all_pas         = squeeze(pa_roiXparam(ii_data,:,ii_param));
        
        % FDR correction on b
        pbs_filtered    = all_pbs;
        pbs_filtered(setdiff(1:n_roi,ROI_selected{ii_data})) = NaN;
        [~,ps_idx]      = sort(pbs_filtered(:),'ascend');
        porder          = [];
        porder(ps_idx)  = 1:numel(ps_idx);
        n_non_nan_value = sum(~isnan(pbs_filtered(:)));
        all_pbs_FDR(:,ii_param)     = pbs_filtered .* n_non_nan_value ./ reshape(porder,size(pbs_filtered));
        
        % Detection of ROI covarying with RT (above and beyond the
        % tasks-effect):
        b_selected_rois  = find((all_pbs_FDR(:,ii_param)<0.05));
        
        % FDR correction on "a" on B selected ROIs
        pas_filtered    = all_pas;
        pas_filtered(~(all_pbs_FDR(:,ii_param)<0.05)) = NaN;
        [~,ps_idx]      = sort(pas_filtered(:),'ascend');
        porder          = [];
        porder(ps_idx)  = 1:numel(ps_idx);
        n_non_nan_value = sum(~isnan(pas_filtered(:)));
        all_pas_bFDR_filtered_FDR(:,ii_param)     = pas_filtered .* n_non_nan_value ./ reshape(porder,size(pas_filtered));
        
        % FDR correction on "a" on all ROIs
        pas_Non_filtered    = all_pas;
        [~,ps_idx]      = sort(pas_Non_filtered(:),'ascend');
        porder          = [];
        porder(ps_idx)  = 1:numel(ps_idx);
        n_non_nan_value = sum(~isnan(pas_Non_filtered(:)));
        all_pas_FDR(:,ii_param)     = pas_Non_filtered .* n_non_nan_value ./ reshape(porder,size(pas_Non_filtered));
        
        % results storing
        pa_roiXparam_FDRbFiltered(ii_data,:,ii_param) = all_pas_bFDR_filtered_FDR(:,ii_param);
        pa_roiXparam_FDR(ii_data,:,ii_param) = all_pas_FDR(:,ii_param);
        pb_roiXparam_FDR(ii_data,:,ii_param) = all_pbs_FDR(:,ii_param);
    end
    
    % Detection of ROI :
    ab_selected_rois{ii_data}  = find(any(all_pas_bFDR_filtered_FDR<0.05,2));
    
    % text display of the results
    if ~isempty(ab_selected_rois{ii_data})
        for selected_roi = ab_selected_rois{ii_data}'
            % Collect coefficients for effect-size comparison
            [~,selected_idx_a{ii_data,selected_roi}] = min(all_pas_bFDR_filtered_FDR(selected_roi,:));
            [~,selected_idx_b{ii_data,selected_roi}] = min(all_pbs_FDR(selected_roi,:));
            
            roi_as = cat(2,Shape_a{ii_data,:,selected_roi});
            roi_bs = cat(2,Shape_b{ii_data,:,selected_roi});
            
            switch ii_data
                case 1, time_name = 'Percentage time';
                case 2, time_name = 'S-locked time';
                case 3, time_name = 'R-locked time';
            end
            roi_txt = {['Selected ROI (',time_name,') ', strrep(all_rois_names{selected_roi},'_',' '),': ']};
            
            disp('===================')
            disp(roi_txt{1})
            
            n_roi_subj = size(roi_bs,2);
            for ii_param = 1:8
                a_mean      = mean(roi_as(ii_param,:));
                a_sem       = std(roi_as(ii_param,:))/sqrt(n_roi_subj);
                
                is_a_sig = all_pas_bFDR_filtered_FDR(selected_roi,ii_param)<0.05;
                is_b_sig = all_pbs_FDR(selected_roi,ii_param)<0.05;
                
                [emphasize_a_start,emphasize_a_stop,emphasize_b_start,emphasize_b_stop] = deal('');
                if is_a_sig
                    emphasize_a_start = '<strong>';
                    emphasize_a_stop = '</strong>';
                end
                if is_b_sig
                    emphasize_b_start = '<strong>';
                    emphasize_b_stop = '</strong>';
                end
                
                a_part_txt = [emphasize_a_start,'X->M Shape-param ', param_name{ii_param},': ', ...
                    num2str(a_mean,2),'+-',num2str(a_sem,2), ...
                    ' (p-FDR:',num2str(all_pas_bFDR_filtered_FDR(selected_roi,ii_param),2),...
                    ', p-uncor:',num2str(pa_roiXparam(ii_data,selected_roi,ii_param),2), ...
                    ', n:',num2str(n_roi_subj),')',emphasize_a_stop];
                b_mean     = mean(roi_bs(ii_param,:));
                b_sem      = std(roi_bs(ii_param,:))/sqrt(n_roi_subj);
                
                b_part_txt = [emphasize_b_start,'M->Y|X : ', ...
                    num2str(b_mean,2),'+-',num2str(b_sem,2), ...
                    ' (p-FDR:',num2str(all_pbs_FDR(selected_roi,ii_param),2),...
                    ', p-uncor:',num2str(pb_roiXparam(ii_data,selected_roi,ii_param),2), ...
                    ', n:',num2str(n_roi_subj),')',emphasize_b_stop];
                
                disp([a_part_txt,'   ', b_part_txt])
                roi_txt{end+1} = [a_part_txt,'   ', b_part_txt];
            end
            stat_report{ii_data,selected_roi} = roi_txt;
        end
    end
end

Table_ComplTime_pValues_path_Diffi_to_ShapeParamSEEG     = array2table(squeeze(pa_roiXparam(1,:,:)),'RowNames',all_rois_names,'VariableNames',param_name);
Table_ComplTime_pValues_path_ShapeParamSEEG_to_RT        = array2table(squeeze(pb_roiXparam(1,:,:)),'RowNames',all_rois_names,'VariableNames',param_name);
Table_ComplTime_pValuesFDR_path_ShapeParamSEEG_to_RT            = array2table(squeeze(pb_roiXparam_FDR(1,:,:)),'RowNames',all_rois_names,'VariableNames',param_name);
Table_ComplTime_pValuesFDR_path_Diffi_to_ShapeParamSEEG         = array2table(squeeze(pa_roiXparam_FDR(1,:,:)),'RowNames',all_rois_names,'VariableNames',param_name);
Table_ComplTime_pValuesFDRfiltered_path_Diffi_to_ShapeParamSEEG = array2table(squeeze(pa_roiXparam_FDRbFiltered(1,:,:)),'RowNames',all_rois_names,'VariableNames',param_name);

Table_SLockTime_pValues_path_Diffi_to_ShapeParamSEEG     = array2table(squeeze(pa_roiXparam(2,:,:)),'RowNames',all_rois_names,'VariableNames',param_name);
Table_SLockTime_pValues_path_ShapeParamSEEG_to_RT        = array2table(squeeze(pb_roiXparam(2,:,:)),'RowNames',all_rois_names,'VariableNames',param_name);
Table_SLockTime_pValuesFDR_path_ShapeParamSEEG_to_RT            = array2table(squeeze(pb_roiXparam_FDR(2,:,:)),'RowNames',all_rois_names,'VariableNames',param_name);
Table_SLockTime_pValuesFDR_path_Diffi_to_ShapeParamSEEG         = array2table(squeeze(pa_roiXparam_FDR(2,:,:)),'RowNames',all_rois_names,'VariableNames',param_name);
Table_SLockTime_pValuesFDRfiltered_path_Diffi_to_ShapeParamSEEG = array2table(squeeze(pa_roiXparam_FDRbFiltered(2,:,:)),'RowNames',all_rois_names,'VariableNames',param_name);

Table_RLockTime_pValues_path_Diffi_to_ShapeParamSEEG     = array2table(squeeze(pa_roiXparam(3,:,:)),'RowNames',all_rois_names,'VariableNames',param_name);
Table_RLockTime_pValues_path_ShapeParamSEEG_to_RT        = array2table(squeeze(pb_roiXparam(3,:,:)),'RowNames',all_rois_names,'VariableNames',param_name);
Table_RLockTime_pValuesFDR_path_ShapeParamSEEG_to_RT            = array2table(squeeze(pb_roiXparam_FDR(3,:,:)),'RowNames',all_rois_names,'VariableNames',param_name);
Table_RLockTime_pValuesFDR_path_Diffi_to_ShapeParamSEEG         = array2table(squeeze(pa_roiXparam_FDR(3,:,:)),'RowNames',all_rois_names,'VariableNames',param_name);
Table_RLockTime_pValuesFDRfiltered_path_Diffi_to_ShapeParamSEEG = array2table(squeeze(pa_roiXparam_FDRbFiltered(3,:,:)),'RowNames',all_rois_names,'VariableNames',param_name);

Shape_parameter_name_and_order = param_name;
Shape_parameter_description    = {'c: Peak time','A: Peak amplitude','w1: Integration window size','p1: Integration concavity','b1: Starting baseline','w2: Depletion window size','p2: Depletion concavity','b2: Finishing baseline'};


Table_ComplTime_pValues_path_Diffi_to_ShapeParamSEEG.ROIs            = all_rois_names';
Table_ComplTime_pValues_path_ShapeParamSEEG_to_RT.ROIs               = all_rois_names';
Table_ComplTime_pValuesFDR_path_ShapeParamSEEG_to_RT.ROIs            = all_rois_names';
Table_ComplTime_pValuesFDR_path_Diffi_to_ShapeParamSEEG.ROIs         = all_rois_names';
Table_ComplTime_pValuesFDRfiltered_path_Diffi_to_ShapeParamSEEG.ROIs = all_rois_names';


% Uncomment the following commands to save the data
%     save('ShapeMediationPvalues_Subj-by-ROI-by-Param.mat',...
%         'Table_ComplTime_pValues_path_Diffi_to_ShapeParamSEEG', 'Table_SLockTime_pValues_path_Diffi_to_ShapeParamSEEG', 'Table_RLockTime_pValues_path_Diffi_to_ShapeParamSEEG',...
%         'Table_ComplTime_pValues_path_ShapeParamSEEG_to_RT'   , 'Table_SLockTime_pValues_path_ShapeParamSEEG_to_RT',    'Table_RLockTime_pValues_path_ShapeParamSEEG_to_RT',...
%         'Table_ComplTime_pValuesFDR_path_ShapeParamSEEG_to_RT', 'Table_SLockTime_pValuesFDR_path_ShapeParamSEEG_to_RT', 'Table_RLockTime_pValuesFDR_path_ShapeParamSEEG_to_RT',...
%         'Table_ComplTime_pValuesFDR_path_Diffi_to_ShapeParamSEEG', 'Table_SLockTime_pValuesFDR_path_Diffi_to_ShapeParamSEEG', 'Table_RLockTime_pValuesFDR_path_Diffi_to_ShapeParamSEEG',...
%         'Table_ComplTime_pValuesFDRfiltered_path_Diffi_to_ShapeParamSEEG', 'Table_SLockTime_pValuesFDRfiltered_path_Diffi_to_ShapeParamSEEG', 'Table_RLockTime_pValuesFDRfiltered_path_Diffi_to_ShapeParamSEEG',...
%         'Shape_parameter_name_and_order','Shape_parameter_description');
%
% writetable(Table_ComplTime_pValues_path_Diffi_to_ShapeParamSEEG,'Table_ComplTime_pValues_path_Diffi_to_ShapeParamSEEG.csv')
% writetable(Table_ComplTime_pValues_path_ShapeParamSEEG_to_RT ,'Table_ComplTime_pValues_path_ShapeParamSEEG_to_RT.csv')
% writetable(Table_ComplTime_pValuesFDR_path_ShapeParamSEEG_to_RT ,'Table_ComplTime_pValuesFDR_path_ShapeParamSEEG_to_RT.csv')
% writetable(Table_ComplTime_pValuesFDR_path_Diffi_to_ShapeParamSEEG ,'Table_ComplTime_pValuesFDR_path_Diffi_to_ShapeParamSEEG.csv')
% writetable(Table_ComplTime_pValuesFDRfiltered_path_Diffi_to_ShapeParamSEEG ,'Table_ComplTime_pValuesFDRfiltered_path_Diffi_to_ShapeParamSEEG.csv')
%
% disp("P values saved.")

%% Show result in the data
% Here we visualize the effect directly on the data. It produces the graph:
% |             Stat report                | sEEG deviation RT< vs RT> in easy trials          |
% |        sEEG easy vs hard trial         | sEEG deviation RT< vs RT> in hard trials          |
% |            Mean model R2s              |                                                   |
% | Relative parameters easy vs hard trial | Relative parameters easy vs % hard vs RT< vs RT > |
% Caption: The above graph visualise the effect of the mediation. Plain
% lines and shaded areas correspond to sEEG mean signal and sem across the group,
% dotted lines to mean model fit. On the top left corner are the mediation
% statistics. The middle left panel show the mean across easy (blue) and
% hard (red) trials. The right pannels show the signal deviation depending
% on their RT, and compared to each subject's mean signal in easy
% (top) or hard (bottom) trials. For each trial types separately, RT were
% median-splitted in two category "RT low" and "RT high" for visualisation.
% The bottom pannels show the shape parameters, z-scored regardless of
% conditions and separated either by trial-types (left) or by RT as well (right).

% To explore the data, regardless of the mediation analysis results, you
% could just fill the variable "ROI_to_visualize" with a name from the
% array 'all_rois_names' and relaunch this section (CTRL+Enter).
% For instance                ROI_to_visualize = 'R_MTCc';
% Otherwise let it empty:     ROI_to_visualize = [];
ROI_to_visualize = [];
if ~isempty(ROI_to_visualize), [~,ROI_to_visualize] = ismember(ROI_to_visualize, all_rois_names); end

REORIENT_AND_REMOVE_START_BASELINE_WITHOUT_CENTERING = false;
% Set it to true and NORMALIZU_VISUALISATION to false, to compare the effect of centering on the signal's peak

all_neural_data = {all_neural_data_percent_time, all_neural_data_stim_lock, all_neural_data_resp_lock};


% Uncomment the following two lines to force the visualisation of the
% completion-time's ROIs onto S-locked and R-locked time.
ab_selected_rois{2}=ab_selected_rois{1};
ab_selected_rois{3}=ab_selected_rois{1};


% Uncomment the following three lines to force the visualisation of all
% regions
% ab_selected_rois{1}= (1:83)'
% ab_selected_rois{2}= (1:83)'
% ab_selected_rois{3}= (1:83)'

removeXeffect = 0; % Put to "1" to visualise the second part of the mediation effect directly on the data
[r2,r2_easy,r2_hard] = deal({});
for ii_data=1:3
    switch ii_data
        case 1; prefix = 'P-time';
        case 2; prefix = 'S-time';
        case 3; prefix = 'R-time';
    end
    figsavefolder = [prefix, 'All_ROI_SignalCondX'];
    mkdir(figsavefolder)
    roi_list = ab_selected_rois{ii_data};
    if ~isempty(ROI_to_visualize), roi_list = ROI_to_visualize; end
    if ~isempty(roi_list)
        for selected_roi = roi_list'
            % Collect neural and fitted data
            %             all_subj_neur = all_neural_data_best_elec{ii_data}(:,selected_roi);
            %             all_subj_simu = subj_roi_simu(ii_data,:,selected_roi);
            
            % Prepare all variables to be visualised
            [   all_simu_mean_easy,         all_simu_mean_hard,...
                all_neur_mean_easy,         all_neur_mean_hard,...
                all_P_mean, all_w2RT_scaled_easy, all_w2RT_scaled_hard,...
                all_w2RT_scaled_RT_low_easy, all_w2RT_scaled_RT_low_hard, ...
                all_w2RT_scaled_RT_hig_easy, all_w2RT_scaled_RT_hig_hard, ...
                all_RT_and_P_easy, all_RT_and_P_hard, ...
                all_P_mean_easy,            all_P_mean_hard,...
                all_R2_mean_easy,           all_R2_mean_hard, all_R2_means, ...
                all_simu_mean_dev_easy_RT_low,  all_simu_mean_dev_hard_RT_low,...
                all_simu_mean_dev_easy_RT_high, all_simu_mean_dev_hard_RT_high,...
                all_neur_mean_dev_easy_RT_low,  all_neur_mean_dev_hard_RT_low,...
                all_neur_mean_dev_easy_RT_high, all_neur_mean_dev_hard_RT_high,...
                all_P_mean_easy_RT_low,     all_P_mean_hard_RT_low,...
                all_P_mean_easy_RT_high,    all_P_mean_hard_RT_high,...
                all_R2_mean_easy_RT_low,     all_R2_mean_hard_RT_low,...
                all_R2_mean_easy_RT_high,    all_R2_mean_hard_RT_high] = deal(cell(n_subj,1));
            for ii_subj=1:n_subj
                
                if ~isempty(subj_roi_simu{ii_data,ii_subj,selected_roi})
                    [simus_per_elec, neur_per_elec, Ps_per_elec, kept_RT_per_elec] = deal(cell(size(subj_roi_simu{ii_data,ii_subj,selected_roi},1),1));
                    for ii_elec=1:size(subj_roi_simu{ii_data,ii_subj,selected_roi},1)
                        
                        %FYI:  X = [cst, zscore([trial_number, trial_type])]; Y = zscore(RT); bad_idx = RT>2.5
                        [X,RT,bad_idx]=deal(XYs_and_bad_index{ii_subj,:});
                        
                        % Prepare all the indexes to separate each condition
                        % and behaviors
                        trial_type  = X(:,3)>0; % >0 trial are easy-trial
                        RT_low_easy = RT <= median(RT(trial_type));   % to be used with 'trial_type'
                        RT_low_hard = RT <= median(RT(~trial_type));  % to be used with '~trial_type'
                        
                        % Collect subject data
                        neur             = squeeze(all_neural_data{ii_data}{ii_subj,selected_roi}(:,ii_elec,:))';
                        neur(bad_idx,:)  = [];
                        
                        simus = cat(1,subj_roi_simu{ii_data,ii_subj,selected_roi}{ii_elec,:});
                        simus(bad_idx,:) = [];
                        
                        Ps = cat(2,subj_roi_P{ii_data,ii_subj,selected_roi}{ii_elec,:})';
                        Ps(bad_idx,:)    = [];
                        
                        if ~NORMALIZU_VISUALISATION && REORIENT_AND_REMOVE_START_BASELINE_WITHOUT_CENTERING
                            simus = (simus-Ps(:,5));
                            neur  = (neur-Ps(:,5));
                        end
                        
                        nancorr=@(x,y) corr(x(~isnan(x))',y(~isnan(y))'); % handy shortcut
                        if NORMALIZU_VISUALISATION
                            switch ii_data
                                % Detect the index of the peak
                                case 1 % Completion time
                                    center_idxs = round(Ps(:,1)*100 - min_percentage_before_stim*1 + 1);
                                case {2,3} % Stim- or Resp-locked
                                    center_idxs = round(Ps(:,1)*100 - round(min_time_before_locking/diff(trial_time(1:2))) + 1);
                            end
                            for ii_t = 1:size(simus,1)
                                simu_ori       = simus(ii_t,:);
                                neur_ori       = neur(ii_t,:);
                                
                                % Adjustment of the reference point
                                % 1) correct for out-of-time-window overflow
                                reference_point = max(1,min(center_idxs(ii_t),size(simus,2)));
                                % 2) correct for the peak not beeing precisely
                                % on the rounded time point:
                                % by 2.1) getting the candidate indexes:
                                search_idx = unique( max(1,min( reference_point + (-2:2), size(simus,2)) ) );
                                % % then 2.2) detecting the maximal
                                % absolute value (after correcting for
                                % the starting baseline (see shape model))
                                [~, idx_extrem] = nanmax( abs( simus(ii_t, search_idx) - Ps(ii_t,5) ) );
                                reference_point = search_idx(idx_extrem);
                                
                                % Now we need to detect which part of
                                % the fitted shape response can be
                                % projected around the center of the
                                % visualisation window.
                                % If the reference_point is below the
                                % middle, than we can take at most 'reference_point'
                                % data points before the peak, otherwise
                                % the maximal number of point is half the window's size.
                                left_window_size  = min( reference_point,          round( numel(simu_ori)/2)      ) - 1;
                                % The same logic apply in reverse for the window after the peak.
                                right_window_size = min( round(numel(simu_ori)/2), numel(simu_ori)-reference_point) - 1;
                                if right_window_size==-1, right_window_size=0; end
                                
                                % Creating the index that will be taken
                                % around the peak-time and placed at
                                % the center of the visualized window
                                relative_idx      = -left_window_size:right_window_size;
                                
                                % Now replace neural and fitted data
                                % with their peak centered version
                                simus(ii_t,:)  = simus(ii_t,:)*nan; % put nan everywhere else
                                neur(ii_t,:)   = neur(ii_t,:) *nan;
                                
                                simus(ii_t,round(end/2)+relative_idx) = simu_ori(reference_point + relative_idx);
                                neur( ii_t,round(end/2)+relative_idx) = neur_ori(reference_point + relative_idx);
                                
                                % Normalize with the starting baseline
                                simus(ii_t,:) = simus(ii_t,:) - Ps(ii_t,5);
                                neur(ii_t,:)  = neur(ii_t,:)  - Ps(ii_t,5);
                                
                                % Re-orient the signal according to the
                                % peak sign (relative to the baseline)
                                peak_sign = sign(Ps(ii_t,2) - Ps(ii_t,5));
                                
                                simus(ii_t,:) = simus(ii_t,:)/peak_sign;
                                neur(ii_t,:)  = neur(ii_t,:)/peak_sign;
                                
                            end
                            
                            
                        end % End of normalisation
                        % store data per elec
                        simus_per_elec{ii_elec} = simus;
                        neur_per_elec{ii_elec}  = neur;
                        Ps_per_elec{ii_elec}    = Ps;
                        kept_RT_per_elec{ii_elec}        = RT(~bad_idx);
                        
                    end
                    % average across elec
                    
                    simus     = cat(3,simus_per_elec{:});
                    neur      = cat(3,neur_per_elec{:});
                    Ps        = cat(3,Ps_per_elec{:});
                    kept_RT   = cat(3,kept_RT_per_elec{:});
                    
                    % Some rare trials will be poorly fitted
                    % and result in bad normalisation
                    % (typically when the signal looks like
                    % pure noise).
                    % While this is not a problem for the
                    % statistics, they need to be remove from
                    % the visualisation.
                    bad_norma = find( any( nanstd(simus,0,2)*25 < nanstd(neur,0,2) ,3));
                    % If the normalised signal varies 25 times more strongly than the fitted signal, there is a problem
                    if any(bad_norma)
                        simus(bad_norma,:,:)     = [];
                        neur(bad_norma,:,:)      = [];
                        Ps(bad_norma,:,:)        = [];
                        trial_type(bad_norma)    = [];
                        RT_low_easy(bad_norma)   = [];
                        RT_low_hard(bad_norma)   = [];
                        kept_RT(bad_norma,:)     = [];
                    end
                    
                    % average across electrodes
                    simus =  squeeze(nanmean(simus,3));
                    neur  =  squeeze(nanmean(neur,3));
                    Ps    =  squeeze(nanmean(Ps,3));
                    kept_RT =  squeeze(nanmean(kept_RT,2));
                    
                    %% Collecting every variable whose mean will be visualised
                    unz_PS = squeeze(nanmean(Ps,3));
                    Ps = zscore(Ps);
                    
                    % handi shortcuts
                    mean_fun = @mean;
                    if NORMALIZU_VISUALISATION, mean_fun = @nanmean; end
                    get_mean_R2 = @(input1,input2,idxs) mean_fun( arrayfun(@(idx) nancorr( input1( idx,:), input2( idx,:) ),find(idxs)).^2 );
                    
                    easy_idx        =  trial_type;
                    hard_idx        = ~trial_type;
                    easy_RTlow_idx  =  trial_type & RT_low_easy;
                    easy_RThigh_idx =  trial_type & ~RT_low_easy;
                    hard_RTlow_idx  = ~trial_type & RT_low_hard;
                    hard_RThigh_idx = ~trial_type & ~RT_low_hard;
                    
                    % X -> M effects
                    % Neural data
                    all_neur_mean_easy{ii_subj} = mean_fun(neur( easy_idx,:),1);
                    all_neur_mean_hard{ii_subj} = mean_fun(neur( hard_idx,:),1);
                    
                    % Fitted data
                    all_simu_mean_easy{ii_subj} = mean_fun(simus( easy_idx,:),1);
                    all_simu_mean_hard{ii_subj} = mean_fun(simus( hard_idx,:),1);
                    
                    % Shape parameters
                    all_P_mean_easy{ii_subj}    = mean(Ps( easy_idx,:),1);
                    all_P_mean_hard{ii_subj}    = mean(Ps( hard_idx,:),1);
                    
                    % suplemental shape parameters info
                    all_P_mean{ii_subj}         = mean(unz_PS,1);
                    if ii_data==1
                        all_w2RT_scaled_easy{ii_subj}    = mean(posi_fun(unz_PS( easy_idx,6)).*kept_RT(easy_idx),1);
                        all_w2RT_scaled_hard{ii_subj}    = mean(posi_fun(unz_PS( hard_idx,6)).*kept_RT(hard_idx),1);
                        
                        all_w2RT_scaled_RT_low_easy{ii_subj}    = mean(posi_fun(unz_PS( easy_RTlow_idx,6))  .*kept_RT(easy_RTlow_idx),1);
                        all_w2RT_scaled_RT_low_hard{ii_subj}    = mean(posi_fun(unz_PS( hard_RTlow_idx,6))  .*kept_RT(hard_RTlow_idx),1);
                        all_w2RT_scaled_RT_hig_easy{ii_subj}    = mean(posi_fun(unz_PS( easy_RThigh_idx ,6)).*kept_RT(easy_RThigh_idx),1);
                        all_w2RT_scaled_RT_hig_hard{ii_subj}    = mean(posi_fun(unz_PS( hard_RThigh_idx,6)) .*kept_RT(hard_RThigh_idx),1);
                        all_RT_and_P_easy{ii_subj} = [kept_RT(easy_idx),unz_PS( easy_idx,:)];
                        all_RT_and_P_hard{ii_subj} = [kept_RT(hard_idx),unz_PS( hard_idx,:)];
                    else
                        all_w2RT_scaled_easy{ii_subj}    = mean(posi_fun(unz_PS( easy_idx,6)),1);
                        all_w2RT_scaled_hard{ii_subj}    = mean(posi_fun(unz_PS( hard_idx,6)),1);
                        
                        all_w2RT_scaled_RT_low_easy{ii_subj}    = mean(posi_fun(unz_PS( easy_RTlow_idx,6))  ,1);
                        all_w2RT_scaled_RT_low_hard{ii_subj}    = mean(posi_fun(unz_PS( hard_RTlow_idx,6))  ,1);
                        all_w2RT_scaled_RT_hig_easy{ii_subj}    = mean(posi_fun(unz_PS( easy_RThigh_idx ,6)),1);
                        all_w2RT_scaled_RT_hig_hard{ii_subj}    = mean(posi_fun(unz_PS( hard_RThigh_idx,6)) ,1);
                        all_RT_and_P_easy{ii_subj} = [kept_RT(easy_idx),unz_PS( easy_idx,:)];
                        all_RT_and_P_hard{ii_subj} = [kept_RT(hard_idx),unz_PS( hard_idx,:)];
                    end
                    
                    % Fit quality
                    all_R2_means{ii_subj}       = get_mean_R2(simus, neur, easy_idx+ hard_idx);
                    all_R2_mean_easy{ii_subj}   = get_mean_R2(simus, neur, easy_idx);
                    all_R2_mean_hard{ii_subj}   = get_mean_R2(simus, neur, hard_idx);
                    
                    % M|X -> Y|X effects
                    % neural data
                    all_neur_mean_dev_easy_RT_low{ii_subj} = mean_fun(neur(  easy_RTlow_idx,:)  - removeXeffect*all_neur_mean_easy{ii_subj},1);
                    all_neur_mean_dev_easy_RT_high{ii_subj}= mean_fun(neur(  easy_RThigh_idx,:) - removeXeffect*all_neur_mean_easy{ii_subj},1);
                    
                    all_neur_mean_dev_hard_RT_low{ii_subj} = mean_fun(neur(  hard_RTlow_idx,:)  - removeXeffect*all_neur_mean_hard{ii_subj},1);
                    all_neur_mean_dev_hard_RT_high{ii_subj}= mean_fun(neur(  hard_RThigh_idx,:) - removeXeffect*all_neur_mean_hard{ii_subj},1);
                    
                    % fitted data
                    all_simu_mean_dev_easy_RT_low{ii_subj} = mean_fun(simus( easy_RTlow_idx,:)  - removeXeffect*all_simu_mean_easy{ii_subj} ,1);
                    all_simu_mean_dev_easy_RT_high{ii_subj}= mean_fun(simus( easy_RThigh_idx,:) - removeXeffect*all_simu_mean_easy{ii_subj} ,1);
                    
                    all_simu_mean_dev_hard_RT_low{ii_subj} = mean_fun(simus( hard_RTlow_idx,:)  - removeXeffect*all_simu_mean_hard{ii_subj} ,1);
                    all_simu_mean_dev_hard_RT_high{ii_subj}= mean_fun(simus( hard_RThigh_idx,:) - removeXeffect*all_simu_mean_hard{ii_subj} ,1);
                    
                    % Parameters
                    all_P_mean_easy_RT_low{ii_subj}    = mean_fun(Ps(    easy_RTlow_idx,:),1);
                    all_P_mean_easy_RT_high{ii_subj}   = mean_fun(Ps(    easy_RThigh_idx,:),1);
                    
                    all_P_mean_hard_RT_low{ii_subj}    = mean_fun(Ps(    hard_RTlow_idx,:),1);
                    all_P_mean_hard_RT_high{ii_subj}   = mean_fun(Ps(    hard_RThigh_idx,:),1);
                    
                    % Fit quality
                    all_R2_mean_easy_RT_low{ii_subj}   = get_mean_R2(simus, neur, easy_RTlow_idx);
                    all_R2_mean_easy_RT_high{ii_subj}  = get_mean_R2(simus, neur, easy_RThigh_idx);
                    
                    all_R2_mean_hard_RT_low{ii_subj}   = get_mean_R2(simus, neur, hard_RTlow_idx);
                    all_R2_mean_hard_RT_high{ii_subj}  = get_mean_R2(simus, neur, hard_RThigh_idx);
                    
                    % Uncomment the next line to see the fit quality (mean & std of trials' R2 per subject)
                    %                     disp(num2str([ii_subj,mean(arrayfun(@(ii) nancorr(simus(ii,:), neur(ii,:)),1:size(simus,1))),std(arrayfun(@(ii) nancorr(simus(ii,:), neur(ii,:)),1:size(simus,1)))],2))
                end
            end
            
            r2_easy{ii_data,selected_roi} = all_R2_mean_easy;
            r2_hard{ii_data,selected_roi} = all_R2_mean_hard;
            r2{ii_data,selected_roi} = all_R2_means;
            continue
            %% Data agregation
            nansem = @(x) nanstd(x)/sqrt(size(x,1)); % handy shortcut
            
            % Neural data
            neur_mean_easy        = nanmean(cat(1, all_neur_mean_easy{:}));
            neur_sem_easy         = nansem( cat(1, all_neur_mean_easy{:}));
            neur_mean_hard        = nanmean(cat(1, all_neur_mean_hard{:}));
            neur_sem_hard         = nansem( cat(1, all_neur_mean_hard{:}));
            
            neur_mean_easy_RT_low = nanmean(cat(1, all_neur_mean_dev_easy_RT_low{:}));
            neur_sem_easy_RT_low  = nansem( cat(1, all_neur_mean_dev_easy_RT_low{:}));
            neur_mean_hard_RT_low = nanmean(cat(1, all_neur_mean_dev_hard_RT_low{:}));
            neur_sem_hard_RT_low  = nansem( cat(1, all_neur_mean_dev_hard_RT_low{:}));
            
            neur_mean_easy_RT_high = nanmean(cat(1, all_neur_mean_dev_easy_RT_high{:}));
            neur_sem_easy_RT_high  = nansem( cat(1, all_neur_mean_dev_easy_RT_high{:}));
            neur_mean_hard_RT_high = nanmean(cat(1, all_neur_mean_dev_hard_RT_high{:}));
            neur_sem_hard_RT_high  = nansem( cat(1, all_neur_mean_dev_hard_RT_high{:}));
            
            % Fitted data
            simu_mean_easy         = nanmean(cat(1, all_simu_mean_easy{:}));
            simu_sem_easy          = nansem( cat(1, all_simu_mean_easy{:}));
            simu_mean_hard         = nanmean(cat(1, all_simu_mean_hard{:}));
            simu_sem_hard          = nansem( cat(1, all_simu_mean_hard{:}));
            
            simu_mean_easy_RT_low  = nanmean(cat(1, all_simu_mean_dev_easy_RT_low{:})) ;
            simu_sem_easy_RT_low   = nansem( cat(1, all_simu_mean_dev_easy_RT_low{:}));
            simu_mean_hard_RT_low  = nanmean(cat(1, all_simu_mean_dev_hard_RT_low{:}));
            simu_sem_hard_RT_low   = nansem( cat(1, all_simu_mean_dev_hard_RT_low{:}));
            
            simu_mean_easy_RT_high = nanmean(cat(1, all_simu_mean_dev_easy_RT_high{:}));
            simu_sem_easy_RT_high  = nansem( cat(1, all_simu_mean_dev_easy_RT_high{:}));
            simu_mean_hard_RT_high = nanmean(cat(1, all_simu_mean_dev_hard_RT_high{:}));
            simu_sem_hard_RT_high  = nansem( cat(1, all_simu_mean_dev_hard_RT_high{:}));
            
            % Shape parameters
            P_mean_easy            = nanmean(cat(1, all_P_mean_easy{:}));
            P_sem_easy             = nansem( cat(1, all_P_mean_easy{:}));
            P_mean_hard            = nanmean(cat(1, all_P_mean_hard{:}));
            P_sem_hard             = nansem( cat(1, all_P_mean_hard{:}));
            
            P_mean_easy_RT_low     = nanmean(cat(1, all_P_mean_easy_RT_low{:}));
            P_sem_easy_RT_low      = nansem( cat(1, all_P_mean_easy_RT_low{:}));
            P_mean_hard_RT_low     = nanmean(cat(1, all_P_mean_hard_RT_low{:}));
            P_sem_hard_RT_low      = nansem( cat(1, all_P_mean_hard_RT_low{:}));
            
            P_mean_easy_RT_high    = nanmean(cat(1, all_P_mean_easy_RT_high{:}));
            P_sem_easy_RT_high     = nansem( cat(1, all_P_mean_easy_RT_high{:}));
            P_mean_hard_RT_high    = nanmean(cat(1, all_P_mean_hard_RT_high{:}));
            P_sem_hard_RT_high     = nansem( cat(1, all_P_mean_hard_RT_high{:}));
            
            
            % suplemental shape parameters info
            P_mean =  (cat(1, all_P_mean{:}));
            w2RT_scaled_easy = cat(1, all_w2RT_scaled_easy{:});
            w2RT_scaled_hard = cat(1, all_w2RT_scaled_hard{:});
            
            w2RT_scaled_RT_low_easy = cat(1, all_w2RT_scaled_RT_low_easy{:});
            w2RT_scaled_RT_low_hard = cat(1, all_w2RT_scaled_RT_low_hard{:});
            w2RT_scaled_RT_hig_easy = cat(1, all_w2RT_scaled_RT_hig_easy{:});
            w2RT_scaled_RT_hig_hard = cat(1, all_w2RT_scaled_RT_hig_hard{:});
            
            % Fit quality
            R2_mean_easy           = nanmean(cat(1, all_R2_mean_easy{:}));
            R2_sem_easy            = nansem( cat(1, all_R2_mean_easy{:}));
            R2_mean_hard           = nanmean(cat(1, all_R2_mean_hard{:}));
            R2_sem_hard            = nansem( cat(1, all_R2_mean_hard{:}));
            
            R2_mean_easy_RT_low    = nanmean(cat(1, all_R2_mean_easy_RT_low{:}));
            R2_sem_easy_RT_low     = nansem( cat(1, all_R2_mean_easy_RT_low{:}));
            R2_mean_hard_RT_low    = nanmean(cat(1, all_R2_mean_hard_RT_low{:}));
            R2_sem_hard_RT_low     = nansem( cat(1, all_R2_mean_hard_RT_low{:}));
            
            R2_mean_easy_RT_high   = nanmean(cat(1, all_R2_mean_easy_RT_high{:}));
            R2_sem_easy_RT_high    = nansem( cat(1, all_R2_mean_easy_RT_high{:}));
            R2_mean_hard_RT_high   = nanmean(cat(1, all_R2_mean_hard_RT_high{:}));
            R2_sem_hard_RT_high    = nansem( cat(1, all_R2_mean_hard_RT_high{:}));
            
            n_subj_here            = sum(~cellfun(@isempty,all_R2_mean_easy));
            
            %% Plot
            try
                h_fig               = figure();
                h_fig.Color         = 'w';
                h_fig.Position(2)   = h_fig.Position(2) + h_fig.Position(4) - 850;
                h_fig.Position(3:4) = [1200, 850];
                clf;
                % This part uses the shaded error bar of https://fr.mathworks.com/matlabcentral/fileexchange/26311-raacampbell-shadederrorbar
                
                % Prepare the figure's labels
                switch ii_data
                    case 1 % Completion time
                        time_list     = percentage_time/100;
                        if NORMALIZU_VISUALISATION
                            my_xlim       = time_list([round(end*0.3),round(end*0.7)]);
                            my_xticks     = time_list(round(end/2))+linspace(-0.6,0.6,13);
                            tmp           = 100*(my_xticks-time_list(round(end/2)));
                            tmp(abs(tmp) < 1e-4) = 0;
                            my_xticks_lbl = [arrayfun(@(x){['c',num2str(x,2),'%']} ,tmp(tmp<0)), {'c'},arrayfun(@(x){['c+',num2str(x,2),'%']} ,tmp(tmp>0))];
                            
                        else
                            my_xlim = time_list([1, end]);
                            my_xticks     = [-10, 0, linspace(25,200,8)]/100;
                            my_xticks_lbl = arrayfun(@(x) {[num2str(x*100),'%']},my_xticks);
                        end
                    case {2,3}
                        time_list     = trial_time;
                        if NORMALIZU_VISUALISATION
                            my_xlim       = time_list([round(end*0.2),round(end*0.8)]);
                            my_xticks     = [-1.5,-1,-0.5,0,0.5,1,1.5]; % relative values
                            my_xticks_lbl = [arrayfun(@(x){['c',num2str(x,2),'s']} ,my_xticks(my_xticks<0)), {'c'},arrayfun(@(x){['c+',num2str(x,2),'s']} ,my_xticks(my_xticks>0))];
                            my_xticks     = time_list(round(end/2))+ my_xticks; % true values
                        else
                            my_xlim = time_list([1, end]);
                            my_xticks     = [-2,-1,0,1,2,3];
                            my_xticks_lbl = arrayfun(@(x) {num2str(x)},my_xticks);
                        end
                end
                
                % Stat report
                subplot(8,2,1);
                my_ax = gca;my_ax.Position(1:2) = [0.05,0.90];
                roi_name = strrep(all_rois_names{selected_roi},'_',' ');
                if isempty(stat_report{ii_data,selected_roi}), text(0,0,[{roi_name},{'Shape Mediation statistics : NOT SIGNIFICANT'}],'Fontsize',7)
                else,                                          text(0,0,[{roi_name},{'Shape Mediation statistics :'}, ...
                        strrep(strrep(     stat_report{ii_data,selected_roi}, '<strong>','{\bf '),'</strong>',' }')
                        ],'Fontsize',7)
                end
                axis off
                
                % Plot of Easy vs Hard trials
                subplot(8,2,[3,5]); hold on
                h1 = shadedErrorBar(time_list, neur_mean_easy, neur_sem_easy,      {'color',color_easy,     'Linewidth',2}, 1);
                h2 = shadedErrorBar(time_list, neur_mean_hard, neur_sem_hard,      {'color',color_hard,     'Linewidth',2}, 1);
                h3 = plot(          time_list, simu_mean_easy,'o-', 'MarkerSize',2, 'color',color_easy*0.8, 'Linewidth',1, 'MarkerFaceColor',color_easy*0.8);
                h4 = plot(          time_list, simu_mean_hard,'o-', 'MarkerSize',2, 'color',color_hard*0.8, 'Linewidth',1, 'MarkerFaceColor',color_hard*0.8);
                
                xticks(my_xticks)
                xticklabels(my_xticks_lbl)
                axis tight
                
                if ~NORMALIZU_VISUALISATION
                    plot([0,0],ylim(),'k')                            % Stim or Response line
                    if ii_data == 1, plot([1,1],ylim(),'-k');         % Response line
                    end
                    ylabel('sEEG high-Gamma')
                else
                    plot([1,1]*time_list(round(end/2)),ylim(),'--k'); % peak line
                    ylabel('Start&Peak-normalized data')
                end
                switch ii_data
                    case 1 , xlabel('Completion time')
                    case 2 , xlabel('Stim-locked time')
                    case 3 , xlabel('Resp-locked time')
                end
                xlim(my_xlim)
                
                legend([h1.mainLine,h3,h2.mainLine,h4],{'Easy','Model','Hard','Model'},'box','off','Location','northwest')
                title('sEEG easy vs hard');
                grid on
                
                subplot(8,2,7); hold on
                text(0,0.25,{['Shape Fit Quality (n subj =',num2str(n_subj_here),'):'],...
                    ['R2 easy:               ',num2str(R2_mean_easy,2),' (+-',num2str(R2_sem_easy,2),')',...
                    '         R2 hard:              ', num2str(R2_mean_hard,2),' (+-',num2str(R2_sem_hard,2),')'], ...
                    ['R2 easy RT low:   ',num2str(R2_mean_easy_RT_low,2),' (+-',num2str(R2_sem_easy_RT_low,2),')',...
                    '       R2 hard RT low:  ',num2str(R2_mean_hard_RT_low,2),' (+-',num2str(R2_sem_hard_RT_low,2),')'], ...
                    ['R2 easy RT high: ',num2str(R2_mean_easy_RT_high,2),' (+-',num2str(R2_sem_easy_RT_high,2),')',...
                    '      R2 hard RT high: ',num2str(R2_mean_hard_RT_high,2),' (+-',num2str(R2_sem_hard_RT_high,2),')'], ...
                    },'Fontsize',7)
                axis off
                
                % Plot of RT low vs high in easy trials
                subplot(8,2,[2,4]); hold on
                h1 = shadedErrorBar( time_list, neur_mean_easy_RT_low ,neur_sem_easy_RT_low,   {'color',color_easyRTlow,     'Linewidth',2}, 1);
                h2 = shadedErrorBar( time_list, neur_mean_easy_RT_high,neur_sem_easy_RT_high,  {'color',color_easyRThigh,    'Linewidth',2}, 1);
                h3 = plot(           time_list, simu_mean_easy_RT_low,    'o-', 'MarkerSize',2, 'color',color_easyRTlow*0.8,  'Linewidth',1, 'MarkerFaceColor',color_easyRTlow *0.8);
                h4 = plot(           time_list, simu_mean_easy_RT_high,   'o-', 'MarkerSize',2, 'color',color_easyRThigh*0.8, 'Linewidth',1, 'MarkerFaceColor',color_easyRThigh*0.8);
                
                xticks(my_xticks)
                xticklabels(my_xticks_lbl)
                
                axis tight
                
                if ~NORMALIZU_VISUALISATION
                    plot([0,0],ylim(),'k')                            % Stim or Response line
                    if ii_data == 1, plot([1,1],ylim(),'-k');         % Response line
                    end
                    ylabel('sEEG high-Gamma')
                else
                    plot([1,1]*time_list(round(end/2)),ylim(),'--k'); % peak line
                    ylabel('Start&Peak-normalized data')
                end
                xlim(my_xlim)
                if removeXeffect, plot(xlim(),[0,0],'k'); end % Line of null deviation
                title('sEEG deviation from mean Easy-trials');
                
                legend([h1.mainLine,h3,h2.mainLine,h4],{'Easy RT-','Model','Easy RT+','Model'},'box','off','Location','southwest')
                text(max(xlim())*1.01, max(ylim())/2, '<-- Above mean', 'Rotation',-90, 'VerticalAlignment','bottom', 'HorizontalAlignment','center');
                text(max(xlim())*1.01, min(ylim())/2, 'Below mean -->', 'Rotation',-90, 'VerticalAlignment','bottom', 'HorizontalAlignment','center');
                grid on
                
                % Plot of RT low vs high in hard trials
                subplot(8,2,[6,8]); hold on
                h1 = shadedErrorBar(time_list, neur_mean_hard_RT_low, neur_sem_hard_RT_low,  {'color',color_hardRTlow, 'Linewidth',2}, 1);
                h2 = shadedErrorBar(time_list, neur_mean_hard_RT_high,neur_sem_hard_RT_high, {'color',color_hardRThigh,'Linewidth',2}, 1);
                h3 = plot(          time_list, simu_mean_hard_RT_low, 'o-','MarkerSize',2,    'color',color_hardRTlow*0.8 ,  'Linewidth',2, 'MarkerFaceColor',color_hardRTlow *0.8);
                h4 = plot(          time_list, simu_mean_hard_RT_high,'o-','MarkerSize',2,    'color',color_hardRThigh*0.8,  'Linewidth',2, 'MarkerFaceColor',color_hardRThigh*0.8);
                xticks(my_xticks)
                xticklabels(my_xticks_lbl)
                
                axis tight
                if ~NORMALIZU_VISUALISATION
                    plot([0,0],ylim(),'k')                            % Stim or Response line
                    if ii_data == 1, plot([1,1],ylim(),'-k');         % Response line
                    end
                    ylabel('sEEG high-Gamma')
                else
                    plot([1,1]*time_list(round(end/2)),ylim(),'--k'); % peak line
                    ylabel('Start&Peak-normalized data')
                end
                switch ii_data
                    case 1 , xlabel('Completion time')
                    case 2 , xlabel('Stim-locked time')
                    case 3 , xlabel('Resp-locked time')
                end
                xlim(my_xlim)
                if removeXeffect, plot(xlim(),[0,0],'k'); end % Line of null deviation
                
                title('sEEG deviation from mean Hard-trials', 'VerticalAlignment','top');
                legend([h1.mainLine,h3,h2.mainLine,h4], {'Hard RT-','Model','Hard RT+','Model'}, 'box','off', 'Location','southwest')
                text(max(xlim())*1.01, max(ylim())/2, '<-- Above mean', 'Rotation',-90, 'VerticalAlignment','bottom', 'HorizontalAlignment','center');
                text(max(xlim())*1.01, min(ylim())/2, 'Below mean -->', 'Rotation',-90, 'VerticalAlignment','bottom', 'HorizontalAlignment','center');
                grid on
                
                % New figure
                set(figure(),'color','w')
                % Plot of Easy vs Hard parameters
                
                subplot(4,2,[9 11]-8); hold on
                h1 = bar( (1:8)-0.025, P_mean_easy,0.3);             h1.FaceColor = color_easy;
                errorbar( (1:8)-0.025, P_mean_easy,P_sem_hard,'.k');
                h2 = bar( (1:8)+0.325, P_mean_hard,0.3);             h2.FaceColor = color_hard;
                errorbar( (1:8)+0.325, P_mean_hard,P_sem_hard,'.k')
                
                xticks((1:8)+0.3/2)
                xticklabels(param_name)
                ylabel('Zscored parameters')
                legend([h1,h2],{'Mean Easy','Mean Hard'},'box','off','Location','southwest')
                grid on
                
                % Plot the parameters separated according to easy vs hard and RT low vs high
                subplot(4,2,[10 12]-8); hold on
                h1 = bar( (1:8) - 0.025,        P_mean_easy_RT_low,0.15);      h1.FaceColor=color_easyRTlow;
                errorbar( (1:8) - 0.025,        P_mean_easy_RT_low, P_sem_easy_RT_low,'.k')
                h2 = bar( (1:8) + 0.125,        P_mean_easy_RT_high,0.15);     h2.FaceColor= color_easyRThigh;
                errorbar( (1:8) + 0.125,        P_mean_easy_RT_high, P_sem_easy_RT_high,'.k')
                
                h3 = bar( (1:8)+0.15 *2 +0.025, P_mean_hard_RT_low,0.15);      h3.FaceColor= color_hardRTlow;
                errorbar( (1:8)+0.15 *2 +0.025, P_mean_hard_RT_low, P_sem_hard_RT_low,'.k')
                h4 = bar( (1:8)+0.15 *3 +0.025, P_mean_hard_RT_high,0.15);     h4.FaceColor=color_hardRThigh;
                errorbar( (1:8)+0.15 *3 +0.025, P_mean_hard_RT_high, P_sem_hard_RT_high,'.k')
                
                xticks((1:8)+0.225)
                xticklabels(param_name)
                ylabel('Zscored parameters')
                legend([h1,h2,h3,h4],{'Easy RT-','Easy RT+','Hard RT-','Hard RT+'},'box','off','Location','southwest')
                grid on
                
                
                % Plot param value
                if 1
                    P_mean =  (cat(1, all_P_mean{:}));
                    
                    left_idx = [1 3 4 6 7];
                    right_idx = [2 5 8];
                    subplot(4,4,9); hold on
                    h1 = bar( 1:numel(left_idx), mean(P_mean(:,left_idx)),0.3*2);       h1.FaceColor = [1 1 1]*0.33;
                    errorbar( 1:numel(left_idx), mean(P_mean(:,left_idx)),nansem(P_mean(:,left_idx)),'.k');
                    ylabel('Parameters')
                    xticks(1:numel(left_idx) )
                    xticklabels(param_name(left_idx))
                    title('Population values')
                    grid on
                    
                    subplot(4,4,10); hold on
                    
                    h1 = bar( (1:numel(right_idx)), mean(P_mean(:,right_idx)),0.3*2);   h1.FaceColor = [1 1 1]*0.33;
                    errorbar( (1:numel(right_idx)), mean(P_mean(:,right_idx)),nansem(P_mean(:,right_idx)),'.k');
                    ylabel({'Amplitude parameters'})
                    xticks(1:numel(right_idx))
                    xticklabels(param_name(right_idx))
                    title('Population values')
                    grid on
                    
                end
                
                % plot w2 rescaled
                subplot(4,2,16-8); hold on
                
                
                h1 = bar( 1, mean(w2RT_scaled_easy),0.3*3);             h1.FaceColor = color_easy;
                errorbar( 1, mean(w2RT_scaled_easy), nansem(w2RT_scaled_easy),'.k');
                h2 = bar( 2, mean(w2RT_scaled_hard),0.3*3);             h2.FaceColor = color_hard ;
                errorbar( 2, mean(w2RT_scaled_hard), nansem(w2RT_scaled_hard),'.k');
                
                
                h11 = bar( 4-0.025*2, mean(w2RT_scaled_RT_low_easy),0.15*2);             h11.FaceColor = color_easyRTlow;
                errorbar( 4-0.025*2, mean(w2RT_scaled_RT_low_easy), nansem(w2RT_scaled_RT_low_easy),'.k');
                h12 = bar( 4+0.125*2, mean(w2RT_scaled_RT_hig_easy),0.15*2);             h12.FaceColor = color_easyRThigh;
                errorbar( 4+0.125*2, mean(w2RT_scaled_RT_hig_easy), nansem(w2RT_scaled_RT_hig_easy),'.k');
                
                h21 = bar( 5-0.025*2, mean(w2RT_scaled_RT_low_hard),0.15*2);             h21.FaceColor = color_hardRTlow;
                errorbar( 5-0.025*2, mean(w2RT_scaled_RT_low_hard), nansem(w2RT_scaled_RT_low_hard),'.k');
                h22 = bar( 5+0.125*2, mean(w2RT_scaled_RT_hig_hard),0.15*2);             h22.FaceColor = color_hardRThigh;
                errorbar( 5+0.125*2, mean(w2RT_scaled_RT_hig_hard), nansem(w2RT_scaled_RT_hig_hard),'.k');
                
                title('Neural return speed (in sec.)')
                xlim([0,11])
                legend([h1,h11,h12,h2,h21,h22],{['(easy) w2*RT: ',num2str(mean(w2RT_scaled_easy),2),'s'],
                    ['(easy RT-) w2*RT: ',num2str(mean(w2RT_scaled_RT_low_easy),2),'s'],
                    ['(easy RT+) w2*RT: ',num2str(mean(w2RT_scaled_RT_hig_easy),2),'s'],
                    ['(hard) w2*RT: ',num2str(mean(w2RT_scaled_hard),2),'s'],
                    ['(hard RT-) w2*RT: ',num2str(mean(w2RT_scaled_RT_low_hard),2),'s'],
                    ['(hard RT+) w2*RT: ',num2str(mean(w2RT_scaled_RT_hig_hard),2),'s']},...
                    'box','off','Location','northeast')
                xticks([1.5,4.05,5.05])
                xticklabels({'Global','Easy','Hard'})
                ylabel('w2')
                
                subplot(4,2,14-8)
                hold on
                [RT_easy_low_param,RT_easy_inter_param,RT_easy_inter2_param ,RT_hard_low_param,RT_hard_inter_param,RT_hard_inter2_param,RT_easy_high_param ,RT_hard_high_param] = deal(NaN(n_subj,1));
                for ii_subj=1:n_subj
                    if ~isempty(all_RT_and_P_easy{ii_subj})
                        param_per_trial = all_RT_and_P_easy{ii_subj}(:,1+6);
                        param_low_idx  = param_per_trial < prctile(param_per_trial,25);
                        param_inter_idx  = (param_per_trial > prctile(param_per_trial,25)) & (param_per_trial < prctile(param_per_trial,50));
                        param_inter2_idx  = (param_per_trial > prctile(param_per_trial,50)) & (param_per_trial < prctile(param_per_trial,75));
                        param_high_idx = param_per_trial > prctile(param_per_trial,75);
                        RT_easy_low_param(ii_subj) = mean(all_RT_and_P_easy{ii_subj}(param_low_idx,1));
                        RT_easy_inter_param(ii_subj) = mean(all_RT_and_P_easy{ii_subj}(param_inter_idx,1));
                        RT_easy_inter2_param(ii_subj) = mean(all_RT_and_P_easy{ii_subj}(param_inter2_idx,1));
                        RT_easy_high_param(ii_subj) = mean(all_RT_and_P_easy{ii_subj}(param_high_idx,1));
                        
                        param_per_trial = all_RT_and_P_hard{ii_subj}(:,1+6);
                        param_low_idx  = param_per_trial < prctile(param_per_trial,25);
                        param_inter_idx  = (param_per_trial > prctile(param_per_trial,25)) & (param_per_trial < prctile(param_per_trial,50));
                        param_inter2_idx  = (param_per_trial > prctile(param_per_trial,50)) & (param_per_trial < prctile(param_per_trial,75));
                        param_high_idx = param_per_trial > prctile(param_per_trial,75);
                        RT_hard_low_param(ii_subj) = mean(all_RT_and_P_hard{ii_subj}(param_low_idx,1));
                        RT_hard_inter_param(ii_subj) = mean(all_RT_and_P_hard{ii_subj}(param_inter_idx,1));
                        RT_hard_inter2_param(ii_subj) = mean(all_RT_and_P_hard{ii_subj}(param_inter2_idx,1));
                        RT_hard_high_param(ii_subj) = mean(all_RT_and_P_hard{ii_subj}(param_high_idx,1));
                    end
                end
                my_bar_data = [RT_easy_low_param,RT_easy_inter_param,RT_easy_inter2_param,RT_easy_high_param,...
                    RT_hard_low_param,RT_hard_inter_param,RT_hard_inter2_param,RT_hard_high_param];
                xs_1 = 1:4;
                xs_2 = 6:9;
                b1 = bar(xs_1,nanmean(my_bar_data(:,1:4)),'FaceColor',color_easy);
                errorbar(xs_1,nanmean(my_bar_data(:,1:4)),...
                    nanstd(my_bar_data(:,1:4)) / sqrt(sum(~isnan(RT_easy_low_param))),'.k')
                b2 = bar(xs_2,nanmean(my_bar_data(:,5:8)),'FaceColor',color_hard);
                errorbar(xs_2,nanmean(my_bar_data(:,5:8)),...
                    nanstd(my_bar_data(:,5:8)) / sqrt(sum(~isnan(RT_easy_low_param))),'.k')
                xticks([1:4 6:9])
                xticklabels({'0-25%','25-50%','50-75%','75-100%','0-25%','25-50%','50-75%','75-100%'})
                legend([b1,b2],{'Easy','Hard'},'box','off','location','northwest')
                title('RT per percentile of w_2')
                grid on
                xlabel('w_2''s quartiles')
                set(gca,'xticklabelrotation',40)
                ylabel('RT (log)')
                set(gca,'yscale','log')
                
                
                drawnow
                % Uncomment the following line to save figures as they are generated
                %             saveas(gcf,fullfile(figsavefolder,[all_rois_names{selected_roi},'.png']))
            end
        end
    end
end


%% EXTRA: saving variable of interest as table
if 0
    [Mean_value_a,Mean_value_b,sem_value_a,sem_value_b]  = deal(NaN(8,n_roi));
    for ii_param=1:8
        Mean_value_a(ii_param,:)  =  nanmean(cellfun(@(x)get_param(x,ii_param),squeeze(Shape_a(1,:,:))));
        sem_value_a(ii_param,:)  =  nanstd(cellfun(@(x)get_param(x,ii_param),squeeze(Shape_a(1,:,:)))) ./ sqrt(sum(~isnan(cellfun(@(x)get_param(x,ii_param),squeeze(Shape_a(1,:,:))))));
        Mean_value_b(ii_param,:)  =  nanmean(cellfun(@(x)get_param(x,ii_param),squeeze(Shape_b(1,:,:))));
        sem_value_b(ii_param,:)  =  nanstd(cellfun(@(x)get_param(x,ii_param),squeeze(Shape_b(1,:,:)))) ./ sqrt(sum(~isnan(cellfun(@(x)get_param(x,ii_param),squeeze(Shape_b(1,:,:))))));
    end
    table_a = rows2vars(array2table(arrayfun(@(x,y) {[num2str(x,2) ' +-' num2str(y,2)]},Mean_value_a,sem_value_a), 'VariableNames',all_rois_names, 'RowNames', param_name));
    table_b = rows2vars(array2table(arrayfun(@(x,y) {[num2str(x,2) ' +-' num2str(y,2)]},Mean_value_b,sem_value_b), 'VariableNames',all_rois_names, 'RowNames', param_name));
    writetable(table_a,'table_mean_sem_path_diff_SEEGshape.csv')
    writetable(table_b,'table_mean_sem_path_SEEGshape_RT.csv')
end

if 0 % set to true to extract parameters into files
    Subj_labels  = arrayfun(@(ii){['Subj n°',num2str(ii)]},1:n_subj);
    
    Table_ComplTime_path_Diffi_to_ShapeParamSEEG = array2table(squeeze( Shape_a(1,:,:)), 'VariableNames',all_rois_names, 'RowNames',Subj_labels);
    Table_SLockTime_path_Diffi_to_ShapeParamSEEG = array2table(squeeze( Shape_a(2,:,:)), 'VariableNames',all_rois_names, 'RowNames',Subj_labels);
    Table_RLockTime_path_Diffi_to_ShapeParamSEEG = array2table(squeeze( Shape_a(3,:,:)), 'VariableNames',all_rois_names, 'RowNames',Subj_labels);
    
    Table_ComplTime_path_ShapeParamSEEG_to_RT = array2table(squeeze( Shape_b(1,:,:)), 'VariableNames',all_rois_names, 'RowNames',Subj_labels);
    Table_SLockTime_path_ShapeParamSEEG_to_RT = array2table(squeeze( Shape_b(2,:,:)), 'VariableNames',all_rois_names, 'RowNames',Subj_labels);
    Table_RLockTime_path_ShapeParamSEEG_to_RT = array2table(squeeze( Shape_b(3,:,:)), 'VariableNames',all_rois_names, 'RowNames',Subj_labels);
    
    concatenated_P = cell(size(subj_roi_P));
    for ii_data=1:size(concatenated_P,1)
        for ii_subj=1:size(concatenated_P,2)
            for ii_roi=1:size(concatenated_P,3)
                if ~isempty(subj_roi_P{ii_data,ii_subj,ii_roi})
                    concatenated_P{ii_data,ii_subj,ii_roi} = cat(2,subj_roi_P{ii_data,ii_subj,ii_roi}{:});
                end
            end
        end
    end
    
    [   Matrix_ComplTime_path_Diffi_to_ShapeParamSEEG, Matrix_ComplTime_path_ShapeParamSEEG_to_RT,...
        Matrix_SLockTime_path_Diffi_to_ShapeParamSEEG, Matrix_SLockTime_path_ShapeParamSEEG_to_RT,...
        Matrix_RLockTime_path_Diffi_to_ShapeParamSEEG, Matrix_RLockTime_path_ShapeParamSEEG_to_RT] = deal(nan(n_subj,n_roi,8));
    for ii_subj=1:size(concatenated_P,2)
        for ii_roi=1:size(concatenated_P,3)
            if ~isempty(Shape_a{1,ii_subj,ii_roi})
                Matrix_ComplTime_path_Diffi_to_ShapeParamSEEG(ii_subj,ii_roi,:) = Shape_a{1,ii_subj,ii_roi};
                Matrix_ComplTime_path_ShapeParamSEEG_to_RT(ii_subj,ii_roi,:)    = Shape_b{1,ii_subj,ii_roi};
            end
            if ~isempty(Shape_a{2,ii_subj,ii_roi})
                Matrix_SLockTime_path_Diffi_to_ShapeParamSEEG(ii_subj,ii_roi,:) = Shape_a{2,ii_subj,ii_roi};
                Matrix_SLockTime_path_ShapeParamSEEG_to_RT(ii_subj,ii_roi,:)    = Shape_b{2,ii_subj,ii_roi};
            end
            if ~isempty(Shape_a{3,ii_subj,ii_roi})
                Matrix_RLockTime_path_Diffi_to_ShapeParamSEEG(ii_subj,ii_roi,:) = Shape_a{3,ii_subj,ii_roi};
                Matrix_RLockTime_path_ShapeParamSEEG_to_RT(ii_subj,ii_roi,:)    = Shape_b{3,ii_subj,ii_roi};
            end
        end
    end
    
    
    Shape_parameter_name_and_order = param_name;
    Shape_parameter_description    = {'c: Peak time','A: Peak amplitude','w1: Integration window size','p1: Integration concavity','b1: Starting baseline','w2: Depletion window size','p2: Depletion concavity','b2: Finishing baseline'};
    ROI_names_and_order = all_rois_names;
    
    save('ShapeMediationCoefficients_Subj-by-ROI-by-Param.mat',...
        'Table_ComplTime_path_Diffi_to_ShapeParamSEEG', 'Table_SLockTime_path_Diffi_to_ShapeParamSEEG', 'Table_RLockTime_path_Diffi_to_ShapeParamSEEG',...
        'Table_ComplTime_path_ShapeParamSEEG_to_RT',    'Table_SLockTime_path_ShapeParamSEEG_to_RT',    'Table_RLockTime_path_ShapeParamSEEG_to_RT',...
        'Shape_parameter_name_and_order','Shape_parameter_description');
    
    save('ShapeMediationCoefficients_MATRIX_Subj-by-ROI-by-Param.mat',...
        'Matrix_ComplTime_path_Diffi_to_ShapeParamSEEG', 'Matrix_SLockTime_path_Diffi_to_ShapeParamSEEG', 'Matrix_RLockTime_path_Diffi_to_ShapeParamSEEG',...
        'Matrix_ComplTime_path_ShapeParamSEEG_to_RT',    'Matrix_SLockTime_path_ShapeParamSEEG_to_RT',    'Matrix_RLockTime_path_ShapeParamSEEG_to_RT',...
        'ROI_names_and_order','Shape_parameter_name_and_order','Shape_parameter_description');
    
    Table_ComplTime_ShapeParamSEEG = array2table(squeeze( concatenated_P(1,:,:)), 'VariableNames',all_rois_names, 'RowNames',Subj_labels);
    Table_SlockTime_ShapeParamSEEG = array2table(squeeze( concatenated_P(1,:,:)), 'VariableNames',all_rois_names, 'RowNames',Subj_labels);
    Table_RlockTime_ShapeParamSEEG = array2table(squeeze( concatenated_P(1,:,:)), 'VariableNames',all_rois_names, 'RowNames',Subj_labels);
    
    save('ComplTime_ShapeParamSEEG_Subj-by-ROI-by-Param-by-trials.mat', 'Table_ComplTime_ShapeParamSEEG','Shape_parameter_name_and_order','Shape_parameter_description');
    save('SlockTime_ShapeParamSEEG_Subj-by-ROI-by-Param-by-trials.mat', 'Table_SlockTime_ShapeParamSEEG','Shape_parameter_name_and_order','Shape_parameter_description');
    save('RlockTime_ShapeParamSEEG_Subj-by-ROI-by-Param-by-trials.mat', 'Table_RlockTime_ShapeParamSEEG','Shape_parameter_name_and_order','Shape_parameter_description');
    
    
    for ii_param=1:numel(param_name)
        Table_ComplTime_path_Diffi_to_ShapeParamSEEG = array2table(cellfun(@(x)get_param(x,ii_param),squeeze( Shape_a(1,:,:))), 'VariableNames',all_rois_names, 'RowNames',Subj_labels);
        Table_ComplTime_path_ShapeParamSEEG_to_RT    = array2table(cellfun(@mean,squeeze( Shape_b(1,:,:))), 'VariableNames',all_rois_names, 'RowNames',Subj_labels);
        writetable(Table_ComplTime_path_Diffi_to_ShapeParamSEEG,  [strrep(param_name{ii_param},' ',''), '_Table_ComplTime_path_Diffi_to_ShapeParamSEEG.csv'])
        writetable(Table_ComplTime_path_ShapeParamSEEG_to_RT, [strrep(param_name{ii_param},' ',''), '_Table_ComplTime_path_ShapeParamSEEG_to_RT.csv'])
    end
end

%% Miscelaneous
function out=get_param(x,i)
if ~isempty(x)
    out= x(i);
else
    out= NaN;
end
end

function out = posi_fun(x) % taken from "g_shape_quick.m"
out = log(1+exp(50.*x)) ./50 + 1e-5;
if isinf(out), out = x; end
end
