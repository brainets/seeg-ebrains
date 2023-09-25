%% Global settings
REPROCESS_DATA             = false; % set to false to load previously preprocessed data
REFIT_ELECTRODES           = true; % set to false to load previously fited data
min_percentage_before_stim = -10; % to capture the pre-stim dynamic
max_percentage_after_stim  = 200; % to have the same window before and after the response

ROI_responsiveness_p_value_threshold = 1; % Threshold to select which ROI will be fitted and analysed

% misc variables
percentage_time = min_percentage_before_stim:max_percentage_after_stim;

%% prepare folder path
script_name = mfilename('fullpath');
script_folder = strsplit(script_name,filesep);
script_folder = strjoin(script_folder(1:end-1),filesep);

event_folder = fullfile(script_folder,'analyses','events');
epoch_folder = fullfile(script_folder,'analyses','epochs',...
    'sample_stim','f50f150-sm0');
anat_folder  = fullfile(epoch_folder,'anat');
data_folder  = fullfile(epoch_folder,'data');
data_folder_Rlocked  = fullfile(script_folder,'analyses','epochs','sample_resp','f50f150-sm0',...
    'data');
assert(prod(cellfun(@(x) isfolder(x), {event_folder,epoch_folder,anat_folder,data_folder,data_folder_Rlocked}))==1, ...
    ['Folders are not recovered correctly. In the same folder than this script there should be an ''analyses'' with the '...
    'folders: ''events'', ''epoch'' as produced by Etienne''s preprocessing scripts.'])

%% PREPROC Completeness check
% Check that every subject has the required files
event_file_list = dir(fullfile(event_folder,'*.xlsx'));
anat_file_list  = dir(fullfile(anat_folder,'*.xlsx'));
data_file_list  = dir(fullfile(data_folder,'*.nc'));
dataR_file_list = dir(fullfile(data_folder_Rlocked,'*.nc'));

subj_list1 = arrayfun(@(file) str2num(file.name(end-7:end-5)),event_file_list);
subj_list2 = arrayfun(@(file) str2num(file.name(end-7:end-5)),anat_file_list);
subj_list3 = arrayfun(@(file) str2num(file.name(9:11)),data_file_list);
subj_list4 = arrayfun(@(file) str2num(file.name(9:11)),dataR_file_list);

good_subj = intersect(intersect(intersect( ...
    arrayfun(@num2str,subj_list1,'UniformOutput',false), ...
    arrayfun(@num2str,subj_list2,'UniformOutput',false)) , ...
    arrayfun(@num2str,subj_list3,'UniformOutput',false)) , ...
    arrayfun(@num2str,subj_list4,'UniformOutput',false));

good_subj = sort(cellfun(@str2num,good_subj));

subject_info = struct();
for ii=1:numel(good_subj)
    ii1 = find(subj_list1==good_subj(ii),1,'first');
    ii2 = find(subj_list2==good_subj(ii),1,'first');
    ii3 = find(subj_list3==good_subj(ii),1,'first');
    ii4 = find(subj_list4==good_subj(ii),1,'first');
    subject_info(ii).num = ii;
    subject_info(ii).files.events = fullfile(event_file_list(ii1).folder, event_file_list(ii1).name);
    subject_info(ii).files.anat   = fullfile( anat_file_list(ii2).folder,  anat_file_list(ii2).name);
    subject_info(ii).files.dataS  = fullfile( data_file_list(ii3).folder,  data_file_list(ii3).name);
    subject_info(ii).files.dataR  = fullfile( dataR_file_list(ii4).folder,  dataR_file_list(ii4).name);
end

%% PREPROC Load data
% Load the behavioral and neural data in a matlab-convenient structures
for ii = 1:numel(subject_info)
    event_fn = subject_info(ii).files.events;
    anat_fn  = subject_info(ii).files.anat;
    data_fn  = subject_info(ii).files.dataS;
    dataR_fn = subject_info(ii).files.dataR;
    
    warning off % to get rid of the back-compatibility message
    [~,~,event_sheet]  =  xlsread(event_fn);
    [~,~,anat_sheet]   =  xlsread(anat_fn);
    warning on
    
    data = struct();
    
    % load event xls
    varnames = event_sheet(1,2:end);
    for jj=1:numel(varnames)
        var_name = strrep(varnames{jj},' ','_');
        data.event.(var_name) = event_sheet(2:end,jj+1);
    end
    
    % load anat xls
    varnames = anat_sheet(1,2:end);
    for jj=1:numel(varnames)
        var_name = strrep(varnames{jj},' ','_');
        data.event.(var_name) = anat_sheet(2:end,jj+1);
    end
    
    % load stim-locked data
    data_info = ncinfo(data_fn);
    varnames = arrayfun(@(x) {x.Name}, data_info.Variables);
    for jj=1:numel(varnames)
        var_name = strrep(varnames{jj},' ','_');
        data.neural.(var_name) = ncread(data_fn, varnames{jj});
    end
    
    % load resp-locked data (only the part that differs from the
    % stim-locked data)
    data2_info = ncinfo(dataR_fn);
    data.neural.f50f150_R = ncread(dataR_fn, 'f50f150');
    
    subject_info(ii).data = data;
end

% reformat
for ii = 1:numel(subject_info)
    neural = subject_info(ii).data.neural;
    
    neural.roi           = cellstr(neural.roi')';
    neural.contacts      = cellstr(neural.contacts');
    neural.T_orientation = cellstr(neural.T_orientation');
    neural.trials        = cellstr(neural.trials');
    neural.difficulty    = cellstr(neural.difficulty');
    
    subject_info(ii).data.neural = neural;
end

%% PREPROC Reformating
% Reformat data according to each time frame (especially for the
% 'Completion-time' time frame)
all_rois_names = arrayfun(@(x) {x.data.neural.roi} ,subject_info);
all_rois_names = unique(cat(2,all_rois_names{:}));
n_roi          = numel(all_rois_names);
n_subj         = numel(subject_info);

trial_time       = subject_info(1).data.neural.times;
save('trial_time.mat','trial_time','percentage_time'); % For easy access in SEEG_Mediation_Analysis_Response_Shape.m
if REPROCESS_DATA
    [~, tzero_index] = min(abs(trial_time)); % extract the index of the stimulus time point
    time_window      = numel(trial_time);
    
    [all_neural_data_stim_lock, all_neural_data_resp_lock, ...
        tmp_all_percent_time, tmp_all_neural_data_percent_time] = ...
        deal(cell(n_subj,n_roi));
    for ii_subj=1:n_subj
        neural    = subject_info(ii_subj).data.neural;
        subj_rois =  neural.roi;
        
        for ii_roi = 1:n_roi
            % extract the electrodes that are in the current ROI:
            subj_idx_electrodes_in_roi = find( ismember(subj_rois, all_rois_names{ii_roi}));
            
            if ~isempty(subj_idx_electrodes_in_roi)
                subj_neural_data_stim_locked = neural.f50f150(:,subj_idx_electrodes_in_roi,:);
                subj_neural_data_resp_locked = neural.f50f150_R(:,subj_idx_electrodes_in_roi,:);
                
                % Resample neural data on the percentage of
                % trial-completion time frame:
                RTs               = cell2mat(subject_info(ii_subj).data.event.response_time);
                RTs_as_time_index = tzero_index + round(RTs ./ diff(trial_time(1:2)));
                
                for ii_elec = 1:numel(subj_idx_electrodes_in_roi)
                    for trial = 1:numel(RTs)
                        full_neural_data = squeeze(subj_neural_data_stim_locked(:,ii_elec,trial));
                        full_trial_time  = trial_time;
                        % extended_idx corresponds to the index of the Resp-locked data in the Stim-locked time frame:
                        extended_idx     = RTs_as_time_index(trial)- tzero_index + (1:numel(trial_time));
                        
                        full_neural_data( extended_idx ) = subj_neural_data_resp_locked(:,ii_elec,trial);
                        full_trial_time(  extended_idx ) = RTs(trial) + trial_time;
                        
                        % upsample by 100 then downsample by RT to get a
                        % resampling of 100/RT-th
                        extended_neural = interp(full_neural_data,100);
                        extended_time   = interp(full_trial_time,100);
                        
                        step_prct       = ceil(RTs(trial)/100/diff(extended_time(1:2))); %downsampling index
                        
                        reduced_neural  = extended_neural(1:step_prct:end);
                        reduced_time    = extended_time(  1:step_prct:end);
                        
                        percent_time    = round(reduced_time/RTs(trial)*100);
                        
                        [~,first_occu_idx] = unique(percent_time); % idx used to remove the duplicate created by the rounding
                        tmp_all_percent_time{ii_subj,ii_roi}{ii_elec,trial}             = percent_time(first_occu_idx);
                        tmp_all_neural_data_percent_time{ii_subj,ii_roi}{ii_elec,trial} = reduced_neural(first_occu_idx);
                    end
                end
                
                all_neural_data_stim_lock{ii_subj,ii_roi} = subj_neural_data_stim_locked;
                all_neural_data_resp_lock{ii_subj,ii_roi} = subj_neural_data_resp_locked;
            end
        end
        disp([num2str(ii_subj),'''s data time-reformated'])
    end
    
    %% PREPROC cut completion time after RT and shortly before Stim
    n_prct_time = max_percentage_after_stim - min_percentage_before_stim + 1; % number of desired time step in the new time-format
    
    [all_percent_time, all_neural_data_percent_time] = deal(cell(n_subj,n_roi));
    for ii_subj=1:n_subj
        subj_rois =  subject_info(ii_subj).data.neural.roi;
        n_trial = numel(subject_info(ii_subj).data.event.difficulty);
        
        for ii_roi = 1:n_roi
            subj_idx_electrodes_in_roi = find( ismember(subj_rois, all_rois_names{ii_roi}));
            
            if ~isempty(subj_idx_electrodes_in_roi)
                all_neural_data_percent_time{ii_subj,ii_roi} = NaN(...
                    max_percentage_after_stim-min_percentage_before_stim+3,...
                    numel(subj_idx_electrodes_in_roi),...
                    n_trial); %+3 is a just a security margin that will soon be removed
                for trial = 1:n_trial
                    
                    for ii_elec = 1:numel(subj_idx_electrodes_in_roi)
                        idx1 = find( tmp_all_percent_time{ii_subj,ii_roi}{ii_elec,trial} >= min_percentage_before_stim, 1,'first');
                        idx2 = find( tmp_all_percent_time{ii_subj,ii_roi}{ii_elec,trial} <= max_percentage_after_stim , 1,'last');
                        tmp = tmp_all_neural_data_percent_time{ii_subj,ii_roi}{ii_elec,trial}(idx1:idx2);
                        all_neural_data_percent_time{ii_subj,ii_roi}(1:numel(tmp),ii_elec,trial) = tmp;
                    end
                end
                % remove extra nan time columns
                bad_columns_time = any(isnan(all_neural_data_percent_time{ii_subj,ii_roi}),[2,3]);
                all_neural_data_percent_time{ii_subj,ii_roi}(bad_columns_time,:,:) = [];
            end
        end
    end
    save('formated_neural_data.mat', 'all_neural_data_percent_time', 'all_neural_data_stim_lock', 'all_neural_data_resp_lock')
else
    load('formated_neural_data.mat')
end

%% PREPROC Save X & Y & failed-trial index
XYs_and_bad_index = cell(n_subj,3);
for ii_subj=1:n_subj
    trial_type   = subject_info(ii_subj).data.event.difficulty;
    trial_type   = cellfun(@(x) isequal(x,'easy'),trial_type);
    RT           = cell2mat(subject_info(ii_subj).data.event.response_time);
    trial_number = (1:numel(RT))';
    
    bad_idx               = RT>2.5;
    bad_idx               = RT>Inf;
    trial_type(bad_idx)   = [];
    RT(bad_idx)           = [];
    trial_number(bad_idx) = [];
    
    % The design matrix correct the effect of the trial-type for the trial-mean, the trial number (fatigue or signal drift).
    X = [ones(numel(RT),1),trial_number,trial_type];
    Y                            = RT;
    
    XYs_and_bad_index{ii_subj,1} = X;
    XYs_and_bad_index{ii_subj,2} = Y;
    XYs_and_bad_index{ii_subj,3} = bad_idx;
end
save('subject_XY_&_rmvd_trials.mat','XYs_and_bad_index');

%% ANALYSIS Region selection based on how they respond to either the Stimulus or Response
[~,idx_zero_trial_time]      = min(abs(trial_time));
[~,idx_zero_trial_time_prct] = min(abs(percentage_time));
[~,idx_RT_trial_time_prct]   = min(abs(percentage_time-100));

% Assess how much the signal change close to the stimuli, the response or both
time_window_bf         = 0.2; % seconds before stim or after response
time_window_after      = 1;   % seconds after stim or before response
% for the percentage time this quantity is normalized by the mean RT (across trials and subjects)
overall_mean_RT        = mean( arrayfun(@(subj) mean( cell2mat(subj.data.event.response_time) ), subject_info) );
time_window_bf_prct    = round(time_window_bf   /overall_mean_RT*100);
time_window_after_prct = round(time_window_after/overall_mean_RT*100);

neural_activity_list         = {all_neural_data_percent_time,all_neural_data_stim_lock, all_neural_data_resp_lock};
[delta_stim,delta_resp]      = deal(cell(3,1));
for ii_data=1:3
    for ii_subj=1:n_subj
        neural    = subject_info(ii_subj).data.neural;
        subj_rois =  neural.roi;
        
        n_trial = numel(subject_info(ii_subj).data.event.difficulty);
        RTs     = subject_info(ii_subj).data.event.response_time;
        
        for ii_roi = 1:n_roi
            subj_idx_electrodes_in_roi = find( ismember(subj_rois, all_rois_names{ii_roi}));
            
            if ~isempty(subj_idx_electrodes_in_roi)
                data = neural_activity_list{ii_data}{ii_subj,ii_roi};
                
                for trial = 1:size(data,3)
                    % prepare indexes determining the window of averaging
                    switch ii_data
                        case 1 % prct time
                            stim_idx = idx_zero_trial_time_prct;
                            RTidx    = idx_RT_trial_time_prct;
                            dt       = 1;
                        case {2;3} % Stim or Resp - locked
                            stim_idx  = idx_zero_trial_time;
                            [~,RTidx] = min(abs(RTs{trial}-trial_time));
                            dt        = diff(trial_time(1:2));
                    end
                    
                    % Perform the averaging. When the time-windows extended
                    % over the recorded signal it is shrinked to the first
                    % (or last) data point recorded
                    switch ii_data
                        case 1 % prct time
                            trial_data_before_stim = data(max(1,           round( stim_idx - time_window_bf_prct   /dt )):stim_idx, :,trial);
                            trial_data_after_stim  = data(stim_idx:min(end,round( stim_idx + time_window_after_prct/dt ))         , :,trial);
                            
                            trial_data_before_resp = data(max(1,           round( stim_idx - time_window_after_prct/dt)):RTidx, :,trial);
                            trial_data_after_resp  = data(stim_idx:min(end,round( stim_idx + time_window_bf_prct   /dt))      , :,trial);
                        case 2 % S-locked
                            trial_data_before_stim = data(max(1,           round( stim_idx - time_window_bf   /dt)):stim_idx, :,trial);
                            trial_data_after_stim  = data(stim_idx:min(end,round( stim_idx + time_window_after/dt))         , :,trial);
                            
                            trial_data_before_resp = NaN;
                            trial_data_after_resp  = NaN;
                        case 3 % r-locked
                            trial_data_before_stim = NaN;
                            trial_data_after_stim  = NaN;
                            
                            trial_data_before_resp = data(max(1,           round( RTidx - time_window_after/dt)):RTidx, :,trial);
                            trial_data_after_resp  = data(stim_idx:min(end,round( RTidx + time_window_bf   /dt))      , :,trial);
                    end
                    
                    delta_stim{ii_data}{ii_subj,ii_roi}(:,trial) = mean(trial_data_before_stim,1) - mean(trial_data_after_stim,1);
                    delta_resp{ii_data}{ii_subj,ii_roi}(:,trial) = mean(trial_data_before_resp,1) - mean(trial_data_after_resp,1);
                end
            end
        end
    end
end

% Extract the most responsive electrode per ROI
[delta_stim_per_roi, delta_resp_per_roi,...
    delta_stim_best_elec_idx,delta_resp_best_elec_idx] = deal(cell(n_roi,3));
for ii_data=1:3
    for ii_roi = 1:n_roi
        % compute the trial-average of signal difference for each electrodes:
        tmp_S = cellfun(@(elec_X_trial) { max( abs( mean( elec_X_trial ,2) ) )}, delta_stim{ii_data}(:,ii_roi));
        delta_stim_per_roi{ii_roi,ii_data}       = cat(1,tmp_S{:});
        
        % and extract the most responsive electrode
        delta_stim_best_elec_idx{ii_roi,ii_data} = cellfun(@(elec_X_trial) ...
            { find( abs(mean(elec_X_trial,2)) == max(abs(mean(elec_X_trial,2))) )}, delta_stim{ii_data}(:,ii_roi));
        
        % repeat the same operation for the signal difference around the response
        tmp_R = cellfun(@(elec_X_trial) {max(abs(mean(elec_X_trial,2)))}, delta_resp{ii_data}(:,ii_roi));
        delta_resp_per_roi{ii_roi,ii_data} = cat(1,tmp_R{:});
        delta_resp_best_elec_idx{ii_roi,ii_data} = cellfun(@(elec_X_trial) {find(abs(mean(elec_X_trial,2)) == max(abs(mean(elec_X_trial,2))))}, delta_resp{ii_data}(:,ii_roi));
    end
end

% Compute uncorrected p-values for ROI-responsiveness across participant
[delta_stim_pval_per_roi, delta_resp_pval_per_roi] = deal(NaN(n_roi,3));
for ii_data=1:3
    for ii_roi = 1:n_roi
        [~,delta_stim_pval_per_roi(ii_roi,ii_data)] = ...
            ttest(delta_stim_per_roi{ii_roi,ii_data},0,'Tail','right');
        [~,delta_resp_pval_per_roi(ii_roi,ii_data)] = ...
            ttest(delta_resp_per_roi{ii_roi,ii_data},0,'Tail','right');
    end
end

% FDR correction see https://en.wikipedia.org/wiki/False_discovery_rate#Benjamini%E2%80%93Hochberg_procedure
% First extract the uncorrected p-value and correct them for 1) the number
% of difference tested and 2) the unsigned test
p_diff = nanmin(delta_stim_pval_per_roi,delta_resp_pval_per_roi).*[2,1,1]*2; % * [2,1,1] for the number of diff tested, *2 for the unsigned to signed test
p_diff_FDR = p_diff*NaN;
ROI_selected = cell(3,1);
for ii_data=1:3
    % FDR correction: each p-value is scaled by the ratio of the number of
    % test and its rank among the other p-values).
    % It leads to a strong correction of the best p-values, and weak for
    % the worst one.
    [~,ps_idx]      = sort(p_diff(:,ii_data),'ascend');
    porder          = [];
    porder(ps_idx)  = 1:numel(ps_idx);
    n_non_nan_value = sum(~isnan(p_diff(:,ii_data)));
    p_diff_FDR(:,ii_data) = p_diff(:,ii_data) .* n_non_nan_value ./ porder';
    
    ROI_selected{ii_data} = find( p_diff_FDR(:,ii_data)<ROI_responsiveness_p_value_threshold);
end

disp(['Selected ROI (prct time):  ', strjoin(all_rois_names(ROI_selected{1}),' ; ')])
disp(['Selected ROI (S-lkd time): ', strjoin(all_rois_names(ROI_selected{2}),' ; ')])
disp(['Selected ROI (R_lkd time): ', strjoin(all_rois_names(ROI_selected{3}),' ; ')])

save('ROI_selection.mat','ROI_selected','all_rois_names','p_diff')

%% ANALYSIS trial level fit a of the shape model
neural_activity_list = {all_neural_data_percent_time,all_neural_data_stim_lock, all_neural_data_resp_lock};

% Compute how much electrodes will be fitted
n_trial_to_fit = 0;
for ii_data = 1:numel(neural_activity_list)
    all_neural_data = neural_activity_list{ii_data};
    for ii_roi =1:n_roi
        if ismember(ii_roi,ROI_selected{ii_data})
            for ii_subj=1:n_subj
                if ~isempty(all_neural_data{ii_subj,ii_roi})
                    all_trials_data = all_neural_data{ii_subj,ii_roi};
                    n_trial_to_fit  = n_trial_to_fit + size(all_trials_data,2)*size(all_trials_data,3);
                end
            end
        end
    end
end

% fit the shape model
[subj_roi_P,subj_roi_simu] = deal(cell(numel(neural_activity_list),n_subj,n_roi));
idx_perct   = 1; % index of the completion-time time-frame
elec_count  = 0; % electrode fitted so far
% function wrap to perform a simulation of the shape model:
simu_shape  = @(x,P) arrayfun(@(ii_t) g_shape_quick([],P,x(ii_t),struct()), 1 :numel(x));
all_neural_data_best_elec = cell(3,1);
if REFIT_ELECTRODES
    
    for ii_data = 1:numel(neural_activity_list)
        t_data = tic;
        all_neural_data = neural_activity_list{ii_data};
        
        %% Extract best electrodes data
        if 0
            for ii_roi =1:n_roi
                for ii_subj=1:n_subj
                    switch ii_data
                        case 1 % completion time
                            electrod_idx1 = delta_stim_best_elec_idx{ii_roi,ii_data}{ii_subj};
                            electrod_idx2 = delta_resp_best_elec_idx{ii_roi,ii_data}{ii_subj};
                            
                            if abs(mean(delta_stim{ii_data}{ii_subj,ii_roi}(electrod_idx1,:))) > ...
                                    abs(mean(delta_resp{ii_data}{ii_subj,ii_roi}(electrod_idx2,:)))
                                electrod_idx = electrod_idx1;
                            else
                                electrod_idx = electrod_idx2;
                            end
                            
                        case 2 % Stim-locked time
                            electrod_idx = delta_stim_best_elec_idx{ii_roi,ii_data}{ii_subj};
                            
                        case 3 % Resp-locked time
                            electrod_idx = delta_resp_best_elec_idx{ii_roi,ii_data}{ii_subj};
                    end
                    all_neural_data{ii_subj,ii_roi} = squeeze( all_neural_data{ii_subj,ii_roi}(:,electrod_idx,:) )';
                end
            end
            all_neural_data_best_elec{ii_data} = all_neural_data;
        end
        
        %% subject dynamic's fit
        % The fit is performed at four levels.
        % 1) the grand-average (across subjects, rois and trials) is fited
        % 2) the parameter of the grand-average fit are used as priors/seed for the
        % ROI-average fit (across subjects and trials in a given ROI)
        % 3) the parameter of the ROI-average fit are used as priors/seed for the
        % best-electrode-average fit (across trials)
        % 4) the parameter of the best-electrode-average fit are used as priors/seed for the
        % trial-level fit.
        % Each fit is performed on zscored to data for best performance.
        % However, The zscoring is then reversed directly on the parameter at the
        % trial-level fit.
        
        if ii_data == idx_perct, time_list = percentage_time/100;
        else,                    time_list = trial_time;
        end
        
        % Compute the grand-average:
        avg_roi_subj_trial = cellfun(@(subj_roi_data){nanmean(subj_roi_data,[2,3])'},all_neural_data);
        avg_roi_subj_trial(cellfun(@isempty,avg_roi_subj_trial)) = [];
        avg_roi_subj_trial = cat(1,avg_roi_subj_trial{:});
        avg_roi_subj_trial =  mean(avg_roi_subj_trial,1);
        
        % zscore it to ease the fit
        data = avg_roi_subj_trial; mu = mean(data); sig = std(data);
        
        % Fit the data
        avg_roi_subj_trial_P = g_shape_quick(time_list,(data-mu)/sig);
        
        % Uncomment the next line to visualise the fit
        %         figure;plot((data-mu)/sig); hold on; plot(simu_shape(time_list,avg_roi_subj_trial_P),'--');
        
        for ii_roi =ROI_selected{ii_data}'
            % Compute the ROI-average:
            avg_subj_trial = cellfun(@(subj_roi_data){nanmean(subj_roi_data,[2,3])'},all_neural_data(:,ii_roi));
            avg_subj_trial(cellfun(@isempty,avg_subj_trial)) = [];
            avg_subj_trial = cat(1,avg_subj_trial{:});
            avg_subj_trial =  mean(avg_subj_trial,1);
            
            data =avg_subj_trial; mu = mean(data); sig = std(data);
            
            % Fit the data by adjusting the parameters from the grand-average fit:
            prev_P           = avg_roi_subj_trial_P;
            avg_subj_trial_P = g_shape_quick(time_list,(data-mu)/sig,prev_P);
            
            % Uncomment the next line to visualise the fit
            %                 figure;plot((data-mu)/sig); hold on; plot(simu_shape(time_list,avg_subj_trial_P),'--');
            
            for ii_subj=1:n_subj
                if ~isempty(all_neural_data{ii_subj,ii_roi})
                    % Compute the electrode-wide-average:
                    %                     all_trials_data = squeeze(all_neural_data{ii_subj,ii_roi});
                    %                     avg_trial       = nanmean(all_trials_data);
                    avg_elec_trial = cellfun(@(subj_roi_data){nanmean(subj_roi_data,[2,3])'},all_neural_data(ii_subj,ii_roi)); avg_elec_trial =  avg_elec_trial{1};
                    
                    data =avg_elec_trial; mu = mean(data); sig = std(data);
                    
                    % Fit the data by adjusting the parameters from the ROI-average fit:
                    prev_P = avg_subj_trial_P;
                    avg_elec_trial_P = g_shape_quick(time_list,(data-mu)/sig,prev_P);
                    
                    % Uncomment the next line to visualise all fits
                    %                     figure;plot((data-mu)/sig);hold on;plot([simu_shape(time_list,avg_roi_subj_trial_P);simu_shape(time_list,avg_subj_trial_P);simu_shape(time_list,avg_elec_trial_P)]','--');legend({'Trial mean','RoixSubjxTrial avg','SubjxTrial avg','Trialxelec avg'})
                    
                    for ii_elec = 1:size(all_neural_data{ii_subj,ii_roi},2)
                        avg_trial = nanmean(all_neural_data{ii_subj,ii_roi}(:,ii_elec,:),3)';
                        
                        data =avg_trial; mu = mean(data); sig = std(data);
                        
                        % Fit the data by adjusting the parameters from the ROI-average fit:
                        prev_P = avg_elec_trial_P;
                        avg_trial_P = g_shape_quick(time_list,(data-mu)/sig,prev_P);
                        
                        % Uncomment the next line to visualise all fits
                        %                     figure;plot((data-mu)/sig); hold on; plot([simu_shape(time_list,avg_roi_subj_trial_P);simu_shape(time_list,avg_subj_trial_P);simu_shape(time_list,avg_elec_trial_P);simu_shape(time_list,avg_trial_P)]','--');                    legend({'Trialxelec mean','RoixSubjxTrial avg','SubjxTrial avg','ElecxTrial avg','Trial avg'})
                        
                        tsubj = tic;
                        for ii_trial = 1:size(all_neural_data{ii_subj,ii_roi}(:,ii_elec,:),3)
                            % Collect data for the trial-level fit
                            %                         data = all_trials_data(ii_trial,:);
                            data = all_neural_data{ii_subj,ii_roi}(:,ii_elec,ii_trial);
                            mu = mean(data);
                            sig = std(data);
                            
                            % Fit the data by adjusting the parameters from the best-electrode-average fit:
                            P = g_shape_quick(time_list,(data-mu)/sig,avg_trial_P);
                            
                            % Uncomment the next line to visualise all fits
                            %                         figure;plot((data-mu)/sig); hold on; plot([simu_shape(time_list,avg_roi_subj_trial_P);simu_shape(time_list,avg_subj_trial_P);simu_shape(time_list,avg_elec_trial_P);simu_shape(time_list,avg_trial_P);simu_shape(time_list,P)]','--','LineWidth',2);legend({'Trial','RoixSubjxTrial avg','SubjxTrial avg','SubjxTrial avg','Trial avg','Trial'})
                            
                            % Reverse the zscoring by adjusting the amplitude parameters:
                            P(2) = P(2)*sig + mu; % A: the peak amplitude
                            P(5) = P(5)*sig + mu; % b1: the starting baseline
                            P(8) = P(8)*sig + mu; % b2: the finishing baseline
                            
                            simu = simu_shape(time_list,P);
                            
                            % Uncomment the next line to visualise the un-zscored fit
                            %                     figure;  hold on; plot(data);plot(simu);
                            
                            subj_roi_P{ii_data,ii_subj,ii_roi}{ii_elec,ii_trial}    = P;
                            subj_roi_simu{ii_data,ii_subj,ii_roi}{ii_elec,ii_trial} = simu;
                            
                            elec_count =  elec_count +1;
                        end
                        % Uncomment the next line to visualize the mean of fitted data
                        %                     figure;plot(mean(all_neural_data{ii_subj,ii_roi},[2,3])'); hold on; plot(mean(cat(1,subj_roi_simu{ii_data,ii_subj,ii_roi}{:})),'--')
                        
                    end
                    t_to_process_subj= toc(tsubj);
                    disp([num2str(t_to_process_subj),' for the electrode ',num2str(elec_count),'/',num2str(n_trial_to_fit)]);
                end
                disp([ ' Data type n°',num2str(ii_data),'/',num2str(numel(neural_activity_list)), ...
                    ' ROI n°',num2str(ii_roi),'/',num2str(n_roi), ...
                    ' Subj n°',num2str(ii_subj),'/',num2str(n_subj)])
            end
        end
        t_to_process_data= toc(tsubj);
    end
    save('formated_neural_data_best_elec.mat','all_neural_data_best_elec')
    save('paramshapeAllElec.mat','subj_roi_P','subj_roi_simu')
else
    previous_fit  = load('paramshape.mat');
    subj_roi_P    = previous_fit.subj_roi_P;
    subj_roi_simu = previous_fit.previous_fit;
end

%% EOF
disp(['EOF: ', script_name])

