function [stats_data] = MrM_stats_erp_baseline_long(cfg, data)

% stats_erp_baseline performs statistics against the baseline period.
% Data should be the output of ft_preprocessing (and all trials the same lenght)
%
% Use as: stats_erp_baseline(cfg,data)
% The configuration can have the following parameters
% cfg.trialinfo     = double refering to the condition defined in trialinfo
%                     if this is not specified all trials will be used.
% cfg.channel       = cells containing selected channel (or 'all')
% cfg.time          = time period over which the stats will be performed
%                     [sec sec] or 'all'
% cfg.erp_baseline  = time period from onset used as baseline
%                     [sec sec]
% cfg.iteration     = number of random permutation
% cfg.erp_prob      = statistic threshold for power (two tails like [2.5 97.5])
% cfg.distribution  = 'yes' save the random distribution
% cfg.mode          = define how the baseline "value" is defined and how the test is done :
%
%                     'random within and across (distribution)' **
%                     picked different random baseline time points for each trials
%                     then average across trials
%                     and looped to build surrogate "baseline" distribution
%
%                     'random within and across (paired)'
%                     picked different random baseline time points for each trials
%                     then permute randomly with post onset time point
%                     compute the average of the difference 
%                     and looped to build surrogate "paired-permutation" distribution
%
%                     'random across (distribution)'
%                     picked same random baseline time point for each trial
%                     then average across trials
%                     and loop to build surrogate "baseline" distribution
%
%                     'random across (paired)' ***
%                     picked same random baseline time point for each trial
%                     then permute randomly with post onset time point
%                     compute the average of the difference
%                     and looped to build the surrogate "paired-permutation" ditribution
%
%                     'mean'
%                     mean baseline permuted randomly with each post baseline time point
%                     looped to build the surrogate "paired-permutation" distribution
%                     based on the average of the difference
%
% cfg.z_score       = 'yes'
%                     compute z-statistic in addition to the mean
% cfg.z_correction  = 'yes'
%                     compute correction for multiple comparison based on 'z-score max'

%channels selection
if strcmp(cfg.channel, 'all')
    id_ch = 1:length(data.label);
elseif iscell(cfg.channel)
    for c=1:length(cfg.channel)
        id_ch = strmatch(cfg.channel(c), data.label, 'exact');
    end
elseif ~isfield (cfg, 'channel')
    error('no channel was specified');
else
    error('channel field wrongly specified');
end

%time index
switch isfield(cfg, 'time')
    case 1
        if length(cfg.time) == 2
            id_time =[find(round(data.time{1}*1000) == round(cfg.time(1)*1000)) find(round(data.time{1}*1000) == round(cfg.time(2)*1000))];
        elseif strcmp(cfg.time, 'all');
            id_time = [data.time{1}(1) data.time{1}(end)];
        else
            error('time period wrongly stated');
        end
    case 0
        error('no time of interest defined');
end

%p values
switch isfield(cfg, 'erp_prob')
    case 1
        if length(cfg.erp_prob)==2
        elseif length(cfg.erp_prob)==1
            error('p threshold for power should contains 2 values (two-tail)');
        else
            error('p threshold for power wrongly stated');
        end
    case 0
        error('no p threshold defined for power');
end

%Prepare the output
stats_data.cfg       = cfg;

%data selection
stats_data.time  = data.time{1}(id_time(1):id_time(2));
stats_data.dimord        = 'chan_time';
stats_data.label = data.label(id_ch);
%select the trials
if isfield(cfg, 'trialinfo')
    trialS       = find(data.trialinfo(:,1)==cfg.trialinfo);
else
    trialS       = 1:size(data.trial);
end
%select the data and compute mean and standard error
stats_data.data  = reshape(cell2mat(data.trial), size(data.trial{1}, 1), size(data.trial{1}, 2), size(data.trial,2));
stats_data.data  = stats_data.data(id_ch,id_time(1):id_time(2),trialS);
stats_data.erp   = squeeze(mean(stats_data.data,3));
stats_data.ste   = std(stats_data.data,0,3)/sqrt(size(stats_data.data,3));

%Baseline index
switch isfield(cfg, 'erp_baseline')
    case 1
        if isnumeric(cfg.erp_baseline) && length(cfg.erp_baseline) ==2
            if cfg.erp_baseline(1)>= cfg.time(1) && cfg.erp_baseline(2)<= cfg.time(2)
                id_Bl = [find(round(stats_data.time*1000) == round(cfg.erp_baseline(1)*1000)) find(round(stats_data.time*1000) == round(cfg.erp_baseline(2)*1000))];
            else
                error('baseline period should be included into the selected time period');
            end
        else
            error('baseline wrongly stated');
        end
    case 0
        error('no baseline period defined');
end

clear data;
%dum variables
nb_ch            = size(stats_data.data,1);
nb_tf            = size(stats_data.data,2);
nb_trial         = size(stats_data.data,3);

% start the permutation process
if any(findstr(cfg.mode, 'distribution)'))
    Rd_mean = nan(nb_ch,nb_tf,cfg.iteration);
    Rd_std  = nan(nb_ch,nb_tf,cfg.iteration);
    for lp = 1:cfg.iteration;
        if strcmp(cfg.mode, 'random within and across (distribution)')
            % different random time points averaged across trial
            % (different for each trial, random across loops)

            % build the random draws index matrix
            pick_Rd_Index = round(rand(nb_trial, nb_tf) * (id_Bl(2) - id_Bl(1)) + id_Bl(1));
            % build the random draws value matrix
            for tf = 1:nb_tf
                for tr = 1:nb_trial
                    pick_Rd_Bl(:,tf, tr) = stats_data.data(:,pick_Rd_Index(tr,tf),tr);
                end
            end

        elseif strcmp(cfg.mode, 'random across (distribution)')
            % same random time point averaged across trial
            % (identical for each trial, random across loops)

            % build the random draws index matrix
            pick_Rd_Index = round(rand(nb_tf,1) * (id_Bl(2) - id_Bl(1)) + id_Bl(1));
            % build the random draws value matrix
            for tf = 1:nb_tf
                pick_Rd_Bl(:,tf,:) = stats_data.data(:,pick_Rd_Index(tf),:);
            end
        end

        %compute the average Baseline value (across trial)
        Rd_mean(:,:, lp) = mean(pick_Rd_Bl, 3);
        Rd_std(:,:,lp)   = std(pick_Rd_Bl,0,3);
        disp(['Randomisation #' num2str(lp) ' over ' num2str(cfg.iteration)]);
    end;

    %define percentile for each time point, each electrode
    thrshld_erp = prctile(Rd_mean, cfg.erp_prob , 3);
    %create stats mask
    stats_data.erp_mask = (stats_data.erp <= thrshld_erp(:,:,1)) | (stats_data.erp >= thrshld_erp(:,:,2));
    %compute probabilities
    stats_data.erp_p = nan(nb_ch,nb_tf);
    for ch = 1:nb_ch
        for tf = 1:nb_tf
            %find on which side of the subrogate disctribution fell the mean
            if stats_data.erp(ch,tf) < mean(Rd_mean(ch,tf,:),3);
                stats_data.erp_p(ch,tf) = length(find(Rd_mean(ch,tf,:)<stats_data.erp(ch,tf)))/size(Rd_mean,3);
            elseif stats_data.erp(ch,tf) > mean(Rd_mean(ch,tf,:),3);
                stats_data.erp_p(ch,tf) = length(find(Rd_mean(ch,tf,:)>stats_data.erp(ch,tf)))/size(Rd_mean,3);
            end
        end
    end

elseif any(strfind(cfg.mode, 'paired)'))
    for lp = 1:cfg.iteration;
        if strcmp(cfg.mode, 'random within and across (paired)')
            % different random baseline time points across trials
            % are paired with postonset time points for permutation
            % (same for each time point, different for each trial and freq random across loops ?)

            % build the random draws index matrix
            pick_Rd_Index = round(rand(nb_trial,nb_tf) * (id_Bl(2) - id_Bl(1)) + id_Bl(1));
            for tf = 1:nb_tf
                for tr = 1:nb_trial
                    pick_Rd_Bl(:,tf, tr) = stats_data.data(:,pick_Rd_Index(tr, tf),tr);
                end
            end
            
        elseif strcmp(cfg.mode, 'random across (paired)')
            % same random baseline time point across trials
            % are paired with time points for permutation
            % (identical for each trial, random across loops)

            % build the random draws index matrix
            pick_Rd_Index = round(rand(nb_tf,1) * (id_Bl(2) - id_Bl(1)) + id_Bl(1));
            % make the random draws value matrix
            for tf = 1:nb_tf
                pick_Rd_Bl(:,tf,:) = stats_data.data(:,pick_Rd_Index(tf), :);
            end
        end

        %compute the observed average of the differences
        stats_data.Bl_meandiff(:,:,lp) = mean(stats_data.data - pick_Rd_Bl, 3);
        stats_data.Bl_stddiff(:,:,lp)  = std(stats_data.data - pick_Rd_Bl, 0, 3);

        % build the mixing paired matrices
        Rdperm_tr               = [1:nb_trial]'.*round(rand(nb_trial,1));
        Rdperm_tr(Rdperm_tr==0) = [];
        erp1_Rd                 = stats_data.data;
        erp2_Rd                 = pick_Rd_Bl;
        erp1_Rd(:,:,Rdperm_tr)  = pick_Rd_Bl(:,:, Rdperm_tr);
        erp2_Rd(:,:,Rdperm_tr)  = stats_data.data(:, :, Rdperm_tr);
        %compute the average of the differences
        diff_Rd = erp1_Rd - erp2_Rd;
        Rd_mean(:,:,lp)  = mean(diff_Rd,3);
        Rd_std(:,:,lp)   = std(diff_Rd,0,3);
        disp(['Randomisation #' num2str(lp) ' over ' num2str(cfg.iteration)]);
    end
        %define percentile for each time point, each electrode
        thrshld_erp = prctile(Rd_mean, cfg.erp_prob , 3);
        %create stats mask
        stats_data.erp_mask = (mean(stats_data.Bl_meandiff,3) <= thrshld_erp(:,:,1) | mean(stats_data.Bl_meandiff, 3) >= thrshld_erp(:,:,2));

        %compute the probabilities
        stats_data.erp_p = nan(nb_ch,nb_tf);
        for ch = 1:nb_ch
            for tf = 1:nb_tf
                %define the tail depending on mean side for the amplitude
                if mean(stats_data.Bl_meandiff(ch,tf,:),3) < mean(Rd_mean(ch,tf,:),3);
                    stats_data.erp_p(ch,tf) = length(find(Rd_mean(ch,tf,:)<mean(stats_data.Bl_meandiff(ch,tf,:),3)))/size(Rd_mean,3);
                elseif mean(stats_data.Bl_meandiff(ch,tf,:),3) > mean(Rd_mean(ch,tf,:),3);
                    stats_data.erp_p(ch,tf) = length(find(Rd_mean(ch,tf,:)>mean(stats_data.Bl_meandiff(ch,tf,:),3)))/size(Rd_mean,3);
                end
            end
        end
        
elseif strcmp(cfg.mode, 'mean')
%paired-test
%permutations between the mean baseline and each time point

    % make the Baseline data matrix
    Bl_mean = repmat(mean(stats_data.data(:,id_Bl(1):id_Bl(2),:),2),[1 nb_tf 1]);
    %compute the observed average of the differences
    stats_data.Bl_meandiff = stats_data.data - Bl_mean;
    %start the randomisation process
    for lp = 1:cfg.iteration;
        % build the mixing paired matrices
        Rdperm_tr               = [1:nb_trial]'.*round(rand(nb_trial,1));
        Rdperm_tr(Rdperm_tr==0) = [];
        erp1_Rd                 = stats_data.data;
        erp2_Rd                 = Bl_mean;
        erp1_Rd(:,:,Rdperm_tr)  = Bl_mean(:, :, Rdperm_tr);
        erp2_Rd(:,:,Rdperm_tr)  = stats_data.data(:, :, Rdperm_tr);
        %compute the average of the difference
        diff_Rd          = erp1_Rd - erp2_Rd;
        Rd_mean(:,:,lp)  = mean(diff_Rd,3);
        Rd_std(:,:,lp)   = std(diff_Rd,0,3);
        disp(['Randomisation #' num2str(lp) ' over ' num2str(cfg.iteration)]);
    end

    %define percentiles for each time point, each electrode
    thrshld_erp = prctile(Rd_mean, cfg.erp_prob , 3);
    %create stats mask
    stats_data.erp_mask = (squeeze(mean(stats_data.Bl_meandiff,3)) <= thrshld_erp(:,:,1) | squeeze(mean(stats_data.Bl_meandiff,3)) >= thrshld_erp(:,:,2));
    %compute probabilities
    stats_data.erp_p = nan(nb_ch,nb_tf);
    for ch = 1:nb_ch
        for tf = 1:nb_tf
            %define the tail depending on mean side for the amplitude
            if mean(stats_data.Bl_meandiff(ch,tf,:),3) < mean(Rd_mean(ch,tf,:),3);
                stats_data.erp_p(ch,tf) = length(find(Rd_mean(ch,tf,:)<mean(stats_data.Bl_meandiff(ch,tf,:),3)))/size(Rd_mean,3);
            elseif stats_data.Bl_meandiff(ch,tf,:) > mean(Rd_mean(ch,tf,:),3);
                stats_data.erp_p(ch,tf) = length(find(Rd_mean(ch,tf,:)>mean(stats_data.Bl_meandiff(ch,tf,:),3)))/size(Rd_mean,3);
            end
        end
    end
else
    error('you have to specify the baseline mode')
end

%________________________________________________
%compute the z-score and corresponding statistics
if strcmp(cfg.z_score, 'yes')
    if any(strfind(cfg.mode, 'distribution)'))
        Rd_mean_pop(:,:,1)                    = stats_data.erp;
        Rd_mean_pop(:,:,2:size(Rd_mean,3)+1)  = Rd_mean;
        Rd_mean_z(:,:,1)      = (mean(stats_data.data,3) - mean(Rd_mean_pop,3))./std(stats_data.data,0,3);
        for t=1:size(Rd_mean,3)
             Rd_mean_z(:,:,t+1) = (Rd_mean(:,:,t) - mean(Rd_mean_pop,3))./Rd_std(:,:,t);
        end

    elseif any(strfind(cfg.mode, 'paired)')) || strcmp(cfg.mode, 'mean')
        Rd_mean_pop(:,:,1)                   = mean(stats_data.Bl_meandiff,3);
        Rd_mean_pop(:,:,2:size(Rd_mean,3)+1) = Rd_mean;
        Rd_mean_z(:,:,1)      = (mean(stats_data.Bl_meandiff,3) - mean(Rd_mean_pop,3))./std(stats_data.Bl_meandiff,0,3);
        for t=1:size(Rd_mean,3)
            Rd_mean_z(:,:,t+1)  = (Rd_mean(:,:,t) - mean(Rd_mean_pop,3))./Rd_std(:,:,t);
        end
    end
    %define percentiles for each time point, each electrode
    thrshld_erp_z = prctile(Rd_mean_z, cfg.erp_prob , 3);
    %create stats mask
    stats_data.erp_mask_z = (squeeze(Rd_mean_z(:,:,1)) <= thrshld_erp_z(:,:,1) | squeeze(Rd_mean_z(:,:,1)) >= thrshld_erp_z(:,:,2));
    %compute probabilities
    stats_data.erp_p_z = nan(nb_ch,nb_tf);
    for ch = 1:nb_ch
        for tf = 1:nb_tf
            %define the tail depending on mean side for the amplitude
            if Rd_mean_z(ch,tf,1) < mean(Rd_mean_z(ch,tf,:),3);
                stats_data.erp_p_z(ch,tf) = length(find(Rd_mean_z(ch,tf,:)<Rd_mean_z(ch,tf,1)))/size(Rd_mean_z,3);
            elseif Rd_mean_z(ch,tf,1) > mean(Rd_mean_z(ch,tf,:),3);
                stats_data.erp_p_z(ch,tf) = length(find(Rd_mean_z(ch,tf,:)>Rd_mean_z(ch,tf,1)))/size(Rd_mean_z,3);
            end
        end
    end
else
end

%_________________________________________________
%z-score max method for multicomparison correction
if strcmp(cfg.z_correction, 'yes')
    % correction made over time dimension
     [Z_max Id_max] = max(abs(Rd_mean_z),[],2);
     Rd_mean_z_max = Rd_mean_z(Id_max);
    %define percentile for each electrode
    thrshld_erp_z_corr    = prctile(Rd_mean_z_max, cfg.erp_prob , 3);
    %create stat mask
    stats_data.erp_mask_z_corr =  (Rd_mean_z(:,:,1) <= repmat(thrshld_erp_z_corr(:,:,1),[1 nb_tf 1])) ...
                                | (Rd_mean_z(:,:,1) >= repmat(thrshld_erp_z_corr(:,:,2),[1 nb_tf 1]));

    %compute probabilities
    stats_data.erp_p_z_corr    = nan(nb_ch,nb_tf);
    for ch = 1:nb_ch
        for tf = 1:nb_tf
            %define the tail depending on mean side for the amplitude
            if Rd_mean_z(ch,tf,1) < mean(Rd_mean_z_max(ch,1,:),3);
                stats_data.erp_p_z_corr(ch,tf) = length(find(Rd_mean_z_max(ch,1,:)<Rd_mean_z(ch,tf,1)))/size(Rd_mean_z_max,3);
            elseif Rd_mean_z(ch,tf,1) > mean(Rd_mean_z_max(ch,1,:),3);
                stats_data.erp_p_z_corr(ch,tf) = length(find(Rd_mean_z_max(ch,1,:)>Rd_mean_z(ch,tf,1)))/size(Rd_mean_z_max,3);
            end
        end
    end
end

%________________________
%save distribution or not
switch isfield(cfg, 'distribution')
    case 1
        if strcmp(cfg.distribution, 'yes');
            stats_data.distribution = Rd_mean;
            if strcmp(cfg.z_score, 'yes')
               stats_data.distribution_z = Rd_mean_z;
                if strcmp(cfg.z_correction, 'yes')
                stats_data.distribution_z_max = Rd_mean_z_max;    
                end                
            end
            disp('random distributions were saved');
        elseif strcmp(cfg.distribution, 'no');
            disp('random distributions were not saved');
        end
    case 0
        disp('random distributions were not saved');
end

stats_data.cfg = cfg;
disp('Permutations finished!');
disp(['Channel tested: ' num2str(id_ch)]);
disp(['over ' num2str(cfg.time(1)) ' to ' num2str(cfg.time(2)) ' time period (sec)']);
disp(['with baseline being from ' num2str(cfg.erp_baseline(1)) ' to ' num2str(cfg.erp_baseline(2)) ' ms']);
disp(['p thresholds were ' num2str(cfg.erp_prob(1)) ' & ' num2str(cfg.erp_prob(2))]);

%% PLOT DISTRIBUTION
% figure; hold on
% period = [ ];
% Dt = find(stats_data.time == period(:));
%
% data_plot = squeeze(mean(diff_Rd(Dt(:), :),1) ,2));
% steps_edges = linspace(min(data_plot), max(data_plot), 101);
% bar(steps_edges, histc(data_plot, steps_edges), 'k');
% dat=squeeze(mean(stats_data.erp(1, Dt(:)),1), 2));
% line([dat dat], [0 max(histc(data_plot, steps_edges))], 'Color', 'r');


end