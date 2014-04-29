function [stats_data] = MrM_stats_erp_baseline(cfg, data)

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
% cfg.erp_baseline  = selected baseline time period from onset
%                     [sec sec]
% cfg.iteration     = number of random permutation
% cfg.erp_prob      = statistic threshold for power (two tails like [2.5 97.5])
% cfg.distribution  = 'yes' save the random distribution
% cfg.mode          = 'average' or 'random' define how the baseline "value" is defined,
%                      either being the period average or a random time point drawn within it 
% cfg.erp_baseline  = 'no' or time period used to baseline correct the erp data only
%                     [sec sec] or 'all'
%                     the single trials are not baseline corrected
%                     (since statisctics are about the baseline!)

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

%erp baseline corrected
if strcmp(cfg.erp_baseline, 'no')
elseif isnumeric(cfg.erp_baseline) && length(cfg.erp_baseline) ==2
        if cfg.erp_baseline(1)>= cfg.time(1) && cfg.erp_baseline(2)<= cfg.time(2)
        id_erp_Bl         = [find(round(stats_data.time*1000) == round(cfg.erp_baseline(1)*1000)) find(round(stats_data.time*1000) == round(cfg.erp_baseline(2)*1000))];
        stats_data.erp_Bl = stats_data.erp - repmat(mean(stats_data.erp(:,id_erp_Bl(1):id_erp_Bl(2)),2), [1 size(stats_data.erp, 2)]);
        else
            error('baseline period should be included into the selected time period');
        end
elseif strcmp(cfg.erp_baseline, 'all')
    id_erp_Bl         = [stats_data.time(1) stats_data.time(end)];
    stats_data.erp_Bl = stats_data.erp - repmat(mean(stats_data.erp(:,id_erp_Bl(1):id_erp_Bl(2)),2), [1 size(stats_data.erp, 2)]);
else
    error('erp baseline not specified correctly');
end;

clear data;
%dum variables
nb_ch            = size(stats_data.data,1);
nb_tf            = size(stats_data.data,2);
nb_trial         = size(stats_data.data,3);

%start the permutation process
if strcmp(cfg.mode, 'random')
% random time points averaged across trial
% (different for each trial, random across loops)
    for lp = 1:cfg.iteration;
    % build the random draws index matrix
    pick_Rd_Index = round(rand(nb_trial,1) * (id_Bl(2) - id_Bl(1)) + id_Bl(1));
    % build the random draws value matrix
    for tr = 1:nb_trial
        pick_Rd_Bl(:,:, tr) = stats_data.data(:,repmat(pick_Rd_Index(tr),[1 nb_tf 1]),tr);
    end
    %compute the average Baseline value (across trial)
    Rd_mean(:,:, lp) = squeeze(mean(pick_Rd_Bl, 3));
    disp(['Randomisation #' num2str(lp) ' over ' num2str(cfg.iteration)]);
end;

% random time point averaged acorss trial
% (identical for each trial, random across loops)
% elseif strcmp(cfg.mode, 'random')
%     for lp = 1:cfg.iteration;
%         % build the random draws index matrix
%         pick_Rd_Index = round(rand(1,1) * (id_Bl(2) - id_Bl(1)) + id_Bl(1));
%         % build the random draws value matrix
%         for tr = 1:nb_trial
%             pick_Rd_Bl(:,:, tr) = stats_data.data(:,repmat(pick_Rd_Index,[1 nb_tf 1]),tr);
%         end
%         %compute the average Baseline value (across trial)
%         Rd_mean(:,:, lp) = squeeze(mean(pick_Rd_Bl, 3));
%         disp(['Randomisation #' num2str(lp) ' over ' num2str(cfg.iteration)]);
%     end;

elseif strcmp(cfg.mode, 'average')
    Bl_mean = repmat(mean(stats_data.data(:,id_Bl(1):id_Bl(2),:),2), [1 nb_tf 1]);
    % build the mixing matrix
    Rd_pool = cat(3, stats_data.data, Bl_mean);
    
    %start the randomisation process
for lp = 1:cfg.iteration;
    %shuffle the trials
    Rd_pool = Rd_pool(:,:, randperm(size(Rd_pool,3)));
    %compute the average
    erp1_Rd = squeeze(mean(Rd_pool(:,:,1:nb_trial),3));
    erp2_Rd = squeeze(mean(Rd_pool(:,:,1+nb_trial:end),3));
    %compute the difference
    Rd_mean(:,:,lp)=erp1_Rd-erp2_Rd;
    disp(['Randomisation #' num2str(lp) ' over ' num2str(cfg.iteration)]);
end
else
        error('you have to specify the baseline mode')
end

%define percentile for each time point and frequencies
thrshld_erp = prctile(Rd_mean, cfg.erp_prob , 3);

%create stat mask
stats_data.erp_mask = (squeeze(stats_data.erp) <= thrshld_erp(:,:,1) | squeeze(stats_data.erp) >= thrshld_erp(:,:,2));

%calculate the probabilities
stats_data.erp_p = nan(nb_ch,nb_tf);
for ch = 1:nb_ch
    for tf = 1:nb_tf
        %define the tail depending on mean side for the amplitude
        if stats_data.erp(ch,tf) < mean(Rd_mean(ch,tf,:),3);
            stats_data.erp_p(ch,tf) = length(find(Rd_mean(ch,tf,:)<stats_data.erp(ch,tf)))/size(Rd_mean,3);
        elseif stats_data.erp(ch,tf) > mean(Rd_mean(ch,tf,:),3);
            stats_data.erp_p(ch,tf) = length(find(Rd_mean(ch,tf,:)>stats_data.erp(ch,tf)))/size(Rd_mean,3);
        end
    end
end
switch isfield(cfg, 'distribution')
    case 1
        if strcmp(cfg.distribution, 'yes');
            stats_data.distribution = Rd_mean;
            disp('random distributions were saved');
        elseif strcmp(cfg.distribution, 'no');
            disp('random distributions were not saved');
        end
    case 0
        disp('random distributions were not saved');
end
bug
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