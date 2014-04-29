function [stats_data] = MrM_stats_erp_unpaired_sum(cfg, data)

% stats_erp_unpaired_sum performs statistics in a additive model context.
% (like in multisensory integration: AV vs. (A+V)
% Data should be the output of ft_preprocessing (and all trials the same lenght)
%
% Use as: stats_erp_unpaired_sum(cfg,data)
% The configuration can have the following parameters
% cfg.trialinfo_dat0 = double refering to the condition defined in trialinfo
%                      that will be compared to the sum of the two others
% cfg.trialinfo_dat1 = double refering to the condition defined in trialinfo
%                      that will be summed
% cfg.trialinfo_dat2 = double refering to the condition defined in trialinfo
%                      that will be summed
% cfg.channel        = cells containing selected channel (or 'all')
% cfg.time           = time period over which the stats will be performed
%                     [sec sec] or 'all'
% cfg.iteration      = number of random permutation
% cfg.erp_prob       = statistic threshold for power (two tails like [2.5 97.5])
% cfg.distribution   = 'yes' save the random distribution
% cfg.baseline       = 'no' or time period used to baseline correct the data
%                     [sec sec] or 'all'
% cfg.method         = 'surrogate sum' or 'paired'
%                      distribution based on sums only or on differences

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
stats_data.cfg     = cfg;
stats_data.data0   = cfg.trialinfo_dat0;
stats_data.data1   = cfg.trialinfo_dat1;
stats_data.data2   = cfg.trialinfo_dat2;
%data selection
stats_data.time    = data.time{1}(id_time(1):id_time(2));
stats_data.dimord  = 'chan_time';
stats_data.label   = data.label(id_ch);

%baseline correction index
if strcmp(cfg.baseline, 'no')
elseif isnumeric(cfg.baseline) && length(cfg.baseline) ==2
    if cfg.baseline(1)>= cfg.time(1) && cfg.baseline(2)<= cfg.time(2)
        id_Bl = [find(round(stats_data.time*1000) == round(cfg.baseline(1)*1000)) find(round(stats_data.time*1000) == round(cfg.baseline(2)*1000))];
    else
        error('baseline period should be included into the selected time period');
    end
elseif strcmp(cfg.baseline, 'all')
    id_Bl = [stats_data.time(1) stats_data.time(end)];
else
    error('baseline not specified correctly');
end;

%select the trials corresponding to each condition
if isfield(cfg, 'trialinfo_dat0') && isfield(cfg, 'trialinfo_dat1') && isfield(cfg, 'trialinfo_dat2')
    trialS_dat0         = find(data.trialinfo(:,1)==cfg.trialinfo_dat0);
    trialS_dat1         = find(data.trialinfo(:,1)==cfg.trialinfo_dat1);
    trialS_dat2         = find(data.trialinfo(:,1)==cfg.trialinfo_dat2);
else
    error('conditions must be specified using trial infos from the preprocessed data');
end
%select the data and compute mean, standard error and baseline correction
stats_data.data0    = reshape(cell2mat(data.trial), size(data.trial{1}, 1), size(data.trial{1}, 2), size(data.trial,2));
stats_data.data0    = stats_data.data0(id_ch,id_time(1):id_time(2),trialS_dat0);
stats_data.data1    = reshape(cell2mat(data.trial), size(data.trial{1}, 1), size(data.trial{1}, 2), size(data.trial,2));
stats_data.data1    = stats_data.data1(id_ch,id_time(1):id_time(2),trialS_dat1);
stats_data.data2    = reshape(cell2mat(data.trial), size(data.trial{1}, 1), size(data.trial{1}, 2), size(data.trial,2));
stats_data.data2    = stats_data.data2(id_ch,id_time(1):id_time(2),trialS_dat2);
stats_data.erp0     = squeeze(mean(stats_data.data0,3));
stats_data.erp1     = squeeze(mean(stats_data.data1,3));
stats_data.erp2     = squeeze(mean(stats_data.data2,3));
stats_data.ste0     = std(stats_data.data0,0,3)/sqrt(size(stats_data.data0,3));
stats_data.ste1     = std(stats_data.data1,0,3)/sqrt(size(stats_data.data1,3));
stats_data.ste2     = std(stats_data.data2,0,3)/sqrt(size(stats_data.data2,3));
stats_data.sum      = stats_data.erp1 + stats_data.erp2;
stats_data.sum_ste  = stats_data.ste1 + stats_data.ste2;
if ~strcmp(cfg.baseline, 'no')
    stats_data.erp0_Bl  = stats_data.erp0 - repmat(mean(stats_data.erp0(:,id_Bl(1):id_Bl(2)),2), [1 size(stats_data.erp0, 2)]);
    stats_data.erp1_Bl  = stats_data.erp1 - repmat(mean(stats_data.erp1(:,id_Bl(1):id_Bl(2)),2), [1 size(stats_data.erp1, 2)]);
    stats_data.erp2_Bl  = stats_data.erp2 - repmat(mean(stats_data.erp2(:,id_Bl(1):id_Bl(2)),2), [1 size(stats_data.erp2, 2)]);
    stats_data.sum_Bl   = stats_data.erp1_Bl + stats_data.erp2_Bl;
    stats_data.diff     = stats_data.erp0_Bl - stats_data.sum_Bl;
end;

clear data;
%dum variables
nb_ch               = length(id_ch);
nb_tf               = length(stats_data.time);
nb_trial_dat0       = size(stats_data.data0,3);
nb_trial_dat1       = size(stats_data.data1,3);
nb_trial_dat2       = size(stats_data.data2,3);

%define the combination matrix of trials sum
Rd = unidrnd(nb_trial_dat1*nb_trial_dat2, nb_trial_dat0, cfg.iteration);
%find the respectives indices
[Rd_1, Rd_2] = ind2sub([nb_trial_dat1, nb_trial_dat2],Rd);

%start the randomisation process
for lp = 1:cfg.iteration;
    for nb_0 = 1:nb_trial_dat0
        Rd_sum(:,:,nb_0) = stats_data.data1(:,:,Rd_1(nb_0,lp)) + stats_data.data2(:,:,Rd_2(nb_0,lp));
    end
    %average across the summed trials
    erp_Rd_sum(:,:,lp)    = squeeze(mean(Rd_sum, 3));
    %substract the baseline3
    if ~strcmp(cfg.baseline, 'no')
        erp_Rd_sum_Bl(:,:,lp) = erp_Rd_sum(:,:,lp) - repmat(mean(erp_Rd_sum(:,id_Bl(1):id_Bl(2),lp),2), [1 size(erp_Rd_sum,2)]);
    end
    disp(['Randomisation #' num2str(lp) ' over ' num2str(cfg.iteration)]);
end

%define percentile for each time point and frequencies and create stat mask
if ~strcmp(cfg.baseline, 'no')
    thrshld_erp = prctile(erp_Rd_sum_Bl, cfg.erp_prob , 3);
    stats_data.erp_mask = (squeeze(stats_data.erp0_Bl) <= thrshld_erp(:,:,1) | squeeze(stats_data.erp0_Bl) >= thrshld_erp(:,:,2));
else
    thrshld_erp = prctile(erp_Rd_sum, cfg.erp_prob , 3);
    stats_data.erp_mask = (squeeze(stats_data.erp0) <= thrshld_erp(:,:,1) | squeeze(stats_data.erp0) >= thrshld_erp(:,:,2));
end

%calculate the probabilities
stats_data.erp_p = nan(nb_ch,nb_tf);
for ch = 1:nb_ch
    for tf = 1:nb_tf
        %define the tail depending on mean side for the amplitude
        if ~strcmp(cfg.baseline, 'no')
            if stats_data.erp0_Bl(ch,tf) < mean(erp_Rd_sum_Bl(ch,tf,:),3);
                stats_data.erp_p(ch,tf) = length(find(erp_Rd_sum_Bl(ch,tf,:)<stats_data.erp0_Bl(ch,tf)))/size(erp_Rd_sum_Bl,3);
            elseif stats_data.erp0_Bl(ch,tf) > mean(erp_Rd_sum_Bl(ch,tf,:),3);
                stats_data.erp_p(ch,tf) = length(find(erp_Rd_sum_Bl(ch,tf,:)>stats_data.erp0_Bl(ch,tf)))/size(erp_Rd_sum_Bl,3);
            else
                stats_data.erp_p(ch,tf) = NaN;   
            end
        else
            if stats_data.erp0(ch,tf) < mean(erp_Rd_sum(ch,tf,:),3);
                stats_data.erp_p(ch,tf) = length(find(erp_Rd_sum(ch,tf,:)<stats_data.erp0(ch,tf)))/size(erp_Rd_sum,3);
            elseif stats_data.erp0(ch,tf) > mean(erp_Rd_sum(ch,tf,:),3);
                stats_data.erp_p(ch,tf) = length(find(erp_Rd_sum(ch,tf,:)>stats_data.erp0(ch,tf)))/size(erp_Rd_sum,3);
            else
            stats_data.erp_p(ch,tf) = NaN;    
            end
        end
    end
end
%eventually save the distribution
switch isfield(cfg, 'distribution')
    case 1
        if strcmp(cfg.distribution, 'yes');
            stats_data.distribution = erp_Rd_sum;
            disp('random distributions were saved');
        elseif strcmp(cfg.distribution, 'no');
            disp('random distributions were not saved');
        end
    case 0
        disp('random distributions were not saved');
end

disp('Permutations finished!');
disp(['Channel tested: ' num2str(id_ch)]);
disp(['over ' num2str(cfg.time(1)) ' to ' num2str(cfg.time(2)) ' time period (sec)']);
disp(['with baseline being from ' num2str(cfg.baseline(1)) ' to ' num2str(cfg.baseline(2)) ' ms']);
disp(['p thresholds were ' num2str(cfg.erp_prob(1)) ' & ' num2str(cfg.erp_prob(2))]);
disp(['statistics compared ' num2str(cfg.trialinfo_dat0) ' against the sum of ' num2str(cfg.trialinfo_dat1) ' and ' num2str(cfg.trialinfo_dat2)]);

%% PLOT DISTRIBUTION
% figure; hold on
% period = [ ];
% ch = [];
% Dt = find(stats_data.time == period(:));
%
% data_plot = squeeze(mean(distribution(ch, Dt(:), :),1) ,2));
% steps_edges = linspace(min(data_plot), max(data_plot), 101);
% bar(steps_edges, histc(data_plot, steps_edges), 'k');
% dat=squeeze(mean(stats_data.erp(ch, Dt(:)), 2));
% line([dat dat], [0 max(histc(data_plot, steps_edges))], 'Color', 'r');

end