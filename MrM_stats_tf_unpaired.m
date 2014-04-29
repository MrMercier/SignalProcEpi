function [stats_data] = MrM_stats_tf_unpaired(cfg, data1, data2)

% stats_tf_unpaired performs statistics between two conditions in the time frequency space.
% (at the single trial level using nonparametric 'permutation' tests)
% Data should be the output of ft_timefrequencyanalysis with output: 'fourier' (complex)
% dist used in Fiebelkorn et al., 2012
% pbi used in Busch et al., 2009
%
% the mode specifies if the difference is computed after or before taking the absolute value:
% absdiff / not sensitive to phase difference
% diffabs / sensitive to phase difference
%
% Use as: stats_tf_baseline(cfg,data1, data2)
% The configuration can have the following parameters
% cfg.channel       = selected channel
% cfg.time          = time period over which the stats will be performed
%                     [sec sec] or 'all'
% cfg.freq          = selected freq
%                     [Hz Hz] or 'all'
% cfg.iteration     = number of random permutation
% cfg.pow_prob      = statistic threshold for power (two tail, like [2.5 97.5])
% cfg.pci_prob      = statistic threshold for pci (one or two tail depending on the mode)
% cfg.dist_prob     = statistic threshold for pci distance
% cfg.pbi_prob      = statistic threshold for phase bifurcation index
% cfg.distribution  = 'yes' save the random distributions
% cfg.mode          = 'absdiff' or 'diffabs' or 'complexdiff'
% cfg.baseline      = 'no' or time period used to baseline correct the data
%                     [sec sec] or 'all' Only applied on power
% cfg.baselinetype  = 'relative' or 'absolute' or 'relchange' ('complexdiff' not  ...)

%check p values regarding mode set
switch isfield(cfg, 'pow_prob')
    case 1
        if length(cfg.pow_prob)==2
        elseif length(cfg.pow_prob)==1
            error('p threshold for power should contains 2 values (two-tail)');
        else
            error('p threshold for power wrongly stated');
        end
    case 0
        error('no p threshold defined for power');
end
switch isfield(cfg, 'pci_prob')
    case 1
        if length(cfg.pci_prob)==2 && strcmp(cfg.mode,'diffabs')
        elseif length(cfg.pci_prob)==2 && strcmp(cfg.mode,'absdiff')
            error('p threshold for pci should contains 1 values (one-tail) because absdiff mode was specified');
        elseif length(cfg.pci_prob)==1 &&strcmp(cfg.mode,'absdiff')
        elseif length(cfg.pci_prob)==1 && strcmp(cfg.mode, 'diffabs')
            error('p threshold for pci should contains 2 values (two-tails) because diffabs mode was specified');
        else
            error('p threshold for pci wrongly stated');
        end
    case 0
        error('no p threshold defined for pci');
end
%time index
switch isfield(cfg, 'time')
    case 1
        if length(cfg.time) == 2 && isnumeric(cfg.time)
            id_time_dat1 =[find(round(data1.time*1000) == round(cfg.time(1)*1000)) find(round(data1.time*1000) == round(cfg.time(2)*1000))];
            id_time_dat2 =[find(round(data2.time*1000) == round(cfg.time(1)*1000)) find(round(data2.time*1000) == round(cfg.time(2)*1000))];
        elseif strcmp (cfg.time, 'all')
            id_time_dat1 =[1 length(data1.time)];
            id_time_dat2 =[1 length(data2.time)];
        else
            error('time period wrongly stated');
        end
    case 0
        error('no time of interest defined');
end
%frequency index
switch isfield(cfg, 'freq')
    case 1
        if length(cfg.freq) == 2 && isnumeric(cfg.time)
            id_freq_dat1 = [find(data1.freq == cfg.freq(1)) find(data1.freq == cfg.freq(2))];
            id_freq_dat2 = [find(data2.freq == cfg.freq(1)) find(data2.freq == cfg.freq(2))];
        elseif strcmp (cfg.freq, 'all')
            id_freq_dat1 =[1 length(data1.freq)];
            id_freq_dat2 =[1 length(data2.freq)];
        else
            error('frequencies wrongly stated');
        end
    case 0
        error('no frequencies of interest defined');
end
%check that time and freq dimensions of the conditions are identical
if ~isequal(data1.time(id_time_dat1(1):id_time_dat1(2)), data2.time(id_time_dat2(1):id_time_dat2(2)));
    error('selected time period do not match between condition');
end
if ~isequal(data1.freq(id_freq_dat1(1):id_freq_dat1(2)), data2.freq(id_freq_dat2(1):id_freq_dat2(2)));
    error('selected frequencies do not match between condition');
end
%check is there is NaN in the selected data matrix
if any(any(any(any(isnan(data1.fourierspctrm(:,cfg.channel,id_freq_dat1(1):id_freq_dat1(2),id_time_dat1(1):id_time_dat1(2))))))) ...
        && any(any(any(any(isnan(data2.fourierspctrm(:,cfg.channel,id_freq_dat2(1):id_freq_dat2(2),id_time_dat2(1):id_time_dat2(2)))))))
    error('There is NaN in the selected data');
end

%Prepare the output
stats_data.cfg     = cfg;
stats_data.data1   = inputname(2);
stats_data.data2   = inputname(3);
%data selection
stats_data.freq    = data1.freq(id_freq_dat1(1):id_freq_dat1(2));
stats_data.time    = data1.time(id_time_dat1(1):id_time_dat1(2));
stats_data.dimord  = 'chan_freq_time';
stats_data.label   = {'1'};

%baseline correction index
if strcmp(cfg.baseline, 'no')
else
    if isnumeric(cfg.baseline) && length(cfg.baseline) ==2
    id_Bl = [find(round(stats_data.time*1000) == round(cfg.baseline(1)*1000)) find(round(stats_data.time*1000) == round(cfg.baseline(2)*1000))];
    elseif strcmp(cfg.baseline, 'all')
    id_Bl = [stats_data.time(1) stats_data.time(end)];
    else
    error('baseline not specified correctly');
    end;
    if id_Bl(1)< id_time_dat1(1) || id_Bl(2) > id_time_dat2(2)
    error('baseline period should be included into the selected time period');
    end     
end

%compute average Power
stats_data.dat1_pow = squeeze(mean(abs(data1.fourierspctrm) .^2,1));
stats_data.dat2_pow = squeeze(mean(abs(data2.fourierspctrm) .^2,1));
%baseline (only done on Power)
if ~strcmp(cfg.baseline, 'no')
    %compute the average power during the baseline period
    dat1_pow_Bl = nanmean(stats_data.dat1_pow(:,:,id_Bl(1):id_Bl(2)),3);
    dat2_pow_Bl = nanmean(stats_data.dat2_pow(:,:,id_Bl(1):id_Bl(2)),3);
    %compute baseline correction
    if strcmp(cfg.baselinetype,'absolute')
        stats_data.dat1_pow = stats_data.dat1_pow - repmat(dat1_pow_Bl, [1 1 size(stats_data.dat1_pow,3)]);
        stats_data.dat2_pow = stats_data.dat2_pow - repmat(dat2_pow_Bl, [1 1 size(stats_data.dat2_pow,3)]);
    elseif strcmp(cfg.baselinetype,'relative')
        stats_data.dat1_pow = stats_data.dat1_pow ./ repmat(dat1_pow_Bl, [1 1 size(stats_data.dat1_pow,3)]);
        stats_data.dat2_pow = stats_data.dat2_pow ./ repmat(dat2_pow_Bl, [1 1 size(stats_data.dat2_pow,3)]);
    elseif strcmp(cfg.baselinetype,'relchange')
        stats_data.dat1_pow = (stats_data.dat1_pow - repmat(dat1_pow_Bl, [1 1 size(stats_data.dat1_pow,3)]))./ repmat(dat1_pow_Bl, [1 1 size(stats_data.dat1_pow,3)]);
        stats_data.dat2_pow = (stats_data.dat2_pow - repmat(dat2_pow_Bl, [1 1 size(stats_data.dat2_pow,3)]))./ repmat(dat2_pow_Bl, [1 1 size(stats_data.dat2_pow,3)]);    
    end
end
%compute Power difference
stats_data.diff_pow = stats_data.dat1_pow - stats_data.dat2_pow;

%compute pci and the difference
switch cfg.mode
    case 'diffabs'
        % first compute the absolute, then take the difference
        % (difference of the norm not sensitive to phase differences)
        stats_data.dat1_pci = squeeze(abs(mean(data1.fourierspctrm ./abs(data1.fourierspctrm),1)));
        stats_data.dat2_pci = squeeze(abs(mean(data2.fourierspctrm ./abs(data2.fourierspctrm),1)));
        stats_data.diff_pci = stats_data.dat1_pci - stats_data.dat2_pci;
    case 'absdiff'
        % first compute the difference, then take the absolute
        % (difference of the normalise mean angle is sensitive to phase differences)
        stats_data.dat1_pci = squeeze(mean(data1.fourierspctrm ./abs(data1.fourierspctrm),1));
        stats_data.dat2_pci = squeeze(mean(data2.fourierspctrm ./abs(data2.fourierspctrm),1));
        stats_data.diff_pci = abs(stats_data.dat1_pci - stats_data.dat2_pci);
            % still compute PCI normally to be able to access it later
            stats_data.dat1_pci = abs(stats_data.dat1_pci);
            stats_data.dat2_pci = abs(stats_data.dat2_pci);
    otherwise
        error('incorrect mode specification');
end

%compute pbi
%should be corrected by the number of trials
if isfield(cfg,'pbi_prob') && length(cfg.pbi_prob)==2 && isnumeric(cfg.pbi_prob)
    alldata.fourierspctrm = cat(1,data1.fourierspctrm,data2.fourierspctrm);
    alldata.pci           = squeeze(abs(mean(alldata.fourierspctrm ./abs(alldata.fourierspctrm),1)));
    stats_data.pbi        = (stats_data.dat1_pci - alldata.pci) .* (stats_data.dat2_pci - alldata.pci);
elseif isfield(cfg,'pbi_prob')
    error('phase bifurcation index p-value wrongly stated');
end

%compute the distance between the two vectors
if isfield(cfg,'dist_prob') && length(cfg.pbi_prob)==2 && isnumeric(cfg.pbi_prob)
        stats_data.dat1_pci = squeeze(mean(data1.fourierspctrm ./abs(data1.fourierspctrm),1));
        stats_data.dat2_pci = squeeze(mean(data2.fourierspctrm ./abs(data2.fourierspctrm),1));
        stats_data.dist_pci = sqrt((real(stats_data.dat1_pci)-real(stats_data.dat2_pci)).^2 + (imag(stats_data.dat1_pci)-imag(stats_data.dat2_pci)).^2);
            % still compute PCI normally to be able to access it later
            stats_data.dat1_pci = abs(stats_data.dat1_pci);
            stats_data.dat2_pci = abs(stats_data.dat2_pci);
elseif isfield(cfg,'dist_prob')
    error('phase bifurcation index p-value wrongly stated');
end

% get the selected subset (channel, frequencies and time points)
stats_data.dat1_pow = stats_data.dat1_pow(cfg.channel,id_freq_dat1(1):id_freq_dat1(2),id_time_dat1(1):id_time_dat1(2));
stats_data.dat2_pow = stats_data.dat2_pow(cfg.channel,id_freq_dat2(1):id_freq_dat2(2),id_time_dat2(1):id_time_dat2(2));
stats_data.diff_pow = stats_data.diff_pow(cfg.channel,id_freq_dat1(1):id_freq_dat1(2),id_time_dat1(1):id_time_dat1(2));

stats_data.dat1_pci = stats_data.dat1_pci(cfg.channel,id_freq_dat1(1):id_freq_dat1(2),id_time_dat1(1):id_time_dat1(2));
stats_data.dat2_pci = stats_data.dat2_pci(cfg.channel,id_freq_dat2(1):id_freq_dat2(2),id_time_dat2(1):id_time_dat2(2));
stats_data.diff_pci = stats_data.diff_pci(cfg.channel,id_freq_dat1(1):id_freq_dat1(2),id_time_dat1(1):id_time_dat1(2));

dat1_fourierspctrm  = data1.fourierspctrm(:,cfg.channel,id_freq_dat1(1):id_freq_dat1(2),id_time_dat1(1):id_time_dat1(2));
dat2_fourierspctrm  = data2.fourierspctrm(:,cfg.channel,id_freq_dat2(1):id_freq_dat2(2),id_time_dat2(1):id_time_dat2(2));

if isfield(cfg,'pbi_prob')
    stats_data.pbi = stats_data.pbi(cfg.channel,id_freq_dat1(1):id_freq_dat1(2),id_time_dat1(1):id_time_dat1(2));
end

if isfield(cfg,'dist_prob')
    stats_data.dist_pci = stats_data.dist_pci(cfg.channel,id_freq_dat1(1):id_freq_dat1(2),id_time_dat1(1):id_time_dat1(2));    
end

%dum variable
nb_trial_dat1 = size(dat1_fourierspctrm,1);

clear data1;
clear data2;
% clear dat1_fourierspctrm;
% clear dat2_fourierspctrm;
clear dat1_pow_Bl;
clear dat2_pow_Bl;

%start the permutation process
for lp = 1:cfg.iteration;
    %build a matrix with both condition
    Rd_pool = vertcat(squeeze(dat1_fourierspctrm),squeeze(dat2_fourierspctrm));    
    %shuffle the trials
    Rd_pool = Rd_pool(randperm(size(Rd_pool,1)),:,:);
    %compute average Power
    dat1_Rd_pow = squeeze(mean(abs(Rd_pool(1:nb_trial_dat1,:,:)) .^2,1));
    dat2_Rd_pow = squeeze(mean(abs(Rd_pool(1+nb_trial_dat1:end,:,:)) .^2,1));
    %baseline Power
    if ~strcmp(cfg.baseline, 'no')
        %compute the average power during the baseline period
        dat1_Rd_pow_Bl = squeeze(nanmean(dat1_Rd_pow(:,id_Bl(1):id_Bl(2)),2));
        dat2_Rd_pow_Bl = squeeze(nanmean(dat2_Rd_pow(:,id_Bl(1):id_Bl(2)),2));
        %compute baseline correction
        if strcmp(cfg.baselinetype,'absolute')
            dat1_Rd_pow = dat1_Rd_pow - repmat(dat1_Rd_pow_Bl, [1 size(dat1_Rd_pow,2)]);
            dat2_Rd_pow = dat2_Rd_pow - repmat(dat2_Rd_pow_Bl, [1 size(dat2_Rd_pow,2)]);
        elseif strcmp(cfg.baselinetype,'relative')
            dat1_Rd_pow = dat1_Rd_pow ./ repmat(dat1_Rd_pow_Bl, [1 size(dat1_Rd_pow,2)]);
            dat2_Rd_pow = dat2_Rd_pow ./ repmat(dat2_Rd_pow_Bl, [1 size(dat2_Rd_pow,2)]);
        elseif strcmp(cfg.baselinetype,'relchange')
            dat1_Rd_pow = (dat1_Rd_pow - repmat(dat1_Rd_pow_Bl, [1 size(dat1_Rd_pow,2)]))./ repmat(dat1_Rd_pow_Bl, [1 size(dat1_Rd_pow,2)]);
            dat2_Rd_pow = (dat2_Rd_pow - repmat(dat2_Rd_pow_Bl, [1 size(dat2_Rd_pow,2)]))./ repmat(dat2_Rd_pow_Bl, [1 size(dat2_Rd_pow,2)]);            
        end
    end
    %compute Power difference
    diff_Rd_pow(:,:,lp) = dat1_Rd_pow - dat2_Rd_pow;
    %compute PCI (ITC/PLF) differences
    switch cfg.mode
        case 'diffabs'
            % first compute the absolute, then take the difference
            % (not sensitive to phase differences)
            dat1_Rd_pci        = abs(mean(Rd_pool(1:nb_trial_dat1,:,:) ./abs(Rd_pool(1:nb_trial_dat1,:,:)),1));
            dat2_Rd_pci        = abs(mean(Rd_pool(1+nb_trial_dat1:end,:,:) ./abs(Rd_pool(1+nb_trial_dat1:end,:,:)),1));
            diff_Rd_pci(:,:,lp) = dat1_Rd_pci - dat2_Rd_pci;
            
        case 'absdiff'
            % first compute the difference, then take the absolute
            % (sensitive to phase differences)
            dat1_Rd_pci        = mean(Rd_pool(1:nb_trial_dat1,:,:) ./abs(Rd_pool(1:nb_trial_dat1,:,:)),1);
            dat2_Rd_pci        = mean(Rd_pool(1+nb_trial_dat1:end,:,:) ./abs(Rd_pool(1+nb_trial_dat1:end,:,:)),1);
            diff_Rd_pci(:,:,lp) = abs(dat1_Rd_pci - dat2_Rd_pci);
        otherwise
            error('incorrect mode specification');
    end
    %compute phase bifurcation index
    if isfield(cfg,'pbi_prob')
        dat1_Rd_pci        = abs(mean(Rd_pool(1:nb_trial_dat1,:,:) ./abs(Rd_pool(1:nb_trial_dat1,:,:)),1));
        dat2_Rd_pci        = abs(mean(Rd_pool(1+nb_trial_dat1:end,:,:) ./abs(Rd_pool(1+nb_trial_dat1:end,:,:)),1));
        all_Rd_pci         = stats_data.pbi;
        Rd_pbi(:,:,lp)     = (dat1_Rd_pci - all_Rd_pci) .* (dat2_Rd_pci - all_Rd_pci);
    end
    
    if isfield(cfg,'dist_prob')
        dat1_Rd_pci         = mean(Rd_pool(1:nb_trial_dat1,:,:) ./abs(Rd_pool(1:nb_trial_dat1,:,:)),1);
        dat2_Rd_pci         = mean(Rd_pool(1+nb_trial_dat1:end,:,:) ./abs(Rd_pool(1+nb_trial_dat1:end,:,:)),1);
        Rd_dist_pci(:,:,lp) = sqrt((real(dat1_Rd_pci)-real(dat2_Rd_pci)).^2 + (imag(dat1_Rd_pci)-imag(dat2_Rd_pci)).^2);
    end
    
    disp(['Randomization #' num2str(lp) ' over ' num2str(cfg.iteration)]);
end

%define percentile for each time point and frequencies
thrshld_pow = prctile(diff_Rd_pow, cfg.pow_prob , 3);
thrshld_pci = prctile(diff_Rd_pci, cfg.pci_prob , 3);
if isfield(cfg,'pbi_prob')
    thrshld_pbi = prctile(Rd_pbi, cfg.pbi_prob , 3);
end
if isfield(cfg,'dist_prob')
    thrshld_dist_pci = prctile(Rd_dist_pci, cfg.dist_prob , 3);
end

%compare the actual data with threshold percentile
%create stat mask for pow
stats_data.pow_mask(1,:, :) = (squeeze(stats_data.diff_pow) <= thrshld_pow(:,:,1) | squeeze(stats_data.diff_pow) >= thrshld_pow(:,:,2));
%create stat mask for pci
if length(cfg.pci_prob)==2
    stats_data.pci_mask(1,:, :) = (squeeze(stats_data.diff_pci) <= thrshld_pci(:,:,1) | squeeze(stats_data.diff_pci) >= thrshld_pci(:,:,2));
elseif length(cfg.pci_prob)==1
    stats_data.pci_mask(1,:, :) = squeeze(stats_data.diff_pci) >= thrshld_pci(:,:,1);
end
%create stat mask for pbi
if isfield(cfg,'pbi_prob')
   stats_data.pbi_mask(1,:, :) = (squeeze(stats_data.pbi) <= thrshld_pbi(:,:,1) | squeeze(stats_data.pbi) >= thrshld_pbi(:,:,2)); 
end
%create stat mask for pci distance
if isfield(cfg,'dist_prob')
   stats_data.dist_pci_mask(1,:, :) = (squeeze(stats_data.dist_pci) <= thrshld_dist_pci(:,:,1) | squeeze(stats_data.dist_pci) >= thrshld_dist_pci(:,:,2)); 
end
%compute probabilities
stats_data.pow_p = nan(1,id_freq_dat1(2),id_time_dat1(2));
stats_data.pci_p = nan(1,id_freq_dat1(2),id_time_dat1(2));
for fq = 1:length(stats_data.freq)
    for tf = 1:length(stats_data.time)
        %define the tail depending on mean side for the power
        if stats_data.diff_pow(1,fq, tf) < mean(diff_Rd_pow(fq,tf,:),3);
            stats_data.pow_p(1,fq,tf) = length(find(diff_Rd_pow(fq,tf,:)<stats_data.diff_pow(1,fq, tf)))/size(diff_Rd_pow,3);
        elseif stats_data.diff_pow(1,fq, tf) > mean(diff_Rd_pow(fq,tf,:),3);
            stats_data.pow_p(1,fq,tf) = length(find(diff_Rd_pow(fq,tf,:)>stats_data.diff_pow(1,fq, tf)))/size(diff_Rd_pow,3);
        end
        %define the tail for pci
        if length(cfg.pci_prob)==2
            if stats_data.diff_pci(1,fq, tf) < mean(diff_Rd_pci(fq,tf,:),3);
                stats_data.pci_p(1,fq,tf) = length(find(diff_Rd_pci(fq,tf,:)<stats_data.diff_pci(1,fq, tf)))/size(diff_Rd_pci,3);
            elseif stats_data.diff_pci(1,fq, tf) > mean(diff_Rd_pci(fq,tf,:),3);
                stats_data.pci_p(1,fq,tf) = length(find(diff_Rd_pci(fq,tf,:)>stats_data.diff_pci(1,fq, tf)))/size(diff_Rd_pci,3);
            end
        elseif length(cfg.pci_prob)==1
            stats_data.pci_p(1,fq,tf) = length(find(diff_Rd_pci(fq,tf,:)>stats_data.diff_pci(1,fq, tf)))/size(diff_Rd_pci,3);
        end
        %define the tail for pbi
        if isfield(cfg,'pbi_prob')
            %define the tail depending on mean side of the observed phase bifurcation index
            if stats_data.pbi(1,fq, tf) < mean(Rd_pbi(fq,tf,:),3);
               stats_data.pbi_p(1,fq,tf) = length(find(Rd_pbi(fq,tf,:)<stats_data.pbi(1,fq, tf)))/size(Rd_pbi,3);
            elseif stats_data.pbi(1,fq, tf) > mean(Rd_pbi(fq,tf,:),3);
               stats_data.pbi_p(1,fq,tf) = length(find(Rd_pbi(fq,tf,:)>stats_data.pbi(1,fq, tf)))/size(Rd_pbi,3);
            end    
        end
        %define the tail for pci distance
        if isfield(cfg,'dist_prob')
            %define the tail depending on mean side of the observed pci distance
            if stats_data.dist_pci(1,fq, tf) < mean(Rd_dist_pci(fq,tf,:),3);
               stats_data.dist_pci_p(1,fq,tf) = length(find(Rd_dist_pci(fq,tf,:)<stats_data.dist_pci(1,fq, tf)))/size(Rd_dist_pci,3);
            elseif stats_data.dist_pci(1,fq, tf) > mean(Rd_dist_pci(fq,tf,:),3);
               stats_data.dist_pci_p(1,fq,tf) = length(find(Rd_dist_pci(fq,tf,:)>stats_data.dist_pci(1,fq, tf)))/size(Rd_dist_pci,3);
            end    
        end
    end
end

switch isfield(cfg, 'distribution')
    case 1
        if strcmp(cfg.distribution, 'yes');
            stats_data.distribution.pow = diff_Rd_pow;
            stats_data.distribution.pci = diff_Rd_pci;
            if isfield(cfg,'pbi_prob')
                stats_data.distribution.pbi = Rd_pbi;
            end
            if isfield(cfg,'dist_prob')
                stats_data.distribution.dist = Rd_dist_pci;
            end
            disp('random distributions were saved');
        elseif strcmp(cfg.distribution, 'no');
            disp('random distributions were not saved');
        end
    case 0
        disp('random distributions were not saved');
end

disp('Unpaired Permutation test performed');

end

% %% check distribution
% figure; hold on
% fq = 16;
% tf = 151;
% steps_edges = linspace(min(Rd_pbi(fq,tf,:)), max(Rd_pbi(fq,tf,:)),101);
% bar(steps_edges, squeeze(histc(Rd_pbi(fq,tf,:), steps_edges)), 'k');
% Obs = 	stats_data.pbi(1,fq,tf);
% line([Obs Obs], [0 max(histc(Rd_pbi(fq,tf,:), steps_edges))], 'Color', 'r');