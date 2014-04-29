function [stats_data] = MrM_stats_tf_baseline(cfg, data)

% stats_tf_baseline performs statistics against the baseline period in the time frequency space.
% Data should be the output of ft_timefrequencyanalysis with output: 'fourier' (complex)
%
% Use as: stats_tf_baseline(cfg,data)
% The configuration can have the following parameters
% cfg.channel       = selected channel
% cfg.time          = time period over which the stats will be performed
%                     [sec sec] or 'all'
% cfg.freq          = selected freq
%                     [Hz Hz] or 'all'
% cfg.baseline      = selected baseline time period from onset
%                     [sec sec]
% cfg.iteration     = number of random permutation
% cfg.pow_prob      = statistic threshold for power (two tail like [2.5 97.5])
% cfg.pci_prob      = statistic threshold for pci (one tail like [95])
% cfg.distribution  = 'yes' save the random distribution
% cfg.mode          = define how the baseline "value" is defined and how the test is done :
%                     (This also change the way ppi stats are computed
%                     HAS TO BE CHECKED
%
%                     'random within and across (distribution)' **
%                     picked different random baseline time points for each trials
%                     then average across trials
%                     and looped to build surrogate "baseline" distribution
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
% cfg.ppi_Tref      = time point of reference for phase preservation index
%                     [sec]
% cfg.ppi_prob      = statistic threshold for ppi (one tail like [95])

%time index
switch isfield(cfg, 'time')
    case 1
        if length(cfg.time) == 2
            id_time =[find(round(data.time*1000) == round(cfg.time(1)*1000)) find(round(data.time*1000) == round(cfg.time(2)*1000))];
        elseif strcmp(cfg.time, 'all');
            id_time = [1 length(data.time)];
        else
            error('time period wrongly stated');
        end
    case 0
        error('no time of interest defined');
end
%Baseline index
switch isfield(cfg, 'baseline')
    case 1
        if isnumeric(cfg.baseline) && length(cfg.baseline) ==2
            if round(1000*cfg.baseline(1))>= round(1000*data.time(id_time(1))) && round(1000*cfg.baseline(2))<= round(1000*data.time(id_time(2)))
                id_Bl = [find(round(data.time*1000) == round(cfg.baseline(1)*1000)) find(round(data.time*1000) == round(cfg.baseline(2)*1000))];
            else
                error('baseline period should be included into the selected time period');
            end
        else
            error('baseline wrongly stated');
        end
    case 0
        error('no baseline period defined');
end
%frequency index
switch isfield(cfg, 'freq')
    case 1
        if length(cfg.freq) == 2
            id_freq = [find(data.freq == cfg.freq(1)) find(data.freq == cfg.freq(2))];
        elseif strcmp(cfg.freq, 'all');
            id_freq = [1 length(data.freq)];
        else
            error('frequencies wrongly stated');
        end
    case 0
        error('no frequency of interest defined');
end
%p values
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
        if length(cfg.pci_prob)==2
            error('p threshold for pci should contains 1 values (one-tail)');
        elseif length(cfg.pci_prob)==1
        else
            error('p threshold for pci wrongly stated');
            
        end
    case 0
        error('no p threshold defined for pci');
end
switch isfield(cfg, 'ppi_prob')
    case 1
        if length(cfg.ppi_prob)==2
            error('p threshold for PPI should contains 1 values (one-tail)');
        elseif length(cfg.ppi_prob)==1
        else
            error('p threshold for PPI wrongly stated');
        end
    case 0
end

%time point of reference for phase concentration index
switch isfield(cfg, 'ppi_Tref')
    case 1
        id_Tref = [find(round(data.time*1000) == round(cfg.ppi_Tref*1000))];
    case 0
        error('time frame of reference for phase preservation index (ppi) should be included in the selected period of interest');
end

%Prepare the output
stats_data.cfg       = cfg;
%data selection
stats_data.fourierspctrm = data.fourierspctrm(:, cfg.channel,id_freq(1):id_freq(2),id_time(1):id_time(2));
stats_data.freq          = data.freq(id_freq(1):id_freq(2));
stats_data.time          = data.time(id_time(1):id_time(2));
stats_data.dimord        = 'chan_freq_time';
stats_data.label         = {'1'};
clear data;

%Baseline index recomputed taking into account the period cropped
id_Bl = [find(round(stats_data.time*1000) == round(cfg.baseline(1)*1000)) find(round(stats_data.time*1000) == round(cfg.baseline(2)*1000))];
%dum variable
nb_trial         = size(stats_data.fourierspctrm,1);
nb_freq          = size(stats_data.fourierspctrm,3);
nb_tf            = size(stats_data.fourierspctrm,4);

%compute Power and PCI (ITC / pci)
stats_data.powspctrm     = mean(abs(squeeze(stats_data.fourierspctrm)) .^2,1);
stats_data.pcispctrm     = abs(mean(squeeze(stats_data.fourierspctrm)./abs(squeeze(stats_data.fourierspctrm)),1));
%phase preservation index
if isfield(cfg, 'ppi_Tref')
    % substraction in the complex plan
    % stats_data.ppispctrm = abs(mean( ...
    % repmat(squeeze(stats_data.fourierspctrm(:,:,:,id_Tref))./abs(squeeze(stats_data.fourierspctrm(:,:,:,id_Tref))), [ 1 1 nb_tf]) ...
    %      - squeeze(stats_data.fourierspctrm)./abs(squeeze(stats_data.fourierspctrm)) ...
    %       ,1));
    %'direct' angle difference
        stats_data.ppispctrm = abs(mean(exp(1i*( ...
            repmat(angle(squeeze(stats_data.fourierspctrm(:,:,:,id_Tref))), [ 1 1 nb_tf]) ...
            - angle(squeeze(stats_data.fourierspctrm)) ...
            )),1));
end
%set the mat for randomisation
pick_Rd_Bl  = NaN(nb_trial, nb_freq, nb_tf);
Rd_pow      = NaN(nb_freq,nb_tf, cfg.iteration);
Rd_pci      = NaN(nb_freq,nb_tf, cfg.iteration);
if isfield(cfg, 'ppi_Tref')
    Rd_ppi_mean = NaN(nb_freq,nb_tf, cfg.iteration);
end

if any(strfind(cfg.mode, 'distribution)'))
    for lp = 1:cfg.iteration;
        if strcmp(cfg.mode, 'random within and across (distribution)')
            % build the random draws index matrix
            pick_Rd_Index = round(rand(nb_trial,1) * (id_Bl(2) - id_Bl(1)) + id_Bl(1));
            % build the random draws value matrix
            for tr = 1:nb_trial
                pick_Rd_Bl(tr, :, :) = squeeze(stats_data.fourierspctrm(tr, 1, :, repmat(pick_Rd_Index(tr),[nb_tf 1])));
            end
        elseif strcmp(cfg.mode, 'random across (distribution)')
            % build the random draws index matrix
            pick_Rd_Index = round(rand(1) * (id_Bl(2) - id_Bl(1)) + id_Bl(1));
            % build the random draws value matrix
            pick_Rd_Bl    = squeeze(stats_data.fourierspctrm(:, 1, :, repmat(pick_Rd_Index,[nb_tf 1])));
        else
            error('you have to specify if the averaged baseline time points have the same latency or not within each iteration')
        end
        %compute the average Baseline value (across trial)
        Rd_pow(:, :, lp) = squeeze(mean(abs(pick_Rd_Bl) .^2, 1));
        Rd_pci(:, :, lp) = squeeze(abs(mean(pick_Rd_Bl./abs(pick_Rd_Bl),1)));
        if isfield(cfg, 'ppi_Tref')
            %prepare shuffled pool
            Tref_tmp    = repmat(stats_data.fourierspctrm(:,1,:,id_Tref), [1 1 1 nb_tf]);
            Rd_pool     = cat(1, squeeze(stats_data.fourierspctrm),squeeze(Tref_tmp));
            % permutations for ppi
            Rd_pool     = Rd_pool(randperm(size(Rd_pool,1)),:,:);
            % substraction in the complex plan
            % Rd_ppi_mean(:,:,lp)    = squeeze(abs(mean( ...
            %    repmat(Rd_pool(1+nb_trial:end,:,id_Tref)./abs(Rd_pool(1+nb_trial:end,:,id_Tref)),[1 1 nb_tf]) ...
            %         - Rd_pool(1:nb_trial,:,:)./abs(Rd_pool(1:nb_trial,:,:)) ...
            %   ,1)));
            %'direct' angle difference
            Rd_ppi_mean(:,:,lp) = abs(mean(exp(1i*( ...
                                  repmat(angle(squeeze(Rd_pool(1:nb_trial,:,id_Tref))), [ 1 1 nb_tf]) ...
                                       - angle(squeeze(Rd_pool(1+nb_trial:end,:,:))) ...
                                  )),1));
        end
        disp(['Randomisation #' num2str(lp) ' over ' num2str(cfg.iteration)]);
    end;
    
    %define percentile for each time point and frequencies
    thrshld_pow = prctile(Rd_pow, cfg.pow_prob , 3);
    thrshld_pci = prctile(Rd_pci, cfg.pci_prob , 3);
    if isfield(cfg, 'ppi_Tref')
        thrshld_ppi = prctile(Rd_ppi_mean, cfg.ppi_prob , 3);
    end
    
    %create stat masks
    stats_data.pow_mask(1,:, :) = (squeeze(stats_data.powspctrm) <= thrshld_pow(:,:,1) | squeeze(stats_data.powspctrm) >= thrshld_pow(:,:,2));
    stats_data.pci_mask(1,:, :) = squeeze(stats_data.pcispctrm) >= thrshld_pci(:,:,1);
    if isfield(cfg, 'ppi_Tref')
        stats_data.ppi_mask(1,:, :) = squeeze(stats_data.ppispctrm) >= thrshld_ppi(:,:,1);
    end
    %calculate the probabilities
    stats_data.pow_p = nan(1,nb_freq,nb_tf);
    stats_data.pci_p = nan(1,nb_freq,nb_tf);
    if isfield(cfg, 'ppi_Tref')
        stats_data.ppi_p = nan(1,nb_freq,nb_tf);
    end
    for fq = 1:nb_freq
        for tf = 1:nb_tf
            %define the tail depending on mean side for the power
            if stats_data.powspctrm(1,fq, tf) < mean(Rd_pow(fq,tf,:),3);
                stats_data.pow_p(1,fq,tf) = length(find(Rd_pow(fq,tf,:)<stats_data.powspctrm(1,fq, tf)))/size(Rd_pow,3);
            elseif stats_data.powspctrm(1,fq, tf) > mean(Rd_pow(fq,tf,:),3);
                stats_data.pow_p(1,fq,tf) = length(find(Rd_pow(fq,tf,:)>stats_data.powspctrm(1,fq, tf)))/size(Rd_pow,3);
            else
                stats_data.pow_p(1,fq,tf) = 0.5;
            end
            %one tail for the pci
            stats_data.pci_p(1,fq,tf) = length(find(Rd_pci(fq,tf,:)>stats_data.pcispctrm(1,fq, tf)))/size(Rd_pci,3);
            %one tail for the ppi
            if isfield(cfg, 'ppi_Tref')
                stats_data.ppi_p(1,fq,tf) = length(find(Rd_ppi_mean(fq,tf,:)>stats_data.ppispctrm(1,fq, tf)))/size(Rd_ppi_mean,3);
            end
        end
    end
    
elseif any(strfind(cfg.mode, 'paired'))
    for lp = 1:cfg.iteration;
        
        % build the random draws index matrix within baseline
        pick_Rd_Index = round(rand(nb_freq,nb_tf) * (id_Bl(2) - id_Bl(1)) + id_Bl(1));
        
        % build the random draws Baseline value matrix
        for fq = 1:nb_freq
            for tf = 1:nb_tf
                pick_Rd_Bl(:,fq,tf)    = squeeze(stats_data.fourierspctrm(:, 1, fq, pick_Rd_Index(fq,tf)));
            end
        end
       
        %compute the observed average of the differences for pow and pci
        stats_data.powspctrm        = mean(abs(squeeze(stats_data.fourierspctrm)) .^2,1);
        stats_data.pcispctrm        = abs(mean(squeeze(stats_data.fourierspctrm)./abs(squeeze(stats_data.fourierspctrm)),1));
        
        stats_data.diff_pow(:,:,lp) = stats_data.powspctrm - mean(abs(pick_Rd_Bl) .^2,1);
        stats_data.diff_pci(:,:,lp) = stats_data.pcispctrm - abs(mean(pick_Rd_Bl./abs(pick_Rd_Bl),1));
        
        % build the mixing paired matrices
        Rdperm_tr               = [1:nb_trial]'.*round(rand(nb_trial,1));
        Rdperm_tr(Rdperm_tr==0) = [];
        pool1_Rd                 = squeeze(stats_data.fourierspctrm);
        pool2_Rd                 = pick_Rd_Bl;
        pool1_Rd(Rdperm_tr,:,:)  = pick_Rd_Bl(Rdperm_tr,:,:);
        pool2_Rd(Rdperm_tr,:,:)  = stats_data.fourierspctrm(Rdperm_tr,1,:,:);
        
        %compute the average of the differences
        Rd_pow(:,:,lp)  = squeeze(mean(abs(pool1_Rd) .^2,1))-squeeze(mean(abs(pool2_Rd) .^2,1));
        Rd_pci(:,:,lp)  = squeeze(abs(mean(pool1_Rd ./abs(pool1_Rd),1))) - squeeze(abs(mean(pool2_Rd ./abs(pool2_Rd),1)));
 
        if isfield(cfg, 'ppi_Tref')        
        % build the mixing paired matrices
        Rdperm_tr               = [1:nb_trial]'.*round(rand(nb_trial,1));
        Rdperm_tr(Rdperm_tr==0) = [];
        pool1_Rd                 = stats_data.fourierspctrm(:,1,:,:);
        pool2_Rd                 = repmat(stats_data.fourierspctrm(:,1,:,id_Tref), [ 1 1 nb_tf]);
        pool1_Rd(Rdperm_tr,:,:)  = repmat(stats_data.fourierspctrm(Rdperm_tr,1,:,id_Tref), [ 1 1 nb_tf]);
        pool2_Rd(Rdperm_tr,:,:)  = stats_data.fourierspctrm(Rdperm_tr,1,:,:);
        
        %compute the average of the differences
        Rd_ppi(:,:,lp) = abs(mean(exp(1i*(angle(pool1_Rd) - angle(pool2_Rd))),1));
        end
        
        disp(['Randomisation #' num2str(lp) ' over ' num2str(cfg.iteration)]);
    end
    
    %define percentile for each time point and frequencies
    thrshld_pow = prctile(Rd_pow, cfg.pow_prob , 3);
    thrshld_pci = prctile(Rd_pci, cfg.pci_prob , 3);
    if isfield(cfg, 'ppi_Tref')
        thrshld_ppi = prctile(Rd_ppi, cfg.ppi_prob , 3);
    end
    
    %compare the actual data with threshold percentile
    %create stat mask for pow and PCI
    stats_data.pow_mask(1,:, :) = (squeeze(mean(stats_data.diff_pow,3)) <= thrshld_pow(:,:,1) | squeeze(mean(stats_data.diff_pow,3)) >= thrshld_pow(:,:,2));
    stats_data.pci_mask(1,:, :) = squeeze(mean(stats_data.diff_pci,3)) >= thrshld_pci(:,:,1);
    if isfield(cfg, 'ppi_Tref')
        stats_data.ppi_mask(1,:, :) = squeeze(stats_data.ppispctrm) >= thrshld_ppi(:,:,1);
    end
   
    %compute probabilities
    stats_data.pow_p = nan(1,id_freq(2),id_time(2));
    stats_data.pci_p = nan(1,id_freq(2),id_time(2));
    if isfield(cfg, 'ppi_Tref')
        stats_data.ppi_p = nan(1,nb_freq,nb_tf);
    end    
    for fq = 1:length(stats_data.freq)
        for tf = 1:length(stats_data.time)
            %define the tail depending on mean side for the power
            if mean(stats_data.diff_pow(fq, tf,:),3) < mean(Rd_pow(fq,tf,:),3);
                stats_data.pow_p(1,fq,tf) = length(find(Rd_pow(fq,tf,:)<mean(stats_data.diff_pow(fq, tf,:),3)))/size(Rd_pow,3);
            elseif mean(stats_data.diff_pow(fq, tf,:),3) > mean(Rd_pow(fq,tf,:),3);
                stats_data.pow_p(1,fq,tf) = length(find(Rd_pow(fq,tf,:)>mean(stats_data.diff_pow(fq, tf,:),3)))/size(Rd_pow,3);
            end
            %for pci
            stats_data.pci_p(1,fq,tf) = length(find(Rd_pci(fq,tf,:)>mean(stats_data.diff_pci(fq, tf,:),3)))/size(Rd_pci,3);
            %one tail for the ppi
            if isfield(cfg, 'ppi_Tref')
            stats_data.ppi_p(1,fq,tf) = length(find(Rd_ppi(fq,tf,:)>stats_data.ppispctrm(1,fq, tf)))/size(Rd_ppi,3);
            end
        end
    end
end

switch isfield(cfg, 'distribution')
    case 1
        if strcmp(cfg.distribution, 'yes');
            stats_data.distribution.pow = Rd_pow;
            stats_data.distribution.pci = Rd_pci;
            if isfield(cfg, 'ppi_Tref')
                stats_data.distribution.ppi = Rd_ppi_mean;
            end
            disp('random distributions were saved');
        elseif strcmp(cfg.distribution, 'no');
            disp('random distributions were not saved');
        end
    case 0
        disp('random distributions were not saved');
end

disp('Permutations finished!');
disp(['Channel tested was ' num2str(cfg.channel)]);
% disp(['over ' num2str(data.time(id_time(1))) ' to ' num2str(data.time(id_time(2))) ' time period (sec)']);
% disp(['from ' num2str(data.freq(id_freq(1))) ' to ' num2str(data.freq(id_freq(2))) ' Hz']);
% disp(['with baseline being from ' num2str(cfg.baseline(1)) ' to ' num2str(cfg.baseline(2)) ' ms']);
% disp(['p thresholds were ' num2str(cfg.pow_prob(1)) ' & ' num2str(cfg.pow_prob(2)) ' for power and ' num2str(cfg.pci_prob) ' for pci']);
% if isfield(cfg, 'ppi_Tref')
%    disp(['ppi stats computed']);
% end

end